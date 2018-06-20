/*
 * Copyright (C) 2011-2015 Karlsruhe Institute of Technology
 *
 * This file is part of Ufo.
 *
 * This library is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "ufo-complex-mult-task.h"
#include "stdio.h"


struct _UfoComplexMultTaskPrivate {
    cl_kernel k_mult;
    UfoProfiler *profiler;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoComplexMultTask, ufo_complex_mult_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_COMPLEX_MULT_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_COMPLEX_MULT_TASK, UfoComplexMultTaskPrivate))

UfoNode *
ufo_complex_mult_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_COMPLEX_MULT_TASK, NULL));
}

static void
ufo_complex_mult_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoComplexMultTaskPrivate *priv;

    priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (task);
    priv->profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));
    priv->k_mult = ufo_resources_get_kernel (resources, "mult.cl", "mult3d", NULL, error);

    if (priv->k_mult)
        UFO_RESOURCES_CHECK_CLERR (clRetainKernel (priv->k_mult));
}

static void
ufo_complex_mult_task_get_requisition (UfoTask *task,
                                       UfoBuffer **inputs,
                                       UfoRequisition *requisition,
                                       GError **error)
{
    ufo_buffer_get_requisition(inputs[1], requisition);

    UfoRequisition req0;
    ufo_buffer_get_requisition (inputs[0], &req0);
    if (req0.n_dims == 3) {
        requisition->dims[2] *= req0.dims[2];
    }
}

static guint
ufo_complex_mult_task_get_num_inputs (UfoTask *task)
{
    return 2;
}

static guint
ufo_complex_mult_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return -1;
}

static UfoTaskMode
ufo_complex_mult_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static gboolean
ufo_complex_mult_task_process (UfoTask *task,
                               UfoBuffer **inputs,
                               UfoBuffer *output,
                               UfoRequisition *requisition)
{
    UfoComplexMultTaskPrivate *priv;
    UfoGpuNode *node;
    UfoProfiler *profiler;
    UfoRequisition req0;
    UfoRequisition req1;
    cl_command_queue cmd_queue;

    priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));
    cmd_queue = ufo_gpu_node_get_cmd_queue (node);
    profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));

    ufo_buffer_get_requisition (inputs[0], &req0);
    ufo_buffer_get_requisition (inputs[1], &req1);

    g_assert_cmpint(req0.dims[0], ==, req1.dims[0]);
    g_assert_cmpint(req0.dims[1], ==, req1.dims[1]);

    cl_mem in_mem0 = ufo_buffer_get_device_array (inputs[0], cmd_queue);
    cl_mem in_mem1 = ufo_buffer_get_device_array (inputs[1], cmd_queue);
    cl_mem out_mem = ufo_buffer_get_device_array (output, cmd_queue);

    // input/output are complex values
    gsize global_work_size[3];
    global_work_size[0] = req1.dims[0] / 2;
    global_work_size[1] = req1.dims[1];
    global_work_size[2] = req1.n_dims == 3 ? req1.dims[2] : 1;

    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->k_mult, 0, sizeof (cl_mem), (gpointer) &(in_mem0)));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->k_mult, 1, sizeof (cl_mem), (gpointer) &(in_mem1)));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->k_mult, 2, sizeof (cl_mem), (gpointer) &(out_mem)));

    ufo_profiler_call (profiler, cmd_queue, priv->k_mult, 3, global_work_size, NULL);

    return TRUE;
}

static void
ufo_complex_mult_task_finalize (GObject *object)
{
    UfoComplexMultTaskPrivate *priv;

    priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (object);

    if (priv->k_mult)
        UFO_RESOURCES_CHECK_CLERR (clReleaseKernel (priv->k_mult));

    G_OBJECT_CLASS (ufo_complex_mult_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_complex_mult_task_setup;
    iface->get_num_inputs = ufo_complex_mult_task_get_num_inputs;
    iface->get_num_dimensions = ufo_complex_mult_task_get_num_dimensions;
    iface->get_mode = ufo_complex_mult_task_get_mode;
    iface->get_requisition = ufo_complex_mult_task_get_requisition;
    iface->process = ufo_complex_mult_task_process;
}

static void
ufo_complex_mult_task_class_init (UfoComplexMultTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->finalize = ufo_complex_mult_task_finalize;

    g_type_class_add_private (oclass, sizeof(UfoComplexMultTaskPrivate));
}

static void
ufo_complex_mult_task_init(UfoComplexMultTask *self)
{
    self->priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE(self);
}
