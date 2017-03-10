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
    cl_context context;
    cl_command_queue cmd_queue;
    cl_kernel k_mult;
    UfoProfiler *profiler;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoComplexMultTask, ufo_complex_mult_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_COMPLEX_MULT_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_COMPLEX_MULT_TASK, UfoComplexMultTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

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
    UfoGpuNode *node;

    priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));
    
    priv->profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));
    priv->context = ufo_resources_get_context (resources);
    priv->cmd_queue = ufo_gpu_node_get_cmd_queue (node);

    priv->k_mult  = ufo_resources_get_kernel (resources, "mult.cl", "mult3d", error);

    UFO_RESOURCES_CHECK_CLERR (clRetainContext (priv->context));
}

static void
ufo_complex_mult_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition(inputs[1], requisition);

    UfoRequisition req0;
    ufo_buffer_get_requisition(inputs[0], &req0);
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
    UfoRequisition req0;
    UfoRequisition req1;
    cl_event event;

    priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (task);

    ufo_buffer_get_requisition(inputs[0], &req0);
    ufo_buffer_get_requisition(inputs[1], &req1);

    g_assert_cmpint(req0.dims[0], ==, req1.dims[0]);
    g_assert_cmpint(req0.dims[1], ==, req1.dims[1]);

    cl_mem in_mem0 = ufo_buffer_get_device_array(inputs[0], priv->cmd_queue);
    cl_mem in_mem1 = ufo_buffer_get_device_array(inputs[1], priv->cmd_queue);
    cl_mem out_mem = ufo_buffer_get_device_array(output, priv->cmd_queue);

    // input/output are complex values
    gsize global_work_size[3];
    global_work_size[0] = req1.dims[0]/2;
    global_work_size[1] = req1.dims[1];
    global_work_size[2] = req1.n_dims == 3 ? req1.dims[2] : 1;

    UFO_RESOURCES_CHECK_CLERR ( clSetKernelArg ( 
            priv->k_mult, 0, sizeof (cl_mem), (gpointer) &( in_mem0 )));
    UFO_RESOURCES_CHECK_CLERR ( clSetKernelArg ( 
            priv->k_mult, 1, sizeof (cl_mem), (gpointer) &( in_mem1 )));
    UFO_RESOURCES_CHECK_CLERR ( clSetKernelArg ( 
            priv->k_mult, 2, sizeof (cl_mem), (gpointer) &( out_mem )));

    UFO_RESOURCES_CHECK_CLERR (
           clEnqueueNDRangeKernel ( priv->cmd_queue, priv->k_mult,
                                    3, NULL, global_work_size, NULL,
                                    0, NULL, &event));
    
    UFO_RESOURCES_CHECK_CLERR (clReleaseEvent (event));
    return TRUE;
}

static void
ufo_complex_mult_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    //UfoComplexMultTaskPrivate *priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_complex_mult_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    //UfoComplexMultTaskPrivate *priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_complex_mult_task_finalize (GObject *object)
{
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

    oclass->set_property = ufo_complex_mult_task_set_property;
    oclass->get_property = ufo_complex_mult_task_get_property;
    oclass->finalize = ufo_complex_mult_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoComplexMultTaskPrivate));
}

static void
ufo_complex_mult_task_init(UfoComplexMultTask *self)
{
    self->priv = UFO_COMPLEX_MULT_TASK_GET_PRIVATE(self);
}
