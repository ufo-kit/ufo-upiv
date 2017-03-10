/*
 * Copyright (C) 2011-2013 Karlsruhe Institute of Technology
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

#ifdef HAVE_AMD
#include <clFFT.h>
#else
#include "oclFFT.h"
#endif

#include <math.h>
#include <assert.h>

#include "ufo-hessian-analysis-task.h"

/**
 * SECTION:ufo-hessian-analysis-task
 * @Short_description: Write TIFF files
 * @Title: hessian_analysis
 *
 */

struct _UfoHessianAnalysisTaskPrivate {
    gint width;
    gint height;
    gint n_samples;

    cl_context context;
    cl_command_queue cmd_queue;
    cl_kernel kernel_det;
    cl_kernel kernel_eigval;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoHessianAnalysisTask, ufo_hessian_analysis_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_HESSIAN_ANALYSIS_TASK, UfoHessianAnalysisTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_hessian_analysis_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_HESSIAN_ANALYSIS_TASK, NULL));
}

static void
ufo_hessian_analysis_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoHessianAnalysisTaskPrivate *priv;
    UfoGpuNode *node;

    priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));

    priv->context = ufo_resources_get_context (resources);
    priv->cmd_queue = ufo_gpu_node_get_cmd_queue (node);

    UFO_RESOURCES_CHECK_CLERR (clRetainContext (priv->context));

    priv->kernel_det  = ufo_resources_get_kernel (resources, "hessian.cl", "hessian_det", error);
    if (priv->kernel_det != NULL) {
        UFO_RESOURCES_CHECK_CLERR (clRetainKernel (priv->kernel_det));
    }

    priv->kernel_eigval  = ufo_resources_get_kernel (resources, "hessian.cl", "hessian_eigval", error);
    if (priv->kernel_eigval != NULL) {
        UFO_RESOURCES_CHECK_CLERR (clRetainKernel (priv->kernel_eigval));
    }
}

static void
ufo_hessian_analysis_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    UfoRequisition req0;
    ufo_buffer_get_requisition (inputs[0], &req0);

    UfoHessianAnalysisTaskPrivate *priv;
    priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (task);
    priv->n_samples = (gint) req0.dims[2] / 3;

    //printf("hessian analysis \t%d %d %d %d\n", (int) req0.n_dims, (int) req0.dims[0], (int) req0.dims[1], (int) req0.dims[2]);

    requisition->n_dims = 3;
    requisition->dims[0] = req0.dims[0];
    requisition->dims[1] = req0.dims[1];
    requisition->dims[2] = (gsize) priv->n_samples;

    priv->width = req0.dims[0];
    priv->height = req0.dims[1];
}

static guint
ufo_hessian_analysis_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_hessian_analysis_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 3;
}

static UfoTaskMode
ufo_hessian_analysis_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static gboolean
ufo_hessian_analysis_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    cl_event event;
    cl_mem in_mem, out_mem;
    gsize global_work_size[2];

    UfoHessianAnalysisTaskPrivate *priv;
    priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (task);

    in_mem = ufo_buffer_get_device_array(inputs[0], priv->cmd_queue);
    out_mem = ufo_buffer_get_device_array(output, priv->cmd_queue);

    global_work_size[0] = (gsize) priv->width;
    global_work_size[1] = (gsize) priv->height;
    //printf("%d %d\n", global_work_size[0], global_work_size[1]);

//#define EIGEN
#ifndef EIGEN
    UFO_RESOURCES_CHECK_CLERR (
               clSetKernelArg ( priv->kernel_det, 0, 
                                sizeof (cl_mem), (gpointer) &out_mem));
    UFO_RESOURCES_CHECK_CLERR (
               clSetKernelArg ( priv->kernel_det, 1, 
                                sizeof (cl_mem), (gpointer) &in_mem));
    UFO_RESOURCES_CHECK_CLERR (
       clEnqueueNDRangeKernel ( priv->cmd_queue, priv->kernel_det,
                                2, NULL, global_work_size, NULL,
                                0, NULL, &event));
#else
    UFO_RESOURCES_CHECK_CLERR (
               clSetKernelArg ( priv->kernel_eigval, 0, 
                                sizeof (cl_mem), (gpointer) &out_mem));
    UFO_RESOURCES_CHECK_CLERR (
               clSetKernelArg ( priv->kernel_eigval, 1, 
                                sizeof (cl_mem), (gpointer) &in_mem));
    UFO_RESOURCES_CHECK_CLERR (
       clEnqueueNDRangeKernel ( priv->cmd_queue, priv->kernel_eigval,
                                2, NULL, global_work_size, NULL,
                                0, NULL, &event));
#endif


    UFO_RESOURCES_CHECK_CLERR (clReleaseEvent (event));

    return TRUE;
}

static void
ufo_hessian_analysis_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    // UfoHessianAnalysisTaskPrivate *priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hessian_analysis_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    // UfoHessianAnalysisTaskPrivate *priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hessian_analysis_task_finalize (GObject *object)
{

    UfoHessianAnalysisTaskPrivate *priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE (object);

    if (priv->kernel_det) {
        UFO_RESOURCES_CHECK_CLERR (clReleaseKernel (priv->kernel_det));
        priv->kernel_det = NULL;
    }

    if (priv->kernel_eigval) {
        UFO_RESOURCES_CHECK_CLERR (clReleaseKernel (priv->kernel_eigval));
        priv->kernel_eigval = NULL;
    }

    G_OBJECT_CLASS (ufo_hessian_analysis_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_hessian_analysis_task_setup;
    iface->get_num_inputs = ufo_hessian_analysis_task_get_num_inputs;
    iface->get_num_dimensions = ufo_hessian_analysis_task_get_num_dimensions;
    iface->get_mode = ufo_hessian_analysis_task_get_mode;
    iface->get_requisition = ufo_hessian_analysis_task_get_requisition;
    iface->process = ufo_hessian_analysis_task_process;
}

static void
ufo_hessian_analysis_task_class_init (UfoHessianAnalysisTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_hessian_analysis_task_set_property;
    oclass->get_property = ufo_hessian_analysis_task_get_property;
    oclass->finalize = ufo_hessian_analysis_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoHessianAnalysisTaskPrivate));
}

static void
ufo_hessian_analysis_task_init(UfoHessianAnalysisTask *self)
{
    self->priv = UFO_HESSIAN_ANALYSIS_TASK_GET_PRIVATE(self);
}
