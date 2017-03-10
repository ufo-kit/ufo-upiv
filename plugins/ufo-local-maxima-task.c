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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ufo-local-maxima-task.h"

struct _UfoLocalMaximaTaskPrivate {
    cl_kernel locmax_kernel;
    gfloat sigma;
};

static void ufo_task_interface_init (UfoTaskIface *iface);
static gfloat array_mean(gfloat *arr, unsigned n);
static gfloat array_std(gfloat *arr, gfloat mean, unsigned n);

G_DEFINE_TYPE_WITH_CODE (UfoLocalMaximaTask, ufo_local_maxima_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_LOCAL_MAXIMA_TASK, UfoLocalMaximaTaskPrivate))

enum {
    PROP_0,
    PROP_SIGMA,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_local_maxima_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_LOCAL_MAXIMA_TASK, NULL));
}

static void
ufo_local_maxima_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoLocalMaximaTaskPrivate* priv;

    priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE(task);
    priv->locmax_kernel = ufo_resources_get_kernel(resources, "localmax.cl", NULL, error);

    if (priv->locmax_kernel != NULL)
    {
        UFO_RESOURCES_CHECK_CLERR(clRetainKernel (priv->locmax_kernel));
    }
}

static void
ufo_local_maxima_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition(inputs[0], requisition);
}

static guint
ufo_local_maxima_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_local_maxima_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_local_maxima_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static gboolean
ufo_local_maxima_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoLocalMaximaTaskPrivate *priv;
    UfoGpuNode *node;
    UfoProfiler *profiler;
    cl_command_queue cmd_queue;

    priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));
    cmd_queue = ufo_gpu_node_get_cmd_queue (node);
    profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));

    gfloat sigma = priv->sigma;
    gsize img_size = ufo_buffer_get_size(inputs[0]);
    gsize img_size_p = img_size / sizeof (gfloat);
    // g_warning ("LocalMaxima: input image size = %u", img_size_p);
    // img_size: image size in bytes
    // img_size_p: image size in pixels

    gfloat * in_mem = ufo_buffer_get_host_array (inputs[0], NULL);

    gfloat mean = array_mean(in_mem, img_size_p);
    gfloat std = array_std(in_mem, mean, img_size_p);
    gfloat threshold = mean + sigma * std;
    // g_warning ("LocalMaxima: mean, std = %6.2e %6.2e", mean, std);

    cl_mem in_mem_gpu = ufo_buffer_get_device_array(inputs[0], cmd_queue);
    cl_mem out_mem_gpu = ufo_buffer_get_device_array(output,cmd_queue);

    UFO_RESOURCES_CHECK_CLERR(clSetKernelArg(priv->locmax_kernel,0,sizeof(cl_mem), &in_mem_gpu));
    UFO_RESOURCES_CHECK_CLERR(clSetKernelArg(priv->locmax_kernel,1,sizeof(cl_mem), &out_mem_gpu));
    UFO_RESOURCES_CHECK_CLERR(clSetKernelArg(priv->locmax_kernel,2,sizeof(cl_float), &threshold));

    ufo_profiler_call(profiler,cmd_queue, priv->locmax_kernel,1,&img_size_p,NULL);

    return TRUE;
}

static void
ufo_local_maxima_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoLocalMaximaTaskPrivate *priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_SIGMA:
            priv->sigma = g_value_get_float (value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_local_maxima_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoLocalMaximaTaskPrivate *priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_SIGMA:
            g_value_set_float (value, priv->sigma);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_local_maxima_task_finalize (GObject *object)
{

    UfoLocalMaximaTaskPrivate* priv;
    priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE(object);

    if(priv->locmax_kernel)
    {
        UFO_RESOURCES_CHECK_CLERR(clReleaseKernel(priv->locmax_kernel));
        priv->locmax_kernel = NULL;
    }
       
    G_OBJECT_CLASS (ufo_local_maxima_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_local_maxima_task_setup;
    iface->get_num_inputs = ufo_local_maxima_task_get_num_inputs;
    iface->get_num_dimensions = ufo_local_maxima_task_get_num_dimensions;
    iface->get_mode = ufo_local_maxima_task_get_mode;
    iface->get_requisition = ufo_local_maxima_task_get_requisition;
    iface->process = ufo_local_maxima_task_process;
}

static void
ufo_local_maxima_task_class_init (UfoLocalMaximaTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_local_maxima_task_set_property;
    oclass->get_property = ufo_local_maxima_task_get_property;
    oclass->finalize = ufo_local_maxima_task_finalize;

    properties[PROP_SIGMA] =
        g_param_spec_float ("sigma",
            "threshold sigma",
            "set threshold at <sigma> above mean",
            0.0f, 10.0f, 5.0f,
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoLocalMaximaTaskPrivate));
}

static void
ufo_local_maxima_task_init(UfoLocalMaximaTask *self)
{
    self->priv = UFO_LOCAL_MAXIMA_TASK_GET_PRIVATE(self);
    self->priv->sigma = 5.0f;
}

//////////////////////////////////////////////////////////////////////////////
// internal helper functions
//////////////////////////////////////////////////////////////////////////////
static gfloat array_mean(gfloat *arr, unsigned n) {
    gfloat res = 0.0f;
    for (unsigned i = 0; i < n; i++) res += arr[i];
    return res/n;
}

static gfloat array_std(gfloat *arr, gfloat mean, unsigned n) {
    gfloat res = 0.0f;
    for (unsigned i = 0; i < n; i++) res += ( arr[i] - mean ) * ( arr[i] - mean );
    return sqrt(res/(n-1));
}
