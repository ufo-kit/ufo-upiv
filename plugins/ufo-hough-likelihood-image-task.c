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

#include "ufo-hough-likelihood-image-task.h"

struct _UfoHoughLikelihoodImageTaskPrivate {
    gint mask_num_ones;
    gint masksize;
    gint masksize_h;
    gfloat maskinnersize;
    cl_context context;
    cl_kernel kernel;
    cl_mem mask_mem;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoHoughLikelihoodImageTask, ufo_hough_likelihood_image_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_HOUGH_LIKELIHOOD_IMAGE_TASK, UfoHoughLikelihoodImageTaskPrivate))

enum {
    PROP_0,
    PROP_MASKSIZE,
    PROP_MASKINNERSIZE,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_hough_likelihood_image_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_HOUGH_LIKELIHOOD_IMAGE_TASK, NULL));
}

static void
ufo_hough_likelihood_image_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoHoughLikelihoodImageTaskPrivate *priv;
    cl_int err;
    int i, j, i0, j0;
    int *mask;
    float maskinnersize_h;

    priv = UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE (task);
    priv->context = ufo_resources_get_context (resources);
    priv->kernel = ufo_resources_get_kernel (resources, "hough.cl", "likelihood_image", error);

    if (priv->kernel)
        UFO_RESOURCES_CHECK_CLERR(clRetainKernel(priv->kernel));
    if (priv->context)
        UFO_RESOURCES_CHECK_CLERR(clRetainContext(priv->context));

    mask = g_malloc0 (priv->masksize * priv->masksize * sizeof (int));
    priv->masksize_h = (priv->masksize - 1) / 2;
    maskinnersize_h = priv->maskinnersize / 2;
    priv->mask_num_ones = 0;
    for (i = 0, i0 = - priv->masksize_h; i < priv->masksize; i++, i0++)
        for (j = 0, j0 = - priv->masksize_h; j < priv->masksize; j++, j0++)
            if ((i0*i0 + j0*j0) >= maskinnersize_h * maskinnersize_h)
            {
                mask[i + j*priv->masksize] = 1;
                priv->mask_num_ones += 1;
            }

    // debug mask
    /*
     *int *mask0 = mask;
     *for (i = 0; i < priv->masksize; i++)
     *{
     *    for (j = 0; j < priv->masksize; j++, mask0++)
     *        g_printf ("%d, ", *mask0);
     *    g_printf("\n");
     *}
     */

    if (priv->mask_mem)
        UFO_RESOURCES_CHECK_CLERR (clReleaseMemObject (priv->mask_mem));
    
    priv->mask_mem = clCreateBuffer (priv->context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            priv->masksize * priv->masksize * sizeof (*mask), mask, &err);
    UFO_RESOURCES_CHECK_CLERR (err);
}

static void
ufo_hough_likelihood_image_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition (inputs[0], requisition);
}

static guint
ufo_hough_likelihood_image_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_hough_likelihood_image_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return -1;
}

static UfoTaskMode
ufo_hough_likelihood_image_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static gboolean
ufo_hough_likelihood_image_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoHoughLikelihoodImageTaskPrivate *priv;
    UfoGpuNode *node;
    UfoProfiler *profiler;
    cl_command_queue cmd_queue;

    priv = UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));
    profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));
    cmd_queue = ufo_gpu_node_get_cmd_queue (node);

    gsize g_work_size[] = { requisition->dims[0], requisition->dims[1], 1};
    if (requisition->n_dims == 3)
        g_work_size[2] = requisition->dims[2];
    gsize l_work_size[] = { 32, 16, 1};

    gsize shift = 6;
    gsize l_mem_size = sizeof(float) * (2*shift + l_work_size[0]) * (2*shift + l_work_size[1]);

    cl_mem in_img = ufo_buffer_get_device_image (inputs[0], cmd_queue);
    cl_mem out_mem = ufo_buffer_get_device_array (output, cmd_queue);

    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 0, sizeof(cl_mem), &in_img));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 1, sizeof(cl_mem), &out_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 2, sizeof(cl_mem), &priv->mask_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 3, sizeof(int), &priv->masksize_h));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 4, sizeof(int), &priv->mask_num_ones));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->kernel, 5, l_mem_size, NULL));

    ufo_profiler_call (profiler, cmd_queue, priv->kernel, 3, g_work_size, l_work_size);

    return TRUE;
}

static void
ufo_hough_likelihood_image_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoHoughLikelihoodImageTaskPrivate *priv = UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_MASKSIZE:
            priv->masksize = g_value_get_int (value);
            if (priv->masksize % 2 != 1)
            {
                g_warning ("masksize (%d) must be odd number", priv->masksize);
                priv->masksize -= 1;
            }

            break;
        case PROP_MASKINNERSIZE:
            priv->maskinnersize = g_value_get_float (value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hough_likelihood_image_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoHoughLikelihoodImageTaskPrivate *priv = UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_MASKSIZE:
            g_value_set_int (value, priv->masksize);
            break;
        case PROP_MASKINNERSIZE:
            g_value_set_float (value, priv->maskinnersize);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hough_likelihood_image_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_hough_likelihood_image_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_hough_likelihood_image_task_setup;
    iface->get_num_inputs = ufo_hough_likelihood_image_task_get_num_inputs;
    iface->get_num_dimensions = ufo_hough_likelihood_image_task_get_num_dimensions;
    iface->get_mode = ufo_hough_likelihood_image_task_get_mode;
    iface->get_requisition = ufo_hough_likelihood_image_task_get_requisition;
    iface->process = ufo_hough_likelihood_image_task_process;
}

static void
ufo_hough_likelihood_image_task_class_init (UfoHoughLikelihoodImageTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_hough_likelihood_image_task_set_property;
    oclass->get_property = ufo_hough_likelihood_image_task_get_property;
    oclass->finalize = ufo_hough_likelihood_image_task_finalize;

    properties[PROP_MASKSIZE] =
        g_param_spec_int ("masksize",
            "Mask size",
            "Mask size",
            3, 21, 9,
            G_PARAM_READWRITE);

    properties[PROP_MASKINNERSIZE] =
        g_param_spec_float ("maskinnersize",
            "Mask inner size",
            "Mask inner size",
            1.f, 10.f, 1.f,
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoHoughLikelihoodImageTaskPrivate));
}

static void
ufo_hough_likelihood_image_task_init(UfoHoughLikelihoodImageTask *self)
{
    self->priv = UFO_HOUGH_LIKELIHOOD_IMAGE_TASK_GET_PRIVATE(self);
    self->priv->masksize = 9;
    self->priv->maskinnersize = 2.0;
}
