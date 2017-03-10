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

#include <math.h>
#include <assert.h>

#include "ufo-hessian-kernel-task.h"

/**
 * SECTION:ufo-hessian-kernel-task
 * @Short_description: Write TIFF files
 * @Title: hessian_kernel
 *
 */

struct _UfoHessianKernelTaskPrivate {
    gsize width;
    gsize height;
    gfloat sigma;
    gboolean flag;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoHessianKernelTask, ufo_hessian_kernel_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_HESSIAN_KERNEL_TASK, UfoHessianKernelTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    PROP_WIDTH,
    PROP_HEIGHT,
    PROP_SIGMA,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_hessian_kernel_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_HESSIAN_KERNEL_TASK, NULL));
}

static guint32
pow2round(guint32 x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

static void
ufo_hessian_kernel_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoHessianKernelTaskPrivate *priv;

    priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (task);

    priv->width = pow2round(priv->width);
    priv->height = pow2round(priv->height);
}

static void
ufo_hessian_kernel_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    UfoHessianKernelTaskPrivate *priv;

    priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (task);

    requisition->n_dims = 3;
    requisition->dims[0] = priv->width;
    requisition->dims[1] = priv->height;
    requisition->dims[2] = 3;
}

static guint
ufo_hessian_kernel_task_get_num_inputs (UfoTask *task)
{
    UfoHessianKernelTaskPrivate *priv;
    priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (task);
    priv->flag = FALSE;
    return 0;
}

static guint
ufo_hessian_kernel_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 0;
}

static UfoTaskMode
ufo_hessian_kernel_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_GENERATOR | UFO_TASK_MODE_CPU;
}

static gboolean
ufo_hessian_kernel_task_generate (UfoTask *task,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    int isig, i, j, s0, s1, s2, x, y;
    int width, height, img_size, kernel_ext;
    float sigma, sig2, sig2inv, fct1, fct2, gauss, hxx, hyy, hxy;
    float PIinv = 0.31830988618f;
    float *res;

    UfoHessianKernelTaskPrivate *priv;
    priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (task);

    if (priv->flag) {
        return FALSE;
    }

    width = priv->width;
    height = priv->height;
    img_size = width * height;

    UfoRequisition req;
    ufo_buffer_get_requisition(output, &req);

    res = ufo_buffer_get_host_array(output, NULL);

    isig = 0;
    sigma = priv->sigma;
    kernel_ext = 3 * sigma;

    sig2 = sigma * sigma;
    sig2inv = 1.0f / sig2;
    fct1 = 0.5f * sig2inv * sig2inv * PIinv;
    fct2 = fct1 * sig2inv;

    // store Hessian kernel in wrap order 
    for (i = -kernel_ext; i < kernel_ext + 1; i++) { // height direction
        y = (i >= 0) ? i : height + i;
        for (j = -kernel_ext; j < kernel_ext + 1; j++) { // width direction
            x = (j >= 0) ? j : width + j;
            s0 = 3 * isig * img_size + y * width + x;
            s1 = (3 * isig + 1) * img_size + y * width + x;
            s2 = (3 * isig + 2) * img_size + y * width + x;
            
            gauss = exp(-0.5f*sig2inv*(i*i + j*j));
            hxx = fct1 * ( i * i * sig2inv - 1 ) * gauss;
            hyy = fct1 * ( j * j * sig2inv - 1 ) * gauss;
            hxy = fct2 * ( i * j ) * gauss;

            /*float tmp1, tmp2, lmd1, lmd2; */
            /*float s = isig * img_size + y * width + x;*/

            /*tmp1 = 0.5 * (hxx  + hyy);*/
            /*tmp2 = 0.5 * sqrt(4 * hxy * hxy + (hxx - hyy)*(hxx - hyy));*/

            /*lmd1 = tmp1 + tmp2;*/
            /*lmd2 = tmp1 - tmp2;*/

            /*
             *if (lmd1 < 0 && lmd2 < 0) {
             *    lmd1 = -1.0f * lmd1;
             *    lmd2 = -1.0f * lmd2;
             *    res[s] = (lmd1 > lmd2) ? lmd2 / lmd1: lmd1 / lmd2;
             *} else {
             *    res[s] = 0.0f;
             *    printf("negative value found\n");
             *}
             */

            res[s0] = hxx;
            res[s1] = hyy;
            res[s2] = hxy;

            /*res[s] = fabs( lmd1 / lmd2 ) / (gfloat) (sigma *sigma);*/

            //res[s] = pow(sigma, 6) *  fabs(hxx*hyy - hxy*hxy);
            //res[s] = pow(sigma, 6) *  hxx;
        }
    }

    priv->flag = TRUE;

    return TRUE;
}

static void
ufo_hessian_kernel_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoHessianKernelTaskPrivate *priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        case PROP_WIDTH:
            priv->width =  g_value_get_uint(value);
            break;
        case PROP_HEIGHT:
            priv->height =  g_value_get_uint(value);
            break;
        case PROP_SIGMA:
            priv->sigma =  g_value_get_float(value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hessian_kernel_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoHessianKernelTaskPrivate *priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        case PROP_WIDTH:
            g_value_set_uint (value, priv->width);
            break;
        case PROP_HEIGHT:
            g_value_set_uint (value, priv->height);
            break;
        case PROP_SIGMA:
            g_value_set_float (value, priv->sigma);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_hessian_kernel_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_hessian_kernel_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_hessian_kernel_task_setup;
    iface->get_num_inputs = ufo_hessian_kernel_task_get_num_inputs;
    iface->get_num_dimensions = ufo_hessian_kernel_task_get_num_dimensions;
    iface->get_mode = ufo_hessian_kernel_task_get_mode;
    iface->get_requisition = ufo_hessian_kernel_task_get_requisition;
    iface->generate = ufo_hessian_kernel_task_generate;
}

static void
ufo_hessian_kernel_task_class_init (UfoHessianKernelTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_hessian_kernel_task_set_property;
    oclass->get_property = ufo_hessian_kernel_task_get_property;
    oclass->finalize = ufo_hessian_kernel_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    properties[PROP_WIDTH] = 
        g_param_spec_uint ("width",
            "Set kernel image width",
            "kernel size width should be equal with image size to be convolved",
            256, 4096, 1024,
            G_PARAM_READWRITE);

    properties[PROP_HEIGHT] = 
        g_param_spec_uint ("height",
            "Set kernel image height",
            "kernel image size should be equal with image size to be convolved",
            256, 4096, 1024,
            G_PARAM_READWRITE);

    properties[PROP_SIGMA] =
        g_param_spec_float ("sigma",
            "Set hessian-kernel width sigma",
            "Set hessian-kernel width sigma",
            0.5, 6, 2,
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoHessianKernelTaskPrivate));
}

static void
ufo_hessian_kernel_task_init(UfoHessianKernelTask *self)
{
    self->priv = UFO_HESSIAN_KERNEL_TASK_GET_PRIVATE(self);
    self->priv->width = 1024;
    self->priv->height = 1024;
    self->priv->sigma = 2;
}
