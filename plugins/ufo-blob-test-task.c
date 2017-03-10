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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ufo-blob-test-task.h"

#define FEPSILON 0.000001

struct _UfoBlobTestTaskPrivate {
    gfloat alpha;
};

static void ufo_task_interface_init (UfoTaskIface *iface);
static gfloat snr(gfloat *img_ptr, gint idx, gint idy, gint s, gint t, gint dimx, gint dimy);

G_DEFINE_TYPE_WITH_CODE (UfoBlobTestTask, ufo_blob_test_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_BLOB_TEST_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_BLOB_TEST_TASK, UfoBlobTestTaskPrivate))

enum {
    PROP_0,
    PROP_THRESHOLD_ALPHA,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_blob_test_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_BLOB_TEST_TASK, NULL));
}

static void
ufo_blob_test_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_blob_test_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition (inputs[0], requisition);
}

static guint
ufo_blob_test_task_get_num_inputs (UfoTask *task)
{
    return 2;
}

static guint
ufo_blob_test_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_blob_test_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR;
}

static gboolean
ufo_blob_test_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoBlobTestTaskPrivate *priv = UFO_BLOB_TEST_TASK_GET_PRIVATE (task);

    float threshold_alpha = priv->alpha;
    unsigned dimx = requisition->dims[0];
    unsigned dimy = requisition->dims[1];
    unsigned size_p = dimx * dimy;
    size_t size = size_p * sizeof (gfloat);

    gfloat * img_mem = ufo_buffer_get_host_array (inputs[0], NULL);
    gfloat * test_mem = ufo_buffer_get_host_array (inputs[1], NULL);
    gfloat * out_mem = ufo_buffer_get_host_array (output, NULL);

    memset (out_mem, 0, size);

    int ct = 0;
    int ct2 = 0;
    for (unsigned ip = 0; ip < size_p; ip++) {
        if (fabs(test_mem[ip] - 1.0f) < FEPSILON) {
            ct++;
            unsigned iy = ip / dimx;
            unsigned ix = ip % dimy; 
            // printf("BlobTestTask: pixel candidate %4u, %4u\n", ix, iy);

            float alpha = snr(img_mem, ix, iy, 2, 3, dimx, dimy);

            if (alpha > threshold_alpha) {
                ct2++;
                out_mem[ip] = alpha;
                // printf("BlobTestTask: detected: %4u %4u %f !!!!!!!\n", ix, iy, alpha);
            } 
        }
    }
    // g_warning ("BlobTestTask: %d out of %d pixels detected as ring cneter", ct2, ct);
    return TRUE;
}

static void
ufo_blob_test_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoBlobTestTaskPrivate *priv = UFO_BLOB_TEST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_THRESHOLD_ALPHA:
            priv->alpha = g_value_get_float (value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_blob_test_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoBlobTestTaskPrivate *priv = UFO_BLOB_TEST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_THRESHOLD_ALPHA:
            g_value_set_float (value, priv->alpha);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_blob_test_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_blob_test_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_blob_test_task_setup;
    iface->get_num_inputs = ufo_blob_test_task_get_num_inputs;
    iface->get_num_dimensions = ufo_blob_test_task_get_num_dimensions;
    iface->get_mode = ufo_blob_test_task_get_mode;
    iface->get_requisition = ufo_blob_test_task_get_requisition;
    iface->process = ufo_blob_test_task_process;
}

static void
ufo_blob_test_task_class_init (UfoBlobTestTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_blob_test_task_set_property;
    oclass->get_property = ufo_blob_test_task_get_property;
    oclass->finalize = ufo_blob_test_task_finalize;

    properties[PROP_THRESHOLD_ALPHA] = 
        g_param_spec_float ("alpha",
                            "alpha threshold",
                            "alpha threshold of signal to noise ratio",
                            0.0, 10.0, 2.0,
                            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoBlobTestTaskPrivate));
}

static void
ufo_blob_test_task_init(UfoBlobTestTask *self)
{
    self->priv = UFO_BLOB_TEST_TASK_GET_PRIVATE(self);
    self->priv->alpha = 2.0f;
}


//////////////////////////////////////////////////////////////////////////////
// internal helper functions
//////////////////////////////////////////////////////////////////////////////
static gfloat array_mean(gfloat *arr, gint n) {
    gfloat res = 0.0f;
    for (gint i = 0; i < n; i++) res += arr[i];
    return res/n;
}

static gfloat array_sigma(gfloat *arr, gfloat mean, gint n) {
    gfloat res = 0.0f;
    for (gint i = 0; i < n; i++) res += ( arr[i] - mean ) * ( arr[i] - mean );
    return sqrt(res/(n-1));
}

static gfloat
snr(gfloat *img_ptr, gint idx, gint idy, gint s, gint t, gint dimx, gint dimy) {
    gint pos;
    gfloat sig[100];
    gfloat bak[400];
    gint mask[21][21];

    // inner size must be smaller than outer size
    // outer size must not be larger than 10
    if (s > t) return -1.0f;
    if (t > 5) return -1.0f;

    // when pixel is close to boundary, discard it
    if (idx - t - 2 < 0) return -1.0f;
    if (idy - t - 2 < 0) return -1.0f;
    if (idx + t + 2>= dimx) return -1.0f;
    if (idy + t + 2>= dimy) return -1.0f;

    // get signal pixels from inner region 
    gint ns = 0;
    for (gint ix = idx - s; ix < idx + s + 1; ix++)
    for (gint iy = idy - s; iy < idy + s + 1; iy++) {
        pos = ix + dimx * iy;
        sig[ns] = img_ptr[pos];
        ns++;
    }

    // get background pixels from outer region
    for (gint i = 0; i < 21; i++)
    for (gint j = 0; j < 21; j++) {
        mask[i][j] = 0;
    }

    gint mid = 10;
    for (gint i = mid-t-2; i < mid+t+3; i++)
    for (gint j = mid-t-2; j < mid+t+3; j++) {
        mask[i][j] = 1;
    }

    for (gint i = mid-t; i < mid+t+1; i++)
    for (gint j = mid-t; j < mid+t+1; j++) {
        mask[i][j] = 0;
    }

    for (gint i = 0; i < 21; i++)                                                                         
    for (gint j = 0; j < 21; j++) {                                                                        
        pos = (i + idx - mid) + (j + idy - mid) * dimx;
        mask[i][j] *= pos;
    }

    gint nb = 0;
    for (gint i = 0; i < 21; i++)                                                                         
    for (gint j = 0; j < 21; j++) {                                                                        
        if ( mask[i][j] != 0 ) {
            pos = mask[i][j];
            bak[nb] = img_ptr[pos];
            nb++;
        }
    }

    gfloat mu_s = array_mean(sig, ns);
    gfloat mu_b = array_mean(bak, nb);

    // gfloat sig_s = array_sigma(sig, mu_s, ns);
    gfloat sig_b = array_sigma(bak, mu_b, nb);

    // gfloat res = (mu_s - mu_b) / sqrt(nb*sig_s*sig_s + ns*sig_b*sig_b);
    gfloat res = (mu_s - mu_b) / sqrt(ns*sig_b*sig_b + ns*sig_b*sig_b/nb);

    // printf("\t%.2f\t%.4e, %.4e, %.4e, %.4e, %4d, %4d\n", res, mu_s-mu_b, mu_b, sig_s, sig_b, ns, nb);
    return res;
}
