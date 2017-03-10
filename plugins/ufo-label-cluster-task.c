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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ufo-label-cluster-task.h"

#define EPSILON 1.0e-8

struct _UfoLabelClusterTaskPrivate {
    gboolean foo;
};

static void ufo_task_interface_init (UfoTaskIface *iface);
static void                                                                     
explore (gfloat *src, gfloat *dst, gfloat threshold, guint label, 
         gint x, gint y, gint dimx, gint dimy,
         gint *depth, gint *minx, gint *maxx, gint *miny, gint *maxy,
         gfloat* minv);

G_DEFINE_TYPE_WITH_CODE (UfoLabelClusterTask, ufo_label_cluster_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_LABEL_CLUSTER_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_LABEL_CLUSTER_TASK, UfoLabelClusterTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_label_cluster_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_LABEL_CLUSTER_TASK, NULL));
}

static void
ufo_label_cluster_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_label_cluster_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    UfoRequisition req_in;
    ufo_buffer_get_requisition (inputs[0], &req_in);
    *requisition = req_in;
}

static guint
ufo_label_cluster_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_label_cluster_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_label_cluster_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_CPU;
}

static gboolean
ufo_label_cluster_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    gfloat threshold = 1.0f;
    int dimx = requisition->dims[0];
    int dimy = requisition->dims[1];
    
    gfloat *in_mem = ufo_buffer_get_host_array (inputs[0], NULL);
    gfloat *out_mem = ufo_buffer_get_host_array (output, NULL);

    gsize out_size = ufo_buffer_get_size (output);
    memset (out_mem, 0, out_size);

    unsigned label = 0;
    for (int ix = 0; ix < (int) dimx; ix++)
    for (int iy = 0; iy < (int) dimx; iy++) {

        unsigned pos = ix + dimx * iy;
        if (fabs(in_mem[pos] - threshold) < EPSILON 
                && fabs(out_mem[pos] - 0.0f) < EPSILON
                && fabs(out_mem[pos] + 1.0f) > EPSILON) {

            /*
             *printf("%2f %+2f (%u, %u, %u)\n", in_mem[pos], out_mem[pos],
             *        fabs(in_mem[pos] - threshold) < EPSILON,
             *        fabs(out_mem[pos]-0.0f) < EPSILON,
             *        fabs(out_mem[pos]+1.0f) > EPSILON );
             */

            int depth = 0;
            int minx = dimx;
            int miny = dimy;
            int maxx = 0;
            int maxy = 0;
            gfloat minv = 0.0f;

            label++;
            explore(in_mem, out_mem, threshold, label, 
                    ix, iy, dimx, dimy, &depth, &minx, &maxx, &miny, &maxy, &minv);
            
            int csize = (maxx-minx+1) * (maxy-miny+1);
            if (csize < 4 || csize > 50 || abs(maxx-minx-maxy+miny) > 2) {
                for (int cx = minx; cx < maxx + 1; cx++)
                for (int cy = miny; cy < maxy + 1; cy++) {
                    unsigned cpos = cx + dimx * cy;
                    out_mem[cpos] = -1.0f;
                }
            } else {
                /*
                 *printf("LabelCluster: %4d (%3d), [%3d:%3d, %3d:%3d], %4d, %6f\n", 
                 *        label, depth, minx, maxx, miny, maxy,
                 *        csize, (float) depth / csize );
                 */
            }
        }
    }

    for (int ix = 0; ix < dimx * dimy; ix++) {
        if (out_mem[ix] > 0) {
            out_mem[ix] = log (out_mem[ix]);
        }
    }
    
    return TRUE;
}

static void
ufo_label_cluster_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoLabelClusterTaskPrivate *priv = UFO_LABEL_CLUSTER_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_label_cluster_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoLabelClusterTaskPrivate *priv = UFO_LABEL_CLUSTER_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_label_cluster_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_label_cluster_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_label_cluster_task_setup;
    iface->get_num_inputs = ufo_label_cluster_task_get_num_inputs;
    iface->get_num_dimensions = ufo_label_cluster_task_get_num_dimensions;
    iface->get_mode = ufo_label_cluster_task_get_mode;
    iface->get_requisition = ufo_label_cluster_task_get_requisition;
    iface->process = ufo_label_cluster_task_process;
}

static void
ufo_label_cluster_task_class_init (UfoLabelClusterTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_label_cluster_task_set_property;
    oclass->get_property = ufo_label_cluster_task_get_property;
    oclass->finalize = ufo_label_cluster_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoLabelClusterTaskPrivate));
}

static void
ufo_label_cluster_task_init(UfoLabelClusterTask *self)
{
    self->priv = UFO_LABEL_CLUSTER_TASK_GET_PRIVATE(self);
}

static void                                                                     
explore (gfloat *src, gfloat *dst, gfloat threshold, guint label, 
         gint x, gint y, gint dimx, gint dimy,
         gint *depth, gint *minx, gint *maxx, gint *miny, gint *maxy,
         gfloat* minv)      
{                                                                               
    if (*depth > 100) return;
    if (x < 0 || y < 0 || x >  dimx - 1 || y >  dimy - 1) return;

    guint pos = x + dimx * y;
    if (src[pos] < threshold) return;
    if (fabs(dst[pos] - label) < EPSILON) return;

    dst[pos] = label;
    (*depth)++;

    /*
     *printf("\tpos = %d (%4d, %4d)", pos, x, y);
     *printf("\tsrc = %5e %5d %5d\n", src[pos], (int) dst[pos], *depth);
     */

    *minx = MIN (x, *minx);
    *miny = MIN (y, *miny);
    *maxx = MAX (x, *maxx);
    *maxy = MAX (y, *maxy);
    *minv = MIN (src[pos], *minv);

    explore(src, dst, threshold, label, x+1, y, dimx, dimy, depth, minx, maxx, miny, maxy, minv);
    explore(src, dst, threshold, label, x-1, y, dimx, dimy, depth, minx, maxx, miny, maxy, minv);
    explore(src, dst, threshold, label, x, y-1, dimx, dimy, depth, minx, maxx, miny, maxy, minv);
    explore(src, dst, threshold, label, x, y+1, dimx, dimy, depth, minx, maxx, miny, maxy, minv);
} 
