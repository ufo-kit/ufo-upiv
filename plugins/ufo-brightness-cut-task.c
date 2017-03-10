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

#include <string.h>
#include <math.h>
#include "ufo-brightness-cut-task.h"


struct _UfoBrightnessCutTaskPrivate {
    float low;
    float high;
};

static void ufo_task_interface_init (UfoTaskIface *iface);
static gfloat array_mean(gfloat *arr, unsigned n);
static gfloat array_std(gfloat *arr, gfloat mean, unsigned n);

G_DEFINE_TYPE_WITH_CODE (UfoBrightnessCutTask, ufo_brightness_cut_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_BRIGHTNESS_CUT_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_BRIGHTNESS_CUT_TASK, UfoBrightnessCutTaskPrivate))

enum {
    PROP_0,
    PROP_LOW,
    PROP_HIGH,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_brightness_cut_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_BRIGHTNESS_CUT_TASK, NULL));
}

static void
ufo_brightness_cut_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_brightness_cut_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition (inputs[0], requisition);
}

static guint
ufo_brightness_cut_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_brightness_cut_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_brightness_cut_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_CPU;
}

static gboolean
ufo_brightness_cut_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoBrightnessCutTaskPrivate *priv = UFO_BRIGHTNESS_CUT_TASK_GET_PRIVATE (task);

    gfloat *in_mem = ufo_buffer_get_host_array (inputs[0], NULL);
    gfloat *out_mem = ufo_buffer_get_host_array (output, NULL);
    unsigned size = requisition->dims[0] * requisition->dims[1];

    gfloat mean = array_mean(in_mem, size);
    gfloat std = array_std(in_mem, mean, size);

    gfloat low = mean + priv->low * std;
    gfloat high = mean + priv->high * std;

    memset (out_mem, 0, size * sizeof(gfloat));
    for (unsigned i = 0; i < size; i++) {
        if (in_mem[i] > low && in_mem[i] < high) {
            out_mem[i] = in_mem[i];
        }
    }

    return TRUE;
}

static void
ufo_brightness_cut_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoBrightnessCutTaskPrivate *priv = UFO_BRIGHTNESS_CUT_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_LOW:
            priv->low = g_value_get_float (value);
            break;
        case PROP_HIGH:
            priv->high = g_value_get_float (value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_brightness_cut_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoBrightnessCutTaskPrivate *priv = UFO_BRIGHTNESS_CUT_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_LOW:
            g_value_set_float (value, priv->low);
            break;
        case PROP_HIGH:
            g_value_set_float (value, priv->high);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_brightness_cut_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_brightness_cut_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_brightness_cut_task_setup;
    iface->get_num_inputs = ufo_brightness_cut_task_get_num_inputs;
    iface->get_num_dimensions = ufo_brightness_cut_task_get_num_dimensions;
    iface->get_mode = ufo_brightness_cut_task_get_mode;
    iface->get_requisition = ufo_brightness_cut_task_get_requisition;
    iface->process = ufo_brightness_cut_task_process;
}

static void
ufo_brightness_cut_task_class_init (UfoBrightnessCutTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_brightness_cut_task_set_property;
    oclass->get_property = ufo_brightness_cut_task_get_property;
    oclass->finalize = ufo_brightness_cut_task_finalize;

    properties[PROP_LOW] =
        g_param_spec_float ("low",
            "low threshold",
            "low threshold",
            0.0, 10.0, 1.0,
            G_PARAM_READWRITE);

    properties[PROP_HIGH] =
        g_param_spec_float ("high",
            "high threshold",
            "high threshold",
            0.0, 10.0, 5.0,
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoBrightnessCutTaskPrivate));
}

static void
ufo_brightness_cut_task_init(UfoBrightnessCutTask *self)
{
    self->priv = UFO_BRIGHTNESS_CUT_TASK_GET_PRIVATE(self);
    self->priv->low = 1.0f;
    self->priv->high = 3.0f;
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
