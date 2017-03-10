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

/*
 *  Filter: Bilateral Filter
 *  Author: Chuan Miao
 *
 *  cf. https://en.wikipedia.org/wiki/Bilateral_fitler
 *
 *  @parm sigr Width of the intencity smoothing kernel (Gaussian function)
 *  @parm sigd Width of the spatial smoothing kernel (Gaussian function)
 *  @parm window Smoothing window size
 */

#include "ufo-bilateral-task.h"
#include <math.h>


struct _UfoBilateralTaskPrivate {
    gboolean foo;
    gfloat sigr; // intensitiy
    gfloat sigd; // spatial
    gint window;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoBilateralTask, ufo_bilateral_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_BILATERAL_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_BILATERAL_TASK, UfoBilateralTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    PROP_SIGR,
    PROP_SIGD,
    PORP_WINDOW,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_bilateral_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_BILATERAL_TASK, NULL));
}

static void
ufo_bilateral_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_bilateral_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition (inputs[0], requisition);
}

static guint
ufo_bilateral_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_bilateral_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_bilateral_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR;
}

static gboolean
ufo_bilateral_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoBilateralTaskPrivate *priv = UFO_BILATERAL_TASK_GET_PRIVATE (task);

    int i, j, k, l, c;
    int lx = requisition->dims[0];
    int ly = requisition->dims[1];
    int window = priv->window;
    float sigr = priv->sigr;
    float sigd = priv->sigd;

    gfloat *in_mem = ufo_buffer_get_host_array (inputs[0], NULL);
    gfloat *out_mem = ufo_buffer_get_host_array (output, NULL);
    
    float i0;
    int ct = 0;
    gfloat *wr = g_malloc0 (window*window*sizeof(gfloat));
    gfloat *wd = g_malloc0 (window*window*sizeof(gfloat));
    gfloat *wp = g_malloc0 (window*window*sizeof(gfloat));
    gfloat wh, wg;

    for (i = 0; i < lx; i++) {
        for (j = 0; j < ly; j++) {
            i0 = in_mem[j*lx + i];
            ct = 0;

            for (k = ((i - window/2) <0 ? 0 : (i - window/2)); 
                    k < ((i + window/2 + 1) > lx ? lx : (i + window/2 + 1)); k++)
                for (l = ((j - window/2) <0 ? 0 : (j - window/2)); 
                        l < ((j + window/2 + 1) > ly ? ly : (j + window/2 + 1)); l++) {
                    wr[ct] = in_mem[l*lx + k];
                    wd[ct] = (k-i)*(k-i) + (l-j)*(l-j);
                    ct++;
                }

            wg = 0.0f;
            wh = 0.0f;
            for (c = 0; c < ct; c++) {
                wp[c] = -0.5 * (wr[c] -i0) * (wr[c] - i0) / (sigr*sigr);
                wp[c] += -0.5 * wd[c] / (sigd*sigd); 
                wp[c] = exp (wp[c]);
                wh += wr[c] * wp[c];
                wg += wp[c];
            }

            out_mem [j*lx+i] = wh / wg;
        }
    }

    g_free (wr);
    g_free (wd);
    g_free (wp);

    return TRUE;
}

static void
ufo_bilateral_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoBilateralTaskPrivate *priv = UFO_BILATERAL_TASK_GET_PRIVATE (object);
    (void) *priv;

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_bilateral_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoBilateralTaskPrivate *priv = UFO_BILATERAL_TASK_GET_PRIVATE (object);
    (void) *priv;

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_bilateral_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_bilateral_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_bilateral_task_setup;
    iface->get_num_inputs = ufo_bilateral_task_get_num_inputs;
    iface->get_num_dimensions = ufo_bilateral_task_get_num_dimensions;
    iface->get_mode = ufo_bilateral_task_get_mode;
    iface->get_requisition = ufo_bilateral_task_get_requisition;
    iface->process = ufo_bilateral_task_process;
}

static void
ufo_bilateral_task_class_init (UfoBilateralTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_bilateral_task_set_property;
    oclass->get_property = ufo_bilateral_task_get_property;
    oclass->finalize = ufo_bilateral_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoBilateralTaskPrivate));
}

static void
ufo_bilateral_task_init(UfoBilateralTask *self)
{
    self->priv = UFO_BILATERAL_TASK_GET_PRIVATE(self);
    self->priv->window = 5;
    self->priv->sigd = 3.0f; //spatial
    self->priv->sigr = 9.0f; //intenstiy
}
