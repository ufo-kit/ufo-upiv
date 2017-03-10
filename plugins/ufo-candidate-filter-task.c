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

#include "ufo-candidate-filter-task.h"
#include "ufo-ring-coordinates.h"
#include <math.h>

struct _UfoCandidateFilterTaskPrivate {
    cl_kernel kernel;
    cl_context context;
    guint ring_start;
    guint ring_step;
    guint ring_end;
    guint ring_current;
    guint scale;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoCandidateFilterTask, ufo_candidate_filter_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_CANDIDATE_FILTER_TASK, UfoCandidateFilterTaskPrivate))

enum {
    PROP_0,
    PROP_RING_START,
    PROP_RING_STEP,
    PROP_RING_END,
    PROP_SCALE,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_candidate_filter_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_CANDIDATE_FILTER_TASK, NULL));
}

static void
ufo_candidate_filter_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoCandidateFilterTaskPrivate *priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE(task);
    
    priv->kernel = ufo_resources_get_kernel(resources, "candidate.cl", NULL, error);

    if (priv->kernel)
        UFO_RESOURCES_CHECK_CLERR(clRetainKernel(priv->kernel));

    priv->context = ufo_resources_get_context(resources);
}

static void
ufo_candidate_filter_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    requisition->n_dims = 1;
    requisition->dims[0] = 1;
}

static guint
ufo_candidate_filter_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_candidate_filter_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return -1;
}

static UfoTaskMode
ufo_candidate_filter_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static int compare_func (gconstpointer a, gconstpointer b)
{
    const UfoRingCoordinate *c = (const UfoRingCoordinate *) a;
    const UfoRingCoordinate *d = (const UfoRingCoordinate *) b;
    if (c->x == d->x) return (int)(c->y - d->y);
    return (int)(c->x - d->x);
}

static GList* sort_candidate (GList *list)
{
    return g_list_sort(list, compare_func);
}

/*
 *static void filter_neighbour (gpointer data, gpointer user_data) {
 *    UfoRingCoordinate *r = (UfoRingCoordinate *) data;
 *    UfoRingCoordinate *s = (UfoRingCoordinate *) user_data;
 *    if (fabs(r->x + r->y - s->x - s->y) < 5.0 && fabs(r->r - s->r) < 2.0)
 *        if (r->contrast < s->contrast)
 *            r->contrast = 0.0;
 *}
 */

static UfoRingCoordinate* get_ring(GList *l) {
    UfoRingCoordinate *r = (UfoRingCoordinate *) l->data;
    return r;
}

static GList* filter_sort_candidate (GList *list)
{
    list = sort_candidate(list);

    // filter local peak in neighbouring fixels, set others to 0
    for (GList *current=list; current; current=current->next) {
        UfoRingCoordinate *r = get_ring(current);
        if (r->contrast == 0.0f) continue;
        for (GList *l = current->next; l; l=l->next) {
            UfoRingCoordinate *s = get_ring(l);
            if (fabs(r->x - s->x) > 5) break;
            if (fabs(r->y - s->y) > 5) continue;
            if (fabs(r->r - s->r) > 5) continue;
            if (r->contrast > s->contrast) {
                s->contrast = 0.0f;
            } else {
                r->contrast = 0.0f;
                break;
            }
        }
    }

    // remove 0's 
    for (GList *current=list; current;) {
        GList *next = current->next;
        UfoRingCoordinate *r = (UfoRingCoordinate*) current->data;
        if (fabs(r->contrast - 0.0f) < 0.000001f) {
            list = g_list_remove_link(list, current);
        }
        current = next;
    }

    return list;
}

static gboolean
ufo_candidate_filter_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoCandidateFilterTaskPrivate *priv;
    UfoRingCoordinate *rings;
    UfoRequisition req;
    int num_cand;
    gfloat *cand_cpu, *out_cpu;
    GList *cand_list = NULL;

    priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE (task);

    cand_cpu = ufo_buffer_get_host_array (inputs[0], NULL);

    rings = (UfoRingCoordinate*) &cand_cpu[2]; 
    num_cand = cand_cpu[0];

    for(int i = 0; i < num_cand; i++) {
        rings[i].r = (priv->ring_start + priv->ring_step * rings[i].r) / priv->scale;
        cand_list = g_list_append(cand_list, (gpointer) &rings[i]);
    }

    // filter candidate neighbours
    cand_list = filter_sort_candidate (cand_list);
    num_cand = g_list_length (cand_list);

    req.n_dims = 1;
    req.dims[0] = 2 + num_cand * sizeof(UfoRingCoordinate) / sizeof (float);
    ufo_buffer_resize(output, &req);

    out_cpu = ufo_buffer_get_host_array(output,NULL);
    out_cpu[0] = num_cand;
    out_cpu[1] = priv->scale;

    rings = (UfoRingCoordinate*) &out_cpu[2];
    for (int i = 0; i < num_cand; i++) 
        rings[i] = * (UfoRingCoordinate *) g_list_nth_data (cand_list, i);

    /*g_message ("number of filtered candidate %d", num_cand);*/
    return TRUE;
}

static void
ufo_candidate_filter_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoCandidateFilterTaskPrivate *priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_RING_START:
            priv->ring_start = g_value_get_uint(value);
            priv->ring_current = priv->ring_start;
            break;
        case PROP_RING_STEP:
            priv->ring_step = g_value_get_uint(value);
            break;
        case PROP_RING_END:
            priv->ring_end = g_value_get_uint(value);
            break;
        case PROP_SCALE:
            priv->scale = g_value_get_uint(value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_candidate_filter_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoCandidateFilterTaskPrivate *priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_RING_START:
            g_value_set_uint (value, priv->ring_start);
            break;
        case PROP_RING_STEP:
            g_value_set_uint (value, priv->ring_step);
            break;
        case PROP_RING_END:
            g_value_set_uint (value, priv->ring_end);
            break;
        case PROP_SCALE:
            g_value_set_uint (value, priv->scale);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_candidate_filter_task_finalize (GObject *object)
{
    UfoCandidateFilterTaskPrivate *priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE (object);

    if(priv->kernel)
        UFO_RESOURCES_CHECK_CLERR(clReleaseKernel(priv->kernel));
    
    G_OBJECT_CLASS (ufo_candidate_filter_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_candidate_filter_task_setup;
    iface->get_num_inputs = ufo_candidate_filter_task_get_num_inputs;
    iface->get_num_dimensions = ufo_candidate_filter_task_get_num_dimensions;
    iface->get_mode = ufo_candidate_filter_task_get_mode;
    iface->get_requisition = ufo_candidate_filter_task_get_requisition;
    iface->process = ufo_candidate_filter_task_process;
}

static void
ufo_candidate_filter_task_class_init (UfoCandidateFilterTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_candidate_filter_task_set_property;
    oclass->get_property = ufo_candidate_filter_task_get_property;
    oclass->finalize = ufo_candidate_filter_task_finalize;

    properties[PROP_RING_START] = 
        g_param_spec_uint ("ring_start",
               "", "",
               1, G_MAXUINT, 5,
               G_PARAM_READWRITE);

    properties[PROP_RING_STEP] = 
        g_param_spec_uint ("ring_step",
               "", "",
               1, G_MAXUINT, 2,
               G_PARAM_READWRITE);

    properties[PROP_RING_END] = 
        g_param_spec_uint ("ring_end",
               "", "",
               1, G_MAXUINT, 5,
               G_PARAM_READWRITE);

    properties[PROP_SCALE] =
        g_param_spec_uint ("scale",
               "Rescale factor", "",
               1, 2, 1,
               G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoCandidateFilterTaskPrivate));
}

static void
ufo_candidate_filter_task_init(UfoCandidateFilterTask *self)
{
    self->priv = UFO_CANDIDATE_FILTER_TASK_GET_PRIVATE(self);
    self->priv->ring_start = 5;
    self->priv->ring_end = 5;
    self->priv->ring_step = 2;
    self->priv->scale = 1;
    self->priv->ring_current = self->priv->ring_start;
}
