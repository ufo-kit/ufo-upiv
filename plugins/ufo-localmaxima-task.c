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
#include <math.h>
#include <string.h>
#include "ufo-localmaxima-task.h"

static gfloat find_max(gfloat *, guint, guint, guint, guint);
static gfloat snr_max(gfloat *, gint, gint, gint, gint, gint, gint);

struct _UfoLocalmaximaTaskPrivate {
    guint dimx;
    guint dimy;
    guint gcount;
    guint max_detection;
    guint ring_start;
    guint ring_end;
    guint ring_step;
    cl_context context;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoLocalmaximaTask, ufo_localmaxima_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_LOCALMAXIMA_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_LOCALMAXIMA_TASK, UfoLocalmaximaTaskPrivate))

enum {
    PROP_0,
    PROP_RING_START,
    PROP_RING_STEP,
    PROP_RING_END,
    PROP_MAX_DETECTION,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_localmaxima_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_LOCALMAXIMA_TASK, NULL));
}

static void
ufo_localmaxima_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoLocalmaximaTaskPrivate *priv;  
    priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE (task); 
    priv->context = ufo_resources_get_context (resources);
    priv->gcount = 0;
}

static void
ufo_localmaxima_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    UfoLocalmaximaTaskPrivate *priv;  
    UfoRequisition req_in; 

    priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE (task); 

    /* input info */
    ufo_buffer_get_requisition (inputs[0], &req_in);
    priv->dimx = req_in.dims[0];
    priv->dimy = req_in.dims[1];
    priv->gcount++;

    /* output requistion */
    requisition->n_dims = 1;
    requisition->dims[0] = 3*priv->max_detection + 1;
}

static guint
ufo_localmaxima_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_localmaxima_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_localmaxima_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_CPU;
}

static gboolean
ufo_localmaxima_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoLocalmaximaTaskPrivate *priv;  
    UfoBuffer * tmp_buffer;
    gfloat * img_ptr, * out_ptr, * tmp_ptr;
    guint idx, idxx, idxy, ct1 = 0, ct2 = 0, ct3 = 0, depth;
    gfloat glb_max, threshold;
    UfoRequisition req_tmp;

    priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE (task); 

    req_tmp.n_dims = 2;
    req_tmp.dims[0] = priv->dimx;
    req_tmp.dims[1] = priv->dimy;
    tmp_buffer = ufo_buffer_new(&req_tmp, priv->context);

    // copy input image to host and get pointer
    // generate output buffer on host
    img_ptr = ufo_buffer_get_host_array(inputs[0], NULL);
    out_ptr = ufo_buffer_get_host_array(output, NULL);
    tmp_ptr = ufo_buffer_get_host_array(tmp_buffer, NULL);

    // find global maxima
    glb_max = 0.0f;
    for (idx = 0; idx < priv->dimx * priv->dimy; idx++) {
        glb_max = (img_ptr[idx] > glb_max) ? img_ptr[idx] : glb_max;
        tmp_ptr[idx] = 0.0f;
    }

    // find local maxima and save them in temp buffer
    threshold = 0.4f;
    for (idxx = 1; idxx < priv->dimx - 1; idxx++)
    for (idxy = 1; idxy < priv->dimy - 1; idxy++) {
        idx = idxx + priv->dimx * idxy;
        if (( img_ptr[idx] > threshold * glb_max ) &&
            ( img_ptr[idx] > find_max(img_ptr, idxx, idxy, priv->dimx, priv->dimy) ))  {
            tmp_ptr[idx] = 1.0f;
            ct1++;
        }
    }

    gint minx, maxx, miny, maxy, blobx , bloby, label = 0;
    gfloat minv, brate;
    char check[8];

    threshold = 0.1f;
    for (idxy = 1; idxy < priv->dimy - 1; idxy++) 
    for (idxx = 1; idxx < priv->dimx - 1; idxx++) {
        idx = idxx + priv->dimx * idxy;
        if (tmp_ptr[idx] > 0.5f) {

#if 0
            label++;
            depth = 0;

            minv = img_ptr[idx];
            minx = priv->dimx; 
            maxx = 0;
            miny = priv->dimy; 
            maxy = 0;

            explore(img_ptr, out_ptr, label,
                    idxx, idxy, priv->dimx, priv->dimy, 
                    threshold*img_ptr[idx], &depth,
                    &minx, &maxx, &miny, &maxy, &minv); 

            blobx = maxx -minx + 1;
            bloby = maxy - miny + 1;
            brate = (float) depth /blobx/bloby;

            if (depth > 7 && depth < 51 && blobx < 7 && bloby < 7 && blobx > 1 && bloby > 1 && brate > 0.5f) {
                    ct2++;
                    strcpy(check, "ok");

            } else {
                tmp_ptr[idx] = 0.0f;
                strcpy(check, "x");
            }
#endif


            gfloat a1 = snr_max(img_ptr, idxx, idxy, 2, 3, priv->dimx, priv->dimy);

            if (a1 > 2.0f) {
                out_ptr[3*ct3+1] = (gfloat) idxx;
                out_ptr[3*ct3+2] = (gfloat) idxy;
                out_ptr[3*ct3+3] = a1;
                strcpy(check, "ok");
                ct3++;
            } else {
                tmp_ptr[idx] = 0.0f;
                strcpy(check, "x");
            }
            
            out_ptr[0] = ct3;

            /*
             *if (tmp_ptr[idx] > -0.5f) {
             *    printf("(%4d, %4d): threshold = %5e, label = %5d", idxx, idxy, threshold*img_ptr[idx], label);
             *    printf("\t %4d %4d %4d %.4f %2.2f %8s\n", depth, blobx, bloby, brate, a1, check);
             *}
             */

        }
    }

    /*
     *printf("ok %4u: number of pixels: %4u %4u %4u, \tglobal max: %8e\n", 
     *        priv->gcount, ct1, ct2, ct3, glb_max);
     */

    return TRUE;
}

static void
ufo_localmaxima_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoLocalmaximaTaskPrivate *priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_RING_START:
            priv->ring_start = g_value_get_uint(value);
            break;
        case PROP_RING_STEP:
            priv->ring_step = g_value_get_uint(value);
            break;
        case PROP_RING_END:
            priv->ring_end = g_value_get_uint(value);
            break;
        case PROP_MAX_DETECTION:
            priv->max_detection = g_value_get_uint(value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_localmaxima_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoLocalmaximaTaskPrivate *priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE (object);

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
        case PROP_MAX_DETECTION:
            g_value_set_uint (value, priv->max_detection);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_localmaxima_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_localmaxima_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_localmaxima_task_setup;
    iface->get_num_inputs = ufo_localmaxima_task_get_num_inputs;
    iface->get_num_dimensions = ufo_localmaxima_task_get_num_dimensions;
    iface->get_mode = ufo_localmaxima_task_get_mode;
    iface->get_requisition = ufo_localmaxima_task_get_requisition;
    iface->process = ufo_localmaxima_task_process;
}

static void
ufo_localmaxima_task_class_init (UfoLocalmaximaTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_localmaxima_task_set_property;
    oclass->get_property = ufo_localmaxima_task_get_property;
    oclass->finalize = ufo_localmaxima_task_finalize;

    properties[PROP_RING_START] =
        g_param_spec_uint ("ring_start",
                           "give starting radius size",
                           "give starting radius size",
                           1, G_MAXUINT, 5,
                           G_PARAM_READWRITE);

    properties[PROP_RING_STEP] =
        g_param_spec_uint ("ring_step",
                           "Gives ring step",
                           "Gives ring step",
                           1, G_MAXUINT, 2,
                           G_PARAM_READWRITE);

    properties[PROP_RING_END] =
        g_param_spec_uint ("ring_end",
                           "give ending radius size",
                           "give ending radius size",
                           1, G_MAXUINT, 5,
                           G_PARAM_READWRITE);

    properties[PROP_MAX_DETECTION] =
        g_param_spec_uint ("max_detection",
                           "max detection per ring size",
                           "max number of detected rings per ring size",
                           1, G_MAXUINT, 100,
                           G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoLocalmaximaTaskPrivate));
}

static void
ufo_localmaxima_task_init(UfoLocalmaximaTask *self)
{
    self->priv = UFO_LOCALMAXIMA_TASK_GET_PRIVATE(self);
    self->priv->ring_start = 5;
    self->priv->ring_end = 5;
    self->priv->ring_step = 2;
    self->priv->max_detection = 100;
}

//////////////////////////////////////////////////////////////////////////////
// Internal helper functions
//////////////////////////////////////////////////////////////////////////////
static gfloat
find_max(gfloat *ptr, guint idxx, guint idxy, guint dimx, guint dimy) {
    gfloat a, b;

    a = ptr[idxx-1  + dimx * (idxy-1)];
    b = a;

    a = ptr[idxx + dimx * (idxy-1)];
    b = (b > a) ? b : a;

    a = ptr[idxx+1 + dimx * (idxy-1)];
    b = (b > a) ? b : a;

    a = ptr[idxx-1 + dimx * idxy];
    b = (b > a) ? b : a;

    a = ptr[idxx+1 + dimx * idxy];
    b = (b > a) ? b : a;

    a = ptr[idxx-1  + dimx * (idxy+1)];
    b = (b > a) ? b : a;

    a = ptr[idxx + dimx * (idxy+1)];
    b = (b > a) ? b : a;

    a = ptr[idxx+1 + dimx * (idxy+1)];
    b = (b > a) ? b : a;

    return b;
}

static gfloat
local_average(gfloat *ptr, gint x, gint y, gint s, gint dimx, gint dimy) {
    gint pos, minx, maxx, miny, maxy;
    gfloat res = 0.0f;
    
    minx = (x - s < 0) ? 0 : x - s;
    maxx = (x + s + 1 > dimx) ? dimx : x + s + 1;
    miny = (y - s < 0) ? 0 : y - s;
    maxy = (y + s + 1 > dimy) ? dimy : y + s + 1;
    
    for (gint ix = minx; ix < maxx; ix++) {
        for (gint iy = miny; iy < maxy; iy++) {
            pos = ix + dimx * iy;
            res += ptr[pos];
        }
    }

    return res / (s+s+1) / (s+s+1);
}

static gfloat
local_signoise_old(gfloat *ptr, gint x, gint y, gint s, gint dimx, gint dimy) {
    gint pos, x0, x1, x2, x3, y0, y1, y2, y3;
    gfloat sig = 0.0f, muback, noiseback = 0.0f;
    gint ct = 0;
    
    x1 = (x - s < 0) ? 0 : x - s;
    x2 = (x + s + 1 > dimx) ? dimx : x + s + 1;
    y1 = (y - s < 0) ? 0 : y - s;
    y2 = (y + s + 1 > dimy) ? dimy : y + s + 1;

    x0 = (x - s - s < 0) ? 0 : x - s - s;
    y0 = (y - s - s < 0) ? 0 : y - s - s;
    x3 = (x + s + s + 1 > dimx) ? dimx : x + s + s + 1;
    y3 = (y + s + s + 1 > dimy) ? dimy : y + s + s + 1;

    ct = 0;
    for (gint ix = x1; ix < x2; ix++) {
        for (gint iy = y1; iy < y2; iy++) {
            pos = ix + dimx * iy;
            sig += ptr[pos];
            ct ++;
        }
    }
    sig = sig/ct;

    ct = 0;
    muback = 0.0f;
    noiseback = 0.0f;
    for (gint ix = x0; ix < x3; ix++) {
        for (gint iy = y0; iy < y3; iy++) {
            if ( ((ix < x1) || (ix >= x2)) && ((iy < y1) || (iy >= y2)) ) {
                pos = ix + dimx * iy;
                muback += ptr[pos];
                ct++;
            }
        }
    }
    muback = muback / ct;
    for (gint ix = x0; ix < x3; ix++) {
        for (gint iy = y0; iy < y3; iy++) {
            if ( ((ix < x1) || (ix >= x2)) && ((iy < y1) || (iy >= y2)) ) {
                noiseback += (ptr[pos]-muback)*(ptr[pos]-muback);
            }
        }
    }
    noiseback = sqrt(noiseback/(ct-1));

    return sig/noiseback;
}

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

/*
 *  @s: inner size
 *  @t: outer size
 */

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
    if (idx - t < 0) return -1.0f;
    if (idy - t < 0) return -1.0f;
    if (idx + t >= dimx) return -1.0f;
    if (idy + t >= dimy) return -1.0f;

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

    gfloat sig_s = array_sigma(sig, mu_s, ns);
    gfloat sig_b = array_sigma(bak, mu_b, nb);

    // gfloat res = (mu_s - mu_b) / sqrt(nb*sig_s*sig_s + ns*sig_b*sig_b);
    gfloat res = (mu_s - mu_b) / sqrt(ns*sig_b*sig_b + ns*sig_b*sig_b/nb);

    // printf("\t%.2f\t%.4e, %.4e, %.4e, %.4e, %4d, %4d\n", res, mu_s-mu_b, mu_b, sig_s, sig_b, ns, nb);
    return res;
}

static gfloat
snr_max(gfloat *img_ptr, gint idx, gint idy, gint s, gint t, gint dimx, gint dimy) {
    gfloat res, snr_array[10] = {0};

    /*
    snr_array[0] = snr(img_ptr, idx, idy, 2, t, dimx, dimy);
    snr_array[1] = snr(img_ptr, idx - 1, idy, 2, t, dimx, dimy);
    snr_array[2] = snr(img_ptr, idx + 1, idy, 2, t, dimx, dimy);
    snr_array[3] = snr(img_ptr, idx, idy - 1, 2, t, dimx, dimy);
    snr_array[4] = snr(img_ptr, idx, idy + 1, 2, t, dimx, dimy);
    */

    snr_array[5] = snr(img_ptr, idx, idy, 1, t, dimx, dimy);
    snr_array[6] = snr(img_ptr, idx - 1, idy, 1, t, dimx, dimy);
    snr_array[7] = snr(img_ptr, idx + 1, idy, 1, t, dimx, dimy);
    snr_array[8] = snr(img_ptr, idx, idy - 1, 1, t, dimx, dimy);
    snr_array[9] = snr(img_ptr, idx, idy + 1, 1, t, dimx, dimy);

    res = snr_array[0];
    for (gint i = 1; i < 10; i++) {
        res = (res > snr_array[i]) ? res : snr_array[i];
    }

    return res;
}


static void                                                                     
explore (gfloat *src, gfloat *dst, gint label, 
         gint x, gint y, gint dimx, gint dimy, 
         gfloat threshold, guint *depth, gint *minx, gint *maxx, 
         gint *miny, gint *maxy, gfloat* minv)      
{                                                                               
    guint pos;

        if (*depth > 50) return;
        if (x < 0 || y < 0 || x > dimx - 1 || y > dimy - 1) return;

        pos = x + dimx * y;
        if (src[pos] < threshold) return;
        if (abs(dst[pos] - label) < 0.1f) return;

        dst[pos] = label;
        (*depth)++;

        /*
         *printf("\tpos = %d (%4d, %4d)", pos, x, y);
         *printf("\tsrc = %5e %5d %5d\n", src[pos], (int) dst[pos], *depth);
         */

        *minx = ( x < *minx ) ? x : *minx;
        *miny = ( y < *miny ) ? y : *miny;
        *maxx = ( x > *maxx ) ? x : *maxx;
        *maxy = ( y > *maxy ) ? y : *maxy;
        *minv = ( src[pos] < *minv ) ? src[pos] : *minv;

        explore(src, dst, label, x+1, y, dimx, dimy, threshold, depth, minx, maxx, miny, maxy, minv);
        explore(src, dst, label, x-1, y, dimx, dimy, threshold, depth, minx, maxx, miny, maxy, minv);
        explore(src, dst, label, x, y-1, dimx, dimy, threshold, depth, minx, maxx, miny, maxy, minv);
        explore(src, dst, label, x, y+1, dimx, dimy, threshold, depth, minx, maxx, miny, maxy, minv);
}             


