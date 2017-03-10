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

#include <math.h>
#include <string.h>
#include <glib/gprintf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "ufo-azimuthal-test-task.h"
#include "ufo-ring-coordinates.h"

#define SHOWMESSAGE
#undef SHOWMESSAGE

#define MAX_THREAD_NUM 100
#define MAX_HIST_LEN 32
#define FIT_DATA_SIZE 5
#define FIT_DATA_SIZE_HALF (FIT_DATA_SIZE/2)

struct _UfoAzimuthalTestTaskPrivate {
    guint radii_range;
    guint displacement;
    gint thread;
    gfloat thld_azimu;
    gfloat thld_likelihood;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoAzimuthalTestTask, ufo_azimuthal_test_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_AZIMUTHAL_TEST_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_AZIMUTHAL_TEST_TASK, UfoAzimuthalTestTaskPrivate))

enum {
    PROP_0,
    PROP_THREAD,
    PROP_AZIMU,
    PROP_LIKELIHOOD,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_azimuthal_test_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_AZIMUTHAL_TEST_TASK, NULL));
}

static void
ufo_azimuthal_test_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_azimuthal_test_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition (inputs[1], requisition);
}

static guint
ufo_azimuthal_test_task_get_num_inputs (UfoTask *task)
{
    return 2;
}

static guint
ufo_azimuthal_test_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return (input == 0) ?  2 : 1;
}

static UfoTaskMode
ufo_azimuthal_test_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR;
}

static int min (int l, int r)
{
    return (l > r) ? r : l;
}

static int max (int l, int r)
{
    return (l > r) ? l : r;
}

static void
get_coords(int *left, int *right, int *top, int *bot, int rad,
        int center_x, int center_y, int img_width, int img_height)
{
    int l = center_x - rad;
    int r = center_x + rad;
    int t = center_y - rad;
    int b = center_y + rad;
    *left = max (l, 0);
    *right = min (r, img_width - 1);
    *top = max (t, 0);
    *bot = min (b, img_height - 1);
}

static void
compute_histogram(float *histogram, float *h, float *img, int center_x, int center_y,
                        int min_r, int max_r, int img_width, int img_height)
{
    int x0, x1, y0, y1;
    int r,i,j,ii,jj,idx;

    get_coords(&x0, &x1, &y0, &y1, max_r, center_x, center_y, img_width, img_height);

    for (i = x0; i < x1; i++) for (j = y0; j < y1; j++) {
        idx = i + j * img_width;
        ii = i - center_x;
        jj = j - center_y;
        if (abs(ii) < 0.7*min_r && abs(jj) < 0.7*min_r) continue;
        r = (int) roundf(sqrt(ii*ii+jj*jj));

        if (r >= min_r && r <= max_r) {
            histogram[r - min_r] +=  img[idx];
            h[r - min_r] += 1;
        }
    }
    
    for (r = 0; r < max_r - min_r + 1; r++) 
        if (h[r] > 0) histogram[r] /= h[r];
        else histogram[r] = 0.0f;

#if 0
    // smooth histogram 
    memset(h, 0, MAX_HIST_LEN * sizeof(*h));
    for (int i = 1; i < max_r - min_r; i++)
        h[i] = 0.5 * histogram[i] + 0.25 * (histogram[i+1] + histogram[i-1]);
    h[0] = 0.5*histogram[0] + 0.5*histogram[1];
    h[max_r-min_r] = 0.5*histogram[max_r-min_r-1] + 0.5*histogram[max_r-min_r];
    memcpy(histogram, h, MAX_HIST_LEN * sizeof(*h));
#endif
}

static float array_mean(float *array, int max) {
    float res = 0;
    float ct = 0;
    for (int i=0; i < max; i++) {
        if (array[i] != 0.0f) {
            res += array[i];
            ct++;
        }
    }
    return res/ct;
}

static float array_std(float *array, float mean, int max) {
    float res = 0;
    float ct = 0;
    for (int i=0; i < max; i++) {
        if (array[i] != 0.0f) {
            res += (array[i]-mean)*(array[i]-mean);
            ct++;
        }
    }
    return sqrt(res/ct);
}

static int search_peaks(float *data, int *peaks, int wlen)
{
    int peak, i, j, k=0;
    for (i = 0; i < MAX_HIST_LEN - wlen - 1; i++) {
        peak = 0;
        for (j = i; j < i + wlen - 1; j++)
            if (data[j] >= data[j+1]) {
                peak = j;
                break;
            }
        for (; j < i + wlen - 1; j++)
            if (data[j] <= data[j+1]) break;
        if (peak != i + wlen/2 || j != i + wlen - 1 ||
                data[j] == 0.0f || data[i] == 0.0f) continue;
        peaks[k++] = peak;
    }
    return k;
}
typedef struct _gaussian_thread_data
{
    int tid;
    UfoAzimuthalTestTaskPrivate *priv;
    UfoRingCoordinate* ring;
    UfoRingCoordinate* winner;
    float* image;   // input image
    int img_width;
    int img_height;
    int tmp_S;
} gaussian_thread_data;

#ifdef SHOWMESSAGE
static void print_histogram (int center_x, int center_y, int min_r, int max_r, float *histogram)
{
    g_message(                                                          
        "(%3d,%3d)[%2d,%2d]: histogram\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e, %8.2e, %8.2e, %8.2e,\
        \n\t\t%8.2e, %8.2e",
        center_x, center_y, min_r, max_r,
        histogram[0], histogram[1], histogram[2], histogram[3], histogram[4],
        histogram[5], histogram[6], histogram[7], histogram[8], histogram[9],
        histogram[10], histogram[11], histogram[12], histogram[13], histogram[14],
        histogram[15], histogram[16], histogram[17], histogram[18], histogram[19],
        histogram[20], histogram[21], histogram[22], histogram[23], histogram[24],
        histogram[25], histogram[26], histogram[27], histogram[28], histogram[29],
        histogram[30], histogram[31]);
}

static void print_peaks(int center_x, int center_y, int *peaks, float* h)
{
    const char *fmt =
        "(%3d,%3d): peaks = %2d,%2d,%2d,%2d,%2d\t%.2f,%.2f,%.2f,%.2f,%.2f";
    float h0 = peaks[0] == 0.0f ? 0.0 : h[peaks[0]];
    float h1 = peaks[1] == 0.0f ? 0.0 : h[peaks[1]];
    float h2 = peaks[2] == 0.0f ? 0.0 : h[peaks[2]];
    float h3 = peaks[3] == 0.0f ? 0.0 : h[peaks[3]];
    float h4 = peaks[4] == 0.0f ? 0.0 : h[peaks[4]];
    g_message(fmt, center_x, center_y,
            peaks[0], peaks[1], peaks[2], peaks[3], peaks[4],
            h0,h1,h2,h3,h4);
}

static void print_ring (UfoRingCoordinate *ring)
{
    g_message("(%5.1f,%6.1f,%6.1f): likelihood = %10.2f, azimu_test = %+8.4f",
            ring->x, ring->y, ring->r, ring->contrast, ring->intensity);
}
#endif

typedef struct _azimu_data {
    int x;
    int y;
    int r;
    float height;
    float snr;
} azimu_data;

static void gaussian_thread(gpointer data, gpointer user_data)
{
    gaussian_thread_data *parm = (gaussian_thread_data*) data;

    float *image = parm->image;
    int img_width = parm->img_width;
    int img_height = parm->img_height;
    UfoRingCoordinate* ring = parm->ring;
    unsigned radii_range = parm->priv->radii_range;
    int displacement = parm->priv->displacement;

    int ct, num_peaks; 
    int num_test = (2*displacement+1)*(2*displacement+1);
    int peaks[MAX_HIST_LEN];
    float h[MAX_HIST_LEN], h2[MAX_HIST_LEN];
    float histogram[num_test*MAX_HIST_LEN];
    azimu_data azimu[num_test];

    /*UfoRingCoordinate best_rings[4];*/

    // decide range of data to be considered in the fitting procedure
    // Since we aim at outer rings, keep the search range on the inner side strict,
    // and relax outer search range. This will help increase the probablility for
    // finding the peak position of the outer ring.
    int min_r = ring->r - radii_range;      
    int max_r = ring->r + radii_range + MAX(4,(int)(0.2*ring->r));  
    min_r = MAX(min_r, (int)ring->r/2);

    if (max_r - min_r + 1 > MAX_HIST_LEN) max_r = min_r + MAX_HIST_LEN - 1;

    // background
    float bk, bks;
    memset (h, 0, MAX_HIST_LEN * sizeof(*h));
    memset (h2, 0, MAX_HIST_LEN * sizeof(*h2));
    compute_histogram(h2, h, image, ring->x, ring->y, max_r, max_r+20, img_width, img_height);
    bk = array_mean(h2, MAX_HIST_LEN);
    bks = array_std(h2, bk, MAX_HIST_LEN);

#ifdef SHOWMESSAGE
    g_message(" ");
    print_ring(ring);
    g_message("%f %f", bk, bks);
#endif

    // compute azimuthal histogram
    ct = -1;
    for (int j = - displacement; j < displacement + 1; j++) {
    for (int k = - displacement; k < displacement + 1; k++) {
        float *hist = histogram + (++ct) * MAX_HIST_LEN;
        int center_x = ring->x + j;
        int center_y = ring->y + k;

        memset (hist, 0, MAX_HIST_LEN * sizeof(*hist));
        memset (h, 0, MAX_HIST_LEN * sizeof(*h));
        memset (peaks, 0, MAX_HIST_LEN * sizeof(*peaks));
        azimu[ct] = (azimu_data) {center_x, center_y, 0, 0.0f, 0.0f};

        compute_histogram(hist, h, image, 
                center_x, center_y, min_r, max_r, img_width, img_height);
        num_peaks = search_peaks(hist, peaks, FIT_DATA_SIZE);

        // highest peak
        float p = -1.0f;
        for (int i = 0; i < num_peaks; i++) {
            if (hist[peaks[i]] > p) {
                p = hist[peaks[i]];
                azimu[ct].r = peaks[i];
            }
        }

        if (azimu[ct].r > 0) {
            for (int i = azimu[ct].r; 
                    i < azimu[ct].r + FIT_DATA_SIZE; i++) {
                azimu[ct].snr += hist[i] - bk;
            }
            azimu[ct].snr /=  (bks * FIT_DATA_SIZE);
            azimu[ct].height = hist[azimu[ct].r];
            azimu[ct].r += min_r;
        }

#ifdef SHOWMESSAGE
        print_histogram(center_x,center_y,min_r,max_r,hist);
        print_histogram(center_x,center_y,max_r,max_r+20,h2);
        print_peaks(center_x, center_y, peaks, hist);
        g_message("(%3d,%3d): %2d %6.2f %6.2f",
            azimu[ct].x, azimu[ct].y, azimu[ct].r, azimu[ct].height, azimu[ct].snr);
#endif
    }}

    ct = 0;
    UfoRingCoordinate ring0 = {0,0,0,ring->contrast,0};
    for (int i = 0; i < num_test; i++) {
        if (abs(azimu[i].r - (int)ring->r) <= 3) {
            ring0.x += azimu[i].x;
            ring0.y += azimu[i].y;
            ring0.r += azimu[i].r;
            ring0.intensity += azimu[i].snr;
            ct++;
        }
    }

    if (ct>0) {
        ring0.x /= ct;
        ring0.y /= ct;
        ring0.r /= ct;
        ring0.intensity /= ct;
        *ring = ring0;
    }

#ifdef SHOWMESSAGE
    print_ring (&ring0);
#endif
}

static gboolean
ufo_azimuthal_test_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoAzimuthalTestTaskPrivate *priv;
    UfoRequisition req;

    priv = UFO_AZIMUTHAL_TEST_TASK_GET_PRIVATE (task);
    float *image = ufo_buffer_get_host_array (inputs[0], NULL);
    float *cand_stream = ufo_buffer_get_host_array (inputs[1], NULL);
    guint num_cand = (guint) *cand_stream;
    guint scale = (guint) cand_stream[1];
    UfoRingCoordinate *cand = (UfoRingCoordinate*) &cand_stream[2];

    ufo_buffer_get_requisition(inputs[0], &req);
    int img_width = req.dims[0];
    int img_height = req.dims[1];

    if (scale == 2)
        for (unsigned i = 0; i < num_cand; i++) {
            cand[i].x *= scale;
            cand[i].y *= scale;
            cand[i].r *= scale;
        }

    GError *err;
    GThreadPool *thread_pool = NULL;
    thread_pool = g_thread_pool_new((GFunc) gaussian_thread, NULL, priv->thread, TRUE,&err);
    gaussian_thread_data thread_data[num_cand];

    for (unsigned i = 0; i < num_cand; i++) {
        thread_data[i].tid = i;
        thread_data[i].priv = priv;
        thread_data[i].ring = &cand[i];
        thread_data[i].image = image;
        thread_data[i].img_width = img_width;
        thread_data[i].img_height = img_height;
        g_thread_pool_push(thread_pool, &thread_data[i], &err); 
    }

    g_thread_pool_free(thread_pool, FALSE, TRUE);

    float *res = ufo_buffer_get_host_array(output, NULL);
    UfoRingCoordinate *rings = (UfoRingCoordinate*) &res[2];

    int num = 0;
    for(unsigned i=0; i < num_cand;i++) {
        if ((cand[i].contrast > priv->thld_likelihood
                && cand[i].intensity  > priv->thld_azimu)
            || (cand[i].contrast < priv->thld_likelihood
                && cand[i].intensity  > 2.0 * priv->thld_azimu)
            || (cand[i].intensity < priv->thld_azimu 
                && cand[i].contrast > 2.0 * priv->thld_likelihood)) {

            rings[num] = cand[i];
            num++;
        }
    }
    res[0] = num;
    res[1] = 1;
    
    /*req.n_dims = 1;*/
    /*req.dims[0] = 2 + num * sizeof(UfoRingCoordinate) / sizeof(float);*/
    /*ufo_buffer_resize(output, &req);*/

    return TRUE;
}

static void
ufo_azimuthal_test_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoAzimuthalTestTaskPrivate *priv = UFO_AZIMUTHAL_TEST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_THREAD:
            priv->thread = g_value_get_uint(value);
            break;
        case PROP_AZIMU:
            priv->thld_azimu = g_value_get_float(value);
            break;
        case PROP_LIKELIHOOD:
            priv->thld_likelihood = g_value_get_float(value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_azimuthal_test_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoAzimuthalTestTaskPrivate *priv = UFO_AZIMUTHAL_TEST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_THREAD:
            g_value_set_uint(value, priv->thread);
            break;
        case PROP_AZIMU:
            g_value_set_float(value, priv->thld_azimu);
            break;
        case PROP_LIKELIHOOD:
            g_value_set_float(value, priv->thld_likelihood);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_azimuthal_test_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_azimuthal_test_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_azimuthal_test_task_setup;
    iface->get_num_inputs = ufo_azimuthal_test_task_get_num_inputs;
    iface->get_num_dimensions = ufo_azimuthal_test_task_get_num_dimensions;
    iface->get_mode = ufo_azimuthal_test_task_get_mode;
    iface->get_requisition = ufo_azimuthal_test_task_get_requisition;
    iface->process = ufo_azimuthal_test_task_process;
}

static void
ufo_azimuthal_test_task_class_init (UfoAzimuthalTestTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_azimuthal_test_task_set_property;
    oclass->get_property = ufo_azimuthal_test_task_get_property;
    oclass->finalize = ufo_azimuthal_test_task_finalize;

    properties[PROP_THREAD] = 
        g_param_spec_uint ("thread",
               "", "",
               1, MAX_THREAD_NUM, 8,
               G_PARAM_READWRITE);

    properties[PROP_AZIMU] = 
        g_param_spec_float ("azimu_thld",
               "", "",
               G_MINFLOAT, G_MAXFLOAT, 5.0f,
               G_PARAM_READWRITE);

    properties[PROP_LIKELIHOOD] = 
        g_param_spec_float ("likelihood_thld",
               "", "",
               G_MINFLOAT, G_MAXFLOAT, 350.0f,
               G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoAzimuthalTestTaskPrivate));
}

static void
ufo_azimuthal_test_task_init(UfoAzimuthalTestTask *self)
{
    self->priv = UFO_AZIMUTHAL_TEST_TASK_GET_PRIVATE(self);
    self->priv->radii_range = 6;
    self->priv->displacement = 1;
    self->priv->thread = 8;
    self->priv->thld_azimu = 5.0f;
    self->priv->thld_likelihood = 350.0f;
}
/*
 *
 *     initialize Gaussian fitting solver
 *    int nd = FIT_DATA_SIZE;
 *    int np = 3;
 *
 *    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
 *    gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, nd, np);
 *    gsl_matrix *parm_covar = gsl_matrix_alloc (np, np);
 *
 *    gsl_multifit_function_fdf f;
 *    f.f = &gaussian_f;
 *    f.df = &gaussian_df;
 *    f.fdf = &gaussian_fdf;
 *    f.n = nd;
 *    f.p = np;
 *
 *
 *    // initialize the list for best rings
 *    GList *br_list = NULL;
 *    for (int i = 0; i < 4; i++) {
 *        best_rings[i] = r0;
 *        br_list = g_list_prepend(br_list, &best_rings[i]);
 *    }
 *
 *    for (int j = - displacement; j < displacement + 1; j++) {
 *    for (int k = - displacement; k < displacement + 1; k++) {
 *        memset (histogram, 0, MAX_HIST_LEN * sizeof(*histogram));
 *        memset (h, 0, MAX_HIST_LEN * sizeof(*h));
 *
 *        int center_x = ring->x + j;
 *        int center_y = ring->y + k;
 *
 *        [>compute_histogram(histogram, h, image, center_x, center_y, <]
 *                    [>min_r, max_r, img_width, img_height);<]
 *
 *
 *#if 1 
 *        memcpy (h, histogram, MAX_HIST_LEN * sizeof(*h));
 *#else
 *        #endif
 *
 *        memset (peaks, 0, MAX_HIST_LEN * sizeof(*peaks));
 *        search_peaks(h, peaks, FIT_DATA_SIZE-1);
 *
 *        // choose the highest peak
 *        float peak_h = 0;
 *        int peak = -1;
 *        int i = -1;
 *        while (peaks[++i] > 0 && i < MAX_HIST_LEN) {
 *            if (histogram[peaks[i]] > peak_h) {
 *                peak = peaks[i];
 *                peak_h = histogram[peak];
 *            }
 *        }
 *
 *
 *        if (peak == -1) continue; // no peak found
 *        if (abs(peak+min_r - ring->r) > 4) continue; // peak off the expectation
 *
 *        // set fitting data
 *        struct fitting_data d = 
 *                {FIT_DATA_SIZE, &histogram[peak-FIT_DATA_SIZE_HALF]};
 *        f.params = &d;
 *
 *        // set initial parameter guess
 *        double parm_init[] = {histogram[peak], FIT_DATA_SIZE_HALF, 1};
 *        gsl_vector_view parm_view = gsl_vector_view_array(parm_init, 3);
 *        gsl_multifit_fdfsolver_set(s, &f, &parm_view.vector);
 *
 *        // fitting
 *        int iter=0;
 *        do{
 *            iter++;
 *            status = gsl_multifit_fdfsolver_iterate(s);
 *            if(status) break;
 *            status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
 *        } while (status == GSL_CONTINUE && iter < 50);
 *
 *        // chi square
 *        gsl_multifit_covar (s->J, 0.0, parm_covar);
 *        float chi = gsl_blas_dnrm2(s->f);
 *        float dof = nd - np;
 *        float c = chi / sqrt(dof);
 *
 *#define FIT(i) gsl_vector_get(s->x, i)
 *#define ERR(i) sqrt(gsl_matrix_get(parm_covar,i,i))
 *        float parm_A = FIT(0);
 *        float parm_sig = FIT(2);
 *        float err_mu = c * ERR(1);
 *        float ratio  = fabs(parm_A/parm_sig);
 *
 *        // discard if error of mu is large
 *        [>if (err_mu > 0.3f) continue;<]
 *            
 *        // maintain a list of best fitted rings, reverse ordered
 *        // replace the first one in list if ratio is higher
 *        r = (UfoRingCoordinate*) br_list->data;
 *        if (r->intensity < ratio) {
 *            r->x = center_x;
 *            r->y = center_y;
 *            r->r = peak + min_r;
 *            r->contrast = ring->contrast;
 *            r->intensity = ratio;
 *        }
 *        br_list = g_list_sort(br_list, cmp_intensity_reverse);
 *
 *        float parm_mu = FIT(1);
 *        float err_A = c * ERR(0);
 *        float err_sig = c * ERR(2);
 *#ifdef SHOWMESSAGE
 *        g_message(
 *            "(%3d,%3d,%3d(%3d)): %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,\tA =%6.3f(%6.3f), mu =%6.3f(%6.3f), sig =%6.3f(%6.3f)  %6.3f",
 *            center_x, center_y, peak+min_r, (int)ring->r,
 *            histogram[peak-3], histogram[peak-2],
 *            histogram[peak-1], histogram[peak-0],
 *            histogram[peak+1], histogram[peak+2],
 *            histogram[peak+3],
 *            parm_A, err_A, parm_mu, err_mu, parm_sig, err_sig, ratio);
 *#endif
 *#ifndef SHOWMESSAGE
 *#endif
 *    }} // displacement loop
 *
 *    br_list = g_list_reverse(br_list);
 *#ifndef SHOWMESSAGE
 *    for (GList *l=br_list; l; l=l->next) {
 *        r = (UfoRingCoordinate*)(l->data);
 *        g_message("(%3d,%3d,%3d): likelihood = %f, azimu_test = %f",
 *                (int)r->x, (int)r->y, (int)r->r, r->contrast, r->intensity);
 *    }
 *#endif
 *
 *    int number = 0;
 *    *ring = r0;
 *    for (GList *l=br_list; l; l=l->next) {
 *        r = (UfoRingCoordinate*)(l->data);
 *        if (r->intensity == 0) break;
 *        number++;
 *        ring->intensity += r->intensity;
 *        ring->x += r->x * r->intensity;
 *        ring->y += r->y * r->intensity;
 *        ring->r += r->r * r->intensity;
 *    }
 *    if (number > 0) {
 *        ring->x /= ring->intensity;
 *        ring->y /= ring->intensity;
 *        ring->r /= ring->intensity;
 *        ring->intensity /= number;
 *    } else {
 *        r0.contrast = ring->contrast;
 *        *ring = r0;
 *    }
 *
 *    r0.contrast = ring->contrast;
 *    ring = r0;
 *
 *    min_r = max_r;
 *    max_r = 0;
 *
 *    br_list = g_list_reverse(br_list);
 *    float max_azimu = ((UfoRingCoordinate*)(br_list->data))->intensity;
 *
 *    [>for (GList *l=br_list; l; l=l->next) {<]
 *        [>r = (UfoRingCoordinate*)(l->data);<]
 *
 *        [>if (r->intensity < 0.9*max_azimu) break;<]
 *        [>number++;<]
 *        
 *        [>// find max_r and min_r from br_list<]
 *        [>max_r = max_r > (int)r->r ? max_r : r->r;<]
 *        [>min_r = min_r < (int)r->r ? min_r : r->r;<]
 *
 *    [>}<]
 *
 *    [>if (max_r - min_r < 3) {<]
 *        [>ring->x /= ring->intensity;<]
 *        [>ring->y /= ring->intensity;<]
 *        [>ring->r /= ring->intensity;<]
 *        [>ring->intensity /= number;<]
 *    [>} else {<]
 *        [>*ring = r0;<]
 *    [>}<]
 *
 *    g_message("(%6.2f,%6.2f,%6.2f): likelihood = %f, azimu_test = %f",
 *            ring->x, ring->y, ring->r, ring->contrast, ring->intensity);
 *    gsl_matrix_free (parm_covar);
 *
 *
 *struct fitting_data 
 *{
 *    size_t n;
 *    float *y;
 *};
 *
 *static int gaussian_f (const gsl_vector *x, void *data, gsl_vector *f)
 *{
 *    size_t n = ((struct fitting_data *) data)->n;
 *    float *y = ((struct fitting_data *) data)->y;
 *
 *    float A = gsl_vector_get(x, 0);
 *    float mu = gsl_vector_get(x, 1);
 *    float sig = gsl_vector_get(x, 2);
 *
 *    for (size_t i = 0; i < n; i++)
 *    {
 *        float Yi = A * exp ( - (i - mu) * (i - mu) / (2.0f * sig * sig));
 *        gsl_vector_set (f, i, Yi - y[i]);
 *    }
 *
 *    return GSL_SUCCESS;
 *}
 *
 *static int gaussian_df (const gsl_vector *x, void *data, gsl_matrix *J)
 *{
 *    size_t n = ((struct fitting_data *) data)->n;
 *
 *    float A = gsl_vector_get(x, 0);
 *    float mu = gsl_vector_get(x, 1);
 *    float sig = gsl_vector_get(x, 2);
 *
 *    for (size_t i = 0; i < n; i++)
 *    {
 *        double t = i;
 *        double e = exp ( - (t - mu)*(t - mu) / (2.0f*sig*sig) );
 *        gsl_matrix_set (J, i, 0, e);
 *        gsl_matrix_set (J, i, 1, A * e * (t-mu) / (sig*sig) );
 *        gsl_matrix_set (J, i, 2, A * e * (t-mu)*(t-mu) / (sig*sig*sig) );
 *    }
 *
 *    return GSL_SUCCESS;
 *}
 *
 *static int gaussian_fdf (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
 *{
 *    gaussian_f (x, data, f);
 *    gaussian_df (x, data, J);
 *    return GSL_SUCCESS;
 *}
 *
 *static gint cmp_intensity_reverse (gconstpointer a, gconstpointer b)
 *{
 *    const UfoRingCoordinate *ra = (const UfoRingCoordinate *) a;
 *    const UfoRingCoordinate *rb = (const UfoRingCoordinate *) b;
 *    return ra->intensity > rb->intensity ? 1 : -1;
 *}
 *
 */

/*
 *static void compute_background(float *mean, float *std, 
 *                        float *img, int center_x, int center_y,
 *                        int min_r, int max_r, int img_width, int img_height)
 *{
 *    int x0, x1, y0, y1;
 *    int i,j,ii,jj,k,idx;
 *    float bk[300], mm=0, ss=0;
 *    get_coords(&x0, &x1, &y0, &y1, max_r, center_x, center_y, img_width, img_height);
 *    k=0;
 *    for (i = x0; i < x1; i++) for (j = y0; j < y1 && k < 300; j++) {
 *    [>for (i = x0; i < x1; i++) for (j = y0; j < y1; j++) {<]
 *        ii = i - center_x;
 *        jj = j - center_y;
 *        if (abs(ii) < min_r && abs(jj) < min_r) continue;
 *        idx = i + j * img_width;
 *        if (img[idx] == 0.0f) continue;
 *        bk[k] = img[idx];
 *        mm += bk[k];
 *        [>g_message("%f", bk[k]);<]
 *        k++;
 *    }
 *    if (k == 0) {
 *        *mean = 0;
 *        *std = 0;
 *    } else {
 *        mm /= k;
 *        for (i = 0; i < 300; i++) ss += (bk[i]-mm)*(bk[i]-mm);
 *        *mean = mm;
 *        *std = sqrt(ss/k);
 *    }
 *    [>g_message("%d, %f %f %f", k, *mean, *std, ss);<]
 *}
 *
 */
