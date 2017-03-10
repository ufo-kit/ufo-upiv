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
#include <math.h>

#include "ufo-piv-contrast-task.h"

struct _UfoPivContrastTaskPrivate {
    cl_kernel sigmoid_kernel;
    cl_context context;
    gfloat c1;
    gfloat c2;
    gfloat c3;
    gfloat c4;
    gfloat gamma;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoPivContrastTask, ufo_piv_contrast_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_PIV_CONTRAST_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_PIV_CONTRAST_TASK, UfoPivContrastTaskPrivate))

enum {
    PROP_0,
    PROP_C1,
    PROP_C2,
    PROP_C3,
    PROP_C4,
    PROP_GAMMA,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_piv_contrast_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_PIV_CONTRAST_TASK, NULL));
}

static void
ufo_piv_contrast_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoPivContrastTaskPrivate *priv;

    priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE (task);
    priv->sigmoid_kernel = ufo_resources_get_kernel (resources, "piv_contrast.cl", "sigmoid", error);
    priv->context = ufo_resources_get_context (resources);

    UFO_RESOURCES_CHECK_CLERR (clRetainContext (priv->context));
    if (priv->sigmoid_kernel != NULL)
        UFO_RESOURCES_CHECK_CLERR (clRetainKernel (priv->sigmoid_kernel));
}

static void
ufo_piv_contrast_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    ufo_buffer_get_requisition(inputs[0], requisition);
}

static guint
ufo_piv_contrast_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_piv_contrast_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_piv_contrast_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static float array_min(float *data, int size) 
{
    float res = data[0];
    for (int i = 1; i < size; i++) 
    {
        if (data[i] < res) res = data[i];
    }
    return res;
}

static float array_max(float *data, int size) 
{
    float res = data[0];
    for (int i = 1; i < size; i++) 
    {
        if (data[i] > res) res = data[i];
    }
    return res;
}

static float array_mean(float *data, int size, float low, float high) 
{
    float res = 0.0f;
    unsigned ct = 0;
    for (int i = 0; i < size; i++) {
        if (data[i] > low && data[i] < high) {
            res += data[i];
            ct++;
        }
    }
    return res/ct;
}

static float array_std(float *data, int size, float mean, float low, float high)
{
    float res = 0.0f;
    unsigned ct = 0;
    for (int i = 0; i < size; i++) {
        if (data[i] > low && data[i] < high) {
            res += (data[i]-mean) * (data[i]-mean);
            ct++;
        }
    }
    return sqrt(res/ct);
}

static int prepare_test(gfloat* test, gfloat* input, int offset_x, int offset_y, 
                         int len_x, int len_y, int req_x, int req_y)
{
    int id, id_0;
    len_x= (len_x + offset_x) > req_x ? req_x - offset_x : len_x;
    len_y = (len_y + offset_y) > req_y ? req_y - offset_y : len_y;
    for (int i = 0; i < len_x; i++) {
        for (int j = 0; j < len_y; j++) {
            id = j * len_x + i;
            id_0 = (j+offset_y) * req_x + (i+offset_x);
            test[id] = input[id_0];
        }
    }
    return len_x*len_y;
}

static void approx_bg(gfloat* test, int size, float* m, float* s)
{
    int ct;
    float low_cut, high_cut, std_old;
    float mean, std;
    
    ct = 0;
    std = 0.0f;
    low_cut = array_min(test, size);
    high_cut = array_max(test, size);

    while (ct < 10)
    {
        ct++;
        std_old = std;
        mean = array_mean(test, size, low_cut, high_cut);
        std = array_std(test, size, mean, low_cut, high_cut);
        //printf("mean = %2.8f std = %2.8f; min = %2.8f max = %2.8f\n", mean, std, minv, maxv);

        high_cut = MIN(high_cut, mean + 2*std);
        low_cut = MAX(low_cut, mean - 2*std);

        if (fabs(std - std_old) / std < 0.05)
            break;
    }

    *m = mean;
    *s = std;
}

static gboolean
ufo_piv_contrast_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoPivContrastTaskPrivate *priv;
    UfoGpuNode *node;
    UfoProfiler *profiler;
    cl_command_queue *cmd_queue;
    cl_mem in_device;
    cl_mem out_device;
    
    unsigned offset_x, offset_y, len_x, len_y, req_x, req_y;
    gfloat *in_mem, *test_mem;
    gfloat mean, std, high;
    int n;
 
    req_x = requisition->dims[0];
    req_y = requisition->dims[1];
    offset_x = req_x / 4;
    offset_y = req_y / 4;
    len_x = req_x / 4;
    len_y = req_y / 4;
   
    in_mem = ufo_buffer_get_host_array(inputs[0], NULL);
    test_mem = (gfloat*) g_malloc0 (len_x * len_y * sizeof(gfloat));

    n = prepare_test(test_mem, in_mem, offset_x, offset_y, len_x, len_y, req_x, req_y);
    approx_bg(test_mem, n, &mean, &std);
    g_free(test_mem);
    //g_message("PivContrast: %f, %f", mean, std);

    gfloat *img = ufo_buffer_get_host_array(inputs[0], NULL);
    high = array_max(img, req_x*req_y);

    priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));
    cmd_queue = ufo_gpu_node_get_cmd_queue (node);
    in_device = ufo_buffer_get_device_array(inputs[0], cmd_queue);
    out_device = ufo_buffer_get_device_array(output, cmd_queue);

    gsize global_work_size = requisition->dims[0] * requisition->dims[1];

    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 0, sizeof (cl_mem), &out_device));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 1, sizeof (cl_mem), &in_device));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 2, sizeof (gfloat), &mean));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 3, sizeof (gfloat), &std));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 4, sizeof (gfloat), &high));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 5, sizeof (gfloat), &priv->gamma));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 6, sizeof (gfloat), &priv->c1));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (priv->sigmoid_kernel, 7, sizeof (gfloat), &priv->c2));

    profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));
    ufo_profiler_call (profiler, cmd_queue, priv->sigmoid_kernel, 1, &global_work_size, NULL);

    return TRUE;
}

static void
ufo_piv_contrast_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoPivContrastTaskPrivate *priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_GAMMA:
            priv->gamma = g_value_get_float(value);
            break;
        case PROP_C1:
            priv->c1 = g_value_get_float(value);
            break;
        case PROP_C2:
            priv->c2 = g_value_get_float(value);
            break;
        case PROP_C3:
            priv->c3 = g_value_get_float(value);
            break;
        case PROP_C4:
            priv->c4 = g_value_get_float(value);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_piv_contrast_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoPivContrastTaskPrivate *priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_GAMMA:
            g_value_set_float (value, priv->gamma);
            break;
        case PROP_C1:
            g_value_set_float (value, priv->c1);
            break;
        case PROP_C2:
            g_value_set_float (value, priv->c2);
            break;
        case PROP_C3:
            g_value_set_float (value, priv->c3);
            break;
        case PROP_C4:
            g_value_set_float (value, priv->c4);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_piv_contrast_task_finalize (GObject *object)
{
    UfoPivContrastTaskPrivate *priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE (object);
    
    if (priv->sigmoid_kernel) {
        UFO_RESOURCES_CHECK_CLERR (clReleaseKernel (priv->sigmoid_kernel));
        priv->sigmoid_kernel = NULL;
    }

    G_OBJECT_CLASS (ufo_piv_contrast_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_piv_contrast_task_setup;
    iface->get_num_inputs = ufo_piv_contrast_task_get_num_inputs;
    iface->get_num_dimensions = ufo_piv_contrast_task_get_num_dimensions;
    iface->get_mode = ufo_piv_contrast_task_get_mode;
    iface->get_requisition = ufo_piv_contrast_task_get_requisition;
    iface->process = ufo_piv_contrast_task_process;
}

static void
ufo_piv_contrast_task_class_init (UfoPivContrastTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_piv_contrast_task_set_property;
    oclass->get_property = ufo_piv_contrast_task_get_property;
    oclass->finalize = ufo_piv_contrast_task_finalize;

    properties[PROP_C1] = g_param_spec_float(
                            "c1",
                            "lower cut",
                            "lower cut relative to the background mean in unit of background standard deviation",
                            -5.0f, 5.0f, 0.5f,
                            G_PARAM_READWRITE);

     properties[PROP_C2] = g_param_spec_float(
                            "c2",
                            "higher cut",
                            "higher cut relative to the background mean in unit of background standard deviation",
                            0.0f, 100.0f, 10.0f,
                            G_PARAM_READWRITE);

      properties[PROP_C3] = g_param_spec_float(
                            "c3",
                            "sigmoid shift",
                            "sigmoid shift relative to the background mean in unit of background standard deviation",
                            0.0f, 20.0f, 3.0f,
                            G_PARAM_READWRITE);
 
      properties[PROP_C4] = g_param_spec_float(
                            "c4",
                            "sigmoid width",
                            "sigmoid shift in unit of background standard deviation",
                            0.0f, 10.0f, 1.5f,
                            G_PARAM_READWRITE);

      properties[PROP_GAMMA] = g_param_spec_float(
                            "gamma",
                            "intensity gamma",
                            "intensity gamma",
                            0.0f, 1.0f, 0.3f,
                            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoPivContrastTaskPrivate));
}

static void
ufo_piv_contrast_task_init(UfoPivContrastTask *self)
{
    self->priv = UFO_PIV_CONTRAST_TASK_GET_PRIVATE(self);
    self->priv->c1 = 0.5f;
    self->priv->c2 = 10.0f;
    self->priv->c3 = 4.0f;
    self->priv->c4 = 2.0f;
    self->priv->gamma = 0.3f;
}
