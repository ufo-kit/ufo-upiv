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

#ifdef HAVE_AMD
#include <clFFT.h>
#else
#include "oclFFT.h"
#endif

#include <math.h>
#include <assert.h>

#include "ufo-fftconvolution-task.h"

/**
 * SECTION:ufo-fftconvolution-task
 * @Short_description: FFT Convolve 2D Image
 * @Title: fftconvolution
 *
 */

struct _UfoFftconvolutionTaskPrivate {
    cl_context context;
    cl_kernel kernel_padzero;
    cl_kernel kernel_fftmult;
    cl_kernel kernel_fftpack;
    cl_kernel kernel_fillzero;
    cl_command_queue cmd_queue;
    gsize width, height, width_p, height_p, batch_size;
    UfoProfiler *profiler;
    clFFT_Plan fft_plan;
    UfoBuffer *fft_buffer_1, *fft_buffer_2;
    cl_mem fft_mem_1, fft_mem_2;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoFftconvolutionTask, ufo_fftconvolution_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_FFTCONVOLUTION_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_FFTCONVOLUTION_TASK, UfoFftconvolutionTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_fftconvolution_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_FFTCONVOLUTION_TASK, NULL));
}

static void
ufo_fftconvolution_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
    UfoFftconvolutionTaskPrivate *priv;
    UfoGpuNode *node;

    priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (task);
    node = UFO_GPU_NODE (ufo_task_node_get_proc_node (UFO_TASK_NODE (task)));

    priv->profiler = ufo_task_node_get_profiler (UFO_TASK_NODE (task));

    priv->kernel_padzero  = ufo_resources_get_kernel (resources, "fft.cl", "fft_spread", error);
    priv->kernel_fftpack  = ufo_resources_get_kernel (resources, "fft.cl", "fft_pack", error);
    priv->kernel_fftmult  = ufo_resources_get_kernel (resources, "mult.cl", "mult", error);
    priv->kernel_fillzero = ufo_resources_get_kernel (resources, "fillzero.cl", "fillzero", error);

    priv->context = ufo_resources_get_context (resources);
    priv->cmd_queue = ufo_gpu_node_get_cmd_queue (node);

    priv->fft_plan = NULL;

    UFO_RESOURCES_CHECK_CLERR (clRetainContext (priv->context));
}

static void
ufo_fftconvolution_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    UfoFftconvolutionTaskPrivate *priv;
    UfoRequisition req0, req1;

    priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (task);
    ufo_buffer_get_requisition (inputs[0], &req0);
    /*printf("input 0  \t%d %d %d\n", req0.n_dims, req0.dims[0], req0.dims[1]);*/
    ufo_buffer_get_requisition (inputs[1], &req1);
    /*printf("input 1  \t%d %d %d %d\n", req1.n_dims, req1.dims[0], req1.dims[1], req1.dims[2]);*/

    priv->width        = req0.dims[0];
    priv->height       = req0.dims[1];
    priv->width_p      = req1.dims[0];
    priv->height_p     = req1.dims[1];

    g_assert_cmpint(priv->width, ==, priv->width_p);
    g_assert_cmpint(priv->height, ==, priv->height_p);

    if (req1.n_dims == 3)
        priv->batch_size  = req1.dims[2];
    else
        priv->batch_size = 1;
    
    requisition->n_dims = 3;
    requisition->dims[0] = priv->width;
    requisition->dims[1] = priv->height;
    requisition->dims[2] = priv->batch_size;

    // create Plan here when it does not exists yet
    if (priv->fft_plan == NULL) {
        cl_int cl_err;
        UfoRequisition req_fft;
        clFFT_Dim3 fft_size;
        clFFT_Dimension dimension = clFFT_2D;

        req_fft.n_dims = 2;
        req_fft.dims[0] = 2*priv->width_p;
        req_fft.dims[1] = priv->height_p;

        fft_size.x = req_fft.dims[0] / 2;
        fft_size.y = req_fft.dims[1];
        fft_size.z = 1;

        priv->fft_buffer_1 = ufo_buffer_new(&req_fft, priv->context);
        priv->fft_buffer_2 = ufo_buffer_new(&req_fft, priv->context);
        priv->fft_mem_1 = ufo_buffer_get_device_array(priv->fft_buffer_1, priv->cmd_queue);
        priv->fft_mem_2 = ufo_buffer_get_device_array(priv->fft_buffer_2, priv->cmd_queue);

        priv->fft_plan = clFFT_CreatePlan (priv->context,
                                     fft_size,
                                     dimension,
                                     clFFT_InterleavedComplexFormat,
                                     &cl_err);
    }
}

static guint
ufo_fftconvolution_task_get_num_inputs (UfoTask *task)
{
    return 2;
}

static guint
ufo_fftconvolution_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    if (input == 1) {
        return 3;
    }
    return 2;
}

static UfoTaskMode
ufo_fftconvolution_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

void ufo_buffer_copy_page(UfoBuffer *src, 
                          UfoBuffer *dst, 
                          const guint n, 
                          cl_command_queue cmd_queue) 
{
    UfoRequisition src_req, dst_req;
    
    g_return_if_fail ( UFO_IS_BUFFER (src) && UFO_IS_BUFFER (dst) );

    cl_mem src_mem  = ufo_buffer_get_device_array(src, cmd_queue);
    cl_mem dst_mem  = ufo_buffer_get_device_array(dst, cmd_queue);

    ufo_buffer_get_requisition(src, &src_req);
    ufo_buffer_get_requisition(dst, &dst_req);

    /*
     *printf("copy page src \t%d %d %d\n", src_req.n_dims, src_req.dims[0], src_req.dims[1]);
     *printf("copy page dst \t%d %d %d\n", dst_req.n_dims, dst_req.dims[0], dst_req.dims[1]);
     */

    const size_t src_origin[] = {0, 0, n};
    const size_t dst_origin[] = {0, 0, 0};
    const size_t region[]     = { 4*src_req.dims[0], src_req.dims[1], 1 };
    size_t src_row_pitch      = src_req.dims[0] * 4;
    size_t src_slice_pitch    = src_req.dims[1];
    size_t dst_row_pitch      = dst_req.dims[0] * 4;
    size_t dst_slice_pitch    = dst_req.dims[1];

    UFO_RESOURCES_CHECK_CLERR ( 
             clEnqueueCopyBufferRect( cmd_queue,
                                      src_mem, dst_mem,
                                      src_origin, dst_origin, region,
                                      src_row_pitch, src_slice_pitch,
                                      dst_row_pitch, dst_slice_pitch,
                                      0,NULL,NULL));
}

void ufo_buffer_past_page(cl_command_queue cmd_queue,
                          UfoBuffer *src, 
                          UfoBuffer *dst, 
                          const guint n) 
{
    UfoRequisition src_req, dst_req;
    
    g_return_if_fail ( UFO_IS_BUFFER (src) && UFO_IS_BUFFER (dst) );

    cl_mem src_mem  = ufo_buffer_get_device_array(src, cmd_queue);
    cl_mem dst_mem  = ufo_buffer_get_device_array(dst, cmd_queue);

    ufo_buffer_get_requisition(src, &src_req);
    ufo_buffer_get_requisition(dst, &dst_req);

    /*
     *printf("paste page src \t%d %d %d\n", src_req.n_dims, src_req.dims[0], src_req.dims[1]);
     *printf("paste page dst \t%d %d %d\n", dst_req.n_dims, dst_req.dims[0], dst_req.dims[1]);
     */

    const size_t src_origin[] = {0, 0, 0};
    const size_t dst_origin[] = {0, 0, n};
    const size_t region[]     = { 4*dst_req.dims[0], dst_req.dims[1], 1 };
    size_t src_row_pitch      = src_req.dims[0] * 4;
    size_t src_slice_pitch    = src_req.dims[1];
    size_t dst_row_pitch      = dst_req.dims[0] * 4;
    size_t dst_slice_pitch    = dst_req.dims[1];

    UFO_RESOURCES_CHECK_CLERR ( 
             clEnqueueCopyBufferRect( cmd_queue,
                                      src_mem, dst_mem,
                                      src_origin, dst_origin, region,
                                      src_row_pitch, src_slice_pitch,
                                      dst_row_pitch, dst_slice_pitch,
                                      0,NULL,NULL));
}

cl_int ufo_buffer_fill_zero(cl_command_queue cmd_queue, cl_kernel *kernel, UfoBuffer *buf)
{
    g_return_if_fail ( UFO_IS_BUFFER (buf) );

    UfoRequisition req;
    ufo_buffer_get_requisition(buf, &req);

    gsize global_work_size[] = {req.dims[0], req.dims[1]};

    cl_mem mem = ufo_buffer_get_device_array(buf, cmd_queue);

    UFO_RESOURCES_CHECK_CLERR (
            clSetKernelArg (*kernel, 0, sizeof (cl_mem), &mem));
    UFO_RESOURCES_CHECK_CLERR (
            clEnqueueNDRangeKernel(cmd_queue, *kernel,
                                   2, NULL,
                                   global_work_size, NULL,
                                   0, NULL, NULL));
}

void ufo_buffer_pad_zero(cl_command_queue cmd_queue, 
                         cl_kernel *kernel, 
                         UfoBuffer *src_buf, 
                         UfoBuffer *dst_buf)
{
    g_return_if_fail ( UFO_IS_BUFFER (src_buf) && UFO_IS_BUFFER (dst_buf) );

    cl_mem src_mem = ufo_buffer_get_device_array(src_buf, cmd_queue);
    cl_mem dst_mem = ufo_buffer_get_device_array(dst_buf, cmd_queue);

    UfoRequisition src_req, dst_req;
    ufo_buffer_get_requisition(src_buf, &src_req);
    ufo_buffer_get_requisition(dst_buf, &dst_req);

    /*
     *printf("pad zeros src \t%d %d %d\n", src_req.n_dims, src_req.dims[0], src_req.dims[1]);
     *printf("pad zeros dst \t%d %d %d\n", dst_req.n_dims, dst_req.dims[0], dst_req.dims[1]);
     */

    gsize global_work_size[] = {dst_req.dims[0]/2, dst_req.dims[1]};

    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg ( 
                                    *kernel, 0, 
                                    sizeof (cl_mem), (gpointer) &dst_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg ( 
                                    *kernel, 1, 
                                    sizeof (cl_mem), (gpointer) &src_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg ( 
                                    *kernel, 2, 
                                    sizeof (cl_int), &src_req.dims[0]));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg ( 
                                    *kernel, 3, 
                                    sizeof (cl_int), &src_req.dims[1]));
    UFO_RESOURCES_CHECK_CLERR (clEnqueueNDRangeKernel (
                                    cmd_queue, *kernel,
                                    2, NULL, global_work_size, NULL,
                                    0, NULL, NULL));
}

void ufo_buffer_fftpack (cl_command_queue cmd_queue, 
                         cl_kernel *kernel, 
                         UfoBuffer *src_buf, 
                         UfoBuffer *dst_buf)
{
    g_return_if_fail ( UFO_IS_BUFFER (src_buf) && UFO_IS_BUFFER (dst_buf) );

    cl_mem src_mem = ufo_buffer_get_device_array(src_buf, cmd_queue);
    cl_mem dst_mem = ufo_buffer_get_device_array(dst_buf, cmd_queue);
    
    UfoRequisition src_req, dst_req;
    ufo_buffer_get_requisition(src_buf, &src_req);
    ufo_buffer_get_requisition(dst_buf, &dst_req);

/*
 *    printf("fft pack \t%d %d %d\n", dst_req.n_dims, dst_req.dims[0], dst_req.dims[1]);
 *
 */
    gsize global_work_size[] = {dst_req.dims[0], dst_req.dims[1]};

    gfloat scale = 1.0f/ ((gfloat) src_req.dims[0]*src_req.dims[0]);
    //scale *= 1.0f/ ((gfloat) src_req.dims[0]*src_req.dims[0]);
    //gfloat scale = 1.0f;

    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (*kernel, 0, 
                                    sizeof (cl_mem), (gpointer) &src_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (*kernel, 1,
                                    sizeof (cl_mem), (gpointer) &dst_mem));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (*kernel, 2, 
                                    sizeof (cl_int), &dst_req.dims[0]));
    UFO_RESOURCES_CHECK_CLERR (clSetKernelArg (*kernel, 3, 
                                    sizeof (gfloat), &scale));

    UFO_RESOURCES_CHECK_CLERR (clEnqueueNDRangeKernel (cmd_queue,
                                                       *kernel,
                                                       2, NULL, global_work_size, NULL,
                                                       0, NULL, NULL));
}

static gboolean
ufo_fftconvolution_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    UfoRequisition req_tmp;
    UfoBuffer *tmp_buffer;
    gsize global_work_size[2]; 
    cl_event event;

    UfoFftconvolutionTaskPrivate *priv;
    priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (task);

    global_work_size[0] = priv->width_p;
    global_work_size[1] = priv->height_p;

    req_tmp.n_dims = 2;
    req_tmp.dims[0] = priv->width_p;
    req_tmp.dims[1] = priv->height_p;

    
    /*printf("fft conv req_tmp \t%d %d\n", req_tmp.dims[0], req_tmp.dims[1]);*/
    /*printf("fft conv req_fft \t%d %d\n", req_fft.dims[0], req_fft.dims[1]);*/
    tmp_buffer      = ufo_buffer_new(&req_tmp, priv->context);
    //img_mem  = ufo_buffer_get_device_array(inputs[0], priv->cmd_queue);
    //ker_mem  = ufo_buffer_get_device_array(inputs[1], priv->cmd_queue);
    //out_mem  = ufo_buffer_get_device_array(output, priv->cmd_queue);
    //tmp_mem  = ufo_buffer_get_device_array(tmp_buffer, priv->cmd_queue);
    // fill tmp buffer with zeros 
    ufo_buffer_fill_zero( priv->cmd_queue, &( priv->kernel_fillzero ), tmp_buffer );

    // prepare the input image for FFT
    // pad zero to input image and save it to fft_buffer_1
    ufo_buffer_pad_zero( priv->cmd_queue, &( priv->kernel_padzero ), inputs[0], priv->fft_buffer_1 );

    // FFT on the input image
#ifdef HAVE_AMD
#else

    clFFT_ExecuteInterleaved_Ufo (priv->cmd_queue, priv->fft_plan,
                                  1, clFFT_Forward,
                                  priv->fft_mem_1, priv->fft_mem_1,
                                  1, &event, NULL, priv->profiler);
#endif

    for (guint i=0; i<priv->batch_size; i++)
    {
        /*printf("\n");*/
        // copy buffer
        ufo_buffer_copy_page(inputs[1], tmp_buffer, i, priv->cmd_queue);
        ufo_buffer_pad_zero( priv->cmd_queue, &( priv->kernel_padzero ), tmp_buffer, priv->fft_buffer_2 );
        
#ifdef HAVE_AMD
#else
        clFFT_ExecuteInterleaved_Ufo (priv->cmd_queue, priv->fft_plan,
                                      1, clFFT_Forward,
                                      priv->fft_mem_2, priv->fft_mem_2,
                                      1, &event, NULL, priv->profiler);
#endif
 
        UFO_RESOURCES_CHECK_CLERR (
                   clSetKernelArg ( priv->kernel_fftmult, 0, 
                                    sizeof (cl_mem), (gpointer) &( priv->fft_mem_1 )));
        UFO_RESOURCES_CHECK_CLERR (
                   clSetKernelArg ( priv->kernel_fftmult, 1, 
                                    sizeof (cl_mem), (gpointer) &( priv->fft_mem_2 )));
        UFO_RESOURCES_CHECK_CLERR (
                   clSetKernelArg ( priv->kernel_fftmult, 2, 
                                    sizeof (cl_mem), (gpointer) &( priv->fft_mem_2 )));
        UFO_RESOURCES_CHECK_CLERR (
           clEnqueueNDRangeKernel ( priv->cmd_queue, priv->kernel_fftmult,
                                    2, NULL, global_work_size, NULL,
                                    0, NULL, &event));

        // Inverse FFT
#ifdef HAVE_AMD
#else
        clFFT_ExecuteInterleaved_Ufo (priv->cmd_queue, priv->fft_plan,
                                      1, clFFT_Inverse,
                                      //priv->batch_size, clFFT_Inverse,
                                      priv->fft_mem_2, priv->fft_mem_2,
                                      1, &event, NULL, priv->profiler);
#endif

        ufo_buffer_fftpack( priv->cmd_queue,  &priv->kernel_fftpack, priv->fft_buffer_2, tmp_buffer);

        ufo_buffer_past_page(priv->cmd_queue, tmp_buffer, output, i);
    }

    //ufo_buffer_copy(fft_buffer_2, output);
    //ufo_buffer_copy(inputs[1], output);
    //ufo_buffer_copy(inputs[0], output);
    //ufo_buffer_copy(tmp_buffer, output);

    //UFO_RESOURCES_CHECK_CLERR (clReleaseEvent (event));

    g_object_unref(tmp_buffer);
    return TRUE;
}

static void
ufo_fftconvolution_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    //UfoFftconvolutionTaskPrivate *priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_fftconvolution_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    //UfoFftconvolutionTaskPrivate *priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_fftconvolution_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_fftconvolution_task_parent_class)->finalize (object);
    UfoFftconvolutionTaskPrivate *priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE (object);

    // release kernel and destroy FFT PLAN
    // destroy fft_buffer
    g_object_unref(priv->fft_buffer_1);
    g_object_unref(priv->fft_buffer_2);
    clFFT_DestroyPlan(priv->fft_plan);        
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_fftconvolution_task_setup;
    iface->get_num_inputs = ufo_fftconvolution_task_get_num_inputs;
    iface->get_num_dimensions = ufo_fftconvolution_task_get_num_dimensions;
    iface->get_mode = ufo_fftconvolution_task_get_mode;
    iface->get_requisition = ufo_fftconvolution_task_get_requisition;
    iface->process = ufo_fftconvolution_task_process;
}

static void
ufo_fftconvolution_task_class_init (UfoFftconvolutionTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_fftconvolution_task_set_property;
    oclass->get_property = ufo_fftconvolution_task_get_property;
    oclass->finalize = ufo_fftconvolution_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoFftconvolutionTaskPrivate));
}

static void
ufo_fftconvolution_task_init(UfoFftconvolutionTask *self)
{
    self->priv = UFO_FFTCONVOLUTION_TASK_GET_PRIVATE(self);
}
