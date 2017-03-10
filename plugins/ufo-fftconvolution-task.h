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

#ifndef __UFO_FFTCONVOLUTION_TASK_H
#define __UFO_FFTCONVOLUTION_TASK_H

#include <ufo/ufo.h>

G_BEGIN_DECLS

#define UFO_TYPE_FFTCONVOLUTION_TASK             (ufo_fftconvolution_task_get_type())
#define UFO_FFTCONVOLUTION_TASK(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj), UFO_TYPE_FFTCONVOLUTION_TASK, UfoFftconvolutionTask))
#define UFO_IS_FFTCONVOLUTION_TASK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj), UFO_TYPE_FFTCONVOLUTION_TASK))
#define UFO_FFTCONVOLUTION_TASK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass), UFO_TYPE_FFTCONVOLUTION_TASK, UfoFftconvolutionTaskClass))
#define UFO_IS_FFTCONVOLUTION_TASK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass), UFO_TYPE_FFTCONVOLUTION_TASK))
#define UFO_FFTCONVOLUTION_TASK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS((obj), UFO_TYPE_FFTCONVOLUTION_TASK, UfoFftconvolutionTaskClass))

typedef struct _UfoFftconvolutionTask           UfoFftconvolutionTask;
typedef struct _UfoFftconvolutionTaskClass      UfoFftconvolutionTaskClass;
typedef struct _UfoFftconvolutionTaskPrivate    UfoFftconvolutionTaskPrivate;

/**
 * UfoFftconvolutionTask:
 *
 * [ADD DESCRIPTION HERE]. The contents of the #UfoFftconvolutionTask structure
 * are private and should only be accessed via the provided API.
 */
struct _UfoFftconvolutionTask {
    /*< private >*/
    UfoTaskNode parent_instance;

    UfoFftconvolutionTaskPrivate *priv;
};

/**
 * UfoFftconvolutionTaskClass:
 *
 * #UfoFftconvolutionTask class
 */
struct _UfoFftconvolutionTaskClass {
    /*< private >*/
    UfoTaskNodeClass parent_class;
};

UfoNode  *ufo_fftconvolution_task_new       (void);
GType     ufo_fftconvolution_task_get_type  (void);

void ufo_buffer_copy_page(UfoBuffer *src, UfoBuffer *dst, const guint n, cl_command_queue cmd_queue);
void ufo_buffer_past_page(cl_command_queue cmd_queue, UfoBuffer *src, UfoBuffer *dst, const guint n);
void ufo_buffer_pad_zero(cl_command_queue cmd_queue, cl_kernel *kernel, UfoBuffer *src_buf, UfoBuffer *dst_buf);
cl_int ufo_buffer_fill_zero(cl_command_queue cmd_queue, cl_kernel *kernel, UfoBuffer *buf);
void ufo_buffer_fftpack (cl_command_queue cmd_queue, cl_kernel *kernel, UfoBuffer *src_buf, UfoBuffer *dst_buf);

G_END_DECLS

#endif
