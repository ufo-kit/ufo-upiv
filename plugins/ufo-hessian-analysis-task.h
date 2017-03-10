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

#ifndef __UFO_HESSIAN_ANALYSIS_TASK_H
#define __UFO_HESSIAN_ANALYSIS_TASK_H

#include <ufo/ufo.h>

G_BEGIN_DECLS

#define UFO_TYPE_HESSIAN_ANALYSIS_TASK             (ufo_hessian_analysis_task_get_type())
#define UFO_HESSIAN_ANALYSIS_TASK(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj), UFO_TYPE_HESSIAN_ANALYSIS_TASK, UfoHessianAnalysisTask))
#define UFO_IS_HESSIAN_ANALYSIS_TASK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj), UFO_TYPE_HESSIAN_ANALYSIS_TASK))
#define UFO_HESSIAN_ANALYSIS_TASK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass), UFO_TYPE_HESSIAN_ANALYSIS_TASK, UfoHessianAnalysisTaskClass))
#define UFO_IS_HESSIAN_ANALYSIS_TASK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass), UFO_TYPE_HESSIAN_ANALYSIS_TASK))
#define UFO_HESSIAN_ANALYSIS_TASK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS((obj), UFO_TYPE_HESSIAN_ANALYSIS_TASK, UfoHessianAnalysisTaskClass))

typedef struct _UfoHessianAnalysisTask           UfoHessianAnalysisTask;
typedef struct _UfoHessianAnalysisTaskClass      UfoHessianAnalysisTaskClass;
typedef struct _UfoHessianAnalysisTaskPrivate    UfoHessianAnalysisTaskPrivate;

/**
 * UfoHessianAnalysisTask:
 *
 * [ADD DESCRIPTION HERE]. The contents of the #UfoHessianAnalysisTask structure
 * are private and should only be accessed via the provided API.
 */
struct _UfoHessianAnalysisTask {
    /*< private >*/
    UfoTaskNode parent_instance;

    UfoHessianAnalysisTaskPrivate *priv;
};

/**
 * UfoHessianAnalysisTaskClass:
 *
 * #UfoHessianAnalysisTask class
 */
struct _UfoHessianAnalysisTaskClass {
    /*< private >*/
    UfoTaskNodeClass parent_class;
};

UfoNode  *ufo_hessian_analysis_task_new       (void);
GType     ufo_hessian_analysis_task_get_type  (void);

G_END_DECLS

#endif