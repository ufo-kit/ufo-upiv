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

#ifndef __UFO_BLOB_TEST_TASK_H
#define __UFO_BLOB_TEST_TASK_H

#include <ufo/ufo.h>

G_BEGIN_DECLS

#define UFO_TYPE_BLOB_TEST_TASK             (ufo_blob_test_task_get_type())
#define UFO_BLOB_TEST_TASK(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj), UFO_TYPE_BLOB_TEST_TASK, UfoBlobTestTask))
#define UFO_IS_BLOB_TEST_TASK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj), UFO_TYPE_BLOB_TEST_TASK))
#define UFO_BLOB_TEST_TASK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass), UFO_TYPE_BLOB_TEST_TASK, UfoBlobTestTaskClass))
#define UFO_IS_BLOB_TEST_TASK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass), UFO_TYPE_BLOB_TEST_TASK))
#define UFO_BLOB_TEST_TASK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS((obj), UFO_TYPE_BLOB_TEST_TASK, UfoBlobTestTaskClass))

typedef struct _UfoBlobTestTask           UfoBlobTestTask;
typedef struct _UfoBlobTestTaskClass      UfoBlobTestTaskClass;
typedef struct _UfoBlobTestTaskPrivate    UfoBlobTestTaskPrivate;

struct _UfoBlobTestTask {
    UfoTaskNode parent_instance;

    UfoBlobTestTaskPrivate *priv;
};

struct _UfoBlobTestTaskClass {
    UfoTaskNodeClass parent_class;
};

UfoNode  *ufo_blob_test_task_new       (void);
GType     ufo_blob_test_task_get_type  (void);

G_END_DECLS

#endif