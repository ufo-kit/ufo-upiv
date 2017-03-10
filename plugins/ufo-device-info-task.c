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

#include "ufo-device-info-task.h"

struct _UfoDeviceInfoTaskPrivate {
    gboolean foo;
};

static void ufo_task_interface_init (UfoTaskIface *iface);

G_DEFINE_TYPE_WITH_CODE (UfoDeviceInfoTask, ufo_device_info_task, UFO_TYPE_TASK_NODE,
                         G_IMPLEMENT_INTERFACE (UFO_TYPE_TASK,
                                                ufo_task_interface_init))

#define UFO_DEVICE_INFO_TASK_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE((obj), UFO_TYPE_DEVICE_INFO_TASK, UfoDeviceInfoTaskPrivate))

enum {
    PROP_0,
    PROP_TEST,
    N_PROPERTIES
};

static GParamSpec *properties[N_PROPERTIES] = { NULL, };

UfoNode *
ufo_device_info_task_new (void)
{
    return UFO_NODE (g_object_new (UFO_TYPE_DEVICE_INFO_TASK, NULL));
}

static void
ufo_device_info_task_setup (UfoTask *task,
                       UfoResources *resources,
                       GError **error)
{
}

static void
ufo_device_info_task_get_requisition (UfoTask *task,
                                 UfoBuffer **inputs,
                                 UfoRequisition *requisition)
{
    requisition->n_dims = 1;
    requisition->dims[0] = 1;
}

static guint
ufo_device_info_task_get_num_inputs (UfoTask *task)
{
    return 1;
}

static guint
ufo_device_info_task_get_num_dimensions (UfoTask *task,
                                             guint input)
{
    return 2;
}

static UfoTaskMode
ufo_device_info_task_get_mode (UfoTask *task)
{
    return UFO_TASK_MODE_PROCESSOR | UFO_TASK_MODE_GPU;
}

static gboolean
ufo_device_info_task_process (UfoTask *task,
                         UfoBuffer **inputs,
                         UfoBuffer *output,
                         UfoRequisition *requisition)
{
    cl_platform_id platform_id;
    cl_device_id device_id = 0;
    cl_uint work_item_dim;
    size_t work_item_sizes[3];
    size_t work_group_size;
    cl_ulong local_mem_size;
    cl_uint max_compute_units;

    clGetPlatformIDs(1, &platform_id, NULL);
    clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);

    clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &work_item_dim, NULL);
    clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(work_item_sizes), work_item_sizes, NULL);
    clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL);
    clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);

    g_message ("device id: %u", device_id);
    g_message ("work_item_dim: %d", work_item_dim);
    g_message ("work_item_sizes: %d %d %d", work_item_sizes[0], work_item_sizes[1], work_item_sizes[2]);
    g_message ("work_group_size: %d", work_group_size);
    g_message ("local_mem_size: %lu", local_mem_size);
    g_message ("maximum compute units: %u", max_compute_units);

    return TRUE;
}

static void
ufo_device_info_task_set_property (GObject *object,
                              guint property_id,
                              const GValue *value,
                              GParamSpec *pspec)
{
    UfoDeviceInfoTaskPrivate *priv = UFO_DEVICE_INFO_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_device_info_task_get_property (GObject *object,
                              guint property_id,
                              GValue *value,
                              GParamSpec *pspec)
{
    UfoDeviceInfoTaskPrivate *priv = UFO_DEVICE_INFO_TASK_GET_PRIVATE (object);

    switch (property_id) {
        case PROP_TEST:
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
ufo_device_info_task_finalize (GObject *object)
{
    G_OBJECT_CLASS (ufo_device_info_task_parent_class)->finalize (object);
}

static void
ufo_task_interface_init (UfoTaskIface *iface)
{
    iface->setup = ufo_device_info_task_setup;
    iface->get_num_inputs = ufo_device_info_task_get_num_inputs;
    iface->get_num_dimensions = ufo_device_info_task_get_num_dimensions;
    iface->get_mode = ufo_device_info_task_get_mode;
    iface->get_requisition = ufo_device_info_task_get_requisition;
    iface->process = ufo_device_info_task_process;
}

static void
ufo_device_info_task_class_init (UfoDeviceInfoTaskClass *klass)
{
    GObjectClass *oclass = G_OBJECT_CLASS (klass);

    oclass->set_property = ufo_device_info_task_set_property;
    oclass->get_property = ufo_device_info_task_get_property;
    oclass->finalize = ufo_device_info_task_finalize;

    properties[PROP_TEST] =
        g_param_spec_string ("test",
            "Test property nick",
            "Test property description blurb",
            "",
            G_PARAM_READWRITE);

    for (guint i = PROP_0 + 1; i < N_PROPERTIES; i++)
        g_object_class_install_property (oclass, i, properties[i]);

    g_type_class_add_private (oclass, sizeof(UfoDeviceInfoTaskPrivate));
}

static void
ufo_device_info_task_init(UfoDeviceInfoTask *self)
{
    self->priv = UFO_DEVICE_INFO_TASK_GET_PRIVATE(self);
}
