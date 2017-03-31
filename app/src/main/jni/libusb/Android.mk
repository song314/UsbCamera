# Android build config for libusb, examples and tests
# Copyright © 2012-2013 RealVNC Ltd. <toby.gray@realvnc.com>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

LOCAL_PATH:= $(call my-dir)

include $(CLEAR_VARS)

LOCAL_SRC_FILES := \
	core.c \
	descriptor.c \
	hotplug.c \
	io.c \
	sync.c \
	strerror.c \
	os/android_usbfs.c \
	os/poll_posix.c \
	os/threads_posix.c \
	os/android_netlink.c

LOCAL_C_INCLUDES += \
	$(LOCAL_PATH)/ \
	$(LOCAL_PATH)/os \

LOCAL_EXPORT_C_INCLUDES := \
	$(LOCAL_PATH)/ \

# add some flags
LOCAL_CFLAGS := $(LOCAL_C_INCLUDES:%=-I%)
LOCAL_CFLAGS += -DANDROID_NDK
LOCAL_CFLAGS += -DLOG_NDEBUG
LOCAL_CFLAGS += -DACCESS_RAW_DESCRIPTORS
LOCAL_CFLAGS += -O3 -fstrict-aliasing -fprefetch-loop-arrays
LOCAL_EXPORT_LDLIBS += -llog
LOCAL_ARM_MODE := arm

LOCAL_MODULE := libusb100_static
include $(BUILD_STATIC_LIBRARY)

######################################################################
# libusb.so
######################################################################
include $(CLEAR_VARS)
LOCAL_MODULE_TAGS := optional
LOCAL_EXPORT_LDLIBS += -llog

LOCAL_WHOLE_STATIC_LIBRARIES = libusb100_static

LOCAL_MODULE := libusb100
include $(BUILD_SHARED_LIBRARY)

