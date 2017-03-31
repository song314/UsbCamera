#include $(call all-subdir-makefiles)
PROJ_PATH	:= $(call my-dir)
include $(CLEAR_VARS)
include $(PROJ_PATH)/UVCCamera/Android.mk
include $(PROJ_PATH)/libjpeg/Android.mk
include $(PROJ_PATH)/libusb/Android.mk
include $(PROJ_PATH)/libuvc/Android.mk