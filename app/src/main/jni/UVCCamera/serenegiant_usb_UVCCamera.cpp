/*
 * UVCCamera
 * library and sample to access to UVC web camera on non-rooted Android device
 *
 * Copyright (c) 2014-2015 saki t_saki@serenegiant.com
 *
 * File name: serenegiant_usb_UVCCamera.cpp
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 * All files in the folder are under this Apache License, Version 2.0.
 * Files in the jni/libjpeg, jni/libusb, jin/libuvc, jni/rapidjson folder may have a different license, see the respective files.
*/

#include <jni.h>
#include <android/native_window_jni.h>

#include "libUVCCamera.h"
#include "UVCCamera.h"

/**
 * set the value into the long field
 * @param env: this param should not be null
 * @param bullet_obj: this param should not be null
 * @param field_name
 * @params val
 */
static jlong setField_long(JNIEnv *env, jobject java_obj, const char *field_name, jlong val) {
#if LOCAL_DEBUG
	LOGV("setField_long:");
#endif

	jclass clazz = env->GetObjectClass(java_obj);
	jfieldID field = env->GetFieldID(clazz, field_name, "J");
	if (LIKELY(field))
		env->SetLongField(java_obj, field, val);
	else {
		LOGE("__setField_long:field '%s' not found", field_name);
	}
#ifdef ANDROID_NDK
	env->DeleteLocalRef(clazz);
#endif
	return val;
}

/**
 * @param env: this param should not be null
 * @param bullet_obj: this param should not be null
 */
static jlong __setField_long(JNIEnv *env, jobject java_obj, jclass clazz, const char *field_name, jlong val) {
#if LOCAL_DEBUG
	LOGV("__setField_long:");
#endif

	jfieldID field = env->GetFieldID(clazz, field_name, "J");
	if (LIKELY(field))
		env->SetLongField(java_obj, field, val);
	else {
		LOGE("__setField_long:field '%s' not found", field_name);
	}
	return val;
}

/**
 * @param env: this param should not be null
 * @param bullet_obj: this param should not be null
 */
jint __setField_int(JNIEnv *env, jobject java_obj, jclass clazz, const char *field_name, jint val) {
	LOGV("__setField_int:");

	jfieldID id = env->GetFieldID(clazz, field_name, "I");
	if (LIKELY(id))
		env->SetIntField(java_obj, id, val);
	else {
		LOGE("__setField_int:field '%s' not found", field_name);
		env->ExceptionClear();	// clear java.lang.NoSuchFieldError exception
	}
	return val;
}

/**
 * set the value into int field
 * @param env: this param should not be null
 * @param java_obj: this param should not be null
 * @param field_name
 * @params val
 */
jint setField_int(JNIEnv *env, jobject java_obj, const char *field_name, jint val) {
	LOGV("setField_int:");

	jclass clazz = env->GetObjectClass(java_obj);
	__setField_int(env, java_obj, clazz, field_name, val);
#ifdef ANDROID_NDK
	env->DeleteLocalRef(clazz);
#endif
	return val;
}

JNIEXPORT static ID_TYPE JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeCreate(JNIEnv *env, jobject thiz) {

	ENTER();
	UVCCamera *camera = new UVCCamera();
	setField_long(env, thiz, "mNativePtr", reinterpret_cast<ID_TYPE>(camera));
	RETURN(reinterpret_cast<ID_TYPE>(camera), ID_TYPE)
}

static void JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeDestroy(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	ENTER();
	setField_long(env, thiz, "mNativePtr", 0);
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		SAFE_DELETE(camera);
	}
	EXIT();
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeConnect(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera,
	jint vid, jint pid, jint fd,  jstring usbfs_str) {

	ENTER();
	int result = JNI_ERR;
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	const char *c_usbfs = env->GetStringUTFChars(usbfs_str, JNI_FALSE);
	if (LIKELY(camera && (fd > 0))) {
		 result =  camera->connect(vid, pid, fd, c_usbfs);
	}
	env->ReleaseStringUTFChars(usbfs_str, c_usbfs);
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeRelease(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	ENTER();
	int result = JNI_ERR;
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->release();
	}
	RETURN(result, jint);
}

static jobject JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSupportedSize(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	ENTER();
	jstring result = NULL;
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		char *c_str = camera->getSupportedSize();
		if (LIKELY(c_str)) {
			result = env->NewStringUTF(c_str);
			free(c_str);
		}
	}
	RETURN(result, jobject);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPreviewSize(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint width, jint height, jint mode, jfloat bandwidth) {

	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		return camera->setPreviewSize(width, height, mode, bandwidth);
	}


	RETURN(JNI_ERR, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeStartPreview(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		return camera->startPreview();
	}
	RETURN(JNI_ERR, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeStopPreview(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->stopPreview();
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPreviewDisplay(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jobject jSurface) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		ANativeWindow *preview_window = jSurface ? ANativeWindow_fromSurface(env, jSurface) : NULL;
		result = camera->setPreviewDisplay(preview_window);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetFrameCallback(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jobject jIFrameCallback, jint pixel_format) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		jobject frame_callback_obj = env->NewGlobalRef(jIFrameCallback);
		result = camera->setFrameCallback(env, frame_callback_obj, pixel_format);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetCaptureDisplay(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jobject jSurface) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		ANativeWindow *capture_window = jSurface ? ANativeWindow_fromSurface(env, jSurface) : NULL;
		result = camera->setCaptureDisplay(capture_window);
	}
	RETURN(result, jint);
}

//======================================================================
static jlong JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetCtrlSupports(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jlong result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		uint64_t supports;
		int r = camera->getCtrlSupports(&supports);
		if (!r)
			result = supports;
	}
	RETURN(result, jlong);
}

static jlong JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetProcSupports(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jlong result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		uint64_t supports;
		int r = camera->getProcSupports(&supports);
		if (!r)
			result = supports;
	}
	RETURN(result, jlong);
}

//======================================================================
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetExposureMode(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, int exposureMode) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setExposureMode(exposureMode);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetExposureMode(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getExposureMode();
	}
	RETURN(result, jint);
}

//======================================================================
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetAutoFocus(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jboolean autofocus) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setAutoFocus(autofocus);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetAutoFocus(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getAutoFocus();
	}
	RETURN(result, jint);
}

//======================================================================
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetAutoWhiteBlance(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jboolean autofocus) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setAutoWhiteBlance(autofocus);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetAutoWhiteBlance(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getAutoWhiteBlance();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateBrightnessLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateBrightnessLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mBrightnessMin", min);
			setField_int(env, thiz, "mBrightnessMax", max);
			setField_int(env, thiz, "mBrightnessDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetBrightness(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint brightness) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setBrightness(brightness);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetBrightness(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getBrightness();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateFocusLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateFocusLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mFocusMin", min);
			setField_int(env, thiz, "mFocusMax", max);
			setField_int(env, thiz, "mFocusDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetFocus(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint focus) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setFocus(focus);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetFocus(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getFocus();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateContrastLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateContrastLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mContrastMin", min);
			setField_int(env, thiz, "mContrastMax", max);
			setField_int(env, thiz, "mContrastDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetContrast(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint contrast) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setContrast(contrast);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetContrast(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getContrast();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateSharpnessLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateSharpnessLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mSharpnessMin", min);
			setField_int(env, thiz, "mSharpnessMax", max);
			setField_int(env, thiz, "mSharpnessDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetSharpness(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint sharpness) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setSharpness(sharpness);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSharpness(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getSharpness();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateGainLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateGainLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mGainMin", min);
			setField_int(env, thiz, "mGainMax", max);
			setField_int(env, thiz, "mGainDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetGain(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint gain) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setGain(gain);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetGain(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getGain();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateGammaLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateGammaLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mGammaMin", min);
			setField_int(env, thiz, "mGammaMax", max);
			setField_int(env, thiz, "mGammaDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetGamma(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint gamma) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setGamma(gamma);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetGamma(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getGamma();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateWhiteBlanceLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateWhiteBlanceLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mWhiteBlanceMin", min);
			setField_int(env, thiz, "mWhiteBlanceMax", max);
			setField_int(env, thiz, "mWhiteBlanceDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetWhiteBlance(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint whiteBlance) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setWhiteBlance(whiteBlance);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetWhiteBlance(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getWhiteBlance();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateSaturationLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateSaturationLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mSaturationMin", min);
			setField_int(env, thiz, "mSaturationMax", max);
			setField_int(env, thiz, "mSaturationDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetSaturation(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint saturation) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setSaturation(saturation);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSaturation(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getSaturation();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateHueLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateHueLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mHueMin", min);
			setField_int(env, thiz, "mHueMax", max);
			setField_int(env, thiz, "mHueDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetHue(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint hue) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setHue(hue);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetHue(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getHue();
	}
	RETURN(result, jint);
}

//======================================================================
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPowerlineFrequency(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint frequency) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setPowerlineFrequency(frequency);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetPowerlineFrequency(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getPowerlineFrequency();
	}
	RETURN(result, jint);
}

//======================================================================
// Java mnethod correspond to this function should not be a static mathod
static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateZoomLimit(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {
	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		int min, max, def;
		result = camera->updateZoomLimit(min, max, def);
		if (!result) {
			// Java側へ書き込む
			setField_int(env, thiz, "mZoomMin", min);
			setField_int(env, thiz, "mZoomMax", max);
			setField_int(env, thiz, "mZoomDef", def);
		}
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeSetZoom(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera, jint zoom) {

	jint result = JNI_ERR;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->setZoom(zoom);
	}
	RETURN(result, jint);
}

static jint JNICALL Java_com_dreamguard_usb_camera_UVCCamera_nativeGetZoom(JNIEnv *env, jobject thiz,
	ID_TYPE id_camera) {

	jint result = 0;
	ENTER();
	UVCCamera *camera = reinterpret_cast<UVCCamera *>(id_camera);
	if (LIKELY(camera)) {
		result = camera->getZoom();
	}
	RETURN(result, jint);
}

//**********************************************************************
//
//**********************************************************************
jint registerNativeMethods(JNIEnv* env, const char *class_name, JNINativeMethod *methods, int num_methods) {
	int result = 0;

	jclass clazz = env->FindClass(class_name);
	if (LIKELY(clazz)) {
		int result = env->RegisterNatives(clazz, methods, num_methods);
		if (UNLIKELY(result < 0)) {
			LOGE("registerNativeMethods failed(class=%s)", class_name);
		}
	} else {
		LOGE("registerNativeMethods: class'%s' not found", class_name);
	}
	return result;
}

static JNINativeMethod methods[] = {
	{ "nativeCreate",					"()J", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeCreate },
	{ "nativeDestroy",					"(J)V", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeDestroy },

	{ "nativeConnect",					"(JIIILjava/lang/String;)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeConnect },
	{ "nativeRelease",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeRelease },

	{ "nativeGetSupportedSize",			"(J)Ljava/lang/String;", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSupportedSize },
	{ "nativeSetPreviewSize",			"(JIIIF)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPreviewSize },
	{ "nativeStartPreview",				"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeStartPreview },
	{ "nativeStopPreview",				"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeStopPreview },
	{ "nativeSetPreviewDisplay",		"(JLandroid/view/Surface;)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPreviewDisplay },
	{ "nativeSetFrameCallback",			"(JLcom/dreamguard/usb/camera/IFrameCallback;I)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetFrameCallback },

	{ "nativeSetCaptureDisplay",		"(JLandroid/view/Surface;)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetCaptureDisplay },

	{ "nativeGetCtrlSupports",			"(J)J", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetCtrlSupports },
	{ "nativeGetProcSupports",			"(J)J", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetProcSupports },

	{ "nativeSetExposureMode",			"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetExposureMode },
	{ "nativeGetExposureMode",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetExposureMode },

	{ "nativeSetAutoFocus",				"(JZ)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetAutoFocus },
	{ "nativeGetAutoFocus",				"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetAutoFocus },

	{ "nativeUpdateFocusLimit",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateFocusLimit },
	{ "nativeSetFocus",					"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetFocus },
	{ "nativeGetFocus",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetFocus },

	{ "nativeSetAutoWhiteBlance",		"(JZ)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetAutoWhiteBlance },
	{ "nativeGetAutoWhiteBlance",		"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetAutoWhiteBlance },

	{ "nativeUpdateWhiteBlanceLimit",	"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateWhiteBlanceLimit },
	{ "nativeSetWhiteBlance",			"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetWhiteBlance },
	{ "nativeGetWhiteBlance",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetWhiteBlance },

	{ "nativeUpdateBrightnessLimit",	"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateBrightnessLimit },
	{ "nativeSetBrightness",			"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetBrightness },
	{ "nativeGetBrightness",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetBrightness },

	{ "nativeUpdateContrastLimit",		"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateContrastLimit },
	{ "nativeSetContrast",				"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetContrast },
	{ "nativeGetContrast",				"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetContrast },

	{ "nativeUpdateSharpnessLimit",		"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateSharpnessLimit },
	{ "nativeSetSharpness",				"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetSharpness },
	{ "nativeGetSharpness",				"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSharpness },

	{ "nativeUpdateGainLimit",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateGainLimit },
	{ "nativeSetGain",					"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetGain },
	{ "nativeGetGain",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetGain },

	{ "nativeUpdateGammaLimit",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateGammaLimit },
	{ "nativeSetGamma",					"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetGamma },
	{ "nativeGetGamma",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetGamma },

	{ "nativeUpdateSaturationLimit",	"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateSaturationLimit },
	{ "nativeSetSaturation",			"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetSaturation },
	{ "nativeGetSaturation",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetSaturation },

	{ "nativeUpdateHueLimit",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateHueLimit },
	{ "nativeSetHue",					"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetHue },
	{ "nativeGetHue",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetHue },

	{ "nativeSetPowerlineFrequency",	"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetPowerlineFrequency },
	{ "nativeGetPowerlineFrequency",	"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetPowerlineFrequency },

	{ "nativeUpdateZoomLimit",			"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeUpdateZoomLimit },
	{ "nativeSetZoom",					"(JI)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeSetZoom },
	{ "nativeGetZoom",					"(J)I", (void *) Java_com_dreamguard_usb_camera_UVCCamera_nativeGetZoom },
};

int register_uvccamera(JNIEnv *env) {
	LOGV("register_uvccamera:");
	if (registerNativeMethods(env,
		"com/dreamguard/usb/camera/UVCCamera",
		methods, NUM_ARRAY_ELEMENTS(methods)) < 0) {
		return -1;
	}
    return 0;
}
