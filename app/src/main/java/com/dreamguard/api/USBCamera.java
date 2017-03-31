package com.dreamguard.api;

import android.content.Context;
import android.graphics.Bitmap;
import android.graphics.SurfaceTexture;
import android.hardware.usb.UsbDevice;
import android.media.AudioManager;
import android.media.SoundPool;
import android.os.Environment;
import android.util.Log;
import android.view.Surface;
import android.view.View;
import android.widget.RelativeLayout;

import com.dreamguard.usb.camera.CameraHandler;
import com.dreamguard.usb.camera.UVCCamera;
import com.dreamguard.usb.detect.DeviceFilter;
import com.dreamguard.usb.detect.USBMonitor;
import com.dreamguard.usb.detect.USBStatus;
import com.dreamguard.widget.UVCCameraTextureView;
import com.quinn.usbcameratest.R;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.lang.reflect.Field;
import java.util.List;

/**
 * Created by hailin on 2016/12/10.
 */

public class USBCamera implements IUSBCameraParamsControll {

    private final static String TAG = "USBCamera";
    private static final String DIR_NAME = "720vv_O2";
    /**
     * for accessing USB
     */
    private USBMonitor mUSBMonitor;
    /**
     * Handler to execute camera releated methods sequentially on private thread
     */
    private CameraHandler mHandler;

    private Context context;

    private RelativeLayout tip;

    private Surface mSurfaceTexture;

    private USBStatus usbStatus = USBStatus.DETACHED;

    private final Object mSync = new Object();
    private SoundPool mSoundPool;
    private int mSoundId;

    public void init(Context context, RelativeLayout tip) {
        Log.v(TAG, "init :");
        this.context = context;
        this.tip = tip;
        mUSBMonitor = new USBMonitor(context, mOnDeviceConnectListener);
        mHandler = CameraHandler.createHandler(context);
        mUSBMonitor.register();
        loadSutterSound(context);
    }

    private void loadSutterSound(final Context context) {
        int streamType;
        try {
            final Class<?> audioSystemClass = Class.forName("android.media.AudioSystem");
            final Field sseField = audioSystemClass.getDeclaredField("STREAM_SYSTEM_ENFORCED");
            streamType = sseField.getInt(null);
        } catch (final Exception e) {
            streamType = AudioManager.STREAM_SYSTEM;
        }
        if (mSoundPool != null) {
            try {
                mSoundPool.release();
            } catch (final Exception e) {
            }
            mSoundPool = null;
        }
        mSoundPool = new SoundPool(2, streamType, 0);
//        mSoundId = mSoundPool.load(context, R.raw.usb_in, 1);
        mSoundId = mSoundPool.load(context, R.raw.camera_click, 1);
    }

    public void destroy() {
        Log.v(TAG, "destroy :");
        mUSBMonitor.unregister();
        mUSBMonitor.destroy();
        mHandler = null;
    }

    public void setPreviewSize(int width, int height) {
        CameraHandler.PREVIEW_WIDTH = width;
        CameraHandler.PREVIEW_HEIGHT = height;
        CameraHandler.CAPTURE_WIDTH = width;
        CameraHandler.CAPTURE_HEIGHT = height;
        CameraHandler.RECORD_WIDTH = width;
        CameraHandler.RECORD_HEIGHT = height;


    }

    public void setCameraType(CameraType cameraType) {
        if (cameraType == CameraType.C3D_NORMAL) {
            CameraHandler.is3D = false;
        } else {
            CameraHandler.is3D = true;
        }
    }

    public void setPreviewTexture(Surface surfaceTexture) {
        mSurfaceTexture = surfaceTexture;
    }

    public void startPreview() {

    }

    public boolean open(int id) {
        Log.v(TAG, "open :");
        final List<DeviceFilter> filter = DeviceFilter.getDeviceFilters(context, R.xml.device_filter);
        List<UsbDevice> deviceList = mUSBMonitor.getDeviceList(filter.get(0));
        UsbDevice device = null;
        if (deviceList.size() > id) {
            device = deviceList.get(id);
        }
        if (device != null) {
            Log.v(TAG, "open :" + device.toString());
            mUSBMonitor.requestPermission(device);
            return true;
        } else {
            Log.v(TAG, "open null:");
            return false;
        }

    }

    public void close() {
        mHandler.closeCamera();
    }

    public boolean isCameraOpened() {
        return mHandler.isCameraOpened();
    }

    public UVCCamera getUVCCamera() {
        return mHandler.getUVCCamera();
    }

    public String getSupportedSize() {
        return mHandler.getSupportedSize();
    }

    public void captureStill(UVCCameraTextureView mCameraView) {
        mSoundPool.play(mSoundId, 0.2f, 0.2f, 0, 0, 1.0f);
        File outputFile = null;
        BufferedOutputStream os = null;

        try {
            File tempFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DCIM),DIR_NAME);
            if (!tempFile.exists()) {
                tempFile.mkdirs();
            }
            outputFile = new File(tempFile, System.currentTimeMillis() + ".jpg");
            os = new BufferedOutputStream(new FileOutputStream(outputFile));
            Bitmap bitmap = mCameraView.getBitmap();
            bitmap.compress(Bitmap.CompressFormat.JPEG, 100, os);
            os.flush();
            os.close();
            if (bitmap!=null) {
                bitmap.recycle();
            }
        } catch (Exception e) {
            Log.e("MianActivity", "onFrame Capture still Error Exception");
        }
    }

    public void startRecording() {
        mHandler.startRecording();
    }

    public void stopRecording() {
        mHandler.stopRecording();
    }

    public boolean isRecording() {
        return mHandler.isRecording();
    }

    private final USBMonitor.OnDeviceConnectListener mOnDeviceConnectListener = new USBMonitor.OnDeviceConnectListener() {
        @Override
        public void onAttach(final UsbDevice device) {
            Log.v(TAG, "onAttach:");
//            tip.setVisibility(View.GONE);
            usbStatus = USBStatus.ATTACHED;
        }

        @Override
        public void onConnect(final UsbDevice device, final USBMonitor.UsbControlBlock ctrlBlock, final boolean createNew) {
            Log.v(TAG, "onConnect:");
            usbStatus = USBStatus.CONNECTED;
            mHandler.openCamera(ctrlBlock);
            mHandler.startPreview(mSurfaceTexture);
        }

        @Override
        public void onDisconnect(final UsbDevice device, final USBMonitor.UsbControlBlock ctrlBlock) {
            Log.v(TAG, "onDisconnect:");
            usbStatus = USBStatus.DISCONNECTED;
            if (mHandler != null) {
                mHandler.closeCamera();
            }
        }

        @Override
        public void onDetach(final UsbDevice device) {
            Log.v(TAG, "onDetach:");
//            tip.setVisibility(View.VISIBLE);
            usbStatus = USBStatus.DETACHED;
        }

        @Override
        public void onCancel() {
        }
    };

    @Override
    public void setExposureMode(int exposureMode) {
        mHandler.setExposureMode(exposureMode);
    }

    @Override
    public int getExposureMode() {
        return mHandler.getExposureMode();
    }

    @Override
    public void setAutoFocus(boolean autofocus) {
        mHandler.setAutoFocus(autofocus);
    }

    @Override
    public boolean getAutoFocus() {
        return mHandler.getAutoFocus();
    }

    @Override
    public void resetFocus() {
        mHandler.resetFocus();
    }

    @Override
    public void setFocus(int focus) {
        mHandler.setFocus(focus);
    }

    @Override
    public int getFocus() {
        return mHandler.getFocus();
    }

    @Override
    public void setAutoWhiteBlance(boolean autoWhiteBlance) {
        mHandler.setAutoWhiteBlance(autoWhiteBlance);
    }

    @Override
    public boolean getAutoWhiteBlance() {
        return mHandler.getAutoWhiteBlance();
    }

    @Override
    public void resetWhiteBlance() {
        mHandler.resetWhiteBlance();
    }

    @Override
    public void setWhiteBlance(int whiteBlance) {
        mHandler.setWhiteBlance(whiteBlance);
    }

    @Override
    public int getWhiteBlance() {
        return mHandler.getWhiteBlance();
    }

    @Override
    public void resetBrightness() {
        mHandler.resetBrightness();
    }

    @Override
    public void setBrightness(int brightness) {
        mHandler.setBrightness(brightness);
    }

    @Override
    public int getBrightness() {
        return mHandler.getBrightness();
    }

    @Override
    public void resetContrast() {
        mHandler.resetContrast();
    }

    @Override
    public void setContrast(int contrast) {
        mHandler.setContrast(contrast);
    }

    @Override
    public int getContrast() {
        return mHandler.getContrast();
    }

    @Override
    public void resetSharpness() {
        mHandler.resetSharpness();
    }

    @Override
    public void setSharpness(int sharpness) {
        mHandler.setSharpness(sharpness);
    }

    @Override
    public int getSharpness() {
        return mHandler.getSharpness();
    }

    @Override
    public void resetGain() {
        mHandler.resetGain();
    }

    @Override
    public void setGain(int gain) {
        mHandler.setGain(gain);
    }

    @Override
    public int getGain() {
        return mHandler.getGain();
    }

    @Override
    public void resetGamma() {
        mHandler.resetGamma();
    }

    @Override
    public void setGamma(int gamma) {
        mHandler.setGamma(gamma);
    }

    @Override
    public int getGamma() {
        return mHandler.getGamma();
    }

    @Override
    public void resetSaturation() {
        mHandler.resetSaturation();
    }

    @Override
    public void setSaturation(int saturation) {
        mHandler.setSaturation(saturation);
    }

    @Override
    public int getSaturation() {
        return mHandler.getSaturation();
    }

    @Override
    public void resetHue() {
        mHandler.resetHue();
    }

    @Override
    public void setHue(int hue) {
        mHandler.setHue(hue);
    }

    @Override
    public int getHue() {
        return mHandler.getHue();
    }

    @Override
    public void setPowerlineFrequency(int frequency) {
        mHandler.setPowerlineFrequency(frequency);
    }

    @Override
    public int getPowerlineFrequency() {
        return mHandler.getPowerlineFrequency();
    }

    @Override
    public void resetZoom() {
        mHandler.resetZoom();
    }

    @Override
    public void setZoom(int zoom) {
        mHandler.setZoom(zoom);
    }

    @Override
    public int getZoom() {
        return mHandler.getZoom();
    }
}
