package com.dreamguard.usb.camera;

/**
 * Created by hailin.dai on 12/2/16.
 * email:hailin.dai@wz-tech.com
 */

import android.content.Context;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.ImageFormat;
import android.graphics.Rect;
import android.graphics.YuvImage;
import android.media.AudioManager;
import android.media.MediaScannerConnection;
import android.media.SoundPool;
import android.os.Environment;
import android.os.Handler;
import android.os.Looper;
import android.os.Message;
import android.util.Log;
import android.view.Surface;

import com.dreamguard.api.IUSBCameraParamsControll;
import com.dreamguard.encoder.MediaEncoder;
import com.dreamguard.encoder.MediaMuxerWrapper;
import com.dreamguard.encoder.MediaVideoEncoder;
import com.dreamguard.usb.detect.USBMonitor;
import com.dreamguard.util.ImageProc;
import com.orhanobut.logger.Logger;
import com.quinn.usbcameratest.BuildConfig;
import com.quinn.usbcameratest.R;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.ref.WeakReference;
import java.lang.reflect.Field;
import java.nio.ByteBuffer;

/**
 * Handler class to execute camera releated methods sequentially on private thread
 */
public class CameraHandler extends Handler implements IUSBCameraParamsControll {

    private static final boolean DEBUG = true;
    private static final String TAG = "CameraHandler";
    private static final String DIR_NAME = "720vv_O2";

    /**
     * preview resolution(width)
     * if your camera does not support specific resolution and mode,
     * {@link UVCCamera#setPreviewSize(int, int, int)} throw exception
     */
    public static int PREVIEW_WIDTH = 640;
    /**
     * preview resolution(height)
     * if your camera does not support specific resolution and mode,
     * {@link UVCCamera#setPreviewSize(int, int, int)} throw exception
     */
    public static int PREVIEW_HEIGHT = 480;
    /**
     * preview mode
     * if your camera does not support specific resolution and mode,
     * {@link UVCCamera#setPreviewSize(int, int, int)} throw exception
     * 0:YUYV, other:MJPEG
     */

    public static int RECORD_WIDTH = 640;

    public static int RECORD_HEIGHT = 480;

    public static int CAPTURE_WIDTH = 640;

    public static int CAPTURE_HEIGHT = 480;

    public static boolean is3D = false;

    private static final int PREVIEW_MODE = 1;


    private static final int MSG_OPEN = 0;
    private static final int MSG_CLOSE = 1;
    private static final int MSG_PREVIEW_START = 2;
    private static final int MSG_PREVIEW_STOP = 3;
    private static final int MSG_CAPTURE_STILL = 4;
    private static final int MSG_CAPTURE_START = 5;
    private static final int MSG_CAPTURE_STOP = 6;
    private static final int MSG_MEDIA_UPDATE = 7;
    private static final int MSG_RELEASE = 9;

    private final WeakReference<CameraThread> mWeakThread;

    public static final CameraHandler createHandler(final Context parent) {
        final CameraThread thread = new CameraThread(parent);
        thread.start();
        return thread.getHandler();
    }

    private CameraHandler(final CameraThread thread) {
        mWeakThread = new WeakReference<CameraThread>(thread);
    }

    public boolean isCameraOpened() {
        final CameraThread thread = mWeakThread.get();
        return thread != null ? thread.isCameraOpened() : false;
    }

    public UVCCamera getUVCCamera() {
        final CameraThread thread = mWeakThread.get();
        return thread != null ? thread.getUVCCamera() : null;
    }

    public String getSupportedSize() {
        final CameraThread thread = mWeakThread.get();
        return thread != null ? thread.getSupportedSize() : null;
    }

    public boolean isRecording() {
        final CameraThread thread = mWeakThread.get();
        return thread != null ? thread.isRecording() : false;
    }

    public void openCamera(final USBMonitor.UsbControlBlock ctrlBlock) {
        sendMessage(obtainMessage(MSG_OPEN, ctrlBlock));
    }

    public void closeCamera() {
        stopPreview();
        sendEmptyMessage(MSG_CLOSE);
    }

    public void startPreview(final Surface sureface) {
        if (sureface != null)
            sendMessage(obtainMessage(MSG_PREVIEW_START, sureface));
    }

    public void stopPreview() {
        stopRecording();
        final CameraThread thread = mWeakThread.get();
        if (thread == null) return;
        synchronized (thread.mSync) {
            sendEmptyMessage(MSG_PREVIEW_STOP);
            // wait for actually preview stopped to avoid releasing Surface/SurfaceTexture
            // while preview is still running.
            // therefore this method will take a time to execute
            try {
                thread.mSync.wait();
            } catch (final InterruptedException e) {
            }
        }
    }

    public void captureStill() {
        sendEmptyMessage(MSG_CAPTURE_STILL);
    }

    public void startRecording() {
        sendEmptyMessage(MSG_CAPTURE_START);
    }

    public void stopRecording() {
        sendEmptyMessage(MSG_CAPTURE_STOP);
    }

/*		public void release() {
            sendEmptyMessage(MSG_RELEASE);
		} */

    @Override
    public void handleMessage(final Message msg) {
        final CameraThread thread = mWeakThread.get();
        if (thread == null) return;
        switch (msg.what) {
            case MSG_OPEN:
                thread.handleOpen((USBMonitor.UsbControlBlock) msg.obj);
                break;
            case MSG_CLOSE:
                thread.handleClose();
                break;
            case MSG_PREVIEW_START:
                thread.handleStartPreview((Surface) msg.obj);
                break;
            case MSG_PREVIEW_STOP:
                thread.handleStopPreview();
                break;
            case MSG_CAPTURE_STILL:
                thread.handleCaptureStill();
                break;
            case MSG_CAPTURE_START:
                thread.handleStartRecording();
                break;
            case MSG_CAPTURE_STOP:
                thread.handleStopRecording();
                break;
            case MSG_MEDIA_UPDATE:
                thread.handleUpdateMedia((String) msg.obj);
                break;
            case MSG_RELEASE:
                thread.handleRelease();
                break;
            case MSG_EXPOSUREMODE:
                thread.setExposureMode((Integer) msg.obj);
                break;
            case MSG_AUTOFOCUS:
                thread.setAutoFocus((Boolean) msg.obj);
                break;
            case MSG_FOCUS:
                thread.setFocus((Integer) msg.obj);
                break;
            case MSG_RESET_FOCUS:
                thread.resetFocus();
                break;
            case MSG_AUTOWHITEBLANCE:
                thread.setAutoWhiteBlance((Boolean) msg.obj);
                break;
            case MSG_WHITEBLANCE:
                thread.setWhiteBlance((Integer) msg.obj);
                break;
            case MSG_RESET_WHITEBLANCE:
                thread.resetWhiteBlance();
                break;
            case MSG_BRIGHTNESS:
                thread.setBrightness((Integer) msg.obj);
                break;
            case MSG_RESET_BRIGHTNESS:
                thread.resetBrightness();
                break;
            case MSG_CONTRAST:
                thread.setContrast((Integer) msg.obj);
                break;
            case MSG_RESET_CONTRAST:
                thread.resetContrast();
                break;
            case MSG_SHARPNESS:
                thread.setSharpness((Integer) msg.obj);
                break;
            case MSG_RESET_SHARPNESS:
                thread.resetSharpness();
                break;
            case MSG_GAIN:
                thread.setGain((Integer) msg.obj);
                break;
            case MSG_RESET_GAIN:
                thread.resetGain();
                break;
            case MSG_GAMMA:
                thread.setGamma((Integer) msg.obj);
                break;
            case MSG_RESET_GAMMA:
                thread.resetGamma();
                break;
            case MSG_SATURATION:
                thread.setSaturation((Integer) msg.obj);
                break;
            case MSG_RESET_SATURATION:
                thread.resetSaturation();
                break;
            case MSG_HUE:
                thread.setHue((Integer) msg.obj);
                break;
            case MSG_RESET_HUE:
                thread.resetHue();
                break;
            case MSG_POWERLINEFREQUENCY:
                thread.setPowerlineFrequency((Integer) msg.obj);
                break;
            case MSG_ZOOM:
                thread.setZoom((Integer) msg.obj);
                break;
            case MSG_RESET_ZOOM:
                thread.resetZoom();
                break;
            default:
                throw new RuntimeException("unsupported message:what=" + msg.what);
        }
    }

    @Override
    public void setExposureMode(int exposureMode) {
        sendMessage(obtainMessage(MSG_EXPOSUREMODE, exposureMode));
    }

    @Override
    public int getExposureMode() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getExposureMode();
    }

    @Override
    public void setAutoFocus(boolean autofocus) {
        sendMessage(obtainMessage(MSG_AUTOFOCUS, autofocus));
    }

    @Override
    public boolean getAutoFocus() {
        final CameraThread thread = mWeakThread.get();
        return thread != null && thread.getAutoFocus();
    }

    @Override
    public void resetFocus() {
        sendEmptyMessage(MSG_RESET_FOCUS);
    }

    @Override
    public void setFocus(int focus) {
        sendMessage(obtainMessage(MSG_FOCUS, focus));
    }

    @Override
    public int getFocus() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getFocus();
    }

    @Override
    public void setAutoWhiteBlance(boolean autoWhiteBlance) {
        sendMessage(obtainMessage(MSG_AUTOWHITEBLANCE, autoWhiteBlance));
    }

    @Override
    public boolean getAutoWhiteBlance() {
        final CameraThread thread = mWeakThread.get();
        return thread != null && thread.getAutoWhiteBlance();
    }

    @Override
    public void resetWhiteBlance() {
        sendEmptyMessage(MSG_RESET_WHITEBLANCE);
    }

    @Override
    public void setWhiteBlance(int whiteBlance) {
        sendMessage(obtainMessage(MSG_WHITEBLANCE, whiteBlance));
    }

    @Override
    public int getWhiteBlance() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getWhiteBlance();
    }

    @Override
    public void resetBrightness() {
        sendEmptyMessage(MSG_RESET_BRIGHTNESS);
    }

    @Override
    public void setBrightness(int brightness) {
        sendMessage(obtainMessage(MSG_BRIGHTNESS, brightness));
    }

    @Override
    public int getBrightness() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getBrightness();
    }

    @Override
    public void resetContrast() {
        sendEmptyMessage(MSG_RESET_CONTRAST);
    }

    @Override
    public void setContrast(int contrast) {
        sendMessage(obtainMessage(MSG_CONTRAST, contrast));
    }

    @Override
    public int getContrast() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getContrast();
    }

    @Override
    public void resetSharpness() {
        sendEmptyMessage(MSG_RESET_SHARPNESS);
    }

    @Override
    public void setSharpness(int sharpness) {
        sendMessage(obtainMessage(MSG_SHARPNESS, sharpness));
    }

    @Override
    public int getSharpness() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getSharpness();
    }

    @Override
    public void resetGain() {
        sendEmptyMessage(MSG_RESET_GAIN);
    }

    @Override
    public void setGain(int gain) {
        sendMessage(obtainMessage(MSG_GAIN, gain));
    }

    @Override
    public int getGain() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getGain();
    }

    @Override
    public void resetGamma() {
        sendEmptyMessage(MSG_RESET_GAMMA);
    }

    @Override
    public void setGamma(int gamma) {
        sendMessage(obtainMessage(MSG_GAMMA, gamma));
    }

    @Override
    public int getGamma() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getGamma();
    }

    @Override
    public void resetSaturation() {
        sendEmptyMessage(MSG_RESET_SATURATION);
    }

    @Override
    public void setSaturation(int saturation) {
        sendMessage(obtainMessage(MSG_SATURATION, saturation));
    }

    @Override
    public int getSaturation() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getSaturation();
    }

    @Override
    public void resetHue() {
        sendEmptyMessage(MSG_RESET_HUE);
    }

    @Override
    public void setHue(int hue) {
        sendMessage(obtainMessage(MSG_HUE, hue));
    }

    @Override
    public int getHue() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getHue();
    }

    @Override
    public void setPowerlineFrequency(int frequency) {
        sendMessage(obtainMessage(MSG_POWERLINEFREQUENCY, frequency));
    }

    @Override
    public int getPowerlineFrequency() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getPowerlineFrequency();
    }

    @Override
    public void resetZoom() {
        sendEmptyMessage(MSG_RESET_ZOOM);
    }

    @Override
    public void setZoom(int zoom) {
        sendMessage(obtainMessage(MSG_ZOOM, zoom));
    }

    @Override
    public int getZoom() {
        final CameraThread thread = mWeakThread.get();
        return thread == null ? 0 : thread.getZoom();
    }


    private static final class CameraThread extends Thread implements IUSBCameraParamsControll {
        private static final String TAG_THREAD = "CameraThread";
        private final Object mSync = new Object();
        private final WeakReference<Context> mWeakParent;
        private boolean mIsRecording = false;
        private boolean isCaptureStill = false;
        /**
         * shutter sound
         */
        private SoundPool mSoundPool;
        private int mSoundId;
        private CameraHandler mHandler;
        /**
         * for accessing UVC camera
         */
        private UVCCamera mUVCCamera;
        private String supportedSize;
        /**
         * muxer for audio/video recording
         */

        private MediaMuxerWrapper mMuxer;

        private MediaVideoEncoder videoEncoder;

        private CameraThread(final Context parent) {
            super("CameraThread");
            mWeakParent = new WeakReference<Context>(parent);
            loadSutterSound(parent);
        }

        @Override
        protected void finalize() throws Throwable {
            Log.i(TAG, "CameraThread#finalize");
            super.finalize();
        }

        public CameraHandler getHandler() {
            if (DEBUG) Log.v(TAG_THREAD, "getHandler:");
            synchronized (mSync) {
                if (mHandler == null)
                    try {
                        mSync.wait();
                    } catch (final InterruptedException e) {
                    }
            }
            return mHandler;
        }

        public boolean isCameraOpened() {
            return mUVCCamera != null;
        }

        public boolean isRecording() {
            return mIsRecording;
        }

        public UVCCamera getUVCCamera() {
            return mUVCCamera;
        }

        public String getSupportedSize() {
            return supportedSize;
        }

        public void handleOpen(final USBMonitor.UsbControlBlock ctrlBlock) {
            if (DEBUG) Log.v(TAG_THREAD, "handleOpen:");
//            handleClose();
            mUVCCamera = new UVCCamera();
            mUVCCamera.open(ctrlBlock);
            supportedSize = mUVCCamera.getSupportedSize();
            if (DEBUG) Log.i(TAG, "supportedSize:" + supportedSize);
        }

        public void handleClose() {
            if (DEBUG) Log.v(TAG_THREAD, "handleClose:");
            handleStopRecording();
            if (mUVCCamera != null) {
                mUVCCamera.stopPreview();
                mUVCCamera.destroy();
                mUVCCamera = null;
            }
        }

        public void handleStartPreview(final Surface surface) {
            if (DEBUG) Log.v(TAG_THREAD, "handleStartPreview:");
            if (mUVCCamera == null) return;
            try {
                mUVCCamera.setPreviewSize(PREVIEW_WIDTH, PREVIEW_HEIGHT, PREVIEW_MODE);
                Log.d(TAG, "fallback to MJPEG mode");
            } catch (final IllegalArgumentException e) {
                try {
                    // fallback to YUV mode
                    Log.d(TAG, "fallback to YUV mode");
                    mUVCCamera.setPreviewSize(PREVIEW_WIDTH, PREVIEW_HEIGHT, UVCCamera.DEFAULT_PREVIEW_MODE);
                } catch (final IllegalArgumentException e1) {
                    handleClose();
                }
            }
            if (mUVCCamera != null) {
                if (is3D) {
                    CAPTURE_WIDTH = PREVIEW_WIDTH;
                    RECORD_WIDTH = PREVIEW_WIDTH;
//                    mUVCCamera.setFrameCallback(mIFrameCallback, UVCCamera.PIXEL_FORMAT_NV21_HALF);
                    mUVCCamera.setFrameCallback(mIFrameCallback, UVCCamera.PIXEL_FORMAT_NV21);
                } else {
                    if (BuildConfig.DEBUG) Log.d(TAG, "==========");
                    mUVCCamera.setFrameCallback(mIFrameCallback, UVCCamera.PIXEL_FORMAT_RGBX);

                }
                mUVCCamera.setPreviewDisplay(surface);
                mUVCCamera.startPreview();
            }
        }

        public void handleStopPreview() {
            if (DEBUG) Log.v(TAG_THREAD, "handleStopPreview:");
            if (mUVCCamera != null) {
                mUVCCamera.stopPreview();
            }
            synchronized (mSync) {
                mSync.notifyAll();
            }
        }

        public void handleCaptureStill() {
            isCaptureStill = true;
            if (DEBUG) Log.v(TAG_THREAD, "handleCaptureStill:");
        }

        public void handleStartRecording() {
            if (DEBUG) Log.v(TAG_THREAD, "handleStartRecording:");

            if (mMuxer == null) {
                try {
                    mMuxer = new MediaMuxerWrapper(".mp4");    // if you record audio only, ".m4a" is also OK.
                    if (true) {
                        // for video capturing
                        videoEncoder = new MediaVideoEncoder(mMuxer, mMediaEncoderListener, RECORD_WIDTH, RECORD_HEIGHT);
                    }
                    if (true) {
                        // for audio capturing
//                        new MediaAudioEncoder(mMuxer, mMediaEncoderListener);
                    }
                    mMuxer.prepare();
                    mMuxer.startRecording();
                } catch (final IOException e) {
                }
            }
        }

        public void handleStopRecording() {
            if (DEBUG) Log.v(TAG_THREAD, "handleStopRecording:");
            if (mMuxer != null) {
                mMuxer.stopRecording();
                mMuxer = null;
                // you should not wait here
            }
        }

        /**
         * callback methods from encoder
         */
        private final MediaEncoder.MediaEncoderListener mMediaEncoderListener = new MediaEncoder.MediaEncoderListener() {
            @Override
            public void onPrepared(final MediaEncoder encoder) {
                mIsRecording = true;
            }

            @Override
            public void onStopped(final MediaEncoder encoder) {
                mIsRecording = false;
            }
        };


        public void handleUpdateMedia(final String path) {
            if (DEBUG) Log.v(TAG_THREAD, "handleUpdateMedia:path=" + path);
            final Context parent = mWeakParent.get();
            if (parent != null && parent.getApplicationContext() != null) {
                try {
                    if (DEBUG) Log.i(TAG, "MediaScannerConnection#scanFile");
                    MediaScannerConnection.scanFile(parent.getApplicationContext(), new String[]{path}, null, null);
                } catch (final Exception e) {
                    Log.e(TAG, "handleUpdateMedia:", e);
                }
            } else {
                Log.w(TAG, "MainActivity already destroyed");
                // give up to add this movice to MediaStore now.
                // Seeing this movie on Gallery app etc. will take a lot of time.
                handleRelease();
            }
        }

        public void handleRelease() {
            if (DEBUG) Log.v(TAG_THREAD, "handleRelease:");
            handleClose();
            if (!mIsRecording)
                Looper.myLooper().quit();
        }

        /**
         * prepare and load shutter sound for still image capturing
         */
        @SuppressWarnings("deprecation")
        private void loadSutterSound(final Context context) {
            // get system stream type using refrection
            int streamType;
            try {
                final Class<?> audioSystemClass = Class.forName("android.media.AudioSystem");
                final Field sseField = audioSystemClass.getDeclaredField("STREAM_SYSTEM_ENFORCED");
                streamType = sseField.getInt(null);
            } catch (final Exception e) {
                streamType = AudioManager.STREAM_SYSTEM;    // set appropriate according to your app policy
            }
            if (mSoundPool != null) {
                try {
                    mSoundPool.release();
                } catch (final Exception e) {
                }
                mSoundPool = null;
            }
            // load sutter sound from resource
            mSoundPool = new SoundPool(2, streamType, 0);
            mSoundId = mSoundPool.load(context, R.raw.camera_click, 1);
        }

        private void captureStill(ByteBuffer frame) {
            Log.d(TAG, "onFrame Capture still");
            mSoundPool.play(mSoundId, 0.2f, 0.2f, 0, 0, 1.0f);
            File outputFile = null;
            BufferedOutputStream os = null;

            int rgb[] = new int[CAPTURE_WIDTH * CAPTURE_HEIGHT];

            try {
                File tempFile = new File(Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_DCIM), DIR_NAME);
                if (!tempFile.exists()) {
                    tempFile.mkdirs();
                }
                outputFile = new File(tempFile, System.currentTimeMillis() + ".jpg");
                os = new BufferedOutputStream(new FileOutputStream(outputFile));
                byte buf[] = new byte[CAPTURE_WIDTH * CAPTURE_HEIGHT * 3 / 2];
                frame.get(buf);

//                ImageProc.decodeYUV420SP(rgb, buf, CAPTURE_WIDTH, CAPTURE_HEIGHT);
//                Bitmap bitmap = Bitmap.createBitmap(rgb, CAPTURE_WIDTH, CAPTURE_HEIGHT, Bitmap.Config.ARGB_8888);
                ByteArrayOutputStream out = new ByteArrayOutputStream();
                YuvImage yuvImage = new YuvImage(buf, ImageFormat.NV21, CAPTURE_WIDTH, CAPTURE_HEIGHT, null);
                yuvImage.compressToJpeg(new Rect(0, 0, CAPTURE_WIDTH, CAPTURE_HEIGHT), 50, out);
                byte[] imageBytes = out.toByteArray();
                Bitmap bitmap = BitmapFactory.decodeByteArray(imageBytes, 0, imageBytes.length);
                bitmap.compress(Bitmap.CompressFormat.JPEG, 100, os);
                os.flush();
                os.close();
                mHandler.sendMessage(mHandler.obtainMessage(MSG_MEDIA_UPDATE, outputFile.getPath()));
            } catch (Exception e) {
                Log.e(TAG, "onFrame Capture still Error Exception");
            }

            isCaptureStill = false;
        }

        // if you need frame data as ByteBuffer on Java side, you can use this callback method with UVCCamera#setFrameCallback
        private final IFrameCallback mIFrameCallback = new IFrameCallback() {
            @Override
            public void onFrame(final ByteBuffer frame) {
//                Log.d(TAG,"onFrame");
                if (isCaptureStill) {
//                    captureStill(frame);
                }
                if (isRecording()) {
                    long startTime = System.currentTimeMillis();
                    byte buf[] = new byte[RECORD_WIDTH * RECORD_HEIGHT * 3 / 2];
                    frame.get(buf);
                    videoEncoder.encodeFrame(buf);
                    long endTime = System.currentTimeMillis();
                    Log.i(TAG, Integer.toString((int) (endTime - startTime)) + "ms");
                }
            }
        };


        @Override
        public void run() {
            Looper.prepare();
            synchronized (mSync) {
                mHandler = new CameraHandler(this);
                mSync.notifyAll();
            }
            Looper.loop();
            synchronized (mSync) {
                mHandler = null;
                mSoundPool.release();
                mSoundPool = null;
                mSync.notifyAll();
            }
        }

        @Override
        public void resetContrast() {
            mUVCCamera.resetContrast();
        }

        @Override
        public void setContrast(int contrast) {
            mUVCCamera.setContrast(contrast);
        }

        @Override
        public int getContrast() {
            return mUVCCamera.getContrast();
        }

        @Override
        public void resetSharpness() {
            mUVCCamera.resetSharpness();
        }

        @Override
        public void setSharpness(int sharpness) {
            mUVCCamera.setSharpness(sharpness);
        }

        @Override
        public int getSharpness() {
            return mUVCCamera.getSharpness();
        }

        @Override
        public void resetGain() {
            mUVCCamera.resetGain();
        }

        @Override
        public void setGain(int gain) {
            mUVCCamera.setGain(gain);
        }

        @Override
        public int getGain() {
            return mUVCCamera.getGain();
        }

        @Override
        public void resetGamma() {
            mUVCCamera.resetGamma();
        }

        @Override
        public void setGamma(int gamma) {
            mUVCCamera.setGamma(gamma);
        }

        @Override
        public int getGamma() {
            return mUVCCamera.getGamma();
        }

        @Override
        public void resetSaturation() {
            mUVCCamera.resetSaturation();
        }

        @Override
        public void setSaturation(int saturation) {
            mUVCCamera.setSaturation(saturation);
        }

        @Override
        public int getSaturation() {
            return mUVCCamera.getSaturation();
        }

        @Override
        public void resetHue() {
            mUVCCamera.resetHue();
        }

        @Override
        public void setHue(int hue) {
            mUVCCamera.setHue(hue);
        }

        @Override
        public int getHue() {
            return mUVCCamera.getHue();
        }

        @Override
        public void setPowerlineFrequency(int frequency) {
            mUVCCamera.setPowerlineFrequency(frequency);
        }

        @Override
        public int getPowerlineFrequency() {
            return mUVCCamera.getPowerlineFrequency();
        }

        @Override
        public void resetZoom() {
            mUVCCamera.resetZoom();
        }

        @Override
        public void setZoom(int zoom) {
            mUVCCamera.setZoom(zoom);
        }

        @Override
        public int getZoom() {
            return mUVCCamera.getZoom();
        }

        @Override
        public void setExposureMode(int exposureMode) {

        }

        @Override
        public int getExposureMode() {
            return 0;
        }

        @Override
        public void setAutoFocus(boolean autofocus) {
            mUVCCamera.setAutoFocus(autofocus);
        }

        @Override
        public boolean getAutoFocus() {
            return mUVCCamera.getAutoFocus();
        }

        @Override
        public void resetFocus() {
            mUVCCamera.resetFocus();
        }

        @Override
        public void setFocus(int focus) {
            mUVCCamera.setFocus(focus);
        }

        @Override
        public int getFocus() {
            return mUVCCamera.getFocus();
        }

        @Override
        public void setAutoWhiteBlance(boolean autoWhiteBlance) {
            mUVCCamera.setAutoWhiteBlance(autoWhiteBlance);
        }

        @Override
        public boolean getAutoWhiteBlance() {
            return mUVCCamera.getAutoWhiteBlance();
        }

        @Override
        public void resetWhiteBlance() {
            mUVCCamera.resetWhiteBlance();
        }

        @Override
        public void setWhiteBlance(int whiteBlance) {
            mUVCCamera.setWhiteBlance(whiteBlance);
        }

        @Override
        public int getWhiteBlance() {
            return mUVCCamera.getWhiteBlance();
        }

        @Override
        public void resetBrightness() {
            mUVCCamera.resetBrightness();
        }

        @Override
        public void setBrightness(int brightness) {
            mUVCCamera.setBrightness(brightness);
        }

        @Override
        public int getBrightness() {
            return mUVCCamera.getBrightness();
        }
    }
}