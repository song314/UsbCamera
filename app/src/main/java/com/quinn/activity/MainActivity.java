package com.quinn.activity;

import android.Manifest;
import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.IntentFilter;
import android.content.pm.ActivityInfo;
import android.graphics.Bitmap;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.os.Handler;
import android.os.PowerManager;
import android.support.v7.app.AlertDialog;
import android.support.v7.widget.AppCompatImageButton;
import android.support.v7.widget.AppCompatSeekBar;
import android.support.v7.widget.Toolbar;
import android.util.DisplayMetrics;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.Menu;
import android.view.MenuItem;
import android.view.Surface;
import android.view.View;
import android.view.Window;
import android.view.WindowManager;
import android.widget.RelativeLayout;
import android.widget.SeekBar;
import android.widget.Toast;

import com.asha.vrlib.MDVRLibrary;
import com.dreamguard.api.CameraType;
import com.dreamguard.api.USBCamera;
import com.dreamguard.widget.UVCCameraTextureView;
import com.google.android.apps.muzei.render.GLTextureView;
import com.google.gson.Gson;
import com.orhanobut.logger.Logger;
import com.quinn.model.Format;
import com.quinn.usbcameratest.R;
import com.tbruyelle.rxpermissions.RxPermissions;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.List;

import rx.functions.Action1;

public class MainActivity extends BaseActivity implements View.OnClickListener, SeekBar.OnSeekBarChangeListener {

    private static final String USB_STATE_RECEIVER = "android.hardware.usb.action.USB_STATE";

    private Toolbar toolbar;
    private GLTextureView mCameraView;
    private AppCompatSeekBar brightness, cBrightness, cContrast, cSaturation, cHue, pFocus, pWhiteBlance, pSharpness, pZoom;
    private AppCompatImageButton switchPhotoVideo, photo, video, autoFocus, autoWhiteBlance, colors, params;
    private RelativeLayout cameraParent, tip;
    private USBCamera camera;
    String resolution[];

    //["1472x736","2176x1088","3008x1504"]    1280,720    640,480
    private static int PREVIEW_WIDTH = 1472;
    private static int PREVIEW_HEIGHT =736;
//    private static int PREVIEW_WIDTH = 1280;
//    private static int PREVIEW_HEIGHT = 720;

    private Handler mHandler;

    private boolean isPhotoVideo;
    private RxPermissions rxPermissions;
    private PowerManager.WakeLock mWakeLock;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        PowerManager pm = (PowerManager) getSystemService(Context.POWER_SERVICE);
        mWakeLock = pm.newWakeLock(PowerManager.SCREEN_BRIGHT_WAKE_LOCK, MainActivity.class.getSimpleName());
        mHandler = new Handler();
        rxPermissions = new RxPermissions(this);

        initView();
        initUVCCamera();
    }

    private void initView() {
        toolbar = (Toolbar) findViewById(R.id.toolBar);
        setSupportActionBar(toolbar);

        mCameraView = (GLTextureView) findViewById(R.id.uvc_camera);

        switchPhotoVideo = (AppCompatImageButton) findViewById(R.id.switchPhotoVideo);
        photo = (AppCompatImageButton) findViewById(R.id.photo);
        video = (AppCompatImageButton) findViewById(R.id.video);
        autoFocus = (AppCompatImageButton) findViewById(R.id.autoFocus);
        autoWhiteBlance = (AppCompatImageButton) findViewById(R.id.autoWhiteBlance);
        colors = (AppCompatImageButton) findViewById(R.id.colors);
        params = (AppCompatImageButton) findViewById(R.id.params);
        brightness = (AppCompatSeekBar) findViewById(R.id.seekBarBrightness);

        cameraParent = (RelativeLayout) findViewById(R.id.cameraParent);
        tip = (RelativeLayout) findViewById(R.id.tip);

        photo.setOnClickListener(this);
        switchPhotoVideo.setOnClickListener(this);
        video.setOnClickListener(this);
        autoFocus.setOnClickListener(this);
        autoWhiteBlance.setOnClickListener(this);
        colors.setOnClickListener(this);
        params.setOnClickListener(this);

        brightness.setOnSeekBarChangeListener(this);
    }

    private void initUVCCamera() {
        camera = new USBCamera();
        camera.init(this, tip);
        camera.setCameraType(CameraType.C3D_SBS);
    }

    private void openCamera(final Surface surface) {
//        tip.setVisibility(View.GONE);
//        mCameraView.setAspectRatio(PREVIEW_WIDTH / (float) PREVIEW_HEIGHT);

        mHandler.postDelayed(new Runnable() {
            @Override
            public void run() {
                if (!camera.isCameraOpened()) {
                    boolean ret = camera.open(0);
                    if (ret) {
//                        tip.setVisibility(View.GONE);
                        camera.setPreviewSize(PREVIEW_WIDTH, PREVIEW_HEIGHT);
                        camera.setPreviewTexture(surface);
                        camera.startPreview();

                        mHandler.postDelayed(new Runnable() {
                            @Override
                            public void run() {
                                try {
//                                    getCameraParams();
                                    Format formats = new Gson().fromJson(camera.getSupportedSize(), Format.class);
                                    List<String> sizeList = formats.getFormats().get(0).getSize();
                                    resolution = sizeList.toArray(new String[sizeList.size()]);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        }, 1500);
                    }
                }
            }
        }, 100);
    }

    private MDVRLibrary mVRLibrary;
    private void initVRLibrary(){
        // new instance
        mVRLibrary = MDVRLibrary.with(this)
                .displayMode(MDVRLibrary.DISPLAY_MODE_NORMAL)
                .interactiveMode(MDVRLibrary.INTERACTIVE_MODE_MOTION)
                .asVideo(new MDVRLibrary.IOnSurfaceReadyCallback() {
                    @Override
                    public void onSurfaceReady(Surface surface) {
                        openCamera(surface);
                    }
                })
                .build(R.id.uvc_camera);
    }

    @Override
    protected void onResume() {
        super.onResume();
        initVRLibrary();
        mWakeLock.acquire();
    }

    @Override
    protected void onPause() {
        super.onPause();
        mWakeLock.release();
        camera.close();
    }

    private void getCameraParams() {
        brightness.setProgress(camera.getBrightness());

        if (camera.getAutoFocus()) {
            autoFocus.setImageResource(R.drawable.ic_focus_auto);
        } else {
            autoFocus.setImageResource(R.drawable.ic_focus_normal);
        }
        if (camera.getAutoWhiteBlance()) {
            autoWhiteBlance.setImageResource(R.drawable.ic_whiteblance_auto);
        } else {
            autoWhiteBlance.setImageResource(R.drawable.ic_whiteblance_normal);
        }
    }

    @Override
    protected void onStop() {
        super.onStop();
        camera.close();
    }

    @Override
    protected void onDestroy() {
        super.onDestroy();
        camera.destroy();
        mVRLibrary.onDestroy();
        if (mHandler != null) {
            mHandler = null;
        }
    }

    @Override
    public void onClick(View v) {
        switch (v.getId()) {
            case R.id.switchPhotoVideo:
                if (isPhotoVideo) {
                    photo.setVisibility(View.VISIBLE);
                    video.setVisibility(View.GONE);
                    switchPhotoVideo.setImageResource(R.drawable.ic_video_switch);
                } else {
                    photo.setVisibility(View.GONE);
                    video.setVisibility(View.VISIBLE);
                    switchPhotoVideo.setImageResource(R.drawable.ic_camera_switch);
                }
                isPhotoVideo = !isPhotoVideo;
                break;
            case R.id.photo:
                if (!camera.isCameraOpened()) return;
                rxPermissions.request(Manifest.permission.WRITE_EXTERNAL_STORAGE).subscribe(new Action1<Boolean>() {
                    @Override
                    public void call(Boolean aBoolean) {
                        if (aBoolean) {
//                            camera.captureStill(mCameraView);
                        } else {
                            Toast.makeText(MainActivity.this, "无法获得权限,请到设置中开启权限!", Toast.LENGTH_SHORT).show();
                        }
                    }
                });
                break;
            case R.id.video:
                if (!camera.isCameraOpened()) return;
                rxPermissions.request(Manifest.permission.WRITE_EXTERNAL_STORAGE, Manifest.permission.RECORD_AUDIO).subscribe(new Action1<Boolean>() {
                    @Override
                    public void call(Boolean aBoolean) {
                        if (aBoolean) {
                            if (camera.isRecording()) {
                                camera.stopRecording();
                                video.setImageResource(R.drawable.controll_video);
                            } else {
                                camera.startRecording();
                                video.setImageResource(R.drawable.controll_video_open);
                            }
                        } else {
                            Toast.makeText(MainActivity.this, "无法获得权限,请到设置中开启权限!", Toast.LENGTH_SHORT).show();
                        }
                    }
                });
                break;
//            case R.id.autoFocus:
//                if (!camera.isCameraOpened()) return;
//                if (camera.getAutoFocus()) {
//                    autoFocus.setImageResource(R.drawable.ic_focus_normal);
//                } else {
//                    autoFocus.setImageResource(R.drawable.ic_focus_auto);
//                }
//                camera.setAutoFocus(!camera.getAutoFocus());
//                break;
//            case R.id.autoWhiteBlance:
//                if (!camera.isCameraOpened()) return;
//                if (camera.getAutoWhiteBlance()) {
//                    autoWhiteBlance.setImageResource(R.drawable.ic_whiteblance_normal);
//                } else {
//                    autoWhiteBlance.setImageResource(R.drawable.ic_whiteblance_auto);
//                }
//                camera.setAutoWhiteBlance(!camera.getAutoWhiteBlance());
//                break;
            case R.id.colors:
                showColorsDialog();
                break;
            case R.id.params:
                showParamsDialog();
//                Intent intent = new Intent(Intent.ACTION_GET_CONTENT);
//                intent.setDataAndType(Uri.fromFile(new File("/sdcard/Movis/")),"*/*");//设置类型，我这里是任意类型，任意后缀的可以这样写。
//                intent.addCategory(Intent.CATEGORY_OPENABLE);
//                startActivityForResult(intent,1);
                break;
        }
    }

    private void showColorsDialog() {
        View view = LayoutInflater.from(this).inflate(R.layout.dialog_colors, null);
        cBrightness = (AppCompatSeekBar) view.findViewById(R.id.color_brightness);
        cBrightness.setOnSeekBarChangeListener(this);
        cContrast = (AppCompatSeekBar) view.findViewById(R.id.color_contrast);
        cContrast.setOnSeekBarChangeListener(this);
        cSaturation = (AppCompatSeekBar) view.findViewById(R.id.color_saturation);
        cSaturation.setOnSeekBarChangeListener(this);
        cHue = (AppCompatSeekBar) view.findViewById(R.id.color_hue);
        cHue.setOnSeekBarChangeListener(this);

        cBrightness.setProgress(camera.getBrightness());
        cContrast.setProgress(camera.getContrast());
        cSaturation.setProgress(camera.getSaturation());
        cHue.setProgress(camera.getHue());

        AlertDialog.Builder dialog = new AlertDialog.Builder(this, R.style.MyDialog);
        dialog.setView(view);
        AlertDialog alertDialog = dialog.create();
        alertDialog.show();
        Window window = alertDialog.getWindow();
        WindowManager.LayoutParams params = window.getAttributes();
        DisplayMetrics outMetrics = new DisplayMetrics();
        getWindowManager().getDefaultDisplay().getMetrics(outMetrics);
        params.width = (int) (outMetrics.widthPixels * 0.8);
        window.setAttributes(params);
    }

    private void showParamsDialog() {
        View view = LayoutInflater.from(this).inflate(R.layout.dialog_params, null);
        pFocus = (AppCompatSeekBar) view.findViewById(R.id.param_focus);
        pFocus.setOnSeekBarChangeListener(this);
        pWhiteBlance = (AppCompatSeekBar) view.findViewById(R.id.param_wihteblance);
        pWhiteBlance.setOnSeekBarChangeListener(this);
        pSharpness = (AppCompatSeekBar) view.findViewById(R.id.param_sharpness);
        pSharpness.setOnSeekBarChangeListener(this);
        pZoom = (AppCompatSeekBar) view.findViewById(R.id.param_zoom);
        pZoom.setOnSeekBarChangeListener(this);

        pFocus.setProgress(camera.getFocus());
        pWhiteBlance.setProgress(camera.getWhiteBlance());
        pSharpness.setProgress(camera.getSharpness());
        pZoom.setProgress(camera.getZoom());

        AlertDialog.Builder dialog = new AlertDialog.Builder(this, R.style.MyDialog);
        dialog.setView(view);
        AlertDialog alertDialog = dialog.create();
        alertDialog.show();
        Window window = alertDialog.getWindow();
        WindowManager.LayoutParams params = window.getAttributes();
        DisplayMetrics outMetrics = new DisplayMetrics();
        getWindowManager().getDefaultDisplay().getMetrics(outMetrics);
        params.width = (int) (outMetrics.widthPixels * 0.8);
        window.setAttributes(params);
    }

    private void showSelectResolutionDialog() {
        final AlertDialog.Builder builder = new AlertDialog.Builder(this);
        builder.setTitle("分辨率选择");
//        Logger.d(resolution);
        if (camera.isCameraOpened())
            builder.setSingleChoiceItems(resolution, -1, new DialogInterface.OnClickListener() {
                @Override
                public void onClick(DialogInterface dialogInterface, int i) {
                    String size = resolution[i];
                    PREVIEW_WIDTH = Integer.parseInt(size.substring(0, size.indexOf("x")));
                    PREVIEW_HEIGHT = Integer.parseInt(size.substring(size.indexOf("x") + 1));
                    camera.close();
//                    openCamera();
                    dialogInterface.dismiss();
                }
            });
        builder.setNegativeButton("取消", new DialogInterface.OnClickListener() {
            @Override
            public void onClick(DialogInterface dialogInterface, int i) {
                dialogInterface.dismiss();
            }
        });
        builder.create().show();
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        getMenuInflater().inflate(R.menu.main, menu);
        return super.onCreateOptionsMenu(menu);
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        switch (item.getItemId()) {
            case R.id.setting:
                startActivity(new Intent(this, PreActivity.class));
                break;
//            case R.id.resolution:
//                showSelectResolutionDialog();
//                break;
        }
        return super.onOptionsItemSelected(item);
    }

    @Override
    public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
        if (fromUser) {
            switch (seekBar.getId()) {
                case R.id.seekBarBrightness:
                    if (camera.isCameraOpened()) {
                        camera.setBrightness(progress);
                    } else {
                        brightness.setProgress(0);
                    }
                    break;
                case R.id.color_brightness:
                    if (camera.isCameraOpened()) {
                        camera.setBrightness(progress);
                    } else {
                        cBrightness.setProgress(0);
                    }
                    break;
                case R.id.color_contrast:
                    if (camera.isCameraOpened()) {
                        camera.setContrast(progress);
                    } else {
                        cContrast.setProgress(0);
                    }
                    break;
                case R.id.color_saturation:
                    if (camera.isCameraOpened()) {
                        camera.setSaturation(progress);
                    } else {
                        cSaturation.setProgress(0);
                    }
                    break;
                case R.id.color_hue:
                    if (camera.isCameraOpened()) {
                        camera.setHue(progress);
                    } else {
                        cHue.setProgress(0);
                    }
                    break;
                case R.id.param_focus:
                    if (camera.isCameraOpened()) {
                        camera.setFocus(progress);
                        if (camera.getAutoFocus()) {
                            autoFocus.setImageResource(R.drawable.ic_focus_normal);
                            camera.setAutoFocus(!camera.getAutoFocus());
                        }
                    } else {
                        pFocus.setProgress(0);
                    }
                    break;
                case R.id.param_wihteblance:
                    if (camera.isCameraOpened()) {
                        camera.setWhiteBlance(progress);
                        if (camera.getAutoWhiteBlance()) {
                            autoWhiteBlance.setImageResource(R.drawable.ic_whiteblance_normal);
                            camera.setAutoWhiteBlance(!camera.getAutoWhiteBlance());
                        }
                    } else {
                        pWhiteBlance.setProgress(0);
                    }
                    break;
                case R.id.param_sharpness:
                    if (camera.isCameraOpened()) {
                        camera.setSharpness(progress);
                    } else {
                        pSharpness.setProgress(0);
                    }
                    break;
                case R.id.param_zoom:
                    if (camera.isCameraOpened()) {
                        camera.setZoom(progress);
                    } else {
                        pZoom.setProgress(0);
                    }
                    break;
            }
        }
    }

    @Override
    public void onStartTrackingTouch(SeekBar seekBar) {
    }

    @Override
    public void onStopTrackingTouch(SeekBar seekBar) {

    }

}
