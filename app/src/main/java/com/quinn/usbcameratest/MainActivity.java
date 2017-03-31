package com.quinn.usbcameratest;

import android.support.v7.app.AppCompatActivity;

public class MainActivity extends AppCompatActivity /*implements View.OnClickListener*/ {
//
//    private UVCCameraTextureView mCameraView;
//    private Button open,phone,video;
//    private USBCamera camera;
//
////    1280,720
////    private static int PREVIEW_WIDTH = 1472;
////    private static int PREVIEW_HEIGHT = 736;
//    private static int PREVIEW_WIDTH = 1280;
//    private static int PREVIEW_HEIGHT = 720;
//
//    @Override
//    protected void onCreate(Bundle savedInstanceState) {
//        super.onCreate(savedInstanceState);
//        setContentView(R.layout.activity_main);
//
//        mCameraView = (UVCCameraTextureView) findViewById(R.id.uvc_camera);
//
//        camera = new USBCamera();
//        camera.init(this);
//        camera.setCameraType(CameraType.C3D_SBS);
//        mCameraView.setAspectRatio(PREVIEW_WIDTH / 2 / (float) PREVIEW_HEIGHT);
//
//
//        open = (Button) findViewById(R.id.open);
//        phone = (Button) findViewById(R.id.phone);
//        video = (Button) findViewById(R.id.video);
//        open.setOnClickListener(this);
//        phone.setOnClickListener(this);
//        video.setOnClickListener(this);
//
//    }
//
//    @Override
//    protected void onStop() {
//        super.onStop();
//        camera.close();
//    }
//
//    @Override
//    protected void onDestroy() {
//        super.onDestroy();
//        camera.destroy();
//    }
//
//    @Override
//    public void onClick(View v) {
//        switch (v.getId()){
//            case R.id.open:
//                if (!camera.isCameraOpened()) {
//                    boolean ret = camera.open(0);
//                    if (!ret) {
//                        Toast.makeText(MainActivity.this, "NO_USB_DEVICE", Toast.LENGTH_SHORT).show();
//                    } else {
//                        camera.setPreviewSize(PREVIEW_WIDTH, PREVIEW_HEIGHT);
//                        camera.setPreviewTexture(mCameraView.getSurfaceTexture());
//                        camera.startPreview();
//                    }
//                }
//                break;
//            case R.id.phone:
//                if (camera.isCameraOpened()) {
//                    camera.captureStill();
//                }
//                break;
//            case R.id.video:
//                if (camera.isRecording()) {
//                    camera.stopRecording();
//                }else{
//                    camera.startRecording();
//                }
//                break;
//        }
//    }
}
