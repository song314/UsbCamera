package com.dreamguard.api;

/**
 * Created by sundy on 2017/3/20 0020.
 */

public interface IUSBCameraParamsControll {

    int MSG_EXPOSUREMODE= 10;
    int MSG_RESET_EXPOSUREMODE= 11;
    int MSG_AUTOFOCUS= 12;
    int MSG_RESET_AUTOFOCUS= 13;
    int MSG_FOCUS= 14;
    int MSG_RESET_FOCUS= 15;
    int MSG_AUTOWHITEBLANCE= 16;
    int MSG_RESET_AUTOWHITEBLANCE= 17;
    int MSG_WHITEBLANCE= 18;
    int MSG_RESET_WHITEBLANCE= 19;
    int MSG_BRIGHTNESS= 20;
    int MSG_RESET_BRIGHTNESS= 21;
    int MSG_CONTRAST= 22;
    int MSG_RESET_CONTRAST= 23;
    int MSG_SHARPNESS= 24;
    int MSG_RESET_SHARPNESS= 25;
    int MSG_GAIN= 26;
    int MSG_RESET_GAIN= 27;
    int MSG_GAMMA= 28;
    int MSG_RESET_GAMMA= 29;
    int MSG_SATURATION= 30;
    int MSG_RESET_SATURATION= 31;
    int MSG_HUE= 32;
    int MSG_RESET_HUE= 33;
    int MSG_POWERLINEFREQUENCY= 34;
    int MSG_RESET_POWERLINEFREQUENCY= 35;
    int MSG_ZOOM= 36;
    int MSG_RESET_ZOOM= 37;


    void setExposureMode(int exposureMode);
    int getExposureMode();

    /**
     * 自动对焦
     * @param autofocus
     */
    void setAutoFocus(boolean autofocus);
    boolean getAutoFocus();

    /**
     * 焦距
     */
    void resetFocus();
    void setFocus(int focus);
    int getFocus();

    /**
     * 自动白平衡
     * @param autoWhiteBlance
     */
    void setAutoWhiteBlance(boolean autoWhiteBlance);
    boolean getAutoWhiteBlance();

    /**
     * 白平衡
     */
    void resetWhiteBlance();
    void setWhiteBlance(int whiteBlance);
    int getWhiteBlance();

    /**
     * 亮度
     */
    void resetBrightness();
    void setBrightness(int brightness);
    int getBrightness();

    /**
     * 对比度
     */
    void resetContrast();
    void setContrast(int contrast);
    int getContrast();

    /**
     * 清晰度
     */
    void resetSharpness();
    void setSharpness(int sharpness);
    int getSharpness();

    void resetGain();
    void setGain(int gain);
    int getGain();

    void resetGamma();
    void setGamma(int gamma);
    int getGamma();

    /**
     * 饱和度
     */
    void resetSaturation();
    void setSaturation(int saturation);
    int getSaturation();

    /**
     * 色相/色调
     */
    void resetHue();
    void setHue(int hue);
    int getHue();

    void setPowerlineFrequency(int frequency);
    int getPowerlineFrequency();

    /**
     * 缩放
     */
    void resetZoom();
    void setZoom(int zoom);
    int getZoom();
}
