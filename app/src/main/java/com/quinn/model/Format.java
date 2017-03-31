package com.quinn.model;

import com.google.gson.annotations.SerializedName;

import java.util.List;

/**
 * Created by Administrator on 2017/3/2.
 */

public class Format {

    @SerializedName("formats")
    private List<Formats> formats;

    public List<Formats> getFormats() {
        return formats;
    }

    public void setFormats(List<Formats> formats) {
        this.formats = formats;
    }

    public static class Formats {
        /**
         * index : 1
         * type : 4
         * default : 1
         * size : ["640x480","352x288","320x240","176x144","160x120"]
         */

        @SerializedName("index")
        private int index;
        @SerializedName("type")
        private int type;
        @SerializedName("default")
        private int defaultX;
        @SerializedName("size")
        private List<String> size;

        public int getIndex() {
            return index;
        }

        public void setIndex(int index) {
            this.index = index;
        }

        public int getType() {
            return type;
        }

        public void setType(int type) {
            this.type = type;
        }

        public int getDefaultX() {
            return defaultX;
        }

        public void setDefaultX(int defaultX) {
            this.defaultX = defaultX;
        }

        public List<String> getSize() {
            return size;
        }

        public void setSize(List<String> size) {
            this.size = size;
        }
    }
}
