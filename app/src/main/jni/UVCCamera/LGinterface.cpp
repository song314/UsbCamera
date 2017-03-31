#include <utilbase.h>
#include <stdio.h>
#include   "LGinterface.h"


#include <time.h>
#include <ctime>





void LG::Panointerface( void* data, int width, int height, int channel , Frame *out_frame,int flag_wb) {
//	LOGW("width height %d %d\n",width,height);

	int w = para.correct_width;
	int h = para.correct_height;

	Frame  imgcor1;
	IniteFrame( imgcor1, w, h, 3 );
	Frame imgcor2;
	IniteFrame( imgcor2, w, h, 3 );

	Frame imMerge1;
	IniteFrame( imMerge1, para.merge1_width, para.merge1_height, 3 );
	Frame imMerge2;
	IniteFrame( imMerge2, para.merge2_width, para.merge2_height, 3 );

	int hh = imMerge1.height;
	int ww = imMerge1.width;
	Frame imMer1;
	IniteFrame( imMer1, ( int ) ( ww >> 1 ), hh, 3 );
	Frame imMer2;
	IniteFrame( imMer2, ( int ) ( ww >> 1 ), hh, 3 );

	Frame imsrc;
	IniteFrame( imsrc, width, height, channel );

    clock_t start = clock();
	//ImageCorrect( coordinate1, distance1, data, width, height, channel, imgcor1, coordinate2, distance2, imgcor2, h, w );
    ImageCorrect_one( coordinate1, distance1, data, width, height, channel, imgcor1, h, w );
    ImageCorrect_one( coordinate2, distance2, data, width, height, channel, imgcor2, h, w );
    clock_t finish = clock();
    //LOGW("per frame all: %d", (finish-start)/(CLOCKS_PER_SEC/1000));

    if ( flag_wb )
        imageMerge_wb( imgcor1, imgcor2, para.merge1_x, para.merge1_y, imMerge1, para.cut_width, 0);
    else
        ImMerge( &imgcor1, &imgcor2, para.merge1_x, para.merge1_y, &imMerge1, para.cut_width );

	splitFrame( imMerge1, imMer1, imMer2 );
    if ( flag_wb )
        imageMerge_wb( imMer2, imMer1, para.merge2_x, para.merge2_y, imMerge2, para.cut_width,0 );
    else
        ImMerge( &imMer2, &imMer1, para.merge2_x, para.merge2_y, &imMerge2, para.cut_width );



    uchar *data_ptr = imMerge2.data;
    for(int y = 0; y < para.merge1_height; y ++)
    {
        for(int x = 0; x < para.merge1_width; x++)
        {
            out_frame->data[y*out_frame->width*out_frame->channel+x*out_frame->channel] = data_ptr[y*imMerge2.width*imMerge2.channel+x*imMerge2.channel];
            out_frame->data[y*out_frame->width*out_frame->channel+x*out_frame->channel+1] = data_ptr[y*imMerge2.width*imMerge2.channel+x*imMerge2.channel+1];
            out_frame->data[y*out_frame->width*out_frame->channel+x*out_frame->channel+2] = data_ptr[y*imMerge2.width*imMerge2.channel+x*imMerge2.channel+2];
        }
    }

    DeleteFrame( imgcor1 );
    DeleteFrame( imgcor2 );
    DeleteFrame( imMerge1 );
    DeleteFrame( imMerge2 );
    DeleteFrame( imMer1 );
    DeleteFrame( imMer2 );
    DeleteFrame( imsrc );


}

void LG::LGInit( int w, int h ) {

    //para.imgsrc_W = 1920;
    //para.imgsrc_H = 1080;
    para.imgsrc_W = w;
    para.imgsrc_H = h;
    //para.L_center_x = 482;
    //para.L_center_y = 1080 - 603;
    /*para.L_center_x = 479;
    para.L_center_y = 480;*/
    para.L_center_x = 319;
    para.L_center_y = 320;
    para.L_pitchAngle = 0 * 3.1415926 / 180;
    //para.L_rotationAngle = ( -0.3 - 90 ) * 3.1415926 / 180;
    para.L_rotationAngle = ( 0 - 90 ) * 3.1415926 / 180;
    para.L_flipAngle = 0;

    //para.R_center_x = 1920 - 483;
    //para.R_center_y = 482;
    //para.R_center_x = 1439;
    //para.R_center_y = 480;
    para.R_center_x = 959;
    para.R_center_y = 320;
    para.R_pitchAngle = 0 * 3.1415926 / 180;
    para.R_rotationAngle = ( 0 + 90 ) * 3.1415926 / 180;
    para.R_flipAngle = 0;

    //para.R = 450;
    para.R = 300;
    //para.fov = 190.5882;
    para.fov = 189.69;
    para.r = ( int ) ( para.R * 180 / para.fov );

    para.correct_height = 2 * para.r + 1;
    para.correct_width = 2 * para.R + 1;

    //para.merge1_x = 850;
    /*para.merge1_x = 855;
    para.merge1_y = 0;*/
    para.merge1_x = 570;
    para.merge1_y = 0;
    para.merge1_width = para.correct_width + para.merge1_x;
    para.merge1_height = para.correct_height - abs( para.merge1_y );

    //para.merge2_x = 825;
    //para.merge2_x = 827;      //827, 831
    para.merge2_y = 0;
    para.merge2_x = 551;

    para.merge2_width = 0.5 * para.merge1_width + para.merge2_x;
    para.merge2_height = para.merge1_height - abs( para.merge2_y );

    //para.final_width = para.merge2_width;
    //para.final_height = para.final_width / 2;
    //para.final_height = 850;
    //para.final_width = 1700;
    para.final_height = 600;
    para.final_width = 1200;
    //para.final_height = 512;
    //para.final_width = 1024;

    para.copy_width = para.final_width;
    para.copy_height = para.merge2_height < para.final_height ? para.merge2_height : para.final_height;

    para.a = 0.0;
    //para.b = -0.2853;
    para.b = 0.0;
    para.c = 0.0;
    para.d = 0.0;
    para.e = 0.0;

    para.cut_width = 0;

    coordinate1 = new int[para.correct_height * para.correct_width];
    memset( coordinate1, 0, sizeof( int ) *para.correct_height * para.correct_width );
    coordinate2 = new int[para.correct_height * para.correct_width];
    memset( coordinate2, 0, sizeof( int ) *para.correct_height * para.correct_width );
    distance1 = new int[para.correct_height * para.correct_width * 2];
    memset( distance1, 0, sizeof( int) *para.correct_height * para.correct_width * 2 );
    distance2 = new int[para.correct_height * para.correct_width * 2];
    memset( distance2, 0, sizeof( int ) *para.correct_height * para.correct_width * 2 );

    ComputeCoordinate( coordinate1, distance1, para.correct_width, para.correct_height, para.imgsrc_W, para.imgsrc_H, para.L_pitchAngle, para.L_rotationAngle, para.L_flipAngle,
                       para.L_center_x, para.L_center_y, para );
    ComputeCoordinate( coordinate2, distance2, para.correct_width, para.correct_height, para.imgsrc_W, para.imgsrc_H, para.R_pitchAngle, para.R_rotationAngle, para.R_flipAngle,
                       para.R_center_x, para.R_center_y, para );


}

void  LG::ComputeCoordinate( int* coordinate, int* distance, int w, int h, int W, int H, float fy, float xz, float fz, int center_x, int center_y, Para& fish_para ) {
	double  PI = 3.14159265359;
	double  sita;
	double  alfa, beta;
	double  x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	double  d1, d2;
	double	cos_fi, sin_fi, x1, y1, r_dst, r_src;
	int     X1, Y1, X2, Y2;     //X��    Y��
	float a, b, c, d, e, r, R;
    int d3,d4;
    int lhbit = 32;
	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	r = fish_para.r;
	R = fish_para.R;

	for ( int i = 0; i < h; i++ ) {
		for ( int j = 0; j < w; j++ ) {
			alfa = PI*( i - r ) / h;
			beta = PI*( j - R ) / h;


			x_ccs = r*sin( alfa );
			y_ccs = r*cos( alfa )*sin( beta );
			z_ccs = r*cos( alfa )*cos( beta );


			x_ccs1 = cos( fy )*x_ccs + sin( fy )*z_ccs;
			y_ccs1 = y_ccs;
			z_ccs1 = -sin( fy )*x_ccs + cos( fy )*z_ccs;


			x_ccs2 = cos( xz )*x_ccs1 + sin( xz )*y_ccs1;
			y_ccs2 = -sin( xz )*x_ccs1 + cos( xz )*y_ccs1;
			z_ccs2 = z_ccs1;


			x_ccs3 = x_ccs2;
			y_ccs3 = cos( fz )*y_ccs2 - sin( fz )*z_ccs2;
			z_ccs3 = sin( fz )*y_ccs2 + cos( fz )*z_ccs2;

			if ( z_ccs3>0 )

				sita = asin( sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 ) / r );

			if ( z_ccs3 == 0 )
				sita = PI / 2;

			if ( z_ccs3 < 0 )
				sita = PI / 2 + asin( -z_ccs3 / r );


			cos_fi = x_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );
			sin_fi = y_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );


			y1 = 2 * r*sita*cos_fi / PI;
			x1 = 2 * r*sita*sin_fi / PI;


			r_dst = sqrt( x1*x1 + y1*y1 ) / R;

			r_src = ( a*r_dst*r_dst*r_dst*r_dst + b*r_dst*r_dst*r_dst + c*r_dst*r_dst + ( 1 - a - b - c )*r_dst )*R;


			x1 = r_src*( x1 / ( R*r_dst ) );
			y1 = r_src*( y1 / ( R*r_dst ) );

			x1 = x1 + center_x + d;
			y1 = y1 + center_y + e;

			X1 = floor( x1 );
			Y1 = floor( y1 );

			d1 = x1 - X1;
			d2 = y1 - Y1;

            d3 = d1 * lhbit;
            d4 = d2 * lhbit;
			if ( X1 >= 0 && Y1 >= 0 && X1 < W - 1 && Y1 < H - 1 ) {
				coordinate[i*w + j] = Y1 + X1*pow_15;
			}
			else {
				coordinate[i*w + j] = H - 1;

			}

			distance[2 * i*w + 2 * j] = d3;
			distance[2 * i*w + 2 * j + 1] = d4;
		}
	}
}

void  LG::ImageCorrect( int* coordinate1, int* distance1, void* data, int W, int H, int channel, Frame& img_cor1, int* coordinate2, int* distance2, Frame& img_cor2, int h, int w ) {

    int x1, y1, x2, y2;
    int s1, s2, s3, s4;
    int id = 0;
    int* ptr_coord1 = NULL;
    int* ptr_coord2 = NULL;
    int* ptr_dis1 = NULL;
    int* ptr_dis2 = NULL;
    unsigned char *ptr_Data1 = NULL;
    unsigned char *ptr_Data2 = NULL;
    unsigned char *src_data = (unsigned char *) data;
    int widthStep = W*channel;
    int QH = 31;

    uchar tl, tr, bl, br;

    for ( int i = 0; i < h; i++ ) {
        ptr_coord1 = coordinate1 + i*w;
        ptr_coord2 = coordinate2 + i*w;
        ptr_dis1 = distance1 + 2 * i*w;
        ptr_dis2 = distance2 + 2 * i*w;
        ptr_Data1 = img_cor1.data + i * img_cor1.widthStep;
        ptr_Data2 = img_cor2.data + i * img_cor2.widthStep;
        for ( int j = 0; j < w; j++ ) {
            x1 = ptr_coord1[0] >> 15;
            y1 = ptr_coord1[0]-(x1<<15);
            y2 = y1 + 1;
            x2 = x1 + 1;

           s1 = ( QH - ptr_dis1[0] )*( QH - ptr_dis1[1] );
           s2 = ( QH - ptr_dis1[0] )*ptr_dis1[1];
            s3 = ptr_dis1[0] * ( QH - ptr_dis1[1] );
            s4 = ptr_dis1[0] * ptr_dis1[1];

            ptr_coord1 += 1;
            ptr_dis1 += 2;

            tl = *( src_data + y1*widthStep + channel * x1 );
            tr = *( src_data + y1*widthStep + channel * x2 );
            bl = *( src_data + y2*widthStep + channel * x1 );
            br = *( src_data + y2*widthStep + channel * x2 );

            ptr_Data1[0] =(uchar)((tl * s1 +tr *s3 + bl * s2 + br* s4)>>10);

            tl = *( src_data + y1*widthStep + 3 * x1 +1);
            tr = *( src_data + y1*widthStep + 3 * x2 +1);
            bl = *( src_data + y2*widthStep + 3 * x1 +1);
            br = *( src_data + y2*widthStep + 3 * x2 +1);

            ptr_Data1[1] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

            tl = *( src_data + y1*widthStep + 3 * x1 + 2 );
            tr = *( src_data + y1*widthStep + 3 * x2 + 2 );
            bl = *( src_data + y2*widthStep + 3 * x1 + 2 );
            br = *( src_data + y2*widthStep + 3 * x2 + 2 );

            ptr_Data1[2] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

            ptr_Data1 += 3;

            x1 = ptr_coord2[0]>>15;
            y1 = ptr_coord2[0]-(x1<<15);

            y2 = y1+1;
            x2 = x1+1;

            s1 = ( QH - ptr_dis2[0] )*( QH - ptr_dis2[1] );
            s2 = ( QH - ptr_dis2[0] )*ptr_dis2[1];
            s3 = ptr_dis2[0] * ( QH - ptr_dis2[1] );
            s4 = ptr_dis2[0] * ptr_dis2[1];

            ptr_coord2 += 1;
            ptr_dis2 += 2;

            tl = *( src_data + y1*widthStep + channel * x1 );
            tr = *( src_data + y1*widthStep + channel * x2 );
            bl = *( src_data + y2*widthStep + channel * x1 );
            br = *( src_data + y2*widthStep + channel * x2 );

            ptr_Data2[0] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

            tl = *( src_data + y1*widthStep + 3 * x1 + 1 );
            tr = *( src_data + y1*widthStep + 3 * x2 + 1 );
            bl = *( src_data + y2*widthStep + 3 * x1 + 1 );
            br = *( src_data + y2*widthStep + 3 * x2 + 1 );

            ptr_Data2[1] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

            tl = *( src_data + y1*widthStep + 3 * x1 + 2 );
            tr = *( src_data + y1*widthStep + 3 * x2 + 2 );
            bl = *( src_data + y2*widthStep + 3 * x1 + 2 );
            br = *( src_data + y2*widthStep + 3 * x2 + 2 );

            ptr_Data2[2] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );
            ptr_Data2 += 3;
        }
    }
}
void  LG::ImageCorrect_one( int* coordinate, int* distance, void* data, int W, int H, int channel, Frame& img_cor, int h, int w ) {

    int x1, y1, x2, y2;
    int s1, s2, s3, s4;
    int id = 0;
    int* ptr_coord = NULL;
    int* ptr_dis= NULL;
    unsigned char *ptr_Data = NULL;
    unsigned char *src_data = (unsigned char *) data;
    int widthStep = W*channel;
    int QH = 31;

    uchar tl, tr, bl, br;

    for ( int i = 0; i < h; i++ ) {
        ptr_coord= coordinate + i*w;
        ptr_dis = distance + 2 * i*w;

        ptr_Data = img_cor.data + i * img_cor.widthStep;
        for ( int j = 0; j < w; j++ ) {
            x1 = ptr_coord[0] >> 15;
            y1 = ptr_coord[0]-(x1<<15);
            y2 = y1 + 1;
            x2 = x1 + 1;

            s1 = ( QH - ptr_dis[0] )*( QH - ptr_dis[1] );
            s2 = ( QH - ptr_dis[0] )*ptr_dis[1];
            s3 = ptr_dis[0] * ( QH - ptr_dis[1] );
            s4 = ptr_dis[0] * ptr_dis[1];

            ptr_coord += 1;
            ptr_dis += 2;

            tl = *( src_data + y1*widthStep + channel * x1 );
            tr = *( src_data + y1*widthStep + channel * x2 );
            bl = *( src_data + y2*widthStep + channel * x1 );
            br = *( src_data + y2*widthStep + channel * x2 );

            ptr_Data[0] =(uchar)((tl * s1 +tr *s3 + bl * s2 + br* s4)>>10);

            tl = *( src_data + y1*widthStep + 3 * x1 +1);
            tr = *( src_data + y1*widthStep + 3 * x2 +1);
            bl = *( src_data + y2*widthStep + 3 * x1 +1);
            br = *( src_data + y2*widthStep + 3 * x2 +1);

            ptr_Data[1] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

            tl = *( src_data + y1*widthStep + 3 * x1 + 2 );
            tr = *( src_data + y1*widthStep + 3 * x2 + 2 );
            bl = *( src_data + y2*widthStep + 3 * x1 + 2 );
            br = *( src_data + y2*widthStep + 3 * x2 + 2 );

            ptr_Data[2] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );
            ptr_Data += 3;

        }
    }
}

void  LG::ImMerge( Frame *img_left, Frame *img_right, int merge_x, int merge_y, Frame* img_dst, int cutwidth ) {

	int wl = img_left->width;
	int hl = img_left->height;
	int wr = img_right->width;
	int hr = img_right->height;
	int cuth = abs( hl - hr );


	int cosswidth = img_left->width - merge_x;
	int h, w;
	h = img_left->height;
	w = img_left->width;

	int cw = cosswidth - 2 * cutwidth;


	int i, j;
	double coff_left, coff_right;

	uchar *ptr_dst = NULL;
	uchar *ptr_src_left = NULL;
	uchar *ptr_src_right = NULL;
	for ( i = 0; i < img_dst->height; ++i ) {

		ptr_dst = ( uchar* ) ( img_dst->data ) + i * img_dst->widthStep;
		ptr_src_left = ( uchar* ) ( img_left->data ) + i * img_left->widthStep;
		for ( j = 0; j < merge_x + cutwidth; ++j ) {
			ptr_dst[0] = ptr_src_left[0];
			ptr_dst[1] = ptr_src_left[1];
			ptr_dst[2] = ptr_src_left[2];
			ptr_dst += 3;
			ptr_src_left += 3;
		}

		ptr_dst = ( uchar* ) ( img_dst->data ) + i * img_dst->widthStep + 3 * ( merge_x + cutwidth );
		ptr_src_right = ( uchar* ) ( img_right->data ) + i * img_right->widthStep + 3 * cutwidth;
		for ( int j = merge_x + cutwidth; j < img_left->width - cutwidth; ++j ) {
			coff_left = ( cw - ( j - merge_x - cutwidth ) ) / ( double ) ( cw );
			coff_right = ( j - merge_x - cutwidth ) / ( double ) ( cw );
			ptr_dst[0] = ( uchar ) ( coff_left * ptr_src_left[0] + coff_right * ptr_src_right[0] );
			ptr_dst[1] = ( uchar ) ( coff_left * ptr_src_left[1] + coff_right * ptr_src_right[1] );
			ptr_dst[2] = ( uchar ) ( coff_left * ptr_src_left[2] + coff_right * ptr_src_right[2] );
			ptr_dst += 3;
			ptr_src_left += 3;
			ptr_src_right += 3;
		}

		for ( int j = img_left->width - cutwidth; j < img_dst->width; ++j ) {
			ptr_dst[0] = ptr_src_right[0];
			ptr_dst[1] = ptr_src_right[1];
			ptr_dst[2] = ptr_src_right[2];

			ptr_dst += 3;
			ptr_src_right += 3;

		}
	}

}

void  LG::IniteFrame( Frame& frame, int w, int h, int channel ) {
	frame.width = w;
	frame.height = h;
	frame.channel = channel;
	frame.widthStep = w*channel;
	frame.data = ( unsigned char* ) malloc( w*h*channel*sizeof( unsigned char ) );
	memset( frame.data, 0, w*h*channel*sizeof( unsigned char ) );
}

void LG::DeleteFrame( Frame& frame ) {

	free( frame.data );
}

void LG::splitFrame( Frame& imframe, Frame& im1, Frame& im2 ) {
	int w;
	int h;
	w = im1.width;
	h = im1.height;
	int i = 0;
	int j = 0;
	for ( i = 0; i < h; i++ )
	for ( j = 0; j < w; j++ ) {
		im1.data[i*im1.widthStep + j*im1.channel] = imframe.data[i*imframe.widthStep + j*imframe.channel];
		im1.data[i*im1.widthStep + j*im1.channel + 1] = imframe.data[i*imframe.widthStep + j*imframe.channel + 1];
		im1.data[i*im1.widthStep + j*im1.channel + 2] = imframe.data[i*imframe.widthStep + j*imframe.channel + 2];
		im2.data[i*im2.widthStep + j*im2.channel] = imframe.data[i*imframe.widthStep + ( j + w )*imframe.channel];
		im2.data[i*im2.widthStep + j*im2.channel + 1] = imframe.data[i*imframe.widthStep + ( j + w )*imframe.channel + 1];
		im2.data[i*im2.widthStep + j*im2.channel + 2] = imframe.data[i*imframe.widthStep + ( j + w )*imframe.channel + 2];
	}
}

void  LG::poly( int* X, int* Y, int len, double* p ) {
    double   sgmaY = 0, sgmaX = 0, sgmaXY = 0, sgmaXX = 0;
    for ( int i = 0; i < len; i++ ) {
        sgmaY += Y[i];
        sgmaX += X[i];
        sgmaXY += Y[i] * X[i];
        sgmaXX += X[i] * X[i];
    }
    p[0] = ( double ) ( ( len*sgmaXY - sgmaX*sgmaY ) ) / ( double ) ( len*sgmaXX - sgmaX*sgmaX );
    p[1] = ( double ) ( sgmaXX*sgmaY - sgmaX*sgmaXY ) / ( double ) ( len*sgmaXX - sgmaX*sgmaX );
}

void  LG::imageMerge_wb( Frame& img_left, Frame& img_right, int merge_x, int merge_y, Frame& img_dst, int cutwidth, int wbw ) {
    clock_t start = clock();
    int wl = img_left.width;
    int hl = img_left.height;
    int wr = img_right.width;
    int hr = img_right.height;
    int channel = img_left.channel;

    int cosswidth = img_left.width - merge_x;  //重合区域宽度 重合区域宽度和裁剪无关
    int h, w;
    h = img_left.height;
    w = img_left.width;
    //砍掉部分重叠区域后的 融合区域宽度
    int cw = cosswidth - 2 * cutwidth;


    const int dennum = 10;  //融合区间差值最大值和最小值之间等分区间数目
    int diff_size = h*cw * 4;
    int* diff = ( int* ) malloc( sizeof( int ) *diff_size );//融合区域差值存储内存   数据的范围[0 255]
    memset( diff, 0, sizeof( int ) *diff_size );

    int denw = ( h*cw + 1 );   //定义每个等分区的宽度  每行的第一个位置存的是该区间点的数目 预留h*cw（融合区最大数目）放点的位置（点的数目不会超过融合区点数目）
    int* Loc = ( int* ) malloc( sizeof( int ) *( h*cw + 1 )*dennum );//每个等分差值区域内的点的位置
    memset( Loc, 0, sizeof( int ) *( h*cw + 1 )*dennum );

    int* dp; //融合区域差值内存当前位置指针
    uchar *ptr_l; //左边融合区域位置指针
    uchar *ptr_r; //右边融合区域位置指针

    int max_diff = 0;
    int min_diff = 10000000; //

    for ( int i = 0; i < h; i++ ) {
        ptr_l = img_left.data  + i*img_left.widthStep  + 3 * ( cutwidth + merge_x );   //指向左图融合区域每行的第一列
        ptr_r = img_right.data + i*img_right.widthStep + 3 * cutwidth;    //指向右图融合区域每行的第一列
        dp = diff + i*cw * 4;  //融合区差值内存每行的第一列
        for ( int j = 0; j < cw; j++ ) {
            dp[0] = abs( ( uchar ) ptr_l[0] - ( uchar ) ptr_r[0] );   //B分量绝对差值
            dp[1] = abs( ( uchar ) ptr_l[1] - ( uchar ) ptr_r[1] );   //G分量绝对差值
            dp[2] = abs( ( uchar ) ptr_l[2] - ( uchar ) ptr_r[2] );   //R分量绝对差值
            dp[3] = dp[0] + dp[1] + dp[2];  //B G R分量绝对差值之和

            max_diff = ( max_diff < dp[3] ? dp[3] : max_diff );    //存放最大值
            min_diff = ( min_diff > dp[3] ? dp[3] : min_diff );    //存放最小值
            /*if(max_diff < dp[3])
            max_diff = dp[3] ;
            if (min_diff > dp[3])
            min_diff = dp[3];*/
            ptr_l += 3;    //左融合区当前位置指针移位
            ptr_r += 3;
            dp += 4;
        }
    }

    clock_t finish1 = clock();
    int dis_diff = max_diff - min_diff; //总体区间宽度
    int qu_block = ( int ) ( 1 / double( dennum )*dis_diff );     //每个区间宽度
    int qu1 = 0;
    int qu2 = 0;
    int temp = 0;
    int id;
    for ( int col = 0; col < cw; col++ )
        for ( int row = 0; row < h; row++ )    //融合区内每个点循环
        {
            temp = diff[row*cw * 4 + 4 * col + 3];
            qu1 = min_diff;
            qu2 = qu1 + qu_block;
            for ( int den = 0; den < dennum; den++ )    //判断点属于哪个区域
            {
                if ( temp >= qu1&&temp < qu2 ) {
                    Loc[den*denw] = Loc[den*denw] + 1;  //每行的第一个位置 存的是属于该区间点的数目
                    id = Loc[den*denw];     //取该区间的数目
                    Loc[den*denw + id] = row*cw + col;    //存融合区坐标  具体的坐标范围为[0  h*cw-1]
                    break;
                }
                if ( den == dennum - 1 ) { qu1 = qu2; qu2 = max_diff + 1; } //最后一个区间范围会大一点
                qu1 = qu2;
                qu2 = qu1 + qu_block;
            }
        }
    //找分区数目最大  第二大 的行号
    int maxnum = 0;
    int maxrow = 0;
    int smaxnum = 0;
    int smaxrow = 0;
    int qunum[10] = { 0 }; //每个区间点的数目
    int sum = 0;
    for ( int i = 0; i < 10; i++ ) {
        qunum[i] = Loc[( i*denw )];
        sum += qunum[i];
        if ( qunum[i]>maxnum ) {
            maxnum = qunum[i];
            maxrow = i;
        }
    }

    //找分区数目第二大的行号
    for ( int i = 0; i < 10; i++ ) {
        temp = Loc[( i*denw )];
        if ( temp <maxnum && temp >smaxnum ) {
            smaxnum = temp;
            smaxrow = i;
        }
    }

    clock_t finish2  = clock();

    int* Y_R = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );
    int* X_R = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );
    int* Y_G = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );
    int* X_G = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );
    int* Y_B = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );
    int* X_B = ( int* ) malloc( sizeof( int ) *( maxnum + smaxnum ) );

    int col, row;
    uchar* ptr1;
    uchar* ptr2;
    for ( int i = 0; i < maxnum; i++ ) {
        temp = Loc[maxrow*denw + i + 1];
        row = temp / cw; //分离出行标
        col = temp % cw; //分离出列标
        ptr1 = img_left.data  + row* img_left.widthStep + channel* ( cutwidth + merge_x + col );
        ptr2 = img_right.data + row* img_right.widthStep + channel * ( cutwidth + col );
        Y_R[i] = ( uchar ) ptr1[0];
        Y_G[i] = ( uchar ) ptr1[1];
        Y_B[i] = ( uchar ) ptr1[2];
        X_R[i] = ( uchar ) ptr2[0];
        X_G[i] = ( uchar ) ptr2[1];
        X_B[i] = ( uchar ) ptr2[2];
    }

    for ( int i = maxnum; i < maxnum + smaxnum; i++ ) {
        temp = Loc[smaxrow*denw + ( i - maxnum ) + 1];
        row = temp / cw;
        col = temp % cw;
        ptr1 = img_left.data + row* img_left.widthStep   + channel * ( cutwidth + merge_x + col );
        ptr2 = img_right.data + row* img_right.widthStep + channel * ( cutwidth + col );
        Y_R[i] = ( uchar ) ptr1[0];
        Y_G[i] = ( uchar ) ptr1[1];
        Y_B[i] = ( uchar ) ptr1[2];
        X_R[i] = ( uchar ) ptr2[0];
        X_G[i] = ( uchar ) ptr2[1];
        X_B[i] = ( uchar ) ptr2[2];
    }

    double  p_r[2];
    double  p_g[2];
    double  p_b[2];
    //两个数组之间的一些运算
    //int n = maxnum + smaxnum;
    int n = maxnum;

    clock_t finish3  = clock();


    //最小二乘法求解
    poly( X_R, Y_R, n, p_r );
    poly( X_G, Y_G, n, p_g );
    poly( X_B, Y_B, n, p_b );


    clock_t finish4  = clock();

    uchar* imr_ptr;
    int line = ( int ) ( 0.5*cosswidth );
    double theta = 0;
    double t = 0;
    for ( int i = 0; i < h; i++ ) {
        imr_ptr = img_right.data + i*img_right.widthStep;
        for ( int j = 0; j < line + 10; j++ ) {
            t = ( int ) ( p_r[0] * ( uchar ) imr_ptr[0] + p_r[1] );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[0] = ( int ) t;
            t = ( int ) ( p_g[0] * ( uchar ) imr_ptr[1] + p_g[1] );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[1] = ( int ) t;
            t = ( int ) ( p_b[0] * ( uchar ) imr_ptr[2] + p_b[1] );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[2] = ( int ) t;
            imr_ptr += 3;
        }
        for ( int j = line + 10; j<w; j++ ) {
            theta = 1 - ( double ) ( j - line - 10 ) / ( double ) ( w - line - 10 );
            t = pow( p_r[0], theta ) *( uchar ) imr_ptr[0] + p_r[1] * ( theta );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[0] = ( int ) t;
            t = pow( p_g[0], theta ) *( uchar ) imr_ptr[1] + p_g[1] * ( theta );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[1] = ( int ) t;
            t = pow( p_b[0], theta ) *( uchar ) imr_ptr[2] + p_b[1] * ( theta );
            t = ( t>255 ? 255 : t );
            t = ( t < 0 ? 0 : t );
            imr_ptr[2] = ( int ) t;
            imr_ptr += 3;
        }
    }

    clock_t finish5  = clock();
    /****************************************************************************/
    /****************************************************************************/

    //填充输出矩阵
    for ( int i = 0; i < img_dst.height; i++ ) {
        //左侧非重合区域直接赋左图值
        for ( int j = 0; j < merge_x + cutwidth; j++ ) {
            img_dst.data[i*img_dst.widthStep + channel * j] = img_left.data[i*img_left.widthStep + channel * j];
            img_dst.data[i*img_dst.widthStep + channel * j + 1] = img_left.data[i*img_left.widthStep + channel * j + 1];
            img_dst.data[i*img_dst.widthStep + channel * j + 2] = img_left.data[i*img_left.widthStep + channel * j + 2];
        }
        for ( int j = merge_x + cutwidth; j < img_left.width - cutwidth; j++ ) {
            img_dst.data[i*img_dst.widthStep + channel * j] = ( int ) ( ( ( cw - ( j - merge_x - cutwidth ) ) / ( double ) ( cw ) )*( uchar ) img_left.data[i*img_left.widthStep + channel * j] + ( j - merge_x - cutwidth ) / ( double ) ( cw ) *( uchar ) img_right.data[i*img_right.widthStep + channel * ( j - merge_x )] );
            img_dst.data[i*img_dst.widthStep + channel * j + 1] = ( int ) ( ( ( cw - ( j - merge_x - cutwidth ) ) / ( double ) ( cw ) )*( uchar ) img_left.data[i*img_left.widthStep + channel * j + 1] + ( j - merge_x - cutwidth ) / ( double ) ( cw ) *( uchar ) img_right.data[i*img_right.widthStep + channel * ( j - merge_x ) + 1] );
            img_dst.data[i*img_dst.widthStep + channel * j + 2] = ( int ) ( ( ( cw - ( j - merge_x - cutwidth ) ) / ( double ) ( cw ) )*( uchar ) img_left.data[i*img_left.widthStep + channel * j + 2] + ( j - merge_x - cutwidth ) / ( double ) ( cw ) *( uchar ) img_right.data[i*img_right.widthStep + channel * ( j - merge_x ) + 2] );
        }
        //右侧非重合区域直接赋右图值
        for ( int j = img_left.width - cutwidth; j < img_dst.width; j++ ) {
            img_dst.data[i*img_dst.widthStep + channel * j] = img_right.data[i*img_right.widthStep + channel * ( j - merge_x )];
            img_dst.data[i*img_dst.widthStep + channel * j + 1] = img_right.data[i*img_right.widthStep + channel * ( j - merge_x ) + 1];
            img_dst.data[i*img_dst.widthStep + channel * j + 2] = img_right.data[i*img_right.widthStep + channel * ( j - merge_x ) + 2];
        }
    }

    clock_t finish6  = clock();
    free( diff );
    free( Loc );
    free( Y_R );
    free( X_R );
    free( Y_G );
    free( X_G );
    free( Y_B );
    free( X_B );

    clock_t finish = clock();
//    LOGW("per frame all: %d", (finish-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 1: %d", (finish1-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 2: %d", (finish2-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 3: %d", (finish3-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 4: %d", (finish4-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 5: %d", (finish5-start)/(CLOCKS_PER_SEC/1000));
//    LOGW("per frame 6: %d", (finish6-start)/(CLOCKS_PER_SEC/1000));

}




















