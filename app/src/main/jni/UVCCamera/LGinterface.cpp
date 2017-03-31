
#include   "LGinterface.h"
//#include   <pthread.h>

void  LG::Panointerface( void* data, int width, int height, int nchannel , int flag_wb) {

	 //clock_t start, finish,start1;
	 //start = clock( );
	 //Correct_blk_L( coordtab1, fishblkL, data, width, height, nchannel, imgcor1);
	 //Correct_blk_R( coordtab2, fishblkR, data, width, height, nchannel, imgcor2 );
     //ImageCorrect_L( coordtab1, data, width, height, nchannel, imgcor1);
	 //ImageCorrect_R( coordtab2, data, width, height, nchannel, imgcor2);

      ImageCorrect_L1( coordtab1, data, width, height, nchannel, imgcor1);
	  ImageCorrect_R1( coordtab2, data, width, height, nchannel, imgcor2);

	//CorrectMergePano(coordtab1, coordtab2, data, width, height, nchannel);
    //finish = clock( );
	// printf( "QH correct time is:%d\n", ( finish - start ) / ( CLOCKS_PER_SEC/1000 ) );

//	  start = clock( );
//	  imcor_wb_para( imgcor1, imgcor2, img1_position[1] - 1, img2_position[0] - 1, pr_r,pr_g,pr_b );
//
//	  imcor_wb( imgcor1, imgcor2, img1_position[1] - 1, img2_position[0] - 1, pr_r, pr_g, pr_b );
//	  finish = clock( );

	 // printf( "wb correct time is:%d\n", ( finish - start ) / ( CLOCKS_PER_SEC / 1000 ) );
	 // imcor_wb_para( imgcor2, imgcor1, img2_position[1] - 1, img1_position[0] - 1, pl_r, pl_g, pl_b );
	  
	  //imcor_wb( imgcor2, imgcor1, img2_position[1] - 1, img1_position[0] - 1, pl_r, pl_g, pl_b );
	
	  //pthread_create( )
	//Merge_pano(&imgcor1, &imgcor2);
	 //start1 = clock();
 	   Merge_pano1(&imgcor1, &imgcor2);
       img_resize_nearst( imgpano, imgdst );
	 //finish = clock();
	// printf("Merge_pano1 time is:%d\n", (finish - start1) / (CLOCKS_PER_SEC/1000));
	  //NoMerge_pano(&imgcor1, &imgcor2);
	 // Blent_Smooth();
	 // Blent_Smooth1();

	/*if ( flag_wb )
		imageMerge_wb( imgcor1, imgcor2, para.merge1_x, para.merge1_y, imMerge1, para.cut_width, 0);
	else
	    ImMerge( &imgcor1, &imgcor2, para.merge1_x, para.merge1_y, &imMerge1, para.cut_width );

	splitFrame( imMerge1, imMer1, imMer2 );

	if ( flag_wb )
		imageMerge_wb( imMer2, imMer1, para.merge2_x, para.merge2_y, imMerge2, para.cut_width,0 );
	else
	   ImMerge( &imMer2, &imMer1, para.merge2_x, para.merge2_y, &imMerge2, para.cut_width );*/

 /*   for(int y = 0; y < imgpano.height; y ++)
    {
        for(int x = 0; x < imgpano.width; x++)
        {
			out_frame->data[y*out_frame->widthStep + x*out_frame->channel] = imgpano.data[y*imgpano.widthStep + x*imgpano.channel];
			out_frame->data[y*out_frame->widthStep + x*out_frame->channel + 1] = imgpano.data[y*imgpano.widthStep + x*imgpano.channel + 1];
			out_frame->data[y*out_frame->widthStep + x*out_frame->channel + 2] = imgpano.data[y*imgpano.widthStep + x*imgpano.channel + 2];
			out_frame->data[y*out_frame->widthStep + x*out_frame->channel + 3] = 255;
        }
    }
*/
    /*for(int k = 0; k < 100; k ++)
    {
        LOGW("pixel value: %d", imMerge2.data[200*imMerge2.width*3+k*3]);
    }*/
}
void  LG::LGInit( int w, int h ,int nchannel,int Blent_hw,int dstw,int dsth) {
    InitImpara();
    para.imgsrc_W = w;   //原始鱼眼图的宽高
	para.imgsrc_H = h;
	channel = nchannel;
	blent_hw = Blent_hw;       //融合区半宽

	//渐入渐出系数表
	cor_jc = new float[blent_hw * 2];
	for (int i = 0; i < blent_hw * 2; i++)  
	{cor_jc[i] = (i*1.0)/(blent_hw * 2.0);}

	//根据PT参数计算
	float diff1,diff2;
	diff1 = (impara[1].imcut.bottom - impara[1].imcut.top) / 2; 	//pt参数计算出中心（正图中心）
	diff2 = (impara[1].imcut.right - impara[1].imcut.left) / 2;
	para.R = diff1;
	para.L_center_x = impara[1].imcut.left + diff2;
	para.L_center_y = impara[1].imcut.top + diff1;
	para.fov = impara[1].fov;
    para.r = floor(para.R * 180 / para.fov);

	para.L_pitchAngle    = impara[1].pitchAngle;
	para.L_rotationAngle = impara[1].rotationAngle;
	para.L_flipAngle     = impara[1].flipAngle;

	diff1 = (impara[2].imcut.bottom - impara[2].imcut.top) / 2;
	diff2 = (impara[2].imcut.right  - impara[2].imcut.left) / 2;
	para.R_center_x = impara[2].imcut.left + para.R;
	para.R_center_y = impara[2].imcut.top + para.R;

	para.R_pitchAngle    = impara[2].pitchAngle ;
	para.R_rotationAngle = impara[2].rotationAngle ;
	para.R_flipAngle     = impara[2].flipAngle;

	flip_angle_correct(para.r, para.L_flipAngle, para.R_flipAngle, para.L_fz_xz, para.R_fz_xz, img1_position, img2_position);  //修正翻转角（确保矫正图移到正确的位置），并得到两图拼缝中线位置

	para.pano_height =  2 * para.r ;
	para.pano_width  =  4 * para.r;

	//桶形畸变参数
	para.a = impara[1].a;
	para.b = impara[1].b;
	para.c = impara[1].c;
	para.d = impara[1].d;
	para.e = impara[1].e;
	para.cut_width = 0;


	//块表相关
	blkh = 32;
	blkw = 96;

	//int corw = ceil( ( img1_position[1] + blent_hw )*1.0 / blkw )*blkw;
	//int corh = ceil( ( img1_position[1] + blent_hw )*1.0 / blkw )*blkw;

    int corw = img1_position[1] + blent_hw ;
    int corh = para.pano_height ;

	IniteFrame( imgcor1, corw, para.pano_height, nchannel );
	IniteFrame( imgcor2, corw, para.pano_height, nchannel );
	IniteFrame(imgpano, para.pano_width, para.pano_height, nchannel);
    IniteFrame( imgdst,  dstw, dsth, nchannel );
	
	IniteCoordTab( coordtab1, corw, para.pano_height, 1 );
	IniteCoordTab( coordtab2, corw, para.pano_height, 1 );
		
	//correct_Lz_coordinate((unsigned int*)coordtab1.data, coordtab1.width, coordtab1.height, (para.imgsrc_W) >> 1, para.imgsrc_H, para);
	//correct_Rz_coordinate((unsigned int*)coordtab2.data, coordtab2.width, coordtab1.height, (para.imgsrc_W) >> 1, para.imgsrc_H, para);  //原图宽高位pt中原图宽高
	correct_coordinate_L((unsigned int*)coordtab1.data, coordtab1.width, coordtab1.height, (para.imgsrc_W) >> 1, para.imgsrc_H, para);
	correct_coordinate_R((unsigned int*)coordtab2.data, coordtab2.width, coordtab1.height, (para.imgsrc_W) >> 1, para.imgsrc_H, para);  //原图宽高位pt中原图宽高

//	maxpicnum1 = 0;
//	maxpicnum2 = 0;
//	blk_w_num = ceil( imgcor1.width*1.0 / blkw );
//	blk_h_num = ceil( imgcor1.height*1.0 / blkh );
//	fishblkL = new FishBlkInfo[blk_w_num*blk_h_num];
//	for ( int i = 0; i < blk_w_num*blk_h_num; i++ )
//	{
//		IniteFishBlkInfo( fishblkL[i] );
//	}
//	fishblkR = new FishBlkInfo[blk_w_num*blk_h_num];
//	for ( int i = 0; i < blk_w_num*blk_h_num; i++ ) {
//		IniteFishBlkInfo( fishblkR[i] );
//	}
//	ComputeBlock(this->fishblkL, this->coordtab1,this->maxpicnum1 );
//	blkdata1 = (unsigned int*)malloc(maxpicnum1*sizeof(unsigned int));
//	memset( blkdata1, 0, maxpicnum1*sizeof( unsigned int ) );
//	ComputeBlock( this->fishblkR, this->coordtab2, this->maxpicnum2 );
//	blkdata2 = ( unsigned int* ) malloc( maxpicnum2*sizeof( unsigned int ) );
//	memset( blkdata2, 0, maxpicnum2*sizeof( unsigned int ) );

	//亮度调整相关	
//	  this->cw = 2 * blent_hw;
//	  this->dennum = 10;  //融合区间差值最大值和最小值之间等分区间数目
//	  this->diff_size = coordtab1.height* cw * 4;
//	  diff = ( int* ) malloc( sizeof( int ) *diff_size );//融合区域差值存储内存   数据的范围[0 255]
//	  memset( diff, 0, sizeof( int ) *diff_size );
//	  denw = ( coordtab1.height*cw + 1 );  //定义每个等分区的宽度  每行的第一个位置存的是该区间点的数目 预留h*cw（融合区最大数目）放点的位置（点的数目不会超过融合区点数目）
//	  Loc = ( int* ) malloc( sizeof( int ) *( h*cw + 1 )*dennum ); //每个等分差值区域内的点的位置
//	  memset( Loc, 0, sizeof( int ) *( h*cw + 1 )*dennum );
    
	//对最小二乘算出的系数量化
	//theta 渐入渐出量化范围
	//量化精度
//	 lh_value = 128;
//	//量化范围 （-0.5 1.5）
//	 si = -1.5 * lh_value;
//	 ei = 1.5 * lh_value;
//	 powh = ei - si;
//	 powblk_h = 400;
//	 powblk_w = 2000;
//    // powblk = (float*)malloc(powblk_h*powblk_w*sizeof(float)); //定义块表的大小
//	 powblk = new int[powblk_h*powblk_w];
//	 for ( int i = 0; i < powblk_h*powblk_w; i++ ) { powblk[i] = 0; }
//	 //memset( powblk, 0, sizeof(float)*powblk_h *powblk_w );
//
//	 int QV = 64;
//	for (int i = si; i<=ei; i++) {
//		for ( int j = 0; j <powblk_w; j++ ) {
//			powblk[(i-si)*powblk_w+j]=pow(i*1.0/lh_value,(j+1)*1.0/powblk_w)*QV;
//		}
//	}
//
//	int line = img1_position[0]-1;
//	int theta_QV = 64;  //theta值量化
//	Theta = new int[coordtab1.width - line];
//	for ( int j = 0; j < coordtab1.width - line; j++ )
//		Theta[j] =(int) ( 1 - ( double ) ( j) / ( double ) ( coordtab1.width - line ) )*theta_QV;
//
//	//float temp1 = pow(-5, 0.5);
//	pr_r[0] = 1.0; pr_r[1] = 0.0;
//	pr_g[0] = 1.0; pr_g[1] = 0.0;
//	pr_b[0] = 1.0; pr_b[1] = 0.0;
//
//	pl_r[0] = 1.0; pl_r[1] = 0.0;
//	pl_g[0] = 1.0; pl_g[1] = 0.0;
//	pl_b[0] = 1.0; pl_b[1] = 0.0;
//
//	Y_R = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
//	X_R = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
//	Y_G = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
//	X_G = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
//	Y_B = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
//	X_B = ( int* ) malloc( sizeof( int ) *( 2000 * 100 ) );
    //最终图像插值相关
    insert_lh=256;
    ratio_h=imgpano.height*1.0*insert_lh/imgdst.height;
    ratio_w = imgpano.width*1.0*insert_lh / imgdst.width;
}
void  LG::ComputeCoordinate(int* coordinate, int* distance, int w, int h, int W, int H, float fy, float xz, float fz, int center_x, int center_y, Para& fish_para ) {
	

	float   sita;
	float   alfa, beta;
	float   x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float   d1, d2;
	float	cos_fi, sin_fi, x1, y1, r_dst, r_src;
	int     X1, Y1, X2, Y2;     //X��    Y��
	float   a, b, c, d, e, r, R;
	int     d3;
	int     d4;

	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	r = fish_para.r;
	R = fish_para.R;

	int lhbit = 32;

	for ( int i = 0; i < h; i++ ) {
		for ( int j = 0; j < w; j++ ) {
			alfa = PI*( i - r ) / h;
			beta = PI*( j - R ) / h;

			x_ccs = r*sin(alfa);
			y_ccs = r*cos(alfa)*sin(beta);
			z_ccs = r*cos(alfa)*cos(beta);

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
void  LG::ComputeCoordinate_L( unsigned int* coordinate, int w, int h, Para& fish_para ) {
	
	float   sita;
	float   alfa, beta;
	float   x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float   d1, d2;
	float	cos_fi, sin_fi, x1, y1, r_dst, r_src;
	int     X1, Y1, X2, Y2;     //X��    Y��
	float   a, b, c, d, e, r, R;
	int     d3;
	int     d4;

	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	r = fish_para.r;
	R = fish_para.R;
	int center_x ;
	int center_y ;
	int W,  H;
	float fy,xz,fz;
	W = fish_para.imgsrc_W;
	H = fish_para.imgsrc_H;
	fy = fish_para.L_pitchAngle;
	xz = fish_para.L_rotationAngle;
	fz = 0;
	center_x = fish_para.L_center_x;
	center_y = fish_para.L_center_y;
	
	int lhbit = 32;
	int max_h = 0;
	int max_w = 0;
	for ( int i = 0; i < h; i++ ) {
		for ( int j = 0; j < w; j++ ) {
			alfa = PI*( i - r + 0.5) / h;
			beta = PI*( j - r + 0.5 ) / h;

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
			/*X1 = round( x1 );
			Y1 = round( y1 );*/

			d1 = x1 - X1;
			d2 = y1 - Y1;

			d3 = d1 * lhbit;
			d4 = d2 * lhbit;

			if ( X1 >= 0 && Y1 >= 0 && X1 < W - 1 && Y1 < H - 1 ) {
				max_h = max_h < Y1 ? Y1 : max_h;
				max_w = max_w < X1 ? X1 : max_w;
				coordinate[i*w + j] = Y1*pow_21+X1*pow_10+d3*pow_5+d4;  // 11 11 5 5           
			}
			else {
				coordinate[i*w + j] =0;
			}
		}
	}
//	printf( "L max_h=%d, max_w=%d\n", max_h, max_w );
}
void  LG::ComputeCoordinate_R( unsigned int* coordinate, int w, int h, Para& fish_para ) {

	float   sita;
	float   alfa, beta;
	float   x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float   d1, d2;
	float	cos_fi, sin_fi, x1, y1, r_dst, r_src;
	int     X1, Y1, X2, Y2;     //X��    Y��
	float   a, b, c, d, e, r, R;
	int     d3;
	int     d4;

	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	r = fish_para.r;
	R = fish_para.R;
	int center_x;
	int center_y;
	int W, H;
	float fy, xz, fz;
	W = fish_para.imgsrc_W;
	H = fish_para.imgsrc_H;
	fy = fish_para.R_pitchAngle;
	xz = fish_para.R_rotationAngle;
	fz =0;
	center_x = fish_para.R_center_x;
	center_y = fish_para.R_center_y;
	int  W_half = W >> 1;

	int lhbit = 32;
	int max_h = 0;
	int max_w = 0;
	for ( int i = 0; i < h; i++ ) {
		for ( int j = 0; j < w; j++ ) {
			alfa = PI*( i - r +0.5) / h;
			beta = PI*( j - r +0.5) / h;

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
		/*	X1 = round( x1 );
			Y1 = round( y1 );*/
			d1 = x1 - X1;
			d2 = y1 - Y1;

			d3 = d1 * lhbit;
			d4 = d2 * lhbit;

			if ( X1 >= 0 && Y1 >= 0 && X1 < W - 1 && Y1 < H - 1 ) {
				max_h = max_h < Y1 ? Y1 : max_h;
				max_w = max_w < X1 ? X1 : max_w;
				coordinate[i*w + j] = Y1*pow_21 + (X1-W_half)*pow_10 + d3*pow_5 + d4;  // 11 11 5 5           
			}
			else {
				coordinate[i*w + j] = 0;
			}
		}
	}
//	printf( "R max_h=%d, max_w=%d\n", max_h, max_w );
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
void  LG::ImageCorrect_L(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor) {

	int x1, y1, x2, y2;
	//int s1, s2, s3, s4;
	int id = 0;
	unsigned int* ptr_coord = NULL;

	unsigned char *ptr_Data = NULL;
	unsigned char *src_data = ( unsigned char * ) data;
	int widthStep = W*channel;
	int QH = 31;

	uchar tl, tr, bl, br;

	unsigned int coord_temp=0;
	int x, y;
	int cor1, cor2;

	for ( int i = 0; i < img_cor.height; i++ ) {
		ptr_coord = (unsigned int*)coordtab.data + i*coordtab.widthStep;
		ptr_Data = img_cor.data + i * img_cor.widthStep;

		for ( int j = 0; j < img_cor.width; j++ ) {
			
			coord_temp = ptr_coord[0];
			//按位提取坐标和系数
			//cor2 = coord_temp&( 0x0000001F );
			//cor1 = (coord_temp&( 0x000003E0 ))>>5;
			x1 = (coord_temp&( 0x001FFC00 ))>>10;
			y1 = (coord_temp&( 0xFFE00000 ))>>21;		
			//y2 = y1 + 1;
			//x2 = x1 + 1;
			//s1 = ( QH - cor1 )*( QH - cor2 );
			//s2 = ( QH - cor1 )*cor2;
			//s3 = cor1 * ( QH - cor2 );
			//s4 = cor1 * cor2;
			ptr_coord += coordtab.channel;
			
			tl = *( src_data + y1*widthStep + channel * x1 );
			//tr = *( src_data + y1*widthStep + channel * x2 );
			//bl = *( src_data + y2*widthStep + channel * x1 );
			//br = *( src_data + y2*widthStep + channel * x2 );
			//ptr_Data[0] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

			ptr_Data[0] = tl;
			tl = *(src_data + y1*widthStep + channel * x1 + 1);
			//tr = *(src_data + y1*widthStep + channel * x2 + 1);
			//bl = *(src_data + y2*widthStep + channel * x1 + 1);
			//br = *(src_data + y2*widthStep + channel * x2 + 1);
			//ptr_Data[1] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

			ptr_Data[1] = tl;
			tl = *(src_data + y1*widthStep + channel * x1 + 2);
			//tr = *(src_data + y1*widthStep + channel * x2 + 2);
			//bl = *(src_data + y2*widthStep + channel * x1 + 2);
			//br = *(src_data + y2*widthStep + channel * x2 + 2);
			//ptr_Data[2] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

			ptr_Data[2] = tl;
			ptr_Data += img_cor.channel;
		}
	}	
}
void  LG::ImageCorrect_R(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor) {

	int x1, y1, x2, y2;
	int s1, s2, s3, s4;
	int id = 0;
	unsigned int* ptr_coord = NULL;

	unsigned char *ptr_Data = NULL;
	unsigned char *src_data = ( unsigned char * ) data;
	int widthStep = W*channel;
	int QH = 31;
	int W_half = (W >>1)-1;

	uchar tl, tr, bl, br;

	unsigned int coord_temp = 0;
	int x, y;
	int cor1, cor2;

	for ( int i = 0; i < img_cor.height; i++ ) {
		ptr_coord = (unsigned int*)coordtab.data + i*coordtab.widthStep;
		ptr_Data = img_cor.data + i * img_cor.widthStep;

		for ( int j = 0; j < img_cor.width; j++ ) {
		
			coord_temp = ptr_coord[0];
			//按位提取坐标和系数
			cor2 = coord_temp&( 0x0000001F );
			cor1 = (coord_temp&( 0x000003E0 )) >> 5;
			x1 = ((coord_temp&( 0x001FFC00 )) >> 10)+W_half;
			y1 = (coord_temp&( 0xFFE00000 )) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = ( QH - cor1 )*( QH - cor2 );
			s2 = ( QH - cor1 )*cor2;
			s3 = cor1 * ( QH - cor2 );
			s4 = cor1 * cor2;
			ptr_coord += coordtab.channel;

			tl = *( src_data + y1*widthStep + channel * x1 );
			tr = *( src_data + y1*widthStep + channel * x2 );
			bl = *( src_data + y2*widthStep + channel * x1 );
			br = *( src_data + y2*widthStep + channel * x2 );
			ptr_Data[0] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );

			//ptr_Data[0] = tl;
			tl = *(src_data + y1*widthStep + channel * x1 + 1);
			tr = *(src_data + y1*widthStep + channel * x2 + 1);
			bl = *(src_data + y2*widthStep + channel * x1 + 1);
			br = *(src_data + y2*widthStep + channel * x2 + 1);
			ptr_Data[1] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );
			//ptr_Data[1] = tl;
			tl = *(src_data + y1*widthStep + channel * x1 + 2);
			tr = *(src_data + y1*widthStep + channel * x2 + 2);
			bl = *(src_data + y2*widthStep + channel * x1 + 2);
			br = *(src_data + y2*widthStep + channel * x2 + 2);
			ptr_Data[2] = ( uchar ) ( ( tl * s1 + tr *s3 + bl * s2 + br* s4 ) >> 10 );
			//ptr_Data[2] = tl;
			ptr_Data += img_cor.channel;

		}
	}
	
}
void  LG::IniteFishBlkInfo( FishBlkInfo & fishblk ) {
	fishblk.height = 0;
	fishblk.width = 0;
	fishblk.jo1 = 0;
	fishblk.jo2 = 0;
	fishblk.pb = 0;
	fishblk.x = 0;
	fishblk.y = 0;
	fishblk.pb = 0;
	fishblk.ptr = 0;
	fishblk.pix_num = 0;
}
void  LG::ImageCorrect_L1( CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor ) {
	    int x1, y1, x2, y2;
		//int s1, s2, s3, s4;
		int id = 0;
		unsigned int*   ptr_coord = NULL;
		unsigned int*   ptr_Data = NULL;
		unsigned char*  src_data = (unsigned char *)data;
		int widthStep = W*channel;
		int QH = 31;
		//unsigned int tl, tr, bl, br;
		unsigned int coord_temp = 0;
		int x, y;
		int cor1, cor2;	
		for (int i = 0; i < img_cor.height; i++) {
			ptr_coord = (unsigned int*)coordtab.data + i*coordtab.widthStep;
			ptr_Data =(unsigned int*)img_cor.data + i * img_cor.width; //此处不能用widthstep
			for (int j = 0; j < img_cor.width; j++) {
				//printf("i=%d,j=%d\n", i, j);
				/*if ( i == 310 && j == 0 ) {				
					int kk = 0;
				}*/
				coord_temp = ptr_coord[0];
				//按位提取坐标和系数
				//cor2 = coord_temp&(0x0000001F);
				//cor1 = (coord_temp&(0x000003E0)) >> 5;
				x1 = (coord_temp&(0x001FFC00)) >> 10;
				y1 = (coord_temp&(0xFFE00000)) >> 21;
				//y2 = y1 + 1;
				//x2 = x1 + 1;
				//s1 = (QH - cor1)*(QH - cor2);
				//s2 = (QH - cor1)*cor2;
				//s3 = cor1 * (QH - cor2);
				//s4 = cor1 * cor2;
				ptr_coord += coordtab.channel;	
				//tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
				//tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
				//bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
				//br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);	
				//要采取饱和操作
				//unsigned int temp1, temp2, temp3;
				/*temp1 = *( src_data + y1*widthStep + channel * x1 );
				temp2 = *( src_data + y1*widthStep + channel * x1 + 1 );
				temp3 = *( src_data + y1*widthStep + channel * x1 + 2 );
				temp1 = tl&0x000000ff;
				temp2 = (tl&0x0000ff00) >> 8;
				temp3 = ( tl&0x00ff0000) >> 16;*/
				//temp1 = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);
				//temp2 = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);
				//temp3 = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10);
				/*unsigned int temp_value  = temp1>255 ? 255 : temp1;
				unsigned int temp_value1 = temp2>255 ? 255 : temp2;
				unsigned int temp_value2 = temp3>255 ? 255 : temp3;
				unsigned int value =  ( ( temp_value+temp_value1<<8+temp_value2<<16 )  );
				memcpy( ptr_Data, (uchar*)&value, 3);*/
				/*ptr_Data[0] = temp_value;
				ptr_Data[1] = temp_value1;
				ptr_Data[2] = temp_value2;*/
				//ptr_Data[0] = temp1;
				//ptr_Data[1] = temp2;
				//ptr_Data[2] = temp3;
				//ptr_Data += img_cor.channel;
			  /*	ptr_Data[0] = *(src_data + y1*widthStep + channel * x1);
				    ptr_Data[1] = *(src_data + y1*widthStep + channel * x1 + 1);
				    ptr_Data[2] = *(src_data + y1*widthStep + channel * x1 + 2);
				    ptr_Data += img_cor.channel;*/
				//unsigned int tempval = 0;
				*ptr_Data = *( unsigned int* ) ( src_data + y1*widthStep + ( x1 << 2 ) );		
				ptr_Data  +=1;	
			}
		}
	
}
void  LG::ImageCorrect_R1(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor) {

	int x1, y1, x2, y2;
	//int s1, s2, s3, s4;
	int id = 0;
	unsigned int* ptr_coord = NULL;


	unsigned int *ptr_Data = NULL;
	unsigned char *src_data = ( unsigned char * ) data;
	int widthStep = W*channel;
	int QH = 31;
	int W_half = (W >>1)-1;
	//unsigned int tl, tr, bl, br;

	unsigned int coord_temp = 0;
	int x, y;
	int cor1, cor2;
	for ( int i = 0; i < img_cor.height; i++ ) {
		ptr_coord = (unsigned int*)coordtab.data + i*coordtab.widthStep;
		ptr_Data =  (unsigned int*)img_cor.data + i * img_cor.width;

		for ( int j = 0; j < img_cor.width; j++ ) {

			coord_temp = ptr_coord[0];
			//按位提取坐标和系数
			//cor2 = coord_temp&( 0x0000001F );
		    //cor1 = ( coord_temp&( 0x000003E0 ) ) >> 5;
			x1 = ( ( coord_temp&( 0x001FFC00 ) ) >> 10 ) + W_half;
			y1 = ( coord_temp&( 0xFFE00000 ) ) >> 21;
			//y2 = y1 + 1;
			//x2 = x1 + 1;
		//	s1 = ( QH - cor1 )*( QH - cor2 );
			//s2 = ( QH - cor1 )*cor2;
			//s3 = cor1 * ( QH - cor2 );
			//s4 = cor1 * cor2;
			ptr_coord += coordtab.channel;

			//tl = *( unsigned int* ) ( src_data + y1*widthStep + channel * x1 );
			//tr = *( unsigned int* ) ( src_data + y1*widthStep + channel * x2 );
			//bl = *( unsigned int* ) ( src_data + y2*widthStep + channel * x1 );
			//br = *( unsigned int* ) ( src_data + y2*widthStep + channel * x2 );

			//要采取饱和操作
			//unsigned int temp1, temp2, temp3;
			/*temp1 = *( src_data + y1*widthStep + channel * x1 );
			temp2 = *( src_data + y1*widthStep + channel * x1 + 1 );
			temp3 = *( src_data + y1*widthStep + channel * x1 + 2 );16;*/
			//temp1 = tl&0x000000ff;
			//temp2 = (tl&0x0000ff00) >> 8;
			//temp3 = (tl & 0x00ff0000) >> 16;
			//temp1 = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);
			//temp2 = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);
			//temp3 = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10);
			/*temp1 = ( ( ( tl & 0x000000ff )*s1 + ( tr&( 0x000000ff ) )*s3 + ( bl&( 0x000000ff ) )*s2 + ( br&( 0x000000ff ) )*s4 ) >> 10 );
			temp2 = ( ( ( ( tl & 0x0000ff00 ) >> 8 )*s1 + ( ( tr & 0x0000ff00 ) >> 8 )*s3 + ( ( bl & 0x0000ff00 ) >> 8 )*s2 + ( ( br & 0x0000ff00 ) >> 8 )*s4 ) >> 10 );
			temp3 = ( ( ( ( tl & 0x00ff0000 ) >> 16 )*s1 + ( ( tr & 0x00ff0000 ) >> 16 )*s3 + ( ( bl & 0x00ff0000 ) >> 16 )*s2 + ( ( br & 0x00ff0000 ) >> 16 )*s4 ) >> 10 );*/
			/*unsigned int temp_value  = temp1>255 ? 255 : temp1;
			unsigned int temp_value1 = temp2>255 ? 255 : temp2;
			unsigned int temp_value2 = temp3>255 ? 255 : temp3;
			unsigned int value =  ( ( temp_value+temp_value1<<8+temp_value2<<16 )  );
			memcpy( ptr_Data, (uchar*)&value, 3);*/
			/*ptr_Data[0] = temp_value;
			ptr_Data[1] = temp_value1;
			ptr_Data[2] = temp_value2;*/
			/*ptr_Data[0] = temp1;
			ptr_Data[1] = temp2;
			ptr_Data[2] = temp3;
			ptr_Data += img_cor.channel;*/

		/*	ptr_Data[0] = *(src_data + y1*widthStep + channel * x1);
			ptr_Data[1] = *(src_data + y1*widthStep + channel * x1+1);
			ptr_Data[2] = *(src_data + y1*widthStep + channel * x1+2);
			ptr_Data += img_cor.channel;*/
			*ptr_Data = *(unsigned int*)( src_data + y1*widthStep + channel * x1 );
			ptr_Data++;



		}
	}

}
void  LG::ImMerge(Frame *img_left, Frame *img_right, int merge_x, int merge_y, Frame* img_dst, int cutwidth ) {

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

		ptr_dst = (uchar*) ( img_dst->data ) + i * img_dst->widthStep + 3 * ( merge_x + cutwidth );
		ptr_src_right = (uchar*) (img_right->data) + i * img_right->widthStep + 3 * cutwidth;
		for ( int j = merge_x + cutwidth; j < img_left->width - cutwidth; ++j ) {
			coff_left = (cw - (j - merge_x - cutwidth ) ) / ( double ) (cw);
			coff_right = ( j - merge_x - cutwidth ) / ( double ) ( cw );
			ptr_dst[0] = ( uchar ) ( coff_left * ptr_src_left[0] + coff_right * ptr_src_right[0] );
			ptr_dst[1] = ( uchar ) ( coff_left * ptr_src_left[1] + coff_right * ptr_src_right[1] );
			ptr_dst[2] = ( uchar ) ( coff_left * ptr_src_left[2] + coff_right * ptr_src_right[2] );
			ptr_dst += 3;
			ptr_src_left  += 3;
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
void  LG::IniteCoordTab(CoordinatTable& coordtab, int w, int h, int channel){
	coordtab.width  = w;
	coordtab.height = h;
	coordtab.channel = channel;
	coordtab.widthStep = w*channel;
	coordtab.data = (void*)malloc(w*h*channel*sizeof(int));
	memset(coordtab.data, 0, w*h*channel*sizeof(int));
}
void  LG::DeleteFrame( Frame& frame ) {

	free( frame.data );
}
void  LG::DeleteCoordTab(CoordinatTable& coordtab) {

	free(coordtab.data);
}
void  LG::splitFrame( Frame& imframe, Frame& im1, Frame& im2 ) {
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
void  LG::poly(int* X, int* Y, int len, double* p) {
	double   sgmaY = 0, sgmaX = 0, sgmaXY = 0, sgmaXX = 0;
	for ( int i = 0; i < len; i++ ) {
		sgmaY += Y[i];
		sgmaX += X[i];
		sgmaXY += Y[i] * X[i];
		sgmaXX += X[i] * X[i];
	}
	p[0] = (double) ( len*sgmaXY - sgmaX*sgmaY )/(double) ( len*sgmaXX - sgmaX*sgmaX );
	p[1] = ( double)( sgmaXX*sgmaY - sgmaX*sgmaXY )/(double)(len*sgmaXX - sgmaX*sgmaX );
}

void  LG::imageMerge_wb( Frame& img_left, Frame& img_right, int merge_x, int merge_y, Frame& img_dst, int cutwidth, int wbw ) {
	
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
	//printf("dis_diff is %d",)
	int dis_diff = max_diff - min_diff; //总体区间宽度     
	int qu_block = (int)(1/double(dennum)*dis_diff);     //每个区间宽度
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
	//最小二乘法求解
	poly( X_R, Y_R, n, p_r );
	poly( X_G, Y_G, n, p_g );
	poly( X_B, Y_B, n, p_b );

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
			t = ( t>255 ? 255 : t);
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

	free( diff );
	free( Loc );
	free( Y_R );
	free( X_R );
	free( Y_G );
	free( X_G );
	free( Y_B );
	free( X_B );
}

void  LG::correct_Rz_coordinate(unsigned int* coordinate, int w, int h, int Width, int Hight, Para& fish_para )
{

	//pt中正的图 向右转的到实际图
	float  sita;
	float  alfa, beta;
	float  x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float  d1, d2;
	float	cos_fi, sin_fi, x, y, x1, y1, x2, y2, x3, y3, r_dst, r_src = 0;
	int    X, Y, X1, Y1, X2, Y2;     //X列    Y行

	float  fy, xz, fz, r, R, a, b, c, d, e, center_x, center_y, r1;
	
	fy = fish_para.R_pitchAngle;
	xz = fish_para.R_rotationAngle;
	fz = fish_para.R_flipAngle;
	fz = fish_para.R_fz_xz;

	fy = fy*PI/180;
	xz = xz*PI/180;
	fz = fz*PI/180;

	r = fish_para.r;
	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	R = fish_para.R;

	int W = Width ;
	int H = Hight ;

	center_x = fish_para.R_center_x;
	center_y = fish_para.R_center_y;
	int num1 = 0;
	int num2 = 0;
	int num3 = 0;
	int lhbit = 32;
	int d3, d4;

	for (int i = 1; i <= h; i++)
	{
		for (int j = 1; j <= w; j++)
		{
			alfa = PI*(i - r - 0.5) / h; //纬度
			beta = PI*(j - 2 * r - 0.5) / h; //经度

			//经纬度转化为球面坐标
			x_ccs = r*sin(alfa);
			y_ccs = r*cos(alfa)*sin(beta);
			z_ccs = r*cos(alfa)*cos(beta);

			//翻转角修正
			x_ccs1 = x_ccs;
			y_ccs1 = cos(fz)*y_ccs - sin(fz)*z_ccs;
			z_ccs1 = sin(fz)*y_ccs + cos(fz)*z_ccs;

			//旋转角修正
			x_ccs2 = cos(xz)*x_ccs1 - sin(xz)*y_ccs1;
			y_ccs2 = sin(xz)*x_ccs1 + cos(xz)*y_ccs1;
			z_ccs2 = z_ccs1;

			//俯仰角修正
			x_ccs3 = cos(fy)*x_ccs2 + sin(fy)*z_ccs2;
			y_ccs3 = y_ccs2;
			z_ccs3 = cos(fy)*z_ccs2 - sin(fy)*x_ccs2;

			if (z_ccs3>0)
				// 球面坐标点与Z轴夹角sita
				sita = PI / 2 - asin(z_ccs3 / r);

			if (z_ccs3 == 0)
				sita = PI / 2;

			if (z_ccs3 < 0)
				sita = PI / 2 + asin(-z_ccs3 / r);

			//球面坐标点到XOY平面的垂直投影与X轴夹角fi
			cos_fi = x_ccs3 / sqrt(x_ccs3*x_ccs3 + y_ccs3*y_ccs3);
			sin_fi = y_ccs3 / sqrt(x_ccs3*x_ccs3 + y_ccs3*y_ccs3);

			//球面坐标转化到图像坐标
			y1 = 2 * r*sita*cos_fi / PI; //行标
			x1 = 2 * r*sita*sin_fi / PI; //列标

			//筒形畸变校正  参数待确认
			r_dst = sqrt(x1*x1 + y1*y1) / R;

			if (r_dst < 1.0)
			{

				r_src = (a*r_dst*r_dst*r_dst*r_dst + b*r_dst*r_dst*r_dst + c*r_dst*r_dst + (1 - a - b - c)*r_dst)*R;

				x2 = r_src*sin_fi;  //列标
				y2 = r_src*cos_fi;  //行标

				r1 = sqrt(x2*x2 + y2*y2);

				x3 = x2 + center_x + d;  //对应原鱼眼图的列
				y3 = y2 + center_y + e;

				/*X = floor(x1);
				Y = floor(y1);*/

				X = floor(W - y3);
				Y = floor(x3);

				d1 = W - y3 - X;  //列系数
				d2 = x3 - Y;    //行系数

				d3 = d1 * lhbit;
				d4 = d2 * lhbit;

				if (X > 0 && Y > 0 && X < W  && Y < H && r1<R)
				{
					coordinate[(i-1)*w + j-1] = (Y-1)*pow_21 + (X-1)*pow_10 + d3*pow_5 + d4;  // 11 11 5 5   
					num1++;
				}
				else
				{
					coordinate[(i - 1)*w + j - 1] = 0;       //存储对应原鱼眼图的  行标				
					num2++;
				}

			}

			else{				
				num3++;
			}

		}
	}

	 printf("%d  %d  %d\n", num1, num2, num3); //测试各环节点的数目
}
void  LG::correct_coordinate_R( unsigned int* coordinate, int w, int h, int Width, int Hight, Para& fish_para ) {

	//pt中正的图 向右转的到实际图
	float  sita;
	float  alfa, beta;
	float  x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float  d1, d2;
	float	cos_fi, sin_fi, x, y, x1, y1, x2, y2, x3, y3, r_dst, r_src = 0;
	int    X, Y, X1, Y1, X2, Y2;     //X列    Y行

	float  fy, xz, fz, r, R, a, b, c, d, e, center_x, center_y, r1;

	fy = fish_para.R_pitchAngle;
	xz = fish_para.R_rotationAngle;
	fz = fish_para.R_flipAngle;
	fz = fish_para.R_fz_xz;

	fy = fy*PI / 180;
	xz = xz*PI / 180;
	fz = fz*PI / 180;

	r = fish_para.r;
	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	R = fish_para.R;

	int W = Width;
	int H = Hight;

	center_x = fish_para.R_center_x;
	center_y = fish_para.R_center_y;
	int num1 = 0;
	int num2 = 0;
	int num3 = 0;
	int lhbit = 32;
	int d3, d4;

	for ( int i = 1; i <= h; i++ ) {
		for ( int j = 1; j <= w; j++ ) {
			alfa = PI*( i - r - 0.5 ) / h; //纬度
			beta = PI*( j - 2 * r - 0.5 ) / h; //经度

			//经纬度转化为球面坐标
			x_ccs = r*sin( alfa );
			y_ccs = r*cos( alfa )*sin( beta );
			z_ccs = r*cos( alfa )*cos( beta );

			//翻转角修正
			x_ccs1 = x_ccs;
			y_ccs1 = cos( fz )*y_ccs - sin( fz )*z_ccs;
			z_ccs1 = sin( fz )*y_ccs + cos( fz )*z_ccs;

			//旋转角修正
			x_ccs2 = cos( xz )*x_ccs1 - sin( xz )*y_ccs1;
			y_ccs2 = sin( xz )*x_ccs1 + cos( xz )*y_ccs1;
			z_ccs2 = z_ccs1;

			//俯仰角修正
			x_ccs3 = cos( fy )*x_ccs2 + sin( fy )*z_ccs2;
			y_ccs3 = y_ccs2;
			z_ccs3 = cos( fy )*z_ccs2 - sin( fy )*x_ccs2;

			if ( z_ccs3>0 )
				// 球面坐标点与Z轴夹角sita
				sita = PI / 2 - asin( z_ccs3 / r );

			if ( z_ccs3 == 0 )
				sita = PI / 2;

			if ( z_ccs3 < 0 )
				sita = PI / 2 + asin( -z_ccs3 / r );

			//球面坐标点到XOY平面的垂直投影与X轴夹角fi
			cos_fi = x_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );
			sin_fi = y_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );

			//球面坐标转化到图像坐标
			y1 = 2 * r*sita*cos_fi / PI; //行标
			x1 = 2 * r*sita*sin_fi / PI; //列标

			//筒形畸变校正  参数待确认
			r_dst = sqrt( x1*x1 + y1*y1 ) / R;

			if ( r_dst < 1.0 ) {

				r_src = ( a*r_dst*r_dst*r_dst*r_dst + b*r_dst*r_dst*r_dst + c*r_dst*r_dst + ( 1 - a - b - c )*r_dst )*R;

				x2 = r_src*sin_fi;  //列标
				y2 = r_src*cos_fi;  //行标

				r1 = sqrt( x2*x2 + y2*y2 );

				x3 = x2 + center_x + d;  //对应原鱼眼图的列
				y3 = y2 + center_y + e;

				X = floor(x3);
				Y = floor(y3);

				/*X = floor( W - y3 );
				Y = floor( x3 );*/

				//d1 = W - y3 - X;  //列系数
				//d2 = x3 - Y;    //行系数
				d1 = x3 - X;
				d2 = y3 - Y;

				d3 = d1 * lhbit;
				d4 = d2 * lhbit;

				if ( X > 0 && Y > 0 && X < W  && Y < H && r1<R ) {
					coordinate[( i - 1 )*w + j - 1] = ( Y - 1 )*pow_21 + ( X - 1 )*pow_10 + d3*pow_5 + d4;  // 11 11 5 5   
					num1++;
				}
				else {
					coordinate[( i - 1 )*w + j - 1] = 0;       //存储对应原鱼眼图的  行标				
					num2++;
				}

			}

			else {
				num3++;
			}

		}
	}

	
}

void  LG::correct_Lz_coordinate(unsigned int* coordinate, int w, int h, int W, int H, Para& fish_para)
{
	float   sita;
	float   alfa, beta;
	float   x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float   d1, d2;
	float	cos_fi, sin_fi, x, y, x1, y1, x2, y2, x3, y3, r_dst, r_src = 0;
	int     X, Y, X1, Y1, X2, Y2;     //X列    Y行

	float  fy, xz, fz, r, R, a, b, c, d, e, center_x, center_y, r1;

	fy = fish_para.L_pitchAngle;
	xz = fish_para.L_rotationAngle;
	fz = fish_para.L_flipAngle;
	fz = fish_para.L_fz_xz;

	fy = fy*PI / 180;
	xz = xz*PI / 180;
	fz = fz*PI / 180;

	r = fish_para.r;
	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	R = fish_para.R;

	center_x = fish_para.L_center_x;
	center_y = fish_para.L_center_y;
	int num1 = 0;
	int num2 = 0;
	int num3 = 0;
	int lhbit = 32;
	int d3, d4;
	for (int i = 1; i <= h; i++)
	{
		for (int j = 1; j <= w; j++)
		{
			alfa = PI*(i - r - 0.5) / h; //纬度
			beta = PI*(j - 2 * r - 0.5) / h; //经度

			//经纬度转化为球面坐标
			x_ccs = r*sin(alfa);
			y_ccs = r*cos(alfa)*sin(beta);
			z_ccs = r*cos(alfa)*cos(beta);

			//翻转角修正
			x_ccs1 = x_ccs;
			y_ccs1 = cos(fz)*y_ccs - sin(fz)*z_ccs;
			z_ccs1 = sin(fz)*y_ccs + cos(fz)*z_ccs;

			//旋转角修正
			x_ccs2 = cos(xz)*x_ccs1 - sin(xz)*y_ccs1;
			y_ccs2 = sin(xz)*x_ccs1 + cos(xz)*y_ccs1;
			z_ccs2 = z_ccs1;

			//俯仰角修正
			x_ccs3 = cos(fy)*x_ccs2 + sin(fy)*z_ccs2;
			y_ccs3 = y_ccs2;
			z_ccs3 = cos(fy)*z_ccs2 - sin(fy)*x_ccs2;

			if (z_ccs3>0)
				// 球面坐标点与Z轴夹角sita
				sita = PI / 2 - asin(z_ccs3 / r);

			if (z_ccs3 == 0)
				sita = PI / 2;

			if (z_ccs3 < 0)
				sita = PI / 2 + asin(-z_ccs3 / r);

			//球面坐标点到XOY平面的垂直投影与X轴夹角fi
			cos_fi = x_ccs3 / sqrt(x_ccs3*x_ccs3 + y_ccs3*y_ccs3);
			sin_fi = y_ccs3 / sqrt(x_ccs3*x_ccs3 + y_ccs3*y_ccs3);

			//球面坐标转化到图像坐标
			y1 = 2 * r*sita*cos_fi / PI; //行标
			x1 = 2 * r*sita*sin_fi / PI; //列标

			//筒形畸变校正  参数待确认
			r_dst = sqrt(x1*x1 + y1*y1) / R;

			if (r_dst < 1.0)
			{
				r_src = (a*r_dst*r_dst*r_dst*r_dst + b*r_dst*r_dst*r_dst + c*r_dst*r_dst + (1 - a - b - c)*r_dst)*R;

				x2 = r_src*sin_fi;  //列标
				y2 = r_src*cos_fi;  //行标

				r1 = sqrt(x2*x2 + y2*y2);

				x3 = x2 + center_x + d;  //对应原鱼眼图的列
				y3 = y2 + center_y + e;

				/*X = floor(x1);
				Y = floor(y1);*/

				X = floor(y3);
				Y = floor(H - x3);

				d1 = y3 - X;  //列系数
				d2 = H - x3 - Y;    //行系数

				d3 = d1 * lhbit;
				d4 = d2 * lhbit;

				if (X > 0 && Y > 0 && X < W  && Y < H && r1<R)
				{
					coordinate[(i-1)*w + j-1] = (Y-1)*pow_21 + (X-1)*pow_10 + d3*pow_5 + d4;
					num1++;
				}
				else
				{
					coordinate[(i - 1)*w + j - 1] = 0;
					num2++;
				}
			}

			else{
				
				num3++;
			}
		}
	}
	// printf("%d  %d  %d\n", num1, num2, num3); 测试各环节点的数目
}
void  LG::correct_coordinate_L( unsigned int* coordinate, int w, int h, int W, int H, Para& fish_para ) {
	float   sita;
	float   alfa, beta;
	float   x_ccs, y_ccs, z_ccs, x_ccs1, y_ccs1, z_ccs1, x_ccs2, y_ccs2, z_ccs2, x_ccs3, y_ccs3, z_ccs3;
	float   d1, d2;
	float	cos_fi, sin_fi, x, y, x1, y1, x2, y2, x3, y3, r_dst, r_src = 0;
	int     X, Y, X1, Y1, X2, Y2;     //X列    Y行

	float  fy, xz, fz, r, R, a, b, c, d, e, center_x, center_y, r1;

	fy = fish_para.L_pitchAngle;
	xz = fish_para.L_rotationAngle;
	fz = fish_para.L_flipAngle;
	fz = fish_para.L_fz_xz;

	fy = fy*PI / 180;
	xz = xz*PI / 180;
	fz = fz*PI / 180;

	r = fish_para.r;
	a = fish_para.a;
	b = fish_para.b;
	c = fish_para.c;
	d = fish_para.d;
	e = fish_para.e;
	R = fish_para.R;

	center_x = fish_para.L_center_x;
	center_y = fish_para.L_center_y;
	int num1 = 0;
	int num2 = 0;
	int num3 = 0;
	int lhbit = 32;
	int d3, d4;
	for ( int i = 1; i <= h; i++ ) {
		for ( int j = 1; j <= w; j++ ) {
			alfa = PI*( i - r - 0.5 ) / h; //纬度
			beta = PI*( j - 2 * r - 0.5 ) / h; //经度

			//经纬度转化为球面坐标
			x_ccs = r*sin( alfa );
			y_ccs = r*cos( alfa )*sin( beta );
			z_ccs = r*cos( alfa )*cos( beta );

			//翻转角修正
			x_ccs1 = x_ccs;
			y_ccs1 = cos( fz )*y_ccs - sin( fz )*z_ccs;
			z_ccs1 = sin( fz )*y_ccs + cos( fz )*z_ccs;

			//旋转角修正
			x_ccs2 = cos( xz )*x_ccs1 - sin( xz )*y_ccs1;
			y_ccs2 = sin( xz )*x_ccs1 + cos( xz )*y_ccs1;
			z_ccs2 = z_ccs1;

			//俯仰角修正
			x_ccs3 = cos( fy )*x_ccs2 + sin( fy )*z_ccs2;
			y_ccs3 = y_ccs2;
			z_ccs3 = cos( fy )*z_ccs2 - sin( fy )*x_ccs2;

			if ( z_ccs3>0 )
				// 球面坐标点与Z轴夹角sita
				sita = PI / 2 - asin( z_ccs3 / r );

			if ( z_ccs3 == 0 )
				sita = PI / 2;

			if ( z_ccs3 < 0 )
				sita = PI / 2 + asin( -z_ccs3 / r );

			//球面坐标点到XOY平面的垂直投影与X轴夹角fi
			cos_fi = x_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );
			sin_fi = y_ccs3 / sqrt( x_ccs3*x_ccs3 + y_ccs3*y_ccs3 );

			//球面坐标转化到图像坐标
			y1 = 2 * r*sita*cos_fi / PI; //行标
			x1 = 2 * r*sita*sin_fi / PI; //列标

			//筒形畸变校正  参数待确认
			r_dst = sqrt( x1*x1 + y1*y1 ) / R;

			if ( r_dst < 1.0 ) {
				r_src = ( a*r_dst*r_dst*r_dst*r_dst + b*r_dst*r_dst*r_dst + c*r_dst*r_dst + ( 1 - a - b - c )*r_dst )*R;

				x2 = r_src*sin_fi;  //列标
				y2 = r_src*cos_fi;  //行标

				r1 = sqrt( x2*x2 + y2*y2 );

				x3 = x2 + center_x + e;  //对应原鱼眼图的列
				y3 = y2 + center_y + d;

				/*X = floor(x1);
				Y = floor(y1);*/

				X = floor( x3 );
				Y = floor( y3 );

				d1 = x3 - X;  //列系数
				d2 = y3- Y;    //行系数

				d3 = d1 * lhbit;
				d4 = d2 * lhbit;

				if ( X > 0 && Y > 0 && X < W  && Y < H && r1<R ) {
					coordinate[( i - 1 )*w + j - 1] = ( Y - 1 )*pow_21 + ( X - 1 )*pow_10 + d3*pow_5 + d4;
					num1++;
				}
				else {
					coordinate[( i - 1 )*w + j - 1] = 0;
					num2++;
				}
			}

			else {

				num3++;
			}
		}
	}
	// printf("%d  %d  %d\n", num1, num2, num3); 测试各环节点的数目
}


void  LG::computeimpara(ImPara& impara){

	float diff1;
	float diff2;
	diff1 = impara.imcut.bottom - impara.imcut.top;
	diff2 = impara.imcut.right - impara.imcut.left;
	impara.R = ((abs(diff1)>abs(diff2)) ? abs(diff1) : abs(diff2)) / 2.0;
	impara.r = round(impara.R * 180 / impara.fov);
	impara.centerx = impara.imcut.left + impara.R;
	impara.centery = impara.imcut.top + impara.R;

}

int   LG::MOD(int x, int y){
	int n;
	n = floor(x*1.0 / (1.0*y));
	n = x - n*y;
	return n;
}
float LG::MODf(float x, float y)
{
	int n;
	float result;
	n = floor(x*1.0 / (1.0*y));
	result = x - n*y;
	return result;
}

void  LG::flip_angle_correct(int r, float& fanzhuan1, float& fanzhuan2, float& fanzhuan1_xz, float& fanzhuan2_xz, int img1_position[2], int img2_position[2])
{
	
	int w = 4 * r;
	int img1_out_right, img1_out_left, img2_out_left, img2_out_right;
	int shit1, shit2;
		
	float AngleDiff = (fanzhuan1 + fanzhuan2) / 2;
	//判读中心角度是否在 角度1和角度2之间 否则中心角度存在循环移位
	if ((AngleDiff > fanzhuan1) && (AngleDiff < fanzhuan2))
	{
		img1_out_right = w / 2 + round(w*AngleDiff / 360);
		img1_out_left = MOD(img1_out_right + w / 2, w);
		img2_out_right = img1_out_left;
		img2_out_left = img1_out_right;
	}
	
	else{
		img1_out_left = w / 2 + round(w*AngleDiff / 360);
		img1_out_right = MOD(img1_out_left + w / 2, w);
		img2_out_right = img1_out_left;
		img2_out_left = img1_out_right;
	}

	//有问题 360.0
    fanzhuan1_xz = -(360.0 * (img1_out_left - blent_hw - 1) / w - fanzhuan1);
	fanzhuan1_xz = MODf(fanzhuan1_xz + 180, 360) - 180;  //modf
	fanzhuan2_xz = -(360.0 * (img2_out_left - blent_hw - 1) / w - fanzhuan2);
    fanzhuan2_xz = MODf(fanzhuan2_xz + 180, 360) - 180;

	shit1 = round((w * (fanzhuan1 - fanzhuan1_xz)) / 360);
	shit2 = round((w * (fanzhuan2 - fanzhuan2_xz)) / 360);

	img1_out_left  = MOD(img1_out_left   - shit1, w);
	img1_out_right = MOD(img1_out_right - shit1, w);


	img2_out_left  = MOD(img2_out_left  - shit2, w);
	img2_out_right = MOD(img2_out_right - shit2, w);

	img1_position[0] = img1_out_left;
	img1_position[1] = img1_out_right;
	img2_position[0] = img2_out_left;
	img2_position[1] = img2_out_right;
}

void  LG::Merge_pano(Frame *img_left, Frame *img_right){

	int panow = imgpano.width;
	int panoh = imgpano.height;
	int cnl = img_left->channel;
	int cnr = img_right->channel;
	int cnd = imgpano.channel;
	int position1 = img1_position[0];
	int position2 = img1_position[1];
	uchar* dst_ptr = NULL;
	uchar* iml_ptr = NULL;
	uchar* imr_ptr = NULL;
	float* cor_jc= new float[blent_hw * 2];
	for (int i = 0; i < blent_hw * 2; i++)    {	cor_jc[i] = (i*1.0) / (blent_hw * 2.0);}
	uchar tempr_l, tempg_l, tempb_l,tempr_r,tempg_r,tempb_r;
	int mean_l, mean_r;
	for (int i = 0; i < panoh; i++){
		dst_ptr = imgpano.data + i*imgpano.widthStep;
		iml_ptr = img_left->data + i*img_left->widthStep;
		imr_ptr = img_right->data + i*img_right->widthStep;

		for (int j = 0; j < blent_hw * 2; j++)   //第一个融合区
		{
			tempr_l = iml_ptr[j * cnl + 0];
			tempg_l = iml_ptr[j * cnl + 1];
			tempb_l = iml_ptr[j * cnl + 2];
			tempr_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 0];
			tempg_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 1];
			tempb_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 2];
			mean_l = (tempr_l + tempg_l + tempb_l) / 3;
			mean_r = (tempr_r + tempg_l + tempb_l) / 3;

			dst_ptr[j * cnd + 0] = (uchar)(cor_jc[j] * tempr_l + (1 - cor_jc[j])*tempr_r);
			dst_ptr[j * cnd + 1] = (uchar)(cor_jc[j] * tempg_l + (1 - cor_jc[j])*tempg_r);
			dst_ptr[j * cnd + 2] = (uchar)(cor_jc[j] * tempb_l + (1 - cor_jc[j])*tempb_r);
			
			if (abs(mean_l - mean_r) > 0.9*(mean_l>mean_r?mean_l:mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_l ? tempr_l : tempr_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? tempg_l : tempg_r);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? tempb_l : tempb_r);
			}
			/*dst_ptr[j * 3 + 0] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 0] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 0]);
			dst_ptr[j * 3 + 1] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 1] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 2] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 2]);*/
		
		}
		for (int j = blent_hw * 2; j < position2 - blent_hw; j++){	
			dst_ptr[j * cnd + 0] = iml_ptr[j * cnl + 0];
			dst_ptr[j * cnd + 1] = iml_ptr[j * cnl + 1];
			dst_ptr[j * cnd + 2] = iml_ptr[j * cnl + 2];	
		}

		for (int j = position2 - blent_hw; j < position2 + blent_hw; j++)  //第二个融合区
		{

			tempr_l = iml_ptr[j * cnl + 0];
			tempg_l = iml_ptr[j * cnl + 1];
			tempb_l = iml_ptr[j * cnl + 2];
			tempr_r = imr_ptr[(j - (position2 - blent_hw)) * cnr];
			tempg_r = imr_ptr[(j - (position2 - blent_hw)) * cnr + 1];
			tempb_r = imr_ptr[(j - (position2 - blent_hw)) * cnr + 2];

			dst_ptr[j * cnd + 0] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempr_l + cor_jc[j - (position2 - blent_hw)] * tempr_r);
			dst_ptr[j * cnd + 1] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempg_l + cor_jc[j - (position2 - blent_hw)] * tempg_r);
			dst_ptr[j * cnd + 2] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempb_l + cor_jc[j - (position2 - blent_hw)] * tempb_r);
			mean_l = (tempr_l + tempg_l + tempb_l) / 3;
			mean_r = (tempr_r + tempg_l + tempb_l) / 3;
			if (abs(mean_l - mean_r) > 0.9*(mean_l>mean_r ? mean_l : mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_r ? tempr_l : tempr_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? tempg_l : tempg_r);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? tempb_l : tempb_r);
			}
			/*dst_ptr[j * 3 + 0] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 0] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw))*3]);
			dst_ptr[j * 3 + 1] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 1] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw)) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 2] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw)) * 3 + 2]);	*/	
		}

		for (int j = position2 + blent_hw; j < panow; j++){	
			dst_ptr[j * cnd + 0] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr];
			dst_ptr[j * cnd + 1] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 1];
			dst_ptr[j * cnd + 2] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 2];
		}
	}

	/*IplImage* impano = cvCreateImageHeader({ panow, panoh }, IPL_DEPTH_8U, img_dst->channel);
	cvSetData(impano, img_dst->data, img_dst->widthStep);
	cvSaveImage("./src_pic/pano.jpg", impano);
	cvReleaseImageHeader(&impano);*/

	delete cor_jc;

}

void  LG::Merge_pano1(Frame *img_left, Frame *img_right){

	int panow = imgpano.width;
	int panoh = imgpano.height;
	int cnl = img_left->channel;
	int cnr = img_right->channel;
	int cnd = imgpano.channel;
	int position1 = img1_position[0];
	int position2 = img1_position[1];
	uchar* dst_ptr = NULL;
	uchar* iml_ptr = NULL;
	uchar* imr_ptr = NULL;
	
	uchar tempr_l, tempg_l, tempb_l, tempr_r, tempg_r, tempb_r;
	int mean_l, mean_r;
	for (int i = 0; i < panoh; i++){
		dst_ptr = imgpano.data + i*imgpano.widthStep;
		iml_ptr = img_left->data + i*img_left->widthStep;
		imr_ptr = img_right->data + i*img_right->widthStep;

		for (int j = 0; j < blent_hw * 2; j++)   //第一个融合区
		{
			tempr_l = iml_ptr[j * cnl + 0];
			tempg_l = iml_ptr[j * cnl + 1];
			tempb_l = iml_ptr[j * cnl + 2];
			tempr_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 0];
			tempg_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 1];
			tempb_r = imr_ptr[(img2_position[1] - blent_hw + j) * cnr + 2];
			//mean_l = (tempr_l + tempg_l + tempb_l) / 3;
			//mean_r = (tempr_r + tempg_l + tempb_l) / 3;

			dst_ptr[j * cnd + 0] = (uchar)(cor_jc[j] * tempr_l + (1 - cor_jc[j])*tempr_r);
			dst_ptr[j * cnd + 1] = (uchar)(cor_jc[j] * tempg_l + (1 - cor_jc[j])*tempg_r);
			dst_ptr[j * cnd + 2] = (uchar)(cor_jc[j] * tempb_l + (1 - cor_jc[j])*tempb_r);

			/*if (abs(mean_l - mean_r) > 0.9*(mean_l>mean_r ? mean_l : mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_l ? tempr_l : tempr_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? tempg_l : tempg_r);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? tempb_l : tempb_r);
			}*/
			/*dst_ptr[j * 3 + 0] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 0] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 0]);
			dst_ptr[j * 3 + 1] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 1] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 2] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 2]);*/

		}

		/*for (int j = blent_hw * 2; j < position2 - blent_hw; j++){
			dst_ptr[j * cnd + 0] = iml_ptr[j * cnl + 0];
			dst_ptr[j * cnd + 1] = iml_ptr[j * cnl + 1];
			dst_ptr[j * cnd + 2] = iml_ptr[j * cnl + 2];
		}*/
		memcpy(dst_ptr + blent_hw * 2 * cnd, iml_ptr + blent_hw * 2 * cnl, sizeof(uchar)*channel*(position2 - 3 * blent_hw));//非重叠区域

		for (int j = position2 - blent_hw; j < position2 + blent_hw; j++)  //第二个融合区
		{

			tempr_l = iml_ptr[j * cnl + 0];
			tempg_l = iml_ptr[j * cnl + 1];
			tempb_l = iml_ptr[j * cnl + 2];
			tempr_r = imr_ptr[(j - (position2 - blent_hw)) * cnr];
			tempg_r = imr_ptr[(j - (position2 - blent_hw)) * cnr + 1];
			tempb_r = imr_ptr[(j - (position2 - blent_hw)) * cnr + 2];

			dst_ptr[j * cnd + 0] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempr_l + cor_jc[j - (position2 - blent_hw)] * tempr_r);
			dst_ptr[j * cnd + 1] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempg_l + cor_jc[j - (position2 - blent_hw)] * tempg_r);
			dst_ptr[j * cnd + 2] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * tempb_l + cor_jc[j - (position2 - blent_hw)] * tempb_r);
			/*mean_l = (tempr_l + tempg_l + tempb_l) / 3;
			mean_r = (tempr_r + tempg_l + tempb_l) / 3;
			if (abs(mean_l - mean_r) > 0.9*(mean_l>mean_r ? mean_l : mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_r ? tempr_l : tempr_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? tempg_l : tempg_r);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? tempb_l : tempb_r);
			}*/
			/*dst_ptr[j * 3 + 0] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 0] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw))*3]);
			dst_ptr[j * 3 + 1] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 1] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw)) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)])*iml_ptr[j * 3 + 2] + cor_jc[j - (position2 - blent_hw)] * imr_ptr[(j - (position2 - blent_hw)) * 3 + 2]);	*/
		}

		/*for (int j = position2 + blent_hw; j < panow; j++){
			dst_ptr[j * cnd + 0] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr];
			dst_ptr[j * cnd + 1] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 1];
			dst_ptr[j * cnd + 2] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 2];
		}*/
		//非重叠区域进行内存拷贝
		 memcpy(dst_ptr + (position2 + blent_hw)*cnd, imr_ptr + blent_hw * 2 * cnr, sizeof(uchar)*channel*(panow-position2 - blent_hw));
	}

	/*IplImage* impano = cvCreateImageHeader({ panow, panoh }, IPL_DEPTH_8U, img_dst->channel);
	cvSetData(impano, img_dst->data, img_dst->widthStep);
	cvSaveImage("./src_pic/pano.jpg", impano);
	cvReleaseImageHeader(&impano);*/
}

void  LG::imcor_wb_para( Frame& img_left, Frame& img_right, int L_line, int R_line, double* P_r, double* P_g, double* P_b ) {

	int wl = img_left.width;
	int hl = img_left.height;
	int wr = img_right.width;
	int hr = img_right.height;
	int channel = img_left.channel;
	int h, w;
	h = img_left.height;
	w = img_left.width;

	int* dp; //融合区域差值内存当前位置指针
	uchar *ptr_l; //左边融合区域位置指针
	uchar *ptr_r; //右边融合区域位置指针
	int max_diff = 0;
	int min_diff = 10000000; //
	int temp_r, temp_g, temp_b;

	memset( diff, 0, sizeof( int ) *diff_size );
	memset( Loc, 0, sizeof( int ) *( h*cw + 1 )*dennum );

	for ( int i = 0; i < h; i++ ) {
		ptr_l = img_left.data + i*img_left.widthStep + img_left.channel  * ( L_line - 1 - blent_hw );   //指向左图融合区域每行的第一列
		ptr_r = img_right.data + i*img_right.widthStep;    //指向右图融合区域每行的第一列
		dp = diff + i*cw * 4;  //融合区差值内存每行的第一列
		for ( int j = 0; j < cw; j++ ) {

			dp[0] = abs( ( uchar ) ptr_l[0] - ( uchar ) ptr_r[0] );   //B分量绝对差值
			dp[1] = abs( ( uchar ) ptr_l[1] - ( uchar ) ptr_r[1] );   //G分量绝对差值
			dp[2] = abs( ( uchar ) ptr_l[2] - ( uchar ) ptr_r[2] );   //R分量绝对差值
			dp[3] = dp[0] + dp[1] + dp[2];  //B G R分量绝对差值之和
			max_diff = ( max_diff < dp[3] ? dp[3] : max_diff );    //存放最大值
			min_diff = ( min_diff > dp[3] ? dp[3] : min_diff );    //存放最小值
			
			ptr_l += img_left.channel;    //左融合区当前位置指针移位
			ptr_r += img_right.channel;
			dp += 4;
		}
	}

	//printf("dis_diff is %d",)
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

	int col, row;
	uchar* ptr1;
	uchar* ptr2;
	for ( int i = 0; i < maxnum; i++ ) {
		temp = Loc[maxrow*denw + i + 1];
		row = temp / cw; //分离出行标
		col = temp % cw; //分离出列标
		ptr1 = img_left.data  + row* img_left.widthStep + channel* ( L_line - 1 - blent_hw + col );
		ptr2 = img_right.data + row* img_right.widthStep + channel * col;
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
		ptr1 = img_left.data + row* img_left.widthStep + channel *( L_line - 1 - blent_hw + col );
		ptr2 = img_right.data + row* img_right.widthStep + channel *  col;
		Y_R[i] = ( uchar ) ptr1[0];
		Y_G[i] = ( uchar ) ptr1[1];
		Y_B[i] = ( uchar ) ptr1[2];
		X_R[i] = ( uchar ) ptr2[0];
		X_G[i] = ( uchar ) ptr2[1];
		X_B[i] = ( uchar ) ptr2[2];
	}
	
	//两个数组之间的一些运算
	//int n = maxnum + smaxnum;
	int n = maxnum;
	//最小二乘法求解
	poly( X_R, Y_R, n, P_r );
	poly( X_G, Y_G, n, P_g );
	poly( X_B, Y_B, n, P_b );
}

void  LG::imcor_wb( Frame& img_left, Frame& img_right, int L_line, int R_line, double* P_r, double* P_g, double* P_b ) {

	//L_line 左边缝合线中线位置
	//R_line 右边缝合线中线位置

	int wl = img_left.width;
	int hl = img_left.height;
	int wr = img_right.width;
	int hr = img_right.height;
	int channel = img_left.channel;
	int h, w;
	h = img_left.height;
	w = img_left.width ;

	////砍掉部分重叠区域后的 融合区域宽度	
	//int cw = 2 * blent_hw;
	//const int dennum = 10;  //融合区间差值最大值和最小值之间等分区间数目
	//int diff_size = h*cw * 4;
	//int* diff = (int*)malloc(sizeof(int)*diff_size);//融合区域差值存储内存   数据的范围[0 255]
	//memset(diff, 0, sizeof(int)*diff_size);
	//int denw = (h*cw + 1);  //定义每个等分区的宽度  每行的第一个位置存的是该区间点的数目 预留h*cw（融合区最大数目）放点的位置（点的数目不会超过融合区点数目）
	//int* Loc = (int*)malloc(sizeof(int)*(h*cw + 1)*dennum); //每个等分差值区域内的点的位置
	//memset(Loc, 0, sizeof(int)*(h*cw + 1)*dennum);

	//int* dp; //融合区域差值内存当前位置指针
	//uchar *ptr_l; //左边融合区域位置指针
	//uchar *ptr_r; //右边融合区域位置指针
	//int max_diff = 0;
	//int min_diff = 10000000; //
	//int temp_r, temp_g, temp_b;

	//for (int i = 0; i < h; i++) {
	//	ptr_l = img_left.data  + i*img_left.widthStep  + img_left.channel  * (L_line-1-blent_hw );   //指向左图融合区域每行的第一列
	//	ptr_r = img_right.data + i*img_right.widthStep ;    //指向右图融合区域每行的第一列
	//	dp = diff + i*cw * 4;  //融合区差值内存每行的第一列
	//	for (int j = 0; j < cw; j++) {
	//		
	//		dp[0] = abs((uchar)ptr_l[0] - (uchar)ptr_r[0]);   //B分量绝对差值
	//		dp[1] = abs((uchar)ptr_l[1] - (uchar)ptr_r[1]);   //G分量绝对差值
	//		dp[2] = abs((uchar)ptr_l[2] - (uchar)ptr_r[2]);   //R分量绝对差值
	//		dp[3] = dp[0] + dp[1] + dp[2];  //B G R分量绝对差值之和
	//		max_diff = (max_diff < dp[3] ? dp[3] : max_diff);    //存放最大值
	//		min_diff = (min_diff > dp[3] ? dp[3] : min_diff);    //存放最小值
	//		/*if(max_diff < dp[3])
	//		max_diff = dp[3] ;
	//		if (min_diff > dp[3])
	//		min_diff = dp[3];*/
	//		ptr_l += img_left.channel;    //左融合区当前位置指针移位
	//		ptr_r += img_right.channel;
	//		dp += 4;
	//	}
	//}
	////printf("dis_diff is %d",)
	//int dis_diff = max_diff - min_diff; //总体区间宽度     
	//int qu_block = (int)(1 / double(dennum)*dis_diff);     //每个区间宽度
	//int qu1 = 0;
	//int qu2 = 0;
	//int temp = 0;
	//int id;
	//for (int col = 0; col < cw; col++)
	//for (int row = 0; row < h; row++)    //融合区内每个点循环
	//{
	//	temp = diff[row*cw * 4 + 4 * col + 3];
	//	qu1 = min_diff;
	//	qu2 = qu1 + qu_block;
	//	for (int den = 0; den < dennum; den++)    //判断点属于哪个区域
	//	{
	//		if (temp >= qu1&&temp < qu2) {
	//			Loc[den*denw] = Loc[den*denw] + 1;  //每行的第一个位置 存的是属于该区间点的数目
	//			id = Loc[den*denw];     //取该区间的数目
	//			Loc[den*denw + id] = row*cw + col;    //存融合区坐标  具体的坐标范围为[0  h*cw-1]   
	//			break;
	//		}
	//		if (den == dennum - 1) { qu1 = qu2; qu2 = max_diff + 1; } //最后一个区间范围会大一点
	//		qu1 = qu2;
	//		qu2 = qu1 + qu_block;
	//	}
	//}
	////找分区数目最大  第二大 的行号
	//int maxnum = 0;
	//int maxrow = 0;
	//int smaxnum = 0;
	//int smaxrow = 0;
	//int qunum[10] = { 0 }; //每个区间点的数目
	//int sum = 0;
	//for (int i = 0; i < 10; i++) {
	//	qunum[i] = Loc[(i*denw)];
	//	sum += qunum[i];
	//	if (qunum[i]>maxnum) {
	//		maxnum = qunum[i];
	//		maxrow = i;
	//	}
	//}

	////找分区数目第二大的行号
	//for (int i = 0; i < 10; i++) {
	//	temp = Loc[(i*denw)];
	//	if (temp <maxnum && temp >smaxnum) {
	//		smaxnum = temp;
	//		smaxrow = i;
	//	}
	//}

	//int* Y_R = (int*)malloc(sizeof(int)*(maxnum + smaxnum));
	//int* X_R = (int*)malloc(sizeof(int)*(maxnum + smaxnum));
	//int* Y_G = (int*)malloc(sizeof(int)*(maxnum + smaxnum));
	//int* X_G = (int*)malloc(sizeof(int)*(maxnum + smaxnum));
	//int* Y_B = (int*)malloc(sizeof(int)*(maxnum + smaxnum));
	//int* X_B = (int*)malloc(sizeof(int)*(maxnum + smaxnum));

	//int col, row;
	//uchar* ptr1;
	//uchar* ptr2;
	//for (int i = 0; i < maxnum; i++) {
	//	temp = Loc[maxrow*denw + i + 1];
	//	row = temp / cw; //分离出行标
	//	col = temp % cw; //分离出列标
	//	ptr1 = img_left.data + row* img_left.widthStep + channel* (L_line - 1 - blent_hw+col);
	//	ptr2 = img_right.data + row* img_right.widthStep + channel * col;
	//	Y_R[i] = (uchar)ptr1[0];
	//	Y_G[i] = (uchar)ptr1[1];
	//	Y_B[i] = (uchar)ptr1[2];
	//	X_R[i] = (uchar)ptr2[0];
	//	X_G[i] = (uchar)ptr2[1];
	//	X_B[i] = (uchar)ptr2[2];
	//}

	//for (int i = maxnum; i < maxnum + smaxnum; i++) {
	//	temp = Loc[smaxrow*denw + (i - maxnum) + 1];
	//	row = temp / cw;
	//	col = temp % cw;
	//	ptr1 = img_left.data + row* img_left.widthStep + channel *(L_line - 1 - blent_hw + col);
	//	ptr2 = img_right.data + row* img_right.widthStep + channel *  col;
	//	Y_R[i] = (uchar)ptr1[0];
	//	Y_G[i] = (uchar)ptr1[1];
	//	Y_B[i] = (uchar)ptr1[2];
	//	X_R[i] = (uchar)ptr2[0];
	//	X_G[i] = (uchar)ptr2[1];
	//	X_B[i] = (uchar)ptr2[2];
	//}

	//double  p_r[2];
	//double  p_g[2];
	//double  p_b[2];
	////两个数组之间的一些运算
	////int n = maxnum + smaxnum;
	//int n = maxnum;
	////最小二乘法求解
	//poly(X_R, Y_R, n, p_r);
	//poly(X_G, Y_G, n, p_g);
	//poly(X_B, Y_B, n, p_b);

	int p_r_lh = floor( P_r[1] );
	int p_g_lh = floor( P_g[1] );
	int p_b_lh = floor( P_b[1] );

	//printf( "pr_r_lh=%d ,p_g_lh=%d ,p_b_lh=%d \n", p_r_lh, p_g_lh, p_b_lh );
	//if ( p_r[0]<0 || p_r[0]>1.0 || p_g[0]<0 || p_g[0]>1.0 || p_b[0]<0 || p_b[0]>1.0  )
	//printf( "pr_r=%f , p_g=%f ,  p_b=%f \n", p_r[0], p_g[0], p_b[0] );

	uchar* imr_ptr;
	int line = R_line;
	double theta = 0;
	int t = 0;
	

	for (int i = 0; i < h; i++) {
		imr_ptr = img_right.data + i*img_right.widthStep;
		for (int j = 0; j < line; j++) {
			t = (int)(P_r[0] * (uchar)imr_ptr[0] + P_r[1]);
			t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[0] = (int)t;
			t = (int)(P_g[0] * (uchar)imr_ptr[1] + P_g[1]);
			t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[1] = (int)t;
			t = (int)(P_b[0] * (uchar)imr_ptr[2] + P_b[1]);
			t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[2] = (int)t;
			imr_ptr += img_right.channel;
		}
		//for (int j = line + 10; j<w; j++) {
	      for ( int j = line; j<w; j++ ) {
			theta = 1 - (double)(j-line)/(double)(w-line);
			//t = pow(p_r[0], theta) *(uchar)imr_ptr[0] + p_r[1] * (theta);
			int kk = ( theta) * powblk_w;
			int kt = P_r[0] * lh_value;
			//t =(int)((( powblk[( kt - si )*powblk_w + kk] ) * ( uchar ) imr_ptr[0]) + P_r[1] * ( theta ));
			int tempow = powblk[( kt - si )*powblk_w + kk] * imr_ptr[0];
			t = ((int)(tempow + floor(P_r[1]) * Theta[j-line]))>>6 ;
			t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[0] = (int)t;
			//t = pow(p_g[0], theta) *(uchar)imr_ptr[1] + p_g[1] * (theta);
			//t = 0.5 *( uchar ) imr_ptr[1] + p_g[1] * ( ii[j - line - 10] );
			kt = P_g[0] * lh_value;
			tempow = powblk[( kt - si )*powblk_w + kk] * imr_ptr[1];
			//t = (int)(( powblk[( kt - si )*powblk_w + kk] ) * ( uchar ) imr_ptr[1] + P_g[1] * (theta));
			t = ( ( int ) ( tempow + floor( P_g[1] ) * Theta[j-line] ) ) >> 6;
		    t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[1] = (int)t;
			//t = pow(p_b[0], theta) *(uchar)imr_ptr[2] + p_b[1] * (theta);
			//t = 0.5 *( uchar ) imr_ptr[2] + p_b[1] * ( ii[j - line - 10] );
			kt = P_b[0] * lh_value;
			tempow = powblk[( kt - si )*powblk_w + kk] * imr_ptr[2];
			//t =(int)(( powblk[( kt - si ) * powblk_w + kk] ) * ( uchar ) imr_ptr[2] + P_b[1] * ( theta ));
			t = ( ( int ) ( tempow + floor( P_b[1] ) * Theta[j-line] ) ) >> 6;
			t = (t>255 ? 255 : t);
			t = (t < 0 ? 0 : t);
			imr_ptr[2] = (int)t;
			imr_ptr += img_right.channel;
		}
	}
}

void  LG::CorrectMergePano(CoordinatTable& coordtab1, CoordinatTable& coordtab2, void* data, int W, int H, int channel){

	int x1, y1, x2, y2;
	int s1, s2, s3, s4;
	int id = 0;
	unsigned int* ptr_coordL = NULL;
	unsigned int* ptr_coordR = NULL;
	unsigned char *ptr_Data = NULL;
	unsigned char *src_data = (unsigned char *)data;
	uchar* dst_ptr = NULL;
	int cnd = imgpano.channel;
	int widthStep = W*channel;
	int QH = 31;
	unsigned int tl, tr, bl, br;
	unsigned int L_tl, L_tr, L_bl, L_br;
	unsigned int R_tl, R_tr, R_bl, R_br;
	unsigned int coord_temp = 0;
	int x, y;
	int cor1, cor2;
	int W_half = (W >> 1) - 1;
	//第一个重叠区
	int panoh = imgpano.height;
	int panow = imgpano.width;
	int position1 = img1_position[0];
	int position2 = img1_position[1];

	for (int i = 0; i < panoh; i++){
		ptr_coordL = (unsigned int*)coordtab1.data + i*coordtab1.widthStep;
		ptr_coordR = (unsigned int*)coordtab2.data + i*coordtab2.widthStep;
		dst_ptr = imgpano.data + i*imgpano.widthStep;	
		for (int j = 0; j < blent_hw * 2; j++)   //第一个融合区
		{
			coord_temp = ptr_coordL[j];
			//按位提取坐标和系数
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = (coord_temp&(0x001FFC00)) >> 10 ;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);
			
			unsigned int L_r, L_g, L_b;
			L_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			L_g = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			L_b = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量

			coord_temp = ptr_coordR[(img2_position[1] - blent_hw + j)]; //块表列对应矫正图列
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = ((coord_temp&(0x001FFC00)) >> 10) + W_half;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);

			unsigned int R_r, R_g, R_b;
			R_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			R_g = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			R_b = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量

			//
			int mean_l = (L_r + L_g + L_b) / 3;
			int mean_r = (R_r + R_g + R_b) / 3;

			dst_ptr[j * cnd + 0] = (uchar)(cor_jc[j] * L_r + (1 - cor_jc[j])*R_r);
			dst_ptr[j * cnd + 1] = (uchar)(cor_jc[j] * L_g + (1 - cor_jc[j])*R_g);
			dst_ptr[j * cnd + 2] = (uchar)(cor_jc[j] * L_b + (1 - cor_jc[j])*R_b);

			if (abs(mean_l - mean_r) > 0.9*(mean_l > mean_r ? mean_l : mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_l ? L_r : R_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? L_g : R_g);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? L_b : R_b);
			}
			/*dst_ptr[j * 3 + 0] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 0] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 0]);
			dst_ptr[j * 3 + 1] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 1] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 2] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 2]);*/
		}
	}

	//图一的非重叠区
	for (int i = 0; i < panoh; i++){
		ptr_coordL = (unsigned int*)coordtab1.data + i*coordtab1.widthStep;
		ptr_coordR = (unsigned int*)coordtab2.data + i*coordtab2.widthStep;
		dst_ptr = imgpano.data + i*imgpano.widthStep;
		for (int j = blent_hw * 2; j < position2 - blent_hw; j++){
			
			coord_temp = ptr_coordL[j];
			//按位提取坐标和系数
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = (coord_temp&(0x001FFC00)) >> 10;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);

			unsigned int L_r, L_g, L_b;
			L_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			L_g = ((((tl & 0x0000ff00)>>8)*s1 + ((tr & 0x0000ff00)>>8)*s3 + ((bl & 0x0000ff00)>>8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			L_b = ((((tl & 0x00ff0000)>>16)*s1 + ((tr & 0x00ff0000)>>16)*s3 + ((bl & 0x00ff0000)>>16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量
	
			dst_ptr[j * cnd + 0] = L_r;
			dst_ptr[j * cnd + 1] = L_g;
			dst_ptr[j * cnd + 2] = L_b;
		}
	}


	//第二个重叠区
	for (int i = 0; i < panoh; i++){
		ptr_coordL = (unsigned int*)coordtab1.data + i*coordtab1.widthStep;
		ptr_coordR = (unsigned int*)coordtab2.data + i*coordtab2.widthStep;
		dst_ptr = imgpano.data + i*imgpano.widthStep;
		for (int j = position2 - blent_hw; j < position2 + blent_hw; j++)   //第二个融合区
		{
			coord_temp = ptr_coordL[j];
			//按位提取坐标和系数
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = (coord_temp&(0x001FFC00)) >> 10;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);

			unsigned int L_r, L_g, L_b;
			L_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			L_g = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			L_b = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量

			coord_temp = ptr_coordR[(j - position2 + blent_hw)]; //块表列对应矫正图列
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = ((coord_temp&(0x001FFC00)) >> 10) + W_half;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);

			unsigned int R_r, R_g, R_b;
			R_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			R_g = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			R_b = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量

			//
			int mean_l = (L_r + L_g + L_b) / 3;
			int mean_r = (R_r + R_g + R_b) / 3;

			dst_ptr[j * cnd + 0] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * L_r + cor_jc[j - (position2 - blent_hw)] * R_r);
			dst_ptr[j * cnd + 1] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * L_g + cor_jc[j - (position2 - blent_hw)] * R_g);
			dst_ptr[j * cnd + 2] = (uchar)((1 - cor_jc[j - (position2 - blent_hw)]) * L_b + cor_jc[j - (position2 - blent_hw)] * R_b);

			if (abs(mean_l - mean_r) > 0.9*(mean_l > mean_r ? mean_l : mean_r))
			{
				dst_ptr[j * cnd + 0] = (mean_l > mean_l ? L_r : R_r);
				dst_ptr[j * cnd + 1] = (mean_l > mean_r ? L_g : R_g);
				dst_ptr[j * cnd + 2] = (mean_l > mean_r ? L_b : R_b);
			}
			/*dst_ptr[j * 3 + 0] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 0] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 0]);
			dst_ptr[j * 3 + 1] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 1] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 1]);
			dst_ptr[j * 3 + 2] = (uchar)(cor_jc[j] * iml_ptr[j * 3 + 2] + (1 - cor_jc[j])*imr_ptr[(img2_position[1] - blent_hw + j) * 3 + 2]);*/
		}
	}

	//图二的非重叠区
	for (int i = 0; i < panoh; i++){	
		ptr_coordR = (unsigned int*)coordtab2.data + i*coordtab2.widthStep;
		dst_ptr = imgpano.data + i*imgpano.widthStep;	
		for (int j = position2 + blent_hw; j < panow; j++)  //第二个融合区
		{
			coord_temp = ptr_coordR[(2 * blent_hw + (j - position2 - blent_hw))]; //块表列对应矫正图列
			cor2 = coord_temp&(0x0000001F);
			cor1 = (coord_temp&(0x000003E0)) >> 5;
			x1 = ((coord_temp&(0x001FFC00)) >> 10) + W_half;
			y1 = (coord_temp&(0xFFE00000)) >> 21;
			y2 = y1 + 1;
			x2 = x1 + 1;
			s1 = (QH - cor1)*(QH - cor2);
			s2 = (QH - cor1)*cor2;
			s3 = cor1 * (QH - cor2);
			s4 = cor1 * cor2;

			tl = *(unsigned int*)(src_data + y1*widthStep + channel * x1);
			tr = *(unsigned int*)(src_data + y1*widthStep + channel * x2);
			bl = *(unsigned int*)(src_data + y2*widthStep + channel * x1);
			br = *(unsigned int*)(src_data + y2*widthStep + channel * x2);

			unsigned int R_r, R_g, R_b;
			R_r = (((tl & 0x000000ff)*s1 + (tr&(0x000000ff))*s3 + (bl&(0x000000ff))*s2 + (br&(0x000000ff))*s4) >> 10);  //矫正出的r分量
			R_g = ((((tl & 0x0000ff00) >> 8)*s1 + ((tr & 0x0000ff00) >> 8)*s3 + ((bl & 0x0000ff00) >> 8)*s2 + ((br & 0x0000ff00) >> 8)*s4) >> 10);   //矫正出的g分量
			R_b = ((((tl & 0x00ff0000) >> 16)*s1 + ((tr & 0x00ff0000) >> 16)*s3 + ((bl & 0x00ff0000) >> 16)*s2 + ((br & 0x00ff0000) >> 16)*s4) >> 10); //矫正出的b分量

			dst_ptr[j * cnd + 0] = R_r;
			dst_ptr[j * cnd + 1] = R_g;
			dst_ptr[j * cnd + 2] = R_b;
		}
	}

}

void  LG::NoMerge_pano(Frame *img_left, Frame *img_right){

	int panow = imgpano.width;
	int panoh = imgpano.height;
	int cnl = img_left->channel;
	int cnr = img_right->channel;
	int cnd = imgpano.channel;
	int position1 = img1_position[0];
	int position2 = img1_position[1];
	uchar* dst_ptr = NULL;
	uchar* iml_ptr = NULL;
	uchar* imr_ptr = NULL;

	uchar tempr_l, tempg_l, tempb_l, tempr_r, tempg_r, tempb_r;
	int mean_l, mean_r;
	for (int i = 0; i < panoh; i++){
		dst_ptr = imgpano.data + i*imgpano.widthStep;
		iml_ptr = img_left->data + i*img_left->widthStep;
		imr_ptr = img_right->data + i*img_right->widthStep;

		

		/*for (int j = blent_hw * 2; j < position2 - blent_hw; j++){
		dst_ptr[j * cnd + 0] = iml_ptr[j * cnl + 0];
		dst_ptr[j * cnd + 1] = iml_ptr[j * cnl + 1];
		dst_ptr[j * cnd + 2] = iml_ptr[j * cnl + 2];
		}*/
		memcpy(dst_ptr, imr_ptr + (img2_position[1] - 1-blent_hw)*cnr, blent_hw*channel*sizeof(uchar));

		memcpy(dst_ptr + blent_hw*cnd , iml_ptr + (img1_position[0] - 1)*cnl, sizeof(uchar)*channel*(img1_position[1] - img1_position[0]));//非重叠区域

		memcpy(dst_ptr + (img1_position[1]-1)*cnd , imr_ptr + blent_hw * cnr, sizeof(uchar)*channel*(panow - position2+1));
		

		/*for (int j = position2 + blent_hw; j < panow; j++){
		dst_ptr[j * cnd + 0] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr];
		dst_ptr[j * cnd + 1] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 1];
		dst_ptr[j * cnd + 2] = imr_ptr[(2 * blent_hw + (j - position2 - blent_hw)) * cnr + 2];
		}*/
		//非重叠区域进行内存拷贝
		
	}


}

int   LG::array_min(unsigned int* array, int len ) {

	//确保选出最小值大于0
	int min = array[0];
	for ( int i = 0; i <len; i++ ) {
		if ( array[i] > 0 ) {
			min = array[i];
			break;
		}
	}

	for ( int i = 0; i <len; i++ ) {
		if ( array[i] > 0 && array[i] < min )
			min = array[i];
	}
	return min;
}

int   LG::array_max(unsigned int* array, int len ) {
	   int max = array[0];
	   for ( int i = 0; i <len; i++ ) {
		max = array[i]>max ? array[i] : max;
	   }
	   return max;
}

void  LG::img_resize_nearst( Frame& imgsrc, Frame& imgdst ) {

    int dstw = imgdst.width;
    int dsth = imgdst.height;
    //int srcw = imgsrc.width;
    //int srch = imgsrc.height;
    /*int lh = 256;
    int ratio_h = srch*lh*1.0 / dsth;
    int ratio_w = srcw*lh*1.0 / dstw;*/
    int x, y;
    unsigned int* src_data = ( unsigned int* ) imgsrc.data;
    unsigned int* dst_data = ( unsigned int* ) imgdst.data;

    for ( int i = 0; i < dsth;i++ )
        for ( int j = 0; j < dstw; j++ ) {
            x= j*ratio_w >> 8;
            y = i*ratio_h >> 8;
            dst_data[i*imgdst.width + j] = src_data[y*imgsrc.width + x];
        }
}
void  LG::InitImpara() {

    impara[0].fov = 197.3602677454423 ;
    impara[0].a =-0.1314107345728006;
    impara[0].b =-0.2406867791744577;
    impara[0].c =0.7853876993753323;
    impara[0].d = -1.377323782010086;
    impara[0].e = 7.295944897322002;

    impara[1].fov =  impara[0].fov;
    impara[1].a = impara[0].a;
    impara[1].b =  impara[0].b ;
    impara[1].c =   impara[0].c;
    impara[1].d =   impara[0].d;
    impara[1].e =  impara[0].e;
    impara[1].flipAngle =176.5845364311028;
    impara[1].rotationAngle =-1.317803465390971;
    impara[1].pitchAngle = -0.04384645862756322;
    impara[1].imcut.left=-66;
    impara[1].imcut.right =1014;
    impara[1].imcut.top=-161;
    impara[1].imcut.bottom=919;

    impara[2].fov =  impara[0].fov;
    impara[2].a = impara[0].a;
    impara[2].b =  impara[0].b ;
    impara[2].c =   impara[0].c;
    impara[2].d =   impara[0].d;
    impara[2].e =  impara[0].e;
    impara[2].flipAngle = -4.169671439394874;
    impara[2].rotationAngle =1.392650863303206;
    impara[2].pitchAngle =-2.082830706860889;
    impara[2].imcut.left=-7;
    impara[2].imcut.right=1073;
    impara[2].imcut.top=-224;
    impara[2].imcut.bottom =856;

//    impara[0].fov = 198.7109991134116;
//    impara[0].a =0.01;
//    impara[0].b = 0.01;
//    impara[0].c = 0.01;
//    impara[0].d = -0.2375260356303749;
//    impara[0].e = 25.11666733362256;
//
//    impara[1].fov =  impara[0].fov;
//    impara[1].a = impara[0].a;
//    impara[1].b =  impara[0].b ;
//    impara[1].c =   impara[0].c;
//    impara[1].d =   impara[0].d;
//    impara[1].e =  impara[0].e;
//    impara[1].flipAngle =-163.3348414776905;
//    impara[1].pitchAngle = -4.259024644509847;
//    impara[1].rotationAngle =-3.674955235106407;
//    impara[1].imcut.left=-72;
//    impara[1].imcut.right = 1022;
//    impara[1].imcut.top=-163;
//    impara[1].imcut.bottom=931;
//
//    impara[2].fov =  impara[0].fov;
//    impara[2].a = impara[0].a;
//    impara[2].b =  impara[0].b ;
//    impara[2].c =   impara[0].c;
//    impara[2].d =   impara[0].d;
//    impara[2].e =  impara[0].e;
//    impara[2].flipAngle = 14.41209463799635;
//    impara[2].pitchAngle =4.029655539998219;
//    impara[2].rotationAngle =-4.272706987749789;
//    impara[2].imcut.left=-19;
//    impara[2].imcut.right=1075;
//    impara[2].imcut.top=-233;
//    impara[2].imcut.bottom =861;

}

//
//void  LG::img_pano_merge(void* dst_data,int dstw,int dsth,int channel, void* src_data, int src_w, int src_h, int nchannel){
//
//
//    int x, y;
//    unsigned int* src_data_ptr = ( unsigned int* ) src_data;
//    unsigned int* dst_data_ptr = ( unsigned int* ) dst_data;
//
//    //缩放后的拼接线
//    int  img1_positon_scale[2];
//    img1_position_scale[0]=img
//
//    for ( int i = 0; i < dsth;i++ )   //根据最终缩放后的全景图 循环
//        for ( int j = 0; j < dstw; j++ ) {
//            x= j*ratio_w >> 8;
//            y = i*ratio_h >> 8;
//            dst_data[i*imgdst.width + j] = src_data[y*imgsrc.width + x];
//        }
//
//
//
//
//
//}
//
//















