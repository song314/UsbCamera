#ifndef LG__
#define LG__

#include <math.h>
#include <string.h>
#include <malloc.h>
#include <cstdlib>
#include   <time.h>
#include   <iostream>
#define  pow_15  32768
#define  pow_5   32
#define  pow_10  1024
#define  pow_21  2097152
#define  PI  3.14159265359


typedef struct s_fishEye_Para {
	int     imgsrc_W;
	int     imgsrc_H;
	int		L_center_x;
	int		L_center_y;
	float	L_pitchAngle; //fy
	float	L_rotationAngle;
	float	L_flipAngle;  //fz
	float   L_fz_xz;
	int		R_center_x;
	int		R_center_y;
	float	R_pitchAngle;
	float	R_rotationAngle;
	float	R_flipAngle;
	float   R_fz_xz;
	float	R;
	float	r;
	float	fov;
	int		correct_width;
	int		correct_height;
	int		merge1_x;
	int		merge1_y;
	int		merge1_width;
	int		merge1_height;
	int		merge2_x;
	int		merge2_y;
	int		merge2_width;
	int		merge2_height;
	int    final_width;
	int    final_height;
	int    pano_width;
	int    pano_height;
	int    copy_width;
	int    copy_height;
	float  a, b, c;
	float  d, e;
	int	cut_width;
} Para;

typedef unsigned char  uchar;

typedef struct _Frame {
	unsigned char*  data;
	int width=0;
	int height=0;
	int channel=0;
	int widthStep=0;
} Frame;

typedef struct _CoordinatTable {
	void*  data;
	int width = 0;
	int height = 0;
	int channel = 0;
	int widthStep = 0;
} CoordinatTable;

typedef	struct _imCut{
	float left;
	float right;
	float top;
	float bottom;
}ImCut;

typedef struct _impara{
	int W;
	int H;

	double fov;
	//桶形畸变参数
	double  a;
	double  b;
	double  c;

	//水平和垂直 偏移量
	double d;
	double e;

	//旋转 翻转 俯仰角度
	double  flipAngle; //fanzhuan
	double  pitchAngle;
	double  rotationAngle;

	float  centerx;
	float  centery;

	float R;
	float r;

	//鱼眼图的裁剪区域
	ImCut imcut;

	int inleft_position;
	int outleft_position;
	int outright_position;
	int inright_position;

} ImPara;

typedef struct _FishBlkInfo {
	int x; //外接矩形框起点列坐标
	int y;
	int width;
	int height;
	int pix_num;
	int ptr;
	int jo1;
	int jo2;
	int pb;
	int sp_num;
} FishBlkInfo;


class LG{
public:
	Para para;
	ImPara impara[3]; //每个鱼眼图的参数，从PT工程中获取数据
	CoordinatTable coordtab1;
	CoordinatTable coordtab2;
	int picture_num;  //鱼眼图像个数	
	int blent_hw ;
	int img1_position[2];
	int img2_position[2];
	int channel;
	Frame  imgcor1;	
	Frame imgcor2;
	Frame imgpano;
    Frame  imgdst;
	float* cor_jc;

	//块表相关计算
	int blkw;
	int blkh;
	int blk_w_num;
	int blk_h_num;
	FishBlkInfo* fishblkL;
	FishBlkInfo* fishblkR;
	int maxpicnum1;
	int maxpicnum2;
	unsigned int* blkdata1;
	unsigned int* blkdata2;

	//亮度调整相关
	int lh_value ;
	int si ;
	int ei ;
	int powh ;
	int cw ;
	int powblk_h;
	int powblk_w;
	int* powblk;  //指数运算查表
	int lh_Q;	
	float* cor_wb;  //亮度渐变系数表
	int cor_wb_w;
	int dennum ;    //融合区间差值最大值和最小值之间等分区间数目
	int diff_size;
	int* diff;      //融合区域差值存储内存   数据的范围[0 255]
	int  denw;      //定义每个等分区的宽度  每行的第一个位置存的是该区间点的数目 预留h*cw（融合区最大数目）放点的位置（点的数目不会超过融合区点数目）
	int* Loc;       //每个等分差值区域内的点的位置
	double pr_r[2];
	double pr_g[2];
	double pr_b[2];

	double pl_r[2];
	double pl_g[2];
	double pl_b[2];

	int* Theta;

	int* Y_R ;
	int* X_R ;
	int* Y_G ;
	int* X_G ;
	int* Y_B ;
	int* X_B ;

    //插值相关
    int insert_lh ;
    int ratio_h ;
    int ratio_w ;
	
public:
	void  Panointerface(void* data, int width, int height, int channel,  int flag_wb);
	void  IniteCoordTab(CoordinatTable& coordtab, int w, int h, int channel);
	void  IniteFrame( Frame& frame, int w, int h, int channel );
	void  DeleteFrame( Frame& frame );
	void  DeleteCoordTab(CoordinatTable& coordtab);
	void  ComputeCoordinate( int* coordinate, int* distance, int w, int h, int W, int H, float fy, float xz, float fz, int center_x, int center_y, Para& fish_para );
	void  ImageCorrect( int* coordinate1, int* distance1, void* data, int W, int H, int channel, Frame& img_cor1, int* coordinate2, int* distance2, Frame& img_cor2, int h, int w );
	void  ImMerge( Frame *img_left, Frame *img_right, int merge_x, int merge_y, Frame* img_dst, int cutwidth );
	void  splitFrame( Frame& imframe, Frame& im1, Frame& im2 );
	void  LGInit(int width, int height, int nchannel, int Blent_hw,int dstw,int dsth);
	void  imageMerge_wb( Frame& img_left, Frame& img_right, int merge_x, int merge_y, Frame& img_dst, int cutwidth, int wbw );
	void  poly( int* Y, int* X, int len, double* p );
	void  ComputeCoordinate_L( unsigned int* coordinate, int w, int h, Para& fish_para );
	void  ComputeCoordinate_R( unsigned int* coordinate, int w, int h, Para& fish_para );
	void  ImageCorrect_L(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor);
	void  ImageCorrect_R(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor);
	void  ImageCorrect_L1(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor);
	void  ImageCorrect_R1(CoordinatTable& coordtab, void* data, int W, int H, int channel, Frame& img_cor);
	void  correct_Rz_coordinate(unsigned int* coordinate, int w, int h, int W, int H, Para& fish_para);
	void  correct_Lz_coordinate(unsigned int* coordinate, int w, int h, int Width, int Hight, Para& fish_para);
	
		void  correct_coordinate_R( unsigned int* coordinate, int w, int h, int Width, int Hight, Para& fish_para );
		void  correct_coordinate_L( unsigned int* coordinate, int w, int h, int Width, int Hight, Para& fish_para );
	//关于读PT工程文件相关

	void  computeimpara(ImPara& impara);
    int   MOD(int x, int y);
	float MODf(float x, float y);
	
	void  flip_angle_correct(int r, float& fanzhuan1, float& fanzhuan2, float& fanzhuan1_xz, float& fanzhuan2_xz, int img1_position[2], int img2_position[2]);
	void  Merge_pano(Frame *img_left, Frame *img_right);
	void  Merge_pano1(Frame *img_left, Frame *img_right);
	void  imcor_wb(Frame& img_left, Frame& img_right,int L_line,int R_line);

	void  CorrectMergePano(CoordinatTable& coordtab1, CoordinatTable& coordtab2, void* data, int W, int H, int channel);

	void  NoMerge_pano(Frame *img_left, Frame *img_right);
	
	
	void  IniteFishBlkInfo( FishBlkInfo & fishblk );

	int   array_min(unsigned int* array, int len );
	int   array_max(unsigned int* array, int len );	 
	void  imcor_wb_para( Frame& img_left, Frame& img_right, int L_line, int R_line, double* P_r, double* P_g, double* P_b );
	void  imcor_wb( Frame& img_left, Frame& img_right, int L_line, int R_line, double* P_r, double* P_g, double* P_b );
	void  imcor_wb_optimize( Frame& img_left, Frame& img_right, int L_line, int R_line );
    void  InitImpara();
    void  img_resize_nearst( Frame& imgsrc, Frame& imgdst );
    void  img_pano_merge(void* dst_data,int dstw,int dsth, int channel, void* src_data, int src_w, int src_h, int nchannel);  //最终全景图产生
};

#endif