#ifndef LG__
#define LG__


#include <math.h>
#include <string.h>
#include <malloc.h>
#include <cstdlib>

#define  pow_15  32768
#include <arm_neon.h>

typedef struct s_fishEye_Para {

	int     imgsrc_W;
	int     imgsrc_H;
	int		L_center_x;
	int		L_center_y;
	float	L_pitchAngle;
	float	L_rotationAngle;
	float	L_flipAngle;

	int		R_center_x;
	int		R_center_y;
	double	R_pitchAngle;
	double	R_rotationAngle;
	double	R_flipAngle;

	float	R;
	float	r;
	double	fov;
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

	int final_width;
	int final_height;

	int copy_width;
	int copy_height;

	float a, b, c;
	float d, e;

	int	cut_width;
} Para;

typedef unsigned char  uchar;

typedef struct _Frame {
	unsigned char*  data;
	int width;
	int height;
	int channel;
	int widthStep;
} Frame;


class LG  {
public:
	Para para;
	int *coordinate1;
	int *coordinate2;
	int * distance1;
	int * distance2;

    void  Panointerface( void* data, int width, int height, int channel, Frame *out_frame,int flag_wb);  //�ӿ�
    void  ImageCorrect( int* coordinate1, int* distance1, void* data, int W, int H, int channel, Frame& img_cor1, int* coordinate2, int* distance2, Frame& img_cor2, int h, int w );
    void  ImageCorrect_one( int* coordinate, int* distance, void* data, int W, int H, int channel, Frame& img_cor, int h, int w );
    void  IniteFrame( Frame& frame, int w, int h, int channel );
	void  DeleteFrame( Frame& frame );
	void  ComputeCoordinate( int* coordinate, int* distance, int w, int h, int W, int H, float fy, float xz, float fz, int center_x, int center_y, Para& fish_para );
	//void  ImageCorrect( int* coordinate1, float* distance1, void* data, int W, int H, int channel, Frame& img_cor1, int* coordinate2, float* distance2, Frame& img_cor2, int h, int w );
	void  ImMerge( Frame *img_left, Frame *img_right, int merge_x, int merge_y, Frame* img_dst, int cutwidth );
	void  splitFrame( Frame& imframe, Frame& im1, Frame& im2 );
	void  LGInit(int width, int height);
	void  imageMerge_wb( Frame& img_left, Frame& img_right, int merge_x, int merge_y, Frame& img_dst, int cutwidth, int wbw );
	void  poly( int* Y, int* X, int len, double* p );

};









#endif