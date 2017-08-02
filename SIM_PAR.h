//////////////// base matrix related
// RS-QC-LDPC parameters
#define GF_M	8

#define BASE_MATRIX_C1	2
#define BASE_MATRIX_C2	1
#define BASE_MATRIX_D1	{4,4,4}
#define BASE_MATRIX_P_ALL {3,3,3}
#define BASE_MATRIX_D2	2
#define BASE_MATRIX_LAMBDA 47
#define BASE_MATRIX_P1	3
#define BASE_MATRIX_P2	1
#define BASE_MATRIX_P	(BASE_MATRIX_C1+BASE_MATRIX_C2)
#define BASE_MATRIX_PxL	(BASE_MATRIX_LAMBDA*(BASE_MATRIX_C1+BASE_MATRIX_C2))
#define CPM_SIZE		(BASE_MATRIX_LAMBDA*BASE_MATRIX_P1*BASE_MATRIX_P2)
//#define BASE_CSHIFT	1

//#define BASE_MATRIX_MASKING
//#define BASE_MATRIX_LOCAL_DOWNSAMPLE	2
#define BASE_MATRIX_GLOBAL_DOWNSAMPLE	2

#define BASE_MATRIX_COL_OFFSET	 {0,BASE_MATRIX_LAMBDA-5,BASE_MATRIX_LAMBDA*2-11}

#define N_COL_BLK	(BASE_MATRIX_P*BASE_MATRIX_LAMBDA-11)
#define N_MASK_ROW	108// 64 //10
#define N_MASK_COL	50
//#define RANDOM_ROW_MASK
//#define RANDOM_COL_MASK
#define EVEN_LOCAL_ROW_MASK	
#define EVEN_LOCAL_COL_MASK
#define ROW_MASK_BEGIN 0
#define ROW_MASK_END  CPM_SIZE*2-1



////////// decoding related parmaters
#define MAX_DECODE_ITERATION 20
#define ROW_SHUFFLE
//#define INTERLACED_MINSUM

//#define HARD_DECODE
//#define TWO_BIT_SOFT
//#define FIX_POINT_EN
//#define Q_BIT 2
//#define SOFT_Q_BIT


#define CNU_MAX_LLR_ABS 128
#define MAX_LLR_ABS	CNU_MAX_LLR_ABS*16*MAX_DECODE_ITERATION	//assume column degree <, or = 16

#define CNU_SCALE1	0.5
#define CNU_SCALE2  0.5
#define CNU_DEG_TH	60

//#define DEBUG_DATA
