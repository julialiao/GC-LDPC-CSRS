// QC-LDPC parameters

#define CPM_SIZE		150

// Type-1 GC-LDPC
#define T1_BASE_MATRIX_q	151
#define T1_BASE_MATRIX_l	50
#define T1_BASE_MATRIX_r	3
#define T1_BASE_MATRIX_t	3
#define T1_BASE_MATRIX_s	2
#define T1_BASE_MATRIX_m	4
#define T1_BASE_MATRIX_n	47

#define RANDOM_MN	


// RS-QC-LDPC parameters
/*
#define BASE_MATRIX_C1	4
#define BASE_MATRIX_C2	4
#define BASE_MATRIX_D1	{2,2,2,2,2,2,2,2}//{2,2, 3, 2, 3, 2, 2,2 }
#define BASE_MATRIX_P_ALL {3,5,3,5,3,5,3,5}//{3,5,3,5,3,5,5,5}
#define BASE_MATRIX_D2	2
#define BASE_MATRIX_LAMBDA 13
#define BASE_MATRIX_P1	3
#define BASE_MATRIX_P2	5
#define BASE_MATRIX_P	(BASE_MATRIX_C1+BASE_MATRIX_C2)
#define BASE_MATRIX_PxL	(BASE_MATRIX_LAMBDA*(BASE_MATRIX_C1+BASE_MATRIX_C2))
*/


//#define CYCLIC_PAD_COL		10

//#define N_MASK_ROW	13// 64 //10
//#define N_MASK_COL	2000
//#define RANDOM_ROW_MASK
//#define RANDOM_COL_MASK
#define ROW_MASK_BEGIN 0
#define ROW_MASK_END  CPM_SIZE*2-1

#define BASE_MATRIX_COL_OFFSET	 {0,T1_BASE_MATRIX_n-1,T1_BASE_MATRIX_n*2-2} 
//{0,BASE_MATRIX_LAMBDA-4,BASE_MATRIX_LAMBDA*2-8,BASE_MATRIX_LAMBDA*3-12,BASE_MATRIX_LAMBDA*4-16,BASE_MATRIX_LAMBDA*5-20,BASE_MATRIX_LAMBDA*6-24,BASE_MATRIX_LAMBDA*7-28} 
#define N_COL_BLK	T1_BASE_MATRIX_n*T1_BASE_MATRIX_t -2// 

// min-sum decoding related
//#define INTERLACE_MINSUM

#define MAX_DECODE_ITERATION 20
#define ROW_SHUFFLE
#define CNU_MAX_LLR_ABS 128
#define MAX_LLR_ABS	CNU_MAX_LLR_ABS*16*MAX_DECODE_ITERATION	//assume column degree <, or = 16

#define CNU_SCALE1	0.5
#define CNU_SCALE2  0.5
#define CNU_DEG_TH	60

//#define DEBUG_DATA
