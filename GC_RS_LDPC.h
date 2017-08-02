#ifndef GC_RS_LDPC_H
#define GC_RS_LDPC_H

#include "QC_LDPC.h"
class GC_RS_LDPC: public QC_LDPC{

public:

	

	void parityCheckMatrixGen();
	//decoding
	void interlacedMinsum(vector<float> & llr, vector<bool> &decoded);
	void interlacedMinsum(vector<int> & llr, vector<bool> &decoded);
	
	GC_RS_LDPC();
	~GC_RS_LDPC();

private:

	// encoding
	void GC_RS_BaseMatrixGen();
	
	void recordH();
	int baseJOffset[BASE_MATRIX_P] = BASE_MATRIX_COL_OFFSET;
	int baseHd1[BASE_MATRIX_P]=BASE_MATRIX_D1;
	int baseHP[BASE_MATRIX_P] = BASE_MATRIX_P_ALL;

	
};
#endif
