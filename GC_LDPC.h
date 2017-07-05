#ifndef GC_LDPC_H
#define GC_LDPC_H

#include "QC_LDPC.h"

class GC_LDPC: public QC_LDPC{

public:

	

	void parityCheckMatrixGen();
	void interlacedMinsum(vector<float> & llr, vector<bool> &decoded);
	GC_LDPC();
	~GC_LDPC();

private:

	// encoding
	void Type1BaseMatrixGen();
	
	void recordH();
	int  localBloclkOffset[T1_BASE_MATRIX_t] = BASE_MATRIX_COL_OFFSET;
};

#endif
