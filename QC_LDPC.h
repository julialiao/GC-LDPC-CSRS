#ifndef QC_LDPC_H
#define QC_LDPC_H




#include<vector>
#include <iostream>
#include <numeric>
#include "SIM_PAR.h"
#include<math.h>
#include<random>
#include <stdio.h>
#include <fstream>
using namespace std;
class QC_LDPC{


public:

	int N;
	int M;
	int K;
	int dataLength;
	vector<vector <bool> > H;
	vector<vector <int> > baseH;
	vector<vector <int> > nzIdxG;
	vector<int> columnPermIdx;
	
	void baseMatrixMasking();
	void generatorMatrixGen();
	void encode(vector<bool> & data, vector<bool> & codeword);

	void minSumDecode(vector<float> & llr, vector<bool> & decoded);

	void minSumDecode(vector<int> & llr, vector<bool> & decoded);

	QC_LDPC();
	~QC_LDPC();
	int iterationCnt;


	void initDecoder();
	int nRowBlockH;
	int nColBlockH;
	int nLocalRowBlk;

	void B2H();
	void maskingH();
	

	vector<bool>rowMask;
	vector<bool>colMask;
	vector<bool>maskBits;
	vector< vector<int> >baseHOffsetJ;
	vector <int> cnDegree;
	

	// decoding
	
 	bool checkSumH(vector<bool> & cword);

	void VNU(int vIndex, float llr, int j); 
	void minSumCNU(vector<float*>& sumLLRV2C, vector<float> &mC2V, vector<bool> bitMsk);
	vector< vector  <float> >  memC2V;	
	vector <float>   memSumLLRV2C;	

	void VNU(int vIndex, int llr, int j); 
	void minSumCNU(vector<int*>& sumLLRV2C, vector<int> &mC2V, vector<bool> bitMsk);
	vector< vector  <int> >  memC2V_i;
	vector <int>   memSumLLRV2C_i;

private:

};



#endif
