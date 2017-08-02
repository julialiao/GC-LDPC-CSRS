#include "GC_RS_LDPC.h"

GC_RS_LDPC::GC_RS_LDPC()
{
	nLocalRowBlk = baseHd1[0];
	for(int p = 1; p < BASE_MATRIX_P; ++p)
	{
		nLocalRowBlk+= baseHd1[p];
	}

	 

	nRowBlockH = nLocalRowBlk + BASE_MATRIX_D2;

	nColBlockH = N_COL_BLK;
	
	M = nRowBlockH*CPM_SIZE;
	N = nColBlockH*CPM_SIZE;
	
}


GC_RS_LDPC::~GC_RS_LDPC()
{
}

void GC_RS_LDPC::GC_RS_BaseMatrixGen(){
	
	baseH.resize(nRowBlockH);

	//local parts


	
	int	iOffset = 0;
	int jOffset = 0;
	
	for (int p = 0; p < BASE_MATRIX_P; ++p) {
	
		int currP = baseHP[p];
		jOffset = baseJOffset[p];
		for (int i = 0; i < baseHd1[p]; ++i) {
			vector<int> R0(BASE_MATRIX_LAMBDA,-1);
			baseH[iOffset + i].resize(nColBlockH,-1);
			
			for (int j = 0; j < BASE_MATRIX_LAMBDA; ++j) {
				int idxRot = ((i + 1)*j*currP)%CPM_SIZE;
				baseH[iOffset + i][jOffset +j]= idxRot;	 
				R0[j] =idxRot;

			}
#ifdef BASE_CSHIFT
			for (int j = 0; j < BASE_MATRIX_LAMBDA; ++j) {
				baseH[iOffset + i][jOffset + (i+j)%BASE_MATRIX_LAMBDA] = R0[j] ;
			}
				 
#endif
#ifdef BASE_MATRIX_LOCAL_DOWNSAMPLE
			 
			for (int j = (jOffset+p+i)%BASE_MATRIX_LOCAL_DOWNSAMPLE; j < BASE_MATRIX_LAMBDA; j+= BASE_MATRIX_LOCAL_DOWNSAMPLE) {
				
				baseH[iOffset + i][jOffset +j]= -1;	 
			}	

#endif


		}
		iOffset += baseHd1[p];
		
	}



	// global part
	for (int i = 0; i < BASE_MATRIX_D2; ++i) {
		baseH[iOffset+i].resize(nColBlockH, -1);
		vector<int> R0(nColBlockH,-1);
		for (int j = 0; j < nColBlockH; ++j) {
			int idxRot = (i + 1)*j%CPM_SIZE;
			baseH[iOffset + i][j] = idxRot;
			R0[j] = idxRot;
		}
#ifdef BASE_CSHIFT
		for (int j = 0; j < nColBlockH; ++j) {
			baseH[iOffset + i][(i+j)%nColBlockH] = R0[j] ;
		}
				 
#endif
#ifdef BASE_MATRIX_GLOBAL_DOWNSAMPLE
			 
		for (int j = i; j < nColBlockH; j+= BASE_MATRIX_GLOBAL_DOWNSAMPLE) {
			
			baseH[iOffset + i][j]= -1;	 
		}	

#endif
	}


	baseHOffsetJ.resize(baseH.size());
	for (int i = 0; i < baseH.size(); ++i) {
		baseHOffsetJ[i].resize(nColBlockH);
		int validIdx = 0;
		for (int j = 0; j < nColBlockH; ++j) {
			if (baseH[i][j] >= 0) {
				baseHOffsetJ[i][j] = validIdx;
				++validIdx;
			}
		}
	}

	for(int i=0; i < baseH.size(); ++i){

		for(int j=0; j < baseH[i].size(); ++j){
			cout <<  baseH[i][j] << ",";
		}
		cout << endl;
	}
/*
	for(int j=0; j < baseH[0].size(); ++j){
		int s=0;
		for(int i=0; i < baseH.size();++i){
			s+= baseH[i][j] >=0;
		}
		cout << j << ": " << s << endl;
	}
*/

}

void GC_RS_LDPC::recordH(){
/////////////////////// 
///////////////////////record results

	
	
	cnDegree.resize(baseH.size());
	for(int i=0; i < baseH.size(); ++i){
		int degree = 0;
		for(int j = 0; j < baseH[i].size();++j){
			degree += baseH[i][j]>=0;
		}
		cnDegree[i]=degree;
	}

	char fname[200];
	sprintf(fname,"row_mask_%d_%d.txt",N,M);
	ofstream  ofs; ofs.open(fname); 

	for(int i=0; i <rowMask.size(); ++i){
		if(rowMask[i])		
			ofs << i << endl;
	}
	ofs.close();

	sprintf(fname,"col_mask_%d_%d.txt",N,M);
    ofstream  ofs1; ofs1.open(fname); 

	for(int i=0; i <colMask.size(); ++i){
		if(colMask[i])		
			ofs1 << i << endl;
	}
	ofs1.close();

	sprintf(fname,"simPar_%d_%d.txt",N,M);
	ofstream  ofs2;ofs2.open(fname); 

	ofs2 << "p = " << BASE_MATRIX_P << endl;
	ofs2 << "D1 = "; for(int i=0; i < BASE_MATRIX_P; ++i)  {ofs2 << baseHd1[i] << ",";} ofs2 << endl;
	ofs2 << "D2 = " << BASE_MATRIX_D2 << endl;
	ofs2 << "pxlambda = " << BASE_MATRIX_PxL << endl;
	ofs2 << "lambda = " << BASE_MATRIX_LAMBDA << endl;
	
	ofs2.close();

	sprintf(fname, "row_degree_%d_%d.txt", N, M);
	ofstream  ofsrd; ofsrd.open(fname);
	vector<int> rDeg(nColBlockH+1,0);
	for(int i=0; i <H.size();++i){
		int d=0;
		for(int j=0; j < H[i].size(); ++j){
			d += H[i][j];
		}
		++rDeg[d];
	}

	for (int i = 0; i < rDeg.size(); ++i) {
		ofsrd << i << "\t" << rDeg[i] <<"\t" << (float)rDeg[i]/(float)M << endl;

	}
	ofsrd.close();
	 
	

	int deg2 = 0;
	sprintf(fname, "col_degree_%d_%d.txt", N, M);
	ofstream  ofscd; ofscd.open(fname);
	vector<int> cDeg(nRowBlockH+1, 0);
	for(int j=0; j < H[0].size(); ++j){
		

		int d=0;
		for(int i=0; i <H.size();++i){
			d += H[i][j];
		}
		++cDeg[d];
		if(d<3 && d >0){
			cout << j << ":" << d << endl;
			deg2 ++;		
		}
	}

	for (int i = 0; i < cDeg.size(); ++i) {
		ofscd << i << "\t" << cDeg[i] << "\t" << (float)cDeg[i] / (float)N << endl;
	}
	ofscd.close();
#ifdef N_MASK_COL
	cout << "low degree proportion is " << (float)deg2/(float)(N-N_MASK_COL) << endl;
#else
	cout << "low degree proportion is " << (float)deg2/(float)(N) << endl;
#endif

}

void GC_RS_LDPC::parityCheckMatrixGen(){
	
	GC_RS_BaseMatrixGen();

#ifdef BASE_MATRIX_MASKING
	baseMatrixMasking();
#endif
	
	maskingH();

	B2H();
	
	recordH();
	

}



void GC_RS_LDPC::interlacedMinsum(vector<float> & llr, vector<bool> &decoded){
		
	// initialize

	initDecoder();
	memSumLLRV2C = llr;
	vector<bool> cword(N);
	//iterative decoding

	iterationCnt = MAX_DECODE_ITERATION;

	int CN_Enable = 0;
	
	for(int itr=0; itr < MAX_DECODE_ITERATION; ++itr){

		CN_Enable = itr%2;
		int iOffset = 0;
		int cnIdx = 0;
		for(int p=0;  p< BASE_MATRIX_P;++p){
			
			if(p%2 == CN_Enable){
				
				for(int i=0; i < baseHd1[p]; ++i){
					for(int k=0; k < CPM_SIZE; ++k){
						int vOffset =0;
						vector<float *> sumLLRV2C;
						vector<bool> bitMask;
						for(int j = 0; j < nColBlockH; ++j){
							if(baseH[iOffset+i][j]>=0){
								int vnIdx = vOffset + (baseH[iOffset+i][j] + k)%CPM_SIZE;						 
								sumLLRV2C.push_back(&memSumLLRV2C[vnIdx]);
								bitMask.push_back(colMask[vnIdx]);
							}
							vOffset += CPM_SIZE;
						}

						if(!rowMask[cnIdx]){
							minSumCNU(sumLLRV2C, memC2V[cnIdx],bitMask);
						}
						cnIdx++;
					}// k
				}//i
				iOffset += baseHd1[p];				
			} //if
			else{
				iOffset += baseHd1[p];
				cnIdx += baseHd1[p]*CPM_SIZE;
			}
			
		}// for p
		

		// global check node update
		for(int i=nLocalRowBlk; i <baseH.size(); ++i ){
			for(int k=0; k <CPM_SIZE; ++k){	
				int vOffset = 0;
				
				vector<float *> sumLLRV2C;//(BASE_MATRIX_PxL);
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size(); ++j){
					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						
						sumLLRV2C.push_back(&memSumLLRV2C[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 
						 
					}
					vOffset += CPM_SIZE;	
				}
				 
				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V[cnIdx],bitMask);
				}
				++cnIdx;
			}
			
		}
		
	
		//var. node update (accumulation) and decode
		int vnIdx = 0;
		
		for(int j=0; j <nColBlockH; ++j ){
			for(int k = 0; k < CPM_SIZE; ++k){
#ifndef ROW_SHUFFLE
			if(!colMask[vnIdx])
				VNU(vnIdx,llr[vnIdx],j);	
#endif
 	
			if(!colMask[vnIdx]){
				cword[vnIdx] = (memSumLLRV2C[vnIdx]< 0);	
			}
			else{
				cword[vnIdx] = 0;
			}
 			++vnIdx;
			}
		}
	/*	
		if (itr > 10) {
			cout << 10;
			for (int i = 0; i < N; ++i) {
				if (llr[i] * memSumLLRV2C[i] < 0) {
					cout << i << "," << llr[i] << "," << memSumLLRV2C[i] << endl;
				}
			}
		}*/
		// parity check
		if( checkSumH(cword)==false){
			//cout << "early termination @ " << itr<< endl;
			iterationCnt = itr+1;
			break;
		}

	}

	

	// rearrange the data
	decoded.resize(dataLength);
	int startIdx = N-dataLength;
	for(int j=0; j < dataLength; ++j){
		decoded[j] = cword[columnPermIdx[j+startIdx]];
	}

}


void GC_RS_LDPC::interlacedMinsum(vector<int> & llr, vector<bool> &decoded){
		
	// initialize

	initDecoder();
	memSumLLRV2C_i = llr;
	vector<bool> cword(N);
	//iterative decoding

	iterationCnt = MAX_DECODE_ITERATION;

	int CN_Enable = 0;
	
	for(int itr=0; itr < MAX_DECODE_ITERATION; ++itr){

		CN_Enable = itr%2;
		int iOffset = 0;
		int cnIdx = 0;
		for(int p=0;  p< BASE_MATRIX_P;++p){
			
			if(p%2 == CN_Enable){
				
				for(int i=0; i < baseHd1[p]; ++i){
					for(int k=0; k < CPM_SIZE; ++k){
						int vOffset =0;
						vector<int *> sumLLRV2C;
						vector<bool> bitMask;
						for(int j = 0; j < nColBlockH; ++j){
							if(baseH[iOffset+i][j]>=0){
								int vnIdx = vOffset + (baseH[iOffset+i][j] + k)%CPM_SIZE;						 
								sumLLRV2C.push_back(&memSumLLRV2C_i[vnIdx]);
								bitMask.push_back(colMask[vnIdx]);
							}
							vOffset += CPM_SIZE;
						}

						if(!rowMask[cnIdx]){
							minSumCNU(sumLLRV2C, memC2V_i[cnIdx],bitMask);
						}
						cnIdx++;
					}// k
				}//i
				iOffset += baseHd1[p];				
			} //if
			else{
				iOffset += baseHd1[p];
				cnIdx += baseHd1[p]*CPM_SIZE;
			}
			
		}// for p
		

		// global check node update
		for(int i=nLocalRowBlk; i <baseH.size(); ++i ){
			for(int k=0; k <CPM_SIZE; ++k){	
				int vOffset = 0;
				
				vector<int *> sumLLRV2C;//(BASE_MATRIX_PxL);
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size(); ++j){
					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						
						sumLLRV2C.push_back(&memSumLLRV2C_i[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 
						 
					}
					vOffset += CPM_SIZE;	
				}
				 
				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V_i[cnIdx],bitMask);
				}
				++cnIdx;
			}
			
		}
		
	
		//var. node update (accumulation) and decode
		int vnIdx = 0;
		
		for(int j=0; j <nColBlockH; ++j ){
			for(int k = 0; k < CPM_SIZE; ++k){
#ifndef ROW_SHUFFLE
			if(!colMask[vnIdx])
				VNU(vnIdx,llr[vnIdx],j);	
#endif
 	
			if(!colMask[vnIdx]){
				cword[vnIdx] = (memSumLLRV2C_i[vnIdx]< 0);	
			}
			else{
				cword[vnIdx] = 0;
			}


 			++vnIdx;
			}

		}
	/*	
		if (itr > 10) {
			cout << 10;
			for (int i = 0; i < N; ++i) {
				if (llr[i] * memSumLLRV2C[i] < 0) {
					cout << i << "," << llr[i] << "," << memSumLLRV2C[i] << endl;
				}
			}
		}*/
		// parity check
		if( checkSumH(cword)==false){
			//cout << "early termination @ " << itr<< endl;
			iterationCnt = itr+1;
			break;
		}

	}

	

	// rearrange the data
	decoded.resize(dataLength);
	int startIdx = N-dataLength;
	for(int j=0; j < dataLength; ++j){
		decoded[j] = cword[columnPermIdx[j+startIdx]];
	}

}
