#include "GC_LDPC.h"



GC_LDPC::GC_LDPC()
{
}


GC_LDPC::~GC_LDPC()
{
}

void GC_LDPC::Type1BaseMatrixGen(){

	
	 

	nRowBlockH = T1_BASE_MATRIX_m*T1_BASE_MATRIX_t + T1_BASE_MATRIX_s;
	nColBlockH = N_COL_BLK; //T1_BASE_MATRIX_n*T1_BASE_MATRIX_t;
	nLocalRowBlk = T1_BASE_MATRIX_m*T1_BASE_MATRIX_t ;

	N = nColBlockH*CPM_SIZE;
	M = nRowBlockH*CPM_SIZE;

	vector< vector<int> > R00;
	vector< vector<int>	> R0;

	//construct field element
	int alpha;
	int a;
	//find primitive element alpha
	for(alpha=2; alpha < T1_BASE_MATRIX_q; ++alpha){
		int a = alpha;
		int i;
		for(i=2; i < T1_BASE_MATRIX_q-1; ++ i){
			a = (a*alpha)%T1_BASE_MATRIX_q;
			if(a==1){
			//	cout << "early break " << alpha << endl;
				break;
			}
		}

		if(i==T1_BASE_MATRIX_q-1){
			a = (a*alpha)%T1_BASE_MATRIX_q;

			if(a==1){break;}
		}
	}

	cout << "primitive element alpha is " << alpha << endl;
	vector<int> A(T1_BASE_MATRIX_q);

	A[0] = 1;
	A[1] = alpha;
	for(int i=2; i < T1_BASE_MATRIX_q; ++i){
		A[i] = (A[i-1]*alpha ) % T1_BASE_MATRIX_q;
	}
	
	vector<int> R(T1_BASE_MATRIX_q-1);
//	cout << "The first row of B1 :" << endl; 
	for(int i=0; i < T1_BASE_MATRIX_q-1; ++i){
		int tmp =  A[i] - 1;
		for(int j=0; j < A.size();++j){
			if(tmp==0){
				R[i] = -1;
				break;
			}			
			else if (A[j] == tmp){
				R[i] = j;
				break;
			}
		}
//		cout << R[i] << ",";
	}
//	cout << endl;
	
//	cout << "Elements of R00 are :" << endl; 
	R00.resize(T1_BASE_MATRIX_l);
	
	int startIdx = 0;
	for(int i=0; i < T1_BASE_MATRIX_l; ++i){
		R00[i].resize(T1_BASE_MATRIX_l);
		for(int j=0; j < T1_BASE_MATRIX_l; ++j){
			int idx = (T1_BASE_MATRIX_q-1 + startIdx + j)%(T1_BASE_MATRIX_q-1);
			R00[i][j] = R[idx];
		}
		startIdx--;
	}
	// select m row, n columns from R00 for local
	vector<bool> rowSel(T1_BASE_MATRIX_l,false);
	vector<bool> clmSel(T1_BASE_MATRIX_l*T1_BASE_MATRIX_t,false);
	
#ifdef RANDOM_MN
	std::random_device rd1;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator1(rd1()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution1(0,T1_BASE_MATRIX_l-1);
	int nCol = 0;	
	while(nCol <T1_BASE_MATRIX_n ){
		int t = distribution1(generator1);
		if(!clmSel[t]){
			clmSel[t] = true;
			++nCol;
		}			
	}

	

	std::random_device rd2;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator2(rd2()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution2(0,T1_BASE_MATRIX_l-1);
	int nRow = 0;	
	while(nRow <T1_BASE_MATRIX_m ){
		int t = distribution2(generator2);
		if(!rowSel[t]){
			rowSel[t] = true;
			++nRow;
		}			
	}

#else
	for(int i=0; i <T1_BASE_MATRIX_m; ++i){
		rowSel[i] = true;
	}
	for(int i=0; i < T1_BASE_MATRIX_n; ++i){
		clmSel[i] = true;
	}
#endif

	for(int t=1; t < T1_BASE_MATRIX_t; ++t){
		for(int j=0; j < T1_BASE_MATRIX_l; ++j){
			clmSel[T1_BASE_MATRIX_l*t+j]=clmSel[j];
		}
	}

	startIdx = 0;
	R0.resize(T1_BASE_MATRIX_l-T1_BASE_MATRIX_m);
	int m1=0;

	
	for(int i=0; i < T1_BASE_MATRIX_l; ++i){
		if(!rowSel[i]){
			 
			for(int j=0; j < T1_BASE_MATRIX_q-1; ++j){
				
				
				if(clmSel[j] ){
					int idx = (T1_BASE_MATRIX_q-1 + startIdx + j)%(T1_BASE_MATRIX_q-1);
					
					R0[m1].push_back(R[idx]);
					 
				}
			}
			m1 += (!rowSel[i]); 	
		}
		startIdx--;
		
		//startIdx = (startIdx - 1)%(T1_BASE_MATRIX_q-1);

	}

	for(int i=0; i < R00.size(); ++i){
		for(int j=0; j < R00[0].size();++j){
//			cout << R00[i][j] << ",";
		}
//		cout << endl;
	}

//	cout << "elements of R0 are:" << endl;
	for(int i=0; i < R0.size(); ++i){
		for(int j=0; j < R0[0].size();++j){
//			cout << R0[i][j] << ",";
		}
//		cout << endl;
	}



	baseH.resize(nRowBlockH);
	for(int i=0; i <nRowBlockH; ++i){
		baseH[i].resize(nColBlockH,-1);
	}
	// local
	int n1=0;
	
	int iOffset = 0;
	int jOffset = 0;
	for(int t=0; t < T1_BASE_MATRIX_t; ++t ){
		m1 = 0;
		jOffset = localBloclkOffset[t];
		for(int i=0; i < T1_BASE_MATRIX_l; i++){
			n1 = 0;
			for(int j=0; j < T1_BASE_MATRIX_l; ++j){
				if(rowSel[i] && clmSel[j]){
					baseH[iOffset+m1][jOffset+n1] = R00[i][j];
					n1++;
				}
			
			}
			m1 += (rowSel[i]); 
		
		}
		iOffset += T1_BASE_MATRIX_m;
		//jOffset+= T1_BASE_MATRIX_n;//
		cout << jOffset << endl;
	}

	// select the first s rows for global part
	for(int s = 0; s < T1_BASE_MATRIX_s; ++ s){
		for(int j=0; j < R0[0].size(); ++j){
			
			baseH[iOffset+s][j]=R0[s][j];
		}
	}
#ifdef GLOBAL_ROW_SPLIT
	int orgSize = baseH.size();
	baseH.resize(baseH.size()+1);
	for(int s = 0; s < GLOBAL_ROW_SPLIT;++s){
		baseH[orgSize+s] = baseH[orgSize-GLOBAL_ROW_SPLIT+s];
	       for(int j=1; j < baseH[0].size();j+=2){
	       	baseH[orgSize+s][j] = -1;
		baseH[orgSize-GLOBAL_ROW_SPLIT+s][j-1] = -1;
	       }	
	}
	nRowBlockH += GLOBAL_ROW_SPLIT;
	M += (GLOBAL_ROW_SPLIT*CPM_SIZE);
#endif

	
//	cout << "The base matris is:" << endl; 
	for(int i=0; i < baseH.size();++i){
		for(int j=0; j < baseH[i].size();++j){
//			cout << baseH[i][j] << ",";
		}
//		cout << endl;
	}
}

void GC_LDPC::parityCheckMatrixGen(){
	 
	Type1BaseMatrixGen();

#ifdef BASE_MATRIX_MASKING
	baseMatrixMasking();
#endif
	maskingH();

	B2H();
	
	recordH();

}

void GC_LDPC::recordH(){


	
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

	
	ofs2 << "N = " << N << endl;
	ofs2 << "M = " << M << endl;
	ofs2 << "q = " << T1_BASE_MATRIX_q << endl;
	ofs2 << "t = " << T1_BASE_MATRIX_t << endl;
	ofs2 << "s = " << T1_BASE_MATRIX_s << endl;
	ofs2 << "m = " << T1_BASE_MATRIX_m << endl;
	ofs2 << "n = " << T1_BASE_MATRIX_n << endl;
	ofs2 << "q = " << T1_BASE_MATRIX_q << endl;
	
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

void GC_LDPC::interlacedMinsum(vector<float> & llr, vector<bool> &decoded){
		
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
		for(int p=0;  p< T1_BASE_MATRIX_t;++p){
			
			if(p%2 == CN_Enable){
				
				for(int i=0; i < T1_BASE_MATRIX_m; ++i){
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
				iOffset += T1_BASE_MATRIX_m;				
			} //if
			else{
				iOffset += T1_BASE_MATRIX_m;
				cnIdx += T1_BASE_MATRIX_m*CPM_SIZE;
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
