#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <time.h>
#include "SIM_PAR.h"
#include "GC_RS_LDPC.h"
#include "COMM_BOX.h"

using namespace std;
int main() {

	// initialize random seed
	srand(time(NULL));

	GC_RS_LDPC gcLdpc;
	COMM_BOX commBox;

	vector<float> snrRange;
	vector<float> uBER;
	vector<float> cBER;
	ofstream  ofs; 
	char fname[200];


	 
	gcLdpc.parityCheckMatrixGen();
	gcLdpc.generatorMatrixGen();

	

#ifdef N_MASK_COL
	int N = gcLdpc.N-N_MASK_COL;
	int K = gcLdpc.dataLength - N_MASK_COL;
#else
	int N = gcLdpc.N;
	int K = gcLdpc.dataLength;
#endif

#ifdef ROW_SHUFFLE	
	sprintf(fname,"simResults_%d_%d_sf1.txt",N,K);
	

#else
	sprintf(fname,"simResults_%d_%d_sf0.txt",N,K);
	 
#endif

	ofs.open(fname);
	float rate = (float)K/(float)N;
	cout << "LDPC (N, K) = (" <<  N << "," << K << "), "	<< "Rate = " << rate << endl;
//	ofs << "LDPC (N, K) = (" << gcLdpc.N << "," << gcLdpc.dataLength << "), "	<< "Rate = " << rate << endl;

	vector< float > iteratoinCnt;
	
	double vecSNR[10] =  { 3.6, 3.7,3.8, 3.9, 4, 4.1,4.2,4.3,4.4, 4.5};

	for (int snrIdx = 0; snrIdx <10; snrIdx ++ ) {
		double SNR = vecSNR[snrIdx] + 10*log10(rate);
		cout << "SNR = " << SNR <<endl;
		cout << "Eb/No = " << vecSNR[snrIdx] << endl;
		unsigned uNoErrorBits = 0;
		unsigned cNoErrorBits = 0;
		unsigned pktSim = 1E7;
		unsigned pkNum = 0;		unsigned nErrorPkt = 0;

		int sumItrCnt = 0; 
		int breakFlg = 0;
		for (int pkt = 0; pkt < pktSim; ++pkt) {

			++pkNum;
			// generate test data
			vector<bool> data;
			commBox.GenData(2, gcLdpc.dataLength, data);

			// encoding
			vector<bool> codeword;
			gcLdpc.encode(data, codeword);
	 
			//BPSK 
			vector<float> modSig;
			commBox.BPSK(codeword, modSig);
			//cout << "BPSK" << modSig.size() << endl;

			//AWGN channel
			vector<float> rx;
			commBox.AWGN2(SNR, modSig, rx);
			//cout << "AWGN" << rx.size() << endl;

			// 
			vector<bool> demod;
			commBox.BPSK_DeMod(rx, demod);

			uNoErrorBits += commBox.ComputeBER(demod, codeword);

			
			vector<bool> decodedCodeword;

#ifdef INTERLACED_MINSUM
			gcLdpc.interlacedMinsum(rx, decodedCodeword);
#else
			gcLdpc.minSumDecode(rx, decodedCodeword);
#endif

			sumItrCnt += gcLdpc.iterationCnt;

			unsigned nErr = commBox.ComputeBER(decodedCodeword, data, gcLdpc.maskBits);
			nErrorPkt += (nErr >0);
			cNoErrorBits += nErr;

			//record high error rate pattern
	/*		if(nErr > uNoErrorBits){
				ofstream fd; 
				fd.open("in.txt");
 
				for(int i=0; i< codeword.size(); ++i){
					fd << codeword[i] << "\t" << rx[i] << "\t" <<gcLdpc.memSumLLRV2C[i]  << endl;
					 
				}
				fd.close();
				breakFlg = 1;
				break;
			}

			*/

			 
			if (cNoErrorBits > 1000 && pkNum > 1000){ break; }
		}
		
		iteratoinCnt.push_back((float)sumItrCnt/(float)pkNum);		
		cout << "average decode iteration " << (float)sumItrCnt/(float)pkNum << endl;

		uBER.push_back((float)uNoErrorBits / (float)pkNum / (float)gcLdpc.N);
		cBER.push_back((float)cNoErrorBits / (float)pkNum / (float)K);

		cout << vecSNR[snrIdx] << "\t" << (float)uNoErrorBits / (float)pkNum / (float)gcLdpc.N << "\t" << (float)cNoErrorBits / (float)pkNum / (float)K
			<< "\t" << (float)nErrorPkt / (float)pkNum << "\t"<< pkNum << endl;

		ofs << SNR << "\t" << vecSNR[snrIdx] << "\t" << (float)uNoErrorBits / (float)pkNum / (float)gcLdpc.N << "\t" << (float)cNoErrorBits / (float)pkNum / (float)K
			<< "\t" << (float)nErrorPkt/(float)pkNum 
			<< "\t" << (float)sumItrCnt/(float)pkNum <<   "\t"<< pkNum << endl;

		ofs.flush();
		
		if(breakFlg) break;
	}

	ofs.close();
	return 0;
}
