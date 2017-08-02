#ifndef COMM_BOX_H
#define COMM_BOX_H

#include<random>
#include<vector>
#include<math.h>
#include<iostream> 
#include <time.h>

using namespace std;

class COMM_BOX
{
public:
	COMM_BOX(void);
	~COMM_BOX(void);

	void GenData(unsigned mode, unsigned len, vector<bool> & data);
	void BPSK(vector<bool> &codeword, vector<float> & signal);
	void AWGN(float snr, vector<float>& in, vector<float> &out);
	unsigned ComputeBER(vector<bool> & decoded, vector<bool> & correct );
	unsigned ComputeBER(vector<bool> & decoded, vector<bool> & correct, vector<bool> mask );

	void BPSK_DeMod(vector<float> & rx, vector<bool> &demod );
	void AWGN2(float snr, vector<float>& in, vector<float> &out);

	double rand_normal(double mean, double stddev);

	//const float sigMag = sqrt(0.5);

	void hardInputGen(vector<float>& in, vector<int>& out, int qBit);
	void softQBitGen(vector<float>& in, vector<float>& out, int qBit);
};

#endif
