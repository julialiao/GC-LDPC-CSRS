#include "COMM_BOX.h"

COMM_BOX::COMM_BOX(void)
{
}

COMM_BOX::~COMM_BOX(void)
{
}

void COMM_BOX::GenData(unsigned mode, unsigned len, vector<bool> & data){
	
	if(mode==2){ // random
		data.resize(len);
		std::random_device rd;  //Will be used to obtain a seed for the random number engine
    	std::mt19937 generator(rd()); 
		//std::default_random_engine generator;
  		std::uniform_int_distribution<int> distribution(0,1);
		
		for(int k=0; k <len; ++k){
			data[k] = distribution(generator);
		}
		
	}
	else if(mode == 0){ // all zero
		data.resize(len,false);
	}
	else if(mode == 1){ // all one
		data.resize(len,true);
	}
	else if(mode ==3){
		data.resize(len);
		for(int k=0; k <len; ++k){
			data[k] = k%2;
		}
	}
}
void COMM_BOX::BPSK(vector<bool> & codeword, vector<float> & modSig){

	modSig.resize(codeword.size());
	for(int i=0; i < codeword.size();++i){
		modSig[i] = (float)1.0 - (float)2.0*(float)codeword[i];
		//modSig[i]*=sigMag;
	}
}

void COMM_BOX::AWGN(float snr, vector<float>& in, vector<float> &out){
	// 10*log10(1.0/noisePower) = snr, sqrt(noisePower ) = pow(-snr/10,10)

	srand(time(NULL));
	double noiseStd = pow(10.0,-snr/20.0)*sqrt(0.5);

	//cout << "noiseStd " << noiseStd << endl;	
	out.resize(in.size());

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); 
	//default_random_engine generator;
    normal_distribution<double> distribution(0.0,noiseStd); 
	float np = 0;
	for(int i=0; i < in.size(); ++i){
		 
		out[i] = in[i] + (float)distribution(generator);
		 
	}
	 
}

double COMM_BOX::rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

void COMM_BOX::AWGN2(float snr, vector<float>& in, vector<float> &out){
	// 10*log10(1.0/noisePower) = snr, sqrt(noisePower ) = pow(-snr/10,10)

	srand(time(NULL));
	double noiseStd = pow(10.0,-snr/20.0)*sqrt(0.5);

	//cout << "noiseStd " << noiseStd << endl;	
	out.resize(in.size());

	
	for(int i=0; i < in.size(); ++i){
		
		out[i] = in[i] + rand_normal(0, noiseStd) ;
		 
	}
	 
}

unsigned COMM_BOX::ComputeBER(vector<bool> & decoded, vector<bool> & correct ){
	unsigned s=0;
	
	if(decoded.size() != correct.size()){
		cout << "COMM_BOX::ComputeBER: input output size mismatch" << endl;
		return 0;
	}
	for(int i=0; i < decoded.size(); ++i){
		s+= (unsigned)(decoded[i]^ correct[i]);
		/*if(decoded[i]^ correct[i]){
			cout << i << "," << correct[i] << "," << decoded[i] << endl;
		}*/
	}

	return s;
}
unsigned COMM_BOX::ComputeBER(vector<bool> & decoded, vector<bool> & correct, vector<bool> mask ){
	unsigned s=0;
	
	if(decoded.size() != correct.size()){
		cout << "COMM_BOX::ComputeBER: input output size mismatch" << endl;
		return 0;
	}
	for(int i=0; i < decoded.size(); ++i){
		s+= (unsigned)(decoded[i]^ correct[i]&mask[i]);
		/*if(decoded[i]^ correct[i] && mask[i]){
			cout << i << "," << correct[i] << "," << decoded[i]<< ", mask " << mask[i] << endl;
		}*/
	}

	return s;
}
void COMM_BOX::BPSK_DeMod(vector<float> & rx, vector<bool> &demod ){

	 demod.resize(rx.size());

	for(int i=0; i < rx.size(); ++i){
		demod[i] = (rx[i]>0) ? 0 : 1;
	}
}

void COMM_BOX:: hardInputGen(vector<float>& in,vector<int>& out, int qBit){

	out.resize(in.size());
	if(qBit==1){
		for(int i=0; i < in.size(); ++i){
			out[i] = (in[i]  >0) ? 7 : -7;
		}
	}
	else if(qBit==2){
		for(int i=0; i < in.size(); ++i){
			float tmp = in[i]*8;
			if(tmp > 0){
				out[i] = tmp > 3.245 ? 5 :2;
			}
			else{
				out[i] = tmp < -3.245 ? -5 :-2;
			}
		}
	}
	
}

void COMM_BOX:: softQBitGen(vector<float>& in,vector<float>& out, int qBit){

	out.resize(in.size());
	if(qBit==1){
		for(int i=0; i < in.size(); ++i){
			out[i] = (in[i]  >0) ? 1 : -1;
		}
	}
	else if(qBit==2){
		for(int i=0; i < in.size(); ++i){
			float tmp = in[i]*8;
			if(tmp > 0){
				out[i] = tmp > 3.245 ? 5 :2;
			}
			else{
				out[i] = tmp < -3.245 ? -5 :-2;
			}
		}
	}
	
}


