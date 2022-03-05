#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<ctime>
#include<cfloat>
#include"Spectral_Inter.h"
#include"Spectral_Non.h"
#ifdef _OPENMP
#include<omp.h>
#endif
using namespace std;

void Gradient(long double[14], Spectral_Inter*, Spectral_Inter*, Spectral_Non*, int);
long double PolakRibiere(long double[14], long double[14]);
void Minimize(long double[14], Spectral_Inter*, Spectral_Inter*, Spectral_Non*, int);

long double Uniform(long double, long double);
long double Protected_Uniform(long double, long double, long double, long double);
long double Chi_Square(long double[2], long double[7], int);
long double Print(Spectral_Inter*, Spectral_Inter*, Spectral_Non*, long double[2], long double[7], int); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square

long double Random_Range[14][2] = {{.3,.5},{4.,6.},{2.5,3.5},{3.,6.},{.12,.18},{3.,6.},{7.,15.},{3.5,6.},{2.5,3.5},{3.,6.},{1.65,1.70},{2.5,6.},{1.5,5.},{5.,6.}};
ofstream OutputFile;

const long double Spatial_Ratio[4][7] = {{1.,1.00006,0.99883,0.992039,0.982366,0.970341,0.95766},
				   {.99,0.988286,0.945063,0.879461,0.798659,0.7259,0.654381},
				   {.98,0.954875,0.856416,0.720447,0.573465,0.45867,0.376707},
				   {.97,0.908029,0.715435,0.524036,0.372788,0.246218,0.18}};
const long double Vacuum_Spatial[7] = {13.5965519695970589115, 0.0415680222792943527448, 0.00120126774580634698112, 4.65000126560610055077e-05, 
				  1.96861980208214073368e-06, 8.66667213197805982178e-08, 3.94190650826746653845e-09};
const long double Vacuum_Euclidean[4][2] = {{0.000259568923526295945662, 9.64220418032793835491e-06},
				       {0.00219476484784199273204, 0.000188143003166908006425},
				       {0.00821024781558290065043, 0.00116410367128954708042},
				       {0.0264690797013430510705, 0.00575549843590767333436}};
/*{0.294516,2.87077,2.61272,5.15654,0.102616,3.7999,9.70375,4.78732,2.97949,1.96609,1.77989,3.65803,4.13782,5.66908,0.178552,1.00264,1.00297,0.98933,0.969937,0.96876, 0.956142,0.951705,0.00366396},{0.378568,4.4667,2.39524,4.6254,0.16796,5.3745,9.16091,5.4519,2.9287,2.23292,1.69856,3.23026,2.5263,5.63885,0.189718,1.00872,0.978135,0.920747,0.857007,0.823711,0.734053, 0.64906,0.0031065},{0.398046,5.25466,2.61066,3.1964,0.171603,5.75971,14.9795,3.50283,3.44209,3.15538,1.6936,5.49187,3.68687,5.06835,0.105155,1.00424,0.954654,0.813603,0.597653,0.472872, 0.46059,0.40688,0.0887167},{0.317571,2.23119,3.12586,3.82102,0.120709,3.19412,9.24223,1.25108,2.00384,2.55467,1.66919,3.17136,4.74933,4.69052,0.00354849,1.00179,0.786064,0.715781,0.699182,0.695192, 0.196742,0.179014,0.557705}*/
int main(int argc, char* argv[])
{
	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.294516,2.87077},{3.09946,2.61272,5.15654},{.106597,.102616,3.7999},{10.58,9.70375,4.78732},{3,2.97949,1.96609}},
						{{.4024,.378568,4.4667},{3.125,2.16747,5.26226},{.372,.171732,5.83927},{8.97,14.7462,3.51337},{3,2.83876,4.14535}},
						{{.403047,.43634,4.78343},{3.151,1.66683,5.0539},{.68,.15422,5.69145},{9.07,10.5457,3.25066},{3,2.78677,2.26082}},
						{{.266834,.31775,2.23138},{3.1855,3.12605,3.82121},{.98,.120898,3.19431},{9.52,9.24223,1.25108},{3,2.00384,2.55467}}};
	long double Psi_Prime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,9.70375,4.78732},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.77989,3.65803},{2.45,4.13782,5.66908}},
					       {{1.59,1.71443,4.64476},{2.7,2.71014,2.31858}},
					       {{1.51,1.66247,3.71779},{2.4,2.91303,5.88634}},
					       {{1.36,1.66938,3.17155},{1.98,4.74952,4.69071}}};

	long double Medium_Spatial[7];
	long double Medium_Euclidean[2];
	/*for(int i = 0; i < 6; i++)	//Superceeded by precalculated values, standing by if services required
		Vacuum_Spatial[i] = Spatial((long double)(i)+.25, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[0][0] = Euclidean(1./.388, .194, 0, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[1][0] = Euclidean(1./.516, .258, 0, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[2][0] = Euclidean(1./.640, .320, 0, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[3][0] = Euclidean(1./.800, .400, 0, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[0][1] = Euclidean(1./.388, .194, 3, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[1][1] = Euclidean(1./.516, .258, 3, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[2][1] = Euclidean(1./.640, .320, 3, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[3][1] = Euclidean(1./.800, .400, 3, JPsi_Parameters[0], Psi_Prime_Parameters[0], Non_Parameters[0], true);*/

	Spectral_Inter* JPsi[5];
	Spectral_Inter* Psi_Prime[5];
	Spectral_Non* Non[5];
	for(int i = 0; i < 5; i++)
	{
		JPsi[i] = new Spectral_Inter(JPsi_Parameters[i], i, !bool(i));
		Psi_Prime[i] = new Spectral_Inter(Psi_Prime_Parameters[i], i, !bool(i));
		Non[i] = new Spectral_Non(Non_Parameters[i], i, !bool(i));
	}

	int Temp = atoi(argv[1]);
	long double T;
	switch(Temp)
	{
		case 1:
			T = .194;
			break;
		case 2:
			T = .258;
			break;
		case 3:
			T = .320;
			break;
		case 4:
			T = .400;
			break;
		default:
			return(2);
	}

	char File[70] = "data/Optimiser_Output.";
	if(argc == 3)
	{
		strcat(File,argv[2]);
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File,ios::app);
		if(!OutputFile.is_open())
			return(1);
	}
	else if(argc == 16)
	{
                strcat(File,"API.");
                strcat(File,argv[1]);
                strcat(File,".csv");
                OutputFile.open(File,ios::app);

		JPsi[Temp]->Replace(atof(argv[2]),0,1);
		JPsi[Temp]->Replace(atof(argv[3]),0,2);
		JPsi[Temp]->Replace(atof(argv[4]),1,1);
		JPsi[Temp]->Replace(atof(argv[5]),1,2);
		JPsi[Temp]->Replace(atof(argv[6]),2,1);
		JPsi[Temp]->Replace(atof(argv[7]),2,2);
		JPsi[Temp]->Replace(atof(argv[8]),3,1);
		JPsi[Temp]->Replace(atof(argv[9]),3,2);
		Psi_Prime[Temp]->Replace(atof(argv[8]),3,1);
		Psi_Prime[Temp]->Replace(atof(argv[9]),3,2);
		JPsi[Temp]->Replace(atof(argv[10]),4,1);
		JPsi[Temp]->Replace(atof(argv[11]),4,2);
		Non[Temp]->Replace(atof(argv[12]),0,1);
		Non[Temp]->Replace(atof(argv[13]),0,2);
		Non[Temp]->Replace(atof(argv[14]),1,1);
		Non[Temp]->Replace(atof(argv[15]),1,2);

		Medium_Euclidean[0] = JPsi[Temp]->Euclidean(.5,0)+Psi_Prime[Temp]->Euclidean(.5,0)+Non[Temp]->Euclidean(.5,0);
		Medium_Euclidean[1] = JPsi[Temp]->Euclidean(.5,3)+Psi_Prime[Temp]->Euclidean(.5,3)+Non[Temp]->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = JPsi[Temp]->Spatial((long double)(j)+.25)+Psi_Prime[Temp]->Spatial((long double)(j)+.25)+Non[Temp]->Spatial((long double)(j)+.25);
		cout << setprecision(18) << Print(JPsi[Temp], Psi_Prime[Temp], Non[Temp], Medium_Euclidean, Medium_Spatial, Temp) << endl;
		return(0);
	}
	else if(argc == 31)
	{
		strcat(File,argv[2]);
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File,ios::app);
		if(!OutputFile.is_open())
			return(1);

		Random_Range[0][0] = atof(argv[3]);
		Random_Range[0][1] = atof(argv[4]);
		Random_Range[1][0] = atof(argv[5]);
		Random_Range[1][1] = atof(argv[6]);
		Random_Range[2][0] = atof(argv[7]);
		Random_Range[2][1] = atof(argv[8]);
		Random_Range[3][0] = atof(argv[9]);
		Random_Range[3][1] = atof(argv[10]);
		Random_Range[4][0] = atof(argv[11]);
		Random_Range[4][1] = atof(argv[12]);
		Random_Range[5][0] = atof(argv[13]);
		Random_Range[5][1] = atof(argv[14]);
		Random_Range[6][0] = atof(argv[15]);
		Random_Range[6][1] = atof(argv[16]);
		Random_Range[7][0] = atof(argv[17]);
		Random_Range[7][1] = atof(argv[18]);
		Random_Range[8][0] = atof(argv[19]);
		Random_Range[8][1] = atof(argv[20]);
		Random_Range[9][0] = atof(argv[21]);
		Random_Range[9][1] = atof(argv[22]);
		Random_Range[10][0] = atof(argv[23]);
		Random_Range[10][1] = atof(argv[24]);
		Random_Range[11][0] = atof(argv[25]);
		Random_Range[11][1] = atof(argv[26]);
		Random_Range[12][0] = atof(argv[27]);
		Random_Range[12][1] = atof(argv[28]);
		Random_Range[13][0] = atof(argv[29]);
		Random_Range[13][1] = atof(argv[30]);
	}
	else
	{
		strcat(File,argv[2]);
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File);
		if(!OutputFile.is_open())
			return(1);

		long double start[10];
		long double finish[10];
		long double step[10];
		int Num_Threads = atoi(argv[3]);
		int Thread_Num = atoi(argv[4]);
		int Dims[10];
		int i;

		for(i = 0; i < 10; i++)
		{
			start[i] = atof(argv[5+3*i]);
			finish[i] = atof(argv[6+3*i]);
			step[i] = atof(argv[7+3*i]);
			Dims[i] = int((finish[i]-start[i])/step[i]+1.0000000000001);
			if(Dims[i] < 1)
				return(2);
		}

		long double** Parameter_List = new long double*[Dims[0]*Dims[1]*Dims[2]*Dims[3]*Dims[4]*Dims[5]*Dims[6]*Dims[7]*Dims[8]*Dims[9]];

		i = 0;
		for(long double A = start[0]; A <= finish[0]*1.0000000000001; A += step[0])
			for(long double PA = start[1]; PA <= finish[1]*1.0000000000001; PA += step[1])
				for(long double M = start[2]; M <= finish[2]*1.0000000000001; M += step[2])
					for(long double PM = start[3]; PM <= finish[3]*1.0000000000001; PM += step[3])
						for(long double Gamma = start[4]; Gamma <= finish[4]*1.0000000000001; Gamma += step[4])
							for(long double PGamma = start[5]; PGamma <= finish[5]*1.0000000000001; PGamma += step[5])
				for(long double MQ = start[6]; MQ <= finish[6]*1.0000000000001; MQ += step[6])
					for(long double PMQ = start[7]; PMQ <= finish[7]*1.0000000000001; PMQ += step[7])
						for(long double n = start[8]; n <= finish[8]*1.0000000000001; n += step[8])
							for(long double Pn = start[9]; Pn <= finish[9]*1.0000000000001; Pn += step[9])
							{
								Parameter_List[i] = new long double[10];
								Parameter_List[i][0] = A;
								Parameter_List[i][1] = PA;
								Parameter_List[i][2] = M;
								Parameter_List[i][3] = PM;
								Parameter_List[i][4] = Gamma;
								Parameter_List[i][5] = PGamma;
								Parameter_List[i][6] = MQ;
								Parameter_List[i][7] = PMQ;
								Parameter_List[i][8] = n;
								Parameter_List[i][9] = Pn;
								i++;
							}

		for(i = Thread_Num; i < Dims[0]*Dims[1]*Dims[2]*Dims[3]*Dims[4]*Dims[5]*Dims[6]*Dims[7]*Dims[8]*Dims[9]; i += Num_Threads)
		{

			JPsi[Temp]->Replace(Parameter_List[i][0],0,1);
			JPsi[Temp]->Replace(Parameter_List[i][1],0,2);
			JPsi[Temp]->Replace(Parameter_List[i][2],1,1);
			JPsi[Temp]->Replace(Parameter_List[i][3],1,2);
			JPsi[Temp]->Replace(Parameter_List[i][4],2,1);
			JPsi[Temp]->Replace(Parameter_List[i][5],2,2);
			Non[Temp]->Replace(Parameter_List[i][6],0,1);
			Non[Temp]->Replace(Parameter_List[i][7],0,2);
			Non[Temp]->Replace(Parameter_List[i][8],1,1);
			Non[Temp]->Replace(Parameter_List[i][9],1,2);

			Print(JPsi[Temp], Psi_Prime[Temp], Non[Temp], Medium_Euclidean, Medium_Spatial, Temp);
		}
		return(0);
	}

	long double Best[15], Chi;
	OutputFile << "Random Search seed: " << time(NULL)+3*atoi(argv[1]) << endl;
	srand(time(NULL)+3*atoi(argv[2]));
	time_t start_time = time(NULL);
	if(atoi(argv[2])>100)
	{
		Best[0] = Uniform(Random_Range[0][0],Random_Range[0][1]);
		Best[1] = Uniform(Random_Range[1][0],Random_Range[1][1]);
		Best[2] = Uniform(Random_Range[2][0],Random_Range[2][1]);
		Best[3] = Uniform(Random_Range[3][0],Random_Range[3][1]);
		Best[4] = Uniform(Random_Range[4][0],Random_Range[4][1]);
		Best[5] = Uniform(Random_Range[5][0],Random_Range[5][1]);
		Best[6] = Uniform(Random_Range[6][0],Random_Range[6][1]);
		Best[7] = Uniform(Random_Range[7][0],Random_Range[7][1]);
		Best[8] = Uniform(Random_Range[8][0],Random_Range[8][1]);
		Best[9] = Uniform(Random_Range[9][0],Random_Range[9][1]);
		Best[10] = Uniform(Random_Range[10][0],Random_Range[10][1]);
		Best[11] = Uniform(Random_Range[11][0],Random_Range[11][1]);
		Best[12] = Uniform(Random_Range[12][0],Random_Range[12][1]);
		Best[13] = Uniform(Random_Range[13][0],Random_Range[13][1]);
		JPsi[Temp]->Replace(Best[0],0,1);
		JPsi[Temp]->Replace(Best[1],0,2);
		JPsi[Temp]->Replace(Best[2],1,1);
		JPsi[Temp]->Replace(Best[3],1,2);
		JPsi[Temp]->Replace(Best[4],2,1);
		JPsi[Temp]->Replace(Best[5],2,2);
		JPsi[Temp]->Replace(Best[6],3,1);
		JPsi[Temp]->Replace(Best[7],3,2);
		Psi_Prime[Temp]->Replace(Best[6],3,1);
		Psi_Prime[Temp]->Replace(Best[7],3,2);
		JPsi[Temp]->Replace(Best[8],4,1);
		JPsi[Temp]->Replace(Best[9],4,2);
		Non[Temp]->Replace(Best[10],0,1);
		Non[Temp]->Replace(Best[11],0,2);
		Non[Temp]->Replace(Best[12],1,1);
		Non[Temp]->Replace(Best[13],1,2);
	}
	else
	{
		Best[0] = JPsi[Temp]->Read(0,1);
		Best[1] = JPsi[Temp]->Read(0,2);
		Best[2] = JPsi[Temp]->Read(1,1);
		Best[3] = JPsi[Temp]->Read(1,2);
		Best[4] = JPsi[Temp]->Read(2,1);
		Best[5] = JPsi[Temp]->Read(2,2);
		Best[6] = JPsi[Temp]->Read(3,1);
		Best[7] = JPsi[Temp]->Read(3,2);
		Best[8] = JPsi[Temp]->Read(4,1);
		Best[9] = JPsi[Temp]->Read(4,2);
		Best[10] = Non[Temp]->Read(0,1);
		Best[11] = Non[Temp]->Read(0,2);
		Best[12] = Non[Temp]->Read(1,1);
		Best[13] = Non[Temp]->Read(1,2);
	}
	Medium_Euclidean[0] = JPsi[Temp]->Euclidean(.5,0)+Psi_Prime[Temp]->Euclidean(.5,0)+Non[Temp]->Euclidean(.5,0);
	Medium_Euclidean[1] = JPsi[Temp]->Euclidean(.5,3)+Psi_Prime[Temp]->Euclidean(.5,3)+Non[Temp]->Euclidean(.5,3);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = JPsi[Temp]->Spatial((long double)(j)+.25)+Psi_Prime[Temp]->Spatial((long double)(j)+.25)+Non[Temp]->Spatial((long double)(j)+.25);
	Best[14] = Print(JPsi[Temp], Psi_Prime[Temp], Non[Temp], Medium_Euclidean, Medium_Spatial, Temp);

	while(difftime(time(NULL), start_time) < 9000)
	{
		JPsi[Temp]->Replace(Protected_Uniform(Best[0],Random_Range[0][0],Random_Range[0][1],Best[14]),0,1);
		JPsi[Temp]->Replace(Protected_Uniform(Best[1],Random_Range[1][0],Random_Range[1][1],Best[14]),0,2);
		JPsi[Temp]->Replace(Protected_Uniform(Best[2],Random_Range[2][0],Random_Range[2][1],Best[14]),1,1);
		JPsi[Temp]->Replace(Protected_Uniform(Best[3],Random_Range[3][0],Random_Range[3][1],Best[14]),1,2);
		JPsi[Temp]->Replace(Protected_Uniform(Best[4],Random_Range[4][0],Random_Range[4][1],Best[14]),2,1);
		JPsi[Temp]->Replace(Protected_Uniform(Best[5],Random_Range[5][0],Random_Range[5][1],Best[14]),2,2);
		JPsi[Temp]->Replace(Protected_Uniform(Best[6],Random_Range[6][0],Random_Range[6][1],Best[14]),3,1);
		JPsi[Temp]->Replace(Protected_Uniform(Best[7],Random_Range[7][0],Random_Range[7][1],Best[14]),3,2);
		Psi_Prime[Temp]->Replace(JPsi[Temp]->Read(3,1),3,1);
		Psi_Prime[Temp]->Replace(JPsi[Temp]->Read(3,2),3,2);
		JPsi[Temp]->Replace(Protected_Uniform(Best[8],Random_Range[8][0],Random_Range[8][1],Best[14]),4,1);
		JPsi[Temp]->Replace(Protected_Uniform(Best[9],Random_Range[9][0],Random_Range[9][1],Best[14]),4,2);
		Non[Temp]->Replace(Protected_Uniform(Best[10],Random_Range[10][0],Random_Range[10][1],Best[14]),0,1);
		Non[Temp]->Replace(Protected_Uniform(Best[11],Random_Range[11][0],Random_Range[11][1],Best[14]),0,2);
		Non[Temp]->Replace(Protected_Uniform(Best[12],Random_Range[12][0],Random_Range[12][1],Best[14]),1,1);
		Non[Temp]->Replace(Protected_Uniform(Best[13],Random_Range[13][0],Random_Range[13][1],Best[14]),1,2);
		Medium_Euclidean[0] = JPsi[Temp]->Euclidean(.5,0)+Psi_Prime[Temp]->Euclidean(.5,0)+Non[Temp]->Euclidean(.5,0);
		Medium_Euclidean[1] = JPsi[Temp]->Euclidean(.5,3)+Psi_Prime[Temp]->Euclidean(.5,3)+Non[Temp]->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = JPsi[Temp]->Spatial((long double)(j)+.25)+Psi_Prime[Temp]->Spatial((long double)(j)+.25)+Non[Temp]->Spatial((long double)(j)+.25);
		Chi = Print(JPsi[Temp], Psi_Prime[Temp], Non[Temp], Medium_Euclidean, Medium_Spatial, Temp);

		if((!isnan(Chi) && Chi < Best[14]) || isnan(Best[14]))
		{
			Best[0] = JPsi[Temp]->Read(0,1);
			Best[1] = JPsi[Temp]->Read(0,2);
			Best[2] = JPsi[Temp]->Read(1,1);
			Best[3] = JPsi[Temp]->Read(1,2);
			Best[4] = JPsi[Temp]->Read(2,1);
			Best[5] = JPsi[Temp]->Read(2,2);
			Best[6] = JPsi[Temp]->Read(3,1);
			Best[7] = JPsi[Temp]->Read(3,2);
			Best[8] = JPsi[Temp]->Read(4,1);
			Best[9] = JPsi[Temp]->Read(4,2);
			Best[10] = Non[Temp]->Read(0,1);
			Best[11] = Non[Temp]->Read(0,2);
			Best[12] = Non[Temp]->Read(1,1);
			Best[13] = Non[Temp]->Read(1,2);
			Best[14] = Chi;
		}
	}
	JPsi[Temp]->Replace(Best[0],0,1);
	JPsi[Temp]->Replace(Best[1],0,2);
	JPsi[Temp]->Replace(Best[2],1,1);
	JPsi[Temp]->Replace(Best[3],1,2);
	JPsi[Temp]->Replace(Best[4],2,1);
	JPsi[Temp]->Replace(Best[5],2,2);
	JPsi[Temp]->Replace(Best[6],3,1);
	JPsi[Temp]->Replace(Best[7],3,2);
	Psi_Prime[Temp]->Replace(Best[6],3,1);
	Psi_Prime[Temp]->Replace(Best[7],3,2);
	JPsi[Temp]->Replace(Best[8],4,1);
	JPsi[Temp]->Replace(Best[9],4,2);
	Non[Temp]->Replace(Best[10],0,1);
	Non[Temp]->Replace(Best[11],0,2);
	Non[Temp]->Replace(Best[12],1,1);
	Non[Temp]->Replace(Best[13],1,2);

	long double gradn_1[14], gradn[14], sn_1[14], sn[14];
	long double betan = 0, betan_1;
	start_time = time(NULL);

	Gradient(gradn_1, JPsi[Temp], Psi_Prime[Temp], Non[Temp], Temp);
	for(int i = 0; i < 14; i++)
		sn_1[i] = gradn_1[i];
	Minimize(sn_1, JPsi[Temp], Psi_Prime[Temp], Non[Temp], Temp);
	do
	{
		betan_1 = betan;
		Gradient(gradn, JPsi[Temp], Psi_Prime[Temp], Non[Temp], Temp);
		betan = PolakRibiere(gradn_1,gradn);
		if(betan < 0)
			for(int i = 0; i < 14; i++)
				sn[i] = gradn[i]+betan*sn_1[i];
		else
			for(int i = 0; i < 14; i++)
				sn[i] = gradn[i];
		Minimize(sn, JPsi[Temp], Psi_Prime[Temp], Non[Temp], Temp);
		for(int i = 0; i < 14; i++)
		{
			gradn_1[i] = gradn[i];
			sn_1[i] = sn[i];
		}
	}while((abs(betan)>1e-12 || abs(betan_1)>1e-12) && time(NULL)-start_time < 9000);

	OutputFile.close();

	return(0);
}

void Minimize(long double sn[14], Spectral_Inter* JPsi, Spectral_Inter* Psi_Prime, Spectral_Non* Non, int Temp)
{
	long double JPsi_Local[5][3];
	long double Psi_Prime_Local[5][3];
	long double Non_Local[2][3];
	long double Medium_Euclidean[2];
	long double Medium_Spatial[7];
	long double length[2][2];

	for(int i = 0; i < 5; i++)
		for(int j = 0; j < 3; j++)
		{
			JPsi_Local[i][j] = JPsi->Read(i,j);
			Psi_Prime_Local[i][j] = Psi_Prime->Read(i,j);
			if(i < 2)
				Non_Local[i][j] = Non->Read(i,j);
		}

	length[0][0] = (sn[0]>0)?(Random_Range[0][1]-JPsi_Local[0][1])/sn[0]:(Random_Range[0][0]-JPsi_Local[0][1])/sn[0];
	length[0][1] = (sn[0]>0)?(Random_Range[0][0]-JPsi_Local[0][1])/sn[0]:(Random_Range[0][1]-JPsi_Local[0][1])/sn[0];
	length[1][0] = (sn[1]>0)?(Random_Range[1][1]-JPsi_Local[0][2])/sn[1]:(Random_Range[1][0]-JPsi_Local[0][2])/sn[1];
	length[1][1] = (sn[1]>0)?(Random_Range[1][0]-JPsi_Local[0][2])/sn[1]:(Random_Range[1][1]-JPsi_Local[0][2])/sn[1];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[2]>0)?(Random_Range[2][1]-JPsi_Local[1][1])/sn[2]:(Random_Range[2][0]-JPsi_Local[1][1])/sn[2];
	length[1][1] = (sn[2]>0)?(Random_Range[2][0]-JPsi_Local[1][1])/sn[2]:(Random_Range[2][1]-JPsi_Local[1][1])/sn[2];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[3]>0)?(Random_Range[3][1]-JPsi_Local[1][2])/sn[3]:(Random_Range[3][0]-JPsi_Local[1][2])/sn[3];
	length[1][1] = (sn[3]>0)?(Random_Range[3][0]-JPsi_Local[1][2])/sn[3]:(Random_Range[3][1]-JPsi_Local[1][2])/sn[3];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[4]>0)?(Random_Range[4][1]-JPsi_Local[2][1])/sn[4]:(Random_Range[4][0]-JPsi_Local[2][1])/sn[4];
	length[1][1] = (sn[4]>0)?(Random_Range[4][0]-JPsi_Local[2][1])/sn[4]:(Random_Range[4][1]-JPsi_Local[2][1])/sn[4];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[5]>0)?(Random_Range[5][1]-JPsi_Local[2][2])/sn[5]:(Random_Range[5][0]-JPsi_Local[2][2])/sn[5];
	length[1][1] = (sn[5]>0)?(Random_Range[5][0]-JPsi_Local[2][2])/sn[5]:(Random_Range[5][1]-JPsi_Local[2][2])/sn[5];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	/*length[1][0] = (sn[6]>0)?(Random_Range[6][1]-JPsi_Local[3][1])/sn[6]:(Random_Range[6][0]-JPsi_Local[3][1])/sn[6];
	length[1][1] = (sn[6]>0)?(Random_Range[6][0]-JPsi_Local[3][1])/sn[6]:(Random_Range[6][1]-JPsi_Local[3][1])/sn[6];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[7]>0)?(Random_Range[7][1]-JPsi_Local[3][2])/sn[7]:(Random_Range[7][0]-JPsi_Local[3][2])/sn[7];
	length[1][1] = (sn[7]>0)?(Random_Range[7][0]-JPsi_Local[3][2])/sn[7]:(Random_Range[7][1]-JPsi_Local[3][2])/sn[7];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[8]>0)?(Random_Range[8][1]-JPsi_Local[4][1])/sn[8]:(Random_Range[8][0]-JPsi_Local[4][1])/sn[8];
	length[1][1] = (sn[8]>0)?(Random_Range[8][0]-JPsi_Local[4][1])/sn[8]:(Random_Range[8][1]-JPsi_Local[4][1])/sn[8];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[9]>0)?(Random_Range[9][1]-JPsi_Local[4][2])/sn[9]:(Random_Range[9][0]-JPsi_Local[4][2])/sn[9];
	length[1][1] = (sn[9]>0)?(Random_Range[9][0]-JPsi_Local[4][2])/sn[9]:(Random_Range[9][1]-JPsi_Local[4][2])/sn[9];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];*/
	length[1][0] = (sn[10]>0)?(Random_Range[10][1]-Non_Local[0][1])/sn[10]:(Random_Range[10][0]-Non_Local[0][1])/sn[10];
	length[1][1] = (sn[10]>0)?(Random_Range[10][0]-Non_Local[0][1])/sn[10]:(Random_Range[10][1]-Non_Local[0][1])/sn[10];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[11]>0)?(Random_Range[11][1]-Non_Local[0][2])/sn[11]:(Random_Range[11][0]-Non_Local[0][2])/sn[11];
	length[1][1] = (sn[11]>0)?(Random_Range[11][0]-Non_Local[0][2])/sn[11]:(Random_Range[11][1]-Non_Local[0][2])/sn[11];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[12]>0)?(Random_Range[12][1]-Non_Local[1][1])/sn[12]:(Random_Range[12][0]-Non_Local[1][1])/sn[12];
	length[1][1] = (sn[12]>0)?(Random_Range[12][0]-Non_Local[1][1])/sn[12]:(Random_Range[12][1]-Non_Local[1][1])/sn[12];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[13]>0)?(Random_Range[13][1]-Non_Local[1][2])/sn[13]:(Random_Range[13][0]-Non_Local[1][2])/sn[13];
	length[1][1] = (sn[13]>0)?(Random_Range[13][0]-Non_Local[1][2])/sn[13]:(Random_Range[13][1]-Non_Local[1][2])/sn[13];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];

	long double fz[101][2];
	long double a = 0, c = 1;
	int Min_i;

	OutputFile << "Line search" << endl;
	for(int i = 0; i <= 100; i++)
	{
		fz[i][0] = a+(c-a)*i/100.;
		JPsi->Replace(JPsi_Local[0][1]+length[0][0]*fz[i][0]*sn[0],0,1);
		JPsi->Replace(JPsi_Local[0][2]+length[0][0]*fz[i][0]*sn[0],0,2);
		JPsi->Replace(JPsi_Local[1][1]+length[0][0]*fz[i][0]*sn[0],1,1);
		JPsi->Replace(JPsi_Local[1][2]+length[0][0]*fz[i][0]*sn[0],1,2);
		JPsi->Replace(JPsi_Local[2][1]+length[0][0]*fz[i][0]*sn[0],2,1);
		JPsi->Replace(JPsi_Local[2][2]+length[0][0]*fz[i][0]*sn[0],2,2);
		/*JPsi->Replace(JPsi_Local[3][1]+length[0][0]*fz[i][0]*sn[0],3,1);
		JPsi->Replace(JPsi_Local[3][2]+length[0][0]*fz[i][0]*sn[0],3,2);
		Psi_Prime->Replace(Psi_Prime_Local[3][1]+length[0][0]*fz[i][0]*sn[0],3,1);
		Psi_Prime->Replace(Psi_Prime_Local[3][2]+length[0][0]*fz[i][0]*sn[0],3,2);
		JPsi->Replace(JPsi_Local[4][1]+length[0][0]*fz[i][0]*sn[0],4,1);
		JPsi->Replace(JPsi_Local[4][2]+length[0][0]*fz[i][0]*sn[0],4,2);*/
		Non->Replace(Non_Local[0][1]+length[0][0]*fz[i][0]*sn[0],0,1);
		Non->Replace(Non_Local[0][2]+length[0][0]*fz[i][0]*sn[0],0,2);
		Non->Replace(Non_Local[1][1]+length[0][0]*fz[i][0]*sn[0],1,1);
		Non->Replace(Non_Local[1][2]+length[0][0]*fz[i][0]*sn[0],1,2);
		Medium_Euclidean[0] = JPsi->Euclidean(.5,0)+Psi_Prime->Euclidean(.5,0)+Non->Euclidean(.5,0);
		Medium_Euclidean[1] = JPsi->Euclidean(.5,3)+Psi_Prime->Euclidean(.5,3)+Non->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = JPsi->Spatial((long double)(j)+.25)+Psi_Prime->Spatial((long double)(j)+.25)+Non->Spatial((long double)(j)+.25);
		fz[i][1] = Print(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, Temp);
	}

	for(int i = 1; i < 101; i++)
		Min_i = (fz[i][1]<fz[Min_i][1])?i:Min_i;
	if(Min_i == 0)	//So as to not exceed the limits of the array
		Min_i = 1;
	else if(Min_i == 200)
		Min_i = 99;

	OutputFile << "Brent's Method" << endl;
	const int ITMAX = 100;
	const long double GOLDEN = (3.-sqrt(5.))/2.;
	const long double ZEPS=LDBL_EPSILON*1e-3;
	long double b,d=0,etemp,fu,fv,fw,fx;
	long double p,q,r,tol1,tol2,u,v,w,x,xm;
	long double e=0;
	const long double tol = 3e-16;

	a = fz[Min_i-1][0];
	b = fz[Min_i+1][0];
	x = w = v = fz[Min_i][0];
	fx = fw = fv = fz[Min_i][1];

	for(int inter = 0; inter < ITMAX; inter++)
	{
		xm = (a+b)/2.;
		tol1 = tol*abs(x)+ZEPS;
		tol2 = 2.*tol1;
		if(abs(x-xm) <= (tol2-(b-a)/2.) || abs(fx/fu-1.) < 1e-7)
		{
			JPsi->Replace(JPsi_Local[0][1]+length[0][0]*x*sn[0],0,1);
			JPsi->Replace(JPsi_Local[0][2]+length[0][0]*x*sn[0],0,2);
			JPsi->Replace(JPsi_Local[1][1]+length[0][0]*x*sn[0],1,1);
			JPsi->Replace(JPsi_Local[1][2]+length[0][0]*x*sn[0],1,2);
			JPsi->Replace(JPsi_Local[2][1]+length[0][0]*x*sn[0],2,1);
			JPsi->Replace(JPsi_Local[2][2]+length[0][0]*x*sn[0],2,2);
			/*JPsi->Replace(JPsi_Local[3][1]+length[0][0]*x*sn[0],3,1);
			JPsi->Replace(JPsi_Local[3][2]+length[0][0]*x*sn[0],3,2);
			Psi_Prime->Replace(Psi_Prime_Local[3][1]+length[0][0]*x*sn[0],3,1);
			Psi_Prime->Replace(Psi_Prime_Local[3][2]+length[0][0]*x*sn[0],3,2);
			JPsi->Replace(JPsi_Local[4][1]+length[0][0]*x*sn[0],4,1);
			JPsi->Replace(JPsi_Local[4][2]+length[0][0]*x*sn[0],4,2);*/
			Non->Replace(Non_Local[0][1]+length[0][0]*x*sn[0],0,1);
			Non->Replace(Non_Local[0][2]+length[0][0]*x*sn[0],0,2);
			Non->Replace(Non_Local[1][1]+length[0][0]*x*sn[0],1,1);
			Non->Replace(Non_Local[1][2]+length[0][0]*x*sn[0],1,2);
			Medium_Euclidean[0] = JPsi->Euclidean(.5,0)+Psi_Prime->Euclidean(.5,0)+Non->Euclidean(.5,0);
			Medium_Euclidean[1] = JPsi->Euclidean(.5,3)+Psi_Prime->Euclidean(.5,3)+Non->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[j] = JPsi->Spatial((long double)(j)+.25)+Psi_Prime->Spatial((long double)(j)+.25)+Non->Spatial((long double)(j)+.25);
			Print(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, Temp);
			return;
		}

		if(abs(e) > tol1)
		{
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q-(x-w)*r;
			q = 2.*(q-r);
			if(q > 0) p = -p;
			q = abs(q);
			etemp = e;
			e = d;
			if(abs(p) >= abs(q*etemp/2.) || p <= q*(a-x) || p >= q*(b-x))
				d = GOLDEN*(e=(x>=xm?a-x:b-x));
			else
			{
				d = p/q;
				u = x+d;
				if(u-a < tol2 || b-u < tol2)
					d = tol1*(xm-x)/abs(xm-x);
			}
		}
		else
			d = GOLDEN*(e=(x>=xm?a-x:b-x));

		u = (abs(d)>tol1?x+d:x+tol1*d/abs(d));
		JPsi->Replace(JPsi_Local[0][1]+length[0][0]*u*sn[0],0,1);
		JPsi->Replace(JPsi_Local[0][2]+length[0][0]*u*sn[0],0,2);
		JPsi->Replace(JPsi_Local[1][1]+length[0][0]*u*sn[0],1,1);
		JPsi->Replace(JPsi_Local[1][2]+length[0][0]*u*sn[0],1,2);
		JPsi->Replace(JPsi_Local[2][1]+length[0][0]*u*sn[0],2,1);
		JPsi->Replace(JPsi_Local[2][2]+length[0][0]*u*sn[0],2,2);
		/*JPsi->Replace(JPsi_Local[3][1]+length[0][0]*u*sn[0],3,1);
		JPsi->Replace(JPsi_Local[3][2]+length[0][0]*u*sn[0],3,2);
		Psi_Prime->Replace(Psi_Prime_Local[3][1]+length[0][0]*u*sn[0],3,1);
		Psi_Prime->Replace(Psi_Prime_Local[3][2]+length[0][0]*u*sn[0],3,2);
		JPsi->Replace(JPsi_Local[4][1]+length[0][0]*u*sn[0],4,1);
		JPsi->Replace(JPsi_Local[4][2]+length[0][0]*u*sn[0],4,2);*/
		Non->Replace(Non_Local[0][1]+length[0][0]*u*sn[0],0,1);
		Non->Replace(Non_Local[0][2]+length[0][0]*u*sn[0],0,2);
		Non->Replace(Non_Local[1][1]+length[0][0]*u*sn[0],1,1);
		Non->Replace(Non_Local[1][2]+length[0][0]*u*sn[0],1,2);
		Medium_Euclidean[0] = JPsi->Euclidean(.5,0)+Psi_Prime->Euclidean(.5,0)+Non->Euclidean(.5,0);
		Medium_Euclidean[1] = JPsi->Euclidean(.5,3)+Psi_Prime->Euclidean(.5,3)+Non->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = JPsi->Spatial((long double)(j)+.25)+Psi_Prime->Spatial((long double)(j)+.25)+Non->Spatial((long double)(j)+.25);
		fu = Print(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, Temp);

		if(fu <= fx)
		{
			if(u >= x)
				a = x;
			else
				b = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else
		{
			if(u < x)
				a = u;
			else
				b = u;

			if(fu <= fw || w == x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if(fu <= fv || v == x|| v == w)
			{
				v = u;
				fv = fu;
			}
		}
	}
}

long double PolakRibiere(long double gradn_1[14], long double gradn[14])
{
	long double num = 0, den = 0;

	for(int i = 0; i < 14; i++)
	{
		num += gradn[i]*(gradn[i]-gradn_1[i]);
		den += pow(gradn_1[i],2);
	}

	return(num/den);
}

void Gradient(long double grad[14], Spectral_Inter* JPsi, Spectral_Inter* Psi_Prime, Spectral_Non* Non, int T)
{
	long double f0, f1, h = 1e-5;
	long double Reduce_Euclidean[3][2];
	long double Reduce_Spatial[3][7];
	long double Medium_Euclidean[2];
	long double Medium_Spatial[7];
	long double norm = 0;
	int i, j;

	OutputFile << "Gradient" << endl;
	Reduce_Euclidean[0][0] = JPsi->Euclidean(.5, 0);
	Reduce_Euclidean[1][0] = Psi_Prime->Euclidean(.5, 0);
	Reduce_Euclidean[2][0] = Non->Euclidean(.5, 0);
	Medium_Euclidean[0] = Reduce_Euclidean[0][0]+Reduce_Euclidean[1][0]+Reduce_Euclidean[2][0];
	Reduce_Euclidean[0][1] = JPsi->Euclidean(.5, 3);
	Reduce_Euclidean[1][1] = Psi_Prime->Euclidean(.5, 3);
	Reduce_Euclidean[2][1] = Non->Euclidean(.5, 3);
	Medium_Euclidean[1] = Reduce_Euclidean[0][1]+Reduce_Euclidean[1][1]+Reduce_Euclidean[2][1];
	for(int j = 0; j < 7; j++)
	{
		Reduce_Spatial[0][j] = JPsi->Spatial((long double)(j)+.25);
		Reduce_Spatial[1][j] = Psi_Prime->Spatial((long double)(j)+.25);
		Reduce_Spatial[2][j] = Non->Spatial((long double)(j)+.25);
		Medium_Spatial[j] = Reduce_Spatial[0][j]+Reduce_Spatial[1][j]+Reduce_Spatial[2][j];
	}
	f0 = Chi_Square(Medium_Euclidean, Medium_Spatial, T);
	//cerr << "Gradient " << f0 << endl;

	for(i = 0; i < 5; i++)
	{
		for(j = 1; j < 3; j++)
		{
			JPsi->Add(h,i,j);
			if(i == 3)
			{
				Psi_Prime->Add(h,i,j);
				Medium_Euclidean[0] = JPsi->Euclidean(.5, 0)+Psi_Prime->Euclidean(.5, 0)+Reduce_Euclidean[2][0];
				Medium_Euclidean[1] = JPsi->Euclidean(.5, 3)+Psi_Prime->Euclidean(.5, 3)+Reduce_Euclidean[2][1];
				for(int j = 0; j < 7; j++)
					Medium_Spatial[j] = JPsi->Spatial((long double)(j)+.25)+JPsi->Spatial((long double)(j)+.25)+Reduce_Spatial[2][j];
			}
			else
			{
				Medium_Euclidean[0] = JPsi->Euclidean(.5, 0)+Reduce_Euclidean[1][0]+Reduce_Euclidean[2][0];
				Medium_Euclidean[1] = JPsi->Euclidean(.5, 3)+Reduce_Euclidean[1][1]+Reduce_Euclidean[2][1];
				for(int j = 0; j < 7; j++)
					Medium_Spatial[j] = JPsi->Spatial((long double)(j)+.25)+Reduce_Spatial[1][j]+Reduce_Spatial[2][j];
			}
			f1 = Print(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, T);
			JPsi->Add(-h,i,j);
			if(i == 3)
				Psi_Prime->Add(-h,i,j);
			grad[2*i+j-1] = (f0-f1)/h;
		}
	}

	for(i = 0; i < 2; i++)
	{
		for(j = 1; j < 3; j++)
		{
			Non->Add(h,i,j);
			Medium_Euclidean[0] = Reduce_Euclidean[0][0]+Reduce_Euclidean[1][0]+Non->Euclidean(.5, 0);
			Medium_Euclidean[1] = Reduce_Euclidean[0][1]+Reduce_Euclidean[1][1]+Non->Euclidean(.5, 3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[j] = Reduce_Spatial[0][j]+Reduce_Spatial[1][j]+Non->Spatial((long double)(j)+.25);
			f1 = Print(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, T);
			Non->Add(-h,i,j);
			grad[2*i+j+9] = (f0-f1)/h;
		}
	}
}

long double Print(Spectral_Inter* JPsi, Spectral_Inter* Psi_Prime, Spectral_Non* Non, long double Medium_Euclidean[2], long double Medium_Spatial[7], int T)
{
	for(int j = 0; j < 5; j++)
		OutputFile << JPsi->Read(j,1) << "," << JPsi->Read(j,2) << "," << flush;
	for(int j = 0; j < 2; j++)
		OutputFile << Non->Read(j,1) << "," << Non->Read(j,2) << "," << flush;
	OutputFile << Medium_Euclidean[1]/Vacuum_Euclidean[T-1][1]-Medium_Euclidean[0]/Vacuum_Euclidean[T-1][0] << "," << flush;
	for(int j = 0; j < 7; j++)
	{
		OutputFile << Medium_Spatial[j]/Vacuum_Spatial[j] << "," << flush;
	}
	long double Chi = Chi_Square(Medium_Euclidean, Medium_Spatial, T);
	OutputFile << Chi << endl;
	return(Chi);
}

long double Chi_Square(long double Medium_Euclidean[2], long double Medium_Spatial[7], int T)
{
	long double answer;
	answer = pow(Medium_Euclidean[1]/Vacuum_Euclidean[T-1][1]-Medium_Euclidean[0]/Vacuum_Euclidean[T-1][0]-.2,2)/.2;
	for(int i = 0; i < 7; i++)
	{
		answer += pow(Medium_Spatial[i]/Vacuum_Spatial[i]-Spatial_Ratio[T-1][i],2)/Spatial_Ratio[T-1][i];
		if(i < 6)answer += pow((Medium_Spatial[i+1]/Vacuum_Spatial[i+1]-Medium_Spatial[i]/Vacuum_Spatial[i]>0)?(Medium_Spatial[i+1]/Vacuum_Spatial[i+1]-Medium_Spatial[i]/Vacuum_Spatial[i]):0,.5);
	}
	return(answer);
}
long double Protected_Uniform(long double x0, long double a, long double b, long double chi)
{
	long double Delta = (b-a)*pow(chi,.6)/10.;
	if(x0-Delta < a && x0+Delta > b)
		return(Uniform(a,b));
	else if(x0-Delta < a)
		return(Uniform(a,x0+Delta));
	else if(x0+Delta > b)
		return(Uniform(x0-Delta,b));
	return(x0+Uniform(-Delta,Delta));
}

long double Uniform(long double a, long double b)
{
	return((long double)(rand())/(long double)(RAND_MAX)*(b-a)+a);
}

void mergeSort(long double List[], int a, int b)
{
	int i, j, k;
	long double Temp[(a+b)/2-a+1];

	if(b-a > 1)	//Divide
	{
		mergeSort(List, a, (a+b)/2);
		mergeSort(List, (a+b)/2+1, b);
	}

	for(i = 0; i <= (a+b)/2-a; i++)	//Copy out the lower half array in prep for copy over
		Temp[i] = List[i+a];

	j = 0;
	k = (a+b)/2+1;
	for(i = a; i <= b && j <= (a+b)/2-a && k <= b; i++)	//Merge/Conqure while both half lists have not been exhausted
	{
		if(Temp[j] <= List[k])
		{
			List[i] = Temp[j];
			j++;
		}
		else
		{
			List[i] = List[k];
			k++;
		}
	}
	for(; i <= b && j <= (a+b)/2-a; i++)	//If the Temp list has not been exhausted, complete the job
	{
		List[i] = Temp[j];
		j++;
	}

	return;
}
