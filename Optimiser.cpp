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

long double SpectralNon(long double, long double, long double[2][3], bool);

void Gradient(long double[14], long double[5][3], long double[5][3], long double[2][3], long double[2], long double[7], long double[7], long double);
long double PolakRibiere(long double[14], long double[14]);
void Minimize(long double[14], long double[5][3], long double[5][3], long double[2][3], long double[2], long double[7], long double[7], long double);

long double Uniform(long double, long double);
long double Chi_Square(long double[2], long double[2], long double[7], long double[7], long double[7]);
long double Least_Squares(long double, long double, long double, long double);
long double Print(long double[5][3], long double[5][3], long double[2][3], long double[2], long double[2], long double[7], long double[7], long double[7]); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square

long double Boundary[] = {0.00865, 0.0267, 0.0491, 0.0985, .421, .802, 1.01, 4.85};
long double Random_Range[14][2] = {{.1,.5},{1.,6.},{1.,3.5},{1.,6.},{.02,.18},{1.,6.},{5.,15.},{1.,6.},{1.5,3.5},{1.,6.},{1.59,1.79},{1.,6.},{1.5,5.},{1.,6.}};
ofstream OutputFile;
/*{0.296663, 2.99912, 2.41053, 5.50001, 0.168958, 3.60003, 10.58, 4.57122, 3, 5.25, 1.69802, 4.00006, 3.00002, 5.49997, 0.181899, 1.00337, 1.00995, 0.995169, 0.985732, 0.980595, 0.973869, 0.956278, 0.00181892} T=194 MeV Min
{0.3, 3, 2.45, 5.5, 0.18, 3.6, 10.58, 4.57122, 3, 5.25, 1.7, 5.5, 2.25, 6, 0.169186, 1.00663, 1.03674, 1.00002, 0.997114, 0.999272, 0.946506, 0.939859, 0.00737184} T=194 MeV Grid
{0.445826, 3.86826, 2.58405, 4.01049, 0.535639, 4.20052, 10.3425, 2.64938, 2.81697, 2.12552, 1.69153, 3.42941, 2.26904, 4.4543, 0.159886, 1.00841, 1.06611, 0.933118, 0.803052, 0.778368, 0.74589, 0.664434, 0.0225266} T=258 Min
{0.32, 3, 2.8, 3.6, 0.033, 4.9, 8.97, 1, 3, 1, 1.55, 5.5, 2.25, 3, 0.131497, 1.0083, 0.958267, 0.97533, 0.861404, 0.783244, 0.695317, 0.692414, 0.0298499} T=258 Grid
{0.184891, 3.32152, 1.3123, 5.23093, 0.064861, 5.80902, 5.03399, 5.12164, 1.55583, 5.24328, 1.64566, 5.31292, 3.93386, 5.40511, 0.143769, 1.00805, 0.864955, 0.847608, 0.667474, 0.619708, 0.473687, 0.376542, 0.0332863} T=320 MeV Min
{0.33, 2.4, 2.87, 2.5, 0.084, 5.5, 7.93198, 4.11694, 3.05149, 4.4995, 1.7, 5.5, 3, 3.5, 0.0622976, 1.0064, 0.936273, 0.864206, 0.728055, 0.50826, 0.476024, 0.400429, 0.105599} T=320 MeV Grid
{0.529109, 4.49175, 2.67237, 2.20553, 0.165846, 5.87198, 6.08157, 2.67594, 2.80637, 3.47611, 1.56825, 3.52402, 1.78377, 3.60395, 0.101371, 1.00725, 1.02967, 0.875425, 0.402078, 0.00803392, 0.235885, 0.255533, 0.519549} T=400 MeV Min
{0.19, 2.1, 2.58, 2.1, 0.104, 3.9, 9.52, 1, 3, 1, 1.36, 3.09, 1.98, 5.5, 0.055641, 1.00304, 1.05204, 1.17194, 0.951777, 0.290306, 0.0942507, 0.0513963, 0.972516} T=400 MeV Grid*/
int main(int argc, char* argv[])
{
	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.29663,2.99912},{3.09946,2.41053,5.50001},{.106597,.168958,3.60003},{10.58,10.58,4.57122},{3,3,5.25}},
						{{.4024,.445826,3.86826},{3.125,2.58405,4.01049},{.372,.535639,4.20052},{8.97,10.3425,2.64938},{3,2.81697,2.12552}},
						{{.403047,.184891,3.32152},{3.151,1.3123,5.23093},{.68,.064861,5.80902},{9.07,5.03399,5.12164},{3,1.55583,5.24328}},
						{{.266834,.529109,4.49175},{3.1855,2.67237,2.20553},{.98,.165846,5.87198},{9.52,6.08157,2.67594},{3,2.80637,3.47611}}};
	long double PsiPrime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,10.58,4.57122},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.69802,4.00006},{2.45,3.00002,5.49997}},
					       {{1.59,1.69153,3.42941},{2.7,2.26904,4.4543}},
					       {{1.51,1.64566,5.31292},{2.4,3.93386,5.40511}},
					       {{1.36,1.56825,3.52402},{1.98,1.78377,3.60395}}};

	Spectral_Inter* JPsi[5];
	Spectral_Inter* Psi_Prime[5];
	Spectral_Non* JPsi[5];
	for(int i = 0; i < 5; i++)
	{
		JPsi[i] = new Spectral_Inter(JPsi_Parameters, i, bool(i));
		Psi_Prime[i] = new Spectral_Inter(PsiPrime_Parameters, i, bool(i));
		Non[i] = new Spectral_Non(Non_Parameters, i, bool(i));
	}

	long double Spatial_Ratio[4][7] = {{1.,1.00006,0.99883,0.992039,0.982366,0.970341,0.95766},
					   {.99,0.988286,0.945063,0.879461,0.798659,0.7259,0.654381},
					   {.98,0.954875,0.856416,0.720447,0.573465,0.45867,0.376707},
					   {.97,0.908029,0.715435,0.524036,0.372788,0.246218,0.18}};
	long double Vacuum_Spatial[7] = {13.5965519695970589115, 0.0415680222792943527448, 0.00120126774580634698112, 4.65000126560610055077e-05, 
					  1.96861980208214073368e-06, 8.66667213197805982178e-08, 3.94190650826746653845e-09};
	long double Vacuum_Euclidean[4][2] = {{0.000259568923526295945662, 9.64220418032793835491e-06},
					       {0.00219476484784199273204, 0.000188143003166908006425},
					       {0.00821024781558290065043, 0.00116410367128954708042},
					       {0.0264690797013430510705, 0.00575549843590767333436}};
	long double Medium_Spatial[7];
	long double Medium_Euclidean[2];
	/*for(int i = 0; i < 6; i++)	//Superceeded by precalculated values, standing by if services required
		Vacuum_Spatial[i] = Spatial((long double)(i)+.25, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[0][0] = Euclidean(1./.388, .194, 0, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[1][0] = Euclidean(1./.516, .258, 0, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[2][0] = Euclidean(1./.640, .320, 0, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[3][0] = Euclidean(1./.800, .400, 0, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[0][1] = Euclidean(1./.388, .194, 3, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[1][1] = Euclidean(1./.516, .258, 3, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[2][1] = Euclidean(1./.640, .320, 3, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);
	Vacuum_Euclidean[3][1] = Euclidean(1./.800, .400, 3, JPsi_Parameters[0], PsiPrime_Parameters[0], Non_Parameters[0], true);*/

	char File[70] = "data/Optimiser_Output";
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
			JPsi_Parameters[Temp][0][1] = Parameter_List[i][0];
			JPsi_Parameters[Temp][0][2] = Parameter_List[i][1];
			JPsi_Parameters[Temp][1][1] = Parameter_List[i][2];
			JPsi_Parameters[Temp][1][2] = Parameter_List[i][3];
			JPsi_Parameters[Temp][2][1] = Parameter_List[i][4];
			JPsi_Parameters[Temp][2][2] = Parameter_List[i][5];
			Non_Parameters[Temp][0][1] = Parameter_List[i][6];
			Non_Parameters[Temp][0][2] = Parameter_List[i][7];
			Non_Parameters[Temp][1][1] = Parameter_List[i][8];
			Non_Parameters[Temp][1][2] = Parameter_List[i][9];

			Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
			Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
			Print(JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Medium_Euclidean, Vacuum_Euclidean[Temp-1], Medium_Spatial, Vacuum_Spatial, Spatial_Ratio[Temp-1]);
		}
		return(0);
	}

	long double Best[15], Chi;
	OutputFile << "Random Search seed: " << time(NULL)+3*atoi(argv[1]) << endl;
	srand(time(NULL)+3*atoi(argv[2]));
	time_t start_time = time(NULL);
	if(atoi(argv[2])!=0)
	{
		Best[0] = JPsi_Parameters[Temp][0][1] = Uniform(Random_Range[0][0],Random_Range[0][1]);
		Best[1] = JPsi_Parameters[Temp][0][2] = Uniform(Random_Range[1][0],Random_Range[1][1]);
		Best[2] = JPsi_Parameters[Temp][1][1] = Uniform(Random_Range[2][0],Random_Range[2][1]);
		Best[3] = JPsi_Parameters[Temp][1][2] = Uniform(Random_Range[3][0],Random_Range[3][1]);
		Best[4] = JPsi_Parameters[Temp][2][1] = Uniform(Random_Range[4][0],Random_Range[4][1]);
		Best[5] = JPsi_Parameters[Temp][2][2] = Uniform(Random_Range[5][0],Random_Range[5][1]);
		Best[6] = JPsi_Parameters[Temp][3][1] = Uniform(Random_Range[6][0],Random_Range[6][1]);
		Best[7] = JPsi_Parameters[Temp][3][2] = Uniform(Random_Range[7][0],Random_Range[7][1]);
		Best[8] = JPsi_Parameters[Temp][4][1] = Uniform(Random_Range[8][0],Random_Range[8][1]);
		Best[9] = JPsi_Parameters[Temp][4][2] = Uniform(Random_Range[9][0],Random_Range[9][1]);
		Best[10] = Non_Parameters[Temp][0][1] = Uniform(Random_Range[10][0],Random_Range[10][1]);
		Best[11] = Non_Parameters[Temp][0][2] = Uniform(Random_Range[11][0],Random_Range[11][1]);
		Best[12] = Non_Parameters[Temp][1][1] = Uniform(Random_Range[12][0],Random_Range[12][1]);
		Best[13] = Non_Parameters[Temp][1][2] = Uniform(Random_Range[13][0],Random_Range[13][1]);
	}
	else
	{
		Best[0] = JPsi_Parameters[Temp][0][1];
		Best[1] = JPsi_Parameters[Temp][0][2];
		Best[2] = JPsi_Parameters[Temp][1][1];
		Best[3] = JPsi_Parameters[Temp][1][2];
		Best[4] = JPsi_Parameters[Temp][2][1];
		Best[5] = JPsi_Parameters[Temp][2][2];
		Best[6] = JPsi_Parameters[Temp][3][1];
		Best[7] = JPsi_Parameters[Temp][3][2];
		Best[8] = JPsi_Parameters[Temp][4][1];
		Best[9] = JPsi_Parameters[Temp][4][2];
		Best[10] = Non_Parameters[Temp][0][1];
		Best[11] = Non_Parameters[Temp][0][2];
		Best[12] = Non_Parameters[Temp][1][1];
		Best[13] = Non_Parameters[Temp][1][2];
	}
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
	Best[14] = Print(JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Medium_Euclidean, Vacuum_Euclidean[Temp-1], Medium_Spatial, Vacuum_Spatial, Spatial_Ratio[Temp-1]);

	while(difftime(time(NULL), start_time) < 9000)
	{
		JPsi_Parameters[Temp][0][1] = Uniform(Random_Range[0][0],Random_Range[0][1]);
		JPsi_Parameters[Temp][0][2] = Uniform(Random_Range[1][0],Random_Range[1][1]);
		JPsi_Parameters[Temp][1][1] = Uniform(Random_Range[2][0],Random_Range[2][1]);
		JPsi_Parameters[Temp][1][2] = Uniform(Random_Range[3][0],Random_Range[3][1]);
		JPsi_Parameters[Temp][2][1] = Uniform(Random_Range[4][0],Random_Range[4][1]);
		JPsi_Parameters[Temp][2][2] = Uniform(Random_Range[5][0],Random_Range[5][1]);
		JPsi_Parameters[Temp][3][1] = Uniform(Random_Range[6][0],Random_Range[6][1]);
		JPsi_Parameters[Temp][3][2] = Uniform(Random_Range[7][0],Random_Range[7][1]);
		JPsi_Parameters[Temp][4][1] = Uniform(Random_Range[8][0],Random_Range[8][1]);
		JPsi_Parameters[Temp][4][2] = Uniform(Random_Range[9][0],Random_Range[9][1]);
		Non_Parameters[Temp][0][1] = Uniform(Random_Range[10][0],Random_Range[10][1]);
		Non_Parameters[Temp][0][2] = Uniform(Random_Range[11][0],Random_Range[11][1]);
		Non_Parameters[Temp][1][1] = Uniform(Random_Range[12][0],Random_Range[12][1]);
		Non_Parameters[Temp][1][2] = Uniform(Random_Range[13][0],Random_Range[13][1]);
		Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
		Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);

		Chi = Print(JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Medium_Euclidean, Vacuum_Euclidean[Temp-1], Medium_Spatial, Vacuum_Spatial, Spatial_Ratio[Temp-1]);

		if((!isnan(Chi) && Chi < Best[14]) || isnan(Best[14]))
		{
			Best[0] = JPsi_Parameters[Temp][0][1];
			Best[1] = JPsi_Parameters[Temp][0][2];
			Best[2] = JPsi_Parameters[Temp][1][1];
			Best[3] = JPsi_Parameters[Temp][1][2];
			Best[4] = JPsi_Parameters[Temp][2][1];
			Best[5] = JPsi_Parameters[Temp][2][2];
			Best[6] = JPsi_Parameters[Temp][3][1];
			Best[7] = JPsi_Parameters[Temp][3][2];
			Best[8] = JPsi_Parameters[Temp][4][1];
			Best[9] = JPsi_Parameters[Temp][4][2];
			Best[10] = Non_Parameters[Temp][0][1];
			Best[11] = Non_Parameters[Temp][0][2];
			Best[12] = Non_Parameters[Temp][1][1];
			Best[13] = Non_Parameters[Temp][1][2];
			Best[14] = Chi;
		}
	}
	JPsi_Parameters[Temp][0][1] = Best[0];
	JPsi_Parameters[Temp][0][2] = Best[1];
	JPsi_Parameters[Temp][1][1] = Best[2];
	JPsi_Parameters[Temp][1][2] = Best[3];
	JPsi_Parameters[Temp][2][1] = Best[4];
	JPsi_Parameters[Temp][2][2] = Best[5];
	JPsi_Parameters[Temp][3][1] = Best[6];
	JPsi_Parameters[Temp][3][2] = Best[7];
	JPsi_Parameters[Temp][4][1] = Best[8];
	JPsi_Parameters[Temp][4][2] = Best[9];
	Non_Parameters[Temp][0][1] = Best[10];
	Non_Parameters[Temp][0][2] = Best[11];
	Non_Parameters[Temp][1][1] = Best[12];
	Non_Parameters[Temp][1][2] = Best[13];

	long double gradn_1[14], gradn[14], sn_1[14], sn[14];
	long double betan = 0, betan_1;
	start_time = time(NULL);

	Gradient(gradn_1, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
	for(int i = 0; i < 14; i++)
		sn_1[i] = gradn_1[i];
	Minimize(sn_1, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
	do
	{
		betan_1 = betan;
		Gradient(gradn, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
		betan = PolakRibiere(gradn_1,gradn);
		if(betan < 0)
			for(int i = 0; i < 14; i++)
				sn[i] = gradn[i]+betan*sn_1[i];
		else
			for(int i = 0; i < 14; i++)
				sn[i] = gradn[i];
		Minimize(sn, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
		for(int i = 0; i < 14; i++)
		{
			gradn_1[i] = gradn[i];
			sn_1[i] = sn[i];
		}
	}while((abs(betan)>1e-12 || abs(betan_1)>1e-12) && time(NULL)-start_time < 9000);

	OutputFile.close();

	return(0);
}

void Minimize(long double sn[14], Spectral_Inter JPsi[5], Spectral_Inter PsiPrime[5], Spectral_Non Non[5], long double Vacuum_Euclidean[4][2], long double Vacuum_Spatial[4][7], long double Spatial_Ratio[4][7], int T)
{
	long double JPsi_Local[5][3];
	long double PsiPrime_Local[5][3];
	long double Non_Local[2][3];
	long double Medium_Euclidean[2];
	long double Medium_Spatial[7];
	long double length[2][2];

	for(int i = 0; i < 5; i++)
		for(int j = 0; j < 3; j++)
		{
			JPsi_Local[i][j] = JPsi_Parameters[i][j];
			PsiPrime_Local[i][j] = PsiPrime_Parameters[i][j];
			if(i < 2)
				Non_Local[i][j] = Non_Parameters[i][j];
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
	length[1][0] = (sn[6]>0)?(Random_Range[6][1]-JPsi_Local[3][1])/sn[6]:(Random_Range[6][0]-JPsi_Local[3][1])/sn[6];
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
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
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
		JPsi_Local[0][1] = JPsi_Parameters[0][1]+length[0][0]*fz[i][0]*sn[0];
		JPsi_Local[0][2] = JPsi_Parameters[0][2]+length[0][0]*fz[i][0]*sn[1];
		JPsi_Local[1][1] = JPsi_Parameters[1][1]+length[0][0]*fz[i][0]*sn[2];
		JPsi_Local[1][2] = JPsi_Parameters[1][2]+length[0][0]*fz[i][0]*sn[3];
		JPsi_Local[2][1] = JPsi_Parameters[2][1]+length[0][0]*fz[i][0]*sn[4];
		JPsi_Local[2][2] = JPsi_Parameters[2][2]+length[0][0]*fz[i][0]*sn[5];
		JPsi_Local[3][1] = JPsi_Parameters[3][1]+length[0][0]*fz[i][0]*sn[6];
		JPsi_Local[3][2] = JPsi_Parameters[3][2]+length[0][0]*fz[i][0]*sn[7];
		JPsi_Local[4][1] = JPsi_Parameters[4][1]+length[0][0]*fz[i][0]*sn[8];
		JPsi_Local[4][2] = JPsi_Parameters[4][2]+length[0][0]*fz[i][0]*sn[9];
		Non_Local[0][1] = Non_Parameters[0][1]+length[0][0]*fz[i][0]*sn[10];
		Non_Local[0][2] = Non_Parameters[0][2]+length[0][0]*fz[i][0]*sn[11];
		Non_Local[1][1] = Non_Parameters[1][1]+length[0][0]*fz[i][0]*sn[12];
		Non_Local[1][2] = Non_Parameters[1][2]+length[0][0]*fz[i][0]*sn[13];
		Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Local, PsiPrime_Local, Non_Local, false);
		Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Local, PsiPrime_Local, Non_Local, false);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Local, PsiPrime_Local, Non_Local, false);
		fz[i][1] = Print(JPsi_Local, PsiPrime_Local, Non_Local, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
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
			JPsi_Parameters[0][1] += length[0][0]*x*sn[0];
			JPsi_Parameters[0][2] += length[0][0]*x*sn[1];
			JPsi_Parameters[1][1] += length[0][0]*x*sn[2];
			JPsi_Parameters[1][2] += length[0][0]*x*sn[3];
			JPsi_Parameters[2][1] += length[0][0]*x*sn[4];
			JPsi_Parameters[2][2] += length[0][0]*x*sn[5];
			JPsi_Parameters[3][1] += length[0][0]*x*sn[6];
			JPsi_Parameters[3][2] += length[0][0]*x*sn[7];
			JPsi_Parameters[4][1] += length[0][0]*x*sn[8];
			JPsi_Parameters[4][2] += length[0][0]*x*sn[9];
			Non_Parameters[0][1] += length[0][0]*x*sn[10];
			Non_Parameters[0][2] += length[0][0]*x*sn[11];
			Non_Parameters[1][1] += length[0][0]*x*sn[12];
			Non_Parameters[1][2] += length[0][0]*x*sn[13];
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
		JPsi_Local[0][1] = JPsi_Parameters[0][1]+length[0][0]*u*sn[0];
		JPsi_Local[0][2] = JPsi_Parameters[0][2]+length[0][0]*u*sn[1];
		JPsi_Local[1][1] = JPsi_Parameters[1][1]+length[0][0]*u*sn[2];
		JPsi_Local[1][2] = JPsi_Parameters[1][2]+length[0][0]*u*sn[3];
		JPsi_Local[2][1] = JPsi_Parameters[2][1]+length[0][0]*u*sn[4];
		JPsi_Local[2][2] = JPsi_Parameters[2][2]+length[0][0]*u*sn[5];
		JPsi_Local[3][1] = JPsi_Parameters[3][1]+length[0][0]*u*sn[6];
		JPsi_Local[3][2] = JPsi_Parameters[3][2]+length[0][0]*u*sn[7];
		JPsi_Local[4][1] = JPsi_Parameters[4][1]+length[0][0]*u*sn[8];
		JPsi_Local[4][2] = JPsi_Parameters[4][2]+length[0][0]*u*sn[9];
		Non_Local[0][1] = Non_Parameters[0][1]+length[0][0]*u*sn[10];
		Non_Local[0][2] = Non_Parameters[0][2]+length[0][0]*u*sn[11];
		Non_Local[1][1] = Non_Parameters[1][1]+length[0][0]*u*sn[12];
		Non_Local[1][2] = Non_Parameters[1][2]+length[0][0]*u*sn[13];
		Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Local, PsiPrime_Local, Non_Local, false);
		Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Local, PsiPrime_Local, Non_Local, false);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Local, PsiPrime_Local, Non_Local, false);
		fu = Print(JPsi_Local, PsiPrime_Local, Non_Local, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);

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

void Gradient(long double grad[14], Spectral_Inter JPsi[5], Spectral_Inter PsiPrime[5], Spectral_Non Non[5], long double Medium_Euclidean[4][2], long double Medium_Spatial[4][7], long double Vacuum_Euclidean[4][2], long double Vacuum_Spatial[7], long double Spatial_Ratio[4][7], int T)
{
	long double f0, f1, h = 1e-5;
	long double Reduce_Euclidean[2][2];
	long double Reduce_Spatial[2][7];
	int i, j;
	
	OutputFile << "Gradient" << endl;
	Medium_Euclidean[T][0] = JPsi[T].Euclidean(.5, 0)+PsiPrime[T].Euclidean(.5, 0)+Non[T].Euclidean(.5, 0);
	Medium_Euclidean[T][1] = JPsi[T].Euclidean(.5, 3)+PsiPrime[T].Euclidean(.5, 3)+Non[T].Euclidean(.5, 3);
	Reduce_Euclidean[0][0] = PsiPrime[T].Euclidean(.5, 0)+Non[T].Euclidean(.5, 0);
	Reduce_Euclidean[0][1] = PsiPrime[T].Euclidean(.5, 3)+Non[T].Euclidean(.5, 3);
	Reduce_Euclidean[1][0] = JPsi[T].Euclidean(.5, 0)+PsiPrime[T].Euclidean(.5, 0);
	Reduce_Euclidean[1][1] = JPsi[T].Euclidean(.5, 3)+PsiPrime[T].Euclidean(.5, 3);
	for(int j = 0; j < 7; j++)
	{
		Medium_Spatial[T][j] = JPsi[T].Spatial((long double)(j)+.25)+PsiPrime[T].Spatial((long double)(j)+.25)+Non[T].Spatial((long double)(j)+.25);
		Reduce_Spatial[0][j] = PsiPrime[T].Spatial((long double)(j)+.25)+Non[T].Spatial((long double)(j)+.25);
		Reduce_Spatial[1][j] = JPsi[T].Spatial((long double)(j)+.25)+PsiPrime[T].Spatial((long double)(j)+.25);
	}
	f0 = Chi_Square(JPsi, PsiPrime, Non, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	//cerr << "Gradient " << f0 << endl;

	for(i = 0; i < 5; i++)
		for(j = 1; j < 3; j++)
		{
			JPsi[T].Add(h,i,j);
			Medium_Euclidean[T][0] = JPsi[T].Euclidean(.5, 0)+Reduce_Euclidean[0][0];
			Medium_Euclidean[T][1] = JPsi[T].Euclidean(.5, 3)+Reduce_Euclidean[0][1];
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T][j] = JPsi[T].Spatial((long double)(j)+.25)+Reduce_Spatial[0][j];
			f1 = Print(JPsi, PsiPrime, Non, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio, T);
			JPsi[T].Add(-h,i,j);
			grad[i*2+j-1] = (f0-f1)/h;
			//cerr << "Gradient " << f1 << " " << grad[0] << endl;
		}

	for(i = 0; i < 2; i++)
		for(j = 1; j < 3; j++)
		{
			Non[T].Add(h,i,j);
			Medium_Euclidean[T][0] = Reduce_Euclidean[1][0]+Non[T].Euclidean(.5, 0);
			Medium_Euclidean[T][1] = Reduce_Euclidean[1][1]+Non[T].Euclidean(.5, 3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T][j] = Reduce_Spatial[1][j]+Non[T].Spatial((long double)(j)+.25);
			f1 = Print(JPsi, PsiPrime, Non, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio, T);
			Non[T].Add(-h,i,j);
			grad[i*2+j+9] = (f0-f1)/h;
			//cerr << "Gradient " << f1 << " " << grad[0] << endl;
		}
}

long double Print(Spectral_Inter JPsi[5], Spectral_Inter PsiPrime[5], Spectral_Non Non[5], long double Medium_Euclidean[4][2], long double Vacuum_Euclidean[4][2], long double Medium_Spatial[4][7], long double Vacuum_Spatial[7], long double Spatial_Ratio[4][7], int Temp)
{
	OutputFile << Temp << "," << flush
	JPsi[Temp].Print(OutputFile);
	OutputFile << ",";
	Non[Temp].Print(OutputFile);
	OutputFile << "," << Medium_Euclidean[Temp-1][1]/Vacuum_Euclidean[Temp-1][1]-Medium_Euclidean[Temp-1][0]/Vacuum_Euclidean[Temp-1][0] << "," << flush;
	for(int j = 0; j < 7; j++)
	{
		OutputFile << Medium_Spatial[Temp-1][j]/Vacuum_Spatial[j] << "," << flush;
	}
	long double Chi = Chi_Square(JPsi, PsiPrime, Non, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	OutputFile << Chi << endl;
	return(Chi);
}

long double Chi_Square(Spectral_Inter JPsi[5], Spectral_Inter PsiPrime[5], Spectral_Non Non[5], long double Medium_Euclidean[4][2], long double Vacuum_Euclidean[4][2], long double Medium_Spatial[4][7], long double Vacuum_Spatial[7], long double Spatial_Ratio[4][7])
{
	long double answer = 0;

	for(int i = 0; i < 4; i++)
		answer += pow(Medium_Euclidean[i][1]/Vacuum_Euclidean[i][1]-Medium_Euclidean[i][0]/Vacuum_Euclidean[i][0]-.2,2)/.2;

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 7; j++)
			answer += pow(Medium_Spatial[i][j]/Vacuum_Spatial[j]-Spatial_Ratio[i][j],2)/Spatial_Ratio[i][j];

	for(int i = 0; i < 15; i++)
		answer += Least_Squares(JPsi[1].Read(i), JPsi[2].Read(i), JPsi[3].Read(i), JPsi[4].Read(i));
	for(int i = 0; i < 6; i++)
		answer += Least_Squares(Non[1].Read(i), Non[2].Read(i), Non[3].Read(i), Non[4].Read(i));

	return(answer);
}

long double Least_Squares(long double y1, long double y2, long double y3, long double y4)
{
	return(0.030408*pow(y1,2)+ 0.064712*pow(y2,2)+ 0.066696*pow(y3,2)+ 0.023816*pow(y4,2)- 0.074128*y1*y2- 0.025024*y1*y3- 0.038848*y2*y3+ 0.038336*y1*y4- 0.016448*y2*y4- 0.06952*y3*y4);
}
