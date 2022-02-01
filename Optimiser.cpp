#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<ctime>
#include<cfloat>
#ifdef _OPENMP
#include<omp.h>
#endif
using namespace std;

long double Q(long double, long double, long double, long double);
long double SpectralJPsi(long double, long double, long double[5][3], bool);
long double SpectralPsiPrime(long double, long double, long double[5][3], bool);
long double Width(long double, long double, long double[5][3]);
long double SpectralNon(long double, long double, long double[2][3], bool);

long double Spatial(long double, long double[5][3], long double[5][3], long double[2][3], bool);
long double Spatial_sInt(long double, long double, long double[5][3], long double[5][3], long double[2][3], bool);
long double Spatial_P0Int(long double, long double, long double[5][3], long double[5][3], long double[2][3], bool);
long double SpatialGeneralKernel(long double, long double, long double);
long double SpatialLorentz(long double, long double);
long double SpatialCutoffKernel(long double, long double, long double);

long double Euclidean(long double, long double, long double, long double[5][3], long double[5][3], long double[2][3], bool);
long double EuclideanKernel(long double, long double, long double, long double);

void Gradient(long double[14], long double[5][3], long double[5][3], long double[2][3], long double[2], long double[7], long double[7], long double);
long double PolakRibiere(long double[14], long double[14]);
void Minimize(long double[14], long double[5][3], long double[5][3], long double[2][3], long double[2], long double[7], long double[7], long double);

void mergeSort(long double[], int, int);
long double Uniform(long double, long double);
void Characterize_PsiPrime(long double[5][3], long double, pair<long double, long double>&, bool);
void Characterize_JPsi(long double[5][3], long double, pair<long double, long double>&, bool);
long double Chi_Square(long double[2], long double[2], long double[7], long double[7], long double[7]);
long double Print(long double[5][3], long double[5][3], long double[2][3], long double[2], long double[2], long double[7], long double[7], long double[7]); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square

long double Boundary[] = {0.00865, 0.0267, 0.0491, 0.0985, .421, .802, 1.01, 4.85};
ofstream OutputFile;
/*{0.381905, 5.67612, 2.80054, 4.18969, 0.175255, 3.226, 10.0313, 3.98761, 2.93239, 3.58459, 1.66367, 3.25552, 2.06784, 4.45895, 0.165222, 0.997411, 0.93765, 0.952806, 0.989706, 0.962973, 0.919522, 0.0129831} T=194 Min*/
/*{0.3, 3, 2.5, 5.5, 0.176, 3.6, 10.58, 4.57122, 3, 5.25, 1.655, 3.09, 2.45, 5.5, 0.165909, 1.00145, 0.973704, 0.983587, 0.999133, 0.962439, 0.948636, 0.0069526} T=194 Grid*/
/*{0.33729, 3.92964, 2.62562, 4.05461, 0.11714, 5.27443, 14.3297, 3.53248, 3.28106, 5.43784, 1.59513, 5.71738, 1.84899, 5.45489, 0.137133, 1.01195, 0.9315, 0.813736, 0.737843, 0.734627, 0.654341, 0.0301706} T=258 Min*/
/*{0.3, 3, 2.85, 3.5, 0.04, 5, 8.97, 1, 3, 1, 1.59, 3.09, 2.7, 5.5, 0.100996, 0.921911, 0.936966, 0.869235, 0.810326, 0.716585, 0.598793, 0.0586666} T=258 MeV Grid*/
/*{0.30626, 2.2706, 2.50082, 3.56841, 0.136461, 3.56011, 7.93198, 4.11694, 3.05149, 4.4995, 1.66565, 4.09551, 2.8846, 2.4588, 0.0935383, 0.911683, 0.854309, 0.599931, 0.498286, 0.519157, 0.36636, 0.0969059} T=320 MeV Min
/*{0.35, 2.5, 2.95, 2.5, 0.12, 5.5, 9.07, 1, 3, 1, 1.51, 3.09, 2.4, 5.5, 0.0803167, 0.956631, 0.958191, 0.877978, 0.640787, 0.333339, 0.503843, 0.203221} T=320 MeV Grid*/
/*{0.205794, 3.34307, 3.28467, 3.602, 0.16547, 4.75333, 9.07015, 1.14995, 1.61949, 3.70208, 1.59186, 4.70663, 2.29624, 2.3616, -0.103146, 1.13145, 1.07916, 0.697669, 0.0222384, 0.340961, 0.172499, 1.12332} T=400 MeV Min*/
/*{0.5, 3.5, 3.15, 1.5, 0.06, 5, 9.52, 1, 1, 1, 1.36, 3.09, 1.98, 5.5, 0.0650876, 0.732207, 0.430126, 0.364224, 0.396391, 0.0391483, -0.2501, 1.49091} T=400 Grid*/
int main(int argc, char* argv[])
{
	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.3,3.},{3.09946,2.5,5.5},{.106597,.176,3.6},{10.58,10.58,4.57122},{3,3,5.25}},
						{{.4024,.33729,3.92964},{3.125,2.62562,4.05461},{.372,.11714,5.27443},{8.97,14.3297,3.53248},{3,3.28106,5.43784}},
						{{.403047,.30626,2.2706},{3.151,2.50082,3.56841},{.68,.136461,3.56011},{9.07,7.93198,4.11694},{3,3.05149,4.4995}},
						{{.266834,.205794,3.34307},{3.1855,3.28467,3.602},{.98,.16547,4.75333},{9.52,9.07015,1.14995},{3,1.61949,3.70208}}};
	long double PsiPrime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,4.2,4.57122},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.655,3.09},{2.45,2.45,5.5}},
					       {{1.59,1.59513,5.71738},{2.7,1.84899,5.45489}},
					       {{1.51,1.66565,4.09551},{2.4,2.8846,2.4588}},
					       {{1.36,1.59186,4.70663},{1.98,2.29624,2.3616}}};

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
	else if(argc == 16)
	{
                strcat(File,"API.");
                strcat(File,argv[1]);
                strcat(File,".csv");
                OutputFile.open(File,ios::app);

		JPsi_Parameters[Temp][0][1] = atof(argv[2]);
		JPsi_Parameters[Temp][0][2] = atof(argv[3]);
		JPsi_Parameters[Temp][1][1] = atof(argv[4]);
		JPsi_Parameters[Temp][1][2] = atof(argv[5]);
		JPsi_Parameters[Temp][2][1] = atof(argv[6]);
		JPsi_Parameters[Temp][2][2] = atof(argv[7]);
		JPsi_Parameters[Temp][3][1] = atof(argv[8]);
		JPsi_Parameters[Temp][3][2] = atof(argv[9]);
		JPsi_Parameters[Temp][4][1] = atof(argv[10]);
		JPsi_Parameters[Temp][4][2] = atof(argv[11]);
		Non_Parameters[Temp][0][1] = atof(argv[12]);
		Non_Parameters[Temp][0][2] = atof(argv[13]);
		Non_Parameters[Temp][1][1] = atof(argv[14]);
		Non_Parameters[Temp][1][2] = atof(argv[15]);

		Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
		Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], false);
		cout << setprecision(18) << Print(JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Medium_Euclidean, Vacuum_Euclidean[Temp-1], Medium_Spatial, Vacuum_Spatial, Spatial_Ratio[Temp-1]) << endl;
		return(0);
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

		long double start[6];
		long double finish[6];
		long double step[6];
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
		Best[0] = JPsi_Parameters[Temp][0][1] = Uniform(.1,.5);
		Best[1] = JPsi_Parameters[Temp][0][2] = Uniform(1.,6.);
		Best[2] = JPsi_Parameters[Temp][1][1] = Uniform(2.5,3.5);
		Best[3] = JPsi_Parameters[Temp][1][2] = Uniform(1.,6.);
		Best[4] = JPsi_Parameters[Temp][2][1] = Uniform(.02,.18);
		Best[5] = JPsi_Parameters[Temp][2][2] = Uniform(1.,6.);
		Best[6] = JPsi_Parameters[Temp][3][1] = Uniform(5.,15.);
		Best[7] = JPsi_Parameters[Temp][3][2] = Uniform(1.,6.);
		Best[8] = JPsi_Parameters[Temp][4][1] = Uniform(1.5,3.5);
		Best[9] = JPsi_Parameters[Temp][4][2] = Uniform(1.,6.);
		Best[10] = Non_Parameters[Temp][0][1] = Uniform(1.59,1.79);
		Best[11] = Non_Parameters[Temp][0][2] = Uniform(1.,6.);
		Best[12] = Non_Parameters[Temp][1][1] = Uniform(1.5,5.);
		Best[13] = Non_Parameters[Temp][1][2] = Uniform(1.,6.);
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
		JPsi_Parameters[Temp][0][1] = Uniform(.1,.5);
		JPsi_Parameters[Temp][1][1] = Uniform(2.5,3.5);
		JPsi_Parameters[Temp][2][1] = Uniform(.02,.18);
		JPsi_Parameters[Temp][3][1] = Uniform(5.,15.);
		JPsi_Parameters[Temp][4][1] = Uniform(1.5,3.5);
		JPsi_Parameters[Temp][0][2] = Uniform(1.,6.);
		JPsi_Parameters[Temp][1][2] = Uniform(1.,6.);
		JPsi_Parameters[Temp][2][2] = Uniform(1.,6.);
		JPsi_Parameters[Temp][3][2] = Uniform(1.,6.);
		JPsi_Parameters[Temp][4][2] = Uniform(1.,6.);
		Non_Parameters[Temp][0][1] = Uniform(1.59,1.79);
		Non_Parameters[Temp][1][1] = Uniform(1.5,5.);
		Non_Parameters[Temp][0][2] = Uniform(1.,6.);
		Non_Parameters[Temp][1][2] = Uniform(1.,6.);
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
	for(int i = 0; i < 6; i++)
		sn_1[i] = gradn_1[i];
	Minimize(sn_1, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
	do
	{
		betan_1 = betan;
		Gradient(gradn, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
		betan = PolakRibiere(gradn_1,gradn);
		if(betan < 0)
			for(int i = 0; i < 6; i++)
				sn[i] = gradn[i]+betan*sn_1[i];
		else
			for(int i = 0; i < 6; i++)
				sn[i] = gradn[i];
		Minimize(sn, JPsi_Parameters[Temp], PsiPrime_Parameters[Temp], Non_Parameters[Temp], Vacuum_Euclidean[Temp-1], Vacuum_Spatial, Spatial_Ratio[Temp-1], T);
		for(int i = 0; i < 6; i++)
		{
			gradn_1[i] = gradn[i];
			sn_1[i] = sn[i];
		}
	}while((abs(betan)>1e-12 || abs(betan_1)>1e-12) && time(NULL)-start_time < 9000);

	OutputFile.close();

	return(0);
}

void Minimize(long double sn[14], long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], long double Vacuum_Euclidean[2], long double Vacuum_Spatial[7], long double Spatial_Ratio[7], long double T)
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

	length[0][0] = (sn[0]>0)?(.5-JPsi_Local[0][1])/sn[0]:(.1-JPsi_Local[0][1])/sn[0];
	length[0][1] = (sn[0]>0)?(.1-JPsi_Local[0][1])/sn[0]:(.5-JPsi_Local[0][1])/sn[0];
	length[1][0] = (sn[1]>0)?(6.-JPsi_Local[0][2])/sn[1]:(1.-JPsi_Local[0][2])/sn[1];
	length[1][1] = (sn[1]>0)?(1.-JPsi_Local[0][2])/sn[1]:(6.-JPsi_Local[0][2])/sn[1];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[2]>0)?(3.5-JPsi_Local[1][1])/sn[2]:(2.5-JPsi_Local[1][1])/sn[2];
	length[1][1] = (sn[2]>0)?(2.5-JPsi_Local[1][1])/sn[2]:(3.5-JPsi_Local[1][1])/sn[2];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[3]>0)?(6.-JPsi_Local[1][2])/sn[3]:(1.-JPsi_Local[1][2])/sn[3];
	length[1][1] = (sn[3]>0)?(1.-JPsi_Local[1][2])/sn[3]:(6.-JPsi_Local[1][2])/sn[3];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[4]>0)?(.18-JPsi_Local[2][1])/sn[4]:(.02-JPsi_Local[2][1])/sn[4];
	length[1][1] = (sn[4]>0)?(.02-JPsi_Local[2][1])/sn[4]:(.18-JPsi_Local[2][1])/sn[4];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[5]>0)?(6.-JPsi_Local[2][2])/sn[5]:(1.-JPsi_Local[2][2])/sn[5];
	length[1][1] = (sn[5]>0)?(1.-JPsi_Local[2][2])/sn[5]:(6.-JPsi_Local[2][2])/sn[5];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	/*length[1][0] = (sn[6]>0)?(15.-JPsi_Local[3][1])/sn[6]:(5.-JPsi_Local[3][1])/sn[6];
	length[1][1] = (sn[6]>0)?(5.-JPsi_Local[3][1])/sn[6]:(15.-JPsi_Local[3][1])/sn[6];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[7]>0)?(6.-JPsi_Local[3][2])/sn[7]:(1.-JPsi_Local[3][2])/sn[7];
	length[1][1] = (sn[7]>0)?(1.-JPsi_Local[3][2])/sn[7]:(6.-JPsi_Local[3][2])/sn[7];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[8]>0)?(3.5-JPsi_Local[4][1])/sn[8]:(1.5-JPsi_Local[4][1])/sn[8];
	length[1][1] = (sn[8]>0)?(1.5-JPsi_Local[4][1])/sn[8]:(3.5-JPsi_Local[4][1])/sn[8];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[9]>0)?(6.-JPsi_Local[4][2])/sn[9]:(1.-JPsi_Local[4][2])/sn[9];
	length[1][1] = (sn[9]>0)?(1.-JPsi_Local[4][2])/sn[9]:(6.-JPsi_Local[4][2])/sn[9];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];*/
	length[1][0] = (sn[10]>0)?(1.79-Non_Local[0][1])/sn[10]:(1.59-Non_Local[0][1])/sn[10];
	length[1][1] = (sn[10]>0)?(1.59-Non_Local[0][1])/sn[10]:(1.79-Non_Local[0][1])/sn[10];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[11]>0)?(6.-Non_Local[0][2])/sn[11]:(1.-Non_Local[0][2])/sn[11];
	length[1][1] = (sn[11]>0)?(1.-Non_Local[0][2])/sn[11]:(6.-Non_Local[0][2])/sn[11];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[12]>0)?(5.-Non_Local[1][1])/sn[12]:(1.5-Non_Local[1][1])/sn[12];
	length[1][1] = (sn[12]>0)?(1.5-Non_Local[1][1])/sn[12]:(5.-Non_Local[1][1])/sn[12];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[13]>0)?(6.-Non_Local[1][2])/sn[13]:(1.-Non_Local[1][2])/sn[13];
	length[1][1] = (sn[13]>0)?(1.-Non_Local[1][2])/sn[13]:(6.-Non_Local[1][2])/sn[13];
	length[0][0] = (length[0][0]<length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]>length[1][1])?length[0][1]:length[1][1];

	long double fz[101][2];
	long double a = 0, c = 1;
	int Min_i;

	OutputFile << "Line search" << endl;
	//cout << "Line search: a " << a << " c " << c << endl;
	for(int i = 0; i <= 100; i++)
	{
		fz[i][0] = a+(c-a)*i/100.;
		JPsi_Local[0][1] = JPsi_Parameters[0][1]+length[0][0]*fz[i][0]*sn[0];
		JPsi_Local[0][2] = JPsi_Parameters[0][2]+length[0][0]*fz[i][0]*sn[1];
		JPsi_Local[1][1] = JPsi_Parameters[1][1]+length[0][0]*fz[i][0]*sn[2];
		JPsi_Local[1][2] = JPsi_Parameters[1][2]+length[0][0]*fz[i][0]*sn[3];
		JPsi_Local[2][1] = JPsi_Parameters[2][1]+length[0][0]*fz[i][0]*sn[4];
		JPsi_Local[2][2] = JPsi_Parameters[2][2]+length[0][0]*fz[i][0]*sn[5];
		/*JPsi_Local[3][1] = JPsi_Parameters[3][1]+length[0][0]*fz[i][0]*sn[6];
		JPsi_Local[3][2] = JPsi_Parameters[3][2]+length[0][0]*fz[i][0]*sn[7];
		JPsi_Local[4][1] = JPsi_Parameters[4][1]+length[0][0]*fz[i][0]*sn[8];
		JPsi_Local[4][2] = JPsi_Parameters[4][2]+length[0][0]*fz[i][0]*sn[9];*/
		Non_Local[0][1] = Non_Parameters[0][1]+length[0][0]*fz[i][0]*sn[10];
		Non_Local[0][2] = Non_Parameters[0][2]+length[0][0]*fz[i][0]*sn[11];
		Non_Local[1][1] = Non_Parameters[1][1]+length[0][0]*fz[i][0]*sn[12];
		Non_Local[1][2] = Non_Parameters[1][2]+length[0][0]*fz[i][0]*sn[13];
		Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Local, PsiPrime_Local, Non_Local, false);
		Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Local, PsiPrime_Local, Non_Local, false);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Local, PsiPrime_Local, Non_Local, false);
		fz[i][1] = Print(JPsi_Local, PsiPrime_Local, Non_Local, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
		//cout << "Line Search " << fz[i][0] << " " << fz[i][1] << endl;
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
			/*JPsi_Parameters[3][1] += length[0][0]*x*sn[6];
			JPsi_Parameters[3][2] += length[0][0]*x*sn[7];
			JPsi_Parameters[4][1] += length[0][0]*x*sn[8];
			JPsi_Parameters[4][2] += length[0][0]*x*sn[9];*/
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
		/*JPsi_Local[3][1] = JPsi_Parameters[3][1]+length[0][0]*u*sn[6];
		JPsi_Local[3][2] = JPsi_Parameters[3][2]+length[0][0]*u*sn[7];
		JPsi_Local[4][1] = JPsi_Parameters[4][1]+length[0][0]*u*sn[8];
		JPsi_Local[4][2] = JPsi_Parameters[4][2]+length[0][0]*u*sn[9];*/
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

	for(int i = 0; i < 6; i++)
	{
		num += gradn[i]*(gradn[i]-gradn_1[i]);
		den += pow(gradn_1[i],2);
	}

	return(num/den);
}

void Gradient(long double grad[14], long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], long double Vacuum_Euclidean[2], long double Vacuum_Spatial[7], long double Spatial_Ratio[7], long double T)
{
	long double f0, f1, h = 1e-5;
	long double Medium_Euclidean[2];
	long double Medium_Spatial[7];

	OutputFile << "Gradient" << endl;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f0 = Chi_Square(Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	//cout << "Gradient " << f0 << endl;

	JPsi_Parameters[0][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[0][1] -= h;
	grad[0] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[0] << endl;

	JPsi_Parameters[0][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[0][2] -= h;
	grad[1] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[1] << endl;

	JPsi_Parameters[1][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[1][1] -= h;
	grad[2] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[2] << endl;

	JPsi_Parameters[1][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[1][2] -= h;
	grad[3] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[3] << endl;

	JPsi_Parameters[2][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[2][1] -= h;
	grad[4] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[4] << endl;

	JPsi_Parameters[2][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[2][2] -= h;
	grad[5] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[5] << endl;

	/*JPsi_Parameters[3][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[3][1] -= h;
	grad[6] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[6] << endl;

	JPsi_Parameters[3][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[3][2] -= h;
	grad[7] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[7] << endl;

	JPsi_Parameters[4][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[4][1] -= h;
	grad[8] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[8] << endl;

	JPsi_Parameters[4][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	JPsi_Parameters[4][2] -= h;
	grad[9] = (f0-f1)/h;
	//cout << "Gradient " << f1 << " " << grad[9] << endl;*/

	Non_Parameters[0][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	Non_Parameters[0][1] -= h;
	grad[10] = (f0-f1)/h;
	//cerr << "Gradient " << f1 << " " << grad[10] << endl;

	Non_Parameters[0][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	Non_Parameters[0][2] -= h;
	grad[11] = (f0-f1)/h;
	//cerr << "Gradient " << f1 << " " << grad[11] << endl;

	Non_Parameters[1][1] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	Non_Parameters[1][1] -= h;
	grad[12] = (f0-f1)/h;
	//cerr << "Gradient " << f1 << " " << grad[12] << endl;

	Non_Parameters[1][2] += h;
	Medium_Euclidean[0] = Euclidean(1./(2.*T), T, 0, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	Medium_Euclidean[1] = Euclidean(1./(2.*T), T, 3, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	for(int j = 0; j < 7; j++)
		Medium_Spatial[j] = Spatial((long double)(j)+.25, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, false);
	f1 = Print(JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	Non_Parameters[1][2] -= h;
	grad[13] = (f0-f1)/h;
	//cerr << "Gradient " << f1 << " " << grad[13] << endl;
}

long double Print(long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], long double Medium_Euclidean[2], long double Vacuum_Euclidean[2], long double Medium_Spatial[7], long double Vacuum_Spatial[7], long double Spatial_Ratio[7])
{
	for(int j = 0; j < 5; j++)
		OutputFile << JPsi_Parameters[j][1] << "," << JPsi_Parameters[j][2] << "," << flush;
	for(int j = 0; j < 2; j++)
		OutputFile << Non_Parameters[j][1] << "," << Non_Parameters[j][2] << "," << flush;
	OutputFile << Medium_Euclidean[1]/Vacuum_Euclidean[1]-Medium_Euclidean[0]/Vacuum_Euclidean[0] << "," << flush;
	for(int j = 0; j < 6; j++)
	{
		OutputFile << Medium_Spatial[j]/Vacuum_Spatial[j] << "," << flush;
	}
	long double Chi = Chi_Square(Medium_Euclidean, Vacuum_Euclidean, Medium_Spatial, Vacuum_Spatial, Spatial_Ratio);
	OutputFile << Chi << endl;
	return(Chi);
}

long double Chi_Square(long double Medium_Euclidean[2], long double Vacuum_Euclidean[2], long double Medium_Spatial[7], long double Vacuum_Spatial[7], long double Spatial_Ratio[7])
{
	long double answer;
	answer = pow(Medium_Euclidean[1]/Vacuum_Euclidean[1]-Medium_Euclidean[0]/Vacuum_Euclidean[0]-.2,2)/.2;
	for(int i = 0; i < 6; i++)
		answer += pow(Medium_Spatial[i]/Vacuum_Spatial[i]-Spatial_Ratio[i],2)/Spatial_Ratio[i];
	return(answer);
}

long double Spatial(long double z, long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], bool Vacuum)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	/*long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double Disp[] = {0.27963041316178305, 0.538469310105683, 0.7541667265708494, 0.9061798459386639, 0.9840853600948425};	//Displacement from center for unknown order Gauss-Kronrod integration from Mathematica
	long double w[] = {0.2829874178574912, 0.2728498019125589, 0.24104033922864776, 0.1868007965564926, 0.11523331662247445, 0.0425820367510818};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = 187.*M_PI/(2.*z);	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	long double Intervals = 2.*M_PI/z;
	int i, j, l;		//Counting varibles

	a = 0;
	b = M_PI/(2.*z);

	do
	{
		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 5; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;
			F += w[l+1]*Spatial_sInt(z, x1, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Vacuum);
			F += w[l+1]*Spatial_sInt(z, x2, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Vacuum);
		}
		F += w[0]*Spatial_sInt(z, a/2.+b/2., JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Vacuum);

		Answer += F*(b-a)/2.;
		a = b;
		b += Intervals;
	}while(a < Max);
	Answer += Spatial_P0Int(z, Max, JPsi_Parameters, PsiPrime_Parameters, Non_Parameters, Vacuum);

	return(Answer);
}

long double Spatial_sInt(long double z, long double P, long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], bool Vacuum)
{
	long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	/*long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};
	long double a, b;	//Sub-interval limits of integration
	long double Max = 552.25;	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int Intervals;
	int i, j, l;		//Counting varibles
	long double Stops[60] = {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 3, 6, 9, 12, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 552.25};

	Characterize_JPsi(JPsi_Parameters, P, zero, Vacuum);
	for(i = 0; i < 17; i++)
		Stops[i+25] = pow(zero.first+zero.second*Range[i],2);

	if(Q(P, PsiPrime_Parameters[0][0], PsiPrime_Parameters[0][1], PsiPrime_Parameters[0][2]) != 0)
	{
		Characterize_PsiPrime(PsiPrime_Parameters, P, zero, Vacuum);
		for(i = 0; i < 17; i++)
			Stops[i+42] = pow(zero.first+zero.second*Range[i],2);
		Stops[59] = pow(2.*Q(P, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 60;
	}
	else
	{
		Stops[42] = pow(2.*Q(P, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 43;
	}

	mergeSort(Stops, 0, Intervals-1);

	i = 0;
	while(Stops[i] < 0)
		i++;

	a = b = 0;
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 24; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*SpatialGeneralKernel(x1, P, z)*(SpectralJPsi(x1, P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x1, P, PsiPrime_Parameters, Vacuum)+SpectralNon(x1, P, Non_Parameters, Vacuum));
			F += w[l+1]*SpatialGeneralKernel(x2, P, z)*(SpectralJPsi(x2, P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x2, P, PsiPrime_Parameters, Vacuum)+SpectralNon(x2, P, Non_Parameters, Vacuum));
		}
		F += w[0]*SpatialGeneralKernel(a/2.+b/2., P, z)*(SpectralJPsi(a/2.+b/2., P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(a/2.+b/2., P, PsiPrime_Parameters, Vacuum)+SpectralNon(a/2.+b/2., P, Non_Parameters, Vacuum));

		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Spatial_P0Int(long double z, long double P0, long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], bool Vacuum)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};
	long double a, b;	//Sub-interval limits of integration
	long double Max = 552.25;	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int Intervals;
	int i, j, l;		//Counting varibles
	long double Stops[52] = {3, 6, 9, 12, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 552.25};

	Characterize_JPsi(JPsi_Parameters, P0, zero, Vacuum);
	for(i = 0; i < 17; i++)
		Stops[i+17] = pow(zero.first+zero.second*Range[i],2);

	if(Q(P0, PsiPrime_Parameters[0][0], PsiPrime_Parameters[0][1], PsiPrime_Parameters[0][2]) != 0)
	{
		Characterize_PsiPrime(PsiPrime_Parameters, P0, zero, Vacuum);
		for(i = 0; i < 17; i++)
			Stops[i+34] = pow(zero.first+zero.second*Range[i],2);
		Stops[51] = pow(2.*Q(P0, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 52;
	}
	else
	{
		Stops[34] = pow(2.*Q(P0, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 35;
	}

	mergeSort(Stops, 0, Intervals-1);

	i = 0;
	while(Stops[i] < 0)
		i++;

	a = b = 0;
	i = 0;
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 9; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*SpatialCutoffKernel(x1, P0, z)*(SpectralJPsi(x1, P0, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x1, P0, PsiPrime_Parameters, Vacuum)+SpectralNon(x1, P0, Non_Parameters, Vacuum));
			F += w[l+1]*SpatialCutoffKernel(x2, P0, z)*(SpectralJPsi(x2, P0, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x2, P0, PsiPrime_Parameters, Vacuum)+SpectralNon(x2, P0, Non_Parameters, Vacuum));
		}
		F += w[0]*SpatialCutoffKernel(a/2.+b/2., P0, z)*(SpectralJPsi(a/2.+b/2., P0, JPsi_Parameters, Vacuum)+SpectralPsiPrime(a/2.+b/2., P0, PsiPrime_Parameters, Vacuum)+SpectralNon(a/2.+b/2., P0, Non_Parameters, Vacuum));
		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Euclidean(long double tau, long double T, long double P, long double JPsi_Parameters[5][3], long double PsiPrime_Parameters[5][3], long double Non_Parameters[2][3], bool Vacuum)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};
	long double a, b;	//Sub-interval limits of integration
	long double Max = pow(400.,2);	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int Intervals;
	int i, j, l;		//Counting varibles
	long double Stops[78] = {3, 6, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 604, 704, 804, 904, 1004, 2004, 3004, 4004, 5004, 6004, 7004, 8004, 9004, 10004, 10004, 20004, 30004, 40004, 50004, 60004, 70004, 80004, 90004, 100004, 110004, 120004, 130004, 140004, 150004, 160000};

	Characterize_JPsi(JPsi_Parameters, P, zero, Vacuum);
	for(i = 0; i < 17; i++)
		Stops[i+43] = pow(zero.first+zero.second*Range[i],2);

	if(Q(P, PsiPrime_Parameters[0][0], PsiPrime_Parameters[0][1], PsiPrime_Parameters[0][2]) != 0)
	{
		Characterize_PsiPrime(PsiPrime_Parameters, P, zero, Vacuum);
		for(i = 0; i < 17; i++)
			Stops[i+60] = pow(zero.first+zero.second*Range[i],2);
		Stops[77] = pow(2.*Q(P, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 78;
	}
	else
	{
		Stops[60] = pow(2.*Q(P, Non_Parameters[0][0], Non_Parameters[0][1], Non_Parameters[0][2]),2);
		Intervals = 61;
	}

	mergeSort(Stops, 0, Intervals-1);

	i = 0;
	while(Stops[i] < 0)
		i++;

	a = b = 0;
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 9; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*EuclideanKernel(x1, P, tau, T)*(SpectralJPsi(x1, P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x1, P, PsiPrime_Parameters, Vacuum)+SpectralNon(x1, P, Non_Parameters, Vacuum));
			F += w[l+1]*EuclideanKernel(x2, P, tau, T)*(SpectralJPsi(x2, P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(x2, P, PsiPrime_Parameters, Vacuum)+SpectralNon(x2, P, Non_Parameters, Vacuum));
		}
		F += w[0]*EuclideanKernel(a/2.+b/2., P, tau, T)*(SpectralJPsi(a/2.+b/2., P, JPsi_Parameters, Vacuum)+SpectralPsiPrime(a/2.+b/2., P, PsiPrime_Parameters, Vacuum)+SpectralNon(a/2.+b/2., P, Non_Parameters, Vacuum));

		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

void Characterize_PsiPrime(long double Parameters[5][3], long double P, pair<long double, long double>& zero, bool Vacuum)
{
	long double M_old = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	long double M_new;
	long double h = .001;
	int i = 0;

	M_new = M_old - h*(SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum)-SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum))/(SpectralJPsi(pow(M_old-h,2), P, Parameters, Vacuum)-2.*SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum)+SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum));

	while(abs(M_new/M_old-1.) > 1e-5 && i < 20)
	{
		M_old = M_new;
		M_new = M_old - h*(SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum)-SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum))/(SpectralJPsi(pow(M_old-h,2), P, Parameters, Vacuum)-2.*SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum)+SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum));
		i++;
	}

	zero.first = M_new;
	zero.second = 2.*sqrt(-pow(h,2)*SpectralPsiPrime(pow(M_new,2), P, Parameters, Vacuum)/(2.*(SpectralPsiPrime(pow(M_new-h,2), P, Parameters, Vacuum)-2.*SpectralPsiPrime(pow(M_new,2), P, Parameters, Vacuum)+SpectralPsiPrime(pow(M_new+h,2), P, Parameters, Vacuum))));

	if(isnan(zero.second) || M_new < 2.5 || M_new > 4.5)
	{
		zero.first = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		zero.second = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
	}
}

void Characterize_JPsi(long double Parameters[5][3], long double P, pair<long double, long double>& zero, bool Vacuum)
{
	long double M_old = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	long double M_new;
	long double h = .001;
	int i = 0;

	M_new = M_old - h*(SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum)-SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum))/(SpectralJPsi(pow(M_old-h,2), P, Parameters, Vacuum)-2.*SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum)+SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum));

	while(abs(M_new/M_old-1.) > 1e-5 && i < 20)
	{
		M_old = M_new;
		M_new = M_old - h*(SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum)-SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum))/(SpectralJPsi(pow(M_old-h,2), P, Parameters, Vacuum)-2.*SpectralJPsi(pow(M_old,2), P, Parameters, Vacuum)+SpectralJPsi(pow(M_old+h,2), P, Parameters, Vacuum));
		i++;
	}

	zero.first = M_new;
	zero.second = 2.*sqrt(-pow(h,2)*SpectralJPsi(pow(M_new,2), P, Parameters, Vacuum)/(2.*(SpectralJPsi(pow(M_new-h,2), P, Parameters, Vacuum)-2.*SpectralJPsi(pow(M_new,2), P, Parameters, Vacuum)+SpectralJPsi(pow(M_new+h,2), P, Parameters, Vacuum))));

	if(isnan(zero.second) || M_new < 1.5 || M_new > 4.5)
	{
		zero.first = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		zero.second = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
	}
}

long double SpatialGeneralKernel(long double s, long double P, long double z)
{
	return(cos(P*z)/(s+pow(P,2)));
}

long double SpatialLorentz(long double s, long double z)
{
	return(M_PI*exp(-sqrt(s)*z)/(2.*sqrt(s)));
}

long double SpatialCutoffKernel(long double s, long double P0, long double z)
{
	return(-((z*(-2.*P0*z*(12852.-5868.*pow(sqrt(s)*z,2)+1737.*pow(sqrt(s)*z,4)-78.*pow(sqrt(s)*z,6)+pow(P0*z,8)+pow(sqrt(s)*z,8)+3.*pow(P0*z,4)*(579.+26.*pow(sqrt(s)*z,2)+2.*pow(sqrt(s)*z,4))+pow(P0,6)*(78.*pow(z,6)+4.*s*pow(z,8))+pow(P0,2)*(5868.*pow(z,2)+2418.*s*pow(z,4)-78.*pow(s,2)*pow(z,6)+4.*pow(s,3)*pow(z,8)))*cos(P0*z)+(-5400.+33948.*pow(sqrt(s)*z,2)-16074.*pow(sqrt(s)*z,4)+2301.*pow(sqrt(s)*z,6)-88.*pow(sqrt(s)*z,8)+pow(P0*z,10)+pow(sqrt(s)*z,10)+pow(P0*z,8)*(84.+5.*pow(sqrt(s)*z,2))+pow(P0*z,6)*(2037.+164.*pow(sqrt(s)*z,2)+10.*pow(sqrt(s)*z,4))+pow(P0*z,4)*(10782.+3975.*pow(sqrt(s)*z,2)-12.*pow(sqrt(s)*z,4)+10.*pow(sqrt(s)*z,6))+pow(P0*z,2)*(27756.-4716.*pow(sqrt(s)*z,2)+4239.*pow(sqrt(s)*z,4)-180.*pow(sqrt(s)*z,6)+5.*pow(sqrt(s)*z,8)))*sin(P0*z)))/(pow(P0*z,12)+6.*pow(P0*z,10)*(15.+pow(sqrt(s)*z,2))+3.*pow(P0*z,8)*(819.+90.*pow(sqrt(s)*z,2)+5.*pow(sqrt(s)*z,4))+pow(-36.+216.*pow(sqrt(s)*z,2)-45.*pow(sqrt(s)*z,4)+pow(sqrt(s)*z,6),2)+4.*pow(P0*z,6)*(4878.+1593.*pow(sqrt(s)*z,2)+45.*pow(sqrt(s)*z,4)+5.*pow(sqrt(s)*z,6))+3.*pow(P0*z,4)*(16632.+6120.*pow(sqrt(s)*z,2)+2610.*pow(sqrt(s)*z,4)-60.*pow(sqrt(s)*z,6)+5.*pow(sqrt(s)*z,8))+6.*pow(P0*z,2)*(2592.+12312.*pow(sqrt(s)*z,2)-3060.*pow(sqrt(s)*z,4)+1062.*pow(sqrt(s)*z,6)-45.*pow(sqrt(s)*z,8)+pow(sqrt(s)*z,10)))));
}

long double EuclideanKernel(long double s, long double P, long double tau, long double T)
{
	return(cosh(sqrt(s+pow(P,2))*(tau-1./(2.*T)))/(2.*sqrt(s+pow(P,2))*sinh(sqrt(s+pow(P,2))/(2.*T))));
}

long double SpectralNon(long double s, long double P, long double Parameters[2][3], bool Vacuum)
{
	if(Vacuum)
		P = 0;

	static long double old_P = P;
	static long double M = Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	static long double n = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
		n = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	}

	if(s < pow(2.*M,2))
		return(0);
	return((3.*s)/(8.*pow(M_PI,2))*sqrt(1.-pow(pow(2.*M,2)/s,n)));
}

long double WidthPsiPrime(long double E, long double P, long double Parameters[5][3])
{
	static long double old_P = P;
	static long double M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	static long double a = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
	static long double b = Q(P,Parameters[4][0],Parameters[4][1],Parameters[4][2]);
	static long double E0 = M-2.*Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2]);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		a = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
		b = Q(P,Parameters[4][0],Parameters[4][1],Parameters[4][2]);
		E0 = M-2.*Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2]);
	}

	return(exp(a*(sqrt(pow(sqrt(pow(E,2)+pow(P,2))-sqrt(pow(E0,2)+pow(P,2)),2)+pow(b,2))-sqrt(pow(E,2)+pow(P,2))+sqrt(pow(E0,2)+pow(P,2)))/(2.*(sqrt(pow(sqrt(pow(M,2)+pow(P,2))-sqrt(pow(E0,2)+pow(P,2)),2)+pow(b,2))-sqrt(pow(M,2)+pow(P,2))+sqrt(pow(E0,2)+pow(P,2))))*(M-E0-sqrt(pow(M-E0,2)+pow(b,2)))+a/2.*(E0-M+sqrt(pow(b,2)+pow(M-E0,2)))));
}

long double WidthJPsi(long double E, long double P, long double Parameters[5][3])
{
	static long double old_P = P;
	static long double M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	static long double a = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
	static long double b = Q(P,Parameters[4][0],Parameters[4][1],Parameters[4][2]);
	static long double E0 = M-2.*Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2]);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		a = Q(P,Parameters[3][0],Parameters[3][1],Parameters[3][2]);
		b = Q(P,Parameters[4][0],Parameters[4][1],Parameters[4][2]);
		E0 = M-2.*Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2]);
	}

	return(exp(a*(sqrt(pow(sqrt(pow(E,2)+pow(P,2))-sqrt(pow(E0,2)+pow(P,2)),2)+pow(b,2))-sqrt(pow(E,2)+pow(P,2))+sqrt(pow(E0,2)+pow(P,2)))/(2.*(sqrt(pow(sqrt(pow(M,2)+pow(P,2))-sqrt(pow(E0,2)+pow(P,2)),2)+pow(b,2))-sqrt(pow(M,2)+pow(P,2))+sqrt(pow(E0,2)+pow(P,2))))*(M-E0-sqrt(pow(M-E0,2)+pow(b,2)))+a/2.*(E0-M+sqrt(pow(b,2)+pow(M-E0,2)))));
}

long double SpectralPsiPrime(long double s, long double P, long double Parameters[5][3], bool Vacuum)
{
	if(Vacuum)
		P = 0;

	static long double old_P = P;
	if(0 == Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]))
		return(0);
	static long double M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	static long double A = 2.*M*Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	long double Gamma = Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2])*WidthPsiPrime(sqrt(s),P,Parameters);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		A = 2.*M*Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	}

	return((A*Gamma*M*sqrt((s+pow(P,2))/(pow(M,2)+pow(P,2))))/(M_PI*(pow(s-pow(M,2),2)+pow(Gamma*M*sqrt((s+pow(P,2))/(pow(M,2)+pow(P,2))),2))));
}

long double SpectralJPsi(long double s, long double P, long double Parameters[5][3], bool Vacuum)
{
	if(Vacuum)
		P = 0;

	static long double old_P = P;
	if(0 == Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]))
		return(0);
	static long double M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	static long double A = 2.*M*Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	long double Gamma = Q(P,Parameters[2][0],Parameters[2][1],Parameters[2][2])*WidthJPsi(sqrt(s),P,Parameters);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
		A = 2.*M*Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	}

	return((A*Gamma*M*sqrt((s+pow(P,2))/(pow(M,2)+pow(P,2))))/(M_PI*(pow(s-pow(M,2),2)+pow(Gamma*M*sqrt((s+pow(P,2))/(pow(M,2)+pow(P,2))),2))));
}

long double Q(long double P, long double Q0, long double QV, long double P0)
{
	return((Q0*pow(P0,2)+QV*pow(P,2))/(pow(P0,2)+pow(P,2)));
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
