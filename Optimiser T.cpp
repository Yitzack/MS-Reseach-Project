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
long double Chi_Square(long double[2], long double[7], int);
long double Print(Spectral_Inter*, Spectral_Inter*, Spectral_Non*, long double[2], long double[7], int); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square

long double Random_Range[14][2] = {{.1,.5},{1.,6.},{1.,3.5},{1.,6.},{.02,.18},{1.,6.},{5.,15.},{1.,6.},{1.5,3.5},{1.,6.},{1.59,1.79},{1.,6.},{1.5,5.},{1.,6.}};
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
/*{0.296663, 2.99912, 2.41053, 5.50001, 0.168958, 3.60003, 10.58, 4.57122, 3, 5.25, 1.69802, 4.00006, 3.00002, 5.49997, 0.181899, 1.00337, 1.00995, 0.995169, 0.985732, 0.980595, 0.973869, 0.956278, 0.00181892} T=194 MeV Min
{0.3, 3, 2.45, 5.5, 0.18, 3.6, 10.58, 4.57122, 3, 5.25, 1.7, 5.5, 2.25, 6, 0.169186, 1.00663, 1.03674, 1.00002, 0.997114, 0.999272, 0.946506, 0.939859, 0.00737184} T=194 MeV Grid
{0.445826, 3.86826, 2.58405, 4.01049, 0.535639, 4.20052, 10.3425, 2.64938, 2.81697, 2.12552, 1.69153, 3.42941, 2.26904, 4.4543, 0.159886, 1.00841, 1.06611, 0.933118, 0.803052, 0.778368, 0.74589, 0.664434, 0.0225266} T=258 Min
{0.3, 3, 2.81, 3.6, 0.032, 4.9, 8.97, 1, 3, 1, 1.59, 3.09, 2.7, 5.5, 0.112458, 1.00738, 0.909125, 0.9208, 0.845615, 0.803736, 0.752825, 0.621187, 0.0496041} T=258 Grid
{0.184891, 3.32152, 1.3123, 5.23093, 0.064861, 5.80902, 5.03399, 5.12164, 1.55583, 5.24328, 1.64566, 5.31292, 3.93386, 5.40511, 0.143769, 1.00805, 0.864955, 0.847608, 0.667474, 0.619708, 0.473687, 0.376542, 0.0332863} T=320 MeV Min
{0.33, 2.4, 2.87, 2.5, 0.084, 5.5, 7.93198, 4.11694, 3.05149, 4.4995, 1.7, 5.5, 3, 3.5, 0.0622976, 1.0064, 0.936273, 0.864206, 0.728055, 0.50826, 0.476024, 0.400429, 0.105599} T=320 MeV Grid
{0.529109, 4.49175, 2.67237, 2.20553, 0.165846, 5.87198, 6.08157, 2.67594, 2.80637, 3.47611, 1.56825, 3.52402, 1.78377, 3.60395, 0.101371, 1.00725, 1.02967, 0.875425, 0.402078, 0.00803392, 0.235885, 0.255533, 0.519549} T=400 MeV Min
{0.2, 2, 2.6, 2, 0.1, 3.5, 9.52, 1, 3, 1, 1.36, 3.09, 1.98, 5.5, 0.0643115, 1.00325, 1.02997, 1.1458, 0.92826, 0.243907, 0.122592, 0.591802, 1.729} T=400 MeV Grid*/
int main(int argc, char* argv[])
{
	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.29663,2.99912},{3.09946,2.41053,5.50001},{.106597,.168958,3.60003},{10.58,10.58,4.57122},{3,3,5.25}},
						{{.4024,.445826,3.86826},{3.125,2.58405,4.01049},{.372,.535639,4.20052},{8.97,10.3425,2.64938},{3,2.81697,2.12552}},
						{{.403047,.184891,3.32152},{3.151,1.3123,5.23093},{.68,.064861,5.80902},{9.07,5.03399,5.12164},{3,1.55583,5.24328}},
						{{.266834,.529109,4.49175},{3.1855,2.67237,2.20553},{.98,.165846,5.87198},{9.52,6.08157,2.67594},{3,2.80637,3.47611}}};
	long double Psi_Prime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,10.58,4.57122},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.69802,4.00006},{2.45,3.00002,5.49997}},
					       {{1.59,1.69153,3.42941},{2.7,2.26904,4.4543}},
					       {{1.51,1.64566,5.31292},{2.4,3.93386,5.40511}},
					       {{1.36,1.56825,3.52402},{1.98,1.78377,3.60395}}};

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

	time_t start_time = time(NULL);
	long double gradn_1[14], gradn[14], sn_1[14], sn[14];
	long double betan = 0, betan_1;

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
			grad[i] = (f0-f1)/h;
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
			grad[i] = (f0-f1)/h;
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
		answer += pow(Medium_Spatial[i]/Vacuum_Spatial[i]-Spatial_Ratio[T-1][i],2)/Spatial_Ratio[T-1][i];
	return(answer);
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
