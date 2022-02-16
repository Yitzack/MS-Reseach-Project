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

/*(*Mathematica code payload for finding the best results and turning it into a Frankenstein's monster of different results*)
Length[table=Select[Import["Tables/Optimiser_Output.csv"],Length[#]==98&&NumericQ[#[[98]]]&]]
Table[MinMax[table[[All,i]]],{i,{24,48,72,96,97,98}}]
pos=Flatten[Table[Position[table[[All,i]],Min[table[[All,i]]]][[1]],{i,{24,48,72,96,97,98}}]]
Frankenstein=Join[table[[pos[[1]],1;;24]],table[[pos[[2]],25;;48]],table[[pos[[3]],49;;72]],table[[pos[[4]],73;;96]]]
Total[Table[0.030408*Frankenstein[[i+1]]^2+0.064712*Frankenstein[[i+25]]^2+0.066696*Frankenstein[[i+49]]^2+0.023816*Frankenstein[[i+73]]^2-0.074128*Frankenstein[[i+1]]*Frankenstein[[i+25]]-0.025024*Frankenstein[[i+1]]*Frankenstein[[i+49]]-0.038848*Frankenstein[[i+25]]*Frankenstein[[i+49]]+0.038336*Frankenstein[[i+1]]*Frankenstein[[i+73]]-0.016448*Frankenstein[[i+25]]*Frankenstein[[i+73]]-0.06952*Frankenstein[[i+49]]*Frankenstein[[i+73]],{i,14}]]
Total[Frankenstein[[{24,48,72,96}]]]+%
Join[Frankenstein,{%%,%}]
*/

void Gradient(long double[20], long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], long double[4][2], long double[4][7]);
long double PolakRibiere(long double[20], long double[20]);
void Minimize(long double[20], long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], long double[4][2], long double[4][7]);

long double Chi_Square(Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], long double[4][2], long double[4][7], int);
long double Least_Squares(long double, long double, long double, long double);
long double Print(long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], long double[4][2], long double[4][7]); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square
long double Uniform(long double, long double);
long double Protected_Uniform(long double, long double, long double, long double);
void DataLoad(Spectral_Inter*[5], Spectral_Non*[5], long double[20]);
				   //A, PA, DeltaM1, DeltaM4, PM, Gamma, PGamma, a, Pa, Delta, PDelta, DeltaMQ1, DeltaPMQ4, PMQ, n, Pn
long double Random_Range[16][2] = {{.1,.5},{1.,6.},{-0.59946, 0.40054},{-0.6855, 0.3145},{1.,6.},{.02,.18},{1.,6.},{5.,15.},{1.,6.},{1.5,3.5},{1.,6.},{-.065,.135},{-.05,.43},{1.,6.},{1.5,5.},{1.,6.}};
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

int main(int argc, char* argv[])
{
	//long double Deviation_Points[20] = {A1, A4, PA1, PA4, DeltaM1, DeltaM4, PM1, PM4, Gamma1, Gamma4, PGamma1, PGamma4, DeltaMQ1, DeltaMQ4, PMQ1, PMQ4, n1, n4, Pn1, Pn4};{0.115296, 0.286618, 2.08208, 1.56006, -0.405478, -0.285573, 5.26038, 5.25057, 0.0238598, 0.0716614, 2.50265, 1.90889, -0.0032708, 0.00857658, 4.48369, 3.60523, 3.06481, 3.43242, 4.25095, 5.06179, 1, 0.115296, 2.08208, 2.69398, 5.26038, 0.0238598, 2.50265, 10.58, 4.00927, 3, 1, 1.65173, 4.48369, 3.06481, 4.25095, -0.0925339, 1.01848, 1.25657, 1.4179, 1.51063, 2.34582, 1.23538, -3.89896, 27.5353, 2, 0.168522, 1.9199, 2.75677, 5.25733, 0.0387107, 2.31818, 8.97, 2.61299, 3, 2.22449, 1.59041, 4.21077, 3.17902, 4.50286, 0.0314608, 0.984802, 0.797205, 0.7265, 0.729743, 0.751454, 0.761201, -0.88084, 3.86127, 3, 0.220085, 1.76279, 2.81886, 5.25438, 0.0530976, 2.13948, 9.07, 5.68392, 3, 3.59181, 1.51398, 3.94638, 3.28966, 4.7469, 0.051174, 0.979265, 0.73075, 0.643505, 0.655122, 0.698825, 0.776109, 3.83516, 32.2206, 4, 0.286618, 1.56006, 2.89993, 5.25057, 0.0716614, 1.90889, 9.52, 1.00741, 3, 1.00741, 1.36858, 3.60523, 3.43242, 5.06179, 0.0632431, 0.976795, 0.70255, 0.612515, 0.645419, 0.990207, -4.08444, -0.134695, 77.9266, 141.544}
	long double Deviation_Points[20] = {0.115296, 0.286618, 2.08208, 1.56006, -0.405478, -0.285573, 5.26038, 5.25057, 0.0238598, 0.0716614, 2.50265, 1.90889, -0.0032708, 0.00857658, 4.48369, 3.60523, 3.06481, 3.43242, 4.25095, 5.06179};

	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.317797,2.35534},{3.09946,3.09946,4.55517},{.106597,.106597,2.34302},{10.58,10.58,4.00927},{3,3,1}},
						{{.4024,.4024,4.5535},{3.125,3.125,1.91861},{.372,.372,3.07496},{8.97,8.97,2.61299},{3,3,2.22449}},
						{{.403047,.403047,2.44122},{3.151,3.151,5.30318},{.68,.68,1.92309},{9.07,9.07,5.68392},{3,3,3.59181}},
						{{.266834,.266834,2.10741},{3.1855,3.1855,2.10741},{.98,.98,3.90741},{9.52,9.52,1.00741},{3,3,1.00741}}};
	long double PsiPrime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,10.58,4.00927},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.655,2.4807},{2.45,2.45,1.93449}},
					       {{1.59,1.59,2.41953},{2.7,2.7,1.92468}},
					       {{1.51,1.51,2.63744},{2.4,2.4,3.18298}},
					       {{1.36,1.36,3.09741},{1.98,1.98,5.50741}}};

	JPsi_Parameters[1][0][1] = Deviation_Points[0];
	JPsi_Parameters[2][0][1] = (71.*Deviation_Points[0]+32.*Deviation_Points[1])/103.;
	JPsi_Parameters[3][0][1] = (40.*Deviation_Points[0]+63.*Deviation_Points[1])/103.;
	JPsi_Parameters[4][0][1] = Deviation_Points[1];
	JPsi_Parameters[1][0][2] = Deviation_Points[2];
	JPsi_Parameters[2][0][2] = (71.*Deviation_Points[2]+32.*Deviation_Points[3])/103.;
	JPsi_Parameters[3][0][2] = (40.*Deviation_Points[2]+63.*Deviation_Points[3])/103.;
	JPsi_Parameters[4][0][2] = Deviation_Points[3];
	JPsi_Parameters[1][1][1] = JPsi_Parameters[1][1][0]+Deviation_Points[4];
	JPsi_Parameters[2][1][1] = JPsi_Parameters[2][1][0]+(71.*Deviation_Points[4]+32.*Deviation_Points[5])/103.;
	JPsi_Parameters[3][1][1] = JPsi_Parameters[3][1][0]+(40.*Deviation_Points[4]+63.*Deviation_Points[5])/103.;
	JPsi_Parameters[4][1][1] = JPsi_Parameters[4][1][0]+Deviation_Points[5];
	JPsi_Parameters[1][1][2] = Deviation_Points[6];
	JPsi_Parameters[2][1][2] = (71.*Deviation_Points[6]+32.*Deviation_Points[7])/103.;
	JPsi_Parameters[3][1][2] = (40.*Deviation_Points[6]+63.*Deviation_Points[7])/103.;
	JPsi_Parameters[4][1][2] = Deviation_Points[7];
	JPsi_Parameters[1][2][1] = Deviation_Points[8];
	JPsi_Parameters[2][2][1] = (71.*Deviation_Points[8]+32.*Deviation_Points[9])/103.;
	JPsi_Parameters[3][2][1] = (40.*Deviation_Points[8]+63.*Deviation_Points[9])/103.;
	JPsi_Parameters[4][2][1] = Deviation_Points[9];
	JPsi_Parameters[1][2][2] = Deviation_Points[10];
	JPsi_Parameters[2][2][2] = (71.*Deviation_Points[10]+32.*Deviation_Points[11])/103.;
	JPsi_Parameters[3][2][2] = (40.*Deviation_Points[10]+63.*Deviation_Points[11])/103.;
	JPsi_Parameters[4][2][2] = Deviation_Points[11];
	Non_Parameters[1][0][1] = Non_Parameters[1][0][0]+Deviation_Points[12];
	Non_Parameters[2][0][1] = Non_Parameters[2][0][0]+(71.*Deviation_Points[12]+32.*Deviation_Points[13])/103.;
	Non_Parameters[3][0][1] = Non_Parameters[3][0][0]+(40.*Deviation_Points[12]+63.*Deviation_Points[13])/103.;
	Non_Parameters[4][0][1] = Non_Parameters[4][0][0]+Deviation_Points[13];
	Non_Parameters[1][0][2] = Deviation_Points[14];
	Non_Parameters[2][0][2] = (71.*Deviation_Points[14]+32.*Deviation_Points[15])/103.;
	Non_Parameters[3][0][2] = (40.*Deviation_Points[14]+63.*Deviation_Points[15])/103.;
	Non_Parameters[4][0][2] = Deviation_Points[15];
	Non_Parameters[1][1][1] = Deviation_Points[16];
	Non_Parameters[2][1][1] = (71.*Deviation_Points[16]+32.*Deviation_Points[17])/103.;
	Non_Parameters[3][1][1] = (40.*Deviation_Points[16]+63.*Deviation_Points[17])/103.;
	Non_Parameters[4][1][1] = Deviation_Points[17];
	Non_Parameters[1][1][2] = Deviation_Points[18];
	Non_Parameters[2][1][2] = (71.*Deviation_Points[18]+32.*Deviation_Points[19])/103.;
	Non_Parameters[3][1][2] = (40.*Deviation_Points[18]+63.*Deviation_Points[19])/103.;
	Non_Parameters[4][1][2] = Deviation_Points[19];

	Spectral_Inter* JPsi[5];
	Spectral_Inter* Psi_Prime[5];
	Spectral_Non* Non[5];
	for(int i = 0; i < 5; i++)
	{
		JPsi[i] = new Spectral_Inter(JPsi_Parameters[i], i, bool(i));
		Psi_Prime[i] = new Spectral_Inter(PsiPrime_Parameters[i], i, bool(i));
		Non[i] = new Spectral_Non(Non_Parameters[i], i, bool(i));
	}

	long double Medium_Spatial[4][7];
	long double Medium_Euclidean[4][2];
	bool Cycle = true;
	for(int T = 1; T < 5; T++)
	{
		Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
		Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
	}
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
	if(argc == 2)
	{
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File,ios::app);
		if(!OutputFile.is_open())
			return(1);
	}
	else if(argc == 35)
	{
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File,ios::app);
		if(!OutputFile.is_open())
			return(1);

		Random_Range[0][0] = atof(argv[2]);
		Random_Range[0][1] = atof(argv[3]);
		Random_Range[1][0] = atof(argv[4]);
		Random_Range[1][1] = atof(argv[5]);
		Random_Range[2][0] = atof(argv[6]);
		Random_Range[2][1] = atof(argv[7]);
		Random_Range[3][0] = atof(argv[8]);
		Random_Range[3][1] = atof(argv[9]);
		Random_Range[4][0] = atof(argv[10]);
		Random_Range[4][1] = atof(argv[11]);
		Random_Range[5][0] = atof(argv[12]);
		Random_Range[5][1] = atof(argv[13]);
		Random_Range[6][0] = atof(argv[14]);
		Random_Range[6][1] = atof(argv[15]);
		Random_Range[7][0] = atof(argv[16]);
		Random_Range[7][1] = atof(argv[17]);
		Random_Range[8][0] = atof(argv[18]);
		Random_Range[8][1] = atof(argv[19]);
		Random_Range[9][0] = atof(argv[20]);
		Random_Range[9][1] = atof(argv[21]);
		Random_Range[10][0] = atof(argv[22]);
		Random_Range[10][1] = atof(argv[23]);
		Random_Range[11][0] = atof(argv[24]);
		Random_Range[11][1] = atof(argv[25]);
		Random_Range[12][0] = atof(argv[26]);
		Random_Range[12][1] = atof(argv[27]);
		Random_Range[13][0] = atof(argv[28]);
		Random_Range[13][1] = atof(argv[29]);
		Random_Range[14][0] = atof(argv[30]);
		Random_Range[14][1] = atof(argv[31]);
		Random_Range[15][0] = atof(argv[32]);
		Random_Range[15][1] = atof(argv[33]);
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
			JPsi[1]->Replace(Parameter_List[i][0], 0, 1);
			JPsi[1]->Replace(Parameter_List[i][1], 0, 2);
			JPsi[1]->Replace(Parameter_List[i][2], 1, 1);
			JPsi[1]->Replace(Parameter_List[i][3], 1, 2);
			JPsi[1]->Replace(Parameter_List[i][4], 2, 1);
			JPsi[1]->Replace(Parameter_List[i][5], 2, 2);
			JPsi[1]->Replace(Parameter_List[i][6], 0, 1);
			JPsi[1]->Replace(Parameter_List[i][7], 0, 2);
			JPsi[1]->Replace(Parameter_List[i][8], 1, 1);
			JPsi[1]->Replace(Parameter_List[i][9], 1, 2);

			Medium_Euclidean[0][0] = JPsi[1]->Euclidean(.5, 0)+Psi_Prime[1]->Euclidean(.5,0)+Non[1]->Euclidean(.5,0);
			Medium_Euclidean[0][1] = JPsi[1]->Euclidean(.5, 3)+Psi_Prime[1]->Euclidean(.5,3)+Non[1]->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[0][j] = JPsi[1]->Spatial((long double)(j)+.25)+Psi_Prime[1]->Spatial((long double)(j)+.25)+Non[1]->Spatial((long double)(j)+.25);
			Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
		}
		return(0);
	}

	long double Best[21], Chi;
	OutputFile << "Random Search seed: " << time(NULL)+3*atoi(argv[1]) << endl;
	srand(time(NULL)+3*atoi(argv[1]));
	time_t master_start_time = time(NULL);
	time_t round_start_time;
	OutputFile << "Random Search" << endl;
	int i = 0;
	if(atoi(argv[1])!=0 && Cycle)
	{
//A, PA, DeltaM1, DeltaM4, PM, Gamma, PGamma, a, Pa, Delta, PDelta, DeltaMQ1, DeltaPMQ4, PMQ, n, Pn
		Deviation_Points[0] = Uniform(Random_Range[0][0],Random_Range[0][1]);
		Deviation_Points[1] = Uniform(Random_Range[0][0],Random_Range[0][1]);
		Deviation_Points[2] = Uniform(Random_Range[1][0],Random_Range[1][1]);
		Deviation_Points[3] = Uniform(Random_Range[1][0],Random_Range[1][1]);
		Deviation_Points[4] = Uniform(Random_Range[2][0],Random_Range[2][1]);
		Deviation_Points[5] = Uniform(Random_Range[3][0],Random_Range[3][1]);
		Deviation_Points[6] = Uniform(Random_Range[4][0],Random_Range[4][1]);
		Deviation_Points[7] = Uniform(Random_Range[4][0],Random_Range[4][1]);
		Deviation_Points[8] = Uniform(Random_Range[5][0],Random_Range[5][1]);
		Deviation_Points[9] = Uniform(Random_Range[5][0],Random_Range[5][1]);
		Deviation_Points[10] = Uniform(Random_Range[6][0],Random_Range[6][1]);
		Deviation_Points[11] = Uniform(Random_Range[6][0],Random_Range[6][1]);
		Deviation_Points[12] = Uniform(Random_Range[11][0],Random_Range[11][1]);
		Deviation_Points[13] = Uniform(Random_Range[12][0],Random_Range[12][1]);
		Deviation_Points[14] = Uniform(Random_Range[13][0],Random_Range[13][1]);
		Deviation_Points[15] = Uniform(Random_Range[13][0],Random_Range[13][1]);
		Deviation_Points[16] = Uniform(Random_Range[14][0],Random_Range[14][1]);
		Deviation_Points[17] = Uniform(Random_Range[14][0],Random_Range[14][1]);
		Deviation_Points[18] = Uniform(Random_Range[15][0],Random_Range[15][1]);
		Deviation_Points[19] = Uniform(Random_Range[15][0],Random_Range[15][1]);

		for(int j = 0; j < 20; j++)
			Best[j] = Deviation_Points[j];
	}
	else
	{
		for(int j = 0; j < 20; j++)
			Best[j] = Deviation_Points[j];
	}

	DataLoad(JPsi, Non, Deviation_Points);

	for(int T = 1; T < 5; T++)
	{
		Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
		Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
		for(int j = 0; j < 7; j++)
			Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
	}
	Best[20] = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);

	round_start_time = time(NULL);
	while(difftime(time(NULL), round_start_time) < 18000 || i < 80) //18000 seconds (5 hours) and 80 attempts
	{
		Deviation_Points[0] = Protected_Uniform(Best[0],Random_Range[0][0],Random_Range[0][1],Best[20]);
		Deviation_Points[1] = Protected_Uniform(Best[1],Random_Range[0][0],Random_Range[0][1],Best[20]);
		Deviation_Points[2] = Protected_Uniform(Best[2],Random_Range[1][0],Random_Range[1][1],Best[20]);
		Deviation_Points[3] = Protected_Uniform(Best[3],Random_Range[1][0],Random_Range[1][1],Best[20]);
		Deviation_Points[4] = Protected_Uniform(Best[4],Random_Range[2][0],Random_Range[2][1],Best[20]);
		Deviation_Points[5] = Protected_Uniform(Best[5],Random_Range[3][0],Random_Range[3][1],Best[20]);
		Deviation_Points[6] = Protected_Uniform(Best[6],Random_Range[4][0],Random_Range[4][1],Best[20]);
		Deviation_Points[7] = Protected_Uniform(Best[7],Random_Range[4][0],Random_Range[4][1],Best[20]);
		Deviation_Points[8] = Protected_Uniform(Best[8],Random_Range[5][0],Random_Range[5][1],Best[20]);
		Deviation_Points[9] = Protected_Uniform(Best[9],Random_Range[5][0],Random_Range[5][1],Best[20]);
		Deviation_Points[10] = Protected_Uniform(Best[10],Random_Range[6][0],Random_Range[6][1],Best[20]);
		Deviation_Points[11] = Protected_Uniform(Best[11],Random_Range[6][0],Random_Range[6][1],Best[20]);
		Deviation_Points[12] = Protected_Uniform(Best[12],Random_Range[11][0],Random_Range[11][1],Best[20]);
		Deviation_Points[13] = Protected_Uniform(Best[13],Random_Range[12][0],Random_Range[12][1],Best[20]);
		Deviation_Points[14] = Protected_Uniform(Best[14],Random_Range[13][0],Random_Range[13][1],Best[20]);
		Deviation_Points[15] = Protected_Uniform(Best[15],Random_Range[13][0],Random_Range[13][1],Best[20]);
		Deviation_Points[16] = Protected_Uniform(Best[16],Random_Range[14][0],Random_Range[14][1],Best[20]);
		Deviation_Points[17] = Protected_Uniform(Best[17],Random_Range[14][0],Random_Range[14][1],Best[20]);
		Deviation_Points[18] = Protected_Uniform(Best[18],Random_Range[15][0],Random_Range[15][1],Best[20]);
		Deviation_Points[19] = Protected_Uniform(Best[19],Random_Range[15][0],Random_Range[15][1],Best[20]);

		DataLoad(JPsi, Non, Deviation_Points);

		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
			Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
		}

		Chi = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);

		if((!isnan(Chi) && Chi < Best[20]) || isnan(Best[20]))
		{
			for(int j = 0; j < 20; j++)
				Best[j] = Deviation_Points[j];
			Best[20] = Chi;
			i = 0;
		}
		else
			i++;
	}
	for(int j = 0; j < 20; j++)
		Deviation_Points[j] = Best[j];

	long double gradn_1[20], gradn[20], sn_1[20], sn[20];
	long double betan = 0, betan_1;
	round_start_time = time(NULL);

	Gradient(gradn_1, Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
	for(int i = 0; i < 20; i++)
		sn_1[i] = gradn_1[i];
	Minimize(sn_1, Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
	do
	{
		betan_1 = betan;
		Gradient(gradn, Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
		betan = PolakRibiere(gradn_1,gradn);
		if(betan < 0)
			for(int i = 0; i < 20; i++)
				sn[i] = gradn[i]+betan*sn_1[i];
		else
			for(int i = 0; i < 20; i++)
				sn[i] = gradn[i];
		Minimize(sn, Deviation_Points, JPsi, Psi_Prime, Non,  Medium_Euclidean, Medium_Spatial);
		for(int i = 0; i < 20; i++)
		{
			gradn_1[i] = gradn[i];
			sn_1[i] = sn[i];
		}
	}while((abs(betan)>1e-12 || abs(betan_1)>1e-12) && difftime(time(NULL), round_start_time) < 1800);

	OutputFile.close();

	return(0);
}

void Minimize(long double sn[20], long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], long double Medium_Euclidean[4][2], long double Medium_Spatial[4][7])
{
	long double fz[101][2];
	long double a = 0, c = 1;
	int Min_i;
	long double Local_Parameters[20];
	long double length[2][2];

	for(int i = 0; i < 20; i++)
		Local_Parameters[i] = Deviation_Points[i];

	length[0][0] = (sn[0]>0)?sn[0]/(Random_Range[0][1]-Local_Parameters[0]):sn[0]/(Random_Range[0][0]-Local_Parameters[0]);
	length[0][1] = (sn[0]>0)?sn[0]/(Random_Range[0][0]-Local_Parameters[0]):sn[0]/(Random_Range[0][1]-Local_Parameters[0]);
	length[1][0] = (sn[1]>0)?sn[1]/(Random_Range[0][1]-Local_Parameters[1]):sn[1]/(Random_Range[0][0]-Local_Parameters[1]);
	length[1][1] = (sn[1]>0)?sn[1]/(Random_Range[0][0]-Local_Parameters[1]):sn[1]/(Random_Range[0][1]-Local_Parameters[1]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[2]>0)?sn[2]/(Random_Range[1][1]-Local_Parameters[2]):sn[2]/(Random_Range[1][0]-Local_Parameters[2]);
	length[1][1] = (sn[2]>0)?sn[2]/(Random_Range[1][0]-Local_Parameters[2]):sn[2]/(Random_Range[1][1]-Local_Parameters[2]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[3]>0)?sn[3]/(Random_Range[1][1]-Local_Parameters[3]):sn[3]/(Random_Range[1][0]-Local_Parameters[3]);
	length[1][1] = (sn[3]>0)?sn[3]/(Random_Range[1][0]-Local_Parameters[3]):sn[3]/(Random_Range[1][1]-Local_Parameters[3]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[4]>0)?sn[4]/(Random_Range[2][1]-Local_Parameters[4]):sn[4]/(Random_Range[2][0]-Local_Parameters[4]);
	length[1][1] = (sn[4]>0)?sn[4]/(Random_Range[2][0]-Local_Parameters[4]):sn[4]/(Random_Range[2][1]-Local_Parameters[4]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[5]>0)?sn[5]/(Random_Range[3][1]-Local_Parameters[5]):sn[5]/(Random_Range[3][0]-Local_Parameters[5]);
	length[1][1] = (sn[5]>0)?sn[5]/(Random_Range[3][0]-Local_Parameters[5]):sn[5]/(Random_Range[3][1]-Local_Parameters[5]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[6]>0)?sn[6]/(Random_Range[4][1]-Local_Parameters[6]):sn[6]/(Random_Range[4][0]-Local_Parameters[6]);
	length[1][1] = (sn[6]>0)?sn[6]/(Random_Range[4][0]-Local_Parameters[6]):sn[6]/(Random_Range[4][1]-Local_Parameters[6]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[7]>0)?sn[7]/(Random_Range[4][1]-Local_Parameters[7]):sn[7]/(Random_Range[4][0]-Local_Parameters[7]);
	length[1][1] = (sn[7]>0)?sn[7]/(Random_Range[4][0]-Local_Parameters[7]):sn[7]/(Random_Range[4][1]-Local_Parameters[7]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[8]>0)?sn[8]/(Random_Range[5][1]-Local_Parameters[8]):sn[8]/(Random_Range[5][0]-Local_Parameters[8]);
	length[1][1] = (sn[8]>0)?sn[8]/(Random_Range[5][0]-Local_Parameters[8]):sn[8]/(Random_Range[5][1]-Local_Parameters[8]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[9]>0)?sn[9]/(Random_Range[5][1]-Local_Parameters[9]):sn[9]/(Random_Range[5][0]-Local_Parameters[9]);
	length[1][1] = (sn[9]>0)?sn[9]/(Random_Range[5][0]-Local_Parameters[9]):sn[9]/(Random_Range[5][1]-Local_Parameters[9]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[10]>0)?sn[10]/(Random_Range[6][1]-Local_Parameters[10]):sn[10]/(Random_Range[6][0]-Local_Parameters[10]);
	length[1][1] = (sn[10]>0)?sn[10]/(Random_Range[6][0]-Local_Parameters[10]):sn[10]/(Random_Range[6][1]-Local_Parameters[10]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[11]>0)?sn[11]/(Random_Range[6][1]-Local_Parameters[11]):sn[11]/(Random_Range[6][0]-Local_Parameters[11]);
	length[1][1] = (sn[11]>0)?sn[11]/(Random_Range[6][0]-Local_Parameters[11]):sn[11]/(Random_Range[6][1]-Local_Parameters[11]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[12]>0)?sn[12]/(Random_Range[11][1]-Local_Parameters[12]):sn[12]/(Random_Range[11][0]-Local_Parameters[12]);
	length[1][1] = (sn[12]>0)?sn[12]/(Random_Range[11][0]-Local_Parameters[12]):sn[12]/(Random_Range[11][1]-Local_Parameters[12]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[13]>0)?sn[13]/(Random_Range[12][1]-Local_Parameters[13]):sn[13]/(Random_Range[12][0]-Local_Parameters[13]);
	length[1][1] = (sn[13]>0)?sn[13]/(Random_Range[12][0]-Local_Parameters[13]):sn[13]/(Random_Range[12][1]-Local_Parameters[13]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[14]>0)?sn[14]/(Random_Range[13][1]-Local_Parameters[14]):sn[14]/(Random_Range[13][0]-Local_Parameters[14]);
	length[1][1] = (sn[14]>0)?sn[14]/(Random_Range[13][0]-Local_Parameters[14]):sn[14]/(Random_Range[13][1]-Local_Parameters[14]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[15]>0)?sn[15]/(Random_Range[13][1]-Local_Parameters[15]):sn[15]/(Random_Range[13][0]-Local_Parameters[15]);
	length[1][1] = (sn[15]>0)?sn[15]/(Random_Range[13][0]-Local_Parameters[15]):sn[15]/(Random_Range[13][1]-Local_Parameters[15]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[16]>0)?sn[16]/(Random_Range[14][1]-Local_Parameters[16]):sn[16]/(Random_Range[14][0]-Local_Parameters[16]);
	length[1][1] = (sn[16]>0)?sn[16]/(Random_Range[14][0]-Local_Parameters[16]):sn[16]/(Random_Range[14][1]-Local_Parameters[16]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[17]>0)?sn[17]/(Random_Range[14][1]-Local_Parameters[17]):sn[17]/(Random_Range[14][0]-Local_Parameters[17]);
	length[1][1] = (sn[17]>0)?sn[17]/(Random_Range[14][0]-Local_Parameters[17]):sn[17]/(Random_Range[14][1]-Local_Parameters[17]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[18]>0)?sn[18]/(Random_Range[15][1]-Local_Parameters[18]):sn[18]/(Random_Range[15][0]-Local_Parameters[18]);
	length[1][1] = (sn[18]>0)?sn[18]/(Random_Range[15][0]-Local_Parameters[18]):sn[18]/(Random_Range[15][1]-Local_Parameters[18]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[1][0] = (sn[19]>0)?sn[19]/(Random_Range[15][1]-Local_Parameters[19]):sn[19]/(Random_Range[15][0]-Local_Parameters[19]);
	length[1][1] = (sn[19]>0)?sn[19]/(Random_Range[15][0]-Local_Parameters[19]):sn[19]/(Random_Range[15][1]-Local_Parameters[19]);
	length[0][0] = (length[0][0]>length[1][0])?length[0][0]:length[1][0];
	length[0][1] = (length[0][1]<length[1][1])?length[0][1]:length[1][1];
	length[0][0] = 1./length[0][0];
	length[0][1] = 1./length[0][1];

	OutputFile << "Line search" << endl;
	for(int i = 0; i <= 100; i++)
	{
		fz[i][0] = a+(c-a)*i/100.;
		for(int j = 0; j < 20; j++)
			Deviation_Points[j] = Local_Parameters[j]+length[0][0]*fz[i][0]*sn[j];
		DataLoad(JPsi, Non, Deviation_Points);
		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
			Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
		}
		fz[i][1] = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
	}

	for(int i = 1; i < 101; i++)
		Min_i = (fz[i][1]<fz[Min_i][1])?i:Min_i;
	if(Min_i == 0)	//So as to not exceed the limits of the array
		Min_i = 1;
	else if(Min_i == 100)
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
			for(int j = 0; j < 20; j++)
				Deviation_Points[j] = Local_Parameters[j]+length[0][0]*x*sn[j];
			DataLoad(JPsi, Non, Deviation_Points);
			for(int T = 1; T < 5; T++)
			{
				Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
				Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
				for(int j = 0; j < 7; j++)
					Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
			}
			Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
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
		for(int j = 0; j < 20; j++)
			Deviation_Points[j] = Local_Parameters[j]+length[0][0]*u*sn[j];
		DataLoad(JPsi, Non, Deviation_Points);
		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
			Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
		}
		fu = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);

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

long double PolakRibiere(long double gradn_1[20], long double gradn[20])
{
	long double num = 0, den = 0;

	for(int i = 0; i < 20; i++)
	{
		num += gradn[i]*(gradn[i]-gradn_1[i]);
		den += pow(gradn_1[i],2);
	}

	return(num/den);
}

void Gradient(long double grad[20], long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], long double Medium_Euclidean[4][2], long double Medium_Spatial[4][7])
{
	long double f0, f1, h = 1e-5;
	long double Reduce_Euclidean[4][3][2];
	long double Reduce_Spatial[4][3][7];
	long double norm = 0;
	int i, j;

	OutputFile << "Gradient" << endl;
	for(int T = 1; T < 5; T++)
	{
		Reduce_Euclidean[T-1][0][0] = JPsi[T]->Euclidean(.5, 0);
		Reduce_Euclidean[T-1][1][0] = Psi_Prime[T]->Euclidean(.5, 0);
		Reduce_Euclidean[T-1][2][0] = Non[T]->Euclidean(.5, 0);
		Medium_Euclidean[T-1][0] = Reduce_Euclidean[T-1][0][0]+Reduce_Euclidean[T-1][1][0]+Reduce_Euclidean[T-1][2][0];

		Reduce_Euclidean[T-1][0][1] = JPsi[T]->Euclidean(.5, 3);
		Reduce_Euclidean[T-1][1][1] = Psi_Prime[T]->Euclidean(.5, 3);
		Reduce_Euclidean[T-1][2][1] = Non[T]->Euclidean(.5, 3);
		Medium_Euclidean[T-1][1] = Reduce_Euclidean[T-1][0][1]+Reduce_Euclidean[T-1][1][1]+Reduce_Euclidean[T-1][2][1];
		for(int j = 0; j < 7; j++)
		{
			Reduce_Spatial[T-1][0][j] = JPsi[T]->Spatial((long double)(j)+.25);
			Reduce_Spatial[T-1][1][j] = Psi_Prime[T]->Spatial((long double)(j)+.25);
			Reduce_Spatial[T-1][2][j] = Non[T]->Spatial((long double)(j)+.25);
			Medium_Spatial[T-1][j] = Reduce_Spatial[T-1][0][j]+Reduce_Spatial[T-1][1][j]+Reduce_Spatial[T-1][2][j];
		}
	}
	f0 = Chi_Square(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);
	//cerr << "Gradient " << f0 << endl;

	for(i = 0; i < 12; i++)
	{
		Deviation_Points[i] += h;
		DataLoad(JPsi, Non, Deviation_Points);
		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Reduce_Euclidean[T-1][1][0]+Reduce_Euclidean[T-1][2][0];
			Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Reduce_Euclidean[T-1][1][0]+Reduce_Euclidean[T-1][2][1];
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Reduce_Spatial[T-1][1][j]+Reduce_Spatial[T-1][2][j];
		}
		f1 = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
		Deviation_Points[i] -= h;
		grad[i] = (f0-f1)/h;
	}

	for(i; i < 20; i++)
	{
		Deviation_Points[i] += h;
		DataLoad(JPsi, Non, Deviation_Points);
		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = Reduce_Euclidean[T-1][0][0]+Reduce_Euclidean[T-1][1][0]+Non[T]->Euclidean(.5, 0);
			Medium_Euclidean[T-1][1] = Reduce_Euclidean[T-1][0][1]+Reduce_Euclidean[T-1][1][1]+Non[T]->Euclidean(.5, 3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = Reduce_Spatial[T-1][0][j]+Reduce_Spatial[T-1][1][j]+Non[T]->Spatial((long double)(j)+.25);
		}
		f1 = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial);
		Deviation_Points[i] -= h;
		grad[i] = (f0-f1)/h;
	}
}

long double Print(long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], long double Medium_Euclidean[4][2], long double Medium_Spatial[4][7])
{
	long double chi[5];
	for(int i = 0; i < 5; i++)
		chi[i] = Chi_Square(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, i);

	for(int i = 0; i < 20; i++)
	{
		OutputFile << Deviation_Points[i] << "," << flush;
	}

	for(int i = 1; i < 5; i++)
	{
		OutputFile << i << "," << flush;
		JPsi[i]->Print(OutputFile);
		OutputFile << ",";
		Non[i]->Print(OutputFile);
		OutputFile << "," << Medium_Euclidean[i-1][1]/Vacuum_Euclidean[i-1][1]-Medium_Euclidean[i-1][0]/Vacuum_Euclidean[i-1][0] << "," << flush;
		for(int j = 0; j < 7; j++)
		{
			OutputFile << Medium_Spatial[i-1][j]/Vacuum_Spatial[j] << "," << flush;
		}
		OutputFile << chi[i] << "," << flush;
	}

	OutputFile << chi[0] << endl;
	return(chi[0]);
}

long double Chi_Square(Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], long double Medium_Euclidean[4][2], long double Medium_Spatial[4][7], int Temp)
{
	long double answer = 0;

	if(Temp != 0)	//To report out each temprature for goodness of fit considerations
	{
		answer += pow(Medium_Euclidean[Temp-1][1]/Vacuum_Euclidean[Temp-1][1]-Medium_Euclidean[Temp-1][0]/Vacuum_Euclidean[Temp-1][0]-.2,2)/.2;
		for(int j = 0; j < 7; j++)
			answer += pow(Medium_Spatial[Temp-1][j]/Vacuum_Spatial[j]-Spatial_Ratio[Temp-1][j],2)/Spatial_Ratio[Temp-1][j];
	}
	else	//Total Chi-squared
	{
		for(int i = 0; i < 4; i++)
			answer += pow(Medium_Euclidean[i][1]/Vacuum_Euclidean[i][1]-Medium_Euclidean[i][0]/Vacuum_Euclidean[i][0]-.2,2)/.2;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 7; j++)
				answer += pow(Medium_Spatial[i][j]/Vacuum_Spatial[j]-Spatial_Ratio[i][j],2)/Spatial_Ratio[i][j];
	}

	/*for(int i = 0; i < 5; i++)
		for(int j = 1; j < 2; j++)
			answer += Least_Squares(JPsi[1]->Read(i,j), JPsi[2]->Read(i,j), JPsi[3]->Read(i,j), JPsi[4]->Read(i,j));
	for(int i = 0; i < 2; i++)
		for(int j = 1; j < 2; j++)
			answer += Least_Squares(Non[1]->Read(i,j), Non[2]->Read(i,j), Non[3]->Read(i,j), Non[4]->Read(i,j));*/

	return(answer);
}

void DataLoad(Spectral_Inter* JPsi[5], Spectral_Non* Non[5], long double Deviation_Points[20])
{
	JPsi[1]->Replace(Deviation_Points[0],0,1);
	JPsi[1]->Replace(Deviation_Points[2],0,2);
	JPsi[1]->Replace(3.09946+Deviation_Points[4],1,1);
	JPsi[1]->Replace(Deviation_Points[6],1,2);
	JPsi[1]->Replace(Deviation_Points[8],2,1);
	JPsi[1]->Replace(Deviation_Points[10],2,2);
	Non[1]->Replace(1.655+Deviation_Points[12],0,1);
	Non[1]->Replace(Deviation_Points[14],0,2);
	Non[1]->Replace(Deviation_Points[16],1,1);
	Non[1]->Replace(Deviation_Points[18],1,2);
	JPsi[2]->Replace((71.*Deviation_Points[0]+32.*Deviation_Points[1])/103.,0,1);
	JPsi[2]->Replace((71.*Deviation_Points[2]+32.*Deviation_Points[3])/103.,0,2);
	JPsi[2]->Replace(3.125+(71.*Deviation_Points[4]+32.*Deviation_Points[5])/103.,1,1);
	JPsi[2]->Replace((71.*Deviation_Points[6]+32.*Deviation_Points[7])/103.,1,2);
	JPsi[2]->Replace((71.*Deviation_Points[8]+32.*Deviation_Points[9])/103.,2,1);
	JPsi[2]->Replace((71.*Deviation_Points[10]+32.*Deviation_Points[11])/103.,2,2);
	Non[2]->Replace(1.59+(71.*Deviation_Points[12]+32.*Deviation_Points[13])/103.,0,1);
	Non[2]->Replace((71.*Deviation_Points[14]+32.*Deviation_Points[15])/103.,0,2);
	Non[2]->Replace((71.*Deviation_Points[16]+32.*Deviation_Points[17])/103.,1,1);
	Non[2]->Replace((71.*Deviation_Points[18]+32.*Deviation_Points[19])/103.,1,2);
	JPsi[3]->Replace((40.*Deviation_Points[0]+63.*Deviation_Points[1])/103.,0,1);
	JPsi[3]->Replace((40.*Deviation_Points[2]+63.*Deviation_Points[3])/103.,0,2);
	JPsi[3]->Replace(3.151+(40.*Deviation_Points[4]+63.*Deviation_Points[5])/103.,1,1);
	JPsi[3]->Replace((40.*Deviation_Points[6]+63.*Deviation_Points[7])/103.,1,2);
	JPsi[3]->Replace((40.*Deviation_Points[8]+63.*Deviation_Points[9])/103.,2,1);
	JPsi[3]->Replace((40.*Deviation_Points[10]+63.*Deviation_Points[11])/103.,2,2);
	Non[3]->Replace(1.51+(40.*Deviation_Points[12]+63.*Deviation_Points[13])/103.,0,1);
	Non[3]->Replace((40.*Deviation_Points[14]+63.*Deviation_Points[15])/103.,0,2);
	Non[3]->Replace((40.*Deviation_Points[16]+63.*Deviation_Points[17])/103.,1,1);
	Non[3]->Replace((40.*Deviation_Points[18]+63.*Deviation_Points[19])/103.,1,2);
	JPsi[4]->Replace(Deviation_Points[1],0,1);
	JPsi[4]->Replace(Deviation_Points[3],0,2);
	JPsi[4]->Replace(3.1855+Deviation_Points[5],1,1);
	JPsi[4]->Replace(Deviation_Points[7],1,2);
	JPsi[4]->Replace(Deviation_Points[9],2,1);
	JPsi[4]->Replace(Deviation_Points[11],2,2);
	Non[4]->Replace(1.36+Deviation_Points[13],0,1);
	Non[4]->Replace(Deviation_Points[15],0,2);
	Non[4]->Replace(Deviation_Points[17],1,1);
	Non[4]->Replace(Deviation_Points[19],1,2);
}

long double Protected_Uniform(long double x0, long double a, long double b, long double chi)
{
	long double Delta = (b-a)*pow(chi,.25)/100.;
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

long double Least_Squares(long double y1, long double y2, long double y3, long double y4)
{
//cout << "Least_Squares " << 0.030408*pow(y1,2)+ 0.064712*pow(y2,2)+ 0.066696*pow(y3,2)+ 0.023816*pow(y4,2)- 0.074128*y1*y2- 0.025024*y1*y3- 0.038848*y2*y3+ 0.038336*y1*y4- 0.016448*y2*y4- 0.06952*y3*y4 << " relative error squared " << pow(y1-1.7459529462551269*y2+0.5240664796028492*y3+0.2218864666522773*y4,2)/pow(y1+0.7582559896395411*y2+0.5240664796028492*y3+0.22188646665227693*y4,2)+(0.23740790925719363*pow(y1-1.2188897658510927*y2-0.4114706656143129*y3+0.6303604314654047*y4,2))/pow(y1+0.5938982181771568*y2+0.20048711703627758*y3-0.3071401102422766*y4,2)+pow(y1-0.42904841402337196*y2-1.8134390651085137*y3+1.2424874791318834*y4,2)/pow(y1-0.42904841402337235*y2-1.813439065108513*y3-3.5997495826377293*y4,2)+pow(y1+1.5524296675191809*y2-5.330562659846552*y3+2.7781329923273677*y4,2)/pow(y1+1.5524296675191813*y2+2.0875959079283897*y3+2.7781329923273645*y4,2) << endl;
	return(0.030408*pow(y1,2)+ 0.064712*pow(y2,2)+ 0.066696*pow(y3,2)+ 0.023816*pow(y4,2)- 0.074128*y1*y2- 0.025024*y1*y3- 0.038848*y2*y3+ 0.038336*y1*y4- 0.016448*y2*y4- 0.06952*y3*y4);
	//return(pow(y1-1.7459529462551269*y2+0.5240664796028492*y3+0.2218864666522773*y4,2)/pow(y1+0.7582559896395411*y2+0.5240664796028492*y3+0.22188646665227693*y4,2)+(0.23740790925719363*pow(y1-1.2188897658510927*y2-0.4114706656143129*y3+0.6303604314654047*y4,2))/pow(y1+0.5938982181771568*y2+0.20048711703627758*y3-0.3071401102422766*y4,2)+pow(y1-0.42904841402337196*y2-1.8134390651085137*y3+1.2424874791318834*y4,2)/pow(y1-0.42904841402337235*y2-1.813439065108513*y3-3.5997495826377293*y4,2)+pow(y1+1.5524296675191809*y2-5.330562659846552*y3+2.7781329923273677*y4,2)/pow(y1+1.5524296675191813*y2+2.0875959079283897*y3+2.7781329923273645*y4,2));
}
