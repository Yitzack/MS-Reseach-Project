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

void Gradient(long double[20], long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], pair<long double,long double>[4][2], pair<long double,long double>[4][7]);	//Calculates gradient for conjugate gradient descent
long double PolakRibiere(long double[20], long double[20]);	//Calculation of Polak-Ribiere conjugate gradient dot product thingy
void Minimize(long double[20], long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], pair<long double,long double>[4][2], pair<long double,long double>[4][7]);	//Minimization algorithm switching between Brent's Method and Parabolic Interpolation searches along a line

long double Chi_Square(Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], pair<long double,long double>[4][2], pair<long double,long double>[4][7], int); //Calculates a chi-square like thing for finding best fits
long double Least_Squares(long double, long double, long double, long double);	//Scores the difference betweeen 4 points and linear dependence on temprature
long double Print(long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], pair<long double,long double>[4][2], pair<long double,long double>[4][7], int); //In addition to printing the parameters, Euclidean difference, Spatial correlator, and Chi-Square, it also returns Chi_Square(), basically as an alias for Chi_Square
long double Print(long double[20], Spectral_Inter*[5], Spectral_Inter*[5], Spectral_Non*[5], pair<long double,long double>[4][2], pair<long double,long double>[4][7], pair<long double, long double>[7], int); //Like the print above but it get an extra array of Lorentz spatial correlations, for use in one place because I don't want in the optimiser
long double Uniform(long double, long double);	//Uniform random distrubution between the inputs
long double Protected_Uniform(long double, long double, long double, long double);	//Calls the uniform random above while making sure it doesn't wonder out of an acceptable window
void DataLoad(Spectral_Inter*[5], Spectral_Non*[5], long double[20]);	//Complete override of data so that parameters have certain constraints and are linear

				   //A, PA, DeltaM1, DeltaM4, PM, Gamma, PGamma, a, Pa, Delta, PDelta, DeltaMQ1, DeltaPMQ4, PMQ, n, Pn
long double Random_Range[16][2] = {{-.25,.25},{1.,6.},{-0.59946, 0.40054},{-0.6855, 0.3145},{1.,6.},{.02,.18},{1.,6.},{5.,15.},{1.,6.},{1.5,3.5},{1.,6.},{-.065,.135},{-.05,.43},{1.,6.},{1.5,5.},{1.,6.}};
ofstream OutputFile;

const long double Spatial_Ratio[4][7] = {{1.,1.00006,0.99883,0.992039,0.982366,0.970341,0.95766},	//Target Spatial Correlation Ratios
				   {.99,0.988286,0.945063,0.879461,0.798659,0.7259,0.654381},
				   {.98,0.954875,0.856416,0.720447,0.573465,0.45867,0.376707},
				   {.97,0.908029,0.715435,0.524036,0.372788,0.246218,0.18}};
const pair<long double,long double> Vacuum_Spatial[7] = {pair<long double,long double>(13.5965519874368885,4.02179730192699146e-07),
	pair<long double,long double>(0.0415680226812305554,1.1556052202630816e-08),pair<long double,long double>(0.0012012677483247847,3.32066969108868755e-10),
	pair<long double,long double>(4.6499993403302829e-05,9.54353876931880564e-12),pair<long double,long double>(1.96859078869344768e-06,2.74386862111523716e-13),
	pair<long double,long double>(8.66288828408022226e-08,7.89830973612391187e-15),pair<long double,long double>(3.89562525681864189e-09,2.62903068734599143e-16)};	//Pre-calculated vacuum spatial correlation function
const pair<long double,long double> Vacuum_Euclidean[4][2] = 
	{{pair<long double,long double>(0.000264166718975248739,2.42154874798803876e-07),pair<long double,long double>(9.71945863898214921e-06,7.60261186227071115e-09)},
	{pair<long double,long double>(0.00221555204564226004,1.75479630470614913e-06),pair<long double,long double>(0.000188270934911146417,1.19763944939391319e-07)},
	{pair<long double,long double>(0.00824183567291781835,5.61577871796860943e-06),pair<long double,long double>(0.00116036845723485179,6.04705149983612623e-07)},
	{pair<long double,long double>(0.0264511436195479036,1.47908592529210127e-05),pair<long double,long double>(0.00572333060820024838,2.32707330240952598e-06)}};	//Pre-calculated vacuum Euclidean correlation function at tau=T/2

int main(int argc, char* argv[])
{
	//long double Deviation_Points[20] = {A1, A4, PA1, PA4, DeltaM1, DeltaM4, PM1, PM4, Gamma1, Gamma4, PGamma1, PGamma4, DeltaMQ1, DeltaMQ4, PMQ1, PMQ4, n1, n4, Pn1, Pn4};{0.486081, 0.155585, 1.1562, 3.61439, -0.180868, -0.151585, 4.05393, 4.08588, 0.0470915, 0.0790153, 5.62678, 5.68307, -0.00134413, 0.00410448, 3.02722, 1.65693, 2.135601, 4.283871, 2.885071, 1.396201, 1, 0.486081, 1.1562, 2.91859, 4.05393, 0.0470915, 5.62678, 10.58, 4.00927, 3, 1, 1.65366, 3.02722, 2.135601, 2.885071, -0.0925551, 1.01848, 1.25656, 1.41628, 1.53576, 1.90257, -3.580289, 0.8759678, 23.17668, 2, 0.383403, 1.91991, 2.95323, 4.06386, 0.0570096, 5.64427, 8.97, 2.61299, 3, 2.22449, 1.59035, 2.6015, 2.803022, 2.422511, 0.0314531, 0.984802, 0.797222, 0.726317, 0.724347, 0.6775956, 0.6065823, 0.538415, 0.31551, 3, 0.283933, 2.65975, 2.98804, 4.07347, 0.0666177, 5.66121, 9.07, 5.68392, 3, 3.59181, 1.51199, 2.18908, 3.449591, 1.974402, 0.0511807, 0.979265, 0.730733, 0.640868, 0.629152, 0.518963, -0.1865786, 4.49366, 46.1355, 4, 0.155585, 3.61439, 3.03391, 4.08588, 0.0790153, 5.68307, 9.52, 1.00741, 3, 1.00741, 1.3641, 1.65693, 4.283871, 1.396201, 0.0632518, 0.976794, 0.702607, 0.61441, 0.641127, 0.674582, 0.479097, 3.44028, 59.6975, 129.325}
	//long double Deviation_Points[20] = {0.486081, 0.155585, 1.1562, 3.61439, -0.180868, -0.151585, 4.05393, 4.08588, 0.0470915, 0.0790153, 5.62678, 5.68307, -0.00134413, 0.00410448, 3.02722, 1.65693, 2.135601, 4.283871, 2.885071, 1.396201};
	//long double Deviation_Points[20] = {DeltaA1, DeltaA4, PA1, PA4, DeltaM1, DeltaM4, PM1, PM4, Gamma1, Gamma4, PGamma1, PGamma4, DeltaMQ1, DeltaMQ4, PMQ1, PMQ4, n1, n4, Pn1, Pn4};{-0.0553617, 0.0819349, 3.85981, 4.41183, -0.024432, -0.558869, 3.26666, 2.22165, 0.163213, 0.132118, 5.72388, 5.4357, 0.0530835, 0.367393, 3.15481, 3.49742, 3.62434, 4.68273, 4.34571, 5.737, 1, 0.262436, 3.85981, 3.07503, 3.26666, 0.163213, 5.72388, 10.58, 4.00927, 3, 1, 1.70808, 3.15481, 3.62434, 4.34571, -0.0816804, 1.00253, 0.962636, 1.01247, 1.05659, 1.05622, 1.02796, 1.17931, 0.462785, 2, 0.389694, 4.03131, 2.93453, 2.942, 0.153553, 5.63434, 8.97, 2.61299, 3, 2.22449, 1.74073, 3.26125, 3.95316, 4.77796, 0.0912522, 1.00532, 0.988231, 1.04288, 0.909692, 0.653854, 0.457874, 0.735357, 0.205771, 3, 0.431663, 4.19746, 2.79968, 2.62748, 0.144194, 5.54761, 9.07, 5.68392, 3, 3.59181, 1.75533, 3.36437, 4.27171, 5.19669, 0.0831779, 1.00333, 1.06793, 1.12703, 0.758354, 0.154638, -0.155333, 0.576148, 1.4031, 4, 0.348769, 4.41183, 2.62663, 2.22165, 0.132118, 5.4357, 9.52, 1.00741, 3, 1.00741, 1.72739, 3.49742, 4.68273, 5.737, 0.0194521, 0.999636, 1.09671, 1.23878, 0.88124, 0.296777, 0.469896, 0.128582, 1.0628, 4.97934}
	long double Deviation_Points[20] = {-.00741034, .0827412, 4.00958, 4.99877, -.474984, -2.10508, 5.20832, 5.04324, 0.131089, 0.149966, 5.72388, 5.4357, .120144, 0.164985, 4.09329, 3.64841, 3.29008, 3.95769, 4.3464, 5.50403};

	long double JPsi_Parameters[5][5][3] = {{{.314831,.314831,1.},{3.0969,3.0969,1},{.032,.032,1},{9.34,9.34,1},{1,1,1}},
						{{1.97/(2.*3.09946),.294516,2.87077},{3.09946,2.61272,5.15654},{.106597,.102616,3.7999},{10.58,9.70375,4.78732},{3,2.97949,1.96609}},
						{{.4024,.454156,5.87388},{3.125,2.16747,5.26226},{.372,.171732,5.83927},{8.97,14.7462,3.51337},{3,2.83876,4.14535}},
						{{.403047,.43634,4.78343},{3.151,1.66683,5.0539},{.68,.15422,5.69145},{9.07,10.5457,3.25066},{3,2.78677,2.26082}},
						{{.266834,.348726,4.41179},{3.1855,2.6266,2.22161},{.98,.132075,5.43566},{9.52,9.52,1.00741},{3,3,1.00741}}};
	long double PsiPrime_Parameters[5][5][3] = {{{.150566,.150566,3.1},{3.6861,3.6861,3.6},{.102,.102,3.5},{9.34,9.34,4.57122},{1,1,5.25}},
						    {{.55/(2.*3.785),.55/(2.*3.785),3.1},{3.785,3.785,3.6},{.43,.43,3.5},{10.58,9.70375,4.78732},{1,1,5.25}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}},
						    {{0,0,1},{0,0,1},{0,0,1},{0,0,1},{0,0,1}}};
	long double Non_Parameters[5][2][3] = {{{1.775,1.775,3.09},{2.41182,2.41182,5.5}},
					       {{1.655,1.77989,3.65803},{2.45,4.13782,5.66908}},
					       {{1.59,1.71443,4.64476},{2.7,2.71014,2.31858}},
					       {{1.51,1.66247,3.71779},{2.4,2.91303,5.88634}},
					       {{1.36,1.72735,3.49738},{1.98,4.68269,5.73696}}};

	JPsi_Parameters[1][0][1] = JPsi_Parameters[1][0][0]+Deviation_Points[0];
	JPsi_Parameters[2][0][1] = JPsi_Parameters[2][0][0]+(71.*Deviation_Points[0]+32.*Deviation_Points[1])/103.;
	JPsi_Parameters[3][0][1] = JPsi_Parameters[3][0][0]+(40.*Deviation_Points[0]+63.*Deviation_Points[1])/103.;
	JPsi_Parameters[4][0][1] = JPsi_Parameters[4][0][0]+Deviation_Points[1];
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
		JPsi[i] = new Spectral_Inter(JPsi_Parameters[i], i, !bool(i));
		Psi_Prime[i] = new Spectral_Inter(PsiPrime_Parameters[i], i, !bool(i));
		Non[i] = new Spectral_Non(Non_Parameters[i], i, !bool(i));
	}

	pair<long double,long double> Medium_Spatial[4][7];
	pair<long double,long double> Medium_Euclidean[4][2];
	bool Cycle = true;
	if(argc != 4)
	{
		for(int T = 1; T < 5; T++)
		{
			Medium_Euclidean[T-1][0] = JPsi[T]->Euclidean(.5, 0)+Psi_Prime[T]->Euclidean(.5,0)+Non[T]->Euclidean(.5,0);
			Medium_Euclidean[T-1][1] = JPsi[T]->Euclidean(.5, 3)+Psi_Prime[T]->Euclidean(.5,3)+Non[T]->Euclidean(.5,3);
			for(int j = 0; j < 7; j++)
				Medium_Spatial[T-1][j] = JPsi[T]->Spatial((long double)(j)+.25)+Psi_Prime[T]->Spatial((long double)(j)+.25)+Non[T]->Spatial((long double)(j)+.25);
		}
	}
	/*for(int i = 0; i < 7; i++)	//Superceeded by precalculated values, standing by if services required
		Vacuum_Spatial[i] = JPsi[0]->Spatial_Lorentz((long double)(i)+.25)+Psi_Prime[0]->Spatial_Lorentz((long double)(i)+.25)+Non[0]->Spatial_Lorentz((long double)(i)+.25);
	Vacuum_Euclidean[0][0] = JPsi[0]->Euclidean(.5,0,.194)+Psi_Prime[0]->Euclidean(.5,0,.194)+Non[0]->Euclidean(.5,0,.194);
	Vacuum_Euclidean[1][0] = JPsi[0]->Euclidean(.5,0,.258)+Psi_Prime[0]->Euclidean(.5,0,.258)+Non[0]->Euclidean(.5,0,.258);
	Vacuum_Euclidean[2][0] = JPsi[0]->Euclidean(.5,0,.320)+Psi_Prime[0]->Euclidean(.5,0,.320)+Non[0]->Euclidean(.5,0,.320);
	Vacuum_Euclidean[3][0] = JPsi[0]->Euclidean(.5,0,.400)+Psi_Prime[0]->Euclidean(.5,0,.400)+Non[0]->Euclidean(.5,0,.400);
	Vacuum_Euclidean[0][1] = JPsi[0]->Euclidean(.5,3,.194)+Psi_Prime[0]->Euclidean(.5,3,.194)+Non[0]->Euclidean(.5,3,.194);
	Vacuum_Euclidean[1][1] = JPsi[0]->Euclidean(.5,3,.258)+Psi_Prime[0]->Euclidean(.5,3,.258)+Non[0]->Euclidean(.5,3,.258);
	Vacuum_Euclidean[2][1] = JPsi[0]->Euclidean(.5,3,.320)+Psi_Prime[0]->Euclidean(.5,3,.320)+Non[0]->Euclidean(.5,3,.320);
	Vacuum_Euclidean[3][1] = JPsi[0]->Euclidean(.5,3,.400)+Psi_Prime[0]->Euclidean(.5,3,.400)+Non[0]->Euclidean(.5,3,.400);
	cout << setprecision(LDBL_DIG) << "long double digits:=" << LDBL_DIG << endl;
	for(int i = 0; i < 7; i++)
		cout << Vacuum_Spatial[i].first << "±" << Vacuum_Spatial[i].second << endl;
	cout << endl << Vacuum_Euclidean[0][0].first << "±" << Vacuum_Euclidean[0][0].second << endl;
	cout << Vacuum_Euclidean[1][0].first << "±" << Vacuum_Euclidean[1][0].second << endl;
	cout << Vacuum_Euclidean[2][0].first << "±" << Vacuum_Euclidean[2][0].second << endl;
	cout << Vacuum_Euclidean[3][0].first << "±" << Vacuum_Euclidean[3][0].second << endl << endl;
	cout << Vacuum_Euclidean[0][1].first << "±" << Vacuum_Euclidean[0][1].second << endl;
	cout << Vacuum_Euclidean[1][1].first << "±" << Vacuum_Euclidean[1][1].second << endl;
	cout << Vacuum_Euclidean[2][1].first << "±" << Vacuum_Euclidean[2][1].second << endl;
	cout << Vacuum_Euclidean[3][1].first << "±" << Vacuum_Euclidean[3][1].second << endl;*/

	char File[70] = "data/Optimiser_Output";
	if(argc == 2)	//Default parameter search and chi-square minimization ./Optimizer ProcessID
	{
		strcat(File,".");
		strcat(File,argv[1]);
		strcat(File,".csv");
		OutputFile.open(File,ios::app);
		if(!OutputFile.is_open())
			return(1);
	}
	else if(argc == 35)	//Parameter search with search boundary overwritten ./Optimizer ProcessID <16 sets of 2 parameter boundaries>
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
	Best[20] = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);

	round_start_time = time(NULL);
	while(difftime(time(NULL), round_start_time) < 43200 || i < 80) //18000 seconds (5 hours)
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

		Chi = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);

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

void Minimize(long double sn[20], long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], pair<long double,long double> Medium_Euclidean[4][2], pair<long double,long double> Medium_Spatial[4][7])
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
		fz[i][1] = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);
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
			Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);
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
		fu = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);

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

void Gradient(long double grad[20], long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], pair<long double,long double> Medium_Euclidean[4][2], pair<long double,long double> Medium_Spatial[4][7])
{
	long double f0, f1, h = 1e-5;
	pair<long double,long double> Reduce_Euclidean[4][3][2];
	pair<long double,long double> Reduce_Spatial[4][3][7];
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
		f1 = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);
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
		f1 = Print(Deviation_Points, JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, 0);
		Deviation_Points[i] -= h;
		grad[i] = (f0-f1)/h;
	}
}

long double Print(long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], pair<long double,long double> Medium_Euclidean[4][2], pair<long double,long double> Medium_Spatial[4][7], int T)
{
	long double chi[5];
	if(T == 0)
	{
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
			OutputFile << "," << Medium_Euclidean[i-1][1].first/Vacuum_Euclidean[i-1][1].first-Medium_Euclidean[i-1][0].first/Vacuum_Euclidean[i-1][0].first << "±" << sqrt(pow(Medium_Euclidean[i-1][1].second/Vacuum_Euclidean[i-1][1].first,2)+pow(Medium_Euclidean[i-1][0].second/Vacuum_Euclidean[i-1][0].first,2)+pow(Medium_Euclidean[i-1][1].first*Vacuum_Euclidean[i-1][1].second/pow(Vacuum_Euclidean[i-1][0].first,2),2)+pow(Medium_Euclidean[i-1][0].first*Vacuum_Euclidean[i-1][0].second/pow(Vacuum_Euclidean[i-1][0].first,2),2)) << "," << flush;
			for(int j = 0; j < 7; j++)
			{
				OutputFile << Medium_Spatial[i-1][j].first/Vacuum_Spatial[j].first << "±" << sqrt(pow(Medium_Spatial[i-1][j].second/Vacuum_Spatial[j].first,2)+pow(Medium_Spatial[i-1][j].first*Vacuum_Spatial[j].second/pow(Vacuum_Spatial[j].first,2),2)) << "," << flush;
			}
			OutputFile << chi[i] << "," << flush;
		}

		OutputFile << chi[0] << endl;
	}
	else
	{
		chi[T] = Chi_Square(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, T);

		OutputFile << T << "," << flush;
		JPsi[T]->Print(OutputFile);
		OutputFile << ",";
		Non[T]->Print(OutputFile);
		OutputFile << "," << Medium_Euclidean[T-1][1].first/Vacuum_Euclidean[T-1][1].first-Medium_Euclidean[T-1][0].first/Vacuum_Euclidean[T-1][0].first << "±" << sqrt(pow(Medium_Euclidean[T-1][1].second/Vacuum_Euclidean[T-1][1].first,2)+pow(Medium_Euclidean[T-1][0].second/Vacuum_Euclidean[T-1][0].first,2)+pow(Medium_Euclidean[T-1][1].first*Vacuum_Euclidean[T-1][1].second/pow(Vacuum_Euclidean[T-1][0].first,2),2)+pow(Medium_Euclidean[T-1][0].first*Vacuum_Euclidean[T-1][0].second/pow(Vacuum_Euclidean[T-1][0].first,2),2)) << "," << flush;
		for(int j = 0; j < 7; j++)
		{
			OutputFile << Medium_Spatial[T-1][j].first/Vacuum_Spatial[j].first << "±" << sqrt(pow(Medium_Spatial[T-1][j].second/Vacuum_Spatial[j].first,2)+pow(Medium_Spatial[T-1][j].first*Vacuum_Spatial[j].second/pow(Vacuum_Spatial[j].first,2),2)) << "," << flush;
		}
		OutputFile << chi[T] << endl;
	}
	return(chi[0]);
}

long double Print(long double Deviation_Points[20], Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], pair<long double, long double> Medium_Euclidean[4][2], pair<long double, long double> Medium_Spatial[4][7], pair<long double, long double> Medium_Lorentz[7], int T)
{
	long double chi[5];
	if(T == 0)
	{
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
			OutputFile << "," << Medium_Euclidean[i-1][1].first/Vacuum_Euclidean[i-1][1].first-Medium_Euclidean[i-1][0].first/Vacuum_Euclidean[i-1][0].first << "±" << sqrt(pow(Medium_Euclidean[i-1][1].second/Vacuum_Euclidean[i-1][1].first,2)+pow(Medium_Euclidean[i-1][0].second/Vacuum_Euclidean[i-1][0].first,2)+pow(Medium_Euclidean[i-1][1].first*Vacuum_Euclidean[i-1][1].second/pow(Vacuum_Euclidean[i-1][0].first,2),2)+pow(Medium_Euclidean[i-1][0].first*Vacuum_Euclidean[i-1][0].second/pow(Vacuum_Euclidean[i-1][0].first,2),2)) << "," << flush;
			for(int j = 0; j < 7; j++)
			{
				OutputFile << Medium_Spatial[i-1][j].first/Vacuum_Spatial[j].first << "±" << sqrt(pow(Medium_Spatial[i-1][j].second/Vacuum_Spatial[j].first,2)+pow(Medium_Spatial[i-1][j].first*Vacuum_Spatial[j].second/pow(Vacuum_Spatial[j].first,2),2)) << "," << flush;
			}
			OutputFile << chi[i] << "," << flush;
		}

		OutputFile << chi[0] << endl;
	}
	else
	{
		chi[T] = Chi_Square(JPsi, Psi_Prime, Non, Medium_Euclidean, Medium_Spatial, T);

		OutputFile << T << "," << flush;
		JPsi[T]->Print(OutputFile);
		OutputFile << ",";
		Non[T]->Print(OutputFile);
		OutputFile << "," << Medium_Euclidean[T-1][1].first/Vacuum_Euclidean[T-1][1].first-Medium_Euclidean[T-1][0].first/Vacuum_Euclidean[T-1][0].first << "±" << sqrt(pow(Medium_Euclidean[T-1][1].second/Vacuum_Euclidean[T-1][1].first,2)+pow(Medium_Euclidean[T-1][0].second/Vacuum_Euclidean[T-1][0].first,2)+pow(Medium_Euclidean[T-1][1].first*Vacuum_Euclidean[T-1][1].second/pow(Vacuum_Euclidean[T-1][0].first,2),2)+pow(Medium_Euclidean[T-1][0].first*Vacuum_Euclidean[T-1][0].second/pow(Vacuum_Euclidean[T-1][0].first,2),2)) << "," << flush;
		for(int j = 0; j < 7; j++)
		{
			OutputFile << Medium_Spatial[T-1][j].first/Vacuum_Spatial[j].first << "±" << sqrt(pow(Medium_Spatial[T-1][j].second/Vacuum_Spatial[j].first,2)+pow(Medium_Spatial[T-1][j].first*Vacuum_Spatial[j].second/pow(Vacuum_Spatial[j].first,2),2)) << "," << flush;
		}
		for(int j = 0; j < 7; j++)
		{
			OutputFile << Medium_Lorentz[j].first/Vacuum_Spatial[j].first << "±" << sqrt(pow(Medium_Lorentz[j].second/Vacuum_Spatial[j].first,2)+pow(Medium_Lorentz[j].first*Vacuum_Spatial[j].second/pow(Vacuum_Spatial[j].first,2),2)) << "," << flush;
		}
		OutputFile << chi[T] << endl;
	}
	return(chi[0]);
}

long double Chi_Square(Spectral_Inter* JPsi[5], Spectral_Inter* Psi_Prime[5], Spectral_Non* Non[5], pair<long double, long double> Medium_Euclidean[4][2], pair<long double, long double> Medium_Spatial[4][7], int Temp)
{
	long double answer = 0;

	if(Temp != 0)	//To report out each temprature for goodness of fit considerations
	{
		answer += pow(Medium_Euclidean[Temp-1][1].first/Vacuum_Euclidean[Temp-1][1].first-Medium_Euclidean[Temp-1][0].first/Vacuum_Euclidean[Temp-1][0].first-.2,2)/.2;
		for(int j = 0; j < 7; j++)
		{
			answer += pow(Medium_Spatial[Temp-1][j].first/Vacuum_Spatial[j].first-Spatial_Ratio[Temp-1][j],2)/Spatial_Ratio[Temp-1][j];
			if(j < 6)answer += pow((Medium_Spatial[Temp-1][j+1].first/Vacuum_Spatial[j+1].first-Medium_Spatial[Temp-1][j].first/Vacuum_Spatial[j].first>0)?(Medium_Spatial[Temp-1][j+1].first/Vacuum_Spatial[j+1].first-Medium_Spatial[Temp-1][j].first/Vacuum_Spatial[j].first):0,.5);
		}
	}
	else	//Total Chi-squared
	{
		for(int i = 0; i < 4; i++)
			answer += pow(Medium_Euclidean[i][1].first/Vacuum_Euclidean[i][1].first-Medium_Euclidean[i][0].first/Vacuum_Euclidean[i][0].first-.2,2)/.2;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 7; j++)
			{
				answer += pow(Medium_Spatial[i][j].first/Vacuum_Spatial[j].first-Spatial_Ratio[i][j],2)/Spatial_Ratio[i][j];
				if(j < 6)answer += pow((Medium_Spatial[i][j+1].first/Vacuum_Spatial[j+1].first-Medium_Spatial[i][j].first/Vacuum_Spatial[j].first>0)?(Medium_Spatial[i][j+1].first/Vacuum_Spatial[j+1].first-Medium_Spatial[i][j].first/Vacuum_Spatial[j].first):0,.5);
			}
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
	JPsi[1]->Replace(1.97/(2.*3.09946)+Deviation_Points[0],0,1);
	JPsi[1]->Replace(Deviation_Points[2],0,2);
	JPsi[1]->Replace(3.09946+Deviation_Points[4],1,1);
	JPsi[1]->Replace(Deviation_Points[6],1,2);
	JPsi[1]->Replace(Deviation_Points[8],2,1);
	JPsi[1]->Replace(Deviation_Points[10],2,2);
	Non[1]->Replace(1.655+Deviation_Points[12],0,1);
	Non[1]->Replace(Deviation_Points[14],0,2);
	Non[1]->Replace(Deviation_Points[16],1,1);
	Non[1]->Replace(Deviation_Points[18],1,2);
	JPsi[2]->Replace(.4024+(71.*Deviation_Points[0]+32.*Deviation_Points[1])/103.,0,1);
	JPsi[2]->Replace((71.*Deviation_Points[2]+32.*Deviation_Points[3])/103.,0,2);
	JPsi[2]->Replace(3.125+(71.*Deviation_Points[4]+32.*Deviation_Points[5])/103.,1,1);
	JPsi[2]->Replace((71.*Deviation_Points[6]+32.*Deviation_Points[7])/103.,1,2);
	JPsi[2]->Replace((71.*Deviation_Points[8]+32.*Deviation_Points[9])/103.,2,1);
	JPsi[2]->Replace((71.*Deviation_Points[10]+32.*Deviation_Points[11])/103.,2,2);
	Non[2]->Replace(1.59+(71.*Deviation_Points[12]+32.*Deviation_Points[13])/103.,0,1);
	Non[2]->Replace((71.*Deviation_Points[14]+32.*Deviation_Points[15])/103.,0,2);
	Non[2]->Replace((71.*Deviation_Points[16]+32.*Deviation_Points[17])/103.,1,1);
	Non[2]->Replace((71.*Deviation_Points[18]+32.*Deviation_Points[19])/103.,1,2);
	JPsi[3]->Replace(.403047+(40.*Deviation_Points[0]+63.*Deviation_Points[1])/103.,0,1);
	JPsi[3]->Replace((40.*Deviation_Points[2]+63.*Deviation_Points[3])/103.,0,2);
	JPsi[3]->Replace(3.151+(40.*Deviation_Points[4]+63.*Deviation_Points[5])/103.,1,1);
	JPsi[3]->Replace((40.*Deviation_Points[6]+63.*Deviation_Points[7])/103.,1,2);
	JPsi[3]->Replace((40.*Deviation_Points[8]+63.*Deviation_Points[9])/103.,2,1);
	JPsi[3]->Replace((40.*Deviation_Points[10]+63.*Deviation_Points[11])/103.,2,2);
	Non[3]->Replace(1.51+(40.*Deviation_Points[12]+63.*Deviation_Points[13])/103.,0,1);
	Non[3]->Replace((40.*Deviation_Points[14]+63.*Deviation_Points[15])/103.,0,2);
	Non[3]->Replace((40.*Deviation_Points[16]+63.*Deviation_Points[17])/103.,1,1);
	Non[3]->Replace((40.*Deviation_Points[18]+63.*Deviation_Points[19])/103.,1,2);
	JPsi[4]->Replace(.266834+Deviation_Points[1],0,1);
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
	long double Delta = (b-a)*pow(chi,.6)/100.;
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
