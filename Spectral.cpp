//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<omp.h>
#include"Spectral.h"
using namespace std;

void Plot(long double, int);	//Plots the functions from table

char* Process;

int main(int argc, char* argv[])
{
#ifndef BB	//use option -D BB= to activate BB macro
	char File[25] = "Spectralcc.";  //Name of the file
#else
     	char File[25] = "Spectralbb.";  //Name of the file
#endif
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot;
	/*if(atoi(argv[6]) == 0)	//If starting from the beginning, overwrite
		TPlot.open(File);
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);//*/
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[616][3];
	//long double Par[5] = {-158.90117114622294, 2.643945190802571, 1.8, 0, 0};
	Elements holder;
	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	/*TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = atoi(argv[6]); i <= 788; i++)	//Argv[6] allows to restart where ever
	{
		#pragma omp parallel for
		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			Par[1] = 2.643945190802571;
			Par[2] = 1.8;
			if(argc == 4)
				switch(Temp)
				{
					case 1:
						Par[1] *= exp(-4./35.);
						Par[2] = 1.848;//1.902;
						break;
					case 2:
						Par[1] *= exp(-2./7.);
						Par[2] = 1.719;//1.777;
						break;
					case 3:
						Par[1] *= exp(-4./7.);
						Par[2] = 1.563;//1.652;
						break;
				}
			else
			{
				Par[1] *= atof(argv[4]);
				Par[2] = atof(argv[5]);
			}

			if(j <= 150)
			{
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-pow(j,2);
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-pow(j,2);
				}
			}
			else
			{
				Par[3] = i*.8;
#ifndef BB
                                if(j <= 181)
                                        Par[4] = pow((j-151.)/10.,2);
                                else if(j <= 381)
                                        Par[4] = pow((j-181.)/100.+3.,2);
                                else if(j <= 566)
                                        Par[4] = pow((j-381.)/10.+5.,2);
#else
                                if(j <= 251)
                                        Par[4] = pow((j-151.)/10.,2);
                                else if(j <= 451)
                                        Par[4] = pow((j-251.)/100.+10.,2);
                                else if(j <= 566)
                                        Par[4] = pow((j-451.)/10.+12.,2);
#endif
				else
					Par[4] = 552.25+GaussLa[j-567];
			}

			if(j > 150 && i > 751)
			{
				Table[j][0] = Table[j][1] = Table[j][2] = 0;
			}
			else
			{
				holder = theta_Int(Par, Temp);
				Table[j][0] = holder.store(0);
				Table[j][1] = holder.store(1);
				Table[j][2] = holder.store(2);
			}
		}

		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-pow(j,2);
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-pow(j,2);
				}
			}
			else
			{
				Par[3] = i*.8;
#ifndef BB
                                if(j <= 181)
                                        Par[4] = pow((j-151.)/10.,2);
                                else if(j <= 381)
                                        Par[4] = pow((j-181.)/100.+3.,2);
                                else if(j <= 566)
                                        Par[4] = pow((j-381.)/10.+5.,2);
#else
                                if(j <= 251)
                                        Par[4] = pow((j-151.)/10.,2);
                                else if(j <= 451)
                                        Par[4] = pow((j-251.)/100.+10.,2);
                                else if(j <= 566)
                                        Par[4] = pow((j-451.)/10.+12.,2);
#endif
				else
					Par[4] = 552.25+GaussLa[j-567];
			}

			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		}
		TPlot << endl;
	}

	//TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	long double Par[5] = {-158.90117114622294, 2.643945190802571, 1.8, 0, 114.49};
	Par[3] = 0;	//P
	Par[4] = 20;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 20;	//P
	Par[4] = 20;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 50;	//P
	Par[4] = 20;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 100;	//P
	Par[4] = 20;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 500;	//P
	Par[4] = 20;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 0;	//P
	Par[4] = 14.5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 20;	//P
	Par[4] = 14.5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 50;	//P
	Par[4] = 14.5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 100;	//P
	Par[4] = 14.5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 500;	//P
	Par[4] = 14.5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 0;	//P
	Par[4] = 12;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 20;	//P
	Par[4] = 12;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 50;	//P
	Par[4] = 12;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 100;	//P
	Par[4] = 12;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 500;	//P
	Par[4] = 12;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 0;	//P
	Par[4] = 5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 20;	//P
	Par[4] = 5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 50;	//P
	Par[4] = 5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 100;	//P
	Par[4] = 5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 500;	//P
	Par[4] = 5;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 10;	//P
	Par[4] = -100;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 10;	//P
	Par[4] = -90;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 10;	//P
	Par[4] = -80;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 10;	//P
	Par[4] = -70;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	Par[3] = 10;	//P
	Par[4] = -60;	//s
	holder = theta_Int(Par, 1);
	cout << Par[3] << " " << Par[4] << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;

	/*Par[3] = 100;
	Par[4] = 14.5;
	long double theta = 1.54153347495665988;
	long double k[] = {0.548445522450565892, 0.559228397420848182, 0.565313950153802751, 0.56685156427547692, 0.567286246613545104, 0.568056451968565267, 0.569142408804113881, 0.570516013273339712, 0.572141695054432615, 0.573977351658789836, 0.575975441317856306, 0.578084214951222266, 0.580249056667365309, 0.582413898383508351, 0.584522672016874311, 0.586520761675940781, 0.588356418280298002, 0.589982100061390905, 0.591355704530616737, 0.59244166136616535, 0.593211866721185514, 0.593646549059253697, 0.59385156427547692, 0.594286246613545103, 0.595056451968565267, 0.59614240880411388, 0.597516013273339712, 0.599141695054432615, 0.600977351658789836, 0.602975441317856306, 0.605084214951222266, 0.607249056667365308, 0.609413898383508351, 0.611522672016874311, 0.613520761675940781, 0.615356418280298002, 0.616982100061390905, 0.618355704530616736, 0.61944166136616535, 0.620211866721185513, 0.620646549059253697, 0.620749056667365308, 0.620851564275476919, 0.621286246613545103, 0.622056451968565267, 0.62314240880411388, 0.624516013273339711, 0.626141695054432615, 0.627977351658789835, 0.629975441317856305, 0.632084214951222266, 0.634249056667365308, 0.63641389838350835, 0.638522672016874311, 0.640520761675940781, 0.642356418280298001, 0.643982100061390905, 0.645355704530616736, 0.64644166136616535, 0.647211866721185513, 0.647646549059253697, 0.647851564275476919, 0.648286246613545103, 0.649056451968565266, 0.65014240880411388, 0.651516013273339711, 0.653141695054432614, 0.654977351658789835, 0.656975441317856305, 0.659084214951222265, 0.661249056667365308, 0.66341389838350835, 0.66552267201687431, 0.66752076167594078, 0.669356418280298001, 0.670982100061390904, 0.672355704530616736, 0.673441661366165349, 0.674211866721185513, 0.674646549059253696, 0.676184163180927865, 0.682269715913882434, 0.693052590884164724, 0.708255986581845314, 0.727486449151006953, 0.7502459940863076, 0.775945186547308689, 0.803918441774239268, 0.833441272641362711, 0.863749056667365305, 0.894056840693367899, 0.923579671560491342, 0.951552926787421921, 0.977252119248423011, 1.00001166418372366};
	theta_Int(Par, 1);
	//cout << setprecision(18);

	//for(int i = 0; i < 95; i++)
	//	for(long double omega = 50.05; omega <= 50.072; omega += .00001)
	//		cout << k[i] << " " << omega << " " << Spin_Sum(Par, omega, k[i], theta)*Folding_Integrand(Par,omega,k[i],theta,1) << endl;*/

	return(0);
}
