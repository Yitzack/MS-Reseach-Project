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
	char File[25] = "Spectral.";	//Name of the file
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	//ofstream TPlot(File);
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[811][3];
	long double Par[5] = {-42.96210630522018, 2.1348192815218754, 1.8, 2, 3};
	Elements holder;
	long double GaussLa[] = {23.5292089494940390418, 23.6539325380822080769, 23.8784519114339929046, 24.203043968841429832, 24.62804449030959115901, 25.15388906539884363591, 25.78111923347644653209, 26.51038628120128830529, 27.34245522739668292116, 28.27820943138205453677, 29.31865597642423461728, 30.46493193346708690195, 31.7183116110416122313, 33.08021491185883249065, 34.5522169380215279328, 36.13605901385725832108, 37.83366132857440339499, 39.64713744153402449126, 41.57881094274913343943, 43.63123462273780157763, 45.8072125823387678126, 48.10982580889231094881, 50.54246186610561423232, 53.10884949880154539486, 55.81309915127963456172, 58.65975065392247902555, 61.65382966748456817771, 64.8009149171740471975, 68.10721884062876818128, 71.5796850753673570501, 75.22610731101421216486, 79.05527556274067844963, 83.0771580886221159235, 87.30313029304261238365, 91.74626653908353044698, 96.42171766800947991981, 101.34720759844820215182, 106.54369909859864667464, 112.03630611197943572002, 117.85557619641319288989, 124.03934816696116679177, 130.63554136224855814149, 137.70653122712858723725, 145.33639878660318539969, 153.64381522449526055617, 162.80719756334274304328, 173.12081975792771442406, 185.14877015704720903095, 200.34630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	#pragma omp parallel for
	for(i = 0; i <= 751; i++)
	{
		for(j = iProcess+400; j < 811; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			Par[1] = 2.1348192815218754;
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

			Par[3] = i*.8;

			if(j <= 400)	//Defining s
			{
				if(i < 13)
					Par[4] = pow(i*j/500.,2)-pow(.8*i,2);
				else
					Par[4] = pow(.8*i+j*13./500.-10.4,2)-pow(.8*i,2);
			}
			else if(j <= 436)
				Par[4] = pow((j-400.)/10.,2);
			else if(j <= 576)
				Par[4] = pow(3.6+(j-436.)/100.,2);
			else if(j <= 761)
				Par[4] = pow(5.+(j-576.)/10.,2);
			else
				Par[4] = pow(GaussLa[j-762],2);

			holder = theta_Int(Par, temp);
			Table[i][0] = holder.store(0);
			Table[i][1] = holder.store(1);
			Table[i][2] = holder.store(2);
		}

		for(j = iProcess+400; j < 811; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 400)	//Defining s
			{
				if(i < 13)
					Par[4] = pow(i*j/500.,2)-pow(.8*i,2);
				else
					Par[4] = pow(.8*i+j*13./500.-10.4,2)-pow(.8*i,2);
			}
			else if(j <= 436)
				Par[4] = pow((j-400.)/10.,2);
			else if(j <= 576)
				Par[4] = pow(3.6+(j-436.)/100.,2);
			else if(j <= 761)
				Par[4] = pow(5.+(j-576.)/10.,2);
			else
				Par[4] = pow(GaussLa[j-762],2);
			TPlot << Temp <<  " " << i << " " << j << " " << Par[4] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		}
		TPlot << endl;
	}

	TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	return(0);
}






























