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
	ofstream TPlot(File);
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[561][3];
	long double Par[5] = {-42.96210630522018, 2.1348192815218754, 1.8, 0, 0};
	Elements holder;
	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	#pragma omp parallel for
	for(i = 0; i <= 788; i++)
	{
		for(j = iProcess+151; j <= 561; j+=Total)	//Does the subset of E that has been assigned to this process
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
				if(j <= 187)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 327)
					Par[4] = pow((j-187.)/100.+3.6,2);
				else if(j <= 512)
					Par[4] = pow((j-327.)/10.+5,2);
				else
					Par[4] = 552.25+GaussLa[j-513];
			}

			if(Par[3] > 600.8)
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

		for(j = iProcess; j < 561; j+=Total)	//Does the subset of E that has been assigned to this process
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
				if(j <= 187)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 327)
					Par[4] = pow((j-187.)/100.+3.6,2);
				else if(j <= 512)
					Par[4] = pow((j-327.)/10.+5,2);
				else
					Par[4] = 552.25+GaussLa[j-513];
			}

			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << endl;//Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		}
		TPlot << endl;
	}

	TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	return(0);
}






























