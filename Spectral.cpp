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
	char File[25] = "DeBugAcc.";	//Name of the file
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot(File);
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[863][3];
	long double Par[5] = {-42.96210630522018, 2.1348192815218754, 1.8, 2, 3};

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	#pragma omp parallel for
	for(i = 0; i <= 0; i++)
	{
		for(j = iProcess+400; j < 863; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			complex<long double> TMat;
			long double Par[5] = {-42.96210630522018, 2.1348192815218754, 1.8, 0, 0};
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
			else if(j <= 425)
				Par[4] = pow((j-400.)/10.,2);
			else if(j <= 626)
				Par[4] = pow(2.540308+(j-426.)/200.,2);
			else if(j <= 667)
				Par[4] = pow(3.55+11.*(j-627.)/800.,2);
			else
				Par[4] = pow(4.1+(j-667.)/10.,2);

			Spectral(Table[j], Par, Temp);
		}

		long double Par[6] = {-127.995280691106, 1.4049344847006076, 1.8, 0, 0, 0};
		for(j = iProcess+400; j < 863; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 400)	//Defining s
			{
				if(i < 13)
					Par[4] = pow(i*j/500.,2)-pow(.8*i,2);
				else
					Par[4] = pow(.8*i+j*13./500.-10.4,2)-pow(.8*i,2);
			}
			else if(j <= 425)
				Par[4] = pow((j-400.)/10.,2);
			else if(j <= 626)
				Par[4] = pow(2.540308+(j-426.)/200.,2);
			else if(j <= 667)
				Par[4] = pow(3.55+11.*(j-627.)/800.,2);
			else
				Par[4] = pow(4.1+(j-667.)/10.,2);
			TPlot << Temp <<  " " << i << " " << j << " " << Par[4] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		}
		TPlot << endl;
	}

	TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	return(0);
}






























