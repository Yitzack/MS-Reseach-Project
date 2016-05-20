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
	int Frame = 0;	//Frame count
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	complex<long double> TMat;
	long double Table[863][3];
	long double Par[6] = {-127.995280691106, 1.4049344847006076, 1.8, 0, 0, 0};

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.

	#pragma omp parallel for
	for(i = 0; i <= 751; i++)
	{
		for(j = iProcess; j < 863; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			Par[1] = 1.4049344847006076;
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

			TMat = TMatrix(Par, Temp);
			if(Par[4] < pow(2.*Par[2],2))
				TMat *= complex<long double>(Par[0]);
			else
				TMat *= complex<long double>(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+Par[4]-pow(2*Par[2],2)),2));
			Table[j][1] = TMat.real();
			Table[j][2] = TMat.imag();
			Table[j][0] = Spectral(Par, Temp);
		}

		for(j = iProcess; j < 863; j+=Total)	//Does the subset of E that has been assigned to this process
			TPlot << Temp <<  " " << i << " " << j << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		TPlot << endl;
	}

	TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	return(0);
}

void Plot(long double Temp, int Frame)
{
	cout << "set terminal pngcairo size 1600,792 enhanced" << endl;	//Output to gnuplot to make the terminal output to a png.
	cout << "set output './Frames/Spectral" << Frame << ".png'" << endl;	//Set the output to a file so specified
	/*switch(Temp)
	{
		case 0:	//Vacuum, T=0
			cout << "set title \"Vacuum\"" << endl;
			break;
		case 3:	//Media, T=1.2T_c
			cout << "set title \"T = 1.2T_c\"" << endl;
			break;
		case 1:	//Media, T=1.5T_c
			cout << "set title \"T = 1.5T_c\"" << endl;
			break;
		case 2:	//Media, T=2T_c
			cout << "set title \"T = 2T_c\"" << endl;
			break;
	}*/
	if(Temp == 0)
		cout << "set title \"Vacuum\"" << endl;
	else
		cout << "set title \"T = " << Temp << " T_c\"" << endl;

	cout << "set multiplot" << endl;

	cout << "set size 1,.5; set origin 0,.5" << endl;
	cout << "set xrange[0:10]" << endl;	//Set the x range to between 0 and 5
	cout << "set xlabel \"E (GeV)\"" << endl;
	cout << "set ylabel \"{/Symbol s}\"" << endl;
	cout << "set label \"m_{J/{/Symbol Y}}^2\" at 3.040308,0 point pointtype 1" << endl;
	cout << "set label \"m_{{/Symbol Y}'}^2\" at 3.662752,0 point pointtype 1" << endl;
	cout << "plot \'Spectral Plot.0\' using 3:4 every :::" << Frame << "::" << Frame << " with lines title \"Spectral Func\"" << endl;//Tell gnuplot to plot the function in the file so named

	cout << "set size 1,.5; set origin 0,0" << endl;
	cout << "set xrange[0:10]" << endl;	//Set the x range to between 0 and 5
	cout << "set xlabel \"E (GeV)\"" << endl;
	cout << "set ylabel \"T\"" << endl;
	cout << "set label \"m_{J/{/Symbol Y}}^2\" at 3.040308,0 point pointtype 1" << endl;
	cout << "set label \"m_{{/Symbol Y}'}^2\" at 3.662752,0 point pointtype 1" << endl;
	cout << "plot \'Spectral Plot.0\' using 3:5 every :::" << Frame << "::" << Frame << " with lines title \"T-Matrix\"" << endl;//Tell gnuplot to plot the function in the file so named

	cout << "unset multiplot" << endl;
	return;
}
