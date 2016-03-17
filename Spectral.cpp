//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<omp.h>
#include"Spectral.h"
using namespace std;

//#define DELTAE .00001	//This macro allows looking at a few specific points around the bound state
#undef DELTAE	//This undef macro allows the computation of the normal table for Fourier.cpp

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
#ifdef DELTAE
	const long double E[] = {3.040308-40*DELTAE, 3.040308-39*DELTAE, 3.040308-38*DELTAE,  3.040308-37*DELTAE, 3.040308-36*DELTAE, 3.040308-35*DELTAE,  3.040308-34*DELTAE, 3.040308-33*DELTAE, 3.040308-32*DELTAE,  3.040308-31*DELTAE, 3.040308-30*DELTAE, 3.040308-29*DELTAE,  3.040308-28*DELTAE, 3.040308-27*DELTAE, 3.040308-26*DELTAE,  3.040308-25*DELTAE, 3.040308-24*DELTAE, 3.040308-23*DELTAE,  3.040308-22*DELTAE, 3.040308-21*DELTAE, 3.040308-20*DELTAE,  3.040308-19*DELTAE, 3.040308-18*DELTAE, 3.040308-17*DELTAE,  3.040308-16*DELTAE, 3.040308-15*DELTAE, 3.040308-14*DELTAE,  3.040308-13*DELTAE, 3.040308-12*DELTAE, 3.040308-11*DELTAE,  3.040308-10*DELTAE, 3.040308-9*DELTAE, 3.040308-8*DELTAE,  3.040308-7*DELTAE, 3.040308-6*DELTAE, 3.040308-5*DELTAE,  3.040308-4*DELTAE, 3.040308-3*DELTAE, 3.040308-2*DELTAE,  3.040308-DELTAE, 3.040308, 3.040308+DELTAE, 3.040308+2*DELTAE,  3.040308+3*DELTAE, 3.040308+4*DELTAE, 3.040308+5*DELTAE,  3.040308+6*DELTAE, 3.040308+7*DELTAE, 3.040308+8*DELTAE,  3.040308+9*DELTAE, 3.040308+10*DELTAE, 3.040308+11*DELTAE,  3.040308+12*DELTAE, 3.040308+13*DELTAE, 3.040308+14*DELTAE,  3.040308+15*DELTAE, 3.040308+16*DELTAE, 3.040308+17*DELTAE,  3.040308+18*DELTAE, 3.040308+19*DELTAE, 3.040308+20*DELTAE,  3.040308+21*DELTAE, 3.040308+22*DELTAE, 3.040308+23*DELTAE,  3.040308+24*DELTAE, 3.040308+25*DELTAE, 3.040308+26*DELTAE,  3.040308+27*DELTAE, 3.040308+28*DELTAE, 3.040308+29*DELTAE,  3.040308+30*DELTAE, 3.040308+31*DELTAE, 3.040308+32*DELTAE,  3.040308+33*DELTAE, 3.040308+34*DELTAE, 3.040308+35*DELTAE,  3.040308+36*DELTAE, 3.040308+37*DELTAE, 3.040308+38*DELTAE,  3.040308+39*DELTAE, 3.040308+40*DELTAE};
#else
	const long double E[] = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.540308, 2.545308, 2.550308, 2.555308, 2.560308, 2.565308, 2.570308, 2.575308, 2.580308, 2.585308, 2.590308, 2.595308, 2.600308, 2.605308, 2.610308, 2.615308, 2.620308, 2.625308, 2.630308, 2.635308, 2.640308, 2.645308, 2.650308, 2.655308, 2.660308, 2.665308, 2.670308, 2.675308, 2.680308, 2.685308, 2.690308, 2.695308, 2.700308, 2.705308, 2.710308, 2.715308, 2.720308, 2.725308, 2.730308, 2.735308, 2.740308, 2.745308, 2.750308, 2.755308, 2.760308, 2.765308, 2.770308, 2.775308, 2.780308, 2.785308, 2.790308, 2.795308, 2.800308, 2.805308, 2.810308, 2.815308, 2.820308, 2.825308, 2.830308, 2.835308, 2.840308, 2.845308, 2.850308, 2.855308, 2.860308, 2.865308, 2.870308, 2.875308, 2.880308, 2.885308, 2.890308, 2.895308, 2.900308, 2.905308, 2.910308, 2.915308, 2.920308, 2.925308, 2.930308, 2.935308, 2.940308, 2.945308, 2.950308, 2.955308, 2.960308, 2.965308, 2.970308, 2.975308, 2.980308, 2.985308, 2.990308, 2.995308, 3.000308, 3.005308, 3.010308, 3.015308, 3.020308, 3.025308, 3.030308, 3.035308, 3.040308, 3.045308, 3.050308, 3.055308, 3.060308, 3.065308, 3.070308, 3.075308, 3.080308, 3.085308, 3.090308, 3.095308, 3.100308, 3.105308, 3.110308, 3.115308, 3.120308, 3.125308, 3.130308, 3.135308, 3.140308, 3.145308, 3.150308, 3.155308, 3.160308, 3.165308, 3.170308, 3.175308, 3.180308, 3.185308, 3.190308, 3.195308, 3.200308, 3.205308, 3.210308, 3.215308, 3.220308, 3.225308, 3.230308, 3.235308, 3.240308, 3.245308, 3.250308, 3.255308, 3.260308, 3.265308, 3.270308, 3.275308, 3.280308, 3.285308, 3.290308, 3.295308, 3.300308, 3.305308, 3.310308, 3.315308, 3.320308, 3.325308, 3.330308, 3.335308, 3.340308, 3.345308, 3.350308, 3.355308, 3.360308, 3.365308, 3.370308, 3.375308, 3.380308, 3.385308, 3.390308, 3.395308, 3.400308, 3.405308, 3.410308, 3.415308, 3.420308, 3.425308, 3.430308, 3.435308, 3.440308, 3.445308, 3.450308, 3.455308, 3.460308, 3.465308, 3.470308, 3.475308, 3.480308, 3.485308, 3.490308, 3.495308, 3.500308, 3.505308, 3.510308, 3.515308, 3.520308, 3.525308, 3.530308, 3.535308, 3.540308, 3.55, 3.56375, 3.5775, 3.59125, 3.605, 3.61875, 3.6325, 3.64625, 3.66, 3.67375, 3.6875,  3.70125, 3.715, 3.72875, 3.7425, 3.75625, 3.77, 3.78375, 3.7975, 3.81125, 3.825, 3.83875, 3.8525, 3.86625, 3.88, 3.89375, 3.9075, 3.92125, 3.935, 3.94875, 3.9625, 3.97625, 3.99, 4.00375, 4.0175, 4.03125, 4.045, 4.05875, 4.0725, 4.08625, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21., 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22., 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23., 23.1, 23.2, 23.3, 23.4, 23.5};
#endif
	int Frame = 0;	//Frame count
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	complex<long double> TMat;
	long double Table[462][3];
	long double Par[6] = {-127.995280691106, 1.4049344847006076, 1.8, 0, 0, .032};

	if(Temp != 0)
		Par[5] = 0;

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.

	for(i = 0; i <= 751; i++)
	{
		#pragma omp parallel for
#ifdef DELTAE
		for(j = 81*iProcess/(Total); j < 81*(iProcess+1)/Total; j++)	//Does the subset of E that has been assigned to this process
#else
		for(j = iProcess; j < 462; j+=Total)	//Does the subset of E that has been assigned to this process
#endif
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
			Par[4] = pow(E[j],2);
			TMat = TMatrix(Par, Temp);
			if(Par[4] < pow(2.*Par[2]))
				TMat *= complex<long double>(Par[0]);
			else
				TMat *= complex<long double>(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+Par[4]-pow(2*Par[2],2)),2));
			Table[j][1] = TMat.real();
			Table[j][2] = TMat.imag();
			Table[j][0] = Spectral(Par, Temp);
			//cout << E[j] << " " << Table[j][1] << " " << Table[j][0] << endl;
		}

#ifdef DELTAE
		for(j = 81*iProcess/(Total); j < 81*(iProcess+1)/Total; j++)	//Does the subset of E that has been assigned to this process
#else
		for(j = iProcess; j < 462; j+=Total)	//Does the subset of E that has been assigned to this process
#endif
			TPlot << Temp <<  " " << .8*i << " " << E[j] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
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
