//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<omp.h>
#include"Spectral.h"
using namespace std;

void Plot(long double, int);	//Plots the functions from table
long double RandFloat(long double, long double);
bool Poll(long double[][2], int);

char* Process;

int main(int argc, char* argv[])
{
#ifndef BB	//use option -D BB= to activate BB macro
	char File[30] = "Spectralcc.";  //Name of the file
#else
     	char File[30] = "Spectralbb.";  //Name of the file
#endif
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot;
	if(atoi(argv[6]) == 0)	//If starting from the beginning, overwrite
		TPlot.open(File);
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[616][3];//*/
	long double Par[5] = {-158.90117114622294, 2.643945190802571, 1.8, 0, 0};
	Elements holder;
	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = atoi(argv[6]); i <= 788; i++)	//Argv[6] allows to restart where ever
	{
		#pragma omp parallel for
		for(j = iProcess+151; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
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
				cout << "If you come down this way, you have to alter the curvilinar system used to sample. You will want to use lines of constant s, or sqrt(E^2+P^2). You must do this before you do anything in s<0GeV^2 region" << endl;
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-j*j;
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-j*j;
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

		for(j = iProcess+151; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-j*j;
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-j*j;
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

	/*cout << setprecision(18);
	//#pragma omp parallel for
	for(int i = 0; i <= 751; i+=27)
	{
		long double Par[5] = {-158.90117114622294, 2.643945190802571, 1.8, 0, 13.69};
		Par[3] = i*0.8;
		holder[0] = theta_Int(Par, 0);
		//#pragma omp critical
		{
			cout << Par[3] << " " << Par[4] << " " << holder[0].store(0) << " " << holder[0].store(1) << " " << holder[0].store(2) << endl;
		}
	}//*/

	/*cerr << setprecision(18);
	long double Previous[] = {.3, .5, .9, 3.5, 9.6, 11.8, 16, .2, 1.5, 2.5, 3, 4, 5.5, 7.7, 0.3, 0.08};
	//long double slist[] = {.01, 3.24, 4., 12.96, 25., 100., 552.25};
	//long double Plist[] = {0, 20, 200, 600};
	//long double slist[] = {4., 25.};
	//long double Plist[] = {0, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600};
	long double slist[] = {552.25};
	long double Plist[] = {0, 21.6, 43.2, 64.8, 86.4, 108, 129.6, 151.2, 172.8, 194.4, 216, 237.6, 259.2, 280.8, 302.4, 324, 345.6, 367.2, 388.8, 410.4, 432, 453.6, 475.2, 496.8, 518.4, 540, 561.6, 583.2, 604.8};
	int count;
	int i, j = 0;
	int s_size = 1;
	int P_size = 28;
	long double error[3*s_size*P_size][2];

	for(i = 0; i < 16; i++)
		Boundary[i] = Previous[i];

	for(i = 0; i < P_size; i++)
	{
		for(j = 0; j < s_size; j++)
		{
			Par[3] = Plist[i];
			Par[4] = slist[j];
			holder[i+P_size*j] = theta_Int(Par, 0);
			cout << Par[3] << " " << sqrt(Par[4]) << " " << holder[i+P_size*j].store(0) << " " << holder[i+P_size*j].store(1) << " " << holder[i+P_size*j].store(2) << endl;
		}
	}

	cout << setprecision(4);
	for(i = 0; i < P_size*s_size; i++)
	{
		error[3*i][0] = abs(holder[i].store(0)/holder[int(floor(float(i)/float(P_size)))*P_size].store(0)-1.);
		error[3*i+1][0] = abs(holder[i].store(1)/holder[int(floor(float(i)/float(P_size)))*P_size].store(1)-1.);
		error[3*i+2][0] = abs(holder[i].store(2)/holder[int(floor(float(i)/float(P_size)))*P_size].store(2)-1.);
		cout << error[3*i][0] << " " << error[3*i+1][0] << " " << error[3*i+2][0] << " " << flush;
	}
	cout << setprecision(18);
	for(i = 0; i < 16; i++)
		cout << Previous[i] << " " << flush;
	cout << endl;
	srand(time(NULL)+atoi(argv[1])*30+100*atoi(argv[2]));

	for(int l = 0; l <= 2000; l++)
	{
		count = rand()%10;
		switch(count)
		{
			case 0:
				Boundary[count] = RandFloat(0,Boundary[1]);
				break;
			case 4:
				Boundary[count] = RandFloat(Boundary[3],Boundary[4]+1.);
				break;
			case 5:
				Boundary[count] = RandFloat(0,Boundary[6]);
				break;
			case 9:
				Boundary[count] = RandFloat(Boundary[8],Boundary[9]+1.);
				break;
			case 10:
			case 11:
				Boundary[count] = RandFloat(0,1);
				break;
			default:
				Boundary[count] = RandFloat(Boundary[count-1],Boundary[count+1]);
				break;
		}

		cout << setprecision(6);
		for(i = 0; i < P_size; i++)
		{
			for(j = 0; j < s_size; j++)
			{
				Par[3] = Plist[i];
				Par[4] = slist[j];
				holder[i+P_size*j] = theta_Int(Par, 0);
				if(l%10 == 9)
					cout << Par[3] << " " << sqrt(Par[4]) << " " << holder[i+P_size*j].store(0) << " " << holder[i+P_size*j].store(1) << " " << holder[i+P_size*j].store(2) << endl;
			}
		}

		cout << setprecision(4);
		for(i = 0; i < P_size*s_size; i++)
		{
			error[3*i][1] = abs(holder[i].store(0)/holder[int(floor(float(i)/float(P_size)))*P_size].store(0)-1.);
			error[3*i+1][1] = abs(holder[i].store(1)/holder[int(floor(float(i)/float(P_size)))*P_size].store(1)-1.);
			error[3*i+2][1] = abs(holder[i].store(2)/holder[int(floor(float(i)/float(P_size)))*P_size].store(2)-1.);
			cout << error[3*i][1] << " " << error[3*i+1][1] << " " << error[3*i+2][1] << " " << flush;
		}
		cout << setprecision(18);
		for(i = 0; i < 12; i++)
			cout << Boundary[i] << " " << flush;
		cout << endl;

		if(Poll(error, P_size*s_size))
		{	//reject
			cout << l << "th attempt rejected" << endl;
			Boundary[count] = Previous[count];
		}
		else	//keep
		{
			cout << l << "th attempt accepted" << endl;
			Previous[count] = Boundary[count];
			for(i = 0; i < 3*s_size*P_size; i++)
				error[0][i] = error[1][i];
		}
	}//*/

	return(0);
}

bool Poll(long double error[][2], int elements)
{
	long double count = 0;

	for(int i = 0; i < elements*3; i++)
		if(error[i][0] != 0)//error[i][0] < error[i][1])	//Count reject conditions
			count += error[i][1]/error[i][0]-1.;

	long double rand = RandFloat(0,4);
	cout << rand << " " << count << endl;
	if(rand < count)
		return(true);
	return(false);
}

long double RandFloat(long double a, long double b)
{
	return(a+(b-a)*((long double)rand())/((long double)RAND_MAX));
}
