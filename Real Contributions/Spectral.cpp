//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<chrono>
#include"Spectral.h"
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double Set_Mq(long double, long double, long double);			//Momentum dependence for the quark mass, <number> 0 causes it to be constant
long double Set_Lambda(long double, long double, long double, long double, int);//Momentum dependence for the potential cutoff, <number> 0 causes it to be constant
long double Set_C(long double, long double, long double, long double, long double);//Momentum dependence for the coupling constant, <number> 0 causes it to be constant

int main(int argc, char* argv[])
{
#ifdef BB	//use option -D BB= to activate bottomium macro
	char File[70] = "data/ReSpectralbb";  //Name of the file
#endif
#ifdef CC	//use option -D CC= to activate charmonium macro
	char File[70] = "data/ReSpectralcc";
#endif

#if VERSION == EXP	//use option -D VERSION={Exp,22,24,42} to select one of the potentials
	strcat(File,"Exp.");
#elif VERSION == 22
	strcat(File,"22.");
#elif VERSION == 24
	strcat(File,"24.");
#elif VERSION == 42
	strcat(File,"42.");
#endif

#ifdef HALF	//use option -D HALF= to divide self-energy in half
	strcat(File, "Half.");
#endif

	char* Process = argv[1];
	strcat(File, argv[3]);	//Appends the temprature to the file name
	strcat(File, ".");
	strcat(File, Process);	//Appends the process number to the file name

	bool Restart = Restart_Check(File, argv[4], argv[5], argv[6], argv[9], argv[10]);	//True if scrapping the file contents and restarting, only if all args are already in file header

	ofstream TPlot;
	if(Restart)	//If starting from the beginning, overwrite
	{
		TPlot.open(File);
		TPlot << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[9] << " " << argv[10] << endl;
	}
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);

	int i,j;					//Counters
	int Start, Finish;
	Start = atoi(argv[7]);				//Initial starting point
	Finish = atoi(argv[8]);			//Finish at the point given
	Start = Start_Point(Start, File);		//Go find last point written to file and maybe start from there
	if(Finish < Start && Finish >= atoi(argv[7]))	//Missing point directive issued, go back and get the missed points
		Start = atoi(argv[7]);

	const int iProcess = atoi(argv[1]) % atoi(argv[2]);	//Assigned column(s)
	const int Total = atoi(argv[2]);			//Number of concurent threads
	const int Temp = atoi(argv[3]);			//Temprature enumeration
	long double Par[5];					//Parameters to be used in calculation {Coupling constant, potential cutoff, quark mass, P, s}
	Elements holder;					//Calculated value before distribution to Table

	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration, usefull for doing integrals to infinite s. Not used as designed downstream.

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = Start; i <= Finish; i++)
	{
		for(j = iProcess; j < 616; j+=Total)	//Does the subset of j that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					Par[3] = i/10.+j/10.;
					Par[4] = -i*j/50.-j*j/100.;
				}
				else
				{
					Par[3] = i+j/10.-187.2;
					Par[4] = -j/5.*(i-187.2)-j*j/100.;
				}
			}
			else
			{
				Par[3] = i*.8;
				if(i < 0)
					Par[3] = ((long double)(i%7)/8.-.125-floor((long double)(i)/7.))*.8;
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

			Par[1] = Set_Lambda(atof(argv[5]), Par[3], atof(argv[9]), atof(argv[10]), Temp);
			Par[0] = -Set_C(atof(argv[4]), Par[3], atof(argv[9]), Par[1], atof(argv[10]));
			Par[2] = atof(argv[6]);

			auto Start_Time = chrono::system_clock::now();
			holder = theta_Int(Par, Temp);
			auto End_Time = chrono::system_clock::now();
			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << holder[0] << " " << holder[1] << " " << holder[2] << " " << holder[3] << " " << chrono::duration_cast<chrono::nanoseconds>(End_Time-Start_Time).count()/1000000000. << endl;
		}
		TPlot << endl;
	}

	return(0);
}

int Start_Point(int Start, char File[70])	//Go through and find largest starting point in file and return it, causes it to repeat last line
{
	ifstream TPlot(File);
	char Line[200];
	int Test;

	TPlot.getline(Line, 200);
	if(!TPlot.is_open())
		return(Start);

	TPlot.getline(Line, 200);
	while(!TPlot.eof())
	{
		Test = atoi(Line);
		if(Test > Start)
			Start = Test;
		TPlot.getline(Line,200);
	}

	return(Start);
}

bool Restart_Check(char File[70], char* g, char* Lambda, char* Mq, char* P0, char* fraction)	//Looks for same input values in fist line of file
{
	ifstream InFile(File);

	if(InFile.is_open() == false)
		return(true);

	double g_File;
	double Lambda_File;
	double Mq_File;
	double P0_File;
	double fraction_File;
	InFile >> g_File;
	InFile >> Lambda_File;
	InFile >> Mq_File;
	InFile >> P0_File;
	InFile >> fraction_File;
	InFile.close();

	if(abs(g_File/atof(g)-1.) < .0001 &&
	   abs(Mq_File/atof(Mq)-1.) < .0001 &&
	   (abs(Lambda_File/atof(Lambda)-1.) < .0001 || Lambda_File-atof(Lambda) < .0001) &&
	   abs(P0_File/atof(P0)-1.) < .0001 &&
	   (abs(fraction_File/atof(fraction)-1.) < .0001 || fraction_File-atof(fraction) < .0001))
		return(false);

	InFile.close();
	return(true);
}

long double Set_Mq(long double Mq0, long double P, long double P0)
{
#ifndef BB
	long double Mqf = 1.8;
#else
	long double Mqf = 5.25;
#endif

	return((Mq0*pow(P0,2)+Mqf*pow(P,2))/(pow(P0,2)+pow(P,2)));
}

long double Set_Lambda(long double G0, long double P, long double P0, long double fraction, int T)
{
	long double G = G0*(pow(P0,2)+(1-fraction)*pow(P,2))/(pow(P0,2)+pow(P,2));
	//long double G = G0;
	long double TempList[] = {0,.194,.258,.32,.4};
	long double Temp = TempList[T];

#if VERSION == Exp
	return(sqrt(pow(2.979,2)+pow(G*Temp,2)));
#elif VERSION == 22
	return(sqrt(pow(1.47132,2)+pow(G*Temp,2)));
#elif VERSION == 24
	return(sqrt(pow(2.29444,2)+pow(G*Temp,2)/2));
#elif VERSION == 42
	return(pow(pow(2.7,4)+pow(G*Temp,4),.25));
#endif
}

long double Set_C(long double f0, long double P, long double P0, long double Lambda, long double fraction)
{
	long double f = (f0*pow(P0,2)+(fraction*(1-f0)+f0)*pow(P,2))/(pow(P0,2)+pow(P,2));

#if VERSION == Exp
	return(50.3627*f);
#elif VERSION == 22
	return(116.253*f*pow(1.47132/Lambda,4));
#elif VERSION == 24
	return(65.8549*f*pow(2.29444/Lambda,8));
#elif VERSION == 42
	return(38.4541*f*pow(2.7/Lambda,8));
#endif
}
