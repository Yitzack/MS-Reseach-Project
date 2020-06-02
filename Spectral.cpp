//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#ifdef _OPENMP
#include<omp.h>
#endif
#include"Spectral.h"
using namespace std;

int Start_Point(int, char[30]);
bool Restart_Check(char[30], char*, char*, char*);
long double Set_G(long double, long double, int);
long double Set_Mq(long double, long double, int);

char* Process;

int main(int argc, char* argv[])
{
#ifdef BB	//use option -D BB= to activate BB macro
	char File[70] = "data/Spectralbb.";  //Name of the file
#endif
#ifdef CC
     	char File[70] = "data/Spectralcc.";  //Name of the file
#endif
#ifdef RIEK
     	char File[30] = "SpectralccRiek.";  //Name of the file
#endif
#ifdef SHUAI
     	char File[30] = "SpectralccShuai.";  //Name of the file
#endif
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name

	bool Restart = Restart_Check(File, argv[4], argv[5], argv[6]);	//True if restarting

	ofstream TPlot;
	if(Restart)	//If starting from the beginning, overwrite
	{
		TPlot.open(File);
		TPlot << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << " " << argv[8] << " " << argv[9] << endl;
	}
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);
	int i,j;	//counters
	int Finish, Start;
	if(argc == 12)
		Finish = atoi(argv[11]);
	else
		Finish = 788;
	if(argc >= 11)
		Start = atoi(argv[10]);
	Start = Start_Point(Start, File);
	if(Finish < Start && Finish >= atoi(argv[10]))	//Go back and get the missed point(s)
		Start = atoi(argv[10]);
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	const int Temp = atoi(argv[3]);
	long double T;
	long double Table[616][4];
	long double Par[8] = {-.5306436016791014, 8.699892671305086, 1.8, 0, 0, 0, 0, 0};
	Elements holder;

	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	switch(Temp)
	{
		case 1:
			T = .194;
			break;
		case 2:
			T = .258;
			break;
		case 3:
			T = .32;
			break;
		case 4:
			T = .4;
			break;
	}

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = Start; i <= Finish; i++)	//Argv[6] allows to restart where ever
	{
		#pragma omp parallel for
		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process, calculation loop
		{
			long double ParPrivate[8];
			if(j <= 150)
			{
				if(i <= 208)
				{
					ParPrivate[3] = i/10.+j/10.;
					ParPrivate[4] = -i*j/50.-j*j/100.;
				}
				else
				{
					ParPrivate[3] = i+j/10.-187.2;
					ParPrivate[4] = -j/5.*(i-187.2)-j*j/100.;
				}
			}
			else
			{
				ParPrivate[3] = i*.8;
#ifndef BB
				if(j <= 181)
					ParPrivate[4] = pow((j-151.)/10.,2);
				else if(j <= 381)
					ParPrivate[4] = pow((j-181.)/100.+3.,2);
				else if(j <= 566)
					ParPrivate[4] = pow((j-381.)/10.+5.,2);
#else
				if(j <= 251)
					ParPrivate[4] = pow((j-151.)/10.,2);
				else if(j <= 451)
					ParPrivate[4] = pow((j-251.)/100.+10.,2);
				else if(j <= 566)
					ParPrivate[4] = pow((j-451.)/10.+12.,2);
#endif
				else
					ParPrivate[4] = 552.25+GaussLa[j-567];
			}

			ParPrivate[0] = -atof(argv[4]);//-.5306436016791014;//*Set_G(atof(argv[4]), ParPrivate[3], Temp);
			ParPrivate[1] = atof(argv[5]);//8.699892671305086;
			ParPrivate[2] = atof(argv[6]);//Set_Mq(atof(argv[6]), ParPrivate[3], Temp);
			ParPrivate[5] = atof(argv[7]);//Set_Mq(atof(argv[6]), ParPrivate[3], Temp);
			ParPrivate[6] = atof(argv[8]);//Set_Mq(atof(argv[6]), ParPrivate[3], Temp);
			ParPrivate[7] = atof(argv[9]);//Set_Mq(atof(argv[6]), ParPrivate[3], Temp);

			if(j > 150 && i > 751)
			{
				Table[j][0] = Table[j][1] = Table[j][2] = Table[j][3] = 0;
			}
			else
			{
				holder = theta_Int(ParPrivate, Temp);
				Table[j][0] = holder.store(0);
				Table[j][1] = holder.store(1);
				Table[j][2] = holder.store(2);
				Table[j][3] = holder.store(3);
			}
		}

		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process, output loop
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

			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << " " << Table[j][3] << endl;
		}
		TPlot << endl;
	}

	/*TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	return(0);
}

int Start_Point(int Start, char File[30])
{
	ifstream TPlot(File);
	char Line[200];
	int Test;

	TPlot.getline(Line, 200);
	if(!TPlot.is_open())
		return(Start);

	do
	{
		Test = atoi(Line);
		if(Test > Start)
			Start = Test;
		TPlot.getline(Line,200);
	}while(!TPlot.eof());

	return(Start);
}

bool Restart_Check(char File[30], char* g, char* Lambda, char* Mq)
{
	ifstream InFile(File);

	if(InFile.is_open() == false)
		return(true);

	double g_File;
	double Lambda_File;
	double Mq_File;
	InFile >> g_File;
	InFile >> Lambda_File;
	InFile >> Mq_File;
	InFile.close();

	if(abs(g_File/atof(g)-1.) < .0001 && abs(Mq_File/atof(Mq)-1.) < .0001 && abs(Lambda_File/atof(Lambda)-1.) < .0001)
		return(false);

	InFile.close();
	return(true);
}

long double Set_Mq(long double Mq0, long double P, int Temp)
{
	long double Lambda = 8.699892671305086;
	long double T;

	switch(Temp)
	{
		case 0:
			T = 0;
			break;
		case 1:
			T = .194;
			break;
		case 2:
			T = .258;
			break;
		case 3:
			T = .32;
			break;
		case 4:
			T = .4;
			break;
	}

#ifndef BB
	long double Mqf = 1.8;
#else
	long double Mqf = 5.25;
#endif
	long double Delta_Mq = Mqf-Mq0;

	return(Mqf-Delta_Mq/(1+log(1.+pow(P/Lambda,2))));
}

long double Set_G(long double G0, long double P, int Temp)
{
	long double Lambda = 8.699892671305086;
	long double T;

	switch(Temp)
	{
		case 0:
			T = 0;
			break;
		case 1:
			T = .194;
			break;
		case 2:
			T = .258;
			break;
		case 3:
			T = .32;
			break;
		case 4:
			T = .4;
			break;
	}

	long double Gf = 1.;
	long double Delta_G = Gf-G0;

	return(Gf-Delta_G/(1+log(1.+pow(P/Lambda,2))));
}
