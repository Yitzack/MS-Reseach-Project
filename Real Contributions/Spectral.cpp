//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<string>
#include<chrono>
#include"Around.h"
#include"Spectral.h"
using namespace std;

int Start_Point(int, char[70]);							//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);			//Checks to see if file header matches input parameters and clears it if not
long double Set_Mq(long double, long double, long double);				//Momentum dependence for the quark mass, <number> 0 causes it to be constant
long double Set_Lambda(long double, long double, long double, long double, int);	//Momentum dependence for the potential cutoff, <number> 0 causes it to be constant
long double Set_C(long double, long double, long double, long double, long double);	//Momentum dependence for the coupling constant, <number> 0 causes it to be constant
void Load_File(char*);									//Load File from disk to RAM for ReG and ReG_Err
long double i_k_wrap(long double, long double[], long double);

int main(int argc, char* argv[])
{
	char FileApp[70];
	char FileReserved[70];
#ifdef BB	//use option -D BB= to activate bottomium macro
	char File[70] = "data/ReSpectralbb";	//Name of the file
#endif
#ifdef CC	//use option -D CC= to activate charmonium macro
	char File[70] = "data/ReSpectralcc";
#endif
	strcpy(FileApp,File);
	strcat(FileApp,".");

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
	strcat(FileApp, "Half.");
#endif

	char Number_c[5];
	string Number_s;
	char* Process = argv[1];
	strcat(File, argv[3]);	//Appends the temprature to the file name
	strcat(FileApp, argv[3]);
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

	int i,j,l;					//Counters
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
	Elements<Around> holder;					//Calculated value before distribution to Table

	strcpy(FileReserved, FileApp);
	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = Start; i <= Finish; i++)
	{
		for(j = iProcess; j < 576; j+=Total)	//Does the subset of j that has been assigned to this process
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
				if(j <= 161)
					Par[4] = pow((j-151.)/100.,2);
				else if(j <= 190)
					Par[4] = pow((j-161.)/10.+.1,2);
				else if(j <= 390)
					Par[4] = pow((j-190.)/100.+3.,2);
				else
					Par[4] = pow((j-390.)/10.+5.,2);
#else
				if(j <= 251)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 451)
					Par[4] = pow((j-251.)/100.+10.,2);
				else
					Par[4] = pow((j-451.)/10.+12.,2);
#endif
			}

			Par[1] = Set_Lambda(atof(argv[5]), Par[3], atof(argv[9]), atof(argv[10]), Temp);
			Par[0] = -Set_C(atof(argv[4]), Par[3], atof(argv[9]), Par[1], atof(argv[10]));
			Par[2] = atof(argv[6]);

			strcpy(FileApp, FileReserved);
			strcat(FileApp, ".");
			Number_s = to_string(i);
			for(l = 0; l < Number_s.length(); l++)
				Number_c[l] = Number_s[l];
			Number_c[l] = '\0';
			strcat(FileApp, Number_c);
			strcat(FileApp, ".");
			Number_s = to_string(j);
			for(l = 0; l < Number_s.length(); l++)
				Number_c[l] = Number_s[l];
			Number_c[l] = '\0';
			strcat(FileApp, Number_c);
			strcat(FileApp, ".csv");
			try
			{
				Load_File(FileApp);
			}
			catch(...)
			{
				continue;
			}

			auto Start_Time = chrono::system_clock::now();
			try
			{
				holder = theta_Int(Par, Temp);
			}
			catch(...)
			{
				continue;
			}
			auto End_Time = chrono::system_clock::now();
			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << holder[0] << " " << holder[1] << " " << holder[2] << " " << holder[3] << " " << chrono::duration_cast<chrono::nanoseconds>(End_Time-Start_Time).count()/1000000000. << endl;
		}
		TPlot << endl;
	}

	return(0);
}

void Load_File(char* File_Name)
{
	int xSize;
	int ySize;
	long double Holder[2];
	char Bin;
	//char Full_File_Name[210] = "/run/user/1000/gvfs/sftp:host=ccomp.tamu.edu/home/rfrgroup/isarver/data/ReSpectralcc.Half.1/ReProp";
	char Full_File_Name[210] = "/home/rfrgroup/isarver/data/ReSpectralcc.Half.1/ReProp";
	strncpy(&Full_File_Name[strlen(Full_File_Name)],&File_Name[17],64<210-strlen(Full_File_Name)?64:210-strlen(Full_File_Name));	//64 is the max amount of string in File_Name, 210-strlen() is the space avalible. Hopefully, the File_name will fit in Full_File_Name. This should only be a question on my desktop when "/run/user/1000/gvfs/sftp:host=ccomp.tamu.edu/home/rfrgroup/isarver/data/ReSpectralccProp.0/" is the base directory.
	ifstream File(Full_File_Name);

	if(!File.is_open() || !File.good())
		throw;

	File >> xSize;
	File >> Bin;
	File >> ySize;

	long double** Control = new long double*[xSize];
	long double** Control_Err = new long double*[xSize];

	for(int i = 0; i < xSize; i++)
	{
		Control[i] = new long double[ySize];
		Control_Err[i] = new long double[ySize];

		for(int j = 0; j < ySize; j++)
		{
			File >> Holder[0] >> Bin >> Holder[1];
			Control[i][j] = Holder[0];
			Control_Err[i][j] = Holder[1];
		}
	}

	File.close();

	ReG = Interpolation<long double>(Control, xSize, ySize);
	ReG_Err = Interpolation<long double>(Control_Err, xSize, ySize);

	return;
}

int Start_Point(int Start, char File[70])	//Go through and find largest starting point in file and return it, causes it to repeat last line
{
	ifstream TPlot(File);
	char Line[400];
	int Test;

	TPlot.getline(Line, 400);
	if(!TPlot.is_open())
		return(Start);

	TPlot.getline(Line, 400);
	while(!TPlot.eof())
	{
		Test = atoi(Line);
		if(Test > Start)
			Start = Test;
		TPlot.getline(Line,400);
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
	long double Mqf = 1.9;
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

#if VERSION == 22
	return(sqrt(pow(0.9753762593631373,2)+pow(G*Temp,2)));
#elif VERSION == 24
	return(sqrt(pow(1.7264468603176357,2)+pow(G*Temp,2)/2));
#elif VERSION == 42
	return(pow(pow(2.1236949508354317,4)+pow(G*Temp,4),.25));
#elif VERSION == Exp
	return(sqrt(pow(2.342371511360313,2)+pow(G*Temp,2)));
#endif
}

long double Set_C(long double f0, long double P, long double P0, long double Lambda, long double fraction)
{
	long double f = (f0*pow(P0,2)+(fraction*(1-f0)+f0)*pow(P,2))/(pow(P0,2)+pow(P,2));

#if VERSION == 22
	return(332.7040863772379*f*pow(0.9753762593631373/Lambda,4));
#elif VERSION == 24
	return(138.48957840171963*f*pow(1.7264468603176357/Lambda,8));
#elif VERSION == 42
	return(72.44153811020188*f*pow(2.1236949508354317/Lambda,8));
#elif VERSION == Exp
	return(95.18401144965306*f);
#endif
}
