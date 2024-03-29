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

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double ReG12(long double, long double, long double, long double, long double);
long double ImG12(long double, long double, long double, long double, long double);
void Loop_Out1(long double[], int, char[]);
void Loop_Out2(long double[], int, char[]);

int main(int argc, char* argv[])
{
#ifdef BB	//use option -D BB= to activate bottomium macro
	//char File[130] = "/run/user/1000/gvfs/sftp:host=ccomp.tamu.edu/home/rfrgroup/isarver/data/ReSpectralbb.";  //Name of the file
	char File[130] = "data/ReSpectralbb.";  //Name of the file
#endif
#ifdef CC	//use option -D CC= to activate charmonium macro
	//char File[130] = "/run/user/1000/gvfs/sftp:host=ccomp.tamu.edu/home/rfrgroup/isarver/data/ReSpectralcc.Half.1/ReSpectralcc.";
	char File[130] = "data/ReSpectralcc.Half.1/ReSpectralcc.";
	//char File[130] = "data/ReSpectralcc.";
#endif

#ifdef HALF	//use option -D HALF= to divide self-energy in half
	strcat(File, "Half.");
#endif

	char* Process = argv[1];
	char FileApp[130];
	char Number_c[5];
	string Number_s;
	strcat(File, argv[3]);	//Appends the temprature to the file name

	bool Restart = Restart_Check(File, argv[4], argv[5], argv[6], argv[9], argv[10]);	//True if scrapping the file contents and restarting, only if all args are already in file header

	ofstream TPlot;
	if(Restart)	//If starting from the beginning, overwrite
	{
		TPlot.open(File);
		TPlot << argv[4] << "," << argv[5] << "," << argv[6] << "," << argv[9] << "," << argv[10] << endl;
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

			Par[2] = atof(argv[6]);

			strcpy(FileApp, File);
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
			Loop_Out1(Par, Temp, FileApp);
			Loop_Out2(Par, Temp, FileApp);
		}
	}

	return(0);
}

void Loop_Out1(long double Par[], int Temp, char File[])
{
	long double theta = M_PI/2.;
	long double on_shell, photon, on_shell_0, photon_0, stop;
	bool Manifest[303][101];
	ofstream oTable;
	ifstream iTable(File);
	int i, j;
	char Bin_c[11];
	long double Bin_n[9];
	long double k;

	for(i = 0; i < 202; i++)
	{
		for(int j = 0; j < 101; j++)
		{
			Manifest[i][j] = false;
		}
	}

	cerr << setprecision(18);
	j = 0;
	while(iTable.good())
	{
		iTable >> Bin_n[0] >> Bin_c[0] >> Bin_n[1] >> Bin_c[1] >> Bin_n[2] >> Bin_c[2] >> Bin_n[3];
		iTable.ignore(300,'\n');
		i = Bin_n[0];
		theta = Bin_n[2];
		if(0 <= i && i < 303 && 0 <= theta && theta <= M_PI)
			Manifest[i][int(theta*200./M_PI)] = true;
		if((Bin_n[0] == 0 || Bin_n[0] == 201 || Bin_n[0] == 202 || Bin_n[0] == 302) && (Bin_n[2] == 0 || Bin_n[2] > 1.57))
			cerr << File << "," << j << "," << Bin_n[0] << "," << Bin_n[2] << "," << Bin_n[3] << "," << Dispersion(Par, Temp, 0, Bin_n[1], Bin_n[2]).Value() << endl;
		if(float(rand())/float(RAND_MAX) < .002 && abs(Bin_n[3]/Dispersion(Par, Temp, 0, Bin_n[1], Bin_n[2])-1.) > 1e-7)
			cerr << File << "," << j << "," << Bin_n[0] << "," << Bin_n[2] << "," << Bin_n[3] << "," << Dispersion(Par, Temp, 0, Bin_n[1], Bin_n[2]).Value() << endl;
		j++;
	}
	i = 201;

	iTable.close();
	oTable.open(File, ios::app);

	oTable << setprecision(18);

	for(theta = 0; theta < M_PI*.502; theta += M_PI/200.)
	{
		on_shell = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		photon = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		on_shell_0 = .5*sqrt(Par[4]-pow(2.*Par[2],2));
		photon_0 = .5*sqrt(Par[4]);
		stop = isnan(photon)?50.:photon+50.;

		for(i = 0; i <= 200; i++)
		{
			if(!Manifest[i][int(theta*200./M_PI)])
			{
				k = k_i(i,on_shell,photon,stop,on_shell_0,photon_0);
				if(k < 100 && k >= 0)
					oTable << i << "," << k << "," << theta << "," << Dispersion(Par, Temp, 0, k, theta) << "," << k0_Int(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) << endl;
			}
		}
		if(!Manifest[i][int(theta*200./M_PI)])
		{
			k = k_i(i,on_shell,photon,stop,on_shell_0,photon_0);
			if(k < 100 && k>= 0)
			{
				oTable << i << "," << k << "," << theta << "," << Dispersion(Par, Temp, 0, k, theta) << "," << k0_Int(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) << endl;
			}
		}
	}

	oTable.close();
}

void Loop_Out2(long double Par[], int Temp, char File[])
{
	long double theta;
	long double on_shell, photon, on_shell_0, photon_0, stop;
	bool Manifest[303][101];
	ofstream oTable;
	ifstream iTable(File);
	int i = 201, j;
	long double Min = 201;
	char Bin_c[11];
	long double Bin_n[9];
	long double k;

	for(i = 0; i < 303; i++)
	{
		for(j = 0; j < 101; j++)
		{
			Manifest[i][j] = false;
		}
	}

	while(iTable.good())
	{
		iTable >> Bin_n[0] >> Bin_c[0] >> Bin_n[1] >> Bin_c[1] >> Bin_n[2] >> Bin_c[2] >> Bin_c[3] >> Bin_c[4];
		iTable.ignore(300,'\n');
		i = Bin_n[0];
		theta = Bin_n[2];
		if((('0' <= Bin_c[4] && Bin_c[4] <= '9') || Bin_c[4] == '.' ) && 0 <= i && i < 303 && 0 <= theta && theta <= M_PI)
			Manifest[i][int(theta*200./M_PI)] = true;
	}

	iTable.close();
	oTable.open(File, ios::app);

	for(i = 201; i >= 0; i--)
	{
		for(j = 0; j < 101; j++)
		{
			if(i < Min && !Manifest[i][j])
				Min = i;
		}
	}

	i = Min;
	Min = 100;
	for(theta = 0; theta < M_PI*.502; theta += M_PI/200.)
	{
		on_shell = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		photon = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		on_shell_0 = .5*sqrt(Par[4]-pow(2.*Par[2],2));
		photon_0 = .5*sqrt(Par[4]);
		stop = isnan(photon)?50.:photon+50.;
		if(Min > k_i(i,on_shell,photon,stop,on_shell_0,photon_0))
			Min = k_i(i,on_shell,photon,stop,on_shell_0,photon_0);
	}

	oTable << setprecision(18);
	for(i = 0; i <= 100; i++)
	{
		for(theta = 0; theta < M_PI*.502; theta += M_PI/200.)
		{
			if(!Manifest[i+202][int(theta*200./M_PI)])
			{
				k = Min+i*(100.-Min)/100.;
				oTable << i+202 << "," << k << "," << theta << "," << Dispersion(Par, Temp, 0, k, theta) << "," << k0_Int(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) << endl;
			}
		}
	}

	oTable.close();
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
