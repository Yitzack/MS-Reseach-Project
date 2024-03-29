//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<string>
#include<chrono>
#include<omp.h>
#include"Around.h"
#include"Spectral.cuh"
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
double ReG12(double, double, double, double, double);
double ImG12(double, double, double, double, double);
double k_i(int, double, double, double, double, double);
void Loop_Out1(double[], int, char[]);
void Loop_Out2(double[], int, char[]);
Around Int_Re_Insert(double Par[], int Temp, double k, double theta);
Around Int_Im_Insert(double Par[], int Temp, double k, double theta);

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

	char FileApp[130];
	char Number_c[5];
	string Number_s;
	strcat(File, argv[3]);	//Appends the temprature to the file name

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
	double Par[5];					//Parameters to be used in calculation {Coupling constant, potential cutoff, quark mass, P, s}

	cerr << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for double.
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
					Par[3] = ((double)(i%7)/8.-.125-floor((double)(i)/7.))*.8;
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
			//Loop_Out2(Par, Temp, FileApp);
		}
	}

	return(0);
}

void Loop_Out1(double Par[], int Temp, char File[])
{
	double theta = M_PI/2.;
	double on_shell, photon, on_shell_0, photon_0, stop;
	bool Manifest[303][101];
	//ofstream oTable;
	ifstream iTable(File);
	int i;
	char Bin_c[11];
	double Bin_n[9];
	double k;

	for(i = 0; i < 303; i++)
	{
		for(int j = 0; j < 101; j++)
		{
			Manifest[i][j] = false;
		}
	}

	iTable >> Bin_c[0];
	while(iTable.good())
	{
		iTable >> Bin_c[0] >> Bin_n[0] >> Bin_c[1] >> Bin_n[1] >> Bin_c[2] >> Bin_n[2];
		iTable.ignore(200,'\n');
		i = Bin_n[0];
		theta = Bin_n[2];
		Manifest[i][int(theta*200./M_PI)] = true;
	}

	iTable.close();
	//oTable.open(File, ios::app);

	cout << setprecision(18);

	for(theta = 0; theta < M_PI*.0025; theta += M_PI/200.)
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
					cout << i << "," << k << "," << theta << "," << Int_Re_Insert(Par, Temp, k, theta) << "," << Int_Im_Insert(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) <<  endl;
			}
		}
		if(!Manifest[i][int(theta*200./M_PI)])
		{
			k = k_i(i,on_shell,photon,stop,on_shell_0,photon_0);
			if(k < 100 && k>= 0)
			{
				cout << i << "," << k << "," << theta << "," << Int_Re_Insert(Par, Temp, k, theta) << "," << Int_Im_Insert(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) <<  endl;
			}
		}
	}

	//oTable.close();
}

void Loop_Out2(double Par[], int Temp, char File[])
{
	double theta;
	double on_shell, photon, on_shell_0, photon_0, stop;
	bool Manifest[303][101];
	ofstream oTable;
	ifstream iTable(File);
	int i = 201, j;
	double Min = 201;
	char Bin_c[11];
	double Bin_n[9];
	double k;

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
				oTable << i << "," << k << "," << theta << "," << Int_Re_Insert(Par, Temp, k, theta) << "," << Int_Im_Insert(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta, Temp) << "," << ImG12(Par[2], Par[4], Par[3], k, theta, Temp) <<  endl;
			}
		}
	}

	oTable.close();
}

Around Int_Re_Insert(double Par[], int Temp, double k, double theta)
{
	Dev_Pointer Pointers;
	cudaMalloc((void**)&Pointers.Par, 11*sizeof(double));
	cudaMalloc((void**)&Pointers.q, sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.omega, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Fermi, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.ImSelf, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.ReSelf, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Ordinate, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Limits, 30*sizeof(pair<double,double>));

	cudaStreamCreate(&Pointers.Stream);

	Around Answer = Dispersion(Pointers, Par, Temp, 0, k, theta);

	cudaFree(Pointers.Par);
	cudaFree(Pointers.q);
	cudaFree(Pointers.omega);
	cudaFree(Pointers.Fermi);
	cudaFree(Pointers.ImSelf);
	cudaFree(Pointers.ReSelf);
	cudaFree(Pointers.Ordinate);
	cudaFree(Pointers.Limits);
	cudaStreamDestroy(Pointers.Stream);

	return(Answer);
}

Around Int_Im_Insert(double Par[], int Temp, double k, double theta)
{
	Dev_Pointer Pointers;
	cudaMalloc((void**)&Pointers.Par, 11*sizeof(double));
	cudaMalloc((void**)&Pointers.q, sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.omega, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Fermi, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.ImSelf, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.ReSelf, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Ordinate, 1950*sizeof(pair<double,double>));
	cudaMalloc((void**)&Pointers.Limits, 30*sizeof(pair<double,double>));

	cudaStreamCreate(&Pointers.Stream);

	Around Answer = k0_Int(Pointers, Par, Temp, k, theta);

	cudaFree(Pointers.Par);
	cudaFree(Pointers.q);
	cudaFree(Pointers.omega);
	cudaFree(Pointers.Fermi);
	cudaFree(Pointers.ImSelf);
	cudaFree(Pointers.ReSelf);
	cudaFree(Pointers.Ordinate);
	cudaFree(Pointers.Limits);
	cudaStreamDestroy(Pointers.Stream);

	return(Answer);
}

double k_i(int i, double x1, double x2, double x3, double x1_0, double x2_0)
{
	if(isnan(x2_0) || x2_0 < .5)	//It needs to follow the policy of the smallest x2 or x3 that it can calculate
	{
		return(.1*i);
	}
	else if(isnan(x1_0) || x1_0 < .5)
	{
		double a = -x2*x3/(120.*(x2-x3));
		double b = (-6.*x2+x3)/(600.*(x2-x3));
		return(a*i/(1.+b*i));
	}
	double a = -x1*x3/(120.*(x1-x3));
	double b = (-6.*x1+x3)/(600.*(x1-x3));
	return(a*i/(1.+b*i));
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
