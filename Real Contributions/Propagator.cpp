//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<string>
#include<chrono>
#include"Around.h"
//#include"Spectral.h"
#include"Spectral adaptive.h"
//#include"Spectral dk0 ds.h"
//#include"Spectral dk0 ds adaptive.h"
//#include"Spectral f(k0 onshell) dk0 ds.h"
//#include"Spectral f(k0 onshell) dk0 ds adaptive.h"
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double ReG12(long double, long double, long double, long double, long double);
long double ImG12(long double, long double, long double, long double, long double);
long double k_i(int, long double, long double, long double);
void Loop_Out(long double[], int, char[]);

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
	char FileApp[70];
	char Number_c[5];
	string Number_s;
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

	cerr << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
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
			strcat(FileApp, ".m");
			Loop_Out(Par, Temp, FileApp);
		}
	}

	return(0);
}

void Loop_Out(long double Par[], int Temp, char File[])
{
	long double k, theta;
	long double on_shell, photon, stop;
	ofstream Table(File);
	int i;
	Table << setprecision(18);

	Table << "{" << flush;
	for(theta = 0; theta < M_PI*.505; theta += M_PI/200.)
	{
		on_shell = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		photon = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(sin(theta)*Par[3],2)));
		stop = isnan(photon)?50.:photon+50.;

		for(i = 0; i <= 700; i++)
		{
			k = k_i(i,on_shell,photon,stop);
			Table << "{" << i << "," << k << "," << theta << "," << Dispersion(Par, Temp, 0, k, theta) << "," << k0_Int(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta) << "," << ImG12(Par[2], Par[4], Par[3], k, theta) << "}," << endl;
		}
		k = k_i(i,on_shell,photon,stop);
		Table << "{" << i << "," << k << "," << theta << "," << Dispersion(Par, Temp, 0, k, theta) << "," << k0_Int(Par, Temp, k, theta) << "," << ReG12(Par[2], Par[4], Par[3], k, theta) << "," << ImG12(Par[2], Par[4], Par[3], k, theta) << "}" << flush;
		if(theta != M_PI/2.)
			Table << "," << endl;
	}
	Table << "}" << endl;

	Table.close();
}

long double ReG12(long double M, long double s, long double P, long double k, long double theta)
{
	return((2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)))*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(pow(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2),2)+pow(.14,2))));
}

long double ImG12(long double M, long double s, long double P, long double k, long double theta)
{
	return((2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)))*.14/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(pow(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2),2)+pow(.14,2))));
}

long double k_i(int i, long double x1, long double x2, long double x3)
{
	if(isnan(x2))
	{
		return(.1*i);
	}
	else if(isnan(x1) || x1 < .5)
	{
		long double a = -x2*x3/(120.*(x2-x3));
		long double b = (-6.*x2+x3)/(600.*(x2-x3));
		return(a*i/(1.+b*i));
	}
	long double a = (3.*x1*x2*x3)/(140.*(x1*(x2-6.*x3)+5.*x2*x3));
	long double b = -((3.*(7.*x1*x2-32.*x1*x3+15.*x2*x3))/(1400.*(x1*(x2-6.*x3)+5.*x2*x3)));
	long double c = (7.*x1*x2-12.*x1*x3+5.*x2*x3)/(140000.*(x1*(x2-6.*x3)+5.*x2*x3));
	return(a*i/(1.+b*i+c*pow(i,2)));
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