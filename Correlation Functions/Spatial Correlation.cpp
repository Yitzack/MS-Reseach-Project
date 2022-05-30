#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

long double Basis0(long double);	//BSpline Basis functions
long double Basis1(long double);
long double Basis2(long double);
long double Basis3(long double);
long double Basisn(long double);
void Import(char*, long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259]);

int main(int argc, char* argv[])
{
	long double ImG[151][259], ImGvC[151][259], ImGvL[151][259], ImGvQ[151][259], ImGV[151][259];	//Imaginary contributions to the spectral function
	long double ReGvC[151][259], ReGvL[151][259], ReGvQ[151][259], ReGV[151][259];			//Real contributions to the spectral function

	Import(argv[1], ImG, ImGvC, ImGvL, ImGvQ, ImGV, ReGvC, ReGvL, ReGvQ, ReGV);

	return(0);
}

long double Basis0(long double x)
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3);
	return(0);
}

long double Basis1(long double x)
{
	if(0 <= x && x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(0 <= x && x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

long double Basis2(long double x)
{
	if(0 <= x && x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(0 <= x && x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(0 <= x && x <= 4)
		return(-pow((x-4.)/2.,3);
	return(0);
}

long double Basis3(long double x)
{
	if(0 <= x && x < 2)
		return(pow(x,3)/24.);
	else if(0 <= x && x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(0 <= x && x < 4)
		return(11.pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(0 <= x && x <= 5)
		return(-pow((x-5.),3)/6.;
	return(0);
}

long double Basisn(long double x)
{
	if(2 <= x && x < 3)
		return(pow(x-2.,3)/6.);
	else if(2 <= x && x < 4)
		return(-pow(x,3)/2.+5*pow(x,2)-16.*x+50./3.);
	else if(2 <= x && x < 5)
		return(pow(x,3)/2.-7*pow(x,2)+32.*x-142./3.);
	else if(2 <= x && x <= 6)
		return(-pow((x-6.),3)/6.;
	return(0);
}

void Import(char* File, long double ImG[151][259], long double ImGvC[151][259], long double ImGvL[151][259], long double ImGvQ[151][259], long double ImGV[151][259], long double ReGvC[151][259], long double ReGcL[151][259], long double ReGvQ[151][259], long double ReGV[151][259])
{
	ifstream Input(char File);

	for(int i = 0; i < 259; i++)
		for(int j = 0; j < 151; j++)
			Input >> ImG[j][i] >> ImGvC[j][i] >> ImGvL[j][i] >> ImGvQ[j][i] >> ImGV[j][i] >> ReGvC[j][i] >> ReGvL[j][i] >> ReGvQ[j][i] >> ReGV[j][i];

	return;
}
