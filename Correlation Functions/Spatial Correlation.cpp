#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

long double Interpolation(long double, long double, long double[151][259]);	//(s,P) comes in, f(s,P) comes out
long double i_index(long double, long double);				//Turns (s,P) into index variables
long double j_index(long double, long double);
long double Basis0(long double);						//BSpline Basis functions
long double Basis1(long double);
long double Basis2(long double);
long double Basis3(long double);
long double Basisn(long double);

long double Interacting(long double, long double, long double, long double, long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259]);	//Returns the Interacting Spectral Function
long double Non_interacting(long double, long double, long double[151][259]);						//Returns the Non-interacting Spectral Function

void Import(char*, long double, long double, long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259]);

int main(int argc, char* argv[])
{
	long double ImG[151][259], ImGvC[151][259], ImGvL[151][259], ImGvQ[151][259], ImGV[151][259];	//Imaginary contributions to the spectral function
	long double ReGvC[151][259], ReGvL[151][259], ReGvQ[151][259], ReGV[151][259];			//Real contributions to the spectral function
	long double Coupling_Fraction, Vacuum_Coupling;							//Fraction of Vacuum Coupling Counstant, Vacuum Coupling Constant

	Import(argv[1], Coupling_Fraction, Vacuum_Coupling, ImG, ImGvC, ImGvL, ImGvQ, ImGV, ReGvC, ReGvL, ReGvQ, ReGV);

	cout << "(s,P) = (-10 GeV^2,5 GeV), sigma_Inter(s,P) = " << Interacting(-10, 5, Coupling_Fraction, Vacuum_Coupling, ImGvC, ImGvL, ImGvQ, ImGV, ReGvC, ReGvL, ReGvQ, ReGV) << "sigma_Non(s,P) = " << Non_interacting(s, P, ImG) << endl;

	return(0);
}

long double Interacting(long double s, long double P, long double Coupling_Fraction, long double Vacuum_Coupling, long double ImGvCf[151][259], long double ImGvLf[151][259], long double ImGvQf[151][259], long double ImGVf[151][259], long double ReGvCf[151][259], long double ReGvLf[151][259], long double ReGvQf[151][259], long double ReGVf[151][259])	//Returns the Interacting Spectral Function
{
	long double ImGvC = Interpolation(s,P,ImGvCf);	//Get the interpolation of the contributions
	long double ImGvL = sqrt(1.5)*Interpolation(s,P,ImGvCf);
	long double ImGvQ = Interpolation(s,P,ImGvCf)/2.;
	long double ImGV = Interpolation(s,P,ImGVf);
	long double ReGvC = Interpolation(s,P,ReGvCf);
	long double ReGvL = Interpolation(s,P,ReGvLf);
	long double ReGvQ = Interpolation(s,P,ReGvQf);
	long double ReGV = Interpolation(s,P,ReGVf);

	//Each of the 3 interacting spectral function contributions
	long double sigmaC = 3./M_PI/4.*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvC,2)-pow(ImGvC,2)*ImGV-2.*(ReGV-1.)*ImGvC*ReGvC)/(pow((ReGV-1.,2)+pow(ImGV,2));
	long double sigmaL = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvL,2)-pow(ImGvL,2)*ImGV-2.*(ReGV-1.)*ImGvL*ReGvL)/(pow((ReGV-1.,2)+pow(ImGV,2));
	long double sigmaQ = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvQ,2)-pow(ImGvQ,2)*ImGV-2.*(ReGV-1.)*ImGvQ*ReGvQ)/(pow((ReGV-1.,2)+pow(ImGV,2));

	return(sigmaC+sigmaL+sigmaQ);
}

long double Non_interacting(long double s, long double P, long double ImG[151][259])	//Returns the Non-interacting Spectral Function
{
	return(-3./M_PI*Interpolation(s, P, ImG));
}

long double Interpolation(long double s, long double P, long double f[151][259])
{
	long double i = i_index(s,P);							//Index i
	long double j = j_index(s,P);							//Index j
	long double (*Basisi)(long double)[4] = {Basisn,Basisn,Basisn,Basisn};	//Basis Functions in the i direction
	long double (*Basisj)(long double)[4] = {Basisn,Basisn,Basisn,Basisn};	//Basis Functions in the j direction
	long double zx[4][4], zy[4][4];						//z from the x direction and y direction
	long double answer = 0;							//Result

	if(i < 2 || (i >= 208 && i <= 210))	//Reassign the function pointers
	{
		Basisi[0] = Basis0;
		Basisi[1] = Basis1;
		Basisi[2] = Basis2;
		Basisi[3] = Basis3;
	}
	else if(i < 3 || (i >= 210 && i < 211))
	{
		Basisi[0] = Basis1;
		Basisi[1] = Basis2;
		Basisi[2] = Basis3;
	}
	else if(i < 4 || (i >= 211 && i < 212))
	{
		Basisi[0] = Basis2;
		Basisi[1] = Basis3;
	}
	else if(i < 5 || (i >= 212 && i < 213))
	{
		Basisi[0] = Basis3;
	}
	else if((206 <= i && 208 < i) || 256 <= i)
	{
		Basisi[0] = Basis3;
		Basisi[1] = Basis2;
		Basisi[2] = Basis1;
		Basisi[3] = Basis0;
	}
	else if((205 <= i && 207 < i) || 255 <= i)
	{
		Basisi[1] = Basis3;
		Basisi[2] = Basis2;
		Basisi[3] = Basis1;
	}
	else if((204 <= i && 206 < i) || 254 <= i)
	{
		Basisi[2] = Basis3;
		Basisi[3] = Basis2;
	}
	else if((203 <= i && 205 < i) || 253 <= i)
	{
		Basisi[3] = Basis3;
	}
	if(j < 2)
	{
		Basisi[0] = Basis0;
		Basisj[1] = Basis1;
		Basisj[2] = Basis2;
		Basisj[3] = Basis3;
	}
	else if(j < 3)
	{
		Basisj[0] = Basis1;
		Basisj[1] = Basis2;
		Basisj[2] = Basis3;
	}
	else if(j < 4)
	{
		Basisj[0] = Basis2;
		Basisj[1] = Basis3;
	}
	else if(j < 5)
	{
		Basisj[0] = Basis3;
	}
	else if(148 <= j)
	{
		Basisi[0] = Basis3;
		Basisj[1] = Basis2;
		Basisj[2] = Basis1;
		Basisj[3] = Basis0;
	}
	else if(147 <= j)
	{
		Basisj[1] = Basis3;
		Basisj[2] = Basis2;
		Basisj[3] = Basis1;
	}
	else if(146 <= j)
	{
		Basisj[2] = Basis3;
		Basisj[3] = Basis2;
	}
	else if(145 <= j)
	{
		Basisj[3] = Basis3;
	}

	for(int i_count = 0; i_count < 4; i_count++)	//Evaluate the Basis Functions
		for(int j_count = 0; j_count < 4; j_count++)
		{
			if(Basisi[i_count] != Basisn && i < 5)
				zx[i_count][j_count] = Basisi[i_count](j);
			else if(Basisi[i_count] != Basisn && (i >= 208 && i < 213))
				zx[i_count][j_count] = Basisi[i_count](j-208);
			else if(Basisi[i_count] != Basisn && (203 <= i && 208 < i))
				zx[i_count][j_count] = Basisi[i_count](208-j);
			else if(Basisi[i_count] != Basisn && 253 <= i)
				zx[i_count][j_count] = Basisi[i_count](253-j);
			else
				zx[i_count][j_count] = Basisi[i_count](j-j_count);

			if(Basisi[i_count] != Basisn && i < 5)
				zy[i_count][j_count] = Basisj[j_count](i);
			else if(Basisi[i_count] != Basisn && 145 <= i)
				zy[i_count][j_count] = Basisj[j_count](150-i);
			else
				zy[i_count][j_count] = Basisj[j_count](i-i_count);

			answer += zx[i_count][j_count]*zy[i_count][j_count]*f[int(i)+i_count][int(j)+j_count];
		}

	return(answer);
}

long double i_index(long double s, long double P)
{
	if(s+pow(P,2) <= 432.64)
		return(10.*sqrt(s+pow(P,2)));
	else
		return(187.2+sqrt(s+pow(P,2)));
}

long double j_index(long double s, long double P)
{
	return(10.*(P-sqrt(s+pow(P,2))));
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

void Import(char* File, long double Coupling Fraction, long double Vacuum_Coupling, long double ImG[151][259], long double ImGvC[151][259], long double ImGvL[151][259], long double ImGvQ[151][259], long double ImGV[151][259], long double ReGvC[151][259], long double ReGcL[151][259], long double ReGvQ[151][259], long double ReGV[151][259])
{
	ifstream Input(char File);

	Input >> Coupling_Fraction >> Vacuum_Coupling;

	for(int i = 0; i < 259; i++)
		for(int j = 0; j < 151; j++)
			Input >> ImG[j][i] >> ImGvC[j][i] >> ImGvL[j][i] >> ImGvQ[j][i] >> ImGV[j][i] >> ReGvC[j][i] >> ReGvL[j][i] >> ReGvQ[j][i] >> ReGV[j][i];

	return;
}
