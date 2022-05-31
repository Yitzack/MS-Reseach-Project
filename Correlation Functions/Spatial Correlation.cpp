#include<iostream>
#include<iomanip>
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

long double Interacting(long double, long double, long double, long double, long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259]);	//Returns the Interacting Spectral Function
long double Non_interacting(long double, long double, long double[151][259]);						//Returns the Non-interacting Spectral Function

void Import(char*, long double&, long double&, long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259], long double[151][259]);

int main(int argc, char* argv[])
{
	long double ImG[151][259], ImGvC[151][259], ImGvL[151][259], ImGvQ[151][259], ImGV[151][259];	//Imaginary contributions to the spectral function
	long double ReGvC[151][259], ReGvL[151][259], ReGvQ[151][259], ReGV[151][259];			//Real contributions to the spectral function
	long double Coupling_Fraction, Vacuum_Coupling;							//Fraction of Vacuum Coupling Counstant, Vacuum Coupling Constant

	Import(argv[1], Coupling_Fraction, Vacuum_Coupling, ImG, ImGvC, ImGvL, ImGvQ, ImGV, ReGvC, ReGvL, ReGvQ, ReGV);

	cout << setprecision(18) << "(s,P) = (-2321.210625 GeV^2, 85.225 GeV), sigma_Non(s,P) = " << Non_interacting(-2321.210625, 85.225, ImG) << endl;//", sigma_Inter(s,P) = " << Interacting(-10, 5, Coupling_Fraction, Vacuum_Coupling, ImGvC, ImGvL, ImGvQ, ImGV, ReGvC, ReGvL, ReGvQ, ReGV) << endl;

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
	long double sigmaC = 3./M_PI/4.*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvC,2)-pow(ImGvC,2)*ImGV-2.*(ReGV-1.)*ImGvC*ReGvC)/(pow(ReGV-1.,2)+pow(ImGV,2));
	long double sigmaL = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvL,2)-pow(ImGvL,2)*ImGV-2.*(ReGV-1.)*ImGvL*ReGvL)/(pow(ReGV-1.,2)+pow(ImGV,2));
	long double sigmaQ = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvQ,2)-pow(ImGvQ,2)*ImGV-2.*(ReGV-1.)*ImGvQ*ReGvQ)/(pow(ReGV-1.,2)+pow(ImGV,2));

	return(sigmaC+sigmaL+sigmaQ);
}

long double Non_interacting(long double s, long double P, long double ImG[151][259])	//Returns the Non-interacting Spectral Function
{
	return(-3./M_PI*Interpolation(s, P, ImG));
}

long double Interpolation(long double s, long double P, long double f[151][259])
{
	long double i = i_index(s,P);						  //Index i
	long double j = j_index(s,P);						  //Index j
	int offset_i = -1, offset_j = -1;					  //Index offsets in calling up control points, normally -1, but can be 0 or -2 on the ends
	long double (*Basisi[4])(long double) = {Basisn,Basisn,Basisn,Basisn}; //Basis Functions in the i direction
	long double (*Basisj[4])(long double) = {Basisn,Basisn,Basisn,Basisn}; //Basis Functions in the j direction
	long double zx[4][4], zy[4][4];					  //z from the x direction and y direction
	long double answer = 0;						  //Result

	if(i < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisi[0] = Basis0;
		Basisi[1] = Basis1;
		Basisi[2] = Basis2;
		Basisi[3] = Basis3;
		if(i < 1)
			offset_i = 0;
	}
	else if(i < 3)
	{
		Basisi[0] = Basis1;
		Basisi[1] = Basis2;
		Basisi[2] = Basis3;
	}
	else if(i < 4)
	{
		Basisi[0] = Basis2;
		Basisi[1] = Basis3;
	}
	else if(i < 5)
	{
		Basisi[0] = Basis3;
	}
	else if(148 <= i)
	{
		Basisi[0] = Basis3;
		Basisi[1] = Basis2;
		Basisi[2] = Basis1;
		Basisi[3] = Basis0;
		if(149 <= i)
			offset_i = -2;
	}
	else if(147 <= i)
	{
		Basisi[1] = Basis3;
		Basisi[2] = Basis2;
		Basisi[3] = Basis1;
	}
	else if(146 <= i)
	{
		Basisi[2] = Basis3;
		Basisi[3] = Basis2;
	}
	else if(145 <= i)
	{
		Basisi[3] = Basis3;
	}

	if(j < 2 || (j >= 208 && j <= 210))	//Reassign the function pointers for the y/j direction
	{
		Basisj[0] = Basis0;
		Basisj[1] = Basis1;
		Basisj[2] = Basis2;
		Basisj[3] = Basis3;
		if(j < 1 || (j >= 208 && j < 209))
			offset_j = 0;
	}
	else if(j < 3 || (j >= 210 && j < 211))
	{
		Basisj[0] = Basis1;
		Basisj[1] = Basis2;
		Basisj[2] = Basis3;
	}
	else if(j < 4 || (j >= 211 && j < 212))
	{
		Basisj[0] = Basis2;
		Basisj[1] = Basis3;
	}
	else if(j < 5 || (j >= 212 && j < 213))
	{
		Basisj[0] = Basis3;
	}
	else if((206 <= j && j < 208) || 256 <= j)
	{
		Basisj[0] = Basis3;
		Basisj[1] = Basis2;
		Basisj[2] = Basis1;
		Basisj[3] = Basis0;
		if((207 <= j && j < 208) || 257 <= j)
			offset_j = -2;
	}
	else if((205 <= j && j < 208) || 255 <= j)
	{
		Basisj[1] = Basis3;
		Basisj[2] = Basis2;
		Basisj[3] = Basis1;
	}
	else if((204 <= j && j < 208) || 254 <= j)
	{
		Basisj[2] = Basis3;
		Basisj[3] = Basis2;
	}
	else if((203 <= j && j < 208) || 253 <= j)
	{
		Basisj[3] = Basis3;
	}
cout << j << "," << i << endl;
	for(int i_count = 0; i_count < 4; i_count++)	//Evaluate the Basis Functions
		for(int j_count = 0; j_count < 4; j_count++)
		{
			if(Basisj[j_count] != Basisn && j < 5)
				zy[i_count][j_count] = Basisj[j_count](j);
			else if(Basisj[j_count] != Basisn && (208 <= j && j < 213))
				zy[i_count][j_count] = Basisj[j_count](j-208);
			else if(Basisj[j_count] != Basisn && (203 <= j && j < 208))
				zy[i_count][j_count] = Basisj[j_count](208-j);
			else if(Basisj[j_count] != Basisn && 253 <= j)
				zy[i_count][j_count] = Basisj[j_count](258-j);
			else
				zy[i_count][j_count] = Basisj[j_count](i-int(i)+j_count+2.);	//Some how the fractional part of the arguments were trading dimensions. Dispite trying to be careful and keeping x/i in the first argument and y/j in the second argument, I was ending up with f(int(i)+j-int(j),int(j)+i-int(i)) instead of f(int(i)+i-int(i),int(j)+j-int(j))=f(i,j)

			if(Basisi[i_count] != Basisn && i < 5)
				zx[i_count][j_count] = Basisi[i_count](i);
			else if(Basisi[i_count] != Basisn && 145 <= i)
				zx[i_count][j_count] = Basisi[i_count](150-i);
			else
				zx[i_count][j_count] = Basisi[i_count](j-int(j)+i_count+2.);
cout << i << " " << j << " " << i_count << " " << j_count << " " << zx[i_count][j_count] << " " << zy[i_count][j_count] << " " << int(j)+j_count+offset_j << " " << int(i)+i_count+offset_i << " " << f[int(i)+i_count+offset_i][int(j)+j_count+offset_j] << endl;
			answer += zx[i_count][j_count]*zy[i_count][j_count]*f[int(i)+i_count+offset_i][int(j)+j_count+offset_j];
		}

	return(answer);
}

long double j_index(long double s, long double P)
{
	if(s+pow(P,2) <= 432.64)
		return(10.*sqrt(s+pow(P,2)));
	else
		return(187.2+sqrt(s+pow(P,2)));
}

long double i_index(long double s, long double P)
{
	return(10.*(P-sqrt(s+pow(P,2))));
}

long double Basis0(long double x)
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3));
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
		return(-pow((x-4.)/2.,3));
	return(0);
}

long double Basis3(long double x)
{
	if(0 <= x && x < 2)
		return(pow(x,3)/24.);
	else if(0 <= x && x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(0 <= x && x < 4)
		return(11.*pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(0 <= x && x <= 5)
		return(-pow((x-5.),3)/6.);
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
		return(-pow((x-6.),3)/6.);
	return(0);
}

void Import(char* File, long double& Coupling_Fraction, long double& Vacuum_Coupling, long double ImG[151][259], long double ImGvC[151][259], long double ImGvL[151][259], long double ImGvQ[151][259], long double ImGV[151][259], long double ReGvC[151][259], long double ReGvL[151][259], long double ReGvQ[151][259], long double ReGV[151][259])
{
	ifstream Input(File);

	Input >> Coupling_Fraction >> Vacuum_Coupling;

	for(int i = 0; i < 259; i++)
		for(int j = 0; j < 151; j++)
			Input >> ImG[j][i] >> ImGvC[j][i] >> ImGvL[j][i] >> ImGvQ[j][i] >> ImGV[j][i] >> ReGvC[j][i] >> ReGvL[j][i] >> ReGvQ[j][i] >> ReGV[j][i];

	return;
}
