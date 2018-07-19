//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
#include<cfloat>
#include"Elements.h"
using namespace std;

//Integrals that define results and ancillary functions
Elements theta_Int(long double[5], int);	//Integrates the theta results
Elements k_Int(long double[5], int, long double);	//Integrates the k momentum results
Elements Folding(long double[5], int, long double, long double);	//Folding integral, energy integral
void Characterize_k_Int(long double[5], int, long double, long double[4], long double[4], int&); //Returns the poles of the k integral's integrands
void Newtons_k_Int(long double, long double, long double, long double, long double[3], long double);	//Newton's method of finding roots for Characterize_k_Int
int Newtons_Test_k_Int(long double, long double, long double, long double, long double[3], long double);	//Trys to determine the number of roots Newton's method  is looking for and where they might be found
void Characterize_Folding(long double[5], int, long double, long double, long double[4], long double[4], int&); //Returns the poles of the folding integral's integrands

//Straight Functions everything is built from
long double ImSelf_Energy(long double, long double, long double, int);	//Single quark self energy
long double ReSelf_Energy(long double, long double, long double, int);	//Single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Potential_on(long double[5]);	//On-shell potential for the on-shell T-Matrix
long double Potential1(long double[5], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[5], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double Quark_Spectrum(long double, long double, long double, int);	//Single quark spectral function
long double Spin_Sum(long double[5], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Folding_Integrand(long double[5], long double, long double, long double, int);	//Integrand of the folding integral

#define GAMMA -.015

//long double Par[5] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[5], int Temp)	//Integrates the theta results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {.01,.05,.1,M_PI/2.,M_PI-.1,M_PI-.05,M_PI-.01,M_PI};
	long double x1, x2;	//Abscissa
	Elements F_a, F_b, F_ave;	//Sum of ordinate*weights
	Elements Answer(0,0,0);	//Answer to be returned
	long double a = 0, b;	//Sub-interval limits of integration
	int i, j;	//Counters

	for(i = 0; i < 8; i++)
	{
		b = Range[i];

		F_a.null();
		F_b.null();
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += k_Int(Par, Temp, x1)*sin(x1)*w[j+1];
			F_b += k_Int(Par, Temp, x2)*sin(x2)*w[j+1];
		}
		F_ave = k_Int(Par, Temp, (a+b)/2.)*sin((a+b)/2.)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;
		a = b;
	}

	return(Answer/pow(2.*M_PI,2));
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements k_Int(long double Par[5], int Temp, long double theta)	//Integrates the k momentum results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	Elements F_a, F_b, F_ave;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Answer to be returned
	Elements PartialAnswer;	//Answer for sub-interval for determining completeness
	long double x1, x2;	//Abscissa
	long double a = 0, b = 0;//Sub-interval limits of integration
	int Poles;	//Number of poles
	long double zero[3];	//The real part of the signular pole
	long double gamma[3];	//The distance to the singular, maybe
	long double Min_upper;	//The integral has to go at least this far
	long double Width;	//Length of the next sub-interval
	long double Early = 0;	//Early change from one pole to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;//The next width that will be used in the event of an early change of poles
	int i = 0, j, l;	//Counters
	Elements holder;

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);

	Min_upper = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//This the upper bound that the vacuum calls for, Partial/total will promote higher as needed

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0) //j!=0 is because the loop will find other early termination
		Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];

	do
	{
		if(b == 0 && j != 0)	//First pole is closer than zero-64*gamma to lower limit of integration
			Width = zero[i]+Range[j]*gamma[i]-b;
		else if((i < Poles && b+100 < zero[i]-64.*gamma[i]) || a-100 > Min_upper || b+100 < Min_upper || (i > 0 && a-100 > zero[i-1]+64.*gamma[i-1]))	//Middle of nowhere intervals
			Width = 100;
		else if((i < Poles && b+50 < zero[i]-64.*gamma[i]) || a-50 > Min_upper || b+50 < Min_upper || (i > 0 && a-50 > zero[i-1]+64.*gamma[i-1]))
			Width = 50;
		else if((i < Poles && b+10 < zero[i]-64.*gamma[i]) || a-10 > Min_upper || b+10 < Min_upper || (i > 0 && a-10 > zero[i-1]+64.*gamma[i-1]))
			Width = 10;
		else
			Width = 3;

		if(j == 8 && i < Poles) //Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i])) //Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Poles != 0 && j < 8 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && (b != 0 || Width == 0))
		{//Integrating a pole and the width wan't resolved by the first condition after do, 24 lines up
			Width = gamma[i]*(Range[j+1]-Range[j]);
			j++;
		}

		if(NextWidth != 0 && NextWidth == NextWidth)	//Resolving early termination of integrating one pole for the next
		{
			Width = NextWidth;
			NextWidth = 0;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Early != 0 && b+Width > Early)
		{
			Width = Early-b;
			Early = 0;
			j = 8-j;
			NextWidth = Width*gamma[i+1]/gamma[i];
			i++;
		}

		b += Width;	//Set next b keeping in mind the upper limit of integration

		if(b > Min_upper && a < Min_upper)	//If (a,b) incompasses Min_upper, make b equal to Min_upper
			b = Min_upper;

		F_a.null();
		F_b.null();
		for(l = 0; l < 9; l++)
		{
			x1 = (b+a-Disp[l]*(b-a))/2.; //Actual evaluation points
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F_a += Folding(Par, Temp, x1, theta)*pow(x1,2)*w[l+1]; //Evaluate function at x1
			F_b += Folding(Par, Temp, x2, theta)*pow(x2,2)*w[l+1]; //Evaluate function at x2
		}
		F_ave = Folding(Par, Temp, (a+b)/2., theta)*pow((a+b)/2.,2)*w[0]; //Evaluate function at (a+b)/2.
		PartialAnswer = (F_a+F_ave+F_b)*(b-a)/(2.);
		Answer += PartialAnswer;
		a = b;
	}while(b <= Min_upper || abs(PartialAnswer/Answer) >= .0001);

	return(Answer);
}

void Characterize_k_Int(long double Par[5], int Temp, long double theta, long double zero[3], long double gamma[3], int &Poles) //Returns the poles of the k integral's integrands
{
	long double holder;
	long double previous[2];
	int i, j = 0, l;

	if(sqrt(Par[4]) > 2.*Par[2])	//Find the poles and estimate a width (distance to pole)
		zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
	else
		zero[2] = 0;

	zero[0] = zero[1] = zero[2];	//Start the Newton's search for other k poles
	i = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta);
	if(i)
	{
		do
		{
			previous[0] = zero[0];
			previous[1] = zero[1];
			Newtons_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta);
			j++;
		}while((abs(previous[0]/zero[0]-1.) > .001 || abs(previous[1]/zero[1]-1.) > .001) && ((i >= 2 && j <= 10) || i <= 1));	//Keep going while both poles are not known better than 1MeV

		if(j <= 10)
		{
			zero[0] = abs(zero[0]);	//Through investigation of the results, I found that if both peaks are on one condition,
			zero[1] = abs(zero[1]);	//then the lower one will be negative on the other condition
			Poles = 3;
		}
		else
		{
			zero[0] = zero[2];
			gamma[0] = Folding(Par, Temp, zero[0], theta).Min();
			if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
				gamma[0] = 1e-3;
			Poles = 1;
			return;
		}
	}
	else
	{
		zero[0] = zero[2];
		gamma[0] = Folding(Par, Temp, zero[0], theta).Min();
		if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[0] = 1e-3;
		Poles = 1;
		return;
	}

	for(i = 2; i >= 0; i--)	//Bubble sort
	{
		for(j = 0; j < i; j++)
		{
			if(zero[j] > zero[j+1])
			{
				holder = zero[j+1];
				zero[j+1] = zero[j];
				zero[j] = holder;
			}
			else if(zero[j] == zero[j+1])	//Remove duplicates
			{
				for(l = j; l < 2; l++)
					zero[l] = zero[l+1];
				zero[2] = sqrt(Par[4]+pow(Par[3],2)); //Need to use the biggest finite value of reason or it will be attempted to be sorted to the bottom when it is invalid
				Poles--;
			}
		}
	}

	for(i = 0; i < Poles; i++)
	{
		gamma[i] = Folding(Par, Temp, zero[i], theta).Min();
		if(gamma[i] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[i] = 1e-3;
	}
	return;
}

void Newtons_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = Energy(M,P/2.,k[0],theta)-.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[0],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[0],4))));
	long double fp1 = pow(pow(2.*k[0],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[0],4)*(4.+sqrt(16.+pow(Lambda/k[0],4)))))/(k[0]*sqrt(8.+pow(Lambda/k[0],4)/2.))+(2.*k[0]+P*cos(theta))/(2.*Energy(M,P/2.,k[0],theta));
	long double f2 = -Energy(M,P/2.,-k[1],theta)-.5*(-sqrt(s+pow(P,2))+pow(pow(2.*k[1],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[1],4))));
	long double fp2 = -pow(pow(2.*k[1],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[1],4)*(4.+sqrt(16.+pow(Lambda/k[1],4)))))/(k[1]*sqrt(8.+pow(Lambda/k[1],4)/2.))-(2.*k[1]-P*cos(theta))/(2.*Energy(M,P/2.,-k[1],theta));

	k[0] -= f1/fp1;
	k[1] -= f2/fp2;
	return;
}

int Newtons_Test_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = Energy(M,P/2.,k[2],theta)-.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[2],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[2],4))));
	long double fp1 = pow(pow(2.*k[2],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4)))))/(k[2]*sqrt(8.+pow(Lambda/k[2],4)/2.))+(2.*k[2]+P*cos(theta))/(2.*Energy(M,P/2.,k[2],theta));
	long double limit1 = sqrt(pow(M,2)+pow(P,2)/4.)-.5*sqrt(s+pow(P,2))+pow(2,-1.5)*Lambda;
	long double f2 = -Energy(M,P/2.,-k[2],theta)+.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[2],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[2],4))));
	long double fp2 = -pow(pow(2.*k[2],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4)))))/(k[2]*sqrt(8.+pow(Lambda/k[2],4)/2.))-(2.*k[2]-P*cos(theta))/(2.*Energy(M,P/2.,k[2],theta));
	long double limit2 = -sqrt(pow(M,2)+pow(P,2)/4.)+.5*sqrt(s+pow(P,2))-pow(2,-1.5)*Lambda;
	long double limitp = P*cos(theta)/sqrt(pow(2.*M,2)+pow(P,2));
	int roots = 0;

	if(f1*limit1 < 0 && f2*limit2 < 0)
		roots = 1;
	else if(fp1*limitp < 0 && limit1*limitp < 0)
		roots = 2;
	else if(fp2*limitp < 0 && limit2*limitp < 0)
		roots = 3;
	else if(f1*limit1 < 0 || f2*limit2 < 0)
		cerr << "Unexepected outcome in Newton_Test_k_Int s=" << s << " P=" << P << " theta=" << theta << endl;

	if(roots == 0)
	{
		k[0] = 0;
		k[1] = 0;
	}

	return(roots);
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements Folding(long double Par[5], int Temp, long double k, long double theta)	//Folding integral, energy integral
{
	if(Temp == 0 && sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta) <= Energy(0,Par[3]/2.,k,theta))	//Let's save some time and just return 0, because it is
		return(Elements(0,0,0));
	else if(Par[4]+pow(Par[3],2) <= 0)
		return(Elements(0,0,0));	//Bad data trap and time saver

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	long double a, b;	//Sub-interval limits of integration
	long double Width;	//Next step size
	long double Max;	//Upper limit of integration
	long double Early = 0;	//Early change from one pole to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;//The next width that will be used in the event of an early change of poles
	Elements F_a, F_b, F_ave;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Results to be returned
	long double x1, x2;	//Abscissa
	long double zero[4];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[4];	//Imaginary part of poles
	long double Stutter = k+.5*sqrt(Par[4]+pow(Par[3],2));	//Marks a discontinuty in the potiential
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles

	Characterize_Folding(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about

	if(Temp != 0)	//Assign limits of integration
	{
		a = b = 0;
		Max = sqrt(Par[4]+pow(Par[3],2));
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta);
		Max = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0) //j!=0 is because the loop will find other early termination
		Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];

	do
	{
		if((b == 0 || b == Energy(0,Par[3]/2.,k,theta)) && j != 0)	//First pole is closer than zero-64*gamma to lower limit of integration
			Width = zero[i]+Range[j]*gamma[i]-b;
		else if((i < Poles && b+100 < zero[i]-64.*gamma[i]) || b+100 < Max)	//Middle of nowhere intervals
			Width = 100;
		else if((i < Poles && b+50 < zero[i]-64.*gamma[i]) || b+50 < Max)
			Width = 50;
		else if((i < Poles && b+10 < zero[i]-64.*gamma[i]) || b+10 < Max)
			Width = 10;
		else
			Width = 3;

		if(j == 8 && i < Poles) //Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i])) //Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Poles != 0 && j < 8 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && !(b == 0 || b == Energy(0,Par[3]/2.,k,theta)))
		{//Integrating a pole and the width wan't resolved by the first condition after do, 24 lines up
			Width = gamma[i]*(Range[j+1]-Range[j]);
			j++;
		}

		if(NextWidth != 0 && NextWidth == NextWidth)	//Resolving early termination of integrating one pole for the next
		{
			Width = NextWidth;
			NextWidth = 0;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Early != 0 && b+Width > Early)
		{
			Width = Early-b;
			Early = 0;
			j = 8-j;
			NextWidth = Width*gamma[i+1]/gamma[i];
			i++;
		}

		b += Width;	//Set next b keeping in mind the upper limit of integration
		if(b > Max)
			b = Max;

		reStutter:	//Use of goto statement to split interval in to two parts that puts the discontinutiy on either limit of integration. This is probably the similest way to do this without screwing everything else up. *Not advised*
		if(b == Stutter && b-a < Width)	//Catch the second pass of the reStutter step first as it is easier catch here than after the initiation
		{
			b = a + Width;
			a = Stutter;
		}
		else if(b > Stutter && a < Stutter) //reStutter loop initiated
			b = Stutter;

		F_a.null();
		F_b.null();
		for(l = 0; l < 9; l++)	//Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F_a += Elements(Spin_Sum(Par, x1, k, theta), Potential1(Par,x1,k), Potential2(Par,x1,k))*Folding_Integrand(Par,x1,k,theta,Temp)*w[l+1];
			F_b += Elements(Spin_Sum(Par, x2, k, theta), Potential1(Par,x2,k), Potential2(Par,x2,k))*Folding_Integrand(Par,x2,k,theta,Temp)*w[l+1];
		}
		F_ave = Elements(1, Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*Folding_Integrand(Par,(a+b)/2.,k,theta,Temp)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;

		if(b == Stutter && b-a < Width)	//execute the second step of the reStutter loop
			goto reStutter;

		a = b;
	}while(b < Max);

	return(Answer/M_PI);
}

void Characterize_Folding(long double Par[5], int Temp, long double k, long double theta, long double zero[4], long double gamma[4], int &Poles)
{
	long double Lower, Upper; //Limits of integration in Folding, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(Temp != 0)
	{
		Lower = 0;
		Upper = sqrt(Par[4]+pow(Par[3],2));
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta);
		Upper = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	zero[0] = .5*(sqrt(Par[4]+pow(Par[3],2))+pow(pow(2.*k,4)+pow(Par[1],4),.25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));	//Potential poles, I know exactly where these are at.
	zero[1] = .5*(sqrt(Par[4]+pow(Par[3],2))-pow(pow(2.*k,4)+pow(Par[1],4),.25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[0] = .5*(pow(pow(2.*k,4)+pow(Par[1],4),.25)*sin(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[1] = .5*(pow(pow(2.*k,4)+pow(Par[1],4),.25)*sin(.5*atan(pow(Par[1]/(2.*k),2))));

	holder = GAMMA;
	zero[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));	//Exact vacuum
	zero[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));	//Start of exact vacuum

	if(Temp != 0)	//media estimate
	{
		holder = Self_Energy(Par[2], zero[2], Energy(0,Par[3]/2.,k,theta), Temp)/sqrt(pow(zero[2],2)-pow(Energy(0,Par[3]/2.,k,theta),2));
		zero[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));
		gamma[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));

		holder = Self_Energy(Par[2], zero[3], Energy(0,Par[3]/2.,-k,theta), Temp)/sqrt(pow(zero[3],2)-pow(Energy(0,Par[3]/2.,k,theta),2));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
		gamma[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));
		gamma[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-zero[3];
	}

	for(i = 3; i >= 0; i--)	//Bubble sort
	{
		for(j = 0; j < i; j++)
		{
			if(zero[j] > zero[j+1])
			{
				holder = zero[j+1];
				zero[j+1] = zero[j];
				zero[j] = holder;
				holder = gamma[j+1];
				gamma[j+1] = gamma[j];
				gamma[j] = holder;
			}
		}
	}

	i = j = 0;	//Find the first zero greater than 0
	while(zero[i] < Lower) i++;

	while(zero[i] <= Upper && i < 4)	//Move zeroes up to front of array, count off poles within the limits of integration
	{
		zero[j] = zero[i];
		gamma[j] = gamma[i];
		i++;
		j++;
	}
	Poles = j;

	return;
}

//long double Par[5] = {g, Lambda, M, P, s}
long double ImSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
#ifdef RIEK
	long double Sigma;
	long double Delta;
	long double b1, b2;
	long double M1,M_T;
	long double a;
	long double sigma1, sigma2;
	long double Floor;
	
	switch(Temp)
	{
		case 0:
			if(pow(omega,2)>=pow(k,2))
				return(sqrt(pow(omega,2)-pow(k,2))*GAMMA);
			else
				return(0);
			break;
		case 1:
			Sigma = -1.12031;
			Delta = -1.11823;
			b1 = 3.42062;
			b2 = 2.34212;
			M_T = 1.848;
			M1 = 1.55744;
			a = 2.2094;
			sigma1 = 2.574538196654789;
			sigma2 = 2.5745381961599816;
			Floor = .928942;
			break;
		case 2:
			Sigma = -.096194;
			Delta = 1.16782;
			b1 = 3.27672;
			b2 = 1.87149;
			M_T = 1.719;
			M1 = 1.54196;
			a = .760248;
			sigma1 = 3.41189;
			sigma2 = 1.38003;
			Floor = .818425;
			break;
		case 3:
			Sigma = -.0995278;
			Delta = 1.29695;
			b1 = 3.21382;
			b2 = 1.35479;
			M_T = 1.563;
			M1 = 1.46014;
			a = .752122;
			sigma1 = 3.38156;
			sigma2 = 1.30001;
			Floor = .905761;
			break;
	}
	Floor = 0;
	long double Shift = M-M_T;
	M1 += Shift;

	if(pow(omega,2)>=pow(k,2))
		return(2.*M*Sigma*exp(Delta+(b1-b2)*(omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2))/2.-sqrt(b1*b2*pow((omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2)),2)+pow(Delta+(b1-b2)*(omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2))/2.,2)))*(a*exp(-pow(k/sigma1,2))+(1.-a)*exp(-pow(k/sigma2,2)))+sqrt(pow(omega,2)-pow(k,2))*GAMMA);
	else
		return(2.*M*Sigma*exp(Delta+(b1-b2)*(omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2))/2.-sqrt(b1*b2*pow((omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2)),2)+pow(Delta+(b1-b2)*(omega-sqrt(pow(M1,2)+pow(k,2)))*sqrt(pow(M1,2)+pow(k,2))/2.,2)))*(a*exp(-pow(k/sigma1,2))+(1.-a)*exp(-pow(k/sigma2,2))));
#endif
#if defined(SHUAI) || defined(CC) || defined(BB)
	long double omega0;	//location of central peak
	long double Sigma;	//size of energy dependance
	long double a, b;	//slope of exponential decrease to left and right
	long double knee;	//space to change from left to right side of peak
	long double M_T, Shift;

	switch(Temp)
	{
		case 0:
			if(pow(omega,2)>=pow(k,2))
				return(sqrt(pow(omega,2)-pow(k,2))*GAMMA);
			else
				return(0);
			break;
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .585335/sqrt(pow(k,2)+pow(1.85515,2));
			a = 5.05953/(pow(k,2)+pow(1.11686,2))+6.09943;
			b = -8.08693/(pow(k,2)+pow(2.70494,2))+3.25177;
			omega0 = sqrt(pow(1.49006+Shift,2)+pow(k,2))+.248573;
			knee = 3.84788*pow(k+1.,(long double)-.335162);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .660137/sqrt(pow(k,2)+pow(1.90299,2));
			a = 2.82635/(pow(k,2)+pow(.916643,2))+4.19118;
			b = -83.3834/(pow(k,2)+pow(8.8641,2))+2.93508;
			omega0 = sqrt(pow(1.45524+Shift,2)+pow(k,2))+.247213;
			knee = 3.29189*pow(k+1.,(long double)-.575497);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .670397/sqrt(pow(k,2)+pow(1.96561,2));
			a = 2.42808/(pow(k,2)+pow(.840297,2))+3.42835;
			b = .0167941/(pow(k,2)+pow(.47573,2))+1.70158;
			omega0 = sqrt(pow(1.42617+Shift,2)+pow(k,2))+.258289;
			knee = 3.59947*pow(k+1.,(long double)-.710425);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .592982/sqrt(pow(k,2)+pow(2.06656,2));
			a = 2.09377/(pow(k,2)+pow(.763871,2))+2.65712;
			b = .366499/(pow(k,2)+pow(1.06864,2))+1.35141;
			omega0 = sqrt(pow(1.38555+Shift,2)+pow(k,2))+.253076;
			knee = 3.49204*pow(k+1.,(long double)-.925502);
			break;
		default:
			Sigma = 0;
			a = 1;
			b = 1;
			omega0 = 1;
			knee = 1;
			break;
	}

	long double ImSigma;	//Calculation of the argument to the exponential, these first 2 are approximations to hopefully avoid catastrophic loss of precision
	if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee < -4.)
		ImSigma = a*(omega-omega0+knee/sqrt(a*b));
	else if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee > 4.)
		ImSigma = b*(omega0-omega+knee/sqrt(a*b));
	else
		ImSigma = -.5*((a-b)*omega0-((a+b)*knee)/sqrt(a*b))+(a-b)*omega/2-sqrt(pow(((a+b)/2.)*(omega-omega0+((a-b)*knee)/(sqrt(a*b)*(a+b))),2)+pow(knee,2));

	if(pow(omega,2)>=pow(k,2))
		return(-2.*M*Sigma*exp(ImSigma)+sqrt(pow(omega,2)-pow(k,2))*GAMMA);
	else
		return(-2.*M*Sigma*exp(ImSigma));
#endif
}

long double ReSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
#ifdef RIEK
	long double Sigma;
	long double M1, M2;
	long double M_T;
	long double gamma;
	long double a;
	long double sigma1, sigma2;
	long double Floor;
	
	switch(Temp)
	{
		case 0:
			return(0);
			break;
		case 1:
			/*Sigma = .040845;
			M1 = 1.67925;
			M2 = 1.65299;
			gamma = .337277;*/
			Sigma = .135862;
			M1 = 1.6734;
			M2 = 1.65605;
			gamma = .435771/sqrt(pow(k,2)+pow(1.45558,2))+.000075623;
			M_T = 1.848;
			a = 2.2094;
			sigma1 = 2.574538196654789;
			sigma2 = 2.5745381961599816;
			Floor = .928942;
			break;
		case 2:
			/*Sigma = .0364537;
			M1 = 1.58845;
			M2 = 1.55012;
			gamma = .394468;*/
			Sigma = .100932;
			M1 = 1.58441;
			M2 = 1.55442;
			M_T = 1.719;
			gamma = .508717/sqrt(pow(k,2)+pow(1.42858,2))+.0000832943;
			a = .760248;
			sigma1 = 3.41189;
			sigma2 = 1.38003;
			Floor = .818425;
			break;
		case 3:
			/*Sigma = .050848;
			M1 = 1.54629;
			M2 = 1.47387;
			gamma = .517924;*/
			Sigma = .104446;
			M1 = 1.53819;
			M2 = 1.4777;
			M_T = 1.563;
			gamma = .627196/sqrt(pow(k,2)+pow(1.3015,2))+.000261014;
			a = .752122;
			sigma1 = 3.38156;
			sigma2 = 1.30001;
			Floor = .905761;
			break;
	}
	Floor = 0;
	long double Shift = M-M_T;
	M1 += Shift;
	M2 += Shift;

	return(Sigma*gamma*(omega-sqrt(pow(M1,2)+pow(k,2)))/(pow(omega-sqrt(pow(M2,2)+pow(k,2)),2)+pow(gamma,2))*(a*exp(-pow(k/sigma1,2))+(1.-a)*exp(-pow(k/sigma2,2))));
#endif
#if defined(SHUAI) || defined(CC) || defined(BB)
	long double Sigma;	//Strength
	long double x0, x1;	//Centrality markers
	long double gamma;	//Width
	long double Shift, M_T;

	switch(Temp)
	{
		case 0:
			return(0);
			break;
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .244368/sqrt(pow(k,2)+pow(1.32368,2));
			x0 = sqrt(pow(k,2)+pow(1.5567+Shift,2))+.279476;
			x1 = sqrt(pow(k,2)+pow(1.50202+Shift,2))+.259;
			gamma = .320676/sqrt(pow(k,2)+pow(1.56455,2))+.080032;
			break;
		case 2://258MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .322887/sqrt(pow(k,2)+pow(1.34236,2));
			x0 = sqrt(pow(k,2)+pow(1.54159+Shift,2))+.280535;
			x1 = sqrt(pow(k,2)+pow(1.46598+Shift,2))+.260561;
			gamma = .694901/sqrt(pow(k,2)+pow(2.13185,2))+.0653795;
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .375163/sqrt(pow(k,2)+pow(1.41612,2));
			x0 = sqrt(pow(k,2)+pow(1.45507+Shift,2))+.337448;
			x1 = sqrt(pow(k,2)+pow(1.40846+Shift,2))+.289292;
			gamma = .690491/sqrt(pow(k,2)+pow(1.97525,2))+.141465;
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .370623/sqrt(pow(k,2)+pow(1.53585,2));
			x0 = sqrt(pow(k,2)+pow(1.39619+Shift,2))+.35548;
			x1 = sqrt(pow(k,2)+pow(1.3481+Shift,2))+.296587;
			gamma = .857781/sqrt(pow(k,2)+pow(2.25072,2))+.196022;
			break;
		default:
			Sigma = 0;
			x0 = 1;
			x1 = 1;
			gamma = 1;
			break;
	}

	return(Sigma*(omega-x0)/(pow(omega-x1,2)+gamma));
#endif
}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	return(sqrt(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta)));
}

long double Fermi(long double omega, int T)	//Fermi factor
{
	long double Temp;

	switch(T)
	{
		case 0:
			return(0);
			break;
		case 1:
			Temp = .2352;
			break;
		case 2:
			Temp = .294;
			break;
		case 3:
			Temp = .392;
			break;
		default:
			return(0);
	}

	return(1./(1.+exp(omega/Temp)));
}

long double Potential_on(long double Par[5])	//On-shell potential for the on-shell T-Matrix
{
	return(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+abs(Par[4]-4.*pow(Par[2],2))),2));
}

long double Potential1(long double Par[5], long double omega, long double k)	//Potiential for the numerator of the boson spectrum
{
	return(pow(Par[1],2)/(pow(Par[1],2)+abs(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2))));
}

long double Potential2(long double Par[5], long double omega, long double k)	//Potiential for the denominator of the T-Matrix and boson spectrum
{
	return(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+abs(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2))),2));
}

long double Quark_Spectrum(long double omega, long double k, long double M, int Temp)	//Single quark spectral function
{
	return(ImSelf_Energy(M, omega, k, Temp)/(pow(pow(omega,2)-pow(k,2)-pow(M,2)-2.*M*ReSelf_Energy(M, omega, k, Temp),2)+pow(ImSelf_Energy(M, omega, k, Temp),2)));
}

long double Spin_Sum(long double Par[5], long double omega, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
{
	return((-pow(Par[3],2)/4.+pow(k,2)+Par[3]*k*cos(theta)+2.*pow(Par[2],2)+omega*(sqrt(Par[4]+pow(Par[3],2))-omega))/pow(Par[2],2));
}

long double Folding_Integrand(long double Par[5], long double omega, long double k, long double theta, int Temp)	//Integrand of the folding integral
{
	return(-pow(Par[2],2)*Quark_Spectrum(omega, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*Quark_Spectrum(sqrt(Par[4]+pow(Par[3],2))-omega, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(1.-Fermi(omega, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))-omega, Temp)));
}
