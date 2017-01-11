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
long double Binary_Pole(long double[5]);	//Binary search for a theta pole
void Characterize_k_Int(long double[5], int, long double, long double[4], long double[4], int&);	//Returns the poles of the k integral's integrands
long double Width(long double[4], long double[4], int, int[4]);	//Returns the width that is approiate for the V/on-shell integral
void Newtons_k_Int(long double, long double, long double, long double, long double[3], long double);	//Newton's method of finding roots for Characterize_k_Int
int Newtons_Test_k_Int(long double, long double, long double, long double, long double[3], long double);	//Trys to determine the number of roots Newton's method  is looking for and where they might be found
void Characterize_Folding(long double[5], int, long double, long double, long double[4], long double[4], long double[2], int&, int[4]);	//Returns the poles of the folding integral's integrands

//Straight Functions everything is built from
long double ImSelf_Energy(long double, long double, long double, int);	//Imaginary single quark self energy
long double ReSelf_Energy(long double, long double, long double, int);	//Real single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Potential_on(long double[5]);	//On-shell potential for the on-shell T-Matrix
long double Potential1(long double[5], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[5], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double Quark_Spectrum(long double, long double, long double, int);	//Single quark spectral function
long double Spin_Sum(long double[5], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Folding_Integrand(long double[5], long double, long double, long double, int);	//Integrand of the folding integral

#define GAMMA -.015
long double Boundary[] = {.07,.25,.4,.65,.9,1.25,1.5,2,2.5,.06,.67,1,1.5,3.1,4.5,5.25,6,.07,.15,.2,.7,1,1.5};

//long double Par[5] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[5], int Temp)	//Integrates the theta results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Range[] = {-M_PI,Boundary[17],Boundary[18],Boundary[19],Boundary[20],Boundary[21],Boundary[22],M_PI/2.};
	long double x1, x2;	//Abscissa
	Elements F;	//Sum of ordinate*weights
	Elements Answer(0,0,0);	//Answer to be returned
	long double a = 0, b;	//Sub-interval limits of integration
	int i, j;	//Counters
	Elements holder;
	//ofstream Table("theta Table", ios::app);
	//Table << setprecision(18);

	if(false)//Range[0] > 0)	//Found a pole in theta
	{
		for(i = 0; i < 7; i++)	//Bubble sort the single element to where it belongs
		{
			if(Range[i] < Range[i+1])
			{
				long double Holder = Range[i];
				Range[i] = Range[i+1];
				Range[i+1] = Holder;
			}
			else
				break;
		}

		for(i = 0; i < 7; i++)
		{
			b = Range[i];

			F.null();
			for(j = 0; j < 9; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = k_Int(Par, Temp, x1);
				F += holder*sin(x1)*w[j+1];
				//Table << Par[3] << " " << Par[4] << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
				holder = k_Int(Par, Temp, x2);
				F += holder*sin(x2)*w[j+1];
				//Table << Par[3] << " " << Par[4] << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			}
			holder = k_Int(Par, Temp, (a+b)/2.);
			F += holder*sin((a+b)/2.)*w[0];
			//Table << Par[3] << " " << Par[4] << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			Answer += F*(b-a)/2.;
			a = b;
		}
	}

	for(i = 1; i < 7; i++)	//No pole found
	{
		b = Range[i];

		F.null();
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			holder = k_Int(Par, Temp, x1);
			F += holder*sin(x1)*w[j+1];
			//Table << Par[3] << " " << Par[4] << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			holder = k_Int(Par, Temp, x2);
			F += holder*sin(x2)*w[j+1];
			//Table << Par[3] << " " << Par[4] << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		}
		holder = k_Int(Par, Temp, (a+b)/2.);
		F += holder*sin((a+b)/2.)*w[0];
		//Table << Par[3] << " " << Par[4] << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		Answer += F*(b-a)/2.;
		a = b;
	}

	return(Answer/pow(2.*M_PI,2)*2.);
}

long double Binary_Pole(long double Par[5])
{
	long double zero[3];
	long double theta[3] = {0,M_PI/4.,M_PI/2.};
	int Num[3];

	zero[2] = zero[1] = zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta[2]),2)));
	Num[2] = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta[2]);
	if(Num[2] != 0)
		return(-M_PI);	//Nothing found
	zero[1] = zero[1] = zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta[1]),2)));
	Num[1] = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta[1]);
	zero[0] = zero[1] = zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta[0]),2)));
	Num[0] = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta[0]);

	while(theta[2]-theta[0] > M_PI/2000.)	//Binary search to accuracy of M_PI/2000=.00157... (No more than 10 subdivisions)
	{
		if(Num[1] == 0)
		{
			Num[2] = Num[1];
			theta[2] = theta[1];
			theta[1] = (theta[2]+theta[0])/2.;
			zero[1] = zero[1] = zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta[1]),2)));
			Num[1] = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta[1]);
		}
		else
		{
			Num[0] = Num[1];
			theta[0] = theta[1];
			theta[1] = (theta[2]+theta[0])/2.;
			zero[1] = zero[1] = zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta[1]),2)));
			Num[1] = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta[1]);
		}
	}

	return(theta[1]);
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements k_Int(long double Par[5], int Temp, long double theta)	//Integrates the k momentum results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Range[] = {-Boundary[16], -Boundary[15], -Boundary[14], -Boundary[13], -Boundary[12], -Boundary[11], -Boundary[10], -Boundary[9], 0, Boundary[9], Boundary[10], Boundary[11], Boundary[12], Boundary[13], Boundary[14], Boundary[15], Boundary[16]};	//Number of gamma from center
	Elements F_a, F_b, F_ave;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Answer to be returned
	Elements PartialAnswer;	//Answer for sub-interval for determining completeness
	long double x1, x2;	//Abscissa
	long double a = 0, b = 0;//Sub-interval limits of integration
	int Poles;	//Number of poles
	long double zero[5];	//The real part of the signular pole
	long double gamma[5];	//The distance to the singular, maybe
	long double Min_upper;	//The integral has to go at least this far
	long double Width;	//Length of the next sub-interval
	long double Early = 0;	//Early change from one pole to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;//The next width that will be used in the event of an early change of poles
	int i = 0, j, l;	//Counters
	Elements holder;
	//ofstream Table("k Table", ios::app);
	//ofstream Pole_Tab("k Poles", ios::app);
	//Table << setprecision(18);
	//Pole_Tab << setprecision(18);

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);
	//for(i = 0; i < Poles; i++)
		//Pole_Tab << Par[3] << " " << Par[4] << " " << theta << " " << zero[i] << " " << gamma[i] << endl;

	Min_upper = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//This the upper bound that the vacuum calls for, Partial/total will promote higher as needed

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0)	//j!=0 is because the loop will find other early termination
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

		if(j == 16 && i < Poles)	//Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i]))	//Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Poles != 0 && j < 16 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && (b != 0 || Width == 0))
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
			j = 16-j;
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

			holder = Folding(Par, Temp, x1, theta);
			F_a += holder*pow(x1,2)*w[l+1]; //Evaluate function at x1
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			holder = Folding(Par, Temp, x2, theta);
			F_b += holder*pow(x2,2)*w[l+1]; //Evaluate function at x2
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		}
		holder = Folding(Par, Temp, (a+b)/2., theta);
		F_ave = holder*pow((a+b)/2.,2)*w[0]; //Evaluate function at (a+b)/2.
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		PartialAnswer = (F_a+F_ave+F_b)*(b-a)/(2.);
		Answer += PartialAnswer;
		a = b;
	}while(b < Min_upper/* || abs(PartialAnswer/Answer) >= .0001*/);

	return(Answer);
}

void Characterize_k_Int(long double Par[5], int Temp, long double theta, long double zero[5], long double gamma[5], int &Poles) //Returns the poles of the k integral's integrands
{
	long double holder;
	long double previous[2];
	int i, j = 0, l;
	long double Nothing[2];
	long double gamma_fetch[4];
	long double zero_fetch[4];
	int Type[4];

	if(sqrt(Par[4]) > 2.*Par[2])    //Find the double on-shell pole and estimate a width (distance to pole)
	{
		zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
		gamma[2] = zero[2]*.0021/13.38036;
		Poles = 3;
	}
	else
	{
		zero[2] = 0;
		Poles = 0;
	}

	zero[0] = zero[1] = zero[2];	//Start the Newton's search for potiential/on-shell pole
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

			Characterize_Folding(Par, Temp, zero[0], theta, zero_fetch, gamma_fetch, Nothing, l, Type);	//Get the omega-zero data at k[0]
			gamma[0] = Width(zero_fetch, gamma_fetch, l, Type);	//Use the incoming data to figure out the correct width to apply
			Characterize_Folding(Par, Temp, zero[1], theta, zero_fetch, gamma_fetch, Nothing, l, Type);	//Get the omega-zero data at k[1]
			gamma[1] = Width(zero_fetch, gamma_fetch, l, Type);	//Use the incoming data to figure out the correct width to apply
			Poles = 3;
		}
		else
		{
			zero[0] = zero[2];
			gamma[0] = zero[2]*.0021/13.38036;//Folding(Par, Temp, zero[0], theta).Min();
			/*if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
				gamma[0] = 1e-3;
			else if(gamma[0] > abs(2.*Par[2]*GAMMA))
				gamma[0] = abs(2.*Par[2]*GAMMA);*/
			Poles = 1;
		}
	}
	else
	{
		zero[0] = zero[2];
		gamma[0] = zero[2]*.0021/13.38036;//Folding(Par, Temp, zero[0], theta).Min();
		/*if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[0] = 1e-3;
		else if(gamma[0] > abs(2.*Par[2]*GAMMA))
			gamma[0] = abs(2.*Par[2]*GAMMA);*/
		Poles = 1;
	}

	if(pow(Par[2],2)-.5*pow(Par[3]*sin(theta),2)+sqrt(pow(Par[3]*sin(theta),2)/8.*(8.*pow(Par[2],2)+2.*pow(Par[3]*sin(theta),2))) <= Par[4])
	{//On-shell departure from the time-like region
		if(Par[3] == 0)
		{
			zero[Poles] = (Par[4]-pow(Par[2],2))/(sqrt(4.*Par[4]));
			gamma[Poles] = abs(2.*Par[2]*GAMMA);
			Poles++;
		}
		else
		{
			zero[Poles] = abs(pow(Par[2],2)*Par[3]*cos(theta)+sqrt((pow(Par[3],2)+Par[4])*(pow(pow(Par[2],2)-Par[4],2)+(Par[4]-2.*pow(Par[2],2))*pow(Par[3]*sin(theta),2))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2)));
			zero[Poles+1] = abs(pow(Par[2],2)*Par[3]*cos(theta)-sqrt((pow(Par[3],2)+Par[4])*(pow(pow(Par[2],2)-Par[4],2)+(Par[4]-2.*pow(Par[2],2))*pow(Par[3]*sin(theta),2))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2)));
			gamma[Poles] = abs(2.*Par[2]*GAMMA);
			gamma[Poles+1] = abs(2.*Par[2]*GAMMA);
			Poles += 2;
		}
	}

	for(i = Poles-1; i >= 0; i--)	//Bubble sort
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
			else if(zero[j] == zero[j+1])	//Remove duplicates
			{
				for(l = j; l < Poles-1; l++)
				{
					zero[l] = zero[l+1];
					gamma[l] = gamma[l+1];
				}
				zero[2] = sqrt(Par[4]+pow(Par[3],2)); //Need to use the biggest finite value of reason or it will be attempted to be sorted to the bottom when it is invalid
				gamma[2] = 0;
				Poles--;
			}
		}
	}

	/*for(i = 0; i < Poles; i++)
	{
		gamma[i] = Folding(Par, Temp, zero[i], theta).Min();
		if(gamma[i] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[i] = 1e-3;
		if(gamma[i] > abs(2.*Par[2]*GAMMA))
			gamma[i] = abs(2.*Par[2]*GAMMA);
	}*/
	return;
}

long double Width(long double zero[4], long double gamma[4], int Poles, int Type[4])
{
	long double holder;
	int i,j;
	int Sum = 0;

	for(i = Poles-1; i >= 0; i--)	//Bubble sort, resort based on Type instead of zero
	{
		for(j = 0; j < i; j++)
		{
			if(Type[j] > Type[j+1])
			{
				holder = zero[j+1];
				zero[j+1] = zero[j];
				zero[j] = holder;
				holder = gamma[j+1];
				gamma[j+1] = gamma[j];
				gamma[j] = holder;
				holder = Type[j+1];
				Type[j+1] = Type[j];
				Type[j] = holder;
			}
		}
	}

	for(i = 0; i < Poles; i++)	//Sum the pole type by cross category
	{
		switch(Type[i])
		{
			case 0:
			case 3:
				Sum += 1;
				break;
			case 1:
			case 2:
				break;
		}
	}

	if(Poles == 4)	//Should only get here if there is an actual zeroes in the omega at k
	{
		if(abs(zero[0]-zero[3]) < abs(zero[1]-zero[2]))	//V+/omega-
			return(gamma[0]+gamma[3]);
		else	//V-/omega+
			return(gamma[1]+gamma[2]);
	}
	else if(Sum == 2)	//There is only 3 poles in either of the next 2 cases, V+/omega-
		return(gamma[0]+gamma[2]);	//If you ever need to debug this, look at pg 148 in 2016 notes
	else	//The Sum can only be 1, V-/omega+
		return(gamma[0]+gamma[1]);

	return(0);
}

void Newtons_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = -sqrt(2.)*sqrt(s+pow(P,2))+pow(16.*pow(k[0],4)+pow(Lambda,4),(long double).25)*sqrt(1.+4./sqrt(16.+pow(Lambda/k[0],4)))+2.*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[0],theta),4)))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[0],theta),4),(long double).25);
	long double fp1 = (4.*pow(k[0],3)*(4.+sqrt(16.+pow(Lambda/k[0],4))))/(pow(16.*pow(k[0],4)+pow(Lambda,4),(long double).75)*sqrt((pow(Lambda,4)+4.*pow(k[0],4)*(4.+sqrt(16.+pow(Lambda/k[0],4))))/(16.*pow(k[0],4)+pow(Lambda,4))))+((2.*k[0]+P*cos(theta))*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[0],theta),4)))*pow(Energy(M,P/2.,k[0],theta),2))/pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[0],theta),4),(long double).75)+(pow(M*GAMMA,2)*(2.*k[0]+P*cos(theta))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[0],theta),4),(long double).25))/(sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[0],theta),4)))*pow(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[0],theta),4),(long double)1.5)*pow(Energy(M,P/2.,k[0],theta),6));
	long double f2 = -sqrt(2.)*sqrt(s+pow(P,2))+pow(16.*pow(k[1],4)+pow(Lambda,4),(long double).25)*sqrt(1.+4./sqrt(16.+pow(Lambda/k[1],4)))+2.*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[1],theta),4)))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[1],theta),4),(long double).25);
	long double fp2 = (4.*pow(k[1],3)*(4.+sqrt(16.+pow(Lambda/k[1],4))))/(pow(16.*pow(k[1],4)+pow(Lambda,4),(long double).75)*sqrt((pow(Lambda,4)+4.*pow(k[1],4)*(4.+sqrt(16.+pow(Lambda/k[1],4))))/(16.*pow(k[1],4)+pow(Lambda,4))))+((2.*k[1]-P*cos(theta))*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[1],theta),4)))*pow(Energy(M,P/2.,-k[1],theta),2))/pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[1],theta),4),(long double).75)+(pow(M*GAMMA,2)*(2.*k[1]-P*cos(theta))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[1],theta),4),(long double).25))/(sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[1],theta),4)))*pow(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[1],theta),4),(long double)1.5)*pow(Energy(M,P/2.,-k[1],theta),6));
	k[0] -= f1/fp1;
	k[1] -= f2/fp2;
	return;
}

int Newtons_Test_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = -sqrt(2.)*sqrt(s+pow(P,2))+pow(16.*pow(k[2],4)+pow(Lambda,4),(long double).25)*sqrt(1.+4./sqrt(16.+pow(Lambda/k[2],4)))+2.*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[2],theta),4)))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[2],theta),4),(long double).25);
	long double fp1 = (4.*pow(k[2],3)*(4.+sqrt(16.+pow(Lambda/k[2],4))))/(pow(16.*pow(k[2],4)+pow(Lambda,4),(long double).75)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4))))/(16.*pow(k[2],4)+pow(Lambda,4))))+((2.*k[2]+P*cos(theta))*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[2],theta),4)))*pow(Energy(M,P/2.,k[2],theta),2))/pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[2],theta),4),(long double).75)+(pow(M*GAMMA,2)*(2.*k[2]+P*cos(theta))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,k[2],theta),4),(long double).25))/(sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[2],theta),4)))*pow(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,k[2],theta),4),(long double)1.5)*pow(Energy(M,P/2.,k[2],theta),6));
	long double limitp1 = (P*((16.*pow(M*GAMMA,2))/sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2))+pow(4.*pow(M,2)+pow(P,2),2)*(1.+1./sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2))))*cos(theta))/(2.*(4.*pow(M,2)+pow(P,2))*pow(16.*pow(M,4)+pow(P,4)+8.*pow(M,2)*(pow(P,2)+2.*pow(GAMMA,2)),(long double).75)*sqrt(1.+1./sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2))));
	long double f2 = -sqrt(2.)*sqrt(s+pow(P,2))+pow(16.*pow(k[2],4)+pow(Lambda,4),(long double).25)*sqrt(1.+4./sqrt(16.+pow(Lambda/k[2],4)))+2.*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[2],theta),4)))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[2],theta),4),(long double).25);
	long double fp2 = (4.*pow(k[2],3)*(4.+sqrt(16.+pow(Lambda/k[2],4))))/(pow(16.*pow(k[2],4)+pow(Lambda,4),(long double).75)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4))))/(16.*pow(k[2],4)+pow(Lambda,4))))+((2.*k[2]-P*cos(theta))*sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[2],theta),4)))*pow(Energy(M,P/2.,-k[2],theta),2))/pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[2],theta),4),(long double).75)+(pow(M*GAMMA,2)*(2.*k[2]-P*cos(theta))*pow(pow(M*GAMMA,2)+pow(Energy(M,P/2.,-k[2],theta),4),(long double).25))/(sqrt(1.+1./sqrt(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[2],theta),4)))*pow(1.+(pow(M*GAMMA,2))/pow(Energy(M,P/2.,-k[2],theta),4),(long double)1.5)*pow(Energy(M,P/2.,-k[2],theta),6));
	long double limitp2 = (P*(-((16.*pow(M*GAMMA,2))/sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2)))-pow(4.*pow(M,2)+pow(P,2),2)*(1.+1./sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2))))*cos(theta))/(2.*(4.*pow(M,2)+pow(P,2))*pow(16.*pow(M,4)+pow(P,4)+8.*pow(M,2)*(pow(P,2)+2.*pow(GAMMA,2)),(long double).75)*sqrt(1.+1./sqrt(1.+(16.*pow(M*GAMMA,2))/pow(4.*pow(M,2)+pow(P,2),2))));
	long double limit = -sqrt(2)*sqrt(s+pow(P,2))+pow(16.*pow(M,4)+pow(P,4)+8.*pow(M,2)*(pow(P,2)+2.*pow(GAMMA,2)),(long double).25)*sqrt(1.+1./sqrt(1.+16.*pow(M*GAMMA,2)/pow(4.*pow(M,2)+pow(P,2),2)))+Lambda;
	int roots = 0;

	if(limit*f1 < 0)
		roots++;
	if(limit*f2 < 0)
		roots++;
	if(limitp1*fp1 < 0 && limitp1*limit < 0)
		roots+=2;
	if(limitp2*fp2 < 0 && limitp2*limit < 0)
		roots+=2;

	if(roots > 2)
		cerr << "Unexepected outcome in Newtons_Test_k_Int s = " << s << " P = " << P << " theta = " << theta << endl;

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
	if(/*Temp == 0 && */abs(sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta)-Energy(0,Par[3]/2.,k,theta)) < 1e-12)	//Let's save some time and just return 0, because it is
		return(Elements(0,0,0));
	else if(Par[4]+pow(Par[3],2) <= 0)
		return(Elements(0,0,0));	//Bad data trap and time saver

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204}; //Weight of the function at Disp
	long double Range[] = {-Boundary[8], -Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7], Boundary[8]};	//Number of gamma from center
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
	long double End_Points[2];	//The end points of the omega integral
	int Nothing[4];		//Something to send to Characterize_Folding
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	Elements holder;
	//ofstream Table("omega Table", ios::app);
	//ofstream Pole_Tab("omega Poles", ios::app);
	//Table << setprecision(18);
	//Pole_Tab << setprecision(18);

	Characterize_Folding(Par, Temp, k, theta, zero, gamma, End_Points, Poles, Nothing);	//Get the poles that I have to be concerned about
	//for(i = 0; i < Poles; i++)
		//Pole_Tab << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << zero[i] << " " << gamma[i] << endl;

	/*a = b = End_Points[0];
	Max = End_Points[1];//*/
	if(false)//Temp != 0) //This may or may not be a permanent change. Due to this uncertainity, I'm leaving it here.
	{
		a = b = 0;
		Max = sqrt(Par[4]+pow(Par[3],2));
		cout << "If you've come down this way,you must first correct issues of pole intersections between potiential and quark propagators. I have attempted to fuck everything up if you miss this warning." << endl;
		return(Elements(0./0.,0./0.,0./0.));
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta);
		Max = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}//*/

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0) //j!=0 is because the loop will find other early termination
		Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];

	do
	{
		if((b == 0 || b == Energy(0,Par[3]/2.,k,theta)) && j != 0 && zero[i]+Range[j]*gamma[i] != b)	//First pole is closer than zero-64*gamma to lower limit of integration
			Width = zero[i]+Range[j]*gamma[i]-b;
		else if((i < Poles && b+100 < zero[i]-64.*gamma[i]) || b+100 < Max)	//Middle of nowhere intervals
			Width = 100;
		else if((i < Poles && b+50 < zero[i]-64.*gamma[i]) || b+50 < Max)
			Width = 50;
		else if((i < Poles && b+10 < zero[i]-64.*gamma[i]) || b+10 < Max)
			Width = 10;
		else
			Width = 3;

		if(j == 18 && i < Poles)	//Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i]))	//Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1])	//There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];	//Sets early terminantion point
		}
		else if(Poles != 0 && j < 18 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && !(b == 0 || b == Energy(0,Par[3]/2.,k,theta)))
		{//Integrating a pole and the width wan't resolved by the first condition after do, 24 lines up
			Width = gamma[i]*(Range[j+1]-Range[j]);
			j++;
		}

		if(NextWidth != 0 && NextWidth == NextWidth)	//Resolving early termination of integrating one pole for the next
		{
			Width = NextWidth;
			NextWidth = 0;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1])	//There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];	//Sets early terminantion point
		}
		else if(Early != 0 && b+Width > Early)
		{
			Width = Early-b;
			Early = 0;
			j = 18-j;
			NextWidth = Width*gamma[i+1]/gamma[i];
			i++;
		}

		b += Width;	//Set next b keeping in mind the upper limit of integration
		if(b > Max)
			b = Max;

		F_a.null();
		F_b.null();
		for(l = 0; l < 9; l++)	//Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			holder = Elements(Spin_Sum(Par, x1, k, theta), 2.*Potential1(Par,x1,k), Potential2(Par,x1,k))*Folding_Integrand(Par,x1,k,theta,Temp);
			F_a += holder*w[l+1];
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			holder = Elements(Spin_Sum(Par, x2, k, theta), 2.*Potential1(Par,x2,k), Potential2(Par,x2,k))*Folding_Integrand(Par,x2,k,theta,Temp);
			F_b += holder*w[l+1];
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		}
		holder = Elements(Spin_Sum(Par, (a+b)/2., k, theta), 2.*Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*Folding_Integrand(Par,(a+b)/2.,k,theta,Temp);
		F_ave = holder*w[0];
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;

		a = b;
	}while(b < Max);

	return(Answer/M_PI);
}

void Characterize_Folding(long double Par[5], int Temp, long double k, long double theta, long double zero[4], long double gamma[4], long double End_Points[2], int &Poles, int Type[4])
{
	long double Lower, Upper;	//Limits of integration in Folding, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(false)//Temp != 0)
	{
		Lower = 0;
		Upper = sqrt(Par[4]+pow(Par[3],2));
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta);
		Upper = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	/*long double List[int(floor((Upper-Lower)/10.))+2];
	List[0] = abs(Folding_Integrand(Par,Lower,k,theta,Temp));
	for(i = 1; i < floor((Upper-Lower)/10.)+1; i++)
		List[i] = abs(Folding_Integrand(Par,Lower+(i-1)*10.,k,theta,Temp));
	List[i+1] = abs(Folding_Integrand(Par,Upper,k,theta,Temp));
	j = 0;
	for(i = 1; i < floor((Upper-Lower)/10.)+2; i++)
		if(List[j] > List[i])
			j = i;
	for(i = 0; i < floor((Upper-Lower)/10.)+1; i++)
	{
		if(List[i] < List[j]/1.e4 && List[i+1] >= List[j]/1.e4)
			if(i != 0)
				Lower += i*10.;
		else if(List[i] > List[j]/1.e4 && List[i+1] <= List[j]/1.e4)
		{
			Upper = Lower+i*10.+10.;
			break;
		}
	}
	End_Points[0] = Lower;
	End_Points[1] = Upper;//*/

	zero[0] = .5*(sqrt(Par[4]+pow(Par[3],2))+pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));	//Potential poles, I know exactly where these are at.
	zero[1] = .5*(sqrt(Par[4]+pow(Par[3],2))-pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[0] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));
	gamma[1] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));
	Type[0] = 0;	//Labels which pole goes with the other at V/on-shell intersections
	Type[1] = 1;

	holder = GAMMA;
	zero[2] = pow(pow(Energy(Par[2],Par[3]/2.,k,theta),4)+pow(Par[2]*holder,2),(long double).25)*cos(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,k,theta),2)));	//Exact vacuum
	zero[3] = pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),4)+pow(Par[2]*holder,2),(long double).25)*cos(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,-k,theta),2)));

	if(Temp != 0)	//media estimate
	{
		holder = ImSelf_Energy(Par[2], zero[2], Energy(0,Par[3]/2.,k,theta), Temp)/sqrt(pow(zero[2],2)-pow(Energy(0,Par[3]/2.,k,theta),2));
		zero[2] = pow(pow(Energy(Par[2],Par[3]/2.,k,theta),4)+pow(Par[2]*holder,2),(long double).25)*cos(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,k,theta),2)));
		gamma[2] = abs(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),4)+pow(Par[2]*holder,2),(long double).25)*sin(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,k,theta),2))));

		holder = ImSelf_Energy(Par[2], zero[3], Energy(0,Par[3]/2.,-k,theta), Temp)/sqrt(pow(zero[3],2)-pow(Energy(0,Par[3]/2.,-k,theta),2));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),4)+pow(Par[2]*holder,2),(long double).25)*cos(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,-k,theta),2)));
		gamma[3] = abs(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),4)+pow(Par[2]*holder,2),(long double).25)*sin(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,-k,theta),2))));
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = abs(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),4)+pow(Par[2]*holder,2),(long double).25)*sin(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,k,theta),2))));
		gamma[3] = abs(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),4)+pow(Par[2]*holder,2),(long double).25)*sin(.5*atan2(Par[2]*holder,pow(Energy(Par[2],Par[3]/2.,-k,theta),2))));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-zero[3];
	}
	Type[2] = 2;
	Type[3] = 3;


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
				holder = Type[j+1];
				Type[j+1] = Type[j];
				Type[j] = holder;
			}
		}
	}

	i = j = 0;	//Find the first zero greater than 0
	while(zero[i] < Lower) i++;

	while(zero[i] <= Upper && i < 4)	//Move zeroes up to front of array, count off poles within the limits of integration
	{
		zero[j] = zero[i];
		gamma[j] = gamma[i];
		Type[j] = Type[i];
		i++;
		j++;
	}
	Poles = j;

	return;
}

//long double Par[5] = {g, Lambda, M, P, s}
long double ImSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Par[6];	//Momentum dependance parameterization
	long double E_0 = Energy(M,k,0,0);	//location of lorentzian
	long double Sigma;	//size of energy dependance
	long double b1, b2;	//slope of exponential decrease to left and right
	long double Delta;	//concavity or length of transition from left to right

	if(omega < k)
		return(0);

	switch(Temp)
	{
		case 0:
			if(omega>=k)
				return(M*GAMMA);
			else
				return(0);
			break;
		case 1://235.2MeV
			Sigma = -0.1289680072519721;
			b1 = 3.322825262881092;
			b2 = 2.2878310836782014;
			Delta = 1.228601982782018;
			Par[0] = .7359389831810698;
			Par[1] = 7.487501146014314;
			Par[2] = 1.9490238595657456;
			Par[3] = .700215754;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 2://294MeV
			Sigma = -0.09606152620146369;
			b1 = 3.285053276019642;
			b2 = 1.886285913340202;
			Delta = 1.1858269101609233;
			Par[0] = .7409390219065235;
			Par[1] = 7.450458343071824;
			Par[2] = 1.8620618988580635;
			Par[3] = .860810762;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 3://362MeV
			Sigma = -0.09933548776506283;
			b1 = 3.2108770392083246;
			b2 = 1.3694064180118886;
			Delta = 1.3043774341616825;
			Par[0] = .7426375963204489;
			Par[1] = 7.698646415632565;
			Par[2] = 1.771465704769189;
			Par[3] = .608717852;
			Par[4] = 10;
			Par[5] = 3;
			break;
		default:
			Sigma = 0;
			b1 = 0;
			b2 = 0;
			Delta = 0;
			Par[0] = 1;
			Par[1] = 1;
			Par[2] = 1;
			Par[3] = 0;
			Par[4] = 0;
			Par[5] = 0;
			break;
	}
	Par[3] = 0;

	return(2.*M*(Par[0]*exp(-pow(k/Par[1],2))+(1-Par[0])*exp(-pow(k/Par[2],2))+Par[3])*(Sigma*exp(Delta+(b1-b2)*(omega-E_0)*E_0/2.-sqrt(b1*b2*pow((omega-E_0)*E_0,2)+pow(Delta+(b1-b2)*(omega-E_0)*E_0/2.,2))))+M*GAMMA);
}

long double ReSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Par[6];	//Momentum dependance parameterization
	long double E_0 = Energy(M,k,0,0);	//location of lorentzian
	long double a;	//size of energy dependance
	long double gamma;	//width of lorentzian
	long double c;	//zero crossing, might be better as the on-shell energy

	if(omega < k)
		return(0);

	switch(Temp)
	{
		case 0:
			return(0);
			break;
		case 1://235.2MeV
			a = .0412729;
			c = 1.68597;
			gamma = .340028;
			Par[0] = .7359389831810698;
			Par[1] = 7.487501146014314;
			Par[2] = 1.9490238595657456;
			Par[3] = .700215754;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 2://294MeV
			a = .0366063;
			c = 1.5945;
			gamma = .39357;
			Par[0] = .7409390219065235;
			Par[1] = 7.450458343071824;
			Par[2] = 1.8620618988580635;
			Par[3] = .860810762;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 3://362MeV
			a = .05127;
			c = 1.55455;
			gamma = .518487;
			Par[0] = .7426375963204489;
			Par[1] = 7.698646415632565;
			Par[2] = 1.771465704769189;
			Par[3] = .608717852;
			Par[4] = 10;
			Par[5] = 3;
			break;
		default:
			a = 0;
			c = 0;
			gamma = 1;
			Par[0] = 1;
			Par[1] = 1;
			Par[2] = 1;
			Par[3] = 0;
			Par[4] = 0;
			Par[5] = 0;
			break;
	}
	Par[3] = 0;

	return(2.*M*(Par[0]*exp(-pow(k/Par[1],2))+(1-Par[0])*exp(-pow(k/Par[2],2))+Par[3])*(a*(omega-c)/(pow(omega-c,2)+pow(gamma,2))));
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
	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(Par[4]-4.*pow(Par[2],2),2)),2));
}

long double Potential1(long double Par[5], long double omega, long double k)	//Potiential for the numerator of the boson spectrum
{
	return(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2),2)));
}

long double Potential2(long double Par[5], long double omega, long double k)	//Potiential for the denominator of the T-Matrix and boson spectrum
{
	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2),2)),2));
}

long double Quark_Spectrum(long double omega, long double k, long double M, int Temp)	//Single quark spectral function
{
	return(ImSelf_Energy(M, omega, k, Temp)/(pow(pow(omega,2)-pow(k,2)-pow(M,2)-ReSelf_Energy(M, omega, k, Temp),2)+pow(ImSelf_Energy(M, omega, k, Temp),2)));
}

long double Spin_Sum(long double Par[5], long double omega, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
{
	return(4.*(-pow(Par[3],2)/4.+pow(k,2)/*+Par[3]*k*cos(theta)*/+2.*pow(Par[2],2)+omega*(sqrt(Par[4]+pow(Par[3],2))-omega)));
}

long double Folding_Integrand(long double Par[5], long double omega, long double k, long double theta, int Temp)	//Integrand of the folding integral
{
	return(-pow(Par[2],2)*Quark_Spectrum(omega, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*Quark_Spectrum(sqrt(Par[4]+pow(Par[3],2))-omega, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(1.-Fermi(omega, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))-omega, Temp)));
}
