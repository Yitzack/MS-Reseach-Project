//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<complex>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<cfloat>
#include"Elements.h"
using namespace std;

//Integrals that define results and ancillary functions
Elements theta_Int(long double[5], int);	//Integrates the theta results
Elements k_Int(long double[5], int, long double);	//Integrates the k momentum results
Elements Folding(long double[5], int, long double, long double);	//Folding integral, energy integral
long double Newtons_theta(long double, long double, long double, long double);	//Executes a Newton's algorithm search for the maximum of f()
long double D1(long double, long double, long double, long double);	//Finite difference definition, 4th order, 1st derivitive of f
long double D2(long double, long double, long double, long double);	//Finite difference definition, 4th order, 2nd derivitive of f
long double f(long double, long double, long double, long double);	//Analytic integrand of finite P, zero-width result for theta integrand
void Characterize_k_Int(long double[5], int, long double, long double[21], long double[21], int&);	//Returns the poles of the k integral's integrands
long double Upper_Bound(long double, long double, long double, long double, long double);
long double Lower_Bound(long double, long double, long double, long double, long double);
void Characterize_Folding(long double[5], int, long double, long double, long double[6], long double[6], int&);	//Returns the poles of the folding integral's integrands
long double Newton_Method_omega(long double, long double[5], long double, long double, int, long double (*)(long double[5], long double, long double, long double, int));
long double omega_Width(long double, long double[5], long double, long double, int, long double (*)(long double[5], long double, long double, long double, int));

//Straight Functions everything is built from
long double ImSelf_Energy(long double, long double, long double, int);	//Imaginary single quark self energy
long double ReSelf_Energy(long double, long double, long double, int);	//Real single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Potential_on(long double[5]);	//On-shell potential for the on-shell T-Matrix
long double Potential1(long double[5], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[5], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double ImD(long double, long double, long double, int);	//Single quark spectral function
long double ReD(long double, long double, long double, int);	//Single quark spectral function
long double Spin_Sum1(long double[5], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Sum2(long double[5], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double ImFolding_Integrand1(long double[5], long double, long double, long double, int);	//Integrand of the folding integral for positive energy
long double ImFolding_Integrand2(long double[5], long double, long double, long double, int);	//Integrand of the folding integral for negative energy (anti-particle/particle-hole)
long double ReFolding_Integrand1(long double[5], long double, long double, long double, int);	//Integrand of the folding integral for positive energy
long double ReFolding_Integrand2(long double[5], long double, long double, long double, int);	//Integrand of the folding integral for negative energy (anti-particle/particle-hole)

void mergeSort(long double List[], int a, int b)
{
	int i, j, k;
	long double Temp[(a+b)/2-a+1];

	if(b-a > 1)	//Divide
	{
		mergeSort(List, a, (a+b)/2);
		mergeSort(List, (a+b)/2+1, b);
	}

	for(i = 0; i <= (a+b)/2-a; i++)	//Copy out the lower half array in prep for copy over
		Temp[i] = List[i+a];

	j = 0;
	k = (a+b)/2+1;
	for(i = a; i <= b && j <= (a+b)/2-a && k <= b; i++)	//Merge/Conqure while both half lists have not been exhausted
	{
		if(Temp[j] <= List[k])
		{
			List[i] = Temp[j];
			j++;
		}
		else
		{
			List[i] = List[k];
			k++;
		}
	}
	for(; i <= b && j <= (a+b)/2-a; i++)	//If the Temp list has not been exhausted, complete the job
	{
		List[i] = Temp[j];
		j++;
	}

	return;
}

#define GAMMA -.015
long double Boundary[] = {0.00865, 0.0267, 0.0491, 0.0985, .421, .802, 1.01, 4.85, 1.5, 2.5, 3, 4, 5.5, 7.7, 1./17., 0.3, 0.08};

//long double Par[5] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[5], int Temp)
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double x1, x2;	//Abscissa
	if(Par[4] > pow(2.*Par[2],2))
		x1 = Newtons_theta(Par[2], Par[3], Par[4], M_PI/25.);
	else
		x1 = Newtons_theta(Par[2], Par[3], pow(2.*Par[2],2)*1.0001, M_PI/25.);
	if(x1>M_PI/10.)
		x1 = M_PI/10.;
	long double Range[] = {x1*Boundary[14], x1*Boundary[15], x1, x1*(2.-Boundary[15]), x1*(2.-Boundary[15])*(1.-Boundary[16])+M_PI/2.*Boundary[16], M_PI/2., asin(sqrt(-Par[4])/Par[3]),0,0};
	Elements F;	//Sum of ordinate*weights
	Elements Answer = Elements(0,0,0);	//Answer to be returned
	Elements holder;
	long double a = 0, b;	//Sub-interval limits of integration
	int i, j;	//Counters
	//ofstream Table("theta Table", ios::app);

	Range[7] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[7] = acos((pow(Range[7],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(Range[7]*Par[3]));
	Range[8] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[8] = acos((pow(Range[8],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(-Range[8]*Par[3]));

	if(Range[6] != Range[6] || Range[6] < 0)
		Range[6] = M_PI;
	if(Range[7] != Range[7])
		Range[7] = M_PI;
	if(Range[8] != Range[8])
		Range[8] = M_PI;

	mergeSort(Range, 0, 8);

	if(Range[0] == 0)
		i = 1;
	else
		i = 0;

	for(; i < 8 && Range[i] <= M_PI/2. && Par[3] != 0; i++)
	{
		b = Range[i];

		F.null();
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;
			holder = k_Int(Par, Temp, x1);
			//Table << Par[3] << " " << Par[4] << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*sin(x1)*w[j+1];
			holder = k_Int(Par, Temp, x2);
			//Table << Par[3] << " " << Par[4] << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*sin(x2)*w[j+1];
		}
		holder = k_Int(Par, Temp, (a+b)/2.);
		//Table << Par[3] << " " << Par[4] << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		F += holder*sin((a+b)/2.)*w[0];
		Answer += F*(b-a)/2.;
		a = b;
	}

	if(Par[3] == 0)
		return(k_Int(Par, Temp, M_PI/2.)/(2.*pow(M_PI,2)));
	return(Answer/pow(2.*M_PI,2)*2.);
}

long double Newtons_theta(long double M, long double P, long double s, long double Window)
{
	long double theta = 0;
	long double next;
	long double holder;
	int i = 0;

	holder = f(M, P, s, 0);
	next = 0;
	do	//Find the max on mesh
	{
		next += Window;
		if(holder < f(M, P, s, next))
		{
			holder = f(M, P, s, next);
			theta = next;
		}
		else
			break;	//Found an early exit
	}while(next < M_PI/2.);

	holder = theta;
	next = theta-D1(M,P,s,theta)/D2(M,P,s,theta);

	while(abs(theta/next-1.) > .0001 && i < 10)	//Run the actual Newton's algorithm
	{
		theta = next;
		next = theta-D1(M,P,s,theta)/D2(M,P,s,theta);
		i++;
	}

	if(next > holder+Window || next < holder-Window)	//The result is outside the window incated by the search on mesh, try again on tighter mesh
		return(Newtons_theta(M, P, s, Window/10.));

	return(next);
}

long double D1(long double M, long double P, long double s, long double theta)
{
	long double h = .0001;
	return((f(M,P,s,theta-2.*h)/12.-2./3.*f(M,P,s,theta-h)+2./3.*f(M,P,s,theta+h)-f(M,P,s,theta+2.*h)/12.)/h);	//4th order, 1st derivitive
}

long double D2(long double M, long double P, long double s, long double theta)
{
	long double h = .0001;
	return((-f(M,P,s,theta-2.*h)/12.+4./3.*f(M,P,s,theta-h)-2.5*f(M,P,s,theta)+4./3.*f(M,P,s,theta+h)-f(M,P,s,theta+2.*h)/12.)/pow(h,2));	//4th order, 2nd derivitive
}

long double f(long double M, long double P, long double s, long double theta)
{
	return(((-4.*pow(M,2)+s)*pow(pow(P,2)+s,(long double)1.5)*sin(theta))/(8.*abs((sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))*(s+pow(P*sin(theta),2)))/(pow(P,2)+s-P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))))*(s+pow(P*sin(theta),2))*(pow(P,2)+s-P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))));
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements k_Int(long double Par[5], int Temp, long double theta)	//Integrates the k momentum results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Range[] = {-Boundary[13], 0, Boundary[13]};	//Number of gamma from center
	Elements F;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Answer to be returned
	Elements Partial;	//Answer for sub-interval for determining completeness
	Elements holder;
	long double x1, x2;	//Abscissa
	long double a = 0, b = 0;//Sub-interval limits of integration
	int Poles;	//Number of poles
	long double zero[21];	//The real part of the signular pole
	long double gamma[21];	//The distance to the singular, maybe
	int i, j, l;	//Counters
	int Intervals;
	//ofstream Table("k Table", ios::app);
	//ofstream Poles_Table("k Poles", ios::app);

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);
	//for(i = 0; i < Poles; i++)
		//Poles_Table << Par[3] << " " << Par[4] << " "  << theta << " " << zero[i] << " " << gamma[i] << endl;
	long double Stops[Poles*3+8];

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(gamma[i] == gamma[i])	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 3; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Stops as required by poles
				l++;
			}
		else	//At lease insert the central point of the pole so that a certain amount of information isn't lost
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//This the upper bound that the vacuum calls for, Partial/total will promote higher as needed
	if(Stops[l] != Stops[l])
		Stops[l] = .5*sqrt(-Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
	Stops[l+1] = abs(.5*(Par[3]*cos(theta)+sqrt(4.*(pow(Par[3],2)+Par[4]-pow(Par[2],2))-pow(Par[3]*sin(theta),2))));	//On-shells leaving the range -E/2 to E/2
	Stops[l+2] = abs(.5*(Par[3]*cos(theta)-sqrt(4.*(pow(Par[3],2)+Par[4]-pow(Par[2],2))-pow(Par[3]*sin(theta),2))));
	Stops[l+3] = abs(.5*(Par[3]*cos(theta)+sqrt(-4.*pow(Par[2],2)-pow(Par[3]*sin(theta),2))));
	Stops[l+4] = abs(.5*(Par[3]*cos(theta)-sqrt(-4.*pow(Par[2],2)-pow(Par[3]*sin(theta),2))));
	Stops[l+5] = abs((8.*pow(Par[2],2)*Par[3]*cos(theta)+sqrt((pow(Par[3],2)+Par[4])*(pow(2.*Par[2],4)+7.*pow(Par[3],4)+8.*pow(Par[2],2)*(pow(Par[3],2)-Par[4])-8.*pow(Par[3],2)*Par[4]+pow(Par[4],4)+2.*pow(Par[2],2)*(8.*pow(Par[2],2)+7.*pow(Par[3],2)-Par[4])*cos(2.*theta))))/(4.*(pow(Par[3],2)-Par[4]+2.*pow(Par[3],2)*cos(theta))));	//On-shell leaving the range -E/2+k+ to E/2-k-
	Stops[l+6] = abs(2.*pow(Par[3],3)*cos(theta)+2.*Par[3]*Par[4]*cos(theta)+sqrt((pow(Par[3],2)+Par[4])*(16.*pow(Par[2],4)-5.*pow(Par[3],4)-4.*pow(Par[3],2)*Par[4]+pow(Par[4],2)-8.*pow(Par[2],2)*(pow(Par[3],2)+Par[4])+2.*pow(Par[3],2)*(pow(Par[3],2)+Par[4])*cos(2.*theta))))/(4.*(pow(Par[3],2)+Par[4]));
	Stops[l+7] = sqrt(Par[4]+pow(Par[3],2))/2.;	//Potiential peak leaves -E/2 to E/2

	for(i = l; i < l+8; i++)
		if(Stops[i] != Stops[i])
			Stops[i] = -1;

	mergeSort(Stops, 0, l+7);

	i = 0;
	j = 0;
	while(Stops[j] <= 0)
		j++;
	for(; j < l+8; j++)
	{
		if(((i > 0 && Stops[i-1] != Stops[j]) || i == 0) && Stops[j] == Stops[j])
		{
			Stops[i] = Stops[j];
			i++;
		}
		else if(Stops[j] != Stops[j])
			break;
	}
	Intervals = i;

	if(j == 0)
		Intervals = 1;

	i = 0;
	do
	{
		if((i < Intervals && b+100 < Stops[i]) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && b+50 < Stops[i]) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && b+10 < Stops[i]) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && b+3 < Stops[i]) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}
		else
			b += 3;

		F.null();
		for(l = 0; l < 9; l++)
		{
			x1 = (b+a-Disp[l]*(b-a))/2.; //Actual evaluation points
			x2 = (b+a+Disp[l]*(b-a))/2.;

			holder = Folding(Par, Temp, x1, theta);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*pow(x1,2)*w[l+1]; //Evaluate function at x1
			holder = Folding(Par, Temp, x2, theta);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*pow(x2,2)*w[l+1]; //Evaluate function at x2
		}
		holder = Folding(Par, Temp, (a+b)/2., theta);
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		F += holder*pow((a+b)/2.,2)*w[0]; //Evaluate function at (a+b)/2.
		Partial = F*(b-a)/(2.);
		Answer += Partial;
		a = b;
	}while(!(Partial == 0) && (i < Intervals || abs(Partial/Answer) >= .0001) && a <= 20.*sqrt(Par[4]+pow(Par[3],2))); //k bigger than 20E is getting pretty stupid, should be sneaking up on 10^-5 of the answer left

	return(Answer);
}

void Characterize_k_Int(long double Par[5], int Temp, long double theta, long double zero[21], long double gamma[21], int &Poles) //Returns the poles of the k integral's integrands
{
	long double holder;
	int i, j, l;

	Poles = 4;
	zero[0] = .5*Par[3]*abs(cos(theta));	//Near intersection of 2 on-shells
	gamma[0] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[0], theta), zero[0], Temp);
	zero[1] = (pow(Par[2],2)-Par[4]/4.)*(sqrt(Par[4]+pow(Par[3],2))+Par[3]*cos(theta))/(Par[4]-pow(Par[3]*sin(theta),2));	//On-shell momentum intersection with V peak
	zero[2] = (Par[4]/4.-pow(Par[2],2))*(sqrt(Par[4]+pow(Par[3],2))-Par[3]*cos(theta))/(Par[4]-pow(Par[3]*sin(theta),2));
	zero[3] = (pow(Par[2],2)-Par[4]/4.)*(sqrt(Par[4]+pow(Par[3],2))-Par[3]*cos(theta))/(Par[4]-pow(Par[3]*sin(theta),2));
	zero[4] = (Par[4]/4.-pow(Par[2],2))*(sqrt(Par[4]+pow(Par[3],2))+Par[3]*cos(theta))/(Par[4]-pow(Par[3]*sin(theta),2));
	gamma[1] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[1], theta), zero[1], Temp);	//on-shell width at that point, gaussian doesn't add any extra width here
	gamma[2] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[2], theta), zero[2], Temp);
	gamma[3] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[3], theta), zero[3], Temp);
	gamma[4] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[4], theta), zero[4], Temp);

	if((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3]*sin(theta),2)) > 0.)
	{
		zero[Poles] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//on-shell momentum, it's almost like every on-shell peak crosses every other one at this point
		gamma[Poles] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp);
		Poles++;
	}

	for(i = 0; i < Poles; i++)
		zero[i] = abs(zero[i]);

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
		}
	}

	j = 0;
	for(i = 0; i < Poles; i++)
	{
		if(zero[i] > 1000)
			break;

		if(zero[i] == zero[j] && i != j){}
		else if(zero[i] == zero[j] && i == j)
			j++;
		else
		{
			zero[j] = zero[i];
			j++;
		}
	}
	Poles = j;

	return;
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements Folding(long double Par[5], int Temp, long double k, long double theta)	//Folding integral, energy integral
{
	if(Temp == 0 && abs(sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta)-Energy(0,Par[3]/2.,k,theta)) < 1e-12)	//Let's save some time and just return 0, because it is
		return(Elements(0,0,0));
	else if(Par[4]+pow(Par[3],2) < 0)
		return(Elements(0,0,0));	//Bad data trap and time saver

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204}; //Weight of the function at Disp
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};	//Number of gamma from center
	long double a, b;	//Sub-interval limits of integration
	long double Max;	//Upper limit of integration
	Elements F;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Results to be returned
	Elements Partial(0,0,0);//Partial Answer
	Elements holder;
	//long double x1, x2;	//Abscissa
	long double zero[13];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[13];	//Imaginary part of poles
	int Intervals;		//Number of intervals required by poles and discontinuities
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	//ofstream Table("omega Table", ios::app);
	//ofstream Poles_Table("omega Poles", ios::app);

	Characterize_Folding(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	//for(i = 0; i < Poles; i++)
		//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << zero[i] << " " << gamma[i] << endl;
	long double Stops[Poles*17+4];	//Intervals that are required by integrating near poles

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(gamma[i] == gamma[i])	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 17; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Stops as required by poles
				l++;
			}
		else	//At lease insert the central point of the pole so that a certain amount of information isn't lost
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
	Stops[l+1] = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);
	Stops[l+2] = sqrt(Par[4]+pow(Par[3],2))/2.;
	Stops[l+3] = -sqrt(Par[4]+pow(Par[3],2))/2.;
	/*Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << Stops[l] << " " << 1 << endl;
	//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << Stops[l+1] << " " << 1 << endl;
	//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << Stops[l+2] << " " << 1 << endl;
	//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << Stops[l+3] << " " << 1 << endl;*/

	mergeSort(Stops, 0, l+3);

	if(Temp != 0)
	{
		a = b = -sqrt(Par[4]+pow(Par[3],2))/2.;
		Max = sqrt(Par[4]+pow(Par[3],2))/2.;
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
		Max = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);
	}

	i = 0;
	j = 0;
	while(Stops[j] < a)
		j++;
	for(; j < l+4; j++)
	{
		if(((i > 0 && Stops[i-1] != Stops[j]) || i == 0) && Stops[j] <= Max)
		{
			Stops[i] = Stops[j];
			i++;
		}
		else if(Stops[j] > Max)
			break;
	}
	Intervals = i;

	i = 1;
	do
	{
		if((i < Intervals && b+100 < Stops[i]) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && b+50 < Stops[i]) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && b+10 < Stops[i]) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && b+3 < Stops[i]) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}

		F.null();
		#pragma omp parallel for
		for(l = 0; l < 9; l++)	//Integrate the sub-interval
		{
			long double x1 = (b+a-Disp[l]*(b-a))/2.;
			long double x2 = (b+a+Disp[l]*(b-a))/2.;

			holder = (Elements(Spin_Sum1(Par, x1, k, theta), 2.*Potential1(Par,x1,k), Potential2(Par,x1,k))*ImFolding_Integrand1(Par,x1,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -x1, k, theta), 2.*Potential1(Par,-x1,k), Potential2(Par,-x1,k))*ImFolding_Integrand2(Par,-x1,k,theta,Temp)*/);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
			holder = (Elements(Spin_Sum1(Par, x2, k, theta), 2.*Potential1(Par,x2,k), Potential2(Par,x2,k))*ImFolding_Integrand1(Par,x2,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -x2, k, theta), 2.*Potential1(Par,-x2,k), Potential2(Par,-x2,k))*ImFolding_Integrand2(Par,-x2,k,theta,Temp)*/);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
		}
		holder = (Elements(Spin_Sum1(Par, (a+b)/2., k, theta), 2.*Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*ImFolding_Integrand1(Par,(a+b)/2.,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -(a+b)/2., k, theta), 2.*Potential1(Par,-(a+b)/2.,k), Potential2(Par,-(a+b)/2.,k))*ImFolding_Integrand2(Par,-(a+b)/2.,k,theta,Temp)*/);
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		F += holder*w[0];

		Partial = F*(b-a)/(2.);
		Answer += Partial;
		a = b;
	}while((!(Partial == 0) || a < Max) && (i < Intervals || ((abs(Partial/Answer) >= .0001)/* && Temp != 0*/)));

	/*if(Temp != 0)
	{
		a = b = Stops[0];

		do
		{
			if(b+100 < Stops[0])	//Middle of nowhere intervals below Stops[0]
				a -= 100;
			else if(b+50 < Stops[0])
				a -= 50;
			else if(b+10 < Stops[0])
				a -= 10;
			else
				a -= 3;

			F.null();
			#pragma omp parallel for
			for(l = 0; l < 9; l++)	//Integrate the sub-interval
			{
				long double x1 = (b+a-Disp[l]*(b-a))/2.;
				long double x2 = (b+a+Disp[l]*(b-a))/2.;

				holder = (Elements(Spin_Sum1(Par, x1, k, theta), 2.*Potential1(Par,x1,k), Potential2(Par,x1,k))*ImFolding_Integrand1(Par,x1,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -x1, k, theta), 2.*Potential1(Par,-x1,k), Potential2(Par,-x1,k))*ImFolding_Integrand2(Par,-x1,k,theta,Temp));
				//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x1 << " " << holder.store(1) << " " << holder.store(2) << " " << holder.store(2) << endl;
				F += holder*w[l+1];
				holder = (Elements(Spin_Sum1(Par, x2, k, theta), 2.*Potential1(Par,x2,k), Potential2(Par,x2,k))*ImFolding_Integrand1(Par,x2,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -x2, k, theta), 2.*Potential1(Par,-x2,k), Potential2(Par,-x2,k))*ImFolding_Integrand2(Par,-x2,k,theta,Temp));
				//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
				F += holder*w[l+1];
			}
			holder = (Elements(Spin_Sum1(Par, (a+b)/2., k, theta), 2.*Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*ImFolding_Integrand1(Par,(a+b)/2.,k,theta,Temp)/*+Elements(Spin_Sum2(Par, -(a+b)/2., k, theta), 2.*Potential1(Par,-(a+b)/2.,k), Potential2(Par,-(a+b)/2.,k))*ImFolding_Integrand2(Par,-(a+b)/2.,k,theta,Temp));
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[0];

			Partial = F*(b-a)/(2.);
			Answer += Partial;
			b = a;
		}while(abs(Partial/Answer) >= .0001);
	}*/

	return(Answer/M_PI);
}

void Characterize_Folding(long double Par[5], int Temp, long double k, long double theta, long double zero[10], long double gamma[10], int &Poles)
{
	long double Lower, Upper;	//Limits of integration in Folding, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(true)//Temp != 0)
	{
		Lower = -LDBL_MAX;
		Upper = LDBL_MAX;	//Integrate from 0 to E and twice E to infinity (ie, I need all points greater than 0)
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta);
		Upper = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	/*zero[0] = pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));	//Potential poles (Lambda^4), I know exactly where these are at.
	zero[1] = -pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[0] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));
	gamma[1] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));
	zero[2] = -sqrt(Par[4]+pow(Par[3],2))+pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	zero[3] = -sqrt(Par[4]+pow(Par[3],2))-pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[2] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));
	gamma[3] = abs(.5*(pow(pow(2.*k,4)+pow(Par[1],4),(long double).25)*sin(.5*atan(pow(Par[1]/(2.*k),2)))));*/
	zero[0] = k;	//Potential "poles" (gaussian), I know exactly where these are at.
	zero[1] = -k;
	gamma[0] = Par[1];
	gamma[1] = Par[1];
	zero[2] = -sqrt(Par[4]+pow(Par[3],2))+k;
	zero[3] = -sqrt(Par[4]+pow(Par[3],2))-k;
	gamma[2] = Par[1];
	gamma[3] = Par[1];
	zero[4] = 0;
	gamma[4] = .009;

	holder = GAMMA;
	zero[5] = real(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.))))-sqrt(Par[4]+pow(Par[3],2))/2.;	//Exact vacuum
	zero[6] = real(-sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.))))-sqrt(Par[4]+pow(Par[3],2))/2.;
	zero[7] = sqrt(Par[4]+pow(Par[3],2))/2.+real(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.))));
	zero[8] = sqrt(Par[4]+pow(Par[3],2))/2.-real(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.))));

	if(Temp != 0)	//media estimate
	{
		zero[9] = -Newton_Method_omega(-zero[5], Par, k, theta, Temp, ImFolding_Integrand2);	//Negative to counter the fact that I'm doing particle-hole in reverse order to minimize region of interest
		zero[8] = -Newton_Method_omega(-zero[6], Par, k, theta, Temp, ImFolding_Integrand2);
		zero[11] = -Newton_Method_omega(-zero[7], Par, k, theta, Temp, ImFolding_Integrand2);
		zero[12] = -Newton_Method_omega(-zero[8], Par, k, theta, Temp, ImFolding_Integrand2);
		zero[5] = Newton_Method_omega(zero[5], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[6] = Newton_Method_omega(zero[6], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[7] = Newton_Method_omega(zero[7], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[8] = Newton_Method_omega(zero[8], Par, k, theta, Temp, ImFolding_Integrand1);

		gamma[5] = omega_Width(zero[5], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[6] = omega_Width(zero[6], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[7] = omega_Width(zero[7], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[8] = omega_Width(zero[8], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[9] = omega_Width(-zero[9], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[10] = omega_Width(-zero[10], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[11] = omega_Width(-zero[11], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[12] = omega_Width(-zero[12], Par, k, theta, Temp, ImFolding_Integrand2);
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[6] = gamma[5] = abs(imag(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.)))));
		gamma[8] = gamma[7] = abs(imag(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.)))));
	}

	if(Temp == 0)
		i = 8;
	else
		i = 12;

	for(; i >= 0; i--)	//Bubble sort
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

	i = j = 0;	//Find the first zero greater than lower bound
	while(zero[i] < Lower) i++;

	while(zero[i] <= Upper && (i < 13 || (i < 9 && Temp == 0)))	//Move zeroes up to front of array, count off poles within the limits of integration
	{
		zero[j] = zero[i];
		gamma[j] = abs(gamma[i]);
		i++;
		j++;
	}
	Poles = j;

	return;
}

long double Newton_Method_omega(long double omega, long double Par[5], long double k, long double theta, int Temp, long double (*Folding)(long double[5], long double, long double, long double, int))
{
	long double new_omega;
	const long double h = 1e-4;
	int i = 0;

	new_omega = omega - .5*h*(1./Folding(Par, omega+h, k, theta, Temp)-1./Folding(Par, omega-h, k, theta, Temp))/((1./Folding(Par, omega-h, k, theta, Temp)-2./Folding(Par, omega, k, theta, Temp)+1./Folding(Par, omega+h, k, theta, Temp)));

	while(abs(1.-new_omega/omega) > 1e-5 && i <= 10)
	{
		omega = new_omega;
		new_omega = omega - .5*h*(1./Folding(Par, omega+h, k, theta, Temp)-1./Folding(Par, omega-h, k, theta, Temp))/((1./Folding(Par, omega-h, k, theta, Temp)-2./Folding(Par, omega, k, theta, Temp)+1./Folding(Par, omega+h, k, theta, Temp)));
		i++;
	}

	return(omega);
}

long double omega_Width(long double zero, long double Par[5], long double k, long double theta, int Temp, long double (*Folding)(long double[5], long double, long double, long double, int))
{
	return(sqrt(abs(2e-10*Folding(Par, zero, k, theta, Temp)/(Folding(Par, zero-1e-5, k, theta, Temp)-2.*Folding(Par, zero, k, theta, Temp)+Folding(Par, zero+1e-5, k, theta, Temp)))));
}

//long double Par[5] = {g, Lambda, M, P, s}
long double ImSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Par[6];	//Momentum dependance parameterization
	long double E_0 = Energy(M,k,0,0);	//location of lorentzian
	long double Sigma;	//size of energy dependance
	long double b1, b2;	//slope of exponential decrease to left and right
	long double Delta;	//concavity or length of transition from left to right

	switch(Temp)
	{
		case 0:
			if(pow(omega,2)>=pow(k,2))
				return(sqrt(pow(omega,2)-pow(k,2))*GAMMA);
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

	if(pow(omega,2)>=pow(k,2))
		return(2.*M*(Par[0]*exp(-pow(k/Par[1],2))+(1-Par[0])*exp(-pow(k/Par[2],2))+Par[3])*(Sigma*exp(Delta+(b1-b2)*(omega-E_0)*E_0/2.-sqrt(b1*b2*pow((omega-E_0)*E_0,2)+pow(Delta+(b1-b2)*(omega-E_0)*E_0/2.,2))))+sqrt(pow(omega,2)-pow(k,2))*GAMMA);
	else
		return(2.*M*(Par[0]*exp(-pow(k/Par[1],2))+(1-Par[0])*exp(-pow(k/Par[2],2))+Par[3])*(Sigma*exp(Delta+(b1-b2)*(omega-E_0)*E_0/2.-sqrt(b1*b2*pow((omega-E_0)*E_0,2)+pow(Delta+(b1-b2)*(omega-E_0)*E_0/2.,2)))));
}

long double ReSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Par[6];	//Momentum dependance parameterization
	long double E_0 = Energy(M,k,0,0);	//location of lorentzian
	long double a;	//size of energy dependance
	long double gamma;	//width of lorentzian
	long double c;	//zero crossing, might be better as the on-shell energy

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
			if(omega >= 0)	//Fermi factor for vacuum
				return(0.);
			else
				return(1.);
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
//	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(Par[4]-4.*pow(Par[2],2),2)),2));
	return(Par[0]/(2.*M_PI*pow(Par[1],2))*exp(-pow((Par[4]-4.*pow(Par[2],2))/Par[1],2)));
}

long double Potential1(long double Par[5], long double k0, long double k)	//Potiential for the numerator of the boson spectrum
{
//	return(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(4.*pow(k0,2)-4.*pow(k,2),2)));
	return(1./sqrt(2.*M_PI)/Par[1]*exp(-pow((4.*pow(k0,2)-4.*pow(k,2))/Par[1],2)/2.));
}

long double Potential2(long double Par[5], long double k0, long double k)	//Potiential for the denominator of the T-Matrix and boson spectrum
{
//	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(4.*pow(k0,2)-4.*pow(k,2),2)),2));
	return(Par[0]/(2.*M_PI*pow(Par[1],2))*exp(-pow((4.*pow(k0,2)-4.*pow(k,2))/Par[1],2)));
}

long double ImD(long double omega, long double k, long double M, int Temp)	//Single quark spectral function
{
	return(ImSelf_Energy(M, omega, k, Temp)/(pow(pow(omega,2)-pow(k,2)-pow(M,2)-ReSelf_Energy(M, omega, k, Temp),2)+pow(ImSelf_Energy(M, omega, k, Temp),2)));
}

long double ReD(long double omega, long double k, long double M, int Temp)	//Single quark spectral function
{
	return((pow(omega,2)-pow(k,2)-pow(M,2)-ReSelf_Energy(M, omega, k, Temp))/(pow(pow(omega,2)-pow(k,2)-pow(M,2)-ReSelf_Energy(M, omega, k, Temp),2)+pow(ImSelf_Energy(M, omega, k, Temp),2)));
}

long double Spin_Sum1(long double Par[5], long double k0, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), strictly pseudoscalar for now
{
	return(2.*(Par[4]/4.+pow(k,2)-pow(k0,2)+pow(Par[2],2)));
}

long double Spin_Sum2(long double Par[5], long double k0, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), strictly pseudoscalar for now
{
	return(2.*(Par[4]/4.+pow(k,2)-pow(k0,2)-pow(Par[2],2)));
}

long double ReFolding_Integrand1(long double Par[5], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	if(Temp == 0)
		return(-pow(Par[2],2)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp));
	return(-pow(Par[2],2)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(1.-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp))-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(1.-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)));
	/*long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204}; //Weight of the function at Disp
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};	//Number of gamma from center
	long double a, b;	//Sub-interval limits of integration
	long double LocalPar[5] = {Par[0],Par[1],Par[2],Par[3],Par[4]};
	long double Min, Max;	//need to store these so that I know over what range I completed the integral
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;		//Results to be returned
	long double Partial = 0;	//Partial Answer
	long double holder;
	long double Width;	//The length of this interval; don't store distance to travel in it
	//long double x1, x2;	//Abscissa
	long double zero[3];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[3];	//Imaginary part of poles
	int Intervals;		//Number of intervals required by poles and discontinuities
	int i, j, l;		//Counting varibles
	//ofstream Table("ReG Table", ios::app);
	//ofstream Poles_Table("ReG Poles", ios::app);

	zero[0] = 4.*pow(k0,2)+4.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(Par[3],2)-8.*abs(k0)*Energy(Par[2],Par[3]/2.,-k,theta);
	zero[1] = 4.*pow(k0,2)+4.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(Par[3],2)+8.*abs(k0)*Energy(Par[2],Par[3]/2.,-k,theta);
	gamma[0] = gamma[1] = abs(imag(sqrt(complex<long double>(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(GAMMA,2)/2.,sqrt(pow(Par[2]*GAMMA,2)-pow(GAMMA,4)/4.)))));
	for(i = 0; i < 2; i++)
		//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << k0 << " " << zero[i] << " " << gamma[i] << endl;
	long double Stops[56];	//Intervals that are required by integrating near poles

	l = 0;
	for(i = 0; i < 2; i++)
	{
		if(gamma[i] == gamma[i])	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 17; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Stops as required by poles
				l++;
			}
		else	//At lease insert the central point of the pole so that a certain amount of information isn't lost
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = 4.*pow(k0+Energy(0,Par[3]/2.,-k,theta),2)-pow(Par[3],2);
	Stops[l+1] = 4.*pow(k0-Energy(0,Par[3]/2.,-k,theta),2)-pow(Par[3],2);
	Stops[l+2] = -pow(Par[3],2);
	Stops[l+3] = 4.*pow(k0,2)-pow(Par[3],2);
	Stops[l+4] = Par[4];	//Inserted here as a cleaver trick to ensure that Par[4] is not an abscissa and cause 0./0.=nan

	mergeSort(Stops, 0, l+4);

	i = 0;
	j = 0;
	while(Stops[j] < -pow(Par[3],2))
		j++;
	for(; j < l+5; j++)
	{
		if((i > 0 && abs(1.-Stops[i-1]/Stops[j])>1e-8) || i == 0)
		{
			Stops[i] = Stops[j];
			i++;
		}
	}
	Intervals = i;

	a = b = Stops[0];

	i = 1;
	Min = a;
	do
	{
		if((i < Intervals && b+100 < Stops[i]) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && b+50 < Stops[i]) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && b+10 < Stops[i]) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && b+3 < Stops[i]) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}
		else
			b += 3;

		Width = b-a;

		F = 0;
		#pragma omp parallel for
		for(l = 0; l < 9; l++)	//Integrate the sub-interval
		{
			long double x1 = (b+a-Disp[l]*(b-a))/2.;
			long double x2 = (b+a+Disp[l]*(b-a))/2.;

			LocalPar[4] = x1;
			holder = (ImFolding_Integrand1(LocalPar, k0, k, theta, Temp)-ImFolding_Integrand1(Par, k0, k, theta, Temp))/(x1-Par[4]);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << k0 << " " << x1 << " " << holder << endl;
			F += holder*w[l+1];
			LocalPar[4] = x2;
			holder = (ImFolding_Integrand1(LocalPar, k0, k, theta, Temp)-ImFolding_Integrand1(Par, k0, k, theta, Temp))/(x2-Par[4]);
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << k0 << " " << x2 << " " << holder << endl;
			F += holder*w[l+1];
		}
		LocalPar[4] = (a+b)/2.;
		holder = (ImFolding_Integrand1(LocalPar, k0, k, theta, Temp)-ImFolding_Integrand1(Par, k0, k, theta, Temp))/((a+b)/2.-Par[4]);
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << k0 << " " << (a+b)/2. << " " << holder << endl;
		F += holder*w[0];

		Partial = F*Width/(2.);
		Answer += Partial;
		a = b;
	}while(i < Intervals || abs(Partial/Answer)*3./Width >= .0001);
	Max = a;

	if(Par[4] != Min)
		Answer += ImFolding_Integrand1(Par, k0, k, theta, Temp)*log((Max-Par[4])/(Par[4]-Min));
	else
	{
		LocalPar[4] = .0001;
		Answer += 2.*ImFolding_Integrand1(LocalPar, k0, k, theta, Temp)*log((Max-LocalPar[4])/(LocalPar[4]-Min));
		LocalPar[4] = .0002;
		Answer -= ImFolding_Integrand1(LocalPar, k0, k, theta, Temp)*log((Max-LocalPar[4])/(LocalPar[4]-Min));
	}

	return(Answer);*/
}

long double ReFolding_Integrand2(long double Par[5], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for negitive energy (anti-particle/particle-hole)
{
	return(-pow(Par[2],2)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp))-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ReD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)));
}

long double ImFolding_Integrand1(long double Par[5], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	return(-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(1.-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)));
}

long double ImFolding_Integrand2(long double Par[5], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for negitive energy (anti-particle/particle-hole)
{
	return(-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*(Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)));
}
