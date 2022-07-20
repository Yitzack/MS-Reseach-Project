//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<cfloat>
#include<complex>
#include"Elements.h"
using namespace std;

//Integrals that define results
Elements theta_Int(long double[], int);						//Theta integral
Elements k_Int(long double[], int, long double);					//k integral
Elements k0_Int(long double[], int, long double, long double);			//k0 integral aka energy integral
long double Dispersion(long double[], int, long double, long double, long double);	//Dispersion relation for turning ImG_12 into ReG_12

//Functions for finding points of interest in the k integral
void Characterize_k_Int(long double[], int, long double, long double[], long double[], int&);	//Returns the poles of the k integral's integrands
bool Newton_Method_k(long double&, long double, long double, long double, long double, long double, long double(*)(long double, long double), long double(*)(long double, long double, long double, long double, long double));	//Returns the k-intesection of two features of interest: potential peak, on-shell peak, light-like boundaries
long double V_Plus(long double, long double);			//Potiential peaks, matches a potential not in use
long double V_Minus(long double, long double);
long double Emm(long double, long double, long double, long double, long double);	//on-shell peaks
long double Epm(long double, long double, long double, long double, long double);
long double mEmp(long double, long double, long double, long double, long double);
long double Emp(long double, long double, long double, long double, long double);
long double mEpp(long double, long double, long double, long double, long double);
long double Epp(long double, long double, long double, long double, long double);
long double Upper_Bound(long double, long double, long double, long double, long double);	//light-like boundaries, where the imaginary self-energy has a known cusp
long double Lower_Bound(long double, long double, long double, long double, long double);

//Functions for finding points of interest in the k0 integral
void Characterize_k0_Int(long double[], int, long double, long double, long double[], long double[], int&);	//Returns the poles of the k0 integral's integrands
long double Newton_Method_k0(long double, long double[], long double, long double, int, long double (*)(long double[], long double, long double, long double, int));	//Returns the k0 of the on-shell peak using Newton's method on 1/f
long double omega_Width(long double, long double[], long double, long double, int, long double (*)(long double[], long double, long double, long double, int));	//Returns the width of on-shell peak using the assumption of a breit-wigner peak 

//Functions for finding points of interest in the dispersion integral
void Characterize_Dispersion(long double[], int, long double, long double, long double, long double[], long double[], int&);
long double sp_Width(long double[], long double, long double, long double, int, long double (*)(long double[], long double, long double, long double, int));	//Breit-Wigner width of the peak

//Functions that return physics for the integrand
void ImSelf_Energy(long double, long double, long double[], int, long double[]);			//Returns the imaginary single quark self-energies for both quarks, contains an alternate T=194 MeV solution
long double ImSelf_Energy(long double, long double, long double, int);				//Returns the imaginary single quark self-energies for one quark, contains an alternate T=194 MeV solution
void ReSelf_Energy(long double, long double, long double[], int, long double[]);					//Returns the real single quark self-energies for both quarks, contains an alternate T=194 MeV solution
void Self_Energy(long double, long double, long double[], int, long double[], long double[]);	//Returns the complex single quark self-energies for both quarks, is a simple Breit-Wigner self-energy and alternate to those above
long double Energy(long double, long double, long double, long double);						//Single quark energy, also used to return total momentum by setting M=0
long double Fermi(long double, int);											//Fermi function
long double Set_Temp(int);												//Decodes 0-4 into numeric temprature for Fermi factor
long double Potential1(long double[], long double, long double);							//One vertex of the potiential
long double Potential2(long double[], long double, long double);							//Two vertices of the potiential
long double Interacting_Linear_Trace(long double[]);									//Linear (linear in sqrt(s)) contribution to the interacting trace
long double Interacting_Quad_Trace(long double[], long double, long double);						//Quadratic contribution to the interacting trace
long double Imk0_Integrand(long double[], long double, long double, long double, int);				//Integrand of the k0 integral for positive energy

auto Start_Time = chrono::system_clock::now();
int Allotment = 90;

//Merge Sort. There are a number of semi-sorted lists that need sorting. This will probably beat quick sort under similar conditions. Don't worry about the second element. It doesn't need sorting as either they don't matter or they will be assigned after sorting.
void mergeSort(long double List[], int a, int b)
{
	int i, j, k;
	long double Temp[(a+b)/2-a+1];

	if(b-a > 1)	//Divide...
	{
		mergeSort(List, a, (a+b)/2);
		mergeSort(List, (a+b)/2+1, b);
	}

	for(i = 0; i <= (a+b)/2-a; i++)	//Copy out the lower half array in prep for copy over
	{
		Temp[i] = List[i+a];
	}

	j = 0;
	k = (a+b)/2+1;
	for(i = a; i <= b && j <= (a+b)/2-a && k <= b; i++)	//... and conqure while both half lists have not been exhausted
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

#ifndef GAMMA	//use option -D GAMMA=<number> to alter single particle vacuum width, default value is 15MeV
#define GAMMA -.015	//Width of single quark propagator
#endif

Elements theta_Int(long double Par[], int Temp)
{
	if(Par[3] == 0)	//Short cut for P=0, theta integral is analytic
		return(k_Int(Par, Temp, M_PI/2.)/pow(2.*M_PI, 2)*2.);

//37th order Gauss-Legendre integration
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight
	long double x1, x2;					//Abscissa
	long double a = 0, b;					//Sub-interval limits of integration
	long double Boundary_theta[] = {1./17., 0.3, 0.08};	//Extra boundary values
	Elements F;						//Sum of ordinate*weights
	Elements Answer(0, 0, 0, 0);				//Answer to be returned
	int i, j;						//Counters

	if(Par[4] > 0 && Par[3] > sqrt(Par[4]/2.)) //Where the maximum of the theta integral ought to land. It might only be correct for BbS reduction, but is close enough for all other cases. Only valid for s>0 and P>sqrt(s/2)
		x1 = asin(sqrt(Par[4]/2.)/Par[3]);
	else
		x1 = M_PI/10.;	//If it isn't valid, value is needed anyways to split up the integral

	//Don't get too close to the pole or details might get lost
	if(x1>M_PI/10.)
		x1 = M_PI/10.;

	//List of boundaries between subintervals
	long double Range[] = {x1*Boundary_theta[0], x1*Boundary_theta[1], x1, x1*(2.-Boundary_theta[1]), x1*(2.-Boundary_theta[1])*(1.-Boundary_theta[2])+M_PI/2.*Boundary_theta[2], M_PI/2., asin(sqrt(-Par[4])/Par[3]), 0, 0};

	//Some kind of intersection, probably between the simultanous on-shell and potential peak, don't rightly remember
	Range[7] = sqrt(4.*pow(Par[3], 4)+8.*pow(Par[3], 2)*Par[4]+4.*pow(Par[4], 2)-pow(Par[1], 4))/pow(256.*pow(Par[3], 4)+512.*pow(Par[3], 2)*Par[4]+256.*pow(Par[4], 2), (long double).25);
	Range[7] = acos((pow(Range[7], 2)+pow(Par[2], 2)-Par[4]-(long double).75*pow(Par[3], 2))/(Range[7]*Par[3]));
	Range[8] = sqrt(4.*pow(Par[3], 4)+8.*pow(Par[3], 2)*Par[4]+4.*pow(Par[4], 2)-pow(Par[1], 4))/pow(256.*pow(Par[3], 4)+512.*pow(Par[3], 2)*Par[4]+256.*pow(Par[4], 2), (long double).25);
	Range[8] = acos((pow(Range[8], 2)+pow(Par[2], 2)-Par[4]-(long double).75*pow(Par[3], 2))/(-Range[8]*Par[3]));

	//Bad data trap for NaN and negative boundaries. These are only ones that can NaN or return negative numbers
	if(isnan(Range[6]) || Range[6] < 0) Range[6] = M_PI;
	if(isnan(Range[7])) Range[7] = M_PI;
	if(isnan(Range[8])) Range[8] = M_PI;

	//Put in asending order
	mergeSort(Range, 0, 8);

	for(i = 0; i < 8 && Range[i] <= M_PI/2.; i++)	//Count through pre-determined intervals
	{
		b = Range[i];	//Upper edge
		F.null();	//Zero out F for a new round

		//Count through points away from center
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F += k_Int(Par, Temp, x1)*sin(x1)*w[j+1];
			F += k_Int(Par, Temp, x2)*sin(x2)*w[j+1];
		}
		F += k_Int(Par, Temp, (a+b)/2.)*sin((a+b)/2.)*w[0];

		Answer += F*(b-a)/2.;	//Add the subinterval to total of the integral

		a = b;	//Upper edge becomes lower edge
	}

	return(Answer/pow(2.*M_PI, 2)*2.);
}

Elements k_Int(long double Par[], int Temp, long double theta)
{
//37th order Gauss-Legendre integration
	long double Disp37[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center
	long double w37[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight
//97th order Gauss-Legendre integration
	long double Disp97[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center
	long double w97[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight
	long double x1, x2;	//Abscissa
	long double a, b;	//Sub-interval limits of integration

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.

	Elements F;			//Sum of ordinates*weights
	Elements Answer(0, 0, 0, 0);	//Answer to be returned
	Elements Partial;		//Answer for sub-interval for determining completeness

	int Poles;		//Number of poles
	long double zero[26];	//The real part of the signular pole
	long double gamma[26];	//The distance to the singular, maybe
	long double Max;
	int i, j, l;		//Counters, would use 'k', but 'k' is occupied by relative 3-momenta in other parts of program
	int Intervals;		//Number of intervals recorded in Stops

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);	//Find the location of the complex poles
	long double Stops[Poles+12];				//List of pre-determined subintervals

	for(l = 0; l < Poles; l++)	//Counting through the poles
	{
		Stops[l] = zero[l];
	}

	//More intervals from features not already considered
	Stops[l] = .5*sqrt(Par[4]*(Par[4]+pow(Par[3], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)));	//k for which quarks are simultanous light-like, highest k needed for vacuum
	if(isnan(Stops[l]))	//If meson is space-like, keep absolute value of it anyways even though it probably does nothing
		Stops[l] = .5*sqrt(-Par[4]*(Par[4]+pow(Par[3], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)));
	Stops[l+1] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]-pow(2.*Par[2], 2)+pow(Par[3]*cos(theta), 2)));	//On-shells leaving the positive energy range
	Stops[l+2] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]-pow(2.*Par[2], 2)+pow(Par[3]*cos(theta), 2)));
	Stops[l+3] = sqrt(4.*pow(Par[3], 4)+8.*pow(Par[3], 2)*Par[4]+4.*pow(Par[4], 2)-pow(Par[1], 4))/pow(256.*pow(Par[3], 4)+512.*pow(Par[3], 2)*Par[4]+256.*pow(Par[4], 2), (long double).25);	//Potiential leaving the positive energy range
	Stops[l+4] = abs((pow(Par[2], 2)*Par[3]*cos(theta)+sqrt((Par[4]+pow(Par[3], 2))*(pow(Par[2], 4)+(Par[4]+pow(Par[3]*sin(theta), 2))*(Par[4]-2.*pow(Par[2], 2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta), 2))));	//On-shell leaving the time-like range
	Stops[l+5] = abs((pow(Par[2], 2)*Par[3]*cos(theta)-sqrt((Par[4]+pow(Par[3], 2))*(pow(Par[2], 4)+(Par[4]+pow(Par[3]*sin(theta), 2))*(Par[4]-2.*pow(Par[2], 2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta), 2))));
	Stops[l+6] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]+pow(Par[3]*cos(theta), 2)));	//Photon point leaving positive energy range. Not sure what photon point
	Stops[l+7] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]+pow(Par[3]*cos(theta), 2)));
	Stops[l+8] = .5*abs(Par[3]*cos(theta)+sqrt(3.*pow(Par[3], 2)+4.*Par[4]+pow(Par[3]*cos(theta), 2)));
	Stops[l+9] = .5*abs(Par[3]*cos(theta)-sqrt(3.*pow(Par[3], 2)+4.*Par[4]+pow(Par[3]*cos(theta), 2)));

	for(i = 0; i < l+10; i++)	//Removes stops that are NaN or bigger than necessary
	{
		if(isnan(Stops[i]))
			Stops[i] = -1;
		else if(isinf(Stops[i]) || Stops[i] > 100)
			Stops[i] = 100;
	}

	mergeSort(Stops, 0, l+9);	//Sort the list of sub-intervals

	i = 0;
	j = 0;
	while(Stops[j] <= 0)	//Skip past negative sub-intervals and form NaN
		j++;
	for(; j < l+10; j++)
	{
		if((i > 0 && Stops[i-1] != Stops[j]) || i == 0)	//Removes duplicates, faster to remove duplicates than to evaluate zero width interval
		{
			Stops[i] = Stops[j];
			i++;
		}
		else if(Stops[j] != Stops[j])
			break;
	}
	Intervals = i;	//Record number of intervals in Stops
	Max = Stops[i-1];

	if(j == 0)
		Intervals = 1;

	a = b = i = 0;
	do
	{
		if((i < Intervals && b+100 < Stops[i] && b-Stops[i-1] > 100) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && 50 < Stops[i]-b && b-Stops[i-1] > 50) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && 10 < Stops[i]-b && b-Stops[i-1] > 10) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && 3 < Stops[i]-b) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}
		else
			b += 3;

		F.null();	//Zero out F for next subinterval

		if(b-a < 1)	//use a higher order when the interval is large
		{
			for(l = 0; l < 9; l++) //Count through points away from center
			{
				x1 = (b+a-Disp37[l]*(b-a))/2.;
				x2 = (b+a+Disp37[l]*(b-a))/2.;

				F += k0_Int(Par, Temp, x1, theta)*pow(x1, 2)*w37[l+1];
				F += k0_Int(Par, Temp, x2, theta)*pow(x2, 2)*w37[l+1];
			}
			F += k0_Int(Par, Temp, (a+b)/2., theta)*pow((a+b)/2., 2)*w37[0];
		}
		else
		{
			for(l = 0; l < 24; l++) //Count through points away from center
			{
				x1 = (b+a-Disp97[l]*(b-a))/2.;
				x2 = (b+a+Disp97[l]*(b-a))/2.;

				F += k0_Int(Par, Temp, x1, theta)*pow(x1, 2)*w97[l+1];
				F += k0_Int(Par, Temp, x2, theta)*pow(x2, 2)*w97[l+1];
			}
			F += k0_Int(Par, Temp, (a+b)/2., theta)*pow((a+b)/2., 2)*w97[0];
		}

		Partial = F*(b-a)/2.;	//Record the subinterval to total of the integral
		Answer += Partial;	//Add the subinterval to total of the integral
		a = b;
	}while(!(Partial == 0) && (i < Intervals || abs(Partial/Answer) >= .0001) && (a <= Max || a <= 20.*sqrt(Par[4]+pow(Par[3], 2)))); //Keep going so long as the last subinterval isn't zero and the intervals haven't been exhausted and the last partial answer for all functions isn't too big compared to the total answer and the highest sub-interval is less than 20E. k bigger than 20E is getting pretty stupid, should be sneaking up on 10^-5 of the answer left

	return(Answer);
}

Elements k0_Int(long double Par[], int Temp, long double k, long double theta)
{
	if(Par[4]+pow(Par[3], 2) < 0)	//Bad data trap and time saver. The point is supposed to zero energy anyways but got evaluated to non-zero
		return(Elements(0, 0, 0, 0));
//9th order Gauss-Legendre integration
	long double Disp9[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.};	//Displacement from center
	long double w9[] = {128./225., (322.+13.*sqrt(70.))/900., (322.-13.*sqrt(70.))/900.};	//Weight
//41st order Gauss-Legendre integration
	long double w41[] = {0.146081133649690427192, 0.144524403989970059064, 0.139887394791073154722, 0.132268938633337461781, 0.121831416053728534195, 0.108797299167148377663, 0.0934444234560338615533, 0.0761001136283793020171, 0.0571344254268572082836, 0.0369537897708524938000, 0.0160172282577743333242};
	long double Disp41[] = {0.145561854160895090937, 0.288021316802401096601, 0.424342120207438783574, 0.551618835887219807059, 0.667138804197412319306, 0.768439963475677908616, 0.853363364583317283647, 0.920099334150400828790, 0.967226838566306294317, 0.993752170620389500260};
//37th order Gauss-Legendre integration
	long double Disp37[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center
	long double w37[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight
//97th order Gauss-Legendre integration
	long double Disp97[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center
	long double w97[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight
	long double a, b;	//Sub-interval limits of integration
	long double x1, x2;	//Abscissa
	long double Max;	//Upper limit of integration

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	long double Boundary_k_k0[] = {.421, 4.85};
	long double Range[] = {-Boundary_k_k0[1], -Boundary_k_k0[0], 0, Boundary_k_k0[0], Boundary_k_k0[1]};	//Number of gamma from center

	Elements F;			//Sum of ordinates*weights
	Elements Answer(0, 0, 0, 0);	//Results to be returned
	Elements Partial;		//Partial sum to determine continuation
	Elements Holder;

	long double zero[12];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[12];	//Imaginary part of poles
	long double Caution = k*Par[3]*cos(theta)/(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta));
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities

	Characterize_k0_Int(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	long double  Stops[Poles*5+6];					//Intervals that are required by integrating near poles

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(!isnan(gamma[i]) && !isnan(zero[i]))	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 5; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Extra subintervals required by poles
				l++;
			}
		else if(!isnan(zero[i]))	//At lease insert the central point of the pole if the width isn't properly measured
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+pow(Par[3], 2))/2.;	//Lower light-like edge
	Stops[l+1] = sqrt(Par[4]+pow(Par[3], 2))/2.-Energy(0, Par[3]/2., -k, theta);	//Upper light-like edge
	Stops[l+2] = -Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+pow(Par[3], 2))/2.;	//Pretty sure this is the negative energy solution of the lower light-like edge
	Stops[l+3] = sqrt(Par[4]+pow(Par[3], 2))/2.+Energy(0, Par[3]/2., -k, theta);	//Pretty sure this is the negative energy solution of the upper light-like edge
	Stops[l+4] = sqrt(Par[4]+pow(Par[3], 2))/2.;					//Upper energy boundary (E/2)
	Stops[l+5] = -sqrt(Par[4]+pow(Par[3], 2))/2.;					//Lower energy boundary (-E/2)
	Stops[l+6] = (Stops[l+2]+Stops[l+3])/2.;

	for(i = 0; i < l+7; i++)
		if(Stops[i] < -sqrt(Par[4]+pow(Par[3], 2))/2.) Stops[i] = -Stops[i];

	mergeSort(Stops, 0, l+6);	//Sort the subintervals

	a = b = -sqrt(Par[4]+pow(Par[3], 2))/2.;	//Lower edge
	Max = sqrt(Par[4]+pow(Par[3], 2))/2.;		//Upper edge before doubling up

	i = 0;
	j = 0;
	while(Stops[j] == Stops[j]+1.)	//Remove subintervals that duplicates
		j++;
	for(; j < l+6; j++)
	{
		if((i > 0 && Stops[i-1] != Stops[j]) || i == 0)	//Remove dublicates
		{
			Stops[i] = Stops[j];
			i++;
		}
	}
	Intervals = i;	//Record the number of intervals


	i = 1;	//The first point should be the lower limit of integration. That's where we start. Next subinterval is what we need to be looking for
	do
	{
		if((i < Intervals && b+100 < Stops[i] && b-Stops[i-1] > 100) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && 50 < Stops[i]-b && b-Stops[i-1] > 50) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && 10 < Stops[i]-b && b-Stops[i-1] > 10) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && 3 < Stops[i]-b) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}

		if(b > Max && a < Max)
			b = Max;	//Be sure E/2 is and sub-interval boundary

		F.null();	//Zero out F for next sub-interval

		//Count through points away from center, use only the least precise quadurature rule needed
		if((abs(a-Caution)<.25 || abs(b-Caution)<.25) && k < 1)	//Within .25 GeV of the double on-shell and relative momentum is less than 1
		{
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp97[j]*(b-a))/2.;
				x2 = (b+a+Disp97[j]*(b-a))/2.;

				F += (Elements(Potential1(Par, x1, k), Interacting_Linear_Trace(Par)*Potential1(Par, x1, k), Interacting_Quad_Trace(Par, x1, k)*Potential1(Par, x1, k), Potential2(Par, x1, k))*Dispersion(Par, Temp, x1, k, theta))*w97[j+1];
				F += (Elements(Potential1(Par, x2, k), Interacting_Linear_Trace(Par)*Potential1(Par, x2, k), Interacting_Quad_Trace(Par, x2, k)*Potential1(Par, x2, k), Potential2(Par, x2, k))*Dispersion(Par, Temp, x2, k, theta))*w97[j+1];
			}
			F += (Elements(Potential1(Par, (a+b)/2., k), Interacting_Linear_Trace(Par)*Potential1(Par, (a+b)/2., k), Interacting_Quad_Trace(Par, (a+b)/2., k)*Potential1(Par, (a+b)/2., k), Potential2(Par, (a+b)/2., k))*Dispersion(Par, Temp, (a+b)/2., k, theta))*w97[0];
		}
		else if(abs(a-Caution)<.5 || abs(b-Caution)<.5 || (.25 < b-a && b-a < 1))	//Within .5 GeV of the double on-shell or the interval is between .25 and 1 GeV (bigger than .25 GeV and 9th order isn't good enough, more than 1 GeV and it doesn't matter)
		{
			for(j = 0; j < 10; j++)
			{
				x1 = (b+a-Disp41[j]*(b-a))/2.;
				x2 = (b+a+Disp41[j]*(b-a))/2.;

				F += (Elements(Potential1(Par, x1, k), Interacting_Linear_Trace(Par)*Potential1(Par, x1, k), Interacting_Quad_Trace(Par, x1, k)*Potential1(Par, x1, k), Potential2(Par, x1, k))*Dispersion(Par, Temp, x1, k, theta))*w41[j+1];
				F += (Elements(Potential1(Par, x2, k), Interacting_Linear_Trace(Par)*Potential1(Par, x2, k), Interacting_Quad_Trace(Par, x2, k)*Potential1(Par, x2, k), Potential2(Par, x2, k))*Dispersion(Par, Temp, x2, k, theta))*w41[j+1];
			}
			F += (Elements(Potential1(Par, (a+b)/2., k), Interacting_Linear_Trace(Par)*Potential1(Par, (a+b)/2., k), Interacting_Quad_Trace(Par, (a+b)/2., k)*Potential1(Par, (a+b)/2., k), Potential2(Par, (a+b)/2., k))*Dispersion(Par, Temp, (a+b)/2., k, theta))*w41[0];
		}
		else
		{
			for(j = 0; j < 2; j++)
			{
				x1 = (b+a-Disp9[j]*(b-a))/2.;
				x2 = (b+a+Disp9[j]*(b-a))/2.;

				F += (Elements(Potential1(Par, x1, k), Interacting_Linear_Trace(Par)*Potential1(Par, x1, k), Interacting_Quad_Trace(Par, x1, k)*Potential1(Par, x1, k), Potential2(Par, x1, k))*Dispersion(Par, Temp, x1, k, theta))*w9[j+1];
				F += (Elements(Potential1(Par, x2, k), Interacting_Linear_Trace(Par)*Potential1(Par, x2, k), Interacting_Quad_Trace(Par, x2, k)*Potential1(Par, x2, k), Potential2(Par, x2, k))*Dispersion(Par, Temp, x2, k, theta))*w9[j+1];
			}
			F += (Elements(Potential1(Par, (a+b)/2., k), Interacting_Linear_Trace(Par)*Potential1(Par, (a+b)/2., k), Interacting_Quad_Trace(Par, (a+b)/2., k)*Potential1(Par, (a+b)/2., k), Potential2(Par, (a+b)/2., k))*Dispersion(Par, Temp, (a+b)/2., k, theta))*w9[0];
		}

		if(a >= Max)	//If you're in this region, the above conditions don't apply
		{
			if(abs(a-Caution)<.5 || abs(b-Caution)<.5 || (.25 < b-a && b-a < 1))	//Within .5 GeV of the double on-shell or the interval is between .25 and 1 GeV (bigger than .25 GeV and 9th order isn't good enough, more than 1 GeV and it doesn't matter)
			{
				for(j = 0; j < 10; j++)
				{
					x1 = -(b+a-Disp41[j]*(b-a))/2.;
					x2 = -(b+a+Disp41[j]*(b-a))/2.;

				F += (Elements(Potential1(Par, x1, k), Interacting_Linear_Trace(Par)*Potential1(Par, x1, k), Interacting_Quad_Trace(Par, x1, k)*Potential1(Par, x1, k), Potential2(Par, x1, k))*Dispersion(Par, Temp, x1, k, theta))*w41[j+1];
				F += (Elements(Potential1(Par, x2, k), Interacting_Linear_Trace(Par)*Potential1(Par, x2, k), Interacting_Quad_Trace(Par, x2, k)*Potential1(Par, x2, k), Potential2(Par, x2, k))*Dispersion(Par, Temp, x2, k, theta))*w41[j+1];
			}
			F += (Elements(Potential1(Par, (a+b)/2., k), Interacting_Linear_Trace(Par)*Potential1(Par, (a+b)/2., k), Interacting_Quad_Trace(Par, (a+b)/2., k)*Potential1(Par, (a+b)/2., k), Potential2(Par, (a+b)/2., k))*Dispersion(Par, Temp, (a+b)/2., k, theta))*w41[0];
			}
			else
			{
				for(j = 0; j < 2; j++)
				{
					x1 = -(b+a-Disp9[j]*(b-a))/2.;
					x2 = -(b+a+Disp9[j]*(b-a))/2.;

				F += (Elements(Potential1(Par, x1, k), Interacting_Linear_Trace(Par)*Potential1(Par, x1, k), Interacting_Quad_Trace(Par, x1, k)*Potential1(Par, x1, k), Potential2(Par, x1, k))*Dispersion(Par, Temp, x1, k, theta))*w9[j+1];
				F += (Elements(Potential1(Par, x2, k), Interacting_Linear_Trace(Par)*Potential1(Par, x2, k), Interacting_Quad_Trace(Par, x2, k)*Potential1(Par, x2, k), Potential2(Par, x2, k))*Dispersion(Par, Temp, x2, k, theta))*w9[j+1];
			}
			F += (Elements(Potential1(Par, (a+b)/2., k), Interacting_Linear_Trace(Par)*Potential1(Par, (a+b)/2., k), Interacting_Quad_Trace(Par, (a+b)/2., k)*Potential1(Par, (a+b)/2., k), Potential2(Par, (a+b)/2., k))*Dispersion(Par, Temp, (a+b)/2., k, theta))*w9[0];
			}
		}

		Partial = F*(b-a)/2.;
		Answer += Partial;		//Add the subinterval to the total
		a = b;
	}while((!(Partial == 0) || a < Max) && (i < Intervals || abs(Partial/Answer) >= .0001));	//Keep going while intervals aren't exhausted and the last interval added is too much

	return(Answer/M_PI);
}

long double Dispersion(long double Par[], int Temp, long double k0, long double k, long double theta)
{
//9th order Gauss-Legendre integration
	long double Disp[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.};	//Displacement from center
	long double w[] = {128./225., (322.+13.*sqrt(70.))/900., (322.-13.*sqrt(70.))/900.};	//Weight
	long double a, b;	//Sub-interval limits of integration
	long double Min;	//Lower limit of integration
	long double Max = 700;	//Upper limit of integration
	long double ParLoc[5] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//Local copy of the parameters as Par[4] corrisponds to s and ParLoc[4] is s'
	long double ImG12 = Imk0_Integrand(Par, k0, k, theta, Temp);		//Holder of the ImG12 that belongs to the other half of G12 that is calculated here

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	long double Boundary_k_k0 = 4.85;
	long double Range[] = {-Boundary_k_k0, 0, Boundary_k_k0};	//Number of gamma from center

	long double F;			//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double Partial;		//Partial results to examine convergance

	long double zero[4];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[4];	//Imaginary part of poles
	int Poles;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities

	Characterize_Dispersion(ParLoc, Temp, k0, k, theta, zero, gamma, Poles);
	long double Stops[Poles*3+7];		//Extra stops to ensure correctness

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(!isnan(gamma[i]))	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 3; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Extra subintervals required by poles
				l++;
			}
		else	//At lease insert the central point of the pole if the width isn't properly measured
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = 4.*pow(k0, 2)-pow(Par[3], 2);	//Both quarks remain energy positive, should be the center of the fermi functions, (1-f-f)=.5
	Stops[l+1] = 4.*pow(k, 2)+4.*pow(k0, 2)+3.*pow(Par[3], 2)+4.*k*Par[3]*cos(theta)-8.*sqrt(pow(k*k0, 2)+pow(k0*Par[3], 2)+k*Par[3]*pow(k0, 2)*cos(theta));	//Light-like quarks
	Stops[l+2] = 4.*pow(k, 2)+4.*pow(k0, 2)+3.*pow(Par[3], 2)-4.*k*Par[3]*cos(theta)+8.*sqrt(pow(k*k0, 2)+pow(k0*Par[3], 2)-k*Par[3]*pow(k0, 2)*cos(theta));
	Stops[l+3] = Par[4];	//Division by zero of dispersion relation
	Stops[l+4] = 4.*(pow(k, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta));
	Stops[l+5] = 4.*(pow(k, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta));

	mergeSort(Stops, 0, l+5);
	Stops[l+6] = Stops[l+5]+100;	//Adds the minimum end point to keep the integration going

	Min = a = b = -pow(Par[3], 2);	//Start from s'=-P^2

	i = 0;
	while(Stops[i] < a)
		i++;
	Intervals = l+5;

	do
	{
		if((i < Intervals && b+100 < Stops[i] && b-Stops[i-1] > 100) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if((i < Intervals && 50 < Stops[i]-b && b-Stops[i-1] > 50) || Stops[Intervals-1] < a-50)
			b += 50;
		else if((i < Intervals && 10 < Stops[i]-b && b-Stops[i-1] > 10) || Stops[Intervals-1] < a-10)
			b += 10;
		else if((i < Intervals && 3 < Stops[i]-b && b-Stops[i-1] > 3) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			while(abs(a-b)<1e-14)
			{
				if(Stops[i]-b < 3 && i < Intervals)
				{
					b = Stops[i];
					i++;
				}
				else
					b += 3;
			}
		}

		F = 0;	//Zero out F for next sub-interval

		/*if(chrono::duration_cast<chrono::seconds>(chrono::system_clock::now()-Start_Time).count() > Allotment)
		{
			throw(chrono::duration_cast<chrono::seconds>(chrono::system_clock::now()-Start_Time).count());
		}*/

		for(l = 0; l < 2; l++) //Count through points away from center
		{
			ParLoc[4] = (b+a-Disp[l]*(b-a))/2.;
			F += (Imk0_Integrand(ParLoc, k0, k, theta, Temp)-ImG12)/(ParLoc[4]-Par[4])*w[l+1];
			ParLoc[4] = (b+a+Disp[l]*(b-a))/2.;
			F += (Imk0_Integrand(ParLoc, k0, k, theta, Temp)-ImG12)/(ParLoc[4]-Par[4])*w[l+1];
		}
		ParLoc[4] = (b+a)/2.;
		F += (Imk0_Integrand(ParLoc, k0, k, theta, Temp)-ImG12)/(ParLoc[4]-Par[4])*w[0];

		Partial = F*(b-a)/2.;
		Answer += Partial;		//Add the subinterval to the total
		a = b;
	}while((a < Max && i < Intervals) || Partial/Answer > 1e-6);	//Keep going while intervals aren't exhausted and upper limit of integration not excceeded or until convergance

	if(abs(ImG12) >= 1e-20)
		return((Answer+ImG12*log(abs((a-Par[4])/(Par[4]-Min))))/M_PI);
	return(Answer/M_PI);
}

void Characterize_k_Int(long double Par[], int Temp, long double theta, long double zero[], long double gamma[], int &Poles)
{
	long double holder;	//Holder for bubble sort at the end
	int i, j;		//Counter

	Poles = 2;	//Two hard coded poles
	zero[0] = .5*Par[3]*abs(cos(theta));	//Supposedly this is near intersection of 2 on-shells. I don't recongize it, lines 430-435 does that
	//gamma[0] = .05;
	zero[1] = Par[2];
	//gamma[1] = Par[1];

	if(Par[4]-pow(2.*Par[2], 2) > 0.)
	{
		zero[Poles] = .5*sqrt((Par[4]-pow(2.*Par[2], 2))*(Par[4]+pow(Par[3], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)));	//relative 3-momentum for which both quarks are on-shell
		//gamma[Poles] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp);	//Twice the self-energy of the quarks. Should be sum of self-energy of the two quarks as I don't think the self-energy of both quarks are equal to each other for P!=0.
		Poles++;
	}

	zero[Poles] = Par[2];	//Use Newtons's method to find intersection of features of interest. This one is the positive pole of the potential and one of the peaks of the propagator, I think its the negative pole of the anti-aligned quark.
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;	//Try again with a different seed in case the first missed
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		//gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp)+sqrt(complex<long double>(pow(2.*zero[Poles], 2), pow(Par[2], 2))).imag();
		Poles++;
	}

	for(i = 0; i < Poles; i++)	//Move negative results to positive result. Should probably drop it, but might have been missed or otherwise caught a feature that is interesting, just with the wrong sign
	{
		zero[i] = abs(zero[i]);
		//gamma[i] = abs(gamma[i]);
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
				/*holder = gamma[j+1];
				gamma[j+1] = gamma[j];
				gamma[j] = holder;*/
			}
		}
	}

	i = 0;
	for(j = 0; j < Poles; j++)
	{
		if(zero[j] > 1000)	//Don't bother with points beyond 1 TeV. All integrals should end well before then
			break;
		if(((i > 0 && zero[i-1] != zero[j]) || i == 0) && !isnan(zero[j]))	//Remove duplicates and NaN
		{
			//zero[i] = zero[j];
			gamma[i] = gamma[j];
			i++;
		}
	}
	Poles = i;

	return;
}

bool Newton_Method_k(long double& k, long double s, long double P, long double theta, long double M, long double Lambda, long double(*V)(long double, long double), long double(*k0)(long double, long double, long double, long double, long double))	//Returns the k-intesection of a potiential and on-shell peak
{
	long double newk;
	const long double h = 1e-4;	//Size of finite difference
	bool Success = true;
	int i = 0;

	newk = k - 2.*h*(V(k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(k-h, Lambda)+V(k+h, Lambda)-k0(s, P, k+h, theta, M));	//First iteration of Newton's method with finite differences for derivatives

	while(abs(1.-newk/k) > 1e-5 && i <= 10)	//Allow up to 12 iteration steps, but stop early if last step was small
	{
		k = newk;
		newk = k - 2.*h*(V(k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(k-h, Lambda)+V(k+h, Lambda)-k0(s, P, k+h, theta, M));
		i++;
	}

	if(abs(1.-newk/k) > 1e-2 || newk < 0 || newk != newk)	//Soundness of value (positive and not NaN) and degree of improvement on last interation (last iteration should have been small)
		Success = false;
	if(abs(1.-V(newk, Lambda)/k0(s, P, newk, theta, M)) > 1e-2)	//Closeness to solution
		Success = false;

	k = newk;

	return(Success);	//Note success or failure of solution find
}

long double V_Plus(long double k, long double Lambda)
{
#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	return(k);
#elif VERSION == 22
	return(k);
#elif VERSION == 24
	return(k);
#elif VERSION == 42
	return(.5*sqrt(complex<long double>(4.*(pow(k, 2)), pow(Lambda, 2))).real());
#endif
}

long double V_Minus(long double k, long double Lambda)
{
#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	return(-k);
#elif VERSION == 22
	return(-k);
#elif VERSION == 24
	return(-k);
#elif VERSION == 42
	return(-.5*sqrt(complex<long double>(4.*(pow(k, 2)), pow(Lambda, 2))).real());
#endif
}

long double Emm(long double s, long double P, long double k, long double theta, long double M)	//peak of the vacuum on-shells
{
	return(sqrt(s+pow(P, 2))/2.-Energy(M, P/2., -k, theta));
}

long double Epm(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P, 2))/2.+Energy(M, P/2., -k, theta));
}

long double mEmp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P, 2))/2.-Energy(M, P/2., k, theta));
}

long double Emp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P, 2))/2.-Energy(M, P/2., k, theta));
}

long double mEpp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P, 2))/2.+Energy(M, P/2., k, theta));
}

long double Epp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P, 2))/2.+Energy(M, P/2., k, theta));
}

long double Upper_Bound(long double s, long double P, long double k, long double theta, long double M)	//Vacuum boundaries
{
	return(sqrt(s+pow(P, 2))/2.-Energy(0, P/2., -k, theta));
}

long double Lower_Bound(long double s, long double P, long double k, long double theta, long double M)
{
	return(Energy(0, P/2., k, theta)-sqrt(s+pow(P, 2))/2.);
}

void Characterize_k0_Int(long double Par[], int Temp, long double k, long double theta, long double zero[10], long double gamma[10], int &Poles)
{
	long double Lower, Upper;	//Limits of integration in k0_Int, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(true)//Temp != 0)
	{
		Lower = -sqrt(Par[4]+pow(Par[3], 2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3], 2))/2.;	//Integrate from -E/2 to E/2
	}
	else
	{
		Lower = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+pow(Par[3], 2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3], 2))/2.-Energy(0, Par[3]/2., -k, theta);
	}

#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
	gamma[0] = Par[1];
	gamma[1] = Par[1];
#elif VERSION == 22
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
	gamma[0] = k-.5*sqrt(abs(pow(2.*k, 2)-pow(Par[1], 2)));
	gamma[1] = k-.5*sqrt(abs(pow(2.*k, 2)-pow(Par[1], 2)));
#elif VERSION == 24
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
	gamma[0] = k-sqrt(abs(pow(k, 2)-pow(Par[1], 2)/2));
	gamma[1] = k-sqrt(abs(pow(k, 2)-pow(Par[1], 2)/2));
#elif VERSION == 42
	zero[0] = .5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).real();	//Potential poles, I know exactly where these are at.
	zero[1] = -.5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).real();
	gamma[0] = abs(.5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).imag());
	gamma[1] = abs(-.5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).imag());
#endif

	zero[2] = .5*(sqrt(Par[4]+pow(Par[3], 2))-real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));	//Exact vacuum poles
	zero[3] = .5*(sqrt(Par[4]+pow(Par[3], 2))+real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	zero[4] = .5*(-sqrt(Par[4]+pow(Par[3], 2))-real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	zero[5] = .5*(-sqrt(Par[4]+pow(Par[3], 2))+real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	zero[6] = k*Par[3]*cos(theta)/(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta));
	zero[7] = cos(theta)*Par[3]/2.*sqrt((Par[4]-pow(2.*Par[2], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)));

	if(Temp != 0)	//media estimate
	{
		zero[2] = Newton_Method_k0(zero[2], Par, k, theta, Temp, Imk0_Integrand);
		zero[3] = Newton_Method_k0(zero[3], Par, k, theta, Temp, Imk0_Integrand);
		zero[4] = Newton_Method_k0(zero[4], Par, k, theta, Temp, Imk0_Integrand);
		zero[5] = Newton_Method_k0(zero[5], Par, k, theta, Temp, Imk0_Integrand);

		gamma[2] = omega_Width(zero[2], Par, k, theta, Temp, Imk0_Integrand);
		gamma[3] = omega_Width(zero[3], Par, k, theta, Temp, Imk0_Integrand);
		gamma[4] = omega_Width(zero[4], Par, k, theta, Temp, Imk0_Integrand);
		gamma[5] = omega_Width(zero[5], Par, k, theta, Temp, Imk0_Integrand);
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = gamma[3] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
		gamma[4] = gamma[5] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	}

	gamma[6] = ImSelf_Energy(Par[2], sqrt(Par[4]+pow(Par[3], 2))/2.+zero[6], Energy(0, Par[3], k, theta), Temp)+ImSelf_Energy(Par[2], sqrt(Par[4]+pow(Par[3], 2))/2.-zero[6], Energy(0, Par[3], -k, theta), Temp)-GAMMA;
	if(!isnan(zero[7]))
		gamma[7] = ImSelf_Energy(Par[2], sqrt(Par[4]+pow(Par[3], 2))/2.+zero[7], Energy(0, Par[3], k, theta), Temp)+ImSelf_Energy(Par[2], sqrt(Par[4]+pow(Par[3], 2))/2.-zero[7], Energy(0, Par[3], -k, theta), Temp);

	j = i = 0;
	while(j < 8)
	{
		while(isnan(zero[j]) && j < 8)	//remove nan
			j++;

		if(zero[j] < Lower)	//move pole from below Lower to above Upper
			zero[i] = abs(zero[j]);
		else
			zero[i] = zero[j];

		gamma[i] = abs(gamma[j]);
		i++;
		j++;
	}

	Poles = i-1;

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

	return;
}

long double Newton_Method_k0(long double k0, long double Par[], long double k, long double theta, int Temp, long double (*Folding)(long double[], long double, long double, long double, int))	//Newton's method for finding poles of f by looking for zeros of 1/f, much more stable to the point of absolute confidence
{
	long double new_k0;
	long double h;	//Finite difference
	long double Exit;
	int i = 0;
	long double Danger[] = {.5*(sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4)))))),.5*(-sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))))};

	if(abs(k0-Danger[0]) < 1e-6 || abs(k0-Danger[1]) < 1e-6)
	{
		h = 1e-10;
		Exit = 1e-9;
	}
	else
	{
		h = 1e-4;
		Exit = 1e-5;
	}

	new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp)-1./Folding(Par, k0-h, k, theta, Temp))/((1./Folding(Par, k0-h, k, theta, Temp)-2./Folding(Par, k0, k, theta, Temp)+1./Folding(Par, k0+h, k, theta, Temp)));	//First iteration of Netwon's method using finite differences

	while(abs(1.-new_k0/k0) > Exit && i < 10 && !isnan(new_k0))	//Allow up to 12 iterations to find the pole
	{
		k0 = new_k0;
		new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp)-1./Folding(Par, k0-h, k, theta, Temp))/((1./Folding(Par, k0-h, k, theta, Temp)-2./Folding(Par, k0, k, theta, Temp)+1./Folding(Par, k0+h, k, theta, Temp)));
		i++;
	}

	return(new_k0);
}

long double omega_Width(long double zero, long double Par[], long double k, long double theta, int Temp, long double (*Folding)(long double[], long double, long double, long double, int))	//Breit-Wigner width of the peak
{
	return(sqrt(abs(2e-10*Folding(Par, zero, k, theta, Temp)/(Folding(Par, zero-1e-5, k, theta, Temp)-2.*Folding(Par, zero, k, theta, Temp)+Folding(Par, zero+1e-5, k, theta, Temp)))));
}

void Characterize_Dispersion(long double Par[], int Temp, long double k0, long double k, long double theta, long double zero[], long double gamma[], int &Poles)
{
	zero[0] = 4.*(pow(k, 2)+pow(k0, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta)-sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)+4.*k*pow(k0, 2)*Par[3]*cos(theta)));
	zero[1] = 4.*(pow(k, 2)+pow(k0, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta)+sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)+4.*k*pow(k0, 2)*Par[3]*cos(theta))); //Both of the possible on-shell s using positive k^mu
	zero[2] = 4.*(pow(k, 2)+pow(k0, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta)-sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)-4.*k*pow(k0, 2)*Par[3]*cos(theta)));
	zero[3] = 4.*(pow(k, 2)+pow(k0, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta)+sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)-4.*k*pow(k0, 2)*Par[3]*cos(theta))); //Both of the possible on-shell s using negative k^mu

	//Calcluate and record the widths of the peaks
	Par[4] = zero[0];
	gamma[0] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[1];
	gamma[1] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[2];
	gamma[2] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[3];
	gamma[3] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));

	Poles = 4;
}

long double sp_Width(long double Par[], long double k0, long double k, long double theta, int Temp, long double (*Folding)(long double[], long double, long double, long double, int))	//Breit-Wigner width of the peak
{
	long double y[3];
	y[1] = Folding(Par, k0, k, theta, Temp);
	Par[4] -= 1e-5;
	y[0] = Folding(Par, k0, k, theta, Temp);
	Par[4] += 2e-5;
	y[2] = Folding(Par, k0, k, theta, Temp);
	return(sqrt(abs(2e-10*y[1]/(y[0]-2.*y[1]+y[2]))));
}

void ImSelf_Energy(long double M, long double omega[], long double k[], int Temp, long double Results[])	//Single quark self energy for both quarks
{
	static long double omega0[2];		//Location of central peak
	static long double Sigma[2];		//Amplitude of energy dependance
	static long double a[2], b[2];	//Slope of exponential decrease to left and right
	static long double knee[2];		//Interval to change from left to right side of peak
	static long double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	static long double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	if(pow(omega[0], 2)>=pow(k[0], 2) && omega[0] >= 0)	//Vacuum width
		Results[0] = sqrt(pow(omega[0], 2)-pow(k[0], 2))*GAMMA;
	else
		Results[0] = 0;
	if(pow(omega[1], 2)>=pow(k[1], 2) && omega[1] >= 0)
		Results[1] = sqrt(pow(omega[1], 2)-pow(k[1], 2))*GAMMA;
	else
		Results[1] = 0;

	if(Temp == 0)
		return;

	if(k[0] != k_old[0] || k[1] != k_old[1])	//If either of the relative momenta have been altered
	{
		k_old[0] = k[0];
		k_old[1] = k[1];
		switch(Temp)
		{
			/*case 1://194MeV, variation of T=194 MeV self-energy not used but more similar to the provisioned self-energy
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(pow(k[0], 2)+pow(1.75236, 2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1], 2)+pow(1.75236, 2))+.0187484;
				a[0] = 4.689/(pow(k[0], 2)+pow(1.18, 2))+4.59495;
				a[1] = 4.689/(pow(k[1], 2)+pow(1.18, 2))+4.59495;
				b[0] = -70400/(pow(k[0]+20, 2)+pow(130, 2))+6.24;
				b[1] = -70400/(pow(k[1]+20, 2)+pow(130, 2))+6.24;
				omega0[0] = sqrt(pow(1.51443+Shift, 2)+pow(k[0], 2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift, 2)+pow(k[1], 2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1., (long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1., (long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(pow(k[0], 2)+pow(1.75236, 2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1], 2)+pow(1.75236, 2))+.0187484;
				a[0] = 12.5349/(pow(k[0], 2)+pow(1.63711, 2))+5.026;
				a[1] = 12.5349/(pow(k[1], 2)+pow(1.63711, 2))+5.026;
				b[0] = -291.579/(pow(k[0]+15.2519, 2)+pow(.0614821, 2))+3.36681;
				b[1] = -291.579/(pow(k[1]+15.2519, 2)+pow(.0614821, 2))+3.36681;
				omega0[0] = sqrt(pow(1.51443+Shift, 2)+pow(k[0], 2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift, 2)+pow(k[1], 2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1., (long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1., (long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;
			case 2://285MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .625855/sqrt(pow(k[0], 2)+pow(1.8429, 2))+.0249334;
				Sigma[1] = .625855/sqrt(pow(k[1], 2)+pow(1.8429, 2))+.0249334;
				a[0] = 3.3971/(pow(k[0], 2)+pow(1.01744, 2))+3.99561;
				a[1] = 3.3971/(pow(k[1], 2)+pow(1.01744, 2))+3.99561;
				b[0] = -65187.5/(pow(k[0]+3.11711, 2)+pow(101.697, 2))+8.15532;
				b[1] = -65187.5/(pow(k[1]+3.11711, 2)+pow(101.697, 2))+8.15532;
				omega0[0] = sqrt(pow(1.5065+Shift, 2)+pow(k[0], 2))+.209135;
				omega0[1] = sqrt(pow(1.5065+Shift, 2)+pow(k[1], 2))+.209135;
				knee[0] = 3.1568*pow(k[0]+1., (long double)-.624827)+.197004*(tanh((k[0]-27.1743)/10.0192)+1);
				knee[1] = 3.1568*pow(k[1]+1., (long double)-.624827)+.197004*(tanh((k[1]-27.1743)/10.0192)+1);
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .587509/sqrt(pow(k[0], 2)+pow(1.84447, 2))+.0309251;
				Sigma[1] = .587509/sqrt(pow(k[1], 2)+pow(1.84447, 2))+.0309251;
				a[0] = 2.44943/(pow(k[0], 2)+pow(.887313, 2))+3.32859;
				a[1] = 2.44943/(pow(k[1], 2)+pow(.887313, 2))+3.32859;
				b[0] = -4439.38/(pow(k[0]-7.23198, 2)+pow(38.9387, 2))+4.55531;
				b[1] = -4439.38/(pow(k[1]-7.23198, 2)+pow(38.9387, 2))+4.55531;
				omega0[0] = sqrt(pow(1.47725+Shift, 2)+pow(k[0], 2))+.219181;
				omega0[1] = sqrt(pow(1.47725+Shift, 2)+pow(k[1], 2))+.219181;
				knee[0] = 3.28564*pow(k[0]+1., (long double)-.721321)+.330483*(tanh((k[0]-22.9096)/10.7139)+1);
				knee[1] = 3.28564*pow(k[1]+1., (long double)-.721321)+.330483*(tanh((k[1]-22.9096)/10.7139)+1);
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .459303/sqrt(pow(k[0], 2)+pow(1.84321, 2))+.0386564;
				Sigma[1] = .459303/sqrt(pow(k[1], 2)+pow(1.84321, 2))+.0386564;
				a[0] = 1.79149/(pow(k[0], 2)+pow(.764836, 2))+2.66209;
				a[1] = 1.79149/(pow(k[1], 2)+pow(.764836, 2))+2.66209;
				b[0] = -1856.16/(pow(k[0]-8.69519, 2)+pow(26.3551, 2))+3.94631;
				b[1] = -1856.16/(pow(k[1]-8.69519, 2)+pow(26.3551, 2))+3.94631;
				omega0[0] = sqrt(pow(1.45428+Shift, 2)+pow(k[0], 2))+.197493;
				omega0[1] = sqrt(pow(1.45428+Shift, 2)+pow(k[1], 2))+.197493;
				knee[0] = 3.06296*pow(k[0]+1., (long double)-.917081)+.394833*(tanh((k[0]-19.5932)/12.0494)+1);
				knee[1] = 3.06296*pow(k[1]+1., (long double)-.917081)+.394833*(tanh((k[1]-19.5932)/12.0494)+1);
				break;
			default:
				omega0[0]=omega0[1]=1.74727;
				Sigma[0]=Sigma[1]=.344006;
				a[0]=a[1]=9.70298;
				b[0]=b[1]=2.11338;
				knee[0]=knee[1]=3.78966;
		}
	}

	long double ImSigma[2];	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
	if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] < -4.)
		ImSigma[0] = a[0]*(omega[0]-omega0[0]+knee[0]/sqrt(a[0]*b[0]));
	else if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] > 4.)
		ImSigma[0] = b[0]*(omega0[0]-omega[0]+knee[0]/sqrt(a[0]*b[0]));
	else	//Lost of precision having been circumvented, the actual value
		ImSigma[0] = -.5*((a[0]-b[0])*omega0[0]-((a[0]+b[0])*knee[0])/sqrt(a[0]*b[0]))+(a[0]-b[0])*omega[0]/2-sqrt(pow(((a[0]+b[0])/2.)*(omega[0]-omega0[0]+((a[0]-b[0])*knee[0])/(sqrt(a[0]*b[0])*(a[0]+b[0]))), 2)+pow(knee[0], 2));

	if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] < -4.)
		ImSigma[1] = a[1]*(omega[1]-omega0[1]+knee[1]/sqrt(a[1]*b[1]));
	else if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] > 4.)
		ImSigma[1] = b[1]*(omega0[1]-omega[1]+knee[1]/sqrt(a[1]*b[1]));
	else
		ImSigma[1] = -.5*((a[1]-b[1])*omega0[1]-((a[1]+b[1])*knee[1])/sqrt(a[1]*b[1]))+(a[1]-b[1])*omega[1]/2-sqrt(pow(((a[1]+b[1])/2.)*(omega[1]-omega0[1]+((a[1]-b[1])*knee[1])/(sqrt(a[1]*b[1])*(a[1]+b[1]))), 2)+pow(knee[1], 2));

#ifdef HALF
	Results[0] += -M*Sigma[0]*exp(ImSigma[0]);	//ImSigma from the in-medium
	Results[1] += -M*Sigma[1]*exp(ImSigma[1]);
#else
	Results[0] += -2.*M*Sigma[0]*exp(ImSigma[0]);
	Results[1] += -2.*M*Sigma[1]*exp(ImSigma[1]);
#endif
	return;
}

long double ImSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{

	long double omega0;	//location of central peak
	long double Sigma;	//size of energy dependance
	long double a, b;	//slope of exponential decrease to left and right
	long double knee;	//space to change from left to right side of peak
	long double M_T, Shift=0;
	long double answer;

	if(pow(omega, 2)>=pow(k, 2) && omega >= 0)
		answer = sqrt(pow(omega, 2)-pow(k, 2))*GAMMA;
	else
		answer = 0;

	if(Temp == 0)
		return(answer);

	switch(Temp)
	{
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k, 2)+pow(1.75236, 2))+.0187484;
			a = 4.689/(pow(k, 2)+pow(1.18, 2))+4.59495;
			b = -70400/(pow(k+20, 2)+pow(130, 2))+6.24;
			omega0 = sqrt(pow(1.51443+Shift, 2)+pow(k, 2))+.232841;
			knee = 3.78956*pow(k+1., (long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;*/
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k, 2)+pow(1.75236, 2))+.0187484;
			a = 12.5349/(pow(k, 2)+pow(1.63711, 2))+5.026;
			b = -291.579/(pow(k+15.2519, 2)+pow(.0614821, 2))+3.36681;
			omega0 = sqrt(pow(1.51443+Shift, 2)+pow(k, 2))+.232841;
			knee = 3.78956*pow(k+1., (long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .625855/sqrt(pow(k, 2)+pow(1.8429, 2))+.0249334;
			a = 3.3971/(pow(k, 2)+pow(1.01744, 2))+3.99561;
			b = -65187.5/(pow(k+3.11711, 2)+pow(101.697, 2))+8.15532;
			omega0 = sqrt(pow(1.5065+Shift, 2)+pow(k, 2))+.209135;
			knee = 3.1568*pow(k+1., (long double)-.624827)+.197004*(tanh((k-27.1743)/10.0192)+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .587509/sqrt(pow(k, 2)+pow(1.84447, 2))+.0309251;
			a = 2.44943/(pow(k, 2)+pow(.887313, 2))+3.32859;
			b = -4439.38/(pow(k-7.23198, 2)+pow(38.9387, 2))+4.55531;
			omega0 = sqrt(pow(1.47725+Shift, 2)+pow(k, 2))+.219181;
			knee = 3.28564*pow(k+1., (long double)-.721321)+.330483*(tanh((k-22.9096)/10.7139)+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .459303/sqrt(pow(k, 2)+pow(1.84321, 2))+.0386564;
			a = 1.79149/(pow(k, 2)+pow(.764836, 2))+2.66209;
			b = -1856.16/(pow(k-8.69519, 2)+pow(26.3551, 2))+3.94631;
			omega0 = sqrt(pow(1.45428+Shift, 2)+pow(k, 2))+.197493;
			knee = 3.06296*pow(k+1., (long double)-.917081)+.394833*(tanh((k-19.5932)/12.0494)+1);
			break;
		case 5://40MeV
			M_T = 1.8;
			Shift = M-M_T;
			Sigma = .00386564;
			a = 6.2;
			b = 2.8;
			omega0 = sqrt(pow(1.53+Shift, 2)+pow(k, 2));
			knee = .56;
			break;
		default:
			omega0=1.74727;
			Sigma=.344006;
			a=9.70298;
			b=2.11338;
			knee=3.78966;
			break;
	}

	long double ImSigma;	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
	if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee < -4.)
		ImSigma = a*(omega-omega0+knee/sqrt(a*b));
	else if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee > 4.)
		ImSigma = b*(omega0-omega+knee/sqrt(a*b));
	else
		ImSigma = -.5*((a-b)*omega0-((a+b)*knee)/sqrt(a*b))+(a-b)*omega/2-sqrt(pow(((a+b)/2.)*(omega-omega0+((a-b)*knee)/(sqrt(a*b)*(a+b))), 2)+pow(knee, 2));

#ifdef HALF
	answer += -M*Sigma*exp(ImSigma);
#else
	answer += -2.*M*Sigma*exp(ImSigma);
#endif

	return(answer);
}

void ReSelf_Energy(long double M, long double omega[], long double k[], int Temp, long double Results[])	//Single quark self energy
{
	static long double Sigma[2];		//Strength
	static long double x0[2], x1[2];	//Centrality markers
	static long double gamma[2];		//Width
	static long double Shift, M_T;
	static long double k_old[2];		//Note on validity of k

	if(Temp == 0 || Temp == 5)
	{
		Results[0] = 0;
		Results[1] = 0;
		return;
	}

	if(k[0] != k_old[0] || k[1] != k_old[1])
	{
		k_old[0] = k[0];
		k_old[1] = k[1];
		switch(Temp)
		{
			/*case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .257498/sqrt(pow(k[0], 2)+pow(1.33201, 2))+.00762638;
				Sigma[1] = .257498/sqrt(pow(k[1], 2)+pow(1.33201, 2))+.00762638;
				x0[0] = sqrt(pow(k[0], 2)+pow(1.54778+Shift, 2))+.276509;
				x0[1] = sqrt(pow(k[1], 2)+pow(1.54778+Shift, 2))+.276509;
				x1[0] = sqrt(pow(k[0], 2)+pow(1.49799+Shift, 2))+.246719;
				x1[1] = sqrt(pow(k[1], 2)+pow(1.49799+Shift, 2))+.246719;
				gamma[0] = .658734/sqrt(pow(k[0], 2)+pow(3.35217, 2))+.0815109;
				gamma[1] = .658734/sqrt(pow(k[1], 2)+pow(3.35217, 2))+.0815109;
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571/sqrt(pow(k[0], 2)+pow(1.17821, 2))+.00762638;
				Sigma[1] = .212571/sqrt(pow(k[1], 2)+pow(1.17821, 2))+.00762638;
				x0[0] = sqrt(pow(k[0], 2)+pow(1.57536+Shift, 2))+.259147;
				x0[1] = sqrt(pow(k[1], 2)+pow(1.57536+Shift, 2))+.259147;
				x1[0] = sqrt(pow(k[0], 2)+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(pow(k[1], 2)+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .336699/sqrt(pow(k[0], 2)+pow(1.87956, 2))+.0651449;
				gamma[1] = .336699/sqrt(pow(k[1], 2)+pow(1.87956, 2))+.0651449;
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972/sqrt(pow(k[0], 2)+pow(1.41483, 2))+.0101423;
				Sigma[1] = .307972/sqrt(pow(k[1], 2)+pow(1.41483, 2))+.0101423;
				x0[0] = sqrt(pow(k[0], 2)+pow(1.56476+Shift, 2))+.251031;
				x0[1] = sqrt(pow(k[1], 2)+pow(1.56476+Shift, 2))+.251031;
				x1[0] = sqrt(pow(k[0], 2)+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(pow(k[1], 2)+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .550628/sqrt(pow(k[0], 2)+pow(2.43968, 2))+.0981269;
				gamma[1] = .550628/sqrt(pow(k[1], 2)+pow(2.43968, 2))+.0981269;
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131/sqrt(pow(k[0], 2)+pow(1.43308, 2))+.0125796;
				Sigma[1] = .339131/sqrt(pow(k[1], 2)+pow(1.43308, 2))+.0125796;
				x0[0] = sqrt(pow(k[0], 2)+pow(1.55034+Shift, 2))+.257788;
				x0[1] = sqrt(pow(k[1], 2)+pow(1.55034+Shift, 2))+.257788;
				x1[0] = sqrt(pow(k[0], 2)+pow(1.46999+Shift, 2))+.231821;
				x1[1] = sqrt(pow(k[1], 2)+pow(1.46999+Shift, 2))+.231821;
				gamma[0] = .615278/sqrt(pow(k[0], 2)+pow(2.22298, 2))+.143376;
				gamma[1] = .615278/sqrt(pow(k[1], 2)+pow(2.22298, 2))+.143376;
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841/sqrt(pow(k[0], 2)+pow(1.42911, 2))+.0157245;
				Sigma[1] = .304841/sqrt(pow(k[1], 2)+pow(1.42911, 2))+.0157245;
				x0[0] = sqrt(pow(k[0], 2)+pow(1.55511+Shift, 2))+.231105;
				x0[1] = sqrt(pow(k[1], 2)+pow(1.55511+Shift, 2))+.231105;
				x1[0] = sqrt(pow(k[0], 2)+pow(1.44714+Shift, 2))+.20956;
				x1[1] = sqrt(pow(k[1], 2)+pow(1.44714+Shift, 2))+.20956;
				gamma[0] = .862629/sqrt(pow(k[0], 2)+pow(2.67193, 2))+.189598;
				gamma[1] = .862629/sqrt(pow(k[1], 2)+pow(2.67193, 2))+.189598;
				break;
			default:
				Sigma[0] = Sigma[1] = .188045;
				x0[0] = x0[1] = 1.83451;
				x1[0] = x1[1] = 1.72447;
				gamma[0] = gamma[1] = .244282;
		}
	}

#ifdef HALF
	Results[0] = Sigma[0]*(omega[0]-x0[0])/(pow(omega[0]-x1[0], 2)+gamma[0])/2.;
	Results[1] = Sigma[1]*(omega[1]-x0[1])/(pow(omega[1]-x1[1], 2)+gamma[1])/2.;
#else
	Results[0] = Sigma[0]*(omega[0]-x0[0])/(pow(omega[0]-x1[0], 2)+gamma[0]);
	Results[1] = Sigma[1]*(omega[1]-x0[1])/(pow(omega[1]-x1[1], 2)+gamma[1]);
#endif
	return;
}

void Self_Energy(long double M, long double omega[], long double k[], int Temp, long double ImSelf[], long double ReSelf[])	//Single quark self energy for both quarks. This one has both imaginary and real parts. It is a simple Breit-Wigner peak and simplier than the other provisioned version
{
	static long double omega0[2];	//location of central peak
	static long double Sigma[2];	//size of energy dependance
	static long double gamma[2];	//space to change from left to right side of peak
	static long double k_old[2];

	if(pow(omega[0], 2)>=pow(k[0], 2))
		ImSelf[0] = sqrt(pow(omega[0], 2)-pow(k[0], 2))*GAMMA;
	else
		ImSelf[0] = 0;
	if(pow(omega[1], 2)>=pow(k[1], 2))
		ImSelf[1] = sqrt(pow(omega[1], 2)-pow(k[1], 2))*GAMMA;
	else
		ImSelf[1] = 0;
	ReSelf[0] = ReSelf[1] = 0;

	if(Temp == 0)
		return;

	if(k[0] != k_old[0] || k[1] != k_old[1])
	{
		k_old[0] = k[0];
		k_old[1] = k[1];
		switch(Temp)
		{
			case 1://194MeV
				Sigma[0] = .840172/sqrt(pow(k[0], 2)+pow(1.45603, 2))+.021257;
				Sigma[1] = .840172/sqrt(pow(k[1], 2)+pow(1.45603, 2))+.021257;
				//omega0[0] = sqrt(pow(M, 2)+pow(k[0], 2));
				//omega0[1] = sqrt(pow(M, 2)+pow(k[1], 2));
				omega0[0] = sqrt(pow(1.99829, 2)+pow(k[0], 2));
				omega0[1] = sqrt(pow(1.99829, 2)+pow(k[1], 2));
				gamma[0] = 1.05035*pow(k[0]+1.3891, (long double)-1.3891)+.01;
				gamma[1] = 1.05035*pow(k[1]+1.3891, (long double)-1.3891)+.01;
				break;
			case 2://285MeV
				Sigma[0] = 1.05337/sqrt(pow(k[0], 2)+pow(1.50861, 2))+.0282696;
				Sigma[1] = 1.05337/sqrt(pow(k[1], 2)+pow(1.50861, 2))+.0282696;
				//omega0[0] = sqrt(pow(M, 2)+pow(k[0], 2));
				//omega0[1] = sqrt(pow(M, 2)+pow(k[1], 2));
				omega0[0] = sqrt(pow(1.97732, 2)+pow(k[0], 2));
				omega0[1] = sqrt(pow(1.97732, 2)+pow(k[1], 2));
				gamma[0] = 1.4624*pow(k[0]+2.64, (long double)-1.41048)+.01;
				gamma[1] = 1.4624*pow(k[1]+2.64, (long double)-1.41048)+.01;
				break;
			case 3://320MeV
				Sigma[0] = 1.14064/sqrt(pow(k[0], 2)+pow(1.54999, 2))+.0350631;
				Sigma[1] = 1.14064/sqrt(pow(k[1], 2)+pow(1.54999, 2))+.0350631;
				//omega0[0] = sqrt(pow(M, 2)+pow(k[0], 2));
				//omega0[1] = sqrt(pow(M, 2)+pow(k[1], 2));
				omega0[0] = sqrt(pow(1.96823, 2)+pow(k[0], 2));
				omega0[1] = sqrt(pow(1.96823, 2)+pow(k[1], 2));
				gamma[0] = 2.07102*pow(k[0]+3.037, (long double)-1.46076)+.01;
				gamma[1] = 2.07102*pow(k[1]+3.037, (long double)-1.46076)+.01;
				break;
			case 4://400MeV
				Sigma[0] = 1.06073/sqrt(pow(k[0], 2)+pow(1.64912, 2))+.0438288;
				Sigma[1] = 1.06073/sqrt(pow(k[1], 2)+pow(1.64912, 2))+.0438288;
				//omega0[0] = sqrt(pow(M, 2)+pow(k[0], 2));
				//omega0[1] = sqrt(pow(M, 2)+pow(k[1], 2));
				omega0[0] = sqrt(pow(1.93309, 2)+pow(k[0], 2));
				omega0[1] = sqrt(pow(1.93309, 2)+pow(k[1], 2));
				gamma[0] = 3.42222*pow(k[0]+3.663, (long double)-1.56165)+.01;
				gamma[1] = 3.42222*pow(k[1]+3.663, (long double)-1.56165)+.01;
				break;
		}
	}

#ifdef HALF
	ImSelf[0] += -M*Sigma[0]*omega[0]*gamma[0]/(M_PI*(pow(omega[0]-omega0[0], 2)+pow(omega[0]*gamma[0], 2)));
	ImSelf[1] += -M*Sigma[1]*omega[1]*gamma[1]/(M_PI*(pow(omega[1]-omega0[1], 2)+pow(omega[1]*gamma[1], 2)));
	ReSelf[0] += Sigma[0]*(omega[0]-omega0[0])/(M_PI*(pow(omega[0]-omega0[0], 2)+pow(omega[0]*gamma[0], 2)))/2.;
	ReSelf[1] += Sigma[1]*(omega[1]-omega0[1])/(M_PI*(pow(omega[1]-omega0[1], 2)+pow(omega[1]*gamma[1], 2)))/2.;
#else
	ImSelf[0] += -2.*M*Sigma[0]*omega[0]*gamma[0]/(M_PI*(pow(omega[0]-omega0[0], 2)+pow(omega[0]*gamma[0], 2)));
	ImSelf[1] += -2.*M*Sigma[1]*omega[1]*gamma[1]/(M_PI*(pow(omega[1]-omega0[1], 2)+pow(omega[1]*gamma[1], 2)));
	ReSelf[0] += Sigma[0]*(omega[0]-omega0[0])/(M_PI*(pow(omega[0]-omega0[0], 2)+pow(omega[0]*gamma[0], 2)));
	ReSelf[1] += Sigma[1]*(omega[1]-omega0[1])/(M_PI*(pow(omega[1]-omega0[1], 2)+pow(omega[1]*gamma[1], 2)));
#endif

	return;
}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	if(pow(M, 2)+pow(P, 2)+pow(k, 2)+2.*P*k*cos(theta) < 0)
		return(0.);
	else
		return(sqrt(pow(M, 2)+pow(P, 2)+pow(k, 2)+2.*P*k*cos(theta)));
}

long double Set_Temp(int T)
{
	const long double Temps[] = {0, .194, .258, .32, .4, .04, .04};
	return(Temps[T]);
}

long double Fermi(long double omega, int T)	//Fermi factor
{
	static long double Temp = Set_Temp(T);

	if(Temp == 0)
	{
		if(omega >= 0)	//Fermi factor for vacuum
			return(0);
		else
			return(1);
	}
	return(1./(1.+exp(omega/Temp)));
}

long double Potential1(long double Par[], long double k0, long double k)	//Single vertex of potiential without coupling constant
{
#if VERSION == Exp
	return(exp(-abs(-4.*pow(k0, 2)+4.*pow(k, 2))/pow(Par[1], 2)));
#elif VERSION == 22
	return(pow(Par[1], 2.)/(pow(Par[1], 2.)+abs(-4.*pow(k0, 2)+4.*pow(k, 2))));
#elif VERSION == 42
	return(pow(Par[1], 4.)/(pow(Par[1], 4.)+pow(-4.*pow(k0, 2)+4.*pow(k, 2), 2)));
#elif VERSION == 24
	return(pow(2.*pow(Par[1], 2.)/(2.*pow(Par[1], 2.)+abs(-4.*pow(k0, 2)+4.*pow(k, 2))), 2));
#endif
}

long double Potential2(long double Par[], long double k0, long double k)	//Two vertices of potential with coupling constant
{
#if VERSION == Exp
	return(Par[0]*exp(-2.*abs(-4.*pow(k0, 2)+4.*pow(k, 2))/pow(Par[1], 2)));
#elif VERSION == 22
	return(Par[0]*pow(pow(Par[1], 2.)/(pow(Par[1], 2.)+abs(-4.*pow(k0, 2)+4.*pow(k, 2))), 2));
#elif VERSION == 42
	return(Par[0]*pow(pow(Par[1], 4.)/(pow(Par[1], 4.)+pow(-4.*pow(k0, 2)+4.*pow(k, 2), 2)), 2));
#elif VERSION == 24
	return(Par[0]*pow(2.*pow(Par[1], 2.)/(2.*pow(Par[1], 2.)+abs(-4.*pow(k0, 2)+4.*pow(k, 2))), 4));
#endif
}

long double Interacting_Linear_Trace(long double Par[])
{
	return(sqrt(3.*Par[4]/(8.*pow(Par[2], 2))));
}

long double Interacting_Quad_Trace(long double Par[], long double k0, long double k)
{
	return(Par[4]/4.-pow(k0, 2)+pow(k, 2))/(2.*pow(Par[2], 2));
}

long double Imk0_Integrand(long double Par[], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	static long double q[2] = {Energy(0, Par[3]/2., k, theta), Energy(0, Par[3]/2., -k, theta)};
	static long double k_old = k;
	long double omega[2] = {sqrt(Par[4]+pow(Par[3], 2))/2.+k0, sqrt(Par[4]+pow(Par[3], 2))/2.-k0};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];

	if(k_old != k)
	{
		k_old = k;
		q[0] = Energy(0, Par[3]/2., k, theta);
		q[1] = Energy(0, Par[3]/2., -k, theta);
	}

	//Self_Energy(Par[2], omega, q, Temp, ImSelf, ReSelf);
	ImSelf_Energy(Par[2], omega, q, Temp, ImSelf);
	ReSelf_Energy(Par[2], omega, q, Temp, ReSelf);

	return(-((4.*ImSelf[0]*ImSelf[1]*pow(Par[2], 2)*(1.-fermi[0]-fermi[1]))/((pow(pow(omega[0], 2)-pow(q[0], 2)-pow(Par[2], 2)-2.*Par[2]*ReSelf[0], 2)+pow(ImSelf[0], 2))*(pow(pow(omega[1], 2)-pow(q[1], 2)-pow(Par[2], 2)-2.*Par[2]*ReSelf[1], 2)+pow(ImSelf[1], 2)))));
}
