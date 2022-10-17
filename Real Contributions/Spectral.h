//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<cfloat>
#include<complex>
#include<queue>
#include"Elements.h"
using namespace std;

//Integrals that define results
Elements<Around> theta_Int(long double[], int);			//Theta integral
Elements<Around> theta_Int(long double[], int, long double, long double, int);
Elements<Around> k_Int(long double[], int, long double);
Elements<Around> k_Int(long double[], int, long double, long double, long double, int, int);
Around Dispersion(long double[], int, long double, long double, long double);	//Dispersion relation for turning ImG_12 into ReG_12
Around Dispersion(long double[], int, long double, long double, long double, long double, long double, Around, int, int);	//Dispersion relation for turning ImG_12 into ReG_12
Around k0_Int(long double[], int, long double, long double);			//k0 integral aka energy integral
Around k0_Int(long double[], int, long double, long double, long double, long double, int, int);			//k0 integral aka energy integral

Elements<long double> Integrand(long double[], long double, long double, int);
long double ReG12(long double, long double, long double, long double, long double);
long double ImG12(long double, long double, long double, long double, long double);
complex<long double> G12Reverse(long double, long double, long double, long double, long double, int);
long double ReG12Reverse(long double, long double, long double, long double, long double, int);
long double ImG12Reverse(long double, long double, long double, long double, long double, int);
complex<long double> G12BbS(long double, long double, long double, long double, long double, int);
complex<long double> SigmaReverse(long double, long double, long double, long double, long double, int);
complex<long double> SigmaBbS(long double, long double, long double, long double, long double, int);
complex<long double> ReBuiltSigma(long double[], int, long double, long double);

//Functions for finding points of interest in the k integral
void Characterize_k_Int(long double[], int, long double, long double[], long double[], int&);	//Returns the poles of the k integral's integrands
bool Newton_Method_k(long double&, long double, long double, long double, long double, long double, long double(*)(long double, long double), long double(*)(long double, long double, long double, long double, long double));			//Returns the k-intesection of two features of interest: potential peak, on-shell peak, light-like boundaries
long double V_Plus(long double, long double);			//Potiential peaks, matches a potential not in use
long double V_Minus(long double, long double);
long double Emm(long double, long double, long double, long double, long double);		//on-shell peaks
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
void ReSelf_Energy(long double, long double, long double[], int, long double[]);			//Returns the real single quark self-energies for both quarks, contains an alternate T=194 MeV solution
long double ReSelf_Energy(long double, long double, long double, int);				//Returns the real single quark self-energies for both quarks, contains an alternate T=194 MeV solution
void Self_Energy(long double, long double, long double[], int, long double[], long double[]);	//Returns the complex single quark self-energies for both quarks, is a simple Breit-Wigner self-energy and alternate to those above
long double Energy(long double, long double, long double, long double);				//Single quark energy, also used to return total momentum by setting M=0
long double Fermi(long double, int);									//Fermi function
long double Set_Temp(int);										//Decodes 0-4 into numeric temprature for Fermi factor
long double Potential1(long double[], long double, long double);					//One vertex of the potiential
long double Potential2(long double[], long double, long double);					//Two vertices of the potiential
long double Interacting_Linear_Trace(long double[]);							//Linear (linear in sqrt(s)) contribution to the interacting trace
long double Interacting_Quad_Trace(long double[], long double, long double);				//Quadratic contribution to the interacting trace
long double Imk0_Integrand(long double[], long double, long double, long double, int);		//Integrand of the k0 integral for positive energy


//Merge Sort. There are a number of semi-sorted lists that need sorting. This will probably beat quick sort under similar conditions.
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
		Temp[i] = List[i+a];

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
#define GAMMA -0.1193456819588533	//Width of single quark propagator
#endif

Elements<Around> theta_Int(long double Par[], int Temp)
{
	//if(Par[3] == 0)	//Short cut for P=0, theta integral is analytic
	//	return(k_Int(Par, Temp, M_PI/2., true)/pow(2.*M_PI,2)*2.);

	long double x1;
	long double a = 0, b;					//Sub-interval limits of integration
	long double Boundary_theta[] = {1./17., 0.3, 0.08};	//Extra boundary values
	Elements<Around> F[2];					//Sum of ordinate*weights
	Elements<Around> Answer(0,0,0,0);			//Answer to be returned
	int i, j;						//Counters

	if(Par[4] > 0 && Par[3] > sqrt(Par[4]/2.)) //Where the maximum of the theta integral ought to land. It might only be correct for BbS reduction, but is close enough for all other cases. Only valid for s>0 and P>sqrt(s/2)
		x1 = asin(sqrt(Par[4]/2.)/Par[3]);
	else
		x1 = M_PI/10.;	//If it isn't valid, value is needed anyways to split up the integral

	//Don't get too close to the pole or details might get lost
	if(x1>M_PI/10.)
		x1 = M_PI/10.;

	//List of boundaries between subintervals
	long double Range[] = {x1*Boundary_theta[0], x1*Boundary_theta[1], x1, x1*(2.-Boundary_theta[1]), x1*(2.-Boundary_theta[1])*(1.-Boundary_theta[2])+M_PI/2.*Boundary_theta[2], M_PI/2., asin(sqrt(-Par[4])/Par[3]),0,0};

	//Some kind of intersection, probably between the simultanous on-shell and potential peak, don't rightly remember
	Range[7] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[7] = acos((pow(Range[7],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(Range[7]*Par[3]));
	Range[8] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[8] = acos((pow(Range[8],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(-Range[8]*Par[3]));

	//Bad data trap for NaN and negative boundaries. These are only ones that can NaN or return negative numbers
	if(isnan(Range[6]) || Range[6] < 0) Range[6] = M_PI;
	if(isnan(Range[7])) Range[7] = M_PI;
	if(isnan(Range[8])) Range[8] = M_PI;

	//Put in asending order
	mergeSort(Range, 0, 8);

	for(i = 0; i < 8 && a < M_PI/2.; i++)	//Count through pre-determined intervals
	{
		b = Range[i];	//Upper edge

		Answer += theta_Int(Par, Temp, a, b, 0);	//Add the subinterval to total of the integral

		a = b;	//Upper edge becomes lower edge
	}

	return(Answer/pow(2.*M_PI,2)*2.);
}

Elements<Around> theta_Int(long double Par[], int Temp, long double a, long double b, int deep)
{
//9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight
	long double Disp9[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	long double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	long double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights*/
/*1st/2nd order Newton-Coates integration
	long double Disp9[] = {1};	//Displacement from center
	long double w9[] = {0,1.};	//1st order Newton-Coates weights
	long double w16[]= {4./3.,1./3.}; //2nd order Newton-Coates weights*/
#if ORDER == 37	//37th order Gauss-Legendre integration
//23th order Gauss-Legendre/37th order Gauss-Kronrod integration
	long double Disp[] = {0.1252334085114689154724414, 0.2485057483204692762677910, 0.3678314989981801937526915, 0.4813394504781570929359436, 0.5873179542866174472967024, 0.6840598954700558939449291, 0.7699026741943046870368938, 0.8435581241611532447921419, 0.9041172563704748566784659, 0.9505377959431212965490602, 0.9815606342467192506905491, 0.9969339225295954269123502};	//Displacement from center
	long double wl[] = {0, 0.2491470458134027850005624, 0, 0.2334925365383548087608499, 0, 0.2031674267230659217490645, 0, 0.1600783285433462263346525, 0, 0.10693932599531843096025472, 0, 0.04717533638651182719461596, 0};	//23rd order Gauss-Legendre weight
	long double wh[] =  {0.12555689390547433530429613, 0.1245841645361560734373125, 0.12162630352394838324609976, 0.1167120535017568262935807, 0.11002260497764407263590740, 0.10164973227906027771568877, 0.091549468295049210528171940, 0.07992027533360170149339261, 0.067250907050839930304940940, 0.05369701760775625122888916, 0.038915230469299477115089632, 0.02303608403898223259108458, 0.0082577114331683957576939224};	//37th order Gauss-Kronrod weight
#elif ORDER == 97
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	long double Disp[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	long double wl[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	long double wh[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight
#endif
	long double x1, x2;			//Abscissa
	Elements<Around> F[2];			//Sum of ordinate*weights
	Elements<Around> Answer(0,0,0,0);	//Answer to be returned
	Elements<Around> Holder;
	int i, j = 0;				//Counters

	F[0].null();	//Zero out F for a new round
	F[1].null();	//Zero out F for a new round

	for(j = 0; j < 5; j++)	//Count through points away from center
	{
		x1 = (b+a-Disp[j]*(b-a))/2.;
		x2 = (b+a+Disp[j]*(b-a))/2.;

		Holder = k_Int(Par, Temp, x1);
		F[0] += Holder*w9[j+1];
		F[1] += Holder*w16[j+1];

		Holder = k_Int(Par, Temp, x2);

		F[0] += Holder*w9[j+1];
		F[1] += Holder*w16[j+1];
	}
	Holder = k_Int(Par, Temp, (a+b)/2.);
	F[0] += Holder*w9[0];
	F[1] += Holder*w16[0];

	if(abs(F[0]-F[1])*2./abs(F[0]+F[1]) > 1 && abs(b-a) > FLT_EPSILON)
		Answer = theta_Int(Par, Temp, a, (a+b)/2., deep+1) + theta_Int(Par, Temp, (a+b)/2., b, deep+1);
	else
	{
		F[0].null();	//Zero out F for a new round
		F[1].null();	//Zero out F for a new round

#if ORDER == 37	//Count through points away from center
		for(j = 0; j < 12; j++)
#elif ORDER == 97
		for(j = 0; j < 32; j++)
#endif
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			Holder = k_Int(Par, Temp, x1);
			F[0] += Holder*wl[j+1];
			F[1] += Holder*wh[j+1];

			Holder = k_Int(Par, Temp, x2);
			F[0] += Holder*wl[j+1];
			F[1] += Holder*wh[j+1];
		}
		Holder = k_Int(Par, Temp, (a+b)/2.);
		F[0] += Holder*wl[0];
		F[1] += Holder*wh[0];

		Answer = Elements<Around>(Around(F[1][0],abs(F[0][0]-F[1][0])),Around(F[1][1],abs(F[0][1]-F[1][1])),Around(F[1][2],abs(F[0][2]-F[1][2])),Around(F[1][3],abs(F[0][3]-F[1][3])))*(b-a)/2.;	//Add the subinterval to total of the integral
	}

	return(Answer);
}

Elements<long double> Integrand(long double Par[], long double k, long double theta, int Temp)
{
	long double k0 = (Energy(Par[2], Par[3]/2., k, theta)-Energy(Par[2], Par[3]/2., -k, theta))/2.;
//	Elements<long double> Holder = Elements<long double>(Potential1(Par, k0, k), Interacting_Linear_Trace(Par)*Potential1(Par, k0, k), Interacting_Quad_Trace(Par, k0, k)*Potential1(Par, k0, k), Potential2(Par, k0, k))*ReG12Reverse(Par[2], Par[4], Par[3], k, theta, Temp);
//cerr << Par[3] << "," << Par[4] << "," << k << "," << theta << "," << Holder[0] << "," << Holder[1] << "," << Holder[2] << "," << Holder[3] << "," << ReG12Reverse(Par[2], Par[4], Par[3], k, theta, Temp) << endl;
	return(Elements<long double>(Potential1(Par, k0, k), Interacting_Linear_Trace(Par)*Potential1(Par, k0, k), Interacting_Quad_Trace(Par, k0, k)*Potential1(Par, k0, k), Potential2(Par, k0, k))*ReG12Reverse(Par[2], Par[4], Par[3], k, theta, Temp)*pow(k,2)*sin(theta));
}

Elements<Around> k_Int(long double Par[], int Temp, long double theta)
{
	long double a, b;	//Sub-interval limits of integration

	Elements<Around> Answer(0, 0, 0, 0);	//Answer to be returned
	Elements<Around> Partial;		//Answer for sub-interval for determining completeness

	int Poles;		//Number of poles
	long double zero[26];	//The real part of the signular pole
	long double gamma[26];	//The distance to the singular, maybe
	long double Max;
	int i, j, l;		//Counters, would use 'k', but 'k' is occupied by relative 3-momenta in other parts of program
	int Intervals;		//Number of intervals recorded in Stops

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);	//Find the location of the complex poles
	long double Stops[Poles+22];				//List of pre-determined subintervals

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
	Stops[l+10] = .5*sqrt((Par[4]+pow(Par[3], 2))*(Par[4]-pow(2.*Par[2], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)))+5.*GAMMA;

	for(i = 1; i <= 10; i++)
	{
		Stops[i+l+10] = Stops[l+10]-i*GAMMA;
	}

	for(i = 0; i < l+21; i++)	//Removes stops that are NaN or bigger than necessary
	{
		if(isnan(Stops[i]))
			Stops[i] = -1;
		else if(isinf(Stops[i]) || Stops[i] > 100)
			Stops[i] = 100;
	}

	mergeSort(Stops, 0, l+20);	//Sort the list of sub-intervals
	Stops[l+21] = 660;

	i = 0;
	j = 0;
	while(Stops[j] <= 0)	//Skip past negative sub-intervals and form NaN
		j++;
	for(; j < l+21; j++)
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
		a = b;
		if(((i < Intervals && b+100 < Stops[i]) && (i > 0 && b-Stops[i-1] > 100)) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if(((i < Intervals && 50 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 50)) || Stops[Intervals-1] < a-50)
			b += 50;
		else if(((i < Intervals && 10 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 10)) || Stops[Intervals-1] < a-10)
			b += 10;
		else if(((i < Intervals && 3 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 3)) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}
		else
			b += 3;

		if(b-a < 1)	//use a higher order when the interval is large
			Partial = k_Int(Par, Temp, theta, a, b, 37, 0);
		else
			Partial = k_Int(Par, Temp, theta, a, b, 97, 0);

		Answer += Partial;	//Add the subinterval to total of the integral
	}while(!(Partial[0] == 0) && (i < Intervals || abs(Partial/Answer)/(b-a) >= .0001) && (a <= Max || a <= 20.*sqrt(Par[4]+pow(Par[3], 2)))); //Keep going so long as the last subinterval isn't zero and the intervals haven't been exhausted and the last partial answer for all functions isn't too big compared to the total answer and the highest sub-interval is less than 20E. k bigger than 20E is getting pretty stupid, should be sneaking up on 10^-5 of the answer left

	return(Answer);
}

Elements<Around> k_Int(long double Par[], int Temp, long double theta, long double a, long double b, int order, int deep)
{
	long double Disp16[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	long double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	long double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
//23th order Gauss-Legendre/37th order Gauss-Kronrod integration
	long double Disp37[] = {0.1252334085114689154724414, 0.2485057483204692762677910, 0.3678314989981801937526915, 0.4813394504781570929359436, 0.5873179542866174472967024, 0.6840598954700558939449291, 0.7699026741943046870368938, 0.8435581241611532447921419, 0.9041172563704748566784659, 0.9505377959431212965490602, 0.9815606342467192506905491, 0.9969339225295954269123502};	//Displacement from center
	long double w23[] = {0, 0.2491470458134027850005624, 0, 0.2334925365383548087608499, 0, 0.2031674267230659217490645, 0, 0.1600783285433462263346525, 0, 0.10693932599531843096025472, 0, 0.04717533638651182719461596, 0};	//23rd order Gauss-Legendre weight
	long double w37[] =  {0.12555689390547433530429613, 0.1245841645361560734373125, 0.12162630352394838324609976, 0.1167120535017568262935807, 0.11002260497764407263590740, 0.10164973227906027771568877, 0.091549468295049210528171940, 0.07992027533360170149339261, 0.067250907050839930304940940, 0.05369701760775625122888916, 0.038915230469299477115089632, 0.02303608403898223259108458, 0.0082577114331683957576939224};	//37th order Gauss-Kronrod weight
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	long double Disp97[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	long double w63[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	long double w97[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight
	long double x1, x2;	//Abscissa
	long double k01, k02;	//On-shell relative energy at the abscissa

	Elements<long double> F[2] = {Elements<long double>(0, 0, 0, 0),Elements<long double>(0, 0, 0, 0)};	//Sum of ordinates*weights
	Elements<Around> Answer(0, 0, 0, 0);	//Answer to be returned
	Elements<long double> Holder;

	switch(order)
	{
	case 16:
		for(int l = 0; l < 5; l++)//for(int l = 0; l < 12; l+=2)// //Count through points away from center
		{
			x1 = (b+a-Disp16[l]*(b-a))/2.;
			x2 = (b+a+Disp16[l]*(b-a))/2.;

			Holder = Integrand(Par, x1, theta, Temp);
			F[0] += Holder*w9[l+1];
			F[1] += Holder*w16[l+1];

			Holder = Integrand(Par, x2, theta, Temp);
			F[0] += Holder*w9[l+1];
			F[1] += Holder*w16[l+1];
		}
		x1 = (a+b)/2.;
		Holder = Integrand(Par, x1, theta, Temp);
		F[0] += Holder*pow((a+b)/2., 2)*w9[0];
		F[1] += Holder*pow((a+b)/2., 2)*w16[0];
		break;
	case 37:
		for(int l = 0; l < 12; l++)//for(int l = 0; l < 12; l+=2)// //Count through points away from center
		{
			x1 = (b+a-Disp37[l]*(b-a))/2.;
			x2 = (b+a+Disp37[l]*(b-a))/2.;

			Holder = Integrand(Par, x1, theta, Temp);
			F[0] += Holder*w23[l+1];
			F[1] += Holder*w37[l+1];
			Holder = Integrand(Par, x2, theta, Temp);
			F[0] += Holder*w23[l+1];
			F[1] += Holder*w37[l+1];
		}
		x1 = (a+b)/2.;
		Holder = Integrand(Par, x1, theta, Temp);
		F[0] += Holder*w23[0];
		F[1] += Holder*w37[0];
		break;
	case 97:
		for(int l = 0; l < 32; l++)//for(int l = 0; l < 32; l+=2)// //Count through points away from center
		{
			x1 = (b+a-Disp97[l]*(b-a))/2.;
			x2 = (b+a+Disp97[l]*(b-a))/2.;

			Holder = Integrand(Par, x1, theta, Temp);
			F[0] += Holder*w63[l+1];
			F[1] += Holder*w97[l+1];
			Holder = Integrand(Par, x2, theta, Temp);
			F[0] += Holder*w63[l+1];
			F[1] += Holder*w97[l+1];
		}
		x1 = (a+b)/2.;
		Holder = Integrand(Par, x1, theta, Temp);
		F[0] += Holder*w63[0];
		F[1] += Holder*w97[0];
		break;
	}

	Answer = Elements<Around>(Around(F[1][0], abs(F[0][0]-F[1][0])), Around(F[1][1], abs(F[0][1]-F[1][1])), Around(F[1][2], abs(F[0][2]-F[1][2])), Around(F[1][3], abs(F[0][3]-F[1][3])))*(b-a)/2.;//F[0]*(b-a)/2.;//	//Record the subinterval to total of the integral
	if((Answer[0].RelErr() > 1e-9 || Answer[0].RelErr() > 1e-9 || Answer[0].RelErr() > 1e-9 || Answer[0].RelErr() > 1e-9) && deep < 10 && abs(b/a-(long double)(1.)) > FLT_EPSILON)
		Answer = k_Int(Par, Temp, theta, a, (a+b)/2., order, deep+1) + k_Int(Par, Temp, theta, (a+b)/2., b, order, deep+1);//*/

	return(Answer);
}

complex<long double> ReBuiltSigma(long double Par[], int Temp, long double k, long double theta)
{
	complex<long double> G12 = complex<long double>(Dispersion(Par, Temp, (Energy(1.9, Par[3]/2., k, 0)-Energy(1.9, Par[3]/2., -k, 0))/2., k, theta).Value(),k0_Int(Par, Temp, k, theta).Value());
	long double omega[2] = {Energy(Par[2], Par[3]/2., k, theta), Energy(Par[2], Par[3]/2., -k, theta)};
	complex<long double> Sigma = (-G12*(pow(omega[1],2)*omega[0]-omega[1]*pow(omega[0],2))+sqrt(-(long double)(2.)*G12*pow(Par[2]*omega[1],2)*omega[0]-(long double)(2.)*G12*omega[1]*pow(Par[2]*omega[0],2)+pow(G12*Par[3]*omega[1]*omega[0],2)+Par[4]*pow(G12*omega[1]*omega[0],2)))/(G12*omega[1]*omega[0]);

	if(Sigma.imag() > 0)
		Sigma = (-G12*(pow(omega[1],2)*omega[0]-omega[1]*pow(omega[0],2))-sqrt(-(long double)(2.)*G12*pow(Par[2]*omega[1],2)*omega[0]-(long double)(2.)*G12*omega[1]*pow(Par[2]*omega[0],2)+pow(G12*Par[3]*omega[1]*omega[0],2)+Par[4]*pow(G12*omega[1]*omega[0],2)))/(G12*omega[1]*omega[0]);

	return(Sigma);
}

complex<long double> SigmaBbS(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {Energy(M, P/2., k, theta), Energy(M, P/2., -k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp);
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp);
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp);
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp);

	return(complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]));
}

complex<long double> SigmaReverse(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {sqrt(s+pow(P,2))-Energy(M, P/2., -k, theta), sqrt(s+pow(P,2))-Energy(M, P/2., k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp)/5.;
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp)/5.;
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp)/2.;
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp)/2.;

	return(complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]));
}

complex<long double> G12BbS(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {Energy(M, P/2., k, theta), Energy(M, P/2., -k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];
	long double Vacuum_Width = 0;

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp);
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp);
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp);
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp);

	if(s >= 0)
		Vacuum_Width = GAMMA*(25.*(sqrt(s)+2.*M*tanh((80.*pow(M,2))/73.)+(2.*M-sqrt(s))*tanh((20.*(-4.*pow(M,2)+s))/73.)))/(76.+50.*M*tanh((80.*pow(M,2))/73.)+(-76.+50.*M)*tanh(23104./9125.-(80.*pow(M,2))/73.));

	return(2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)+complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]),2)+complex<long double>(0,Vacuum_Width))));
}

complex<long double> G12Reverse(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {sqrt(s+pow(P,2))-Energy(M, P/2., -k, theta), sqrt(s+pow(P,2))-Energy(M, P/2., k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];
	long double Vacuum_Width = 0;

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp)/5.;
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp)/5.;
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp)/2.;
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp)/2.;

	if(s >= 0)
		Vacuum_Width = -GAMMA*(25.*(sqrt(s)+2.*M*tanh((80.*pow(M,2))/73.)+(2.*M-sqrt(s))*tanh((20.*(-4.*pow(M,2)+s))/73.)))/(76.+50.*M*tanh((80.*pow(M,2))/73.)+(-76.+50.*M)*tanh(23104./9125.-(80.*pow(M,2))/73.));

	return(2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)+complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]),2)+complex<long double>(0,Vacuum_Width))));
}

long double ReG12Reverse(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {sqrt(s+pow(P,2))-Energy(M, P/2., -k, theta), sqrt(s+pow(P,2))-Energy(M, P/2., k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];
	long double Vacuum_Width = 0;

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp)/5.;
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp)/5.;
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp)/2.;
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp)/2.;

	if(s >= 0)
		Vacuum_Width = -GAMMA*(25.*(sqrt(s)+2.*M*tanh((80.*pow(M,2))/73.)+(2.*M-sqrt(s))*tanh((20.*(-4.*pow(M,2)+s))/73.)))/(76.+50.*M*tanh((80.*pow(M,2))/73.)+(-76.+50.*M)*tanh(23104./9125.-(80.*pow(M,2))/73.));

	return(2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)+complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]),2)+complex<long double>(0,Vacuum_Width)))).real();
}

long double ImG12Reverse(long double M, long double s, long double P, long double k, long double theta, int Temp)
{
	long double q[2] = {Energy(0, P/2., k, theta), Energy(0, P/2., -k, theta)};
	long double omega[2] = {sqrt(s+pow(P,2))-Energy(M, P/2., -k, theta), sqrt(s+pow(P,2))-Energy(M, P/2., k, theta)};
	long double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];
	long double Vacuum_Width = 0;

	ImSelf[0] = ImSelf_Energy(M, omega[0], q[0], Temp)/5.;
	ImSelf[1] = ImSelf_Energy(M, omega[1], q[1], Temp)/5.;
	ReSelf[0] = ReSelf_Energy(M, omega[0], q[0], Temp)/2.;
	ReSelf[1] = ReSelf_Energy(M, omega[1], q[1], Temp)/2.;

	if(s >= 0)
		Vacuum_Width = -GAMMA*(25.*(sqrt(s)+2.*M*tanh((80.*pow(M,2))/73.)+(2.*M-sqrt(s))*tanh((20.*(-4.*pow(M,2)+s))/73.)))/(76.+50.*M*tanh((80.*pow(M,2))/73.)+(-76.+50.*M)*tanh(23104./9125.-(80.*pow(M,2))/73.));

	return(2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)+complex<long double>(ReSelf[0],ImSelf[0])+complex<long double>(ReSelf[1],ImSelf[1]),2)+complex<long double>(0,Vacuum_Width)))).imag();
}

long double ImG12(long double M, long double s, long double P, long double k, long double theta)
{
	return((-2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)))*.14/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(pow(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2),2)+pow(.14,2))));
}

long double ReG12(long double M, long double s, long double P, long double k, long double theta)
{
	return((2.*pow(M,2)*(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta)))*(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2))/(Energy(M,P/2.,k,theta)*Energy(M,P/2.,-k,theta)*(pow(s+pow(P,2)-pow(Energy(M,P/2.,k,theta)+Energy(M,P/2.,-k,theta),2),2)+pow(.14,2))));
}

Around Dispersion(long double Par[], int Temp, long double k0, long double k, long double theta)
{
	long double a, b;	//Sub-interval limits of integration
	long double Min;	//Lower limit of integration
	long double Max = 0;	//Upper limit of principal value integration
	Around ImG12 = k0_Int(Par, Temp, k, theta);		//Holder of the ImG12 that belongs to the other half of G12 that is calculated here

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	long double Boundary_k_k0[] = {0.0491, .421};
	long double Range[] = {-Boundary_k_k0[1], -Boundary_k_k0[0], 0, Boundary_k_k0[0], Boundary_k_k0[1]};	//Number of gamma from center
	long double ParLoc[5] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//Local copy of the parameters as Par[4] corrisponds to s and ParLoc[4] is s'

	Around Answer(0, 0);	//Results to be returned
	Around Partial;	//Partial results to examine convergance
	Around Holder;

	long double zero[5];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[5];	//Imaginary part of poles
	int Poles;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities
	int Principal = 1;	//Principal value integral calculating or not

	Characterize_Dispersion(ParLoc, Temp, k0, k, theta, zero, gamma, Poles);
	long double Stops[Poles*5+9];		//Extra stops to ensure correctness

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(!isnan(gamma[i]))	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 5; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Extra Regions required by poles
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
	Stops[l+6] = -pow(Par[3], 2);
	Stops[l+7] = 4.*pow(k, 2);

	mergeSort(Stops, 0, l+7);
	Stops[l+8] = Stops[l+7]+100;	//Adds the minimum end point to keep the integration going

	Min = a = b = -pow(Par[3], 2);	//Start from s'=-P^2
	if(abs(Min-Par[4])<.1 && Min-.1 > -pow(Par[3], 2))
	{
		Min -= .1;
		a -= .1;
		b -= .1;
	}

	for(i = 0; Stops[i]+3 < Par[4]; i++)	//Replace all stops too close to s with s. Will skip over them in the s' traversal.
		if(abs(b/Par[4]-(long double)(1.)) < FLT_EPSILON)
			Stops[i] = Par[4];

	for(i = 0; Stops[i] < Min+3 || Holder == 0; i++)	//Remove any stops below the minimum of the limit of integration. Faster to illiminate here than by popping
	{
		ParLoc[4] = Stops[i];
		Holder = k0_Int(ParLoc, Temp, k, theta);
		if(!(Holder == 0))
		{
			if(i!=0) i--;
			Min = a = b = Stops[i];
			break;
		}
	}

	Intervals = l+9;

	do
	{
		if(((i < Intervals && b+100 < Stops[i]) && (i > 0 && b-Stops[i-1] > 100)) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if(((i < Intervals && 50 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 50)) || Stops[Intervals-1] < a-50)
			b += 50;
		else if(((i < Intervals && 10 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 10)) || Stops[Intervals-1] < a-10)
			b += 10;
		else if(((i < Intervals && 3 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 3)) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			while(i < Intervals && (b == Min || abs(b/a-(long double)(1.)) < FLT_EPSILON) || b-a == 0)
			{
				b = Stops[i];
				i++;
			}

			if(b-a > 3)
			{
				b = a+3;
				i--;
			}
		}

		if(a > Par[4]+100. && Max == 0)
		{
			Principal = 0;
			Max = a;
		}

		/*if(chrono::duration_cast<chrono::seconds>(chrono::system_clock::now()-Start_Time).count() > Allotment)
		{
			throw(chrono::duration_cast<chrono::seconds>(chrono::system_clock::now()-Start_Time).count());
		}*/

		if(abs(Par[4]-a)<.25 || abs(b-Par[4])<.25)
			Partial = Dispersion(Par, Temp, k0, k, theta, a, b, ImG12*Principal, 97, 0);
		else
			Partial = Dispersion(Par, Temp, k0, k, theta, a, b, ImG12*Principal, 16, 0);
		Answer += Partial;		//Add the Region to the total
		a = b;
	}while((a < 4.*(pow(k, 2)+pow(Par[2], 2))+50. && i < Intervals) || Partial/Answer > 1e-6 || isinf(log(abs((a-Par[4])))));	//Keep going while intervals aren't exhausted and upper limit of integration not excceeded or until convergance

	if(Max == 0)	//Just in case it terminates before getting to s+100
		Max = a;
	if(abs(ImG12) >= 1e-12)
		return((Answer+ImG12*log(abs((Max-Par[4])/(Par[4]-Min))))/M_PI);
	return(Answer/M_PI);
}

Around Dispersion(long double Par[], int Temp, long double k0, long double k, long double theta, long double a, long double b, Around ImG12, int order, int deep)
{
//9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight
	long double Disp9[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	long double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	long double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	long double Disp97[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	long double w63[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	long double w97[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight

	Around F[2] = {0, 0};	//Sum of ordinates*weights
	Around Answer(0, 0);	//Results to be returned
	Around Holder;

	long double ParLoc[5] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//Local copy of the parameters as Par[4] corrisponds to s and ParLoc[4] is s'

	switch(order)
	{
	case 16:
		for(int l = 0; l < 5; l++) //Count through points away from center
		{
			ParLoc[4] = (b+a-Disp9[l]*(b-a))/2.;
			Holder = k0_Int(ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[l+1];
			ParLoc[4] = (b+a+Disp9[l]*(b-a))/2.;
			Holder = k0_Int(ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[l+1];
		}
		ParLoc[4] = (b+a)/2.;
		Holder = k0_Int(ParLoc, Temp, k, theta);
		F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[0];
		F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[0];
		break;
	case 97:
		for(int l = 0; l < 32; l++) //Count through points away from center
		{
			ParLoc[4] = (b+a-Disp97[l]*(b-a))/2.;
			Holder = k0_Int(ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[l+1];
			ParLoc[4] = (b+a+Disp97[l]*(b-a))/2.;
			Holder = k0_Int(ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[l+1];
		}
		ParLoc[4] = (b+a)/2.;
		Holder = k0_Int(ParLoc, Temp, k, theta);
		F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[0];
		F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[0];
		break;
	}

	Answer = Around(F[1], abs(F[0]-F[1]))*(b-a)/2.;
	if(!(Answer.RelErr() < 1e-6 || Answer.Error() < 1e-7) && deep < 4 && abs(b/a-(long double)(1.)) > FLT_EPSILON)
		Answer = Dispersion(Par, Temp, k0, k, theta, a, (a+b)/2., ImG12, order, deep+1) + Dispersion(Par, Temp, k0, k, theta, (a+b)/2., b, ImG12, order, deep+1);

	return(Answer);
}

Around k0_Int(long double Par[], int Temp, long double k, long double theta)
{
	if(Par[4]+pow(Par[3], 2) <= 0) //Can't go below zero energy for the dispersion relation.
		return(Around(0.));

	long double a, b;	//Sub-interval limits of integration
	long double Max;	//Upper limit of integration

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	long double Boundary_k_k0[] = {.421, 4.85};
	long double Range[] = {-Boundary_k_k0[1], -Boundary_k_k0[0], 0, Boundary_k_k0[0], Boundary_k_k0[1]};	//Number of gamma from center

	Around Answer(0, 0);	//Results to be returned
	Around Partial;		//Partial sum to determine continuation

	long double zero[12];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[12];	//Imaginary part of poles
	long double Caution = (Par[4]-pow(2.*Par[2], 2)>0)?(Par[3]*cos(theta)/2.*sqrt((Par[4]-pow(2.*Par[2], 2))/(Par[4]-pow(Par[3]*sin(theta), 2)))):((Energy(Par[2], Par[3]/2., k, theta)-Energy(0, Par[3]/2., -k, theta))/2.);	//Double on-shell in the k0 direction
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities

	Characterize_k0_Int(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	long double Stops[Poles*5+6];					//Intervals that are required by integrating near poles

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(!isnan(gamma[i]) && !isnan(zero[i]))	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 0; j < 5; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Extra Regions required by poles
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
	Stops[l+2] = Energy(0, Par[3]/2., k, theta)+sqrt(Par[4]+pow(Par[3], 2))/2.;	//Pretty sure this is the negative energy solution of the lower light-like edge
	Stops[l+3] = sqrt(Par[4]+pow(Par[3], 2))/2.+Energy(0, Par[3]/2., -k, theta);	//Pretty sure this is the negative energy solution of the upper light-like edge
	Stops[l+4] = sqrt(Par[4]+pow(Par[3], 2))/2.;					//Upper energy boundary (E/2)
	Stops[l+5] = -sqrt(Par[4]+pow(Par[3], 2))/2.;					//Lower energy boundary (-E/2)

	if(Temp != 0)
	{
		a = b = -sqrt(Par[4]+pow(Par[3], 2))/2.;	//Lower edge for non-vacuum
		Max = sqrt(Par[4]+pow(Par[3], 2))/2.;		//Upper edge for non-vacuum
	}
	else
	{
		a = b = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+pow(Par[3], 2))/2.;	//Lower edge for vacuum
		Max = sqrt(Par[4]+pow(Par[3], 2))/2.-Energy(0, Par[3]/2., -k, theta);	//Upper edge for vacuum
		if(a>Max)
		{
			a = b = sqrt(Par[4]+pow(Par[3], 2))/2.-Energy(0, Par[3]/2., -k, theta);	//Lower edge for vacuum
			Max = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+pow(Par[3], 2))/2.;	//Upper edge for vacuum
		}
	}

	for(i = 0; i < l+6; i++)
	{
		if(Stops[i] < a)
			Stops[i] = -Stops[i];
	}

	mergeSort(Stops, 0, l+5);	//Sort the Regions

	i = 0;
	j = 0;
	while(Stops[j] == Stops[j]+1. || isnan(Stops[j]))	//Remove Regions that duplicates or below the lower edge
		j++;
	for(; j < l+6; j++)
	{
		if(((i > 0 && Stops[i-1] != Stops[j]) || i == 0))	//Remove dublicates and intervals above the upper edge
		{
			Stops[i] = Stops[j];
			i++;
		}
	}
	Intervals = i;	//Record the number of intervals

	i = 1;	//The first point should be the lower limit of integration. That's where we start. Next Region is what we need to be looking for
	do
	{
		if(((i < Intervals && b+100 < Stops[i]) && (i > 0 && b-Stops[i-1] > 100)) || Stops[Intervals-1] < a-100)	//Middle of nowhere intervals not specified by Stops
			b += 100;
		else if(((i < Intervals && 50 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 50)) || Stops[Intervals-1] < a-50)
			b += 50;
		else if(((i < Intervals && 10 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 10)) || Stops[Intervals-1] < a-10)
			b += 10;
		else if(((i < Intervals && 3 < Stops[i]-b) && (i > 0 && b-Stops[i-1] > 3)) || Stops[Intervals-1] < a-3)
			b += 3;
		else if(i < Intervals)
		{
			b = Stops[i];
			i++;
		}

		if(b > Max && a < Max)
			b = Max;	//Be sure E/2 is and sub-interval boundary

		if((abs(a-Caution) < 1 || abs(b-Caution) < 1) && Par[4] > pow(2.*Par[2], 2) && abs(k-.5*sqrt((Par[4]-pow(2.*Par[2], 2))*(Par[4]+pow(Par[3], 2))/(Par[4]+pow(Par[3]*sin(theta), 2)))) < 1)
			Partial = k0_Int(Par, Temp, k, theta, a, b, 97, 0);
		else if((abs(a-Caution) < 1 || abs(b-Caution) < 1))
			Partial = k0_Int(Par, Temp, k, theta, a, b, 97, 0);
		else
			Partial = k0_Int(Par, Temp, k, theta, a, b, 97, 0);

		Answer += Partial;		//Add the Region to the total
		a = b;
	}while((i < Intervals || abs(Partial/Answer) >= .0001) && a < Max);	//Keep going while intervals aren't exhausted and upper limit of ntegration not excceeded

	return(Answer/M_PI);
}

Around k0_Int(long double Par[], int Temp, long double k, long double theta, long double a, long double b, int order, int deep)
{
//9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight
	long double Disp9[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	long double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	long double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
//23th order Gauss-Legendre/37th order Gauss-Kronrod integration
	long double Disp37[] = {0.1252334085114689154724414, 0.2485057483204692762677910, 0.3678314989981801937526915, 0.4813394504781570929359436, 0.5873179542866174472967024, 0.6840598954700558939449291, 0.7699026741943046870368938, 0.8435581241611532447921419, 0.9041172563704748566784659, 0.9505377959431212965490602, 0.9815606342467192506905491, 0.9969339225295954269123502};	//Displacement from center
	long double w23[] = {0, 0.2491470458134027850005624, 0, 0.2334925365383548087608499, 0, 0.2031674267230659217490645, 0, 0.1600783285433462263346525, 0, 0.10693932599531843096025472, 0, 0.04717533638651182719461596, 0};	//23rd order Gauss-Legendre weight
	long double w37[] =  {0.12555689390547433530429613, 0.1245841645361560734373125, 0.12162630352394838324609976, 0.1167120535017568262935807, 0.11002260497764407263590740, 0.10164973227906027771568877, 0.091549468295049210528171940, 0.07992027533360170149339261, 0.067250907050839930304940940, 0.05369701760775625122888916, 0.038915230469299477115089632, 0.02303608403898223259108458, 0.0082577114331683957576939224};	//37th order Gauss-Kronrod weight
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	long double Disp97[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	long double w63[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	long double w97[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight
	long double x1, x2;	//Abscissa

	long double F[2] = {0, 0};
	Around Answer;

	switch(order)
	{
	case 97:
		for(int j = 0; j < 32; j++)//for(int j = 0; j < 32; j+=2)//
		{
			x1 = (b+a-Disp97[j]*(b-a))/2.;
			x2 = (b+a+Disp97[j]*(b-a))/2.;

			F[0] += Imk0_Integrand(Par, x1, k, theta, Temp)*w63[j+1];
			F[1] += Imk0_Integrand(Par, x1, k, theta, Temp)*w97[j+1];
			F[0] += Imk0_Integrand(Par, x2, k, theta, Temp)*w63[j+1];
			F[1] += Imk0_Integrand(Par, x2, k, theta, Temp)*w97[j+1];
		}
		F[0] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w63[0];
		F[1] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w97[0];
		break;
	case 37:
		for(int j = 0; j < 12; j++)//for(int j = 0; j < 12; j+=2)//
		{
			x1 = (b+a-Disp37[j]*(b-a))/2.;
			x2 = (b+a+Disp37[j]*(b-a))/2.;

			F[0] += Imk0_Integrand(Par, x1, k, theta, Temp)*w23[j+1];
			F[1] += Imk0_Integrand(Par, x1, k, theta, Temp)*w37[j+1];
			F[0] += Imk0_Integrand(Par, x2, k, theta, Temp)*w23[j+1];
			F[1] += Imk0_Integrand(Par, x2, k, theta, Temp)*w37[j+1];
		}
		F[0] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w23[0];
		F[1] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w37[0];
		break;
	case 16:
		for(int j = 0; j < 5; j++)//for(int j = 1; j < 5; j+=2)//
		{
			x1 = (b+a-Disp9[j]*(b-a))/2.;
			x2 = (b+a+Disp9[j]*(b-a))/2.;

			F[0] += Imk0_Integrand(Par, x1, k, theta, Temp)*w9[j+1];
			F[1] += Imk0_Integrand(Par, x1, k, theta, Temp)*w16[j+1];
			F[0] += Imk0_Integrand(Par, x2, k, theta, Temp)*w9[j+1];
			F[1] += Imk0_Integrand(Par, x2, k, theta, Temp)*w16[j+1];
		}
		F[0] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w9[0];
		F[1] += Imk0_Integrand(Par, (a+b)/2., k, theta, Temp)*w16[0];
		break;
	}

	Answer = Around(F[1], abs(F[0]-F[1]))*(b-a)/2.;//Around(F[0])*(b-a)/2.;//
	/*if(Answer.RelErr() > 1e-8 && deep < 4 && abs(b/a-(long double)(1.)) > FLT_EPSILON)
		Answer = k0_Int(Par, Temp, k, theta, a, (a+b)/2., order, deep+1) + k0_Int(Par, Temp, k, theta, (a+b)/2., b, order, deep+1);//*/

	return(Answer);
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

	if(true)//Temp != 0)	//media estimate
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
	long double Danger[] = {.5*(sqrt(Par[4]+pow(Par[3], 2))+real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)-k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4)))))), .5*(-sqrt(Par[4]+pow(Par[3], 2))-real(sqrt(complex<long double>(4.*(pow(k, 2)+pow(Par[2], 2)+k*Par[3]*cos(theta))+pow(Par[3], 2)-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))))};

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
	zero[4] = .5*(pow(2.*k, 2)+pow(2.*Par[2], 2)-pow(Par[3], 2)+sqrt(pow(pow(Par[3], 2)-pow(2.*k, 2)-pow(2.*Par[2], 2), 2)+pow(2.*Par[3], 2)*(pow(2.*Par[2], 2)+pow(2.*k*sin(theta), 2))));	//Intersection of two on-shells in terms of k, M, P, and theta

	//Calcluate and record the widths of the peaks
	Par[4] = zero[0];
	gamma[0] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[1];
	gamma[1] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[2];
	gamma[2] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[3];
	gamma[3] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));
	Par[4] = zero[4];
	gamma[4] = abs(sp_Width(Par, k0, k, theta, Temp, Imk0_Integrand));

	Poles = 5;
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
	static long double k_old[2] = {-1,-1}; //Previous value of k to know if the parmeters need to recalculated

	Results[0] = 0;
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

long double ReSelf_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Sigma;	//Strength
	long double x0, x1;	//Centrality markers
	long double gamma;	//Width
	long double Shift, M_T;
	long double Results;

	if(Temp == 0 || Temp == 5)
		return(0);

	switch(Temp)
	{
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .257498/sqrt(pow(k, 2)+pow(1.33201, 2))+.00762638;
			x0 = sqrt(pow(k, 2)+pow(1.54778+Shift, 2))+.276509;
			x1 = sqrt(pow(k, 2)+pow(1.49799+Shift, 2))+.246719;
			gamma = .658734/sqrt(pow(k, 2)+pow(3.35217, 2))+.0815109;
			break;*/
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .212571/sqrt(pow(k, 2)+pow(1.17821, 2))+.00762638;
			x0 = sqrt(pow(k, 2)+pow(1.57536+Shift, 2))+.259147;
			x1 = sqrt(pow(k, 2)+pow(1.50194+Shift, 2))+.222526;
			gamma = .336699/sqrt(pow(k, 2)+pow(1.87956, 2))+.0651449;
			break;
		case 2://258MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .307972/sqrt(pow(k, 2)+pow(1.41483, 2))+.0101423;
			x0 = sqrt(pow(k, 2)+pow(1.56476+Shift, 2))+.251031;
			x1 = sqrt(pow(k, 2)+pow(1.50194+Shift, 2))+.222526;
			gamma = .550628/sqrt(pow(k, 2)+pow(2.43968, 2))+.0981269;
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .339131/sqrt(pow(k, 2)+pow(1.43308, 2))+.0125796;
			x0 = sqrt(pow(k, 2)+pow(1.55034+Shift, 2))+.257788;
			x1 = sqrt(pow(k, 2)+pow(1.46999+Shift, 2))+.231821;
			gamma = .615278/sqrt(pow(k, 2)+pow(2.22298, 2))+.143376;
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .304841/sqrt(pow(k, 2)+pow(1.42911, 2))+.0157245;
			x0 = sqrt(pow(k, 2)+pow(1.55511+Shift, 2))+.231105;
			x1 = sqrt(pow(k, 2)+pow(1.44714+Shift, 2))+.20956;
			gamma = .862629/sqrt(pow(k, 2)+pow(2.67193, 2))+.189598;
			break;
		default:
			Sigma = .188045;
			x0 = 1.83451;
			x1 = 1.72447;
			gamma = .244282;
	}

#ifdef HALF
	Results = Sigma*(omega-x0)/(pow(omega-x1, 2)+gamma)/2.;
#else
	Results = Sigma*(omega-x0)/(pow(omega-x1, 2)+gamma);
#endif
	return(Results);
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
