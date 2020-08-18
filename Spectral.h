//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
#include<cfloat>
#include"Elements.h"
using namespace std;

//Integrals that define results and ancillary functions
Elements theta_Int(long double[], int);	//Integrates the theta results
Elements k_Int(long double[], int, long double);	//Integrates the k momentum results
Elements Folding(long double[], int, long double, long double);	//Folding integral, energy integral
void Characterize_k_Int(long double[], int, long double, long double[], long double[], int&);	//Returns the poles of the k integral's integrands
bool Newton_Method_k(long double&, long double, long double, long double, long double, long double, long double(*)(long double, long double, long double, long double), long double(*)(long double, long double, long double, long double, long double));	//Returns the k-intesection of a potiential and on-shell peak
long double V_Plus(long double, long double, long double, long double);	//Potiential peaks
long double V_Minus(long double, long double, long double, long double);
long double Emm(long double, long double, long double, long double, long double);//on-shell peaks
long double Epm(long double, long double, long double, long double, long double);
long double mEmp(long double, long double, long double, long double, long double);
long double Emp(long double, long double, long double, long double, long double);
long double mEpp(long double, long double, long double, long double, long double);
long double Epp(long double, long double, long double, long double, long double);
long double Upper_Bound(long double, long double, long double, long double, long double);
long double Lower_Bound(long double, long double, long double, long double, long double);
void Characterize_Folding(long double[], int, long double, long double, long double[], long double[], int&);	//Returns the poles of the folding integral's integrands
long double Newton_Method_k0(long double, long double[], long double, long double, int, long double (*)(long double[], long double, long double, long double, int));
long double omega_Width(long double, long double[], long double, long double, int, long double (*)(long double[], long double, long double, long double, int));

//Straight Functions everything is built from
long double ImSelf_Energy(long double, long double, long double, int, long double, long double[]);	//Imaginary single quark self energy
long double ReSelf_Energy(long double, long double, long double, int, long double, long double[]);	//Real single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Potential_on(long double[]);	//On-shell potential for the on-shell T-Matrix
long double Potential1(long double[], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double WidthFraction(long double, long double[]);	//Energy dependance of the self-energy width to go to vacuum as quark energy goes high
long double MassFraction(long double, long double[]);	//Energy dependance of the quark mass to go to vacuum as quark energy goes high
long double ImD(long double, long double, long double, int, long double[]);	//Single quark spectral function
long double Spin_Sum1(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Linear(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Quad(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Sum2(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double ImFolding_Integrand1(long double[], long double, long double, long double, int);	//Integrand of the folding integral for positive energy
long double ImFolding_Integrand2(long double[], long double, long double, long double, int);	//Integrand of the folding integral for negative energy (anti-particle/particle-hole)
long double ReFolding_Integrand1(long double[], long double, long double, long double, int);	//Integrand of the folding integral for positive energy
long double ReFolding_Integrand2(long double[], long double, long double, long double, int);	//Integrand of the folding integral for negative energy (anti-particle/particle-hole)

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

#ifndef GAMMA
#define GAMMA -.015
#endif

#define LATTICE_N 24.
#define A_INVERSE 2.8
long double Boundary[] = {0.00865, 0.0267, 0.0491, 0.0985, .421, .802, 1.01, 4.85, 1.5, 2.5, 3, 4, 5.5, 7.7, 1./17., 0.3, 0.08};

//long double Par[] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[], int Temp)
{
	//long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	//long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
	long double x1, x2;	//Abscissa
	if(Par[4] > 0 && Par[3] > sqrt(Par[4]/2.)) //The maximum of the theta integral, valid for all s>0, arcsin(sqrt(s/(2P^2))) or pi/2
		x1 = asin(sqrt(Par[4]/2.)/Par[3]);
	else
		x1 = M_PI/10.;
	if(x1>M_PI/10.)
		x1 = M_PI/10.;
	long double Range[] = {x1*Boundary[14], x1*Boundary[15], x1, x1*(2.-Boundary[15]), x1*(2.-Boundary[15])*(1.-Boundary[16])+M_PI/2.*Boundary[16], M_PI/2., asin(sqrt(-Par[4])/Par[3]),0,0};
	Elements F;	//Sum of ordinate*weights
	Elements Answer = Elements(0,0,0,0,0);	//Answer to be returned
	Elements holder;
	long double a = 0, b;	//Sub-interval limits of integration
	int i, j;	//Counters
	//ofstream Table("theta Table", ios::app);

	if(Par[3] == 0)	//Short cut for P=0, theta integral is analytic
		return(k_Int(Par, Temp, M_PI/2.)/pow(2.*M_PI,2)*2.);

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

	for(i = 0; i < 8 && Range[i] <= M_PI/2.; i++)
	{
		b = Range[i];

		F.null();
		for(j = 0; j < 24; j++)
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

	return(Answer/pow(2.*M_PI,2)*2.);
}

//long double Par[] = {g, Lambda, M, P, s}
Elements k_Int(long double Par[], int Temp, long double theta)	//Integrates the k momentum results
{
	/*if(4.*M_PI*A_INVERSE < Par[3])	//Prevents work when clearly the only possible answer is 0 as all momentum for aligned quark is greater than UV cutoff
	{
		Elements Zero(0,0,0,0,0);
		return(Zero);
	}*/

	//long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	//long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};	//Number of gamma from center
	Elements F;	//Sum of ordinates*weights
	Elements Answer(0,0,0,0,0);	//Answer to be returned
	Elements Partial;	//Answer for sub-interval for determining completeness
	Elements holder;
	long double x1, x2;	//Abscissa
	long double a = 0, b = 0;//Sub-interval limits of integration
	int Poles;	//Number of poles
	long double zero[26];	//The real part of the signular pole
	long double gamma[26];	//The distance to the singular, maybe
	long double UV_End, IR_Resume, IR_Stop;	//UV and IR boundaries in k to stop the inclusion of quark momentum prohibited by lattice space and box size.
	int i, j, l;	//Counters
	int Intervals;
	//ofstream Table("k Table", ios::app);
	//ofstream Poles_Table("k Poles", ios::app);

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);
	//for(i = 0; i < Poles; i++)
		//Poles_Table << Par[3] << " " << Par[4] << " "  << theta << " " << zero[i] << " " << gamma[i] << endl;
	long double Stops[Poles*17+12];

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(zero[i] == .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2))))	//Only true for on-shell k
			for(j = 0; j < 17; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Stops as required by poles
				l++;
			}
		else if(gamma[i] == gamma[i])	//Prevents bad poles from getting in (It would be better to find the source of bad poles and eliminate it)
			for(j = 1; j < 14; j+=4)
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
	Stops[l+1] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]-pow(2.*Par[2],2)+pow(Par[3]*cos(theta),2)));	//On-shells leaving the range 0 to E
	Stops[l+2] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]-pow(2.*Par[2],2)+pow(Par[3]*cos(theta),2)));
	Stops[l+3] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);	//Potiential leaving the range 0 to E
	Stops[l+4] = abs((pow(Par[2],2)*Par[3]*cos(theta)+sqrt((Par[4]+pow(Par[3],2))*(pow(Par[2],4)+(Par[4]+pow(Par[3]*sin(theta),2))*(Par[4]-2.*pow(Par[2],2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2))));	//On-shell leaving the range k+ to E-k-
	Stops[l+5] = abs((pow(Par[2],2)*Par[3]*cos(theta)-sqrt((Par[4]+pow(Par[3],2))*(pow(Par[2],4)+(Par[4]+pow(Par[3]*sin(theta),2))*(Par[4]-2.*pow(Par[2],2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2))));
	Stops[l+6] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]+pow(Par[3]*cos(theta),2)));	//Photon point leaving 0 to E
	Stops[l+7] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]+pow(Par[3]*cos(theta),2)));
	Stops[l+8] = .5*abs(Par[3]*cos(theta)+sqrt(3.*pow(Par[3],2)+4.*Par[4]+pow(Par[3]*cos(theta),2)));
	Stops[l+9] = .5*abs(Par[3]*cos(theta)-sqrt(3.*pow(Par[3],2)+4.*Par[4]+pow(Par[3]*cos(theta),2)));
	UV_End = Stops[l+10] = (-4.*Par[3]*cos(theta)+sqrt(pow(8.*M_PI*A_INVERSE,2)-pow(4.*Par[3]*sin(theta),2)))/8.;
	/*if(pow(16.*M_PI*A_INVERSE/LATTICE_N,2)-pow(4.*Par[3]*sin(theta),2) > 0)
	{
		IR_Stop = Stops[l+11] = (4.*Par[3]*cos(theta)-sqrt(pow(16.*M_PI*A_INVERSE/LATTICE_N,2)-pow(4.*Par[3]*sin(theta),2)))/8.;
		IR_Resume = Stops[l+12] = (4.*Par[3]*cos(theta)+sqrt(pow(16.*M_PI*A_INVERSE/LATTICE_N,2)-pow(4.*Par[3]*sin(theta),2)))/8.;
		if(IR_Stop < 0)	//If this point is less than 0, then we start at the stop
			IR_Stop = 0;
		l += 2;
	}*/

	for(i = 0; i < l+11; i++)	//Removes stops in forbiden regions and points that aren't self-equal
		if(Stops[i] != Stops[i])// || Stops[i] > UV_End || (pow(16.*M_PI*A_INVERSE/LATTICE_N,2)-pow(4.*Par[3]*sin(theta),2) > 0 && Stops[i] > IR_Stop && Stops[i] < IR_Resume))
			Stops[i] = -1;

	mergeSort(Stops, 0, l+10);

	i = 0;
	j = 0;
	while(Stops[j] <= 0)
		j++;
	for(; j < l+11; j++)
	{
		if(((i > 0 && Stops[i-1] != Stops[j]) || i == 0) && Stops[j] == Stops[j])	//Removes Duplicates and nan
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
		/*if(a == IR_Stop)
		{
			a = b = IR_Resume;
			i++;
		}*/

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
		for(l = 0; l < 24; l++)
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
	}while(!(Partial == 0) && (i < Intervals || abs(Partial/Answer) >= .0001) && a <= 20.*sqrt(Par[4]+pow(Par[3],2)));// UV_End); //k bigger than 20E is getting pretty stupid, should be sneaking up on 10^-5 of the answer left

	return(Answer);
}

void Characterize_k_Int(long double Par[], int Temp, long double theta, long double zero[], long double gamma[], int &Poles) //Returns the poles of the k integral's integrands
{
	long double holder;
	int i, j, l;

	Poles = 2;
	zero[0] = .5*Par[3]*abs(cos(theta));	//Near intersection of 2 on-shells
	gamma[0] = .05;
	zero[1] = Par[2];
	gamma[1] = Par[1];

	if((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3]*sin(theta),2)) > 0.)
	{
		zero[Poles] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//on-shell momentum, it's almost like every on-shell peak crosses every other one at this point
		gamma[Poles] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par);
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Temp, Par[5], Par)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	for(i = 0; i < Poles; i++)	//Any point that is less than is zero is replaced with zero width so that it will be removed in the bubble sort.
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

	i = 0;
	for(j = 0; j < Poles; j++)
	{
		if(((i > 0 && zero[i-1] != zero[j]) || i == 0) && zero[j] == zero[j])
		{
			zero[i] = zero[j];
			gamma[i] = gamma[j];
			i++;
		}
		else if(zero[j] > 1000)
			break;
	}
	Poles = i;

	return;
}

bool Newton_Method_k(long double& k, long double s, long double P, long double theta, long double M, long double Lambda, long double(*V)(long double, long double, long double, long double), long double(*k0)(long double, long double, long double, long double, long double))	//Returns the k-intesection of a potiential and on-shell peak
{
	long double newk;
	const long double h = 1e-4;
	bool Success = true;
	int i = 0;

	newk = k - 2.*h*(V(s, M, k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(s, M, k-h, Lambda)+V(s, M, k+h, Lambda)-k0(s, P, k+h, theta, M));

	while(abs(1.-newk/k) > 1e-5 && i <= 10)
	{
		k = newk;
		newk = k - 2.*h*(V(s, M, k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(s, M, k-h, Lambda)+V(s, M, k+h, Lambda)-k0(s, P, k+h, theta, M));
		i++;
	}

	if(abs(1.-newk/k) > 1e-2 || newk < 0 || newk != newk)	//Soundness of number and degree of improvement on last interation
		Success = false;
	if(abs(1.-V(s, M, newk, Lambda)/k0(s, P, newk, theta, M)) > 1e-2)	//Closeness to solution
		Success = false;

	k = newk;

	return(Success);
}

long double V_Plus(long double s, long double M, long double k, long double Lambda)
{
	return(.5*sqrt(complex<long double>(4.*(pow(k,2)+pow(M,2)),pow(Lambda,2))).real());
}

long double V_Minus(long double s, long double M, long double k, long double Lambda)
{
	return(-.5*sqrt(complex<long double>(4.*(pow(k,2)+pow(M,2)),pow(Lambda,2))).real());
}

long double Emm(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.-Energy(M,P/2.,-k,theta));
}

long double Epm(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.+Energy(M,P/2.,-k,theta));
}

long double mEmp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P,2))/2.-Energy(M,P/2.,k,theta));
}

long double Emp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.-Energy(M,P/2.,k,theta));
}

long double mEpp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P,2))/2.+Energy(M,P/2.,k,theta));
}

long double Epp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.+Energy(M,P/2.,k,theta));
}

long double Upper_Bound(long double s, long double P, long double k, long double theta, long double M)	//Vacuum boundaries
{
	return(sqrt(s+pow(P,2))/2.-Energy(0,P/2.,-k,theta));
}

long double Lower_Bound(long double s, long double P, long double k, long double theta, long double M)
{
	return(Energy(0,P/2.,k,theta)-sqrt(s+pow(P,2))/2.);
}

//long double Par[] = {g, Lambda, M, P, s}
Elements Folding(long double Par[], int Temp, long double k, long double theta)	//Folding integral, energy integral
{
	if(Temp == 0 && abs(sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta)-Energy(0,Par[3]/2.,k,theta)) < 1e-12)	//Let's save some time and just return 0, because it is
		return(Elements(0,0,0,0,0));
	else if(Par[4]+pow(Par[3],2) < 0)
		return(Elements(0,0,0,0,0));	//Bad data trap and time saver

	//long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	//long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
	long double Range[] = {-Boundary[7], -Boundary[6], -Boundary[5], -Boundary[4], -Boundary[3], -Boundary[2], -Boundary[1], -Boundary[0], 0, Boundary[0], Boundary[1], Boundary[2], Boundary[3], Boundary[4], Boundary[5], Boundary[6], Boundary[7]};	//Number of gamma from center
	long double a, b;	//Sub-interval limits of integration
	long double Max;	//Upper limit of integration
	Elements F;	//Sum of ordinates*weights
	Elements Answer(0,0,0,0,0);	//Results to be returned
	Elements Partial(0,0,0,0,0);//Partial Answer
	Elements holder;
	long double x1, x2;	//Abscissa
	long double zero[12];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[12];	//Imaginary part of poles
	int Intervals;		//Number of intervals required by poles and discontinuities
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	//ofstream Table("omega Table", ios::app);
	//ofstream Poles_Table("omega Poles", ios::app);

	Characterize_Folding(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	//for(i = 0; i < Poles; i++)
		//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << zero[i] << " " << gamma[i] << endl;
	long double Stops[Poles*17+6];	//Intervals that are required by integrating near poles

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
	Stops[l+2] = -Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
	Stops[l+3] = sqrt(Par[4]+pow(Par[3],2))/2.+Energy(0,Par[3]/2.,-k,theta);
	Stops[l+4] = sqrt(Par[4]+pow(Par[3],2))/2.;
	Stops[l+5] = -sqrt(Par[4]+pow(Par[3],2))/2.;

	mergeSort(Stops, 0, l+5);

	if(Temp != 0)
	{
		a = b = -sqrt(Par[4]+pow(Par[3],2))/2.;
		Max = sqrt(Par[4]+pow(Par[3],2))/2.;//zero[Poles-1]+Boundary[7]*gamma[Poles-1];
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
		Max = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta)/2.;
	}

	i = 0;
	j = 0;
	while(Stops[j] < a || Stops[j] == Stops[j]+1.)
		j++;
	for(; j < l+6; j++)
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

	if(a == a+1.)
		a = b = Stops[0];

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

		if(b > Max)
			b = Max;

		F.null();
		#pragma omp parallel for
		for(l = 0; l < 24; l++)	//Integrate the sub-interval
		{
			long double x1 = (b+a-Disp[l]*(b-a))/2.;
			long double x2 = (b+a+Disp[l]*(b-a))/2.;

			holder = (Elements(Spin_Sum1(Par, x1, k, theta), Potential1(Par,x1,k), Spin_Linear(Par, x1, k, theta)*Potential1(Par,x1,k), Spin_Quad(Par, x1, k, theta)*Potential1(Par,x1,k), Potential2(Par,x1,k))*ImFolding_Integrand1(Par,x1,k,theta,Temp));
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
			holder = (Elements(Spin_Sum1(Par, x2, k, theta), Potential1(Par,x2,k), Spin_Linear(Par, x2, k, theta)*Potential1(Par,x2,k), Spin_Quad(Par, x2, k, theta)*Potential1(Par,x2,k), Potential2(Par,x2,k))*ImFolding_Integrand1(Par,x2,k,theta,Temp));
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
		}
		holder = (Elements(Spin_Sum1(Par, (a+b)/2., k, theta), Potential1(Par,(a+b)/2.,k), Spin_Linear(Par, (a+b)/2., k, theta)*Potential1(Par,(a+b)/2.,k), Spin_Quad(Par, (a+b)/2., k, theta)*Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*ImFolding_Integrand1(Par,(a+b)/2.,k,theta,Temp));
		//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << (a+b)/2. << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
		F += holder*w[0];

		Partial = F*(b-a)/(2.);
		Answer += Partial;
		a = b;
	}while((!(Partial == 0) || a < Max) && (i < Intervals || ((abs(Partial/Answer) >= .0001))));

	return(Answer/M_PI);
}

void Characterize_Folding(long double Par[], int Temp, long double k, long double theta, long double zero[10], long double gamma[10], int &Poles)
{
	long double Lower, Upper;	//Limits of integration in Folding, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(true)//Temp != 0)
	{
		Lower = -sqrt(Par[4]+pow(Par[3],2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3],2))/2.;	//Integrate from 0 to E and twice E to infinity (ie, I need all points greater than 0)
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);
	}

	zero[0] = .5*sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)),pow(Par[1],2))).real();	//Potential poles, I know exactly where these are at.
	zero[1] = -.5*sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)),pow(Par[1],2))).real();
	gamma[0] = abs(.5*sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)),pow(Par[1],2))).imag());
	gamma[1] = abs(-.5*sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)),pow(Par[1],2))).imag());

	zero[2] = .5*(sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));	//Exact vacuum
	zero[3] = .5*(sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[4] = .5*(-sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[5] = .5*(-sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));

	if(Temp != 0)	//media estimate
	{
		zero[2] = Newton_Method_k0(zero[2], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[3] = Newton_Method_k0(zero[3], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[4] = Newton_Method_k0(zero[4], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[5] = Newton_Method_k0(zero[5], Par, k, theta, Temp, ImFolding_Integrand1);
		zero[6] = -Newton_Method_k0(-zero[3], Par, k, theta, Temp, ImFolding_Integrand2);	//Negative to counter the fact that I'm doing particle-hole in reverse order to minimize region of interest
		zero[7] = -Newton_Method_k0(-zero[4], Par, k, theta, Temp, ImFolding_Integrand2);
		zero[8] = -Newton_Method_k0(-zero[5], Par, k, theta, Temp, ImFolding_Integrand2);
		zero[9] = -Newton_Method_k0(-zero[6], Par, k, theta, Temp, ImFolding_Integrand2);

		gamma[2] = omega_Width(zero[2], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[3] = omega_Width(zero[3], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[4] = omega_Width(zero[4], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[5] = omega_Width(zero[5], Par, k, theta, Temp, ImFolding_Integrand1);
		gamma[6] = omega_Width(-zero[6], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[7] = omega_Width(-zero[7], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[8] = omega_Width(-zero[8], Par, k, theta, Temp, ImFolding_Integrand2);
		gamma[9] = omega_Width(-zero[9], Par, k, theta, Temp, ImFolding_Integrand2);
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = gamma[3] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
		gamma[4] = gamma[5] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	}

	if(Temp == 0)
		i = 5;
	else
		i = 9;

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

	while(zero[i] <= Upper && (i < 10 || (i < 6 && Temp == 0)))	//Move zeroes up to front of array, count off poles within the limits of integration
	{
		zero[j] = zero[i];
		gamma[j] = abs(gamma[i]);
		i++;
		j++;
	}
	Poles = j;

	return;
}

long double Newton_Method_k0(long double k0, long double Par[], long double k, long double theta, int Temp, long double (*Folding)(long double[], long double, long double, long double, int))
{
	long double new_k0;
	const long double h = 1e-4;
	int i = 0;

	new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp)-1./Folding(Par, k0-h, k, theta, Temp))/((1./Folding(Par, k0-h, k, theta, Temp)-2./Folding(Par, k0, k, theta, Temp)+1./Folding(Par, k0+h, k, theta, Temp)));

	while(abs(1.-new_k0/k0) > 1e-5 && i <= 10)
	{
		k0 = new_k0;
		new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp)-1./Folding(Par, k0-h, k, theta, Temp))/((1./Folding(Par, k0-h, k, theta, Temp)-2./Folding(Par, k0, k, theta, Temp)+1./Folding(Par, k0+h, k, theta, Temp)));
		i++;
	}

	return(k0);
}

long double omega_Width(long double zero, long double Par[], long double k, long double theta, int Temp, long double (*Folding)(long double[], long double, long double, long double, int))
{
	return(sqrt(abs(2e-10*Folding(Par, zero, k, theta, Temp)/(Folding(Par, zero-1e-5, k, theta, Temp)-2.*Folding(Par, zero, k, theta, Temp)+Folding(Par, zero+1e-5, k, theta, Temp)))));
}

//long double Par[] = {g, Lambda, M, P, s}
long double ImSelf_Energy(long double M, long double omega, long double k, int Temp, long double WidthFraction, long double Par[])	//Single quark self energy
{
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
			//Sigma = .585335/sqrt(pow(k,2)+pow(1.85515,2));
			Sigma = .569969/sqrt(pow(k,2)+pow(1.75236,2))+.0187484;//3.38027e-6;
			//a = 5.05953/(pow(k,2)+pow(1.11686,2))+6.09943;
			a = 12.5349/(pow(k,2)+pow(1.63711,2))+5.026;
			//b = -8.08693/(pow(k,2)+pow(2.70494,2))+3.25177;
			b = -291.579/(pow(k+15.2519,2)+pow(.0614821,2))+3.36681;
			//omega0 = sqrt(pow(1.49006+Shift,2)+pow(k,2))+.248573;
			omega0 = sqrt(pow(1.51443+Shift,2)+pow(k,2))+.232841;
			//knee = 3.84788*pow(k+1.,(long double)-.335162);
			knee = 3.78956*pow(k+1.,(long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			//Sigma = .660137/sqrt(pow(k,2)+pow(1.90299,2));
			Sigma = .625855/sqrt(pow(k,2)+pow(1.8429,2))+.0249334;//.0085596;
			//a = 2.82635/(pow(k,2)+pow(.916643,2))+4.19118;
			a = 3.3971/(pow(k,2)+pow(1.01744,2))+3.99561;
			//b = -83.3834/(pow(k,2)+pow(8.8641,2))+2.93508;
			b = -65187.5/(pow(k+3.11711,2)+pow(101.697,2))+8.15532;
			//omega0 = sqrt(pow(1.45524+Shift,2)+pow(k,2))+.247213;
			omega0 = sqrt(pow(1.5065+Shift,2)+pow(k,2))+.209135;
			//knee = 3.29189*pow(k+1.,(long double)-.575497);
			knee = 3.1568*pow(k+1.,(long double)-.624827)+.197004*(tanh((k-27.1743)/10.0192)+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			//Sigma = .670397/sqrt(pow(k,2)+pow(1.96561,2));
			Sigma = .587509/sqrt(pow(k,2)+pow(1.84447,2))+.0309251;//.0236134;
			//a = 2.42808/(pow(k,2)+pow(.840297,2))+3.42835;
			a = 2.44943/(pow(k,2)+pow(.887313,2))+3.32859;
			//b = .0167941/(pow(k,2)+pow(.47573,2))+1.70158;
			b = -4439.38/(pow(k-7.23198,2)+pow(38.9387,2))+4.55531;
			//omega0 = sqrt(pow(1.42617+Shift,2)+pow(k,2))+.258289;
			omega0 = sqrt(pow(1.47725+Shift,2)+pow(k,2))+.219181;
			//knee = 3.59947*pow(k+1.,(long double)-.710425);
			knee = 3.28564*pow(k+1.,(long double)-.721321)+.330483*(tanh((k-22.9096)/10.7139)+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			//Sigma = .592982/sqrt(pow(k,2)+pow(2.06656,2));
			Sigma = .459303/sqrt(pow(k,2)+pow(1.84321,2))+.0386564;
			//a = 2.09377/(pow(k,2)+pow(.763871,2))+2.65712;
			a = 1.79149/(pow(k,2)+pow(.764836,2))+2.66209;
			//b = .366499/(pow(k,2)+pow(1.06864,2))+1.35141;
			b = -1856.16/(pow(k-8.69519,2)+pow(26.3551,2))+3.94631;
			//omega0 = sqrt(pow(1.38555+Shift,2)+pow(k,2))+.253076;
			omega0 = sqrt(pow(1.45428+Shift,2)+pow(k,2))+.197493;
			//knee = 3.49204*pow(k+1.,(long double)-.925502);
			knee = 3.06296*pow(k+1.,(long double)-.917081)+.394833*(tanh((k-19.5932)/12.0494)+1);
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
		return(-2.*WidthFraction*M*Sigma*exp(ImSigma)+sqrt(pow(omega,2)-pow(k,2))*GAMMA);
	else
		return(-2.*WidthFraction*M*Sigma*exp(ImSigma));
}

long double ReSelf_Energy(long double M, long double omega, long double k, int Temp, long double WidthFraction, long double Par[])	//Single quark self energy
{
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
			//Sigma = .244368/sqrt(pow(k,2)+pow(1.32368,2));
			Sigma = .212571/sqrt(pow(k,2)+pow(1.17821,2))+.00762638;//5.77403e-6;
			//x0 = sqrt(pow(k,2)+pow(1.5567+Shift,2))+.279476;
			x0 = sqrt(pow(k,2)+pow(1.57536+Shift,2))+.259147;
			//x1 = sqrt(pow(k,2)+pow(1.50202+Shift,2))+.259;
			x1 = sqrt(pow(k,2)+pow(1.50194+Shift,2))+.222526;
			//gamma = .320676/sqrt(pow(k,2)+pow(1.56455,2))+.080032;
			gamma = .336699/sqrt(pow(k,2)+pow(1.87956,2))+.0651449;
			break;
		case 2://258MeV
			M_T = 1.69584;
			Shift = M-M_T;
			//Sigma = .322887/sqrt(pow(k,2)+pow(1.34236,2));
			Sigma = .307972/sqrt(pow(k,2)+pow(1.41483,2))+.0101423;//.0000230294;
			//x0 = sqrt(pow(k,2)+pow(1.54159+Shift,2))+.280535;
			x0 = sqrt(pow(k,2)+pow(1.56476+Shift,2))+.251031;
			//x1 = sqrt(pow(k,2)+pow(1.46598+Shift,2))+.260561;
			x1 = sqrt(pow(k,2)+pow(1.50194+Shift,2))+.222526;
			//gamma = .694901/sqrt(pow(k,2)+pow(2.13185,2))+.0653795;
			gamma = .550628/sqrt(pow(k,2)+pow(2.43968,2))+.0981269;
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			//Sigma = .375163/sqrt(pow(k,2)+pow(1.41612,2));
			Sigma = .339131/sqrt(pow(k,2)+pow(1.43308,2))+.0125796;//.00623679;
			//x0 = sqrt(pow(k,2)+pow(1.45507+Shift,2))+.337448;
			x0 = sqrt(pow(k,2)+pow(1.55034+Shift,2))+.257788;
			//x1 = sqrt(pow(k,2)+pow(1.40846+Shift,2))+.289292;
			x1 = sqrt(pow(k,2)+pow(1.46999+Shift,2))+.231821;
			//gamma = .690491/sqrt(pow(k,2)+pow(1.97525,2))+.141465;
			gamma = .615278/sqrt(pow(k,2)+pow(2.22298,2))+.143376;
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			//Sigma = .370623/sqrt(pow(k,2)+pow(1.53585,2));
			Sigma = .304841/sqrt(pow(k,2)+pow(1.42911,2))+.0157245;
			//x0 = sqrt(pow(k,2)+pow(1.39619+Shift,2))+.35548;
			x0 = sqrt(pow(k,2)+pow(1.55511+Shift,2))+.231105;
			//x1 = sqrt(pow(k,2)+pow(1.3481+Shift,2))+.296587;
			x1 = sqrt(pow(k,2)+pow(1.44714+Shift,2))+.20956;
			//gamma = .857781/sqrt(pow(k,2)+pow(2.25072,2))+.196022;
			gamma = .862629/sqrt(pow(k,2)+pow(2.67193,2))+.189598;
			break;
		default:
			Sigma = 0;
			x0 = 1;
			x1 = 1;
			gamma = 1;
			break;
	}

	return(WidthFraction*Sigma*(omega-x0)/(pow(omega-x1,2)+gamma));
}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	if(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta) < 0)
		return(0.);
	else
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
			Temp = .194;
			break;
		case 2:
			Temp = .258;
			break;
		case 3:
			Temp = .32;
			break;
		case 4:
			Temp = .4;
			break;
		default:
			return(0);
	}

	return(1./(1.+exp(omega/Temp)));
}

long double Potential1(long double Par[], long double k0, long double k)        //Potiential for the numerator of the boson spectrum
{
#if VERSION == 22
	return(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))));
#elif VERSION == 42
	return(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)));
#elif VERSION == 24
	return(pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#endif
}

long double Potential2(long double Par[], long double k0, long double k)        //Potiential for the denominator of the T-Matrix and boson spectrum
{
#if VERSION == 22
	return(Par[0]*pow(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#elif VERSION == 42
	return(Par[0]*pow(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)),2));
#elif VERSION == 24
	return(Par[0]*pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),4));
#endif
}

long double ImD(long double omega, long double k, long double M, int Temp, long double Par[])	//Single quark spectral function
{
	return(2.*ImSelf_Energy(M, omega, k, Temp, 1, Par)/(pow(pow(omega,2)-pow(k,2)-pow(M,2)-2.*M*ReSelf_Energy(M, omega, k, Temp, 1, Par),2)+pow(ImSelf_Energy(M, omega, k, Temp, 1, Par),2)));
}

long double Spin_Sum1(long double Par[], long double k0, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), strictly pseudoscalar for now
{
	return((Par[4]/4.+pow(k,2)-pow(k0,2)+pow(Par[2],2))/(pow(Par[2],2)));
}

long double Spin_Linear(long double Par[], long double k0, long double k , long double theta)
{
	if(Par[4] >= 0.)
		return(sqrt(3.*Par[4]/(8.*pow(Par[2],2))));
	else
		return(0.);
}

long double Spin_Quad(long double Par[], long double k0, long double k , long double theta)
{
	return(Par[4]/4.-pow(k0,2)+pow(k,2))/(2.*pow(Par[2],2));
}

long double Spin_Sum2(long double Par[], long double k0, long double k , long double theta)	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly pseudoscalar for now
{
	return((Par[4]/4.+pow(k,2)-pow(k0,2)-pow(Par[2],2))/pow(Par[2],2));
}

long double ImFolding_Integrand1(long double Par[], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	return(-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp, Par)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp, Par)*(1.-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)));
}

long double ImFolding_Integrand2(long double Par[], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for negitive energy (anti-particle/particle-hole)
{
	return(-pow(Par[2],2)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Energy(0, Par[3]/2., k, theta), Par[2], Temp, Par)*ImD(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Energy(0, Par[3]/2., -k, theta), Par[2], Temp, Par)*(Fermi(sqrt(Par[4]+pow(Par[3],2))/2.-k0, Temp)-Fermi(sqrt(Par[4]+pow(Par[3],2))/2.+k0, Temp)));
}
