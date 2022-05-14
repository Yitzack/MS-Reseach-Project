//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
#include<cfloat>
#include<fstream>
#include<iostream>
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
void ImSelf_Energy(long double, long double, long double[], long double[],int, long double[]);	//Both imaginary single quark self energies
void Self_Energy(long double, long double, long double[], long double[],int, long double[], long double[]);	//Both single quark self energies
long double ImSelf_Energy(long double, long double, long double, long double[], int); //Imaginary single quark self energy
void ReSelf_Energy(long double, long double, long double[], int, long double[]);	//Real single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Set_Temp(int);
long double Potential1(long double[], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double Spin_Sum1(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Linear(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Quad(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Spin_Sum2(long double[], long double, long double, long double);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double ImFolding_Integrand(long double[], long double, long double, long double, int);	//Integrand of the folding integral for positive energy

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
char* Process;

//long double Par[] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[], int Temp)
{
#if ORDER == 37
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
#elif ORDER == 97
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
#endif
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
#if ORDER == 37
		for(j = 0; j < 9; j++)
#elif ORDER == 97
		for(j = 0; j < 24; j++)
#endif
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

#if ORDER == 37
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
#elif ORDER == 97
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
#endif
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
		//Poles_Table << Par[3] << " " << Par[4] << " " << theta << " " << zero[i] << " " << gamma[i] << endl;
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
#if ORDER == 37
		for(l = 0; l < 9; l++)
#elif ORDER == 97
		for(l = 0; l < 24; l++)
#endif
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
		gamma[Poles] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp);
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
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
#if ORDER == 37
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp
#elif ORDER == 97
	long double Disp[] = {0.06342068498268678602883,0.1265859972696720510680,0.1892415924618135864853,0.2511351786125772735072,0.3120175321197487622079,0.3716435012622848888637,0.4297729933415765246586,0.4861719414524920421770,0.5406132469917260665582,0.5928776941089007124559,0.6427548324192376640569,0.6900438244251321135048,0.7345542542374026962137,0.7761068943454466350181,0.8145344273598554315395,0.8496821198441657010349,0.8814084455730089100370,0.9095856558280732852130,0.9341002947558101490590,0.9548536586741372335552,0.9717622009015553801400,0.9847578959142130043593,0.9937886619441677907601,0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825,0.06333550929649174859084,0.06295270746519569947440,0.06231641732005726740108,0.06142920097919293629683,0.06029463095315201730311,0.05891727576002726602453,0.05730268153018747548516,0.05545734967480358869043,0.05338871070825896852794,0.05110509433014459067462,0.04861569588782824027765,0.04593053935559585354250,0.04306043698125959798835,0.04001694576637302136861,0.03681232096300068981947,0.03345946679162217434249,0.02997188462058382535069,0.02636361892706601696095,0.02264920158744667649877,0.01884359585308945844445,0.01496214493562465102958,0.01102055103159358049751,0.007035099590086451473451,0.003027278988922905077481};	//Weight of the function at Disp
#endif
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
#if ORDER == 37
		for(l = 0; l < 9; l++) //Integrate the sub-interval
#elif ORDER == 97
		for(l = 0; l < 24; l++) //Integrate the sub-interval
#endif
		{
			long double x1 = (b+a-Disp[l]*(b-a))/2.;
			long double x2 = (b+a+Disp[l]*(b-a))/2.;

			holder = (Elements(Spin_Sum1(Par, x1, k, theta), Potential1(Par,x1,k), Spin_Linear(Par, x1, k, theta)*Potential1(Par,x1,k), Spin_Quad(Par, x1, k, theta)*Potential1(Par,x1,k), Potential2(Par,x1,k))*ImFolding_Integrand(Par,x1,k,theta,Temp));
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x1 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
			holder = (Elements(Spin_Sum1(Par, x2, k, theta), Potential1(Par,x2,k), Spin_Linear(Par, x2, k, theta)*Potential1(Par,x2,k), Spin_Quad(Par, x2, k, theta)*Potential1(Par,x2,k), Potential2(Par,x2,k))*ImFolding_Integrand(Par,x2,k,theta,Temp));
			//Table << Par[3] << " " << Par[4] << " " << theta << " " << k << " " << x2 << " " << holder.store(0) << " " << holder.store(1) << " " << holder.store(2) << endl;
			F += holder*w[l+1];
		}
		holder = (Elements(Spin_Sum1(Par, (a+b)/2., k, theta), Potential1(Par,(a+b)/2.,k), Spin_Linear(Par, (a+b)/2., k, theta)*Potential1(Par,(a+b)/2.,k), Spin_Quad(Par, (a+b)/2., k, theta)*Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*ImFolding_Integrand(Par,(a+b)/2.,k,theta,Temp));
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

	zero[0] = .5*sqrt(complex<long double>(4.*pow(k,2),pow(Par[1],2))).real();	//Potential poles, I know exactly where these are at.
	zero[1] = -.5*sqrt(complex<long double>(4.*pow(k,2),pow(Par[1],2))).real();
	gamma[0] = abs(.5*sqrt(complex<long double>(4.*pow(k,2),pow(Par[1],2))).imag());
	gamma[1] = abs(-.5*sqrt(complex<long double>(4.*pow(k,2),pow(Par[1],2))).imag());

	zero[2] = .5*(sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));	//Exact vacuum
	zero[3] = .5*(sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[4] = .5*(-sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[5] = .5*(-sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));

	if(Temp != 0)	//media estimate
	{
		zero[2] = Newton_Method_k0(zero[2], Par, k, theta, Temp, ImFolding_Integrand);
		zero[3] = Newton_Method_k0(zero[3], Par, k, theta, Temp, ImFolding_Integrand);
		zero[4] = Newton_Method_k0(zero[4], Par, k, theta, Temp, ImFolding_Integrand);
		zero[5] = Newton_Method_k0(zero[5], Par, k, theta, Temp, ImFolding_Integrand);

		gamma[2] = omega_Width(zero[2], Par, k, theta, Temp, ImFolding_Integrand);
		gamma[3] = omega_Width(zero[3], Par, k, theta, Temp, ImFolding_Integrand);
		gamma[4] = omega_Width(zero[4], Par, k, theta, Temp, ImFolding_Integrand);
		gamma[5] = omega_Width(zero[5], Par, k, theta, Temp, ImFolding_Integrand);
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = gamma[3] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
		gamma[4] = gamma[5] = abs(.5*imag(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	}

	for(i = 5; i >= 0; i--)	//Bubble sort
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

	while(zero[i] <= Upper && i < 6)	//Move zeroes up to front of array, count off poles within the limits of integration
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
void ImSelf_Energy(long double M, long double omega[], long double k[], long double Par[], int Temp, long double Results[])	//Single quark self energy
{
	static long double omega0[2];	//location of central peak
	static long double Sigma[2];	//size of energy dependance
	static long double a[2], b[2];	//slope of exponential decrease to left and right
	static long double knee[2];	//space to change from left to right side of peak
	static long double M_T, Shift;
	static long double k_old[2];

	/*static bool Been_Here = false;
	static char File_Name[70] = "Self_E_Modifier.";
	static long double Mod[26];
	if(!Been_Here)
	{
		char TempChar[1] = {Temp+48};
		strncat(File_Name,TempChar,1);
		strncat(File_Name,".",1);
		strcat(File_Name,Process);
		static ifstream File(File_Name);
		for(int i = 0; i < 26; i++)
			File >> Mod[i];
		Been_Here = true;
	}*/

	if(pow(omega[0],2)>=pow(k[0],2))
		Results[0] = sqrt(pow(omega[0],2)-pow(k[0],2))*GAMMA;
		//Results[0] = GAMMA*sqrt(Par[4])/3.0404;
	else
		Results[0] = 0;
	if(pow(omega[1],2)>=pow(k[1],2))
		Results[1] = sqrt(pow(omega[1],2)-pow(k[1],2))*GAMMA;
		//Results[1] = GAMMA*sqrt(Par[4])/3.0404;
	else
		Results[1] = 0;

	if(Temp == 0)
		return;

	if(k[0] != k_old[0] || k[1] != k_old[1])
	{
		k_old[0] = k[0];
		k_old[1] = k[1];
		switch(Temp)
		{
			/*case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(pow(k[0],2)+pow(1.75236,2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1],2)+pow(1.75236,2))+.0187484;
				a[0] = 4.689/(pow(k[0],2)+pow(1.18,2))+4.59495;
				a[1] = 4.689/(pow(k[1],2)+pow(1.18,2))+4.59495;
				b[0] = -70400/(pow(k[0]+20,2)+pow(130,2))+6.24;
				b[1] = -70400/(pow(k[1]+20,2)+pow(130,2))+6.24;
				omega0[0] = sqrt(pow(1.51443+Shift,2)+pow(k[0],2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift,2)+pow(k[1],2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1.,(long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1.,(long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(pow(k[0],2)+pow(1.75236,2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1],2)+pow(1.75236,2))+.0187484;
				a[0] = 12.5349/(pow(k[0],2)+pow(1.63711,2))+5.026;
				a[1] = 12.5349/(pow(k[1],2)+pow(1.63711,2))+5.026;
				b[0] = -291.579/(pow(k[0]+15.2519,2)+pow(.0614821,2))+3.36681;
				b[1] = -291.579/(pow(k[1]+15.2519,2)+pow(.0614821,2))+3.36681;
				omega0[0] = sqrt(pow(1.51443+Shift,2)+pow(k[0],2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift,2)+pow(k[1],2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1.,(long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1.,(long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;
			case 2://285MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .625855/sqrt(pow(k[0],2)+pow(1.8429,2))+.0249334;
				Sigma[1] = .625855/sqrt(pow(k[1],2)+pow(1.8429,2))+.0249334;
				a[0] = 3.3971/(pow(k[0],2)+pow(1.01744,2))+3.99561;
				a[1] = 3.3971/(pow(k[1],2)+pow(1.01744,2))+3.99561;
				b[0] = -65187.5/(pow(k[0]+3.11711,2)+pow(101.697,2))+8.15532;
				b[1] = -65187.5/(pow(k[1]+3.11711,2)+pow(101.697,2))+8.15532;
				omega0[0] = sqrt(pow(1.5065+Shift,2)+pow(k[0],2))+.209135;
				omega0[1] = sqrt(pow(1.5065+Shift,2)+pow(k[1],2))+.209135;
				knee[0] = 3.1568*pow(k[0]+1.,(long double)-.624827)+.197004*(tanh((k[0]-27.1743)/10.0192)+1);
				knee[1] = 3.1568*pow(k[1]+1.,(long double)-.624827)+.197004*(tanh((k[1]-27.1743)/10.0192)+1);
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .587509/sqrt(pow(k[0],2)+pow(1.84447,2))+.0309251;
				Sigma[1] = .587509/sqrt(pow(k[1],2)+pow(1.84447,2))+.0309251;
				a[0] = 2.44943/(pow(k[0],2)+pow(.887313,2))+3.32859;
				a[1] = 2.44943/(pow(k[1],2)+pow(.887313,2))+3.32859;
				b[0] = -4439.38/(pow(k[0]-7.23198,2)+pow(38.9387,2))+4.55531;
				b[1] = -4439.38/(pow(k[1]-7.23198,2)+pow(38.9387,2))+4.55531;
				omega0[0] = sqrt(pow(1.47725+Shift,2)+pow(k[0],2))+.219181;
				omega0[1] = sqrt(pow(1.47725+Shift,2)+pow(k[1],2))+.219181;
				knee[0] = 3.28564*pow(k[0]+1.,(long double)-.721321)+.330483*(tanh((k[0]-22.9096)/10.7139)+1);
				knee[1] = 3.28564*pow(k[1]+1.,(long double)-.721321)+.330483*(tanh((k[1]-22.9096)/10.7139)+1);
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .459303/sqrt(pow(k[0],2)+pow(1.84321,2))+.0386564;
				Sigma[1] = .459303/sqrt(pow(k[1],2)+pow(1.84321,2))+.0386564;
				a[0] = 1.79149/(pow(k[0],2)+pow(.764836,2))+2.66209;
				a[1] = 1.79149/(pow(k[1],2)+pow(.764836,2))+2.66209;
				b[0] = -1856.16/(pow(k[0]-8.69519,2)+pow(26.3551,2))+3.94631;
				b[1] = -1856.16/(pow(k[1]-8.69519,2)+pow(26.3551,2))+3.94631;
				omega0[0] = sqrt(pow(1.45428+Shift,2)+pow(k[0],2))+.197493;
				omega0[1] = sqrt(pow(1.45428+Shift,2)+pow(k[1],2))+.197493;
				knee[0] = 3.06296*pow(k[0]+1.,(long double)-.917081)+.394833*(tanh((k[0]-19.5932)/12.0494)+1);
				knee[1] = 3.06296*pow(k[1]+1.,(long double)-.917081)+.394833*(tanh((k[1]-19.5932)/12.0494)+1);
				break;
			/*case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969*Mod[0]/sqrt(pow(k[0],2)+pow(1.75236*Mod[1],2))+.0187484*Mod[2];
				Sigma[1] = .569969*Mod[0]/sqrt(pow(k[1],2)+pow(1.75236*Mod[1],2))+.0187484*Mod[2];
				a[0] = 12.5349*Mod[3]/(pow(k[0],2)+pow(1.63711*Mod[4],2))+5.026*Mod[5];
				a[1] = 12.5349*Mod[3]/(pow(k[1],2)+pow(1.63711*Mod[4],2))+5.026*Mod[5];
				b[0] = -291.579*Mod[6]/(pow(k[0]+15.2519*Mod[7],2)+pow(.0614821*Mod[8],2))+3.36681*Mod[9];
				b[1] = -291.579*Mod[6]/(pow(k[1]+15.2519*Mod[7],2)+pow(.0614821*Mod[8],2))+3.36681*Mod[9];
				omega0[0] = sqrt(pow(1.51443*Mod[10]+Shift,2)+pow(k[0],2))+.232841*Mod[11];
				omega0[1] = sqrt(pow(1.51443*Mod[10]+Shift,2)+pow(k[1],2))+.232841*Mod[11];
				knee[0] = 3.78956*Mod[12]*pow(k[0]+1.,(long double)-.530289*Mod[13])+.305*Mod[14]*(tanh((k[0]-48.4*Mod[15])/11.1111/Mod[16])+1);
				knee[1] = 3.78956*Mod[12]*pow(k[1]+1.,(long double)-.530289*Mod[13])+.305*Mod[14]*(tanh((k[1]-48.4*Mod[15])/11.1111/Mod[16])+1);
				break;
			case 2://285MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .625855*Mod[0]/sqrt(pow(k[0],2)+pow(1.8429*Mod[1],2))+.0249334*Mod[2];
				Sigma[1] = .625855*Mod[0]/sqrt(pow(k[1],2)+pow(1.8429*Mod[1],2))+.0249334*Mod[2];
				a[0] = 3.3971*Mod[3]/(pow(k[0],2)+pow(1.01744*Mod[4],2))+3.99561*Mod[5];
				a[1] = 3.3971*Mod[3]/(pow(k[1],2)+pow(1.01744*Mod[4],2))+3.99561*Mod[5];
				b[0] = -65187.5*Mod[6]/(pow(k[0]+3.11711*Mod[7],2)+pow(101.697*Mod[8],2))+8.15532*Mod[9];
				b[1] = -65187.5*Mod[6]/(pow(k[1]+3.11711*Mod[7],2)+pow(101.697*Mod[8],2))+8.15532*Mod[9];
				omega0[0] = sqrt(pow(1.5065*Mod[10]+Shift,2)+pow(k[0],2))+.209135*Mod[11];
				omega0[1] = sqrt(pow(1.5065*Mod[10]+Shift,2)+pow(k[1],2))+.209135*Mod[11];
				knee[0] = 3.1568*Mod[12]*pow(k[0]+1.,(long double)-.624827*Mod[13])+.197004*Mod[14]*(tanh((k[0]-27.1743*Mod[15])/10.0192/Mod[16])+1);
				knee[1] = 3.1568*Mod[12]*pow(k[1]+1.,(long double)-.624827*Mod[13])+.197004*Mod[14]*(tanh((k[1]-27.1743*Mod[15])/10.0192/Mod[16])+1);
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .587509*Mod[0]/sqrt(pow(k[0],2)+pow(1.84447*Mod[1],2))+.0309251*Mod[2];
				Sigma[1] = .587509*Mod[0]/sqrt(pow(k[1],2)+pow(1.84447*Mod[1],2))+.0309251*Mod[2];
				a[0] = 2.44943*Mod[3]/(pow(k[0],2)+pow(.887313*Mod[4],2))+3.32859*Mod[5];
				a[1] = 2.44943*Mod[3]/(pow(k[1],2)+pow(.887313*Mod[4],2))+3.32859*Mod[5];
				b[0] = -4439.38*Mod[6]/(pow(k[0]-7.23198*Mod[7],2)+pow(38.9387*Mod[8],2))+4.55531*Mod[9];
				b[1] = -4439.38*Mod[6]/(pow(k[1]-7.23198*Mod[7],2)+pow(38.9387*Mod[8],2))+4.55531*Mod[9];
				omega0[0] = sqrt(pow(1.47725*Mod[10]+Shift,2)+pow(k[0],2))+.219181*Mod[11];
				omega0[1] = sqrt(pow(1.47725*Mod[10]+Shift,2)+pow(k[1],2))+.219181*Mod[11];
				knee[0] = 3.28564*Mod[12]*pow(k[0]+1.,(long double)-.721321*Mod[13])+.330483*Mod[14]*(tanh((k[0]-22.9096*Mod[15])/10.7139/Mod[16])+1);
				knee[1] = 3.28564*Mod[12]*pow(k[1]+1.,(long double)-.721321*Mod[13])+.330483*Mod[14]*(tanh((k[1]-22.9096*Mod[15])/10.7139/Mod[16])+1);
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .459303*Mod[0]/sqrt(pow(k[0],2)+pow(1.84321*Mod[1],2))+.0386564*Mod[2];
				Sigma[1] = .459303*Mod[0]/sqrt(pow(k[1],2)+pow(1.84321*Mod[1],2))+.0386564*Mod[2];
				a[0] = 1.79149*Mod[3]/(pow(k[0],2)+pow(.764836*Mod[4],2))+2.66209*Mod[5];
				a[1] = 1.79149*Mod[3]/(pow(k[1],2)+pow(.764836*Mod[4],2))+2.66209*Mod[5];
				b[0] = -1856.16*Mod[6]/(pow(k[0]-8.69519*Mod[7],2)+pow(26.3551*Mod[8],2))+3.94631*Mod[9];
				b[1] = -1856.16*Mod[6]/(pow(k[1]-8.69519*Mod[7],2)+pow(26.3551*Mod[8],2))+3.94631*Mod[9];
				omega0[0] = sqrt(pow(1.45428*Mod[10]+Shift,2)+pow(k[0],2))+.197493*Mod[11];
				omega0[1] = sqrt(pow(1.45428*Mod[10]+Shift,2)+pow(k[1],2))+.197493*Mod[11];
				knee[0] = 3.06296*Mod[12]*pow(k[0]+1.,(long double)-.917081*Mod[13])+.394833*Mod[14]*(tanh((k[0]-19.5932*Mod[15])/12.0494/Mod[16])+1);
				knee[1] = 3.06296*Mod[12]*pow(k[1]+1.,(long double)-.917081*Mod[13])+.394833*Mod[14]*(tanh((k[1]-19.5932*Mod[15])/12.0494/Mod[16])+1);
				break;*/
			case 5:
				M_T = 1.8;
				Shift = M-M_T;
				Sigma[0] = .00386564;
				Sigma[1] = .00386564;
				a[0] = 6.2;
				a[1] = 6.2;
				b[0] = 2.8;
				b[1] = 2.8;
				omega0[0] = sqrt(pow(1.53+Shift,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.53+Shift,2)+pow(k[1],2));
				knee[0] = .56;
				knee[1] = .56;
				break;
		}
	}

	long double ImSigma[2];	//Calculation of the argument to the exponential, these first 2 are approximations to hopefully avoid catastrophic loss of precision
	if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] < -4.)
		ImSigma[0] = a[0]*(omega[0]-omega0[0]+knee[0]/sqrt(a[0]*b[0]));
	else if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] > 4.)
		ImSigma[0] = b[0]*(omega0[0]-omega[0]+knee[0]/sqrt(a[0]*b[0]));
	else
		ImSigma[0] = -.5*((a[0]-b[0])*omega0[0]-((a[0]+b[0])*knee[0])/sqrt(a[0]*b[0]))+(a[0]-b[0])*omega[0]/2-sqrt(pow(((a[0]+b[0])/2.)*(omega[0]-omega0[0]+((a[0]-b[0])*knee[0])/(sqrt(a[0]*b[0])*(a[0]+b[0]))),2)+pow(knee[0],2));

	if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] < -4.)
		ImSigma[1] = a[1]*(omega[1]-omega0[1]+knee[1]/sqrt(a[1]*b[1]));
	else if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] > 4.)
		ImSigma[1] = b[1]*(omega0[1]-omega[1]+knee[1]/sqrt(a[1]*b[1]));
	else
		ImSigma[1] = -.5*((a[1]-b[1])*omega0[1]-((a[1]+b[1])*knee[1])/sqrt(a[1]*b[1]))+(a[1]-b[1])*omega[1]/2-sqrt(pow(((a[1]+b[1])/2.)*(omega[1]-omega0[1]+((a[1]-b[1])*knee[1])/(sqrt(a[1]*b[1])*(a[1]+b[1]))),2)+pow(knee[1],2));

	Results[0] += -2.*M*Sigma[0]*exp(ImSigma[0]);
	Results[1] += -2.*M*Sigma[1]*exp(ImSigma[1]);

	return;
}

void Self_Energy(long double M, long double omega[], long double k[], long double Par[], int Temp, long double ImSelf[], long double ReSelf[])	//Single quark self energy
{
	static long double omega0[2];	//location of central peak
	static long double Sigma[2];	//size of energy dependance
	static long double gamma[2];	//space to change from left to right side of peak
	static long double k_old[2];

	if(pow(omega[0],2)>=pow(k[0],2))
		ImSelf[0] = sqrt(pow(omega[0],2)-pow(k[0],2))*GAMMA;
	else
		ImSelf[0] = 0;
	if(pow(omega[1],2)>=pow(k[1],2))
		ImSelf[1] = sqrt(pow(omega[1],2)-pow(k[1],2))*GAMMA;
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
				Sigma[0] = .840172/sqrt(pow(k[0],2)+pow(1.45603,2))+.021257;
				Sigma[1] = .840172/sqrt(pow(k[1],2)+pow(1.45603,2))+.021257;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.99829,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.99829,2)+pow(k[1],2));
				gamma[0] = 1.05035*pow(k[0]+1.3891,(long double)-1.3891)+.01;
				gamma[1] = 1.05035*pow(k[1]+1.3891,(long double)-1.3891)+.01;
				break;
			case 2://285MeV
				Sigma[0] = 1.05337/sqrt(pow(k[0],2)+pow(1.50861,2))+.0282696;
				Sigma[1] = 1.05337/sqrt(pow(k[1],2)+pow(1.50861,2))+.0282696;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.97732,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.97732,2)+pow(k[1],2));
				gamma[0] = 1.4624*pow(k[0]+2.64,(long double)-1.41048)+.01;
				gamma[1] = 1.4624*pow(k[1]+2.64,(long double)-1.41048)+.01;
				break;
			case 3://320MeV
				Sigma[0] = 1.14064/sqrt(pow(k[0],2)+pow(1.54999,2))+.0350631;
				Sigma[1] = 1.14064/sqrt(pow(k[1],2)+pow(1.54999,2))+.0350631;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.96823,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.96823,2)+pow(k[1],2));
				gamma[0] = 2.07102*pow(k[0]+3.037,(long double)-1.46076)+.01;
				gamma[1] = 2.07102*pow(k[1]+3.037,(long double)-1.46076)+.01;
				break;
			case 4://400MeV
				Sigma[0] = 1.06073/sqrt(pow(k[0],2)+pow(1.64912,2))+.0438288;
				Sigma[1] = 1.06073/sqrt(pow(k[1],2)+pow(1.64912,2))+.0438288;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.93309,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.93309,2)+pow(k[1],2));
				gamma[0] = 3.42222*pow(k[0]+3.663,(long double)-1.56165)+.01;
				gamma[1] = 3.42222*pow(k[1]+3.663,(long double)-1.56165)+.01;
				break;
		}
	}

	ImSelf[0] += -2.*M*Sigma[0]*omega[0]*gamma[0]/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)));
	ImSelf[1] += -2.*M*Sigma[1]*omega[1]*gamma[1]/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)));
	ReSelf[0] += Sigma[0]*(omega[0]-omega0[0])/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)));
	ReSelf[1] += Sigma[1]*(omega[1]-omega0[1])/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)));

	return;
}

long double ImSelf_Energy(long double M, long double omega, long double k, long double Par[], int Temp)	//Single quark self energy
{
	/*long double omega0;	//location of central peak
	long double Sigma;	//size of energy dependance
	long double gamma;	//space to change from left to right side of peak
	long double ImSelf;

	if(pow(omega,2)>=pow(k,2))
		ImSelf = sqrt(pow(omega,2)-pow(k,2))*GAMMA;
	else
		ImSelf = 0;

	if(Temp == 0)
		return(ImSelf);

	switch(Temp)
	{
		case 1://194MeV
			Sigma = .840172/sqrt(pow(k,2)+pow(1.45603,2))+.021257;
			//omega0 = sqrt(pow(M,2)+pow(k,2));
			omega0 = sqrt(pow(1.99829,2)+pow(k,2));
			gamma = 1.05035*pow(k+2.343,(long double)-1.3891)+.01;
			break;
		case 2://285MeV
			Sigma = 1.05337/sqrt(pow(k,2)+pow(1.50861,2))+.0282696;
			//omega0 = sqrt(pow(M,2)+pow(k,2));
			omega0 = sqrt(pow(1.97732,2)+pow(k,2));
			gamma = 1.4624*pow(k+2.64,(long double)-1.41048)+.01;
			break;
		case 3://320MeV
			Sigma = 1.14064/sqrt(pow(k,2)+pow(1.54999,2))+.0350631;
			//omega0 = sqrt(pow(M,2)+pow(k,2));
			omega0 = sqrt(pow(1.96823,2)+pow(k,2));
			gamma = 2.07102*pow(k+3.037,(long double)-1.46076)+.01;
			break;
		case 4://400MeV
			Sigma = 1.06073/sqrt(pow(k,2)+pow(1.64912,2))+.0438288;
			//omega0 = sqrt(pow(M,2)+pow(k,2));
			omega0 = sqrt(pow(1.93309,2)+pow(k,2));
			gamma = 3.42222*pow(k+3.663,(long double)-1.56165)+.01;
			break;
	}

	ImSelf += -M*Sigma*omega*gamma/(M_PI*(pow(omega-omega0,2)+pow(omega*gamma,2)));
	return(ImSelf);*/

	long double omega0;	//location of central peak
	long double Sigma;	//size of energy dependance
	long double a, b;	//slope of exponential decrease to left and right
	long double knee;	//space to change from left to right side of peak
	long double M_T, Shift;
	long double answer;

	/*static bool Been_Here = false;
	static char File_Name[70] = "Self_E_Modifier.";
	static long double Mod[26];
	if(!Been_Here)
	{
		char TempChar[1] = {Temp+48};
		strncat(File_Name,TempChar,1);
		strncat(File_Name,".",1);
		strcat(File_Name,Process);
		static ifstream File(File_Name);
		for(int i = 0; i < 26; i++)
			File >> Mod[i];
		Been_Here = true;
	}*/

	if(pow(omega,2)>=pow(k,2))
		answer = sqrt(pow(omega,2)-pow(k,2))*GAMMA;
		//answer = GAMMA*sqrt(Par[4])/3.0404;
	else
		answer = 0;

	if(Temp == 0)
		return(answer);

	switch(Temp)
	{
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k,2)+pow(1.75236,2))+.0187484;
			a = 4.689/(pow(k,2)+pow(1.18,2))+4.59495;
			b = -70400/(pow(k+20,2)+pow(130,2))+6.24;
			omega0 = sqrt(pow(1.51443+Shift,2)+pow(k,2))+.232841;
			knee = 3.78956*pow(k+1.,(long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;*/
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k,2)+pow(1.75236,2))+.0187484;
			a = 12.5349/(pow(k,2)+pow(1.63711,2))+5.026;
			b = -291.579/(pow(k+15.2519,2)+pow(.0614821,2))+3.36681;
			omega0 = sqrt(pow(1.51443+Shift,2)+pow(k,2))+.232841;
			knee = 3.78956*pow(k+1.,(long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .625855/sqrt(pow(k,2)+pow(1.8429,2))+.0249334;
			a = 3.3971/(pow(k,2)+pow(1.01744,2))+3.99561;
			b = -65187.5/(pow(k+3.11711,2)+pow(101.697,2))+8.15532;
			omega0 = sqrt(pow(1.5065+Shift,2)+pow(k,2))+.209135;
			knee = 3.1568*pow(k+1.,(long double)-.624827)+.197004*(tanh((k-27.1743)/10.0192)+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .587509/sqrt(pow(k,2)+pow(1.84447,2))+.0309251;
			a = 2.44943/(pow(k,2)+pow(.887313,2))+3.32859;
			b = -4439.38/(pow(k-7.23198,2)+pow(38.9387,2))+4.55531;
			omega0 = sqrt(pow(1.47725+Shift,2)+pow(k,2))+.219181;
			knee = 3.28564*pow(k+1.,(long double)-.721321)+.330483*(tanh((k-22.9096)/10.7139)+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .459303/sqrt(pow(k,2)+pow(1.84321,2))+.0386564;
			a = 1.79149/(pow(k,2)+pow(.764836,2))+2.66209;
			b = -1856.16/(pow(k-8.69519,2)+pow(26.3551,2))+3.94631;
			omega0 = sqrt(pow(1.45428+Shift,2)+pow(k,2))+.197493;
			knee = 3.06296*pow(k+1.,(long double)-.917081)+.394833*(tanh((k-19.5932)/12.0494)+1);
			break;
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969*Mod[0]/sqrt(pow(k,2)+pow(1.75236*Mod[1],2))+.0187484*Mod[2];
			a = 12.5349*Mod[3]/(pow(k,2)+pow(1.63711*Mod[4],2))+5.026*Mod[5];
			b = -291.579*Mod[6]/(pow(k+15.2519*Mod[7],2)+pow(.0614821*Mod[8],2))+3.36681*Mod[9];
			omega0 = sqrt(pow(1.51443*Mod[10]+Shift,2)+pow(k,2))+.232841*Mod[11];
			knee = 3.78956*Mod[12]*pow(k+1.,(long double)-.530289*Mod[13])+.305*Mod[14]*(tanh((k-48.4*Mod[15])/11.1111/Mod[16])+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .625855*Mod[0]/sqrt(pow(k,2)+pow(1.8429*Mod[1],2))+.0249334*Mod[2];
			a = 3.3971*Mod[3]/(pow(k,2)+pow(1.01744*Mod[4],2))+3.99561*Mod[5];
			b = -65187.5*Mod[6]/(pow(k+3.11711*Mod[7],2)+pow(101.697*Mod[8],2))+8.15532*Mod[9];
			omega0 = sqrt(pow(1.5065*Mod[10]+Shift,2)+pow(k,2))+.209135*Mod[11];
			knee = 3.1568*Mod[12]*pow(k+1.,(long double)-.624827*Mod[13])+.197004*Mod[14]*(tanh((k-27.1743*Mod[15])/10.0192/Mod[16])+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .587509*Mod[0]/sqrt(pow(k,2)+pow(1.84447*Mod[1],2))+.0309251*Mod[2];
			a = 2.44943*Mod[3]/(pow(k,2)+pow(.887313*Mod[4],2))+3.32859*Mod[5];
			b = -4439.38*Mod[6]/(pow(k-7.23198*Mod[7],2)+pow(38.9387*Mod[8],2))+4.55531*Mod[9];
			omega0 = sqrt(pow(1.47725*Mod[10]+Shift,2)+pow(k,2))+.219181*Mod[11];
			knee = 3.28564*Mod[12]*pow(k+1.,(long double)-.721321*Mod[13])+.330483*Mod[14]*(tanh((k-22.9096*Mod[15])/10.7139/Mod[16])+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .459303*Mod[0]/sqrt(pow(k,2)+pow(1.84321*Mod[1],2))+.0386564*Mod[2];
			a = 1.79149*Mod[3]/(pow(k,2)+pow(.764836*Mod[4],2))+2.66209*Mod[5];
			b = -1856.16*Mod[6]/(pow(k-8.69519*Mod[7],2)+pow(26.3551*Mod[8],2))+3.94631*Mod[9];
			omega0 = sqrt(pow(1.45428*Mod[10]+Shift,2)+pow(k,2))+.197493*Mod[11];
			knee = 3.06296*Mod[12]*pow(k+1.,(long double)-.917081*Mod[13])+.394833*Mod[14]*(tanh((k-19.5932*Mod[15])/12.0494/Mod[16])+1);
			break;*/
		case 5://40MeV
			M_T = 1.8;
			Shift = M-M_T;
			Sigma = .00386564;
			a = 6.2;
			b = 2.8;
			omega0 = sqrt(pow(1.53+Shift,2)+pow(k,2));
			knee = .56;
			break;
	}

	long double ImSigma;	//Calculation of the argument to the exponential, these first 2 are approximations to hopefully avoid catastrophic loss of precision
	if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee < -4.)
		ImSigma = a*(omega-omega0+knee/sqrt(a*b));
	else if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee > 4.)
		ImSigma = b*(omega0-omega+knee/sqrt(a*b));
	else
		ImSigma = -.5*((a-b)*omega0-((a+b)*knee)/sqrt(a*b))+(a-b)*omega/2-sqrt(pow(((a+b)/2.)*(omega-omega0+((a-b)*knee)/(sqrt(a*b)*(a+b))),2)+pow(knee,2));

	answer += -2.*M*Sigma*exp(ImSigma);

	return(answer);
}

void ReSelf_Energy(long double M, long double omega[], long double k[], int Temp, long double Results[])	//Single quark self energy
{
	static long double Sigma[2];		//Strength
	static long double x0[2], x1[2];	//Centrality markers
	static long double gamma[2];		//Width
	static long double Shift, M_T;
	static long double k_old[2];		//Note on validity of k

	/*static bool Been_Here = false;
	static char File_Name[70] = "Self_E_Modifier.";
	static long double Mod[27];
	if(!Been_Here)
	{
		char TempChar[1] = {Temp+48};
		strncat(File_Name,TempChar,1);
		strncat(File_Name,".",1);
		strcat(File_Name,Process);
		static ifstream File(File_Name);
		for(int i = 0; i < 27; i++)
			File >> Mod[i];
		Been_Here = true;
	}*/

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
				Sigma[0] = .257498/sqrt(pow(k[0],2)+pow(1.33201,2))+.00762638;
				Sigma[1] = .257498/sqrt(pow(k[1],2)+pow(1.33201,2))+.00762638;
				x0[0] = sqrt(pow(k[0],2)+pow(1.54778+Shift,2))+.276509;
				x0[1] = sqrt(pow(k[1],2)+pow(1.54778+Shift,2))+.276509;
				x1[0] = sqrt(pow(k[0],2)+pow(1.49799+Shift,2))+.246719;
				x1[1] = sqrt(pow(k[1],2)+pow(1.49799+Shift,2))+.246719;
				gamma[0] = .658734/sqrt(pow(k[0],2)+pow(3.35217,2))+.0815109;
				gamma[1] = .658734/sqrt(pow(k[1],2)+pow(3.35217,2))+.0815109;
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571/sqrt(pow(k[0],2)+pow(1.17821,2))+.00762638;
				Sigma[1] = .212571/sqrt(pow(k[1],2)+pow(1.17821,2))+.00762638;
				x0[0] = sqrt(pow(k[0],2)+pow(1.57536+Shift,2))+.259147;
				x0[1] = sqrt(pow(k[1],2)+pow(1.57536+Shift,2))+.259147;
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194+Shift,2))+.222526;
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194+Shift,2))+.222526;
				gamma[0] = .336699/sqrt(pow(k[0],2)+pow(1.87956,2))+.0651449;
				gamma[1] = .336699/sqrt(pow(k[1],2)+pow(1.87956,2))+.0651449;
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972/sqrt(pow(k[0],2)+pow(1.41483,2))+.0101423;
				Sigma[1] = .307972/sqrt(pow(k[1],2)+pow(1.41483,2))+.0101423;
				x0[0] = sqrt(pow(k[0],2)+pow(1.56476+Shift,2))+.251031;
				x0[1] = sqrt(pow(k[1],2)+pow(1.56476+Shift,2))+.251031;
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194+Shift,2))+.222526;
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194+Shift,2))+.222526;
				gamma[0] = .550628/sqrt(pow(k[0],2)+pow(2.43968,2))+.0981269;
				gamma[1] = .550628/sqrt(pow(k[1],2)+pow(2.43968,2))+.0981269;
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131/sqrt(pow(k[0],2)+pow(1.43308,2))+.0125796;
				Sigma[1] = .339131/sqrt(pow(k[1],2)+pow(1.43308,2))+.0125796;
				x0[0] = sqrt(pow(k[0],2)+pow(1.55034+Shift,2))+.257788;
				x0[1] = sqrt(pow(k[1],2)+pow(1.55034+Shift,2))+.257788;
				x1[0] = sqrt(pow(k[0],2)+pow(1.46999+Shift,2))+.231821;
				x1[1] = sqrt(pow(k[1],2)+pow(1.46999+Shift,2))+.231821;
				gamma[0] = .615278/sqrt(pow(k[0],2)+pow(2.22298,2))+.143376;
				gamma[1] = .615278/sqrt(pow(k[1],2)+pow(2.22298,2))+.143376;
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841/sqrt(pow(k[0],2)+pow(1.42911,2))+.0157245;
				Sigma[1] = .304841/sqrt(pow(k[1],2)+pow(1.42911,2))+.0157245;
				x0[0] = sqrt(pow(k[0],2)+pow(1.55511+Shift,2))+.231105;
				x0[1] = sqrt(pow(k[1],2)+pow(1.55511+Shift,2))+.231105;
				x1[0] = sqrt(pow(k[0],2)+pow(1.44714+Shift,2))+.20956;
				x1[1] = sqrt(pow(k[1],2)+pow(1.44714+Shift,2))+.20956;
				gamma[0] = .862629/sqrt(pow(k[0],2)+pow(2.67193,2))+.189598;
				gamma[1] = .862629/sqrt(pow(k[1],2)+pow(2.67193,2))+.189598;
				break;
			/*case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571*Mod[17]/sqrt(pow(k[0],2)+pow(1.17821*Mod[18],2))+.00762638*Mod[19];
				Sigma[1] = .212571*Mod[17]/sqrt(pow(k[1],2)+pow(1.17821*Mod[18],2))+.00762638*Mod[19];
				x0[0] = sqrt(pow(k[0],2)+pow(1.57536*Mod[20]+Shift,2))+.259147*Mod[21];
				x0[1] = sqrt(pow(k[1],2)+pow(1.57536*Mod[20]+Shift,2))+.259147*Mod[21];
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194*Mod[22]+Shift,2))+.222526*Mod[23];
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194*Mod[22]+Shift,2))+.222526*Mod[23];
				gamma[0] = .336699*Mod[24]/sqrt(pow(k[0],2)+pow(1.87956*Mod[25],2))+.0651449*Mod[26];
				gamma[1] = .336699*Mod[24]/sqrt(pow(k[1],2)+pow(1.87956*Mod[25],2))+.0651449*Mod[26];
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972*Mod[17]/sqrt(pow(k[0],2)+pow(1.41483*Mod[18],2))+.0101423*Mod[19];
				Sigma[1] = .307972*Mod[17]/sqrt(pow(k[1],2)+pow(1.41483*Mod[18],2))+.0101423*Mod[19];
				x0[0] = sqrt(pow(k[0],2)+pow(1.56476*Mod[20]+Shift,2))+.251031*Mod[21];
				x0[1] = sqrt(pow(k[1],2)+pow(1.56476*Mod[20]+Shift,2))+.251031*Mod[21];
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194*Mod[22]+Shift,2))+.222526*Mod[23];
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194*Mod[22]+Shift,2))+.222526*Mod[23];
				gamma[0] = .550628*Mod[24]/sqrt(pow(k[0],2)+pow(2.43968*Mod[25],2))+.0981269*Mod[26];
				gamma[1] = .550628*Mod[24]/sqrt(pow(k[1],2)+pow(2.43968*Mod[25],2))+.0981269*Mod[26];
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131*Mod[17]/sqrt(pow(k[0],2)+pow(1.43308*Mod[18],2))+.0125796*Mod[19];
				Sigma[1] = .339131*Mod[17]/sqrt(pow(k[1],2)+pow(1.43308*Mod[18],2))+.0125796*Mod[19];
				x0[0] = sqrt(pow(k[0],2)+pow(1.55034*Mod[20]+Shift,2))+.257788*Mod[21];
				x0[1] = sqrt(pow(k[1],2)+pow(1.55034*Mod[20]+Shift,2))+.257788*Mod[21];
				x1[0] = sqrt(pow(k[0],2)+pow(1.46999*Mod[22]+Shift,2))+.231821*Mod[23];
				x1[1] = sqrt(pow(k[1],2)+pow(1.46999*Mod[22]+Shift,2))+.231821*Mod[23];
				gamma[0] = .615278*Mod[24]/sqrt(pow(k[0],2)+pow(2.22298*Mod[25],2))+.143376*Mod[26];
				gamma[1] = .615278*Mod[24]/sqrt(pow(k[1],2)+pow(2.22298*Mod[25],2))+.143376*Mod[26];
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841*Mod[17]/sqrt(pow(k[0],2)+pow(1.42911*Mod[18],2))+.0157245*Mod[19];
				Sigma[1] = .304841*Mod[17]/sqrt(pow(k[1],2)+pow(1.42911*Mod[18],2))+.0157245*Mod[19];
				x0[0] = sqrt(pow(k[0],2)+pow(1.55511*Mod[20]+Shift,2))+.231105*Mod[21];
				x0[1] = sqrt(pow(k[1],2)+pow(1.55511*Mod[20]+Shift,2))+.231105*Mod[21];
				x1[0] = sqrt(pow(k[0],2)+pow(1.44714*Mod[22]+Shift,2))+.20956*Mod[23];
				x1[1] = sqrt(pow(k[1],2)+pow(1.44714*Mod[22]+Shift,2))+.20956*Mod[23];
				gamma[0] = .862629*Mod[24]/sqrt(pow(k[0],2)+pow(2.67193*Mod[25],2))+.189598*Mod[26];
				gamma[1] = .862629*Mod[24]/sqrt(pow(k[1],2)+pow(2.67193*Mod[25],2))+.189598*Mod[26];
				break;*/
		}
	}

	Results[0] = Sigma[0]*(omega[0]-x0[0])/(pow(omega[0]-x1[0],2)+gamma[0]);
	Results[1] = Sigma[1]*(omega[1]-x0[1])/(pow(omega[1]-x1[1],2)+gamma[1]);
	return;
}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	if(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta) < 0)
		return(0.);
	else
		return(sqrt(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta)));
}

long double Set_Temp(int T)
{
	const long double Temps[] = {0,.194,.258,.32,.4,.04,.04};
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

long double Potential1(long double Par[], long double k0, long double k)	//Potiential for the numerator of the boson spectrum
{
#if VERSION == Exp
	return(exp(-abs(-4.*pow(k0,2)+4.*pow(k,2))/pow(Par[1],2)));
#elif VERSION == 22
	return(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))));
#elif VERSION == 42
	return(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)));
#elif VERSION == 24
	return(pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#endif
}

long double Potential2(long double Par[], long double k0, long double k)	//Potiential for the denominator of the T-Matrix and boson spectrum
{
#if VERSION == Exp
	return(Par[0]*exp(-2.*abs(-4.*pow(k0,2)+4.*pow(k,2))/pow(Par[1],2)));
#elif VERSION == 22
	return(Par[0]*pow(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#elif VERSION == 42
	return(Par[0]*pow(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)),2));
#elif VERSION == 24
	return(Par[0]*pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),4));
#endif
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

long double ImFolding_Integrand(long double Par[], long double k0, long double k, long double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	static long double q[2] = {Energy(0, Par[3]/2., k, theta),Energy(0, Par[3]/2., -k, theta)};
	static long double k_old = k;
	long double omega[2] = {sqrt(Par[4]+pow(Par[3],2))/2.+k0,sqrt(Par[4]+pow(Par[3],2))/2.-k0};
	long double fermi[2] = {Fermi(omega[0], Temp),Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];

	if(k_old != k)
	{
		k_old = k;
		q[0] = Energy(0, Par[3]/2., k, theta);
		q[1] = Energy(0, Par[3]/2., -k, theta);
	}

	ImSelf_Energy(Par[2], omega, q, Par, Temp, ImSelf);
	//Self_Energy(Par[2], omega, q, Par, Temp, ImSelf, ReSelf);
	ReSelf_Energy(Par[2], omega, q, Temp, ReSelf);

	return(-((4.*ImSelf[0]*ImSelf[1]*pow(Par[2],2)*(1.-fermi[0]-fermi[1]))/((pow(pow(omega[0],2)-pow(q[0],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[0],2)+pow(ImSelf[0],2))*(pow(pow(omega[1],2)-pow(q[1],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[1],2)+pow(ImSelf[1],2)))));
}
