#include<fstream>
#include<iostream>
#include<cmath>
#include<complex>
#include<cstring>
#include<omp.h>
#include<iomanip>
#include"Spectral.h"
using namespace std;

long double Spectral(long double***, long double, long double, long double);	//The tabulated points and 2 inputs and returns the bicubic interpolation of that input
long double Correlator(long double***, long double, long double);	//Evaluates the spatial correlator
long double Integrate1(long double***, long double, long double);	//1D integral with fixed theta
long double Integrate2(long double***, long double);	//2D integral to evaluate the correlator

char* Process;

int main(int argc, char* argv[])
{
	char* File = new char[25];
	strcpy(File, argv[3]);	//Name of the file
	Process = argv[1];
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot(File);
	TPlot << setprecision(18);
	long double z;		//The position value of the spactial correlator
	long double R_max = 80;	//Largest sqrt(omega^2+p^2)=sqrt(s+2p^2) that will be integrated
	long double holder;
	const int iProcess = atoi(argv[1]);
	const int Total = atoi(argv[2]);

	#pragma omp parallel for private(z, holder)
	for(int i = 290*iProcess/Total; i <= 290*(iProcess+1)/Total; i++)
	{
		z = .3+i*.02;
		holder = Correlator(Table, z, R_max);
		#pragma omp critical
		{
			TPlot << z << " " << holder << endl;
		}
	}//*/

	return(0);
}

long double Correlator(long double*** Table, long double z, long double R_max)
{
	return(Integrate2(Table, z));
}

long double Integrate2(long double*** Table, long double z)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1[24];	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x3[24];
	long double distance[] = {5e-2, 2e-2, 1.5e-2, 1e-2, 2.5e-3, 1e-4, 1e-5, 1e-6};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = 1;
	long double b = 0;
	int i, j;

	for(i = 0; i < 8; i++)
	{
		b = 3.040308 - distance[i];
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Table, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 3.040308;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Table, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	for(i = 7; i >= 0; i--)
	{
		b = 3.040308 + distance[i];
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Table, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 3.6;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Table, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	for(i = 7; i >= 0; i--)
	{
		b = 3.6 + distance[i];
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Table, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 4;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Table, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	b = 10;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Table, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	b = 23.5;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Table, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Table, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Table, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);

	return(Answer);
}

long double Integrate1(long double*** Table, long double E, long double z)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1[24];	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x3[24];
	long double Answer = 0;
	long double stride = 2*M_PI/z;	//Stride of the integral
	long double F_a, F_b, F_ave;
	long double a = 0;
	long double b = 0;
	long double P_Max = 600;
	int i;

//This code is for integating out to a z dependant boundary and may be used.
	for(int j = 0; j < 1000; j++)	//need to start the count off from where it left off in the previous integration block
	{
		b += stride;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Spectral(Table, E, x1[i], z)*4.*E/(x1[i]*x1[i]+E*E)*cos(z*x1[i])*w[i+1];	//Evaluate k integral at x1
			F_b += Spectral(Table, E, x3[i], z)*4.*E/(x3[i]*x3[i]+E*E)*cos(z*x3[i])*w[i+1];	//Evaluate k integral at x3
		}
		F_ave = Spectral(Table, E, a/2.+b/2., z)*4.*E/(pow(a/2.+b/2.,2)+E*E)*cos(z*(a/2.+b/2.));
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}	//For the bulk of the integral where either the result is well approximated by either the finite or zero width analytic result*/

	b += stride/2.;	//Evaluate an extra half stride, should result in a result between the extrems of a quarter and three quarter strides
	F_a = F_b = 0;	//Start integration at 0
	for(i = 0; i < 24; i++)
	{
		x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
		x3[i] = (b+a+Disp[i]*(b-a))/2.;

		F_a += Spectral(Table, E, x1[i], z)*4.*E/(x1[i]*x1[i]+E*E)*cos(z*x1[i])*w[i+1];	//Evaluate k integral at x1
		F_b += Spectral(Table, E, x3[i], z)*4.*E/(x3[i]*x3[i]+E*E)*cos(z*x3[i])*w[i+1];	//Evaluate k integral at x3
	}
	F_ave = Spectral(Table, E, a/2.+b/2., z)*4.*E/(pow(a/2.+b/2.,2)+E*E)*cos(z*(a/2.+b/2.));
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);//*/

	return(Answer);	//return the best estimate of the integral on the interval*/
}

long double Spectral(long double*** Table, long double E, long double p, long double z)
{
	Gamma = .032+.0023855870336295*(p*p+pow(7.3820000424184,2))*exp(-pow(p/11.486249371586,2));
	if(E > 3.6)
		return(arctanh(E)*1.6/M_PI*Gamma/(pow(E-3.040308,2)+Gamma*Gamma)+14.58/(M_PI*M_PI)*sqrt(1-12.96/(E*E)));
	else
		return(arctanh(E)*1.6/M_PI*Gamma/(pow(E-3.040308,2)+Gamma*Gamma));
}
