#include<fstream>
#include<iostream>
#include<cmath>
#include<complex>
#include<cstring>
#include<omp.h>
#include<iomanip>
#include"Spectral.h"
using namespace std;

long double Correlator(long double[6], long double[3], long double[5], long double);	//Evaluates the spatial correlator
long double Integrate1(long double[6], long double[3], long double[5], long double, long double);	//1D integral with fixed theta
long double Integrate2(long double[6], long double[3], long double[5], long double);	//2D integral to evaluate the correlator

long double Epsilon = .034;//3.40672e-4;//.000380625; //.5MeV
char* Process;

int main(int argc, char* argv[])
{
	char* File = new char[100];
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
	long double Par[6] = {atof(argv[4]), atof(argv[5]), atof(argv[6]), 0, 0};	//-127.995280691106, 1.4049344847006076, 1.8: g, Lambda, M, |vec p|, E=sqrt(s)
	long double SelfPPar[3] = {atof(argv[7]), atof(argv[8]), atof(argv[9])};	//a, Sigma1, Sigma2
	long double SelfEPar[5] = {atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atof(argv[14])};	//Sigma, gamma, E_0), atof(a, b

	#pragma omp parallel for private(z, holder, Par)
	for(int i = 290*iProcess/Total; i <= 290*(iProcess+1)/Total; i++)
	{
		z = .3+i*.02;
		holder = Correlator(Par, SelfPPar, SelfEPar, z);
		#pragma omp critical
		{
			TPlot << z << " " << holder << endl;
		}
	}//*/

	return(0);
}

long double Correlator(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double z)
{
	return(Integrate2(Par, SelfPPar, SelfEPar, z));
}

long double Integrate2(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double z)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1[24];	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x3[24];
	long double distance[] = {5e-2, 2e-2, 1.5e-2, 1e-2, 2.5e-3, 1e-4, 1e-5, 1e-6};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = 0;
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

			F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 3.040308;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
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

			F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 3.6;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
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

			F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
			F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	b = 4;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	b = 10;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	a = b;

	b = 23.5;
	F_a = F_b = 0;	//Start integration at 0
	for(j = 0; j < 24; j++)
	{
		x1[j] = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
		x3[j] = (b+a+Disp[j]*(b-a))/2.;

		F_a += Integrate1(Par, SelfPPar, SelfEPar, x1[j], z)*w[j+1];	//Evaluate k integral at x1
		F_b += Integrate1(Par, SelfPPar, SelfEPar, x3[j], z)*w[j+1];	//Evaluate k integral at x3
	}
	F_ave = Integrate1(Par, SelfPPar, SelfEPar, a/2.+b/2., z);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);

	return(Answer);
}

long double Integrate1(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double E, long double z)
{
	Par[4] = E;
	long double Value = Spectral(Par, SelfPPar, SelfEPar, E, 0);
	if(abs(z-.3) <= .05)
		cout << E << " " << Value << endl;
	return(2.*M_PI*exp(-E*z)*Value);	//return the best estimate of the integral on the interval*/
}
