//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
#include<cfloat>
using namespace std;

long double Energy(long double, long double, long double, long double);	//Returns sqrt(M^2+p^2+k^2-pk cos(theta))
long double ReInt(long double[6], long double, long double, int);	//Returns the real part of the integrad from 2p_0 to infinity
long double ImInt(long double[6], long double, long double, int);	//Returns the imaginary part of the integrad from 2p_0 to infinity
long double Integrate1(long double(*)(long double[6], long double, long double, int), long double[6], long double, int);	//Integrates the part that it is told to integrate. It uses the difference in evaluating the trapaziod rule and Simpson's rule to pick the points that it integrate between. Since the Gaussian quadrature is very accuate for using very few points (2n-1 order polynomial accurate with n points), it then returns the Gaussian quadrature for 3 points on that subinterval.
long double Integrate2(long double(*)(long double[6], long double, long double, int), long double[6], int);	//Contains more brains than Integrate1() as it will need to divide the integral into 2 parts and pass the endpoints down for faster times but it uses the same algorithm to acheive its results.
long double Self_Energy(long double, long double, long double, int);	//Returns the Self-Energy
long double Self_P_Depends(int, long double);	//Returns the momentum dependance of the self-energy
long double Self_E_Depends(int, long double, long double, long double);	//Returns the energy dependance of the self-energy that moves with the momentum
long double ReProp(long double[6], long double, long double, int);	//Returns the real part of the propagator
long double ImProp(long double[6], long double, long double, int);	//Returns the imaginary part of the propagator
void Characterize(long double[6], long double, long double, int, long double*&, long double*&, int&);	//Characterizes peaks
bool Minimize(long double[6], long double, long double, int, long double&, long double&);	//Minimizes by binary search
long double PropIntegrand(long double, long double[6], long double, long double, int);	//The integrand for the imaginary part of the propagator
long double Rho(long double, long double[6], long double, long double, int);	//Single particle propagator
long double LawCosines(long double, long double, long double);	//Returns the law of cosines for two vectors with an angle inbetween.
long double Potential(long double[6], long double, long double, int);	//Returns the potential CC*Lambda^2/(M*(Lambda^2-4k^mu k_mu))
complex<long double> TMatrix(long double[6], int);	//Returns the T-matrix for a given M, P, E=sqrt(s), and T
long double Spectral(long double[6], int);	//Returns the spectral function of the T-matrix
long double G_0Int(long double[6], long double, long double, int);	//Returns the integrand for G_0. This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
long double ReDelta_GInt(long double[6], long double, long double, int);	//Returns the real part of the Delta G integrand
long double ImDelta_GInt(long double[6], long double, long double, int);	//Returns the imaginary part of the Delta G integrand
long double Potential1(long double[6], long double, long double, int);	//Returns one of the factors of the potiential that is in Potential() without the coupling constant.
long double Fermi(long double, int);	//Returns the fermi factor
long double Fermi(long double[6], long double, long double, int);
long double Analytic(long double, long double);	//The analytic spectral function for vacuum, which is constant in P
long double ApproxSpectral(long double, long double);
complex<long double> arctan(complex<long double>);
complex<long double> arctanh(complex<long double>);

long double ApproxSpectral(long double E, long double P)
{
	long double Epsilon;
	static long double Table[375][2];
	static bool Run = false;
	ifstream List("./Epsilon List");
	int i;

	if(!Run)
	{
		for(i = 0; i < 375; i++)
		{
			List >> Table[i][0];	//Grab the first column and overwrite the previous fourth column
			List >> Table[i][1] >> Table[i][1];	//Grab the second and then over write with the third column
			if(i != 401)
				List >> Table[i+1][0];	//Grab the fourth column except for the last one as it will go into oblivion which may cause segfault
		}
		Run = true;	//prevent a run of this code
	}

	if(P <= 256.8)	//Last reasonable P for measured epsilon
	{
		i = P/.8;
		Epsilon = Table[i][1]*(1.+i-1.25*P)+Table[i+1][1]*(1.25*P-i)*2.;	//1.25=.8^-1
	}
	else
		Epsilon = .032;

	if(Epsilon < 2.81e-6)	//Anything below this is unreliable and should thus be counted as zero width to decent approximation. This is on the edge of resolution of the numerical integrator (Gamma_FWHM<800eV, smallest interval=1keV), While this is true for analytic case, it fails sooner for in-media.
		Epsilon = .032;	

	return(Analytic(E, Epsilon));
}

complex<long double> arctan(complex<long double>x)
{
	complex<long double> i(0,1);
	complex<long double> one(1,0);
	return(complex<long double>(0,.5)*log((one-i*x)/(one+i*x)));
}

complex<long double> arctanh(complex<long double>x)
{
	complex<long double> one(1,0);
	return(complex<long double>(.5,0)*log((one+x)/(one-x)));
}

long double Analytic(long double E, long double Epsilon)	//This strictly vacuum, the width has no depandance beside the E, ever. The epsilon is .008GeV for the correct natural width.
{//Stupid complex number object doesn't include binary functions of reals and complex numbers
	long double M = 1.8;
	long double CC = -127.995280691106;
	long double Lambda = 1.4049344847006076;
	Epsilon = Epsilon*pow((E*E-1.258884)/7.984588734864,2.5)*pow(1.618884/(.36+E*E),2);

	complex<long double> SMEpsilon(E*E/4.-M*M,Epsilon*E);	//The analytic needs to divide the energy by 2 to match what the integration is doing
	complex<long double> SEpsilon(E*E/4.,Epsilon*E);
	complex<long double> MLambda1(pow(Lambda/2.,2)-M*M,0);
	complex<long double> MLambda2(pow(M/Lambda,2)-.25,0);
	complex<long double> SMLambdaEpsilon((E*E+pow(Lambda,2))/4.-M*M,Epsilon*E);
	complex<long double> Num, Den, Non;

	Num = complex<long double>(pow(Lambda*M/(4.*M_PI),2),0)*(complex<long double>(2,0)*sqrt(MLambda1*SMEpsilon)*arctan(sqrt(-SEpsilon/SMEpsilon))-complex<long double>(Lambda,0)*sqrt(SEpsilon)*arctan(complex<long double>(2,0)*sqrt(MLambda2)))/(sqrt(SEpsilon*-MLambda1)*SMLambdaEpsilon);
	Den = complex<long double>(CC*pow(Lambda*Lambda*M/(8*M_PI),2),0)/pow(SMLambdaEpsilon,2)*(SMLambdaEpsilon/MLambda1+(complex<long double>(pow(Lambda,4)/4.-pow(Lambda*M,2)/2,0)+complex<long double>(2*M*M,0)*SMEpsilon)/(complex<long double>(Lambda,0)*pow(-MLambda1,(long double)(1.5)))*arctan(complex<long double>(2,0)*sqrt(MLambda2))+complex<long double>(2,0)*sqrt(-SMEpsilon/SEpsilon)*arctan(sqrt(-SEpsilon/SMEpsilon)));
	Non = complex<long double>(M*M/(2.*M_PI*M_PI),0)*sqrt(SMEpsilon/SEpsilon)*arctanh(sqrt(SEpsilon/SMEpsilon));

	//cout << E << " " << Num.real() << " " << Num.imag() << " " << Den.real() << " " << Den.imag() << " " << -18./M_PI*Non.imag() << endl;

	return(-18./M_PI*(Non+CC*(pow(Num,(long double)(2.))/(complex<long double>(1.,0)-Den))).imag());
}

long double Self_Energy(long double E, long double P, long double M, int Temp)
{
	long double Ans = Self_E_Depends(Temp, E, P, M)*Self_P_Depends(Temp, P);
	return(Ans);
}

long double Self_P_Depends(int Temp, long double P)
{
	long double Par[6];

	switch(Temp)
	{
		case 0:
			Par[0] = 0;
			Par[1] = 1;
			Par[2] = 1;
			Par[3] = 0;
			Par[4] = 1;
			Par[5] = 1;
			break;
		case 1:
			Par[0] = .7359389831810698;
			Par[1] = 7.487501146014314;
			Par[2] = 1.9490238595657456;
			Par[3] = .1;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 2:
			Par[0] = .7409390219065235;
			Par[1] = 7.450458343071824;
			Par[2] = 1.8620618988580635;
			Par[3] = .1;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 3:
			Par[0] = .7426375963204489;
			Par[1] = 7.698646415632565;
			Par[2] = 1.771465704769189;
			Par[3] = .1;
			Par[4] = 10;
			Par[5] = 3;
			break;
	}

	return(Par[0]*exp(-pow(P/Par[1],2))+(1-Par[0])*exp(-pow(P/Par[2],2))+Par[3]/(1.+exp((Par[4]-P)/Par[5])));
}

long double Self_E_Depends(int Temp, long double E, long double P, long double M)
{
	long double E_0 = Energy(M,P,0,0); //location of lorentzian
	long double Sigma; //size of energy dependance
	long double gamma; //width of lorentzian
	long double a; //Exponential parameters

	switch(Temp)
	{
		case 0:
			Sigma = 0;
			gamma = 1;
			a = 1;
			break;
		case 1://235.2MeV
			Sigma = .14647773684797566;
			gamma = .34024341413053605;
			a = .9105749017583176; //1/a=1.098GeV=4.67T 
			break;
		case 2://294MeV
			Sigma = .13311625968464938;
			gamma = .4114357364423404;
			a = 1.0017017464500062; //1/a=.9983GeV=3.40T
			break;
		case 3://362MeV
			Sigma = .20683830569398817;
			gamma = .5994753305596638;
			a = 1.2684796687619544; //1/a=.788GeV=2.18T
			break;
	}
	if(E > E_0)
		a = 0;

	return(exp(a*(E-E_0))*Sigma*gamma/M_PI*(1/(pow(E+E_0,2)+pow(gamma,2))-1/(pow(E-E_0,2)+pow(gamma,2))));
}

long double Spectral(long double Par[6], int Temp)
{
	long double G_0;	//The imaginary part of G_0
	complex<long double> TMat = TMatrix(Par, Temp);
	complex<long double> Num;	//The integral in the numerator of Delta G
	int N_f = 3;
	int N_c = 3;

	G_0 = Integrate2(G_0Int, Par, Temp);
	Num = complex<long double>(Integrate2(ReDelta_GInt, Par, Temp),0);
	Num += complex<long double>(0,Integrate2(ImDelta_GInt, Par, Temp));
	//cerr << Par[4] << " " << Num.real() << " " << Num.imag() << " " << TMat.real() << " " << TMat.imag() << " " << G_0 << endl;

	return(-2.*N_f*N_c/M_PI*(G_0+(Par[0]*pow(Num,2)*TMat).imag()));
}

long double Fermi(long double E, int T)
{
	long double Temp; //T_c = .196GeV = 196MeV

	switch(T)
	{
		case 0:
			return(0);
			break;
		case 1:
			Temp = .196*1.2;
			break;
		case 2:
			Temp = .196*1.5;
			break;
		case 3:
			Temp = .196*2.;
			break;
	}

	return(1./(exp(E/Temp)+1.)); //Fermi factor
}

long double Fermi(long double Par[6], long double k, long double theta, int T)
{
	long double Temp; //T_c = .196GeV = 196MeV

	switch(T)
	{
		case 0:
			return(0);
			break;
		case 1:
			Temp = .196*1.2;
			break;
		case 2:
			Temp = .196*1.5;
			break;
		case 3:
			Temp = .196*2.;
			break;
	}

	return(1./(exp(Energy(Par[2], Par[3]/2., k, theta)/Temp)+1.)); //Fermi factor
}

long double G_0Int(long double Par[6], long double k, long double theta, int Temp)	//This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
{
	return(ImProp(Par, k, theta, Temp)*sin(theta)*k*k);
}

long double ReDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ReProp(Par, k, theta, Temp)*Potential1(Par, k, theta, Temp));
}

long double ImDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ImProp(Par, k, theta, Temp)*Potential1(Par, k, theta, Temp));
}

long double Potential1(long double Par[6], long double k, long double theta, int Temp)
{
	switch(Temp)
	{
		case 1:
			Par[1] *= exp(-1./60.);
			break;
		case 2:
			Par[1] *= exp(-1./30.);
			break;
		case 3:
			Par[1] *= exp(-.05);
			break;
	}
	return(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.));
}

complex<long double> TMatrix(long double Parameters[6], int Temp)
{
	complex<long double> Int_Holder;	//Holder for the result of the integration, allows it to be calculated once
	Int_Holder = complex<long double>(Integrate2(ReInt, Parameters, Temp), 0);
	Int_Holder += complex<long double>(0, Integrate2(ImInt, Parameters, Temp));
	Int_Holder = complex<long double>(1.,0.)/(complex<long double>(1.,0.)-Int_Holder);	//Integrate once, where it says Parameters[6] in the numerator, I need to put V(p,p') in the event that I didn't get that correct

	return(Int_Holder);
}

//Parameters = g, Lambda, M, |vec p|, E=sqrt(s)
long double ReInt(long double Par[6], long double k, long double theta, int Temp)	//Returns the real part of the integrand
{
	return(ReProp(Par, k, theta, Temp)*Potential(Par, k, theta, Temp)*sin(theta)*k*k);
}

long double ImInt(long double Par[6], long double k, long double theta, int Temp)	//Returns the imaginary part of the integrand
{
	return(ImProp(Par, k, theta, Temp)*Potential(Par, k, theta, Temp)*sin(theta)*k*k);
}

long double ReProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the real part of the propagator
{
	//long double Quasi = (pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp), 2))*2.*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp)), 2));
	//return(Quasi);

	long double Disp[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.}; //Displacement from center for 9th order Gauss-Legendre integration
	long double w[] = {128./225., (322.+13.*sqrt(70))/900., (322.-13.*sqrt(70))/900.}; //Weight of the function at Disp
	long double Range16[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double Range8[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	long double Range4[] = {-1,-.5,0,.5,1};	//Number of gamma from center
	long double LocalPar[] = {Par[0], Par[1], Par[2], Par[3], Par[4], Par[5]};	//The local copy of Par to be sent to ImProp
	long double zero = sqrt(2.*(k*k+Par[2]*Par[2]-Par[3]*Par[3]/4.+Energy(Par[2],Par[3]/2.,k,theta)*Energy(Par[2],Par[3]/2.,-k,theta)));	//2 particle on-shell
	long double gamma = -2.*Self_Energy(zero, LawCosines(Par[3]/2., k, theta), Par[2], Temp);	//These are the widths of the features near 2 Particle on shell
	long double Answer = 0;
	long double a = 0;
	long double b;
	long double F_a, F_b, F_ave;
	long double x1[9], x3[9];
	long double Width;	//Step size for integration
	long double E = zero+64.*gamma+3.;	//Largest feature I can find
	long double f0 = ImProp(Par, k, theta, Temp);	//Par[4]^2 is the location of the division by zero
	int i,j,l = 0;
	int version;

	if(k <= 4 && abs(zero-Par[4]) <= 1)
		version = 16;
	else if((k > 4 && abs(zero-Par[4]) <= 1) || k <= 4)
		version = 8;
	else
		version = 4;

	if(pow(Par[4],2) == 0)
		f0 = 0;

	if(zero-64*gamma-3. < 0)
		a = b = 0;
	else
		a = b = zero-64*gamma-3.;
	i = 0;
	do
	{
		Width = 3.;	//No-man's land

		if(a<zero-64.*gamma && b+Width>=zero-64.*gamma)	//Stutter step before the peak
		{
			Width = zero-64.*gamma-b;
			l = 0;	//Resets l before entering the peak
		}
		else if((a>=zero-64.*gamma && b<=zero+64.*gamma) && l < version)	//Integrating the peak itself
		{
			switch(version)
			{
				case 4:
					Width = gamma*(Range4[l+1]-Range4[l]);
					break;
				case 8:
					Width = gamma*(Range8[l+1]-Range8[l]);
					break;
				case 16:
					Width = gamma*(Range16[l+1]-Range16[l]);
					break;
			}
			l++;	//Leaving the peak prevents illegal space access at Range[l+1]
		}

		b += Width;

		F_a = F_b = 0;
		for(j = 0; j < 2; j++)
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.; //Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			LocalPar[4] = x1[j];
			F_a += 2.*x1[j]*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow(x1[j],2)-pow(Par[4],2))*w[j+1]; //Evaluate function at x1
			LocalPar[4] = x3[j];
			F_b += 2.*x3[j]*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow(x3[j],2)-pow(Par[4],2))*w[j+1]; //Evaluate function at x3
		}
		LocalPar[4] = (a+b)/2.;
		F_ave = (a+b)*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow((a+b)/2.,2)-pow(Par[4],2))*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}while(b < E);

	if(f0 == 0)
		Answer = Answer/M_PI;
	else
		Answer = Answer+f0*log(abs(1.-pow(b/Par[4],2)))/M_PI;

	//cout << Par[4] << " " << Par[3] << " " << k << " " << Answer << " " << Quasi << endl;
	return(Answer);//*/
}

long double ImProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the imaginary part of the propagator
{
	/*long double Quasi = 2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp))*2.*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,-k,theta), LawCosines(Par[3]/2., k, theta), Par[2], Temp)+Self_Energy(sqrt(pow(Par[4],2)+pow(Par[3],2))-Energy(Par[2],Par[3]/2.,k,theta), LawCosines(Par[3]/2., -k, theta), Par[2], Temp)), 2));
	return(Quasi);
	long double PNarrow = pow(Par[2],2)*Rho(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta), Par, -k, theta, Temp)*(1.-Fermi(Energy(Par[2], Par[3]/2., k, theta), Temp)-Fermi(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta), Temp))/(pow(2.*M_PI,2)*Energy(Par[2], Par[3]/2., k, theta));
	long double NNarrow = pow(Par[2],2)*Rho(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., -k, theta), Par, k, theta, Temp)*(1.-Fermi(Energy(Par[2], Par[3]/2., -k, theta), Temp)-Fermi(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., -k, theta), Temp))/(pow(2.*M_PI,2)*Energy(Par[2], Par[3]/2., -k, theta));

	if(PNarrow<0 && NNarrow<0)
		return(Quasi);
	else if(PNarrow<0 || NNarrow<0)	//Effectively exclusive or
	{
		if(PNarrow < 0)
			return(PNarrow);
		else
			return(NNarrow);
	}*/

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	//long double Disp[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.}; //Displacement from center for 9th order Gauss-Legendre integration
	//long double w[] = {128./225., (322.+13.*sqrt(70))/900., (322.-13.*sqrt(70))/900.}; //Weight of the function at Disp
	//long double Range[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	//long double Range[] = {-1,-.5,0,.5,1};	//Number of gamma from center
	long double* zero;	//These are the points that may cause the greatest problems, but only if they are between 0 and Par[4]=E=sqrt s
	long double* gamma;	//These are the widths of the features near zero1 and zero2
	long double Answer = 0;
	long double a = 0;
	long double b;
	long double F_a, F_b, F_ave;
	long double x1[9], x3[9];
	long double Width;	//Step size for integration
	int i1,i2;	//Peak counters
	int j,l;	//Point and interval counters
	int Peaks;
	long double Early = 0;	//Early change from one peak to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;	//The next width that will be used in the event of an early change of peaks
	//long double TotalArea = 0;
	int version = 8;

	Characterize(Par, k, theta, Temp, zero, gamma, Peaks);
	//for(i1 = 0; i1 < Peaks; i1++)
	//	TotalArea += area[i1];

	a = b = 0;
	i1 = l = 0;
	i2 = 1;
	while(zero[i1]+Range[l]*gamma[i1] < 0 && Peaks != 0)	//Moves l up until zero[i]+Range[l]*gamma[i] is greater than 0
		l++;

	if(zero[i1]+64.*gamma[i1] > zero[i2]-64.*gamma[i2] && i2 < Peaks && l != 0)
		Early = zero[i1]+(zero[i2]-zero[i1])/(gamma[i1]+gamma[i2])*gamma[i1];
	else
		Early = 0;

	do
	{
		if(b == 0 && l != 0)	//First peak is closer than 64*gamma to 0
			Width = zero[i1]+Range[l]*gamma[i1];
		else if(b < 2.*Par[2] || b >= sqrt(Par[4]*Par[4]+Par[3]*Par[3])-2.*Par[2])
			Width = 1;	//The Self-energy peaks
		else
			Width = 3;	//No-man's land

		if(l == version && i1 < Peaks)	//Last peak has been integrated and there exists a next peak
		{
			i1++;
			i2++;
			l = 0;
		}
		if(i1 == i2)
			i2++;

		if((a<zero[i1]-64.*gamma[i1] && b+Width>=zero[i1]-64.*gamma[i1]) && Peaks != 0)	//Stutter step before the peak
		{
			Width = zero[i1]-64.*gamma[i1]-b;
			l = 0;	//Resets l before entering the peak
			if(zero[i1]+64.*gamma[i1] > zero[i2]-64.*gamma[i2] && i2 < Peaks)
				Early = zero[i1]+(zero[i2]-zero[i1])/(gamma[i1]+gamma[i2])*gamma[i1];
			else
				Early = 0;
		}
		else if((a>=zero[i1]-64.*gamma[i1] && b<=zero[i1]+64.*gamma[i1]) && Peaks != 0 && l < version && b != 0)	//Integrating the peak itself
		{
			Width = gamma[i1]*(Range[l+1]-Range[l]);
			l++;	//Leaving the peak prevents illegal space access at Range[l+1]
		}

		if(NextWidth != 0 && NextWidth == NextWidth)
		{
			Width = NextWidth;
			NextWidth = 0;
		}

		b += Width;

		if(Early != 0 && b > Early)	//Code for changing peaks early
		{
			b = Early;
			Early = 0;
			l = 16-l;
			NextWidth = (b-a)*gamma[i2]/gamma[i1];
			i1 = i2;
		}

		if(b > sqrt(Par[4]*Par[4]+Par[3]*Par[3]))
			b = sqrt(Par[4]*Par[4]+Par[3]*Par[3]);

		F_a = F_b = 0;
		for(j = 0; j < 9; j++)//2
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.; //Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			F_a += PropIntegrand(x1[j], Par, k, theta, Temp)*w[j+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[j], Par, k, theta, Temp)*w[j+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par, k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}while(b < sqrt(Par[4]*Par[4]+Par[3]*Par[3]));

	delete[] zero;
	delete[] gamma;

	Answer /= pow(2.*M_PI,2);

	//cout << Par[4] << " " << Par[3] << " " << k << " " << Answer << endl;

	return(Answer);
}

void Characterize(long double Par[6], long double k, long double theta, int Temp, long double*& zero, long double*& gamma, int& Peaks)	//Searches for minimum, then gamma by binary search
{
	float DeltaE = .5;
	long double Array[3][int(Par[4]/DeltaE)+1];	//0 is function, 1 is second derivative, 2 is boundaries
	long double Center, Width;
	long double Maxima, MaximaW;	//Tempory storage
	long double Value[4];
	int i, j;
	bool Done;

	Peaks = 0;
	for(i = 0; i < Par[4]/DeltaE; i++)	//Evaluate the function on a mesh and take the second derivative
		Array[0][i] = PropIntegrand(i*DeltaE, Par, k, theta, Temp);
	Array[1][0]=(3*Array[0][0]-4*Array[0][1]+Array[0][2])/(2.*DeltaE);
	Array[1][int(Par[4]/DeltaE)]=(-Array[0][int(Par[4]/DeltaE)]+4*Array[0][int(Par[4]/DeltaE)-1]-3*Array[0][int(Par[4]/DeltaE)-2])/(2.*DeltaE);
	for(i = 1; i < Par[4]/DeltaE-1; i++)
		Array[1][i] = (Array[0][i+1]-Array[0][i-1])/(DeltaE*2.);
	for(i = 1; i < Par[4]/DeltaE; i++)	//Count the number of minima by the number of changes from positive to negative second derivatives
	{
		if(Array[1][i-1]*Array[1][i] <= 0 && Array[1][i-1] < Array[1][i])
		{
			Array[2][Peaks] = i*DeltaE;	//List the regions were there exist at least 1 peak
			Peaks++;
		}
	}/*/
	Array[2][0] = Energy(Par[2], Par[3]/2., k, theta);
	Array[2][1] = sqrt(Par[3]*Par[3]+Par[4]*Par[4])-Energy(Par[2], Par[3]/2., -k, theta);
	Peaks = 2;*/

	j = 0; //Number of peaks located
	for(i = 0; i < Peaks; i++)
	{
		Center = Array[2][i];
		Width = DeltaE;

		Done = Minimize(Par, k, theta, Temp, Center, Width);
		if(Done)
		{
			Array[0][j] = Center;
			j++;
		}
		else	//Assume 2 peaks in region and no more
		{
			Maxima = Center;
			MaximaW = Width;

			Center = Maxima-MaximaW/2.;
			Width = MaximaW/2.;
			if(!Minimize(Par, k, theta, Temp, Center, Width))
				cerr << "#Minimization failure, E = " << Par[4] << " P = " << Par[3] << " k = " << k << " theta = " << theta << " Center = " << Center << " Width = " << Width << endl;
			Array[0][j] = Center;
			j++;

			Center = Maxima+MaximaW/2.;
			Width = MaximaW/2.;
			if(!Minimize(Par, k, theta, Temp, Center, Width))
				cerr << "#Minimization failure, E = " << Par[4] << " P = " << Par[3] << " k = " << k << " theta = " << theta << " Center = " << Center << " Width = " << Width << endl;
			Array[0][j] = Center;
			j++;
		}
	}
	Peaks = j;

	zero = new long double[Peaks];
	gamma = new long double[Peaks];
	for(i = 0; i < Peaks; i++)
	{
		zero[i] = Array[0][i];
		Maxima = Array[0][i]*sqrt(LDBL_EPSILON);
		Value[0] = PropIntegrand(Array[0][i]-Maxima, Par, k, theta, Temp);
		Value[1] = PropIntegrand(Array[0][i], Par, k, theta, Temp);
		Value[2] = PropIntegrand(Array[0][i]+Maxima, Par, k, theta, Temp);
		gamma[i] = -2.*Value[1]*pow(Maxima,2)/(Value[0]-2.*Value[1]+Value[2]);

		if(gamma[i] < 0 || gamma[i] > 10)
		{
			Value[3] = PropIntegrand(Array[0][i]+2.*Maxima, Par, k, theta, Temp);
			gamma[i] = -pow(Maxima,2)+pow(2.*Maxima,2)*Value[2]/(Value[3]-Value[1]);
		}

		if(gamma[i] < 0 || gamma[i] > 10)
		{
			Width = .001;
			Minimize(Par, k, theta, Temp, zero[i], Width);	//Not quite minimum, very near by and needs to be retried

			Maxima = zero[i]*sqrt(LDBL_EPSILON);
			Value[0] = PropIntegrand(zero[i]-Maxima, Par, k, theta, Temp);
			Value[1] = PropIntegrand(zero[i], Par, k, theta, Temp);
			Value[2] = PropIntegrand(zero[i]+Maxima, Par, k, theta, Temp);
			gamma[i] = -2.*Value[1]*pow(Maxima,2)/(Value[0]-2.*Value[1]+Value[2]);
		}

		if(gamma[i] < 0 || gamma[i] > 10)
		{
			Value[3] = PropIntegrand(zero[i]+2.*Maxima, Par, k, theta, Temp);
			gamma[i] = -pow(Maxima,2)+pow(2.*Maxima,2)*Value[2]/(Value[3]-Value[1]);
		}

		if(gamma[i] < 0 || gamma[i] > 10)
			gamma[i] = -Self_Energy(Par[4], LawCosines(Par[3],k,theta), Par[2], Temp);
		else
			gamma[i] = sqrt(gamma[i]);
	}

	return;
}

bool Minimize(long double Par[6], long double k, long double theta, int Temp, long double& Center, long double& Width)
{	//Returns success or failure. Failure indicates multiple minima. Center is the starting point, and returns the answer or center at time of failure. Width is how far to the left and right it started and finished wheather in success or failure. Function proceeds by binary search.
	long double Pos[3] = {Center-Width, Center, Center+Width};
	long double Value[3] = {PropIntegrand(Pos[0], Par, k, theta, Temp), PropIntegrand(Pos[1], Par, k, theta, Temp), PropIntegrand(Pos[2], Par, k, theta, Temp)};
	long double TestPos[2] = {Center-Width/2., Center+Width/2.};
	long double TestValue[2] = {PropIntegrand(TestPos[0], Par, k, theta, Temp), PropIntegrand(TestPos[1], Par, k, theta, Temp)};

	do	//Binary Search for minimum
	{
		if((TestValue[0] < Value[0] && TestValue[0] < Value[1]) && (TestValue[1] < Value[1] && TestValue[1] < Value[2]))	//TestValue is greater than either value on either side indicating multiple minima in range
		{
			Center = Pos[1];
			Width = Pos[2]-Pos[1];
			if(Width < 1e-5 || abs(TestValue[0]/TestValue[1]-1) < 1e-5)
				return(true);
			return(false);
		}
		else if(TestValue[0] < TestValue[1])
		{
			Value[2] = Value[1];
			Pos[2] = Pos[1];
			Value[1] = TestValue[0];
			Pos[1] = TestPos[0];
		}
		else
		{
			Value[0] = Value[1];
			Pos[0] = Pos[1];
			Value[1] = TestValue[1];
			Pos[1] = TestPos[1];
		}
		TestPos[0] = (Pos[0]+Pos[1])/2.;
		TestPos[1] = (Pos[2]+Pos[1])/2.;
		TestValue[0] = PropIntegrand(TestPos[0], Par, k, theta, Temp);
		TestValue[1] = PropIntegrand(TestPos[1], Par, k, theta, Temp);
	}while(Pos[2]-Pos[0] > 1e-5);

	Center = Pos[1];
	Width = Pos[2]-Pos[1];
	return(true);
}

long double PropIntegrand(long double omega, long double Par[6], long double k, long double theta, int Temp)
{
	return(-Par[2]*Par[2]/M_PI*Rho(omega, Par, k, theta, Temp)*Rho(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-omega, Par, -k, theta, Temp)*(1.-Fermi(omega, Temp)-Fermi(sqrt(Par[4]*Par[4]+Par[3]*Par[3])-omega, Temp)));
}

long double Rho(long double omega, long double Par[6], long double k, long double theta, int Temp)
{
	return(Self_Energy(omega, LawCosines(Par[3]/2., k, theta), Par[2], Temp)/(Energy(Par[2], Par[3]/2., k, theta)*(pow(omega-Energy(Par[2], Par[3]/2., k, theta),2)+pow(Self_Energy(omega, LawCosines(Par[3]/2., k, theta), Par[2], Temp),2))));
}

long double Potential(long double Par[6], long double k, long double theta, int Temp)	//Returns the potential CC*(Lambda^2/(M*(Lambda^2-4k^mu k_mu)))^2
{
	switch(Temp)
	{
		case 1:
			Par[1] *= exp(-1./60.);
			break;
		case 2:
			Par[1] *= exp(-1./30.);
			break;
		case 3:
			Par[1] *= exp(-.05);
			break;
	}
	return(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.), 2));
}

long double Energy(long double M, long double P, long double k, long double theta)	//Returns twice the energy sqrt(M^2+(vec P/2+vec k)^2)
{
	return(sqrt(M*M+P*P+k*k-2.*P*k*cos(theta)));
}

long double LawCosines(long double P, long double k, long double theta)	//Returns the law of cosines for two vectors with an angle in between.
{
	return(sqrt(pow(P,2)+pow(k,2)-2.*P*k*cos(theta)));
}

long double Integrate2(long double(*Integrand)(long double[6], long double, long double, int), long double Par[6], int Temp)
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double x1[9];	//These are the two other points required for 5th order Gaussian quadrature for this interval
	long double x3[9];
	long double range[] = {acos(sqrt(-pow(Par[4]/Par[3],2)+sqrt(1.+pow(Par[4]/Par[3],2)+pow(Par[4]/Par[3],4)))),2.*acos(sqrt(-pow(Par[4]/Par[3],2)+sqrt(1.+pow(Par[4]/Par[3],2)+pow(Par[4]/Par[3],4)))),M_PI/2.};
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = 0, b = 0;
	int j = 0;

	if(Par[3] == 0)
		j = 2;

	for(j; j < 3; j++)
	{
		b = range[j];
		F_a = F_b = 0;	//Start integration at 0
		for(int i = 0; i < 9; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrate1(Integrand, Par, x1[i], Temp)*w[i+1];	//Evaluate k integral at x1
			F_b += Integrate1(Integrand, Par, x3[i], Temp)*w[i+1];	//Evaluate k integral at x3
		}
		F_ave = Integrate1(Integrand, Par, (a+b)/2., Temp)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a);
		a = b;
	}

	return(Answer);	//return the best estimate of the integral on the interval*/
}

long double Integrate1(long double(*Integrand)(long double[6], long double, long double, int), long double Par[6], long double theta, int Temp)
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double x1[9];	//These are the two other points required for 5th order Gaussian quadrature for this interval
	long double x3[9];	//x1 is extended for use in Gauss-Laguerre integration
	long double F_a, F_b, F_ave;
	long double a = 0, b = 0;
	long double Answer = 0;
	long double k = 0;	//Values locating the various values of k where the division by zero gets closest to the real number line
	long double Range[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double gamma = 0;	//These are the widths of the features near 2 Particle on shell
	long double Width;	//Step size for integration
	long double E;		//Largest feature I can find
	long double Value;
	int i = 0, j;

	if(Par[4] > 2.*Par[2])	//If above threshold, locate and guestimate the width of the feature near the division by zero
	{
		long double LocalPar[] = {Par[0],Par[1],Par[2],Par[3],2.*Par[2],Par[5]};
		k = .5*sqrt((pow(Par[4],2)-pow(2.*Par[2],2))*(pow(Par[4],2)+pow(Par[3],2))/(pow(Par[4],2)+pow(Par[3]*sin(theta),2)));
		gamma = -ImProp(LocalPar, k, theta, Temp);
	}
	else if(Par[4] > 2.*Par[2]-.1)	//If near but below threshold, make up a width to get the important things covered.
	{
		long double LocalPar[] = {Par[0],Par[1],Par[2],Par[3],2.*Par[2],Par[5]};
		k = 0;
		gamma = -ImProp(LocalPar, k, theta, Temp);
	}

	if(gamma < 1e-3)
		gamma = 1e-3;

	//cerr << theta << " " << k << " " << gamma << endl;
	while(k+Range[i]*gamma < 0 && k != 0)	//Moves l up until zero[i]+Range[l]*gamma[i] is greater than 0
		i++;

	E = k+(11.8571+.07*Par[3]+.00185714*pow(Par[3],2));
	if(k < (11.8571+.07*Par[3]+.00185714*pow(Par[3],2)))
		a = b = 0;
	else
		a = b = k-(11.8571+.07*Par[3]+.00185714*pow(Par[3],2));
	do
	{
		if(b == 0 && i != 0)	//First peak is closer than 64*gamma to 0
			Width = k+Range[i]*gamma;
		else if(b < k-10 || b > k+10)
			Width = 10;
		else
			Width = 3.;	//No-man's land

		if(a<k-64.*gamma && b+Width>=k-64.*gamma)	//Stutter step before the peak
		{
			Width = k-64.*gamma-b;
			i = 0;	//Resets l before entering the peak
		}
		else if((a>=k-64.*gamma && b<=k+64.*gamma) && i < 16 && b != 0)	//Integrating the peak itself
		{
			Width = gamma*(Range[i+1]-Range[i]);
			i++;	//Leaving the peak prevents illegal space access at Range[l+1]
		}

		b += Width;

		F_a = F_b = 0;
		for(j = 0; j < 9; j++)
		{
			x1[j] = (b+a-Disp[j]*(b-a))/2.; //Actual evaluation points
			x3[j] = (b+a+Disp[j]*(b-a))/2.;

			Value = Integrand(Par, x1[j], theta, Temp);
			//if(Integrand == G_0Int)
			//	cout << x1[j] << " " << theta << " " << Value << endl;
			F_a += Value*w[j+1]; //Evaluate function at x1
			Value = Integrand(Par, x3[j], theta, Temp);
			//if(Integrand == G_0Int)
			//	cout << x3[j] << " " << theta << " " << Value << endl;
			F_b += Value*w[j+1]; //Evaluate function at x3
		}
		Value = Integrand(Par, (a+b)/2., theta, Temp);
		//if(Integrand == G_0Int)
		//	cout << (a+b)/2. << " " << theta << " " << Value << endl;
		F_ave = Value*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}while(b < E);

	return(Answer);	//return the best estimate of the integral on the interval*/
}
