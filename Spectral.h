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
long double Integrate1(long double(*)(long double[6], long double, long double, int), long double[6], long double, int, bool);	//Integrates the part that it is told to integrate. It uses the difference in evaluating the trapaziod rule and Simpson's rule to pick the points that it integrate between. Since the Gaussian quadrature is very accuate for using very few points (2n-1 order polynomial accurate with n points), it then returns the Gaussian quadrature for 3 points on that subinterval.
long double Integrate2(long double, long double, long double, long double, long double(*)(long double[6], long double, long double, int), long double[6], int);	//Contains more brains than Integrate1() as it will need to divide the integral into 2 parts and pass the endpoints down for faster times but it uses the same algorithm to acheive its results.
long double Self_Energy(long double, long double, long double, long double, long double, int);	//Returns the Self-Energy
long double Self_P_Depends(int,long double, long double,  long double);	//Returns the momentum dependance of the self-energy
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

long double Self_Energy(long double s, long double P, long double E, long double q, long double M, int Temp)
{
	long double Ans = M*Self_E_Depends(Temp, E, q, M)*Self_P_Depends(Temp, q, s, P)/Energy(M,q,0,0);
	return(Ans);
}

long double Self_P_Depends(int Temp, long double q, long double s, long double P)
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
			Par[3] = .0863733905579399/1.8;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 2:
			Par[0] = .7409390219065235;
			Par[1] = 7.450458343071824;
			Par[2] = 1.8620618988580635;
			Par[3] = .115370834826227/1.8;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 3:
			Par[0] = .7426375963204489;
			Par[1] = 7.698646415632565;
			Par[2] = 1.771465704769189;
			Par[3] = .108198924731183/1.8;
			Par[4] = 10;
			Par[5] = 3;
			break;
	}

	long double f;
	if(P <= 100)
		f = 1;
	else if(P > 200)
		f = 0;
	else
		f = 2.*pow(P/100.-1.,3)-3.*pow(P/100.-1.,2)+1.;

	if(s >= .685971426239)
		return(Par[0]*exp(-pow(q/Par[1],2))+(1-Par[0])*exp(-pow(q/Par[2],2))+Par[3]/(1.+exp((Par[4]-q)/Par[5]))*pow((s-.685971426239)/8.5575013086254,(long double)2.5)*f*pow(9.603472734864/(.36+s),2));
	else
		return(Par[0]*exp(-pow(q/Par[1],2))+(1-Par[0])*exp(-pow(q/Par[2],2)));
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

	G_0 = Integrate2(0, M_PI, Integrate1(G_0Int, Par, 0, Temp, false), Integrate1(G_0Int, Par, M_PI, Temp, false), G_0Int, Par, Temp);
	Num = complex<long double>(Integrate2(0, M_PI, Integrate1(ReDelta_GInt, Par, 0, Temp, false), Integrate1(ReDelta_GInt, Par, M_PI, Temp, false), ReDelta_GInt, Par, Temp),0);
	Num += complex<long double>(0,Integrate2(0, M_PI, Integrate1(ImDelta_GInt, Par, 0, Temp, false), Integrate1(ImDelta_GInt, Par, M_PI, Temp, false), ImDelta_GInt, Par, Temp));

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
	if(Temp != 0)
		return((pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2.,-k, theta),2)-Par[3]*Par[3])/(2.*Par[2]*Par[2])*ImProp(Par, k, theta, Temp)*sin(theta)*k*k);
	else
		return((-k*k*Par[4]*sin(theta)/(4.*M_PI*abs(Par[3]*cos(theta)*(Energy(Par[2], Par[3]/2., -k, theta)-Energy(Par[2], Par[3]/2., k, theta))-2.*k*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2.,-k, theta))))));
}

long double ReDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ReProp(Par, k, theta, Temp)*Potential1(Par, k, theta, Temp));
}

long double ImDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	if(Temp != 0)
		return(k*k*sin(theta)*ImProp(Par, k, theta, Temp)*Potential1(Par, k, theta, Temp));
	else
		return(-pow(Par[2]*k,2)*sin(theta)*Potential1(Par, k, theta, Temp)/(4.*M_PI*abs(Par[3]*cos(theta)*(Energy(Par[2], Par[3]/2., -k, theta)-Energy(Par[2], Par[3]/2., k, theta))-2.*k*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2.,-k, theta)))));
}

long double Potential1(long double Par[6], long double k, long double theta, int Temp)
{
	return(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.));
}

complex<long double> TMatrix(long double Parameters[6], int Temp)
{
	complex<long double> Int_Holder;	//Holder for the result of the integration, allows it to be calculated once
	Int_Holder = complex<long double>(Integrate2(0, M_PI, Integrate1(ReInt, Parameters, 0, Temp, false), Integrate1(ReInt, Parameters, M_PI, Temp, false), ReInt, Parameters, Temp), 0);
	Int_Holder += complex<long double>(0, Integrate2(0, M_PI, Integrate1(ImInt, Parameters, 0, Temp, false), Integrate1(ImInt, Parameters, M_PI, Temp, false), ImInt, Parameters, Temp));
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
	if(Temp != 0)
		return(k*k*sin(theta)*ImProp(Par, k, theta, Temp)*Potential(Par, k, theta, Temp));
	else
		return(-pow(Par[2]*k,2)*sin(theta)*Potential(Par, k, theta, Temp)/(4.*M_PI*abs(Par[3]*cos(theta)*(Energy(Par[2], Par[3]/2., -k, theta)-Energy(Par[2], Par[3]/2., k, theta))-2.*k*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2.,-k, theta)))));
}

long double ReProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the real part of the propagator
{
	if(Temp == 0)
	{
 		if(Par[4] >= .685971426239)
 			return(2.*((Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2), 2)+pow(-Par[5]*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*sqrt(Par[4]), 2))));
 		else
 			return(2.*((Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2), 2))));
 	}
 
 	long double f;
 	if(Par[3] > 200)
 		f = 1;
 	else if(Par[3] <= 100)
 		f = 0;
 	else
 		f = 3.*pow(Par[3]/100.-1.,2)-2.*pow(Par[3]/100.-1.,3);
 
 	if(Par[4] >= .685971426239)
		return(2.*((Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2))*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta)))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)),2)+pow(-Par[5]*f*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*sqrt(Par[4]), 2)));
	else if(Par[4]+pow(Par[3],2) < 0)	//Catches an issue where energy is zero but appears to be negative due to a difference very nearlly equal numbers. Would have resulted in cascading nan.
		return(2.*((-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2))*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta)))/(pow(-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)),2)));
	else
		return(2.*((Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2))*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta)))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)),2)));
	
	long double Disp[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.}; //Displacement from center for 9th order Gauss-Legendre integration
	long double w[] = {128./225., (322.+13.*sqrt(70))/900., (322.-13.*sqrt(70))/900.}; //Weight of the function at Disp
	long double Range16[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double Range8[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	long double Range4[] = {-1,-.5,0,.5,1};	//Number of gamma from center
	long double LocalPar[] = {Par[0], Par[1], Par[2], Par[3], Par[4], Par[5]};	//The local copy of Par to be sent to ImProp
	long double zero = sqrt(2.*(k*k+Par[2]*Par[2]-Par[3]*Par[3]/4.+Energy(Par[2],Par[3]/2.,k,theta)*Energy(Par[2],Par[3]/2.,-k,theta)));	//2 particle on-shell
	long double gamma = -2.*Self_Energy(Par[4],Par[3], zero, LawCosines(Par[3]/2., k, theta), Par[2], Temp);	//These are the widths of the features near 2 Particle on shell
	long double Answer = 0;
	long double a = 0;
	long double b;
	long double F_a, F_b, F_ave;
	long double x1[9], x3[9];
	long double Width;	//Step size for integration
	long double E = zero+64.*gamma+3.;	//Largest feature I can find
	long double f0 = ImProp(Par, k, theta, Temp);	//Par[4] is the location of the division by zero
	long double A;	//Starting position
	int i,j,l = 0;
	int version = 16;

	if(Par[4] == 0)
		f0 = 0;

	if(zero-64*gamma-3. < 0)
		a = b = 0;
	else
		a = b = zero-64*gamma-3.;
	A = a;

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
			F_a += 2.*x1[j]*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow(x1[j],2)-Par[4])*w[j+1]; //Evaluate function at x1
			LocalPar[4] = x3[j];
			F_b += 2.*x3[j]*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow(x3[j],2)-Par[4])*w[j+1]; //Evaluate function at x3
		}
		LocalPar[4] = (a+b)/2.;
		F_ave = (a+b)*(ImProp(LocalPar, k, theta, Temp)-f0)/(pow((a+b)/2.,2)-Par[4])*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}while(b < E);

	if(f0 == 0)
		Answer = Answer/M_PI;
	else if(A < sqrt(Par[4]) && sqrt(Par[4]) < b)
		Answer = Answer/M_PI+f0*log((Par[4] - b*b)/(A*A - Par[4]))/M_PI;
	else
		Answer = Answer/M_PI+f0*log((b*b - Par[4])/(A*A - Par[4]))/M_PI;

	return(Answer);
}

long double ImProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the imaginary part of the propagator
{
	if(Temp == 0)
	{
		if(Par[4] >= .685971426239)
			return((-Par[5]*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*sqrt(Par[4]))*(2.*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2), 2)+pow(-Par[5]*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*sqrt(Par[4]), 2))));
		else
			return(0);
	}

	if(sqrt(Par[4]+Par[3]*Par[3])-.02 < Energy(Par[2], Par[3]/2., k, theta) && sqrt(Par[4]+Par[3]*Par[3])-.02 < Energy(Par[2], Par[3]/2., -k, theta))
	{
		if(Par[3] > 100)
			return(0);
	}
	else if(sqrt(Par[4]+Par[3]*Par[3]) > Energy(Par[2], Par[3]/2., k, theta) && sqrt(Par[4]+Par[3]*Par[3]) > Energy(Par[2], Par[3]/2., -k, theta))
	{
	 	long double f;
	 	if(Par[3] > 200)
	 		f = 1;
	 	else if(Par[3] <= 100)
	 		f = 0;
	 	else
	 		f = 3.*pow(Par[3]/100.-1.,2)-2.*pow(Par[3]/100.-1.,3);
 
		if(Par[4] >= .685971426239)
			return(((4.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp))-2*Par[5]*f*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*Par[4])*pow(Par[2],2)/pow(2.*M_PI,2)*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta)))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)),2)+pow(-Par[5]*f*pow((Par[4]-.685971426239)/8.5575013086254,(long double)2.5)*pow(9.603472734864/(.36+Par[4]),2)*Par[4], 2)));
		else
			return(((4.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)))*pow(Par[2],2)/pow(2.*M_PI,2)*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta)))/(pow(Par[4]+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),LawCosines(Par[3]/2.,-k,theta),Par[2],Temp)+Self_Energy(Par[4],Par[3],sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2.,-k, theta),LawCosines(Par[3]/2.,k,theta),Par[2],Temp)),2)));
	}
	else if(sqrt(Par[4]+Par[3]*Par[3]) > Energy(Par[2], Par[3]/2., k, theta) || sqrt(Par[4]+Par[3]*Par[3]) > Energy(Par[2], Par[3]/2., -k, theta))
	{
		if(sqrt(Par[4]+Par[3]*Par[3]) > Energy(Par[2], Par[3]/2., -k, theta))
			return(pow(Par[2]/(2.*M_PI),2)*Rho(sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., -k, theta),Par, k, theta, Temp)*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))/Energy(Par[2], Par[3]/2., -k, theta));
		else
			return(pow(Par[2]/(2.*M_PI),2)*Rho(sqrt(Par[4]+Par[3]*Par[3])-Energy(Par[2], Par[3]/2., k, theta),Par, -k, theta, Temp)*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))/Energy(Par[2], Par[3]/2., k, theta));
	}
	
	if(abs((Par[4]+pow(Par[3],2))/Par[4]) < 1e-12) //zero energy bad data trap, avoid things of the size sqrt(-1e-17)
		return(0);

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	//long double Disp[] = {sqrt(5.-2.*sqrt(10./7.))/3., sqrt(5.+2.*sqrt(10./7.))/3.}; //Displacement from center for 9th order Gauss-Legendre integration
	//long double w[] = {128./225., (322.+13.*sqrt(70))/900., (322.-13.*sqrt(70))/900.}; //Weight of the function at Disp
	//long double Range[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	//long double Range[] = {-1,-.5,0,.5,1};	//Number of gamma from center
	long double* zero;	//These are the points that may cause the greatest problems, but only if they are between 0 and E=sqrt(s+P^2)=sqrt(Par[4]+P^2)
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
	int version = 8;

	Characterize(Par, k, theta, Temp, zero, gamma, Peaks);

	a = b = 0;
	i1 = l = 0;
	i2 = 1;
	while(zero[i1]+Range[l]*gamma[i1] < 0 && Peaks != 0 && l < version)	//Moves l up until zero[i]+Range[l]*gamma[i] is greater than 0
		l++;

	if(zero[i1]+64.*gamma[i1] > zero[i2]-64.*gamma[i2] && i2 < Peaks && l != 0)
		Early = zero[i1]+(zero[i2]-zero[i1])/(gamma[i1]+gamma[i2])*gamma[i1];
	else
		Early = 0;

	do
	{
		if(b == 0 && l != 0)	//First peak is closer than 64*gamma to 0
			Width = zero[i1]+Range[l]*gamma[i1];
		else if(b < 2.*Par[2] || b >= sqrt(Par[4]+Par[3]*Par[3])-2.*Par[2])
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

		if(b > sqrt(Par[4]+Par[3]*Par[3]))
			b = sqrt(Par[4]+Par[3]*Par[3]);

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
	}while(b < sqrt(Par[4]+Par[3]*Par[3]));

	delete[] zero;
	delete[] gamma;

	Answer /= pow(2.*M_PI,2);

	return(Answer);
}

void Characterize(long double Par[6], long double k, long double theta, int Temp, long double*& zero, long double*& gamma, int& Peaks)	//Searches for minimum, then gamma by binary search
{
	float DeltaE = .5;
	long double Array[3][int(sqrt(Par[4]+Par[3]*Par[3])/DeltaE)+1];	//0 is function, 1 is second derivative, 2 is boundaries
	long double Center, Width;
	long double Maxima, MaximaW;	//Tempory storage
	long double Value[4];
	int i, j;
	bool Done;

	Peaks = 0;
	for(i = 0; i < sqrt(Par[4]+Par[3]*Par[3])/DeltaE; i++)	//Evaluate the function on a mesh and take the second derivative
		Array[0][i] = PropIntegrand(i*DeltaE, Par, k, theta, Temp);
	Array[1][0]=(3*Array[0][0]-4*Array[0][1]+Array[0][2])/(2.*DeltaE);
	Array[1][int(sqrt(Par[4]+Par[3]*Par[3])/DeltaE)]=(-Array[0][int(sqrt(Par[4]+Par[3]*Par[3])/DeltaE)]+4*Array[0][int(sqrt(Par[4]+Par[3]*Par[3])/DeltaE)-1]-3*Array[0][int(sqrt(Par[4]+Par[3]*Par[3])/DeltaE)-2])/(2.*DeltaE);
	for(i = 1; i < sqrt(Par[4]+Par[3]*Par[3])/DeltaE-1; i++)
		Array[1][i] = (Array[0][i+1]-Array[0][i-1])/(DeltaE*2.);
	for(i = 1; i < sqrt(Par[4]+Par[3]*Par[3])/DeltaE; i++)	//Count the number of minima by the number of changes from positive to negative second derivatives
	{
		if(Array[1][i-1]*Array[1][i] <= 0 && Array[1][i-1] < Array[1][i])
		{
			Array[2][Peaks] = i*DeltaE;	//List the regions were there exist at least 1 peak
			Peaks++;
		}
	}

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
				cerr << "#Minimization failure, s = " << Par[4] << " P = " << Par[3] << " k = " << k << " theta = " << theta << " Center = " << Center << " Width = " << Width << endl;
			Array[0][j] = Center;
			j++;

			Center = Maxima+MaximaW/2.;
			Width = MaximaW/2.;
			if(!Minimize(Par, k, theta, Temp, Center, Width))
				cerr << "#Minimization failure, s = " << Par[4] << " P = " << Par[3] << " k = " << k << " theta = " << theta << " Center = " << Center << " Width = " << Width << endl;
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
			gamma[i] = -Self_Energy(Par[4],Par[3],sqrt(Par[4]), LawCosines(Par[3],k,theta), Par[2], Temp);
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
	return(-Par[2]*Par[2]/M_PI*Rho(omega, Par, k, theta, Temp)*Rho(sqrt(Par[4]+Par[3]*Par[3])-omega, Par, -k, theta, Temp)*(1.-Fermi(omega, Temp)-Fermi(sqrt(Par[4]+Par[3]*Par[3])-omega, Temp)));
}

long double Rho(long double omega, long double Par[6], long double k, long double theta, int Temp)
{
	return(Self_Energy(Par[4],Par[3],omega, LawCosines(Par[3]/2., k, theta), Par[2], Temp)/(Energy(Par[2], Par[3]/2., k, theta)*(pow(omega-Energy(Par[2], Par[3]/2., k, theta),2)+pow(Self_Energy(Par[4],Par[3],omega, LawCosines(Par[3]/2., k, theta), Par[2], Temp),2))));
}

long double Potential(long double Par[6], long double k, long double theta, int Temp)	//Returns the potential CC*(Lambda^2/(M*(Lambda^2-4k^mu k_mu)))^2
{
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

long double Integrate2(long double a, long double b, long double F_a, long double F_b, long double(*Integrand)(long double[6], long double, long double, int), long double Par[6], int Temp)
{
	long double F_ave = Integrate1(Integrand, Par, a/2.+b/2., Temp, false);	//Evaluate k integral at (a+b)/2

	if(Temp == 0 && Integrand == G_0Int)
	{
		long double k = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin((a+b)/2.),2)));
		F_ave = G_0Int(Par, k, (a+b)/2., Temp);
		if(k != k)
			F_ave = 0;
	}

	long double Trapazoid = (F_a+F_b)*(b-a)/2.;		//Trapazoid rule
	long double Simpsons = (F_a+F_ave*4.+F_b)*(b-a)/6.;	//Simpson's rule
	if(abs(b-a) > M_PI/8. || (abs(Trapazoid-Simpsons)*2./abs(Trapazoid+Simpsons) > 1 && abs(b-a) > M_PI/100.))	//If difference between measurements is too large and the differnce between the two points is large enough. The accuracy needs to be better than .00005 and the resolution equal to 1e-18 (segfault if too small)
		return(Integrate2(a, a/2.+b/2., F_a, F_ave, Integrand, Par, Temp)+Integrate2(a/2.+b/2., b, F_ave, F_b, Integrand, Par, Temp)); //subdivide the interval and return integral of two sub-intervals
	else	//else
	{
		long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center{sqrt(.6)};//
		long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point{8./9.,5./9.};//
		long double x1[24];	//These are the two other points required for 5th order Gaussian quadrature for this interval
		long double x3[24];

		F_a = F_b = 0;	//Start integration at 0
		for(int i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			if(Temp == 0 && Integrand == G_0Int)
			{
				long double k = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(x1[i]),2)));
				F_a += G_0Int(Par, k, x1[i], 0)*w[i+1];
				k = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(x3[i]),2)));
				F_b += G_0Int(Par, k, x3[i], 0)*w[i+1];

				if(F_a != F_a)
				{
					F_a = 0;
					F_b = 0;
				}
			}
			else
			{
				F_a += Integrate1(Integrand, Par, x1[i], Temp, true)*w[i+1];	//Evaluate k integral at x1
				F_b += Integrate1(Integrand, Par, x3[i], Temp, true)*w[i+1];	//Evaluate k integral at x3
			}
		}

		if(Temp == 0 && Integrand == G_0Int)
		{
			long double k = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin((a+b)/2.),2)));
			F_ave = G_0Int(Par, k, (a+b)/2., Temp)*w[0];
			if(k != k)
				F_ave = 0;
		}
		else
			F_ave = Integrate1(Integrand, Par, (a+b)/2., Temp, true)*w[0];

		return((F_a+F_ave+F_b)*(b-a)/(2.));	//return the best estimate of the integral on the interval*/
	}
}

long double Integrate1(long double(*Integrand)(long double[6], long double, long double, int), long double Par[6], long double theta, int Temp, bool Important)
{
	long double Disp97[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};
	long double w97[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481}; //97th order Gauss-Legendre integration
	long double Disp37[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 37th order Gauss-Legendre integration
	long double w37[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double x1[24];	//These are the two other points required for 5th order Gaussian quadrature for this interval
	long double x3[24];	//x1 is extended for use in Gauss-Laguerre integration
	long double F_a, F_b, F_ave;
	long double a = 0, b = 0;
	long double Answer = 0;
	long double PartialAnswer;
	long double k = 0;	//Values locating the various values of k where the division by zero gets closest to the real number line
	long double Range[] = {-64,-32,-16,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,16,32,64};	//Number of gamma from center
	long double gamma = 0;	//These are the widths of the features near 2 Particle on shell
	long double Width;	//Step size for integration
	long double E;		//Largest feature I can find
	int i = 0, j;

	if(sqrt(Par[4]) > 2.*Par[2])	//If above threshold, locate and guestimate the width of the feature near the division by zero
	{
		long double LocalPar[] = {Par[0],Par[1],Par[2],Par[3],2.*Par[2],Par[5]};
		k = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
		gamma = -ImProp(LocalPar, k, theta, Temp);
	}
	else if(sqrt(Par[4]) > 2.*Par[2]-.1)	//If near but below threshold, make up a width to get the important things covered.
	{
		long double LocalPar[] = {Par[0],Par[1],Par[2],Par[3],2.*Par[2],Par[5]};
		k = 0;
		gamma = -ImProp(LocalPar, k, theta, Temp);
	}

	if(gamma < 1e-3)
		gamma = 1e-3;

	while(k+Range[i]*gamma < 0 && k != 0)	//Moves l up until zero[i]+Range[l]*gamma[i] is greater than 0
		i++;

        E = k+(11.8571+.57*Par[3]+.00185714*pow(Par[3],2));
	a = b = 0;

	do
	{
		if(b == 0 && i != 0)	//First peak is closer than 64*gamma to 0
			Width = k+Range[i]*gamma;
		else if(b < k-100 || b > k+100)
			Width = 50;
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
		if(Width > 10 && Important)
		{
			for(j = 0; j < 24; j++)
			{
				x1[j] = (b+a-Disp97[j]*(b-a))/2.; //Actual evaluation points
				x3[j] = (b+a+Disp97[j]*(b-a))/2.;

				F_a += Integrand(Par, x1[j], theta, Temp)*w97[j+1]; //Evaluate function at x1
				F_b += Integrand(Par, x3[j], theta, Temp)*w97[j+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Par, (a+b)/2., theta, Temp)*w97[0]; //Evaluate the function at the center
		}
		else
		{
			for(j = 0; j < 9; j++)
			{
				x1[j] = (b+a-Disp37[j]*(b-a))/2.; //Actual evaluation points
				x3[j] = (b+a+Disp37[j]*(b-a))/2.;

				F_a += Integrand(Par, x1[j], theta, Temp)*w37[j+1]; //Evaluate function at x1
				F_b += Integrand(Par, x3[j], theta, Temp)*w37[j+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Par, (a+b)/2., theta, Temp)*w37[0]; //Evaluate the function at the center
		}

		PartialAnswer = (F_a+F_ave+F_b)*(b-a)/(2.);
		Answer += PartialAnswer;
		a = b;
	}while(b < E || (abs(PartialAnswer/Answer) >= .0000001 && Important));

	return(Answer);	//return the best estimate of the integral on the interval*/
}
