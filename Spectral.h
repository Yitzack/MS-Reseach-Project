//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
using namespace std;

inline long double Energy(long double, long double, long double, long double); //Returns sqrt(M^2+p^2+k^2-pk cos(theta))
long double ReInt(long double[6], long double [3], long double[5], long double, long double, int); //Returns the real part of the integrad from 2p_0 to infinity
long double ImInt(long double[6], long double [3], long double[5], long double, long double, int); //Returns the imaginary part of the integrad from 2p_0 to infinity
long double Integrate1(long double, long double, long double, long double, long double(*)(long double[6], long double [3], long double[5], long double, long double, int), long double[6], long double [3], long double[5], long double, int); //Integrates the part that it is told to integrate. It uses the difference in evaluating the trapaziod rule and Simpson's rule to pick the points that it integrate between. Since the Gaussian quadrature is very accuate for using very few points (2n-1 order polynomial accurate with n points), it then returns the Gaussian quadrature for 3 points on that subinterval.
long double Integrate2(long double, long double, long double, long double, long double, long double, long double(*)(long double[6], long double [3], long double[5], long double, long double, int), long double[6], long double [3], long double[5], int); //Contains more brains than Integrate1() as it will need to divide the integral into 2 parts and pass the endpoints down for faster times but it uses the same algorithm to acheive its results.
long double Self_Energy(long double[3], long double); //Returns Sigma*(a*Lambda^2/(Lambda^2+P^2)+(1-a)exp(-P^2/sigma)), this is my choice of function for the self-energy
inline long double ReProp(long double[6], long double [3], long double[5], long double, long double, int); //Returns the real part of the propagator
inline long double ImProp(long double[6], long double [3], long double[5], long double, long double, int); //Returns the imaginary part of the propagator
inline long double PropIntegrand(long double, long double, long double, long double, long double, long double, int);
inline long double LawCosines(long double, long double, long double); //Returns the law of cosines for two vectors with an angle inbetween.
inline long double Potential(long double[6], long double, long double); //Returns the potential CC*Lambda^2/(M*(Lambda^2-4k^mu k_mu))
inline long double Common(long double[6], long double [3], long double[5], long double, long double, int); //Returns the common part of propagators
inline long double Self_E_Depends(long double[5], long double); //Contains a function that will give a dependance on E and Temp for the self-energy
complex<long double> TMatrix(long double[6], long double[3], long double[5], long double, int); //Returns the T-matrix
long double Spectral(long double[6], long double[3], long double[5], long double, int); //Returns the spectral function of the T-matrix
long double G_0Int(long double[6], long double [3], long double[5], long double, long double, int); //Returns the integrand for G_0. This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
long double ReDelta_GInt(long double[6], long double [3], long double[5], long double, long double, int); //Returns the real part of the Delta G integrand
long double ImDelta_GInt(long double[6], long double [3], long double[5], long double, long double, int); //Returns the imaginary part of the Delta G integrand
inline long double Potential1(long double[6], long double, long double); //Returns one of the factors of the potiential that is in Potential() without the coupling constant.
long double Fermi(long double[6], long double, long double, int); //Returns the fermi factor where energy is computed
long double FermiProp(long double E, int T)//Returns the fermi factor where the energy is given
long double Analytic(long double, long double); //The analytic spectral function for vacuum, which is constant in P
inline complex<long double> arctan(complex<long double>);
inline complex<long double> arctanh(complex<long double>);
inline complex<long double> arctan(complex<long double>x)
{
	complex<long double> i(0,1);
	complex<long double> one(1,0);
	return(complex<long double>(0,.5)*log((one-i*x)/(one+i*x)));
}

inline complex<long double> arctanh(complex<long double>x)
{
	complex<long double> one(1,0);
	return(complex<long double>(.5,0)*log((one+x)/(one-x)));
}

long double Analytic(long double E, long double Epsilon) //This strictly vacuum, the width has no depandance beside the E, ever. The epsilon is .008GeV for the correct natural width.
{//Stupid complex number object doesn't include binary functions of reals and complex numbers
	long double M = 1.8;
	long double CC = -127.995280691106;
	long double Lambda = 1.4049344847006076;
	Epsilon = Epsilon*pow((E*E-1.258884)/7.984588734864,2.5)*pow(1.618884/(.36+E*E),2);

	complex<long double> SMEpsilon(E*E/4.-M*M,Epsilon*E); //The analytic needs to divide the energy by 2 to match what the integration is doing
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

inline long double Self_Energy(long double Par[3], long double P)
{
	return(Par[0]*exp(-pow(P/Par[1],2))+(1-Par[0])*exp(-pow(P/Par[2],2)));
}

inline long double Self_E_Depends(long double Par[5], long double E)
{
	E /= 2.;
	long double Sigma = Par[0]; //size of energy dependance
	long double gamma = Par[1]; //width of lorentzian
	long double E_0 = Par[2]; //location of lorentzian
	long double a = Par[3], b = Par[4]; //Exponential parameters, length and power
	return(pow(tanh(a*E),2)*Sigma*gamma/M_PI*(1/(pow(E+E_0,2)+pow(gamma,2))-1/(pow(E-E_0,2)+pow(gamma,2))));
}

complex<long double> TMatrix(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double E, int Temp)
{
	complex<long double> Int_Holder; //Holder for the result of the integration, allows it to be calculated once
	long double F_a, a = 0.;
	long double F_b, b = M_PI;
	long double F_c, c = 0.;
	long double F_d, d = 100.;

	F_c = ReInt(Par, SelfPPar, SelfEPar, a, c, Temp); //Inital end points of the boundary
	F_d = ReInt(Par, SelfPPar, SelfEPar, a, d, Temp);
	F_a = Integrate1(c, d, F_c, F_d, ReInt, Par, SelfPPar, SelfEPar, a, Temp);

	F_c = ReInt(Par, SelfPPar, SelfEPar, b, c, Temp); //Inital end points of the boundary
	F_d = ReInt(Par, SelfPPar, SelfEPar, b, d, Temp);
	F_b = Integrate1(c, d, F_c, F_d, ReInt, Par, SelfPPar, SelfEPar, b, Temp);

	Int_Holder = complex<long double>(Integrate2(a, b, F_a, F_b, c, d, ReInt, Par, SelfPPar, SelfEPar, Temp), 0);

	F_c = ImInt(Par, SelfPPar, SelfEPar, a, c, Temp); //Inital end points of the boundary
	F_d = ImInt(Par, SelfPPar, SelfEPar, a, d, Temp);
	F_a = Integrate1(c, d, F_c, F_d, ImInt, Par, SelfPPar, SelfEPar, a, Temp);

	F_c = ImInt(Par, SelfPPar, SelfEPar, b, c, Temp); //Inital end points of the boundary
	F_d = ImInt(Par, SelfPPar, SelfEPar, b, d, Temp);
	F_b = Integrate1(c, d, F_c, F_d, ImInt, Par, SelfPPar, SelfEPar, b, Temp);

	Int_Holder += complex<long double>(0 ,Integrate2(a, b, F_a, F_b, c, d, ImInt, Par, SelfPPar, SelfEPar, Temp));
	Int_Holder = complex<long double>(1.,0.)/(complex<long double>(1.,0.)-Int_Holder); //Integrate once, where it says Parameters[6] in the numerator, I need to put V(p,p') in the event that I didn't get that correct

	return(Int_Holder);
}

long double Spectral(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double E, int Temp)
{
	long double G_0; //The imaginary part of G_0
	complex<long double> TMat;
	complex<long double> Num; //The integral in the numerator of Delta G
	long double F_a, a = 0.;
	long double F_b, b = M_PI;
	long double F_c, c = 0.;
	long double F_d, d = 500.;
	int N_f = 3;
	int N_c = 3;
	int i = 0;
	static long double Done[2][1421];	//Needs to be persistant between function executions
	
	while(i < 1421)	//While I haven't made it to the end of the array
	{
		if(Done[0][i] == E)	//If I find E in the array
			return(Done[1][i]);	//Return the Spectral function value that goes with it
		i++;	//If not found, increament i and try again
	}
	
	TMat = TMatrix(Par, SelfPPar, SelfEPar, E, Temp);
	
	F_c = G_0Int(Par, SelfPPar, SelfEPar, a, c, Temp); //Inital end points of the boundary
	F_d = G_0Int(Par, SelfPPar, SelfEPar, a, d, Temp);
	F_a = Integrate1(c, d, F_c, F_d, G_0Int, Par, SelfPPar, SelfEPar, a, Temp);

	F_c = G_0Int(Par, SelfPPar, SelfEPar, b, c, Temp); //Inital end points of the boundary
	F_d = G_0Int(Par, SelfPPar, SelfEPar, b, d, Temp);
	F_b = Integrate1(c, d, F_c, F_d, G_0Int, Par, SelfPPar, SelfEPar, b, Temp);

	G_0 = Integrate2(a, b, F_a, F_b, c, d, G_0Int, Par, SelfPPar, SelfEPar, Temp);

	F_c = ReDelta_GInt(Par, SelfPPar, SelfEPar, a, c, Temp); //Inital end points of the boundary
	F_d = ReDelta_GInt(Par, SelfPPar, SelfEPar, a, d, Temp);
	F_a = Integrate1(c, d, F_c, F_d, ReDelta_GInt, Par, SelfPPar, SelfEPar, a, Temp);

	F_c = ReDelta_GInt(Par, SelfPPar, SelfEPar, b, c, Temp); //Inital end points of the boundary
	F_d = ReDelta_GInt(Par, SelfPPar, SelfEPar, b, d, Temp);
	F_b = Integrate1(c, d, F_c, F_d, ReDelta_GInt, Par, SelfPPar, SelfEPar, b, Temp);

	Num = complex<long double>(Integrate2(a, b, F_a, F_b, c, d, ReDelta_GInt, Par, SelfPPar, SelfEPar, Temp),0);

	F_c = ImDelta_GInt(Par, SelfPPar, SelfEPar, a, c, Temp); //Inital end points of the boundary
	F_d = ImDelta_GInt(Par, SelfPPar, SelfEPar, a, d, Temp);
	F_a = Integrate1(c, d, F_c, F_d, ImDelta_GInt, Par, SelfPPar, SelfEPar, a, Temp);

	F_c = ImDelta_GInt(Par, SelfPPar, SelfEPar, b, c, Temp); //Inital end points of the boundary
	F_d = ImDelta_GInt(Par, SelfPPar, SelfEPar, b, d, Temp);
	F_b = Integrate1(c, d, F_c, F_d, ImDelta_GInt, Par, SelfPPar, SelfEPar, b, Temp);

	Num += complex<long double>(0,Integrate2(a, b, F_a, F_b, c, d, ImDelta_GInt, Par, SelfPPar, SelfEPar, Temp));

	i = 0;	//Start at 0
	while(i < 1421 && Done[0][i] != 0)	//And go find the first point in the array with no values
		i++;

	if(i < 1421)	//Only do if it will be in the alloted array, otherwise I'll get SigFaulted if I actually do more than 1421 points
	{
		Done[0][i] = E;	//Now that you've found the first 0, store the information here
		Done[1][i] = -2.*N_f*N_c/M_PI*(G_0+(Par[0]*pow(Num,2)*TMat).imag());
	}

	return(-2.*N_f*N_c/M_PI*(G_0+(Par[0]*pow(Num,2)*TMat).imag()));
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

long double FermiProp(long double E, int T)
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

long double G_0Int(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
{
	return(ImProp(Par, SelfPPar, SelfEPar, k, theta, Temp)*sin(theta)*k*k);
}

long double ReDelta_GInt(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ReProp(Par, SelfPPar, SelfEPar, k, theta, Temp)*Common(Par, SelfPPar, SelfEPar, k, theta, Temp)*Potential1(Par, k, theta));
}

long double ImDelta_GInt(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ImProp(Par, SelfPPar, SelfEPar, k, theta, Temp)*Potential1(Par, k, theta));
}

inline long double Potential1(long double Par[6], long double k, long double theta)
{
	return(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.));
}

//Parameters = g, Lambda, M, |vec p|, E=sqrt(s)
long double ReInt(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //Returns the real part of the integrand
{
	return(ReProp(Par, SelfPPar, SelfEPar, k, theta, Temp)*Potential(Par, k, theta)*Common(Par, SelfPPar, SelfEPar, k, theta, Temp)*sin(theta)*k*k*(1.-Fermi(Par, k, theta, Temp)-Fermi(Par, -k, theta, Temp)));
}

long double ImInt(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //Returns the imaginary part of the integrand
{
	return(ImProp(Par, SelfPPar, SelfEPar, k, theta, Temp)*Potential(Par, k, theta)*sin(theta)*k*k);
}

inline long double Common(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //Returns the common part of both propagators
{
	if(Par[4] >= 1.122)
		return(2.*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta)), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta)))-Par[5]*pow((pow(Par[4],2)-1.258884)/7.984588734864,2.5)*pow(1.618884/(.36+pow(Par[4],2)),2)*Par[4], 2)));
	else
		return(2.*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta)), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta))), 2)));
}

inline long double ReProp(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //Returns the real part of the propagator
{
	return(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(SelfEPar, Par[4])*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta)), 2));
}

long double ImProp(long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp) //Returns the imaginary part of the propagator
{
	long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618}; //Dispacement from center for Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481}; //Weight of the function at Disp
	long double zero1 = Energy(Par[2], Par[3]/2., k, theta);	//These are the pionts that may cause the greatest problems, but only if they are between 0 and Par[4]=E=sqrt s
	long double zero2 = Par[4]-Energy(Par[2], Par[3]/2., -k, theta);
	long double Answer = 0;
	long double a = 0;
	long double b;
	long double F_a, F_b, F_ave;
	int i;
	long double x1[24], x3[24];

	if(zero1 > zero2)	//reorder zero1 before zero2
	{
		b = zero1;
		zero1 = zero2;
		zero2 = zero1;
	}
	b = 0;

	if(zero1 > 0 && zero1 < Par[4])	//Integrate from 0 to zero1
	{
		if(zero1 < .9)
			b = 1;
		else
			b = zero1-.1;

		while(b < zero1-.1)
		{
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
				F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b;
			b += 1;
		}

		b = zero1-.1;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;

		b = zero1;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;

		if(zero2-zero1 > .1 && zero2 < Par[4])
		{
			b = zero1+.1;
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
				F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b;
		}
	}

	if(zero2 > 0 && zero2 < Par[4])	//Integrate from 0 or zero1 to zero2, zero2 at this point is going to be bigger than zero1
	{
		if(b < zero2-1.1)
			b += 1;
		else
			b = zero2-.1;

		while(b < zero2-.1)
		{
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
				F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b;
			b += 1;
		}

		b = zero2-.1;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;

		b = zero2;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}
	
	if(Par[4]-zero2 > .1)
	{
		b = zero2+.1;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	//Integrate from 0, zero1, or zero2 to Par[4]=E=sqrt(s)
	while(b < Par[4])
	{
		b += 1;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
			F_a += PropIntegrand(x1[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += PropIntegrand(x3[i], Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = PropIntegrand((a+b)/2., Par[6], SelfPPar[3], SelfEPar[5], k, theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	return(Answer);
}

inline long double PropIntegrand(long double omega, long double Par[6], long double SelfPPar[3], long double SelfEPar[5], long double k, long double theta, int Temp)
{
	return(Par[2]*Par[2]*Self_E_Depends(SelfEPar, omega)*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta))*Self_E_Depends(SelfEPar, Par[4]-omega)*Self_Energy(SelfPPar, LawCosines(Par[3]/2., -k, theta))/(Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta)*(pow(omega-Energy(Par[2], Par[3]/2., k, theta),2)+pow(Self_E_Depends(SelfEPar, omega)*Self_Energy(SelfPPar, LawCosines(Par[3]/2., k, theta)),2))*(pow(Par[4]-omega-Energy(Par[2], Par[3]/2.,-k, theta),2)+pow(Self_E_Depends(SelfEPar, Par[4]-omega)*Self_Energy(SelfPPar, LawCosines(Par[3]/2.,-k, theta)),2))));
}

inline long double Potential(long double Par[6], long double k, long double theta) //Returns the potential CC*(Lambda^2/(M*(Lambda^2-4k^mu k_mu)))^2
{
	return(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.), 2));
}

inline long double Energy(long double M, long double P, long double k, long double theta) //Returns twice the energy sqrt(M^2+(vec P/2+vec k)^2)
{
	return(sqrt(M*M+P*P+k*k-2.*P*k*cos(theta)));
}

inline long double LawCosines(long double P, long double k, long double theta) //Returns the law of cosines for two vectors with an angle in between.
{
	return(sqrt(pow(P,2)+pow(k,2)-2.*P*k*cos(theta)));
}

long double Integrate2(long double a, long double b, long double F_a, long double F_b, long double c, long double d, long double(*Integrand)(long double[6], long double [3], long double[5], long double, long double, int), long double Parameters[6], long double SelfPPar[3], long double SelfEPar[5], int Temp)
{
	long double F_c = Integrand(Parameters, SelfPPar, SelfEPar, c, (a+b)/2., Temp); //Inital end points of the boundary
	long double F_d = Integrand(Parameters, SelfPPar, SelfEPar, d, (a+b)/2., Temp);
	long double F_ave = Integrate1(c, d, F_c, F_d, Integrand, Parameters, SelfPPar, SelfEPar, a/2.+b/2., Temp); //Evaluate k integral at (a+b)/2

	long double Trapazoid = (F_a+F_b)*(b-a)/2.; //Trapazoid rule
	long double Simpsons = (F_a+F_ave*4.+F_b)*(b-a)/6.; //Simpson's rule
	if(abs(Trapazoid-Simpsons)*2./abs(Trapazoid+Simpsons) > 1 && abs(b-a) > 1e-10) //If difference between measurements is too large and the differnce between the two points is large enough. The accuracy needs to be better than .00005 and the resolution equal to 1e-18 (segfault if too small). Possibly stronger requirements now that epsilon is -.000001
		return(Integrate2(a, a/2.+b/2., F_a, F_ave, c, d, Integrand, Parameters, SelfPPar, SelfEPar, Temp)+Integrate2(a/2.+b/2., b, F_ave, F_b, c, d, Integrand, Parameters, SelfPPar, SelfEPar, Temp)); //subdivide the interval and return integral of two sub-intervals
	else	//else
	{
		long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618}; //Dispacement from center{sqrt(.6)};//
		long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481}; //Weight of data point{8./9.,5./9.};//
		long double x1[24]; //These are the two other points required for 5th order Gaussian quadrature for this interval
		long double x3[24];
		
		F_a = F_b = 0; //Start integration at 0
		for(int i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_c = Integrand(Parameters, SelfPPar, SelfEPar, c, x1[i], Temp); //Inital end points of the boundary
			F_d = Integrand(Parameters, SelfPPar, SelfEPar, d, x1[i], Temp);
			F_a += Integrate1(c, d, F_c, F_d, Integrand, Parameters, SelfPPar, SelfEPar, x1[i], Temp)*w[i+1]; //Evaluate k integral at x1

			F_c = Integrand(Parameters, SelfPPar, SelfEPar, c, x3[i], Temp); //Inital end points of the boundary
			F_d = Integrand(Parameters, SelfPPar, SelfEPar, d, x3[i], Temp);
			F_b += Integrate1(c, d, F_c, F_d, Integrand, Parameters, SelfPPar, SelfEPar, x3[i], Temp)*w[i+1]; //Evaluate k integral at x3
		}
		
	return((F_a+w[0]*F_ave+F_b)*(b-a)/(2.)); //return the best estimate of the integral on the interval*/
	}
}

long double Integrate1(long double a, long double b, long double F_a, long double F_b, long double(*Integrand)(long double[6], long double [3], long double[5], long double, long double, int), long double Parameters[6], long double SelfPPar[3], long double SelfEPar[5], long double theta, int Temp)
{
	long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618}; //Dispacement from center for Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481}; //Weight of the function at Disp
	long double DispLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409}; //Displacement from 0 for Gauss-Laguerre integration
	long double wLa[] = {0.07496328305102102808055, 0.1745735743605928864303, 0.2745074833881225250022, 0.3747323102655645620060, 0.4753412526072084401161, 0.5764380939967183636147, 0.6781307242364945406823, 0.7805307978511547593175, 0.8837542316062452388883, 0.9879219194279636096671, 1.0931605619330277996916, 1.1996035979670979427973, 1.3073922479469277349326, 1.416676687469297701993, 1.5276173754408796787012, 1.640386566702889623924, 1.7551700457872174635214, 1.8721691266543402861779, 1.9916029736088098866132, 2.1137113117669909276048, 2.2387576123844772725684, 2.3670328602831611098048, 2.4988600392644108123394, 2.6345995091430390709, 2.7746554982525006307172, 2.9194840027576204632431, 3.0696024758091833914472, 3.2256018156600758204608, 3.3881613374746331979827, 3.5580676615951707296054, 3.7362388067183244743069, 3.9237552950635210172968, 4.1219008467729629867363, 4.3322164077399479741288, 4.5565730632309056055423, 4.7972722621195591678357, 5.057186469320242487569, 5.3399612774797865633198, 5.6503138450512931300331, 5.9944877492232503537552, 6.3809726096501927329094, 6.8216946862388774056326, 7.3340972531892936469048, 7.9450326451948326187906, 8.6987143462393085933469, 9.6750102652900375180015, 11.039313738067347840094, 13.220456867750092021034, 17.982575250664959108273}; //Weight of the function at DispLa
	long double x1[49]; //These are the two other points required for 5th order Gaussian quadrature for this interval
	long double x3[24]; //x1 is extended for use in Gauss-Laguerre integration
	long double F_ave;
	long double Answer = 0;
	int i, j, start = 0;

	if(Parameters[4] > 2.*Parameters[2]) //If above threshold, execute this method designed for divisions by zero closely approching the real number line
	{
		long double k, k_min, k_max; //Values locating the various values of k where the division by zero gets closest to the real number line
		long double distance[] = {1e-1, 7.5e-2, 5e-2, 2.5e-2, 1e-2, 7.5e-3, 5e-3, 2.5e-3, 1e-3, 5e-4, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17}; //magic numbers that indicates the distance from the near division by zero.
		k = .5*sqrt((pow(Parameters[4],2)-pow(2.*Parameters[2],2))*(pow(Parameters[4],2)+pow(Parameters[3],2))/(pow(Parameters[4],2)+pow(Parameters[3]*sin(theta),2)));
		k_min = .5*sqrt((pow(Parameters[4],2)-pow(2.*Parameters[2],2))*(pow(Parameters[4],2)+pow(Parameters[3],2))/(pow(Parameters[4],2)+pow(Parameters[3],2)));
		k_max = .5*sqrt((pow(Parameters[4],2)-pow(2.*Parameters[2],2))*(pow(Parameters[4],2)+pow(Parameters[3],2))/(pow(Parameters[4],2)));
		while(2.*distance[start] > k_min && start < 14) //Finds the starting value that won't over run the 0 lower boundary
		{
			start++;
		}
		
		a = b = 0; //0GeV to near divsion by zero line

		while(b+25 < k_min-2.*distance[start]) //Do the interval 25GeV at a time until k_min-25 is reached, k_min may be out a fair distance
		{
			b += 25;
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		b = k_min-2.*distance[start];
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);

		for(j = start; j < 24; j++)
		{
			a = b; //previous location to set distance from division by zero
			b = k-distance[j];
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
		
				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}
	
		a = b; //last near divsion by zero to division by zero
		b = k;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
			F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	
		for(j = 23; j >= 0; j--)
		{
			a = b; //last value to set distance from division by zero
			b = k+distance[j];
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}
	
		a = b; //near divsion by zero to near division by zero line
		b = k_max+distance[0]*2.;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
			F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	
		while(b < 660) //Do the integration 25GeV at time until 500GeV is reached. k_max may be a fair distance from 660GeV
		{
			a = b; //near divsion by zero line to +100GeV
			b += 25;
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b; //near divsion by zero line+100GeV to 500GeV, or back to 500GeV according to P and if the last integral came up short or when past 500GeV
		b = 660;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
			F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	}
	else
	{
		long double distance[] = {2.5,5,7.5,10,50}; //magic numbers that indicates the distance from k=0GeV
		a = 0; //0GeV to 2.5GeV
		for(j = 0; j < 5; j++)
		{
			b = distance[j]; //New upper boundary
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b; //New lower boundary
		}
		
		while(b < 660) //Do the integration 25GeV at time until 500GeV is reached. k_max may be a fair distance from 660GeV
		{
			a = b; //near divsion by zero line to +100GeV
			b += 25;
			F_a = F_b = 0; //Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
				F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
			}
			F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b; //near divsion by zero line+100GeV to 500GeV, or back to 500GeV according to P and if the last integral came up short or when past 500GeV
		b = 660;
		F_a = F_b = 0; //Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.; //Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*w[i+1]; //Evaluate function at x1
			F_b += Integrand(Parameters, SelfPPar, SelfEPar, x3[i], theta, Temp)*w[i+1]; //Evaluate function at x3
		}
		F_ave = Integrand(Parameters, SelfPPar, SelfEPar, (a+b)/2., theta, Temp)*w[0]; //Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	}

	//Gauss-Laguerre integration, needs to be done by both integral methods
	a = b; //Make b the new lower boundary
	F_a = 0; //Start integration at 0
	for(i = 0; i < 49; i++)
	{
		x1[i] = DispLa[i]+a; //Actual evaluation points

		F_a += Integrand(Parameters, SelfPPar, SelfEPar, x1[i], theta, Temp)*wLa[i]; //Evaluate function at x1
	}
	Answer += F_a;

	return(Answer); //return the best estimate of the integral on the interval*/
}
