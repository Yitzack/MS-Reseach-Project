//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
using namespace std;

inline long double Energy(long double, long double, long double, long double);	//Returns sqrt(M^2+p^2+k^2-pk cos(theta))
long double ReInt(long double[6], long double, long double, int);	//Returns the real part of the integrad from 2p_0 to infinity
long double ImInt(long double[6], long double, long double, int);	//Returns the imaginary part of the integrad from 2p_0 to infinity
long double Integrate1(long double(*)(long double[6], long double, long double, int), long double[6], long double, int);	//Integrates the part that it is told to integrate. It uses the difference in evaluating the trapaziod rule and Simpson's rule to pick the points that it integrate between. Since the Gaussian quadrature is very accuate for using very few points (2n-1 order polynomial accurate with n points), it then returns the Gaussian quadrature for 3 points on that subinterval.
long double Integrate2(long double, long double, long double, long double, long double, long double, long double(*)(long double[6], long double, long double, int), long double[6], int);	//Contains more brains than Integrate1() as it will need to divide the integral into 2 parts and pass the endpoints down for faster times but it uses the same algorithm to acheive its results.
long double Self_Energy(int, long double); //Returns Sigma*(a*Lambda^2/(Lambda^2+P^2)+(1-a)exp(-P^2/sigma)), this is my choice of function for the self-energy
inline long double ReProp(long double[6], long double, long double, int);	//Returns the real part of the propagator
inline long double ImProp(long double[6], long double, long double, int);	//Returns the imaginary part of the propagator
inline long double LawCosines(long double, long double, long double);	//Returns the law of cosines for two vectors with an angle inbetween.
inline long double Potential(long double[6], long double, long double);	//Returns the potential CC*Lambda^2/(M*(Lambda^2-4k^mu k_mu))
inline long double Common(long double[6], long double, long double, int);	//Returns the common part of propagators
inline long double Self_E_Depends(long double, int);	//Contains a function that will give a dependance on E and Temp for the self-energy
complex<long double> TMatrix(long double, long double, long double, int);	//Returns the T-matrix for a given M, P, E=sqrt(s), and T
long double Spectral(long double, long double, long double, int);	//Returns the spectral function of the T-matrix
long double G_0Int(long double[6], long double, long double, int);	//Returns the integrand for G_0. This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
long double ReDelta_GInt(long double[6], long double, long double, int);	//Returns the real part of the Delta G integrand
long double ImDelta_GInt(long double[6], long double, long double, int);	//Returns the imaginary part of the Delta G integrand
inline long double Potential1(long double[6], long double, long double);	//Returns one of the factors of the potiential that is in Potential() without the coupling constant.
long double Fermi(long double[6], long double, long double, int);	//Returns the fermi factor (1-2n_f(E;T))
long double Analytic(long double, long double);	//The analytic spectral function for vacuum, which is constant in P
long double ApproxSpectral(long double, long double);
inline complex<long double> arctan(complex<long double>);
inline complex<long double> arctanh(complex<long double>);

long double Lambda = 1.4049344847006076;//I found it to be 1.405 GeV
long double CC = -127.995280691106;	//I found it to be -127.973 GeV^-2

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

long double Self_Energy(int Temp, long double P) //Returns Sigma*(a*Lambda^2/(Lambda^2+P^2)+(1-a)exp(-P^2/sigma)) as an approximation to the self-energy
{
	long double a;
	long double sigma1;
	long double sigma2;

	switch(Temp)
	{
		case 0:
			return(0);	//Return 0 because we have zero temp, vacuum
			break;
		case 1:
			a = .754419572;
			sigma1 = 7.381325744;
			sigma2 = 1.761645239;
			break;
		case 2:
			a = .775293966;
			sigma1 = 7.580476948;
			sigma2 = 1.650232458;
			break;
		case 3:
			a = .761585379;
			sigma1 = 7.9799853517;
			sigma2 = 1.6664531561;
			break;
	}
	
	return(a*exp(-P*P/(sigma1*sigma1))+(1-a)*exp(-P*P/(sigma2*sigma2)));
}

inline long double Self_E_Depends(long double E, int Temp)
{
	E /= 2.;
	long double Sigma;	//size of energy dependance
	long double gamma;	//width of lorentzian
	long double E_0;	//location of lorentzian
	long double a, b;	//Exponential parameters, length and power

	switch(Temp)
	{
		case 0:
			return(0);	//Return 0 because we have zero temp, vacuum
			break;
		case 1:
			Sigma = .154607;
			gamma = .313225;
			E_0 = 1.6405;
			a = .97956;
			break;
		case 2:
			Sigma = .139149;
			gamma = .375252;
			E_0 = 1.5307;
			a = 1.03312;
			break;
		case 3:
			Sigma = .175974;
			gamma = .545823;
			E_0 = 1.50913;
			a = 1.44961;
			break;
	}


	return(pow(tanh(a*E),2)*Sigma*gamma/M_PI*(1/(pow(E+E_0,2)+pow(gamma,2))-1/(pow(E-E_0,2)+pow(gamma,2))));
}

long double Spectral(long double M, long double P, long double E, int Temp)
{
	long double Par[6] = {CC, Lambda, M, P, E};//Parameters = g, Lambda, M, |vec P|, E=sqrt(s), epsilon
	long double G_0;	//The imaginary part of G_0
	complex<long double> TMat = TMatrix(M, P, E, Temp);
	complex<long double> Num;	//The integral in the numerator of Delta G
	long double F_a, a = 0.;
	long double F_b, b = M_PI;
	long double F_c, c = 0.;
	long double F_d, d = 500.;
	int N_f = 3;
	int N_c = 3;

	F_a = Integrate1(G_0Int, Par, a, Temp);
	F_b = Integrate1(G_0Int, Par, b, Temp);
	G_0 = Integrate2(a, b, F_a, F_b, c, d, G_0Int, Par, Temp);

	F_a = Integrate1(ReDelta_GInt, Par, a, Temp);
	F_b = Integrate1(ReDelta_GInt, Par, b, Temp);
	Num = complex<long double>(Integrate2(a, b, F_a, F_b, c, d, ReDelta_GInt, Par, Temp),0);

	F_a = Integrate1(ImDelta_GInt, Par, a, Temp);
	F_b = Integrate1(ImDelta_GInt, Par, b, Temp);
	Num += complex<long double>(0,Integrate2(a, b, F_a, F_b, c, d, ImDelta_GInt, Par, Temp));

	return(-2.*N_f*N_c/M_PI*(G_0+(Par[0]*pow(Num,2)*TMat).imag()));
}

long double Fermi(long double Par[6], long double k, long double theta, int T)
{
	long double Temp;	//T_c = .196GeV = 196MeV

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

	return(1./(exp(Energy(Par[2], Par[3]/2., k, theta)/Temp)+1.));	//Fermi factor
}

long double G_0Int(long double Par[6], long double k, long double theta, int Temp)	//This argument sturcture is so that I don't have to reinvent the intgrate functions that are known to work
{
	return(ImProp(Par, k, theta, Temp)*Common(Par, k, theta, Temp)*(1.-Fermi(Par, k, theta, Temp)-Fermi(Par, -k, theta, Temp))*sin(theta)*k*k);
}

long double ReDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ReProp(Par, k, theta, Temp)*Common(Par, k, theta, Temp)*Potential1(Par, k, theta)*(1.-Fermi(Par, k, theta, Temp)-Fermi(Par, -k, theta, Temp)));
}

long double ImDelta_GInt(long double Par[6], long double k, long double theta, int Temp)
{
	return(k*k*sin(theta)*ImProp(Par, k, theta, Temp)*Common(Par, k, theta, Temp)*Potential1(Par, k, theta)*(1.-Fermi(Par, k, theta, Temp)-Fermi(Par, -k, theta, Temp)));
}

inline long double Potential1(long double Par[6], long double k, long double theta)
{
	return(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.));
}

complex<long double> TMatrix(long double M, long double P, long double E, int Temp)
{
	complex<long double> Int_Holder;	//Holder for the result of the integration, allows it to be calculated once
	long double Parameters[6] = {CC, Lambda, M, P, E};//Parameters = g, Lambda, M, |vec P|, E=sqrt(s), epsilon
	long double F_a, a = 0.;
	long double F_b, b = M_PI;
	long double F_c, c = 0.;
	long double F_d, d = 100.;

	F_a = Integrate1(ReInt, Parameters, a, Temp);
	F_b = Integrate1(ReInt, Parameters, b, Temp);
	Int_Holder = complex<long double>(Integrate2(a, b, F_a, F_b, c, d, ReInt, Parameters, Temp), 0);

	F_a = Integrate1(ImInt, Parameters, a, Temp);
	F_b = Integrate1(ImInt, Parameters, b, Temp);
	Int_Holder += complex<long double>(0 ,Integrate2(a, b, F_a, F_b, c, d, ImInt, Parameters, Temp));
	Int_Holder = complex<long double>(1.,0.)/(complex<long double>(1.,0.)-Int_Holder);	//Integrate once, where it says Parameters[6] in the numerator, I need to put V(p,p') in the event that I didn't get that correct

	return(Int_Holder);
}

//Parameters = g, Lambda, M, |vec p|, E=sqrt(s)
long double ReInt(long double Par[6], long double k, long double theta, int Temp)	//Returns the real part of the integrand
{
	return(ReProp(Par, k, theta, Temp)*Potential(Par, k, theta)*Common(Par,k,theta,Temp)*sin(theta)*k*k);
}

long double ImInt(long double Par[6], long double k, long double theta, int Temp)	//Returns the imaginary part of the integrand
{
	return(ImProp(Par, k, theta, Temp)*Potential(Par, k, theta)*Common(Par,k,theta,Temp)*sin(theta)*k*k);
}

inline long double Common(long double Par[6], long double k, long double theta, int Temp)	//Returns the common part of both propagators
{
	return(2.*(1.-Fermi(Par, -k, theta, Temp)-Fermi(Par, k, theta, Temp))*pow(Par[2],2)/pow(2.*M_PI,2)*(1./Energy(Par[2], Par[3]/2., -k, theta)+1./Energy(Par[2], Par[3]/2., k, theta))/(pow(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., -k, theta)), 2), 2)+pow(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., -k, theta)))-.0322*Par[4]*pow((pow(Par[4],2)-1.258884)/7.984588734864,2.5)*pow(1.618884/(.36+pow(Par[4],2)),2), 2)));
}

inline long double ReProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the real part of the propagator
{
	return(pow(Par[4],2)+pow(Par[3],2)-pow(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta), 2)+pow(Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., -k, theta)), 2));
}

inline long double ImProp(long double Par[6], long double k, long double theta, int Temp)	//Returns the imaginary part of the propagator
{
	return(2.*(Energy(Par[2], Par[3]/2., k, theta)+Energy(Par[2], Par[3]/2., -k, theta))*(Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., k, theta))+Self_E_Depends(Par[4],Temp)*Self_Energy(Temp, LawCosines(Par[3]/2., -k, theta)))-.0322*Par[4]*pow((pow(Par[4],2)-1.258884)/7.984588734864,2.5)*pow(1.618884/(.36+pow(Par[4],2));
}

inline long double Potential(long double Par[6], long double k, long double theta)	//Returns the potential CC*(Lambda^2/(M*(Lambda^2-4k^mu k_mu)))^2
{
	return(Par[0]*pow(pow(Par[1],2)/(pow(Par[1],2)+2.*(k*k-pow(Par[2],2)+Energy(Par[2], Par[3]/2., k, theta)*Energy(Par[2], Par[3]/2., -k, theta))-pow(Par[3],2)/2.), 2));
}

inline long double Energy(long double M, long double P, long double k, long double theta)	//Returns twice the energy sqrt(M^2+(vec P/2+vec k)^2)
{
	return(sqrt(M*M+P*P+k*k-2.*P*k*cos(theta)));
}

inline long double LawCosines(long double P, long double k, long double theta)	//Returns the law of cosines for two vectors with an angle in between.
{
	return(sqrt(pow(P,2)+pow(k,2)-2.*P*k*cos(theta)));
}

long double Integrate2(long double a, long double b, long double F_a, long double F_b, long double c, long double d, long double(*Integrand)(long double[6], long double, long double, int), long double Parameters[6], int Temp)
{
	long double F_c = Integrand(Parameters, c, (a+b)/2., Temp);	//Inital end points of the boundary
	long double F_d = Integrand(Parameters, d, (a+b)/2., Temp);
	long double F_ave = Integrate1(Integrand, Parameters, a/2.+b/2., Temp);	//Evaluate k integral at (a+b)/2

	long double Trapazoid = (F_a+F_b)*(b-a)/2.;		//Trapazoid rule
	long double Simpsons = (F_a+F_ave*4.+F_b)*(b-a)/6.;	//Simpson's rule
	if(abs(Trapazoid-Simpsons)*2./abs(Trapazoid+Simpsons) > 1 && abs(b-a) > 1e-10 && abs(Simpsons) > 1e-5)	//If difference between measurements is too large and the differnce between the two points is large enough. The accuracy needs to be better than .00005 and the resolution equal to 1e-18 (segfault if too small). Possibly stronger requirements now that epsilon is -.000001
		return(Integrate2(a, a/2.+b/2., F_a, F_ave, c, d, Integrand, Parameters, Temp)+Integrate2(a/2.+b/2., b, F_ave, F_b, c, d, Integrand, Parameters, Temp)); //subdivide the interval and return integral of two sub-intervals
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

			F_c = Integrand(Parameters, c, x1[i], Temp);	//Inital end points of the boundary
			F_d = Integrand(Parameters, d, x1[i], Temp);
			F_a += Integrate1(Integrand, Parameters, x1[i], Temp)*w[i+1];	//Evaluate k integral at x1

			F_c = Integrand(Parameters, c, x3[i], Temp);	//Inital end points of the boundary
			F_d = Integrand(Parameters, d, x3[i], Temp);
			F_b += Integrate1(Integrand, Parameters, x3[i], Temp)*w[i+1];	//Evaluate k integral at x3
		}

		return((F_a+w[0]*F_ave+F_b)*(b-a)/(2.));	//return the best estimate of the integral on the interval*/
	}
}

long double Integrate1(long double(*Integrand)(long double[6], long double, long double, int), long double Parameters[6], long double theta, int Temp)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center for Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp
	long double DispLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration
	long double wLa[] = {0.07496328305102102808055, 0.1745735743605928864303, 0.2745074833881225250022, 0.3747323102655645620060, 0.4753412526072084401161, 0.5764380939967183636147, 0.6781307242364945406823, 0.7805307978511547593175, 0.8837542316062452388883, 0.9879219194279636096671, 1.0931605619330277996916, 1.1996035979670979427973, 1.3073922479469277349326, 1.416676687469297701993, 1.5276173754408796787012, 1.640386566702889623924, 1.7551700457872174635214, 1.8721691266543402861779, 1.9916029736088098866132, 2.1137113117669909276048, 2.2387576123844772725684, 2.3670328602831611098048, 2.4988600392644108123394, 2.6345995091430390709, 2.7746554982525006307172, 2.9194840027576204632431, 3.0696024758091833914472, 3.2256018156600758204608, 3.3881613374746331979827, 3.5580676615951707296054, 3.7362388067183244743069, 3.9237552950635210172968, 4.1219008467729629867363, 4.3322164077399479741288, 4.5565730632309056055423, 4.7972722621195591678357, 5.057186469320242487569, 5.3399612774797865633198, 5.6503138450512931300331, 5.9944877492232503537552, 6.3809726096501927329094, 6.8216946862388774056326, 7.3340972531892936469048, 7.9450326451948326187906, 8.6987143462393085933469, 9.6750102652900375180015, 11.039313738067347840094, 13.220456867750092021034, 17.982575250664959108273};	//Weight of the function at DispLa
	long double x1[49];	//These are the two other points required for 5th order Gaussian quadrature for this interval
	long double x3[24];	//x1 is extended for use in Gauss-Laguerre integration
	long double F_a, F_b, F_ave;
	long double a, b;
	long double Answer = 0;
	int i, j, start = 0;

	if(Parameters[4] > 2.*Parameters[2])	//If above threshold, execute this method designed for divisions by zero closely approching the real number line
	{
		long double k;	//Values locating the various values of k where the division by zero gets closest to the real number line
		long double distance[] = {1e-1, 7.5e-2, 5e-2, 2.5e-2, 1e-2, 7.5e-3, 5e-3, 2.5e-3, 1e-3, 5e-4, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17};	//magic numbers that indicates the distance from the near division by zero.
		k = .5*sqrt((pow(Parameters[4],2)-pow(2.*Parameters[2],2))*(pow(Parameters[4],2)+pow(Parameters[3],2))/(pow(Parameters[4],2)+pow(Parameters[3]*sin(theta),2)));
		while(2.*distance[start] > k && start < 24)	//Finds the starting value that won't over run the 0 lower boundary
		{
			start++;
		}

		a = 0; b = 0;	//0GeV to near divsion by zero line

		while(b+10 < k-2.*distance[start])	//Do the interval 25GeV at a time until k-25 is reached, k may be out a fair distance
		{
			b += 10;
			F_a = F_b = 0;	//Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
				F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
				F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
			}
			F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		b = k-2.*distance[start];
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
			F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
		}
		F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);

		for(j = start; j < 24; j++)
		{
			a = b;	//previous location to set distance from division by zero
			b = k-distance[j];
			F_a = F_b = 0;	//Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
				F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
			}
			F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b;	//last near divsion by zero to division by zero
		b = k;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
			F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
		}
		F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);


		for(j = 23; j >= 0; j--)
		{
			a = b;	//last value to set distance from division by zero
			b = k+distance[j];
			F_a = F_b = 0;	//Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;

				F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
				F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
			}
			F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b;	//near divsion by zero to near division by zero line
		b = k+distance[0]*2.;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
			F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
		}
		F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);

		while(b < 660)	//Do the integration 25GeV at time until 500GeV is reached. k_max may be a fair distance from 660GeV
		{
			a = b;	//near divsion by zero line to +100GeV
			b += 10;
			F_a = F_b = 0;	//Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
				F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
				F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
			}
			F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b;	//near divsion by zero line+100GeV to 500GeV, or back to 500GeV according to P and if the last integral came up short or when past 500GeV
		b = 660;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
			F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
		}
		F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	}
	else
	{
		a = b = 0;
		while(b < 660)	//Do the integration 25GeV at time until 500GeV is reached. k_max may be a fair distance from 660GeV
		{
			a = b;	//near divsion by zero line to +100GeV
			b += 10;
			F_a = F_b = 0;	//Start integration at 0
			for(i = 0; i < 24; i++)
			{
				x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
				x3[i] = (b+a+Disp[i]*(b-a))/2.;
	
				F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
				F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
			}
			F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
			Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		}

		a = b;	//near divsion by zero line+100GeV to 500GeV, or back to 500GeV according to P and if the last integral came up short or when past 500GeV
		b = 660;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1[i] = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x3[i] = (b+a+Disp[i]*(b-a))/2.;

			F_a += Integrand(Parameters, x1[i], theta, Temp)*w[i+1];	//Evaluate function at x1
			F_b += Integrand(Parameters, x3[i], theta, Temp)*w[i+1];	//Evaluate function at x3
		}
		F_ave = Integrand(Parameters, (a+b)/2., theta, Temp)*w[0];	//Evaluate the function at the center
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
	}

	//Gauss-Laguerre integration, needs to be done by both integral methods
	a = b;	//Make b the new lower boundary
	F_a = 0;	//Start integration at 0
	for(i = 0; i < 49; i++)
	{
		x1[i] = DispLa[i]+a;	//Actual evaluation points

		F_a += Integrand(Parameters, x1[i], theta, Temp)*wLa[i];	//Evaluate function at x1
	}
	Answer += F_a;

	return(Answer);	//return the best estimate of the integral on the interval*/
}
