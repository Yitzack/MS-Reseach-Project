//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fstream>
#include<cfloat>
#include"Elements.h"
using namespace std;

void Spectral(long double[3], long double[5], int);	//Function that puts everything in motion

//Integrals that define results and ancillary functions
Elements Dispersion(long double[5], int, Elements);	//Produces the real side of the results
Elements theta_Int(long double[5], int);	//Integrates the theta results
Elements k_Int(long double[5], int, long double);	//Integrates the k momentum results
Elements Folding(long double[5], int, long double, long double);	//Folding integral, energy integral
void Characterize_k_Int(long double[5], int, long double, long double[4], long double[4], int&); //Returns the poles of the k integral's integrands
void Newtons_k_Int(long double, long double, long double, long double, long double[3], long double);	//Newton's method of finding roots for Characterize_k_Int
int Newtons_Test_k_Int(long double, long double, long double, long double, long double[3], long double);	//Trys to determine the number of roots Newton's method  is looking for and where they might be found
void Characterize_Folding(long double[5], int, long double, long double, long double[4], long double[4], int&); //Returns the poles of the folding integral's integrands

//Straight Functions everything is built from
long double Self_Energy(long double, long double, long double, int);	//Single quark self energy
long double Energy(long double, long double, long double, long double);	//Single quark energy, can return momentum if M=0
long double Fermi(long double, int);	//Fermi factor
long double Potential_on(long double[5]);	//On-shell potential for the on-shell T-Matrix
long double Potential1(long double[5], long double, long double);	//Potiential for the numerator of the boson spectrum
long double Potential2(long double[5], long double, long double);	//Potiential for the denominator of the T-Matrix and boson spectrum
long double Quark_Spectrum(long double, long double, long double, int);	//Single quark spectral function
long double Spin_Sum(long double[5]);	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
long double Folding_Integrand(long double[5], long double, long double, long double, int);	//Integrand of the folding integral

#define GAMMA -.3

void Spectral(long double Results[3], long double Par[5], int Temp)	//Function that puts everything in motion
{
	int N_f = 3;
	int N_c = 3;

	Elements Imaginary = theta_Int(Par, Temp);
	Elements Real = Dispersion(Par, Temp, Imaginary);

	long double G_0 = Imaginary.store(0);
	complex<long double> Num(Real.store(1), Imaginary.store(1));
	complex<long double> TMat(Real.store(2), Imaginary.store(2));
	TMat = complex<long double>(Potential_on(Par),0.)/(complex<long double>(1.,0.)-TMat);

	Results[0] = -2.*N_f*N_c/M_PI*(G_0+(Par[0]*pow(Num,2)*TMat).imag());
	Results[1] = TMat.real();
	Results[2] = TMat.imag();
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements Dispersion(long double Par[5], int Temp, Elements Vofs)	//Produces the real side of the results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double DispLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration
	long double wLa[] = {0.07496328305102102808055, 0.1745735743605928864303, 0.2745074833881225250022, 0.3747323102655645620060, 0.4753412526072084401161, 0.5764380939967183636147, 0.6781307242364945406823, 0.7805307978511547593175, 0.8837542316062452388883, 0.9879219194279636096671, 1.0931605619330277996916, 1.1996035979670979427973, 1.3073922479469277349326, 1.416676687469297701993, 1.5276173754408796787012, 1.640386566702889623924, 1.7551700457872174635214, 1.8721691266543402861779, 1.9916029736088098866132, 2.1137113117669909276048, 2.2387576123844772725684, 2.3670328602831611098048, 2.4988600392644108123394, 2.6345995091430390709, 2.7746554982525006307172, 2.9194840027576204632431, 3.0696024758091833914472, 3.2256018156600758204608, 3.3881613374746331979827, 3.5580676615951707296054, 3.7362388067183244743069, 3.9237552950635210172968, 4.1219008467729629867363, 4.3322164077399479741288, 4.5565730632309056055423, 4.7972722621195591678357, 5.057186469320242487569, 5.3399612774797865633198, 5.6503138450512931300331, 5.9944877492232503537552, 6.3809726096501927329094, 6.8216946862388774056326, 7.3340972531892936469048, 7.9450326451948326187906, 8.6987143462393085933469, 9.6750102652900375180015, 11.039313738067347840094, 13.220456867750092021034, 17.982575250664959108273};	//Weight of the function at DispLa
	long double Range[2];	//Various end points of sub-intervals
	long double LocalPar[] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//The local copy of Par to be sent to theta_Int
	Elements Answer(0,0,0);		//Final answer for return
	Elements F_a, F_b, F_ave;	//Sum of ordinate*weights
	long double x1, x2;		//Abscissa
	long double a = -pow(Par[3],2);	//Lower limit of integration for sub-interval -P^2->s,s->2s,2s->M^2,M^2->4M^2 Principle value
	long double b;		//Upper limit of integration for sub-interval Max(2s,4M^2)->inf Gauss-Laugerre
	int i = 0, j;	//Counters

	if(Par[4] < 0)
	{
		Range[0] = 2.*Par[4];
		Range[1] = Par[4];
	}
	else
	{
		Range[0] = Par[4];
		Range[1] = 2.*Par[4];
	}

	if(Temp == 0)
		a = 0;

	while(Range[i] < a) i++;	//May need to discard first end point if 2s<-P^2

	for(i; i < 2; i++)
	{
		b = Range[i];

		F_a.null();
		F_b.null();
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			LocalPar[4] = x1;
			F_a += (theta_Int(LocalPar, Temp)-Vofs)/(Par[5]-x1)*w[j+1];
			LocalPar[4] = x2;
			F_b += (theta_Int(LocalPar, Temp)-Vofs)/(Par[5]-x2)*w[j+1];
		}
		LocalPar[4] = (a+b)/2.;
		F_ave = (theta_Int(LocalPar, Temp)-Vofs)/(Par[5]-(a+b)/2.)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;

		a = b;
	}

	if(Vofs.store(0) != 0 && Temp != 0)	//If Principle value integral is none zero contribution
		Answer += Vofs*log((Par[4]-b)/(Par[4]+pow(Par[3],2)));
	else if(Vofs.store(0) != 0)	//If Principle value integral is none zero contribution
		Answer += Vofs*log((Par[4]-b)/Par[4]);

	F_a.null();
	for(i = 0; i < 49; i++)	//Gauss-Lagerre to get from 4M^2 or 2s to infinity
	{
		x1 = DispLa[i]+a;
		LocalPar[4] = x1;
		F_a += theta_Int(LocalPar, Temp)/(Par[5]-x1)*wLa[i];
	}

	Answer += F_a;

	return(Answer/M_PI);
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements theta_Int(long double Par[5], int Temp)	//Integrates the theta results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {.01,.05,.1,M_PI/2.,M_PI-.1,M_PI-.05,M_PI-.01,M_PI};
	long double x1, x2;	//Abscissa
	Elements F_a, F_b, F_ave;	//Sum of ordinate*weights
	Elements Answer(0,0,0);	//Answer to be returned
	long double a = 0, b;	//Sub-interval limits of integration
	int i, j;	//Counters

	for(i = 0; i < 8; i++)
	{
		b = Range[i];

		F_a.null();
		F_b.null();
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += k_Int(Par, Temp, x1)*sin(x1)*w[j+1];
			F_b += k_Int(Par, Temp, x2)*sin(x2)*w[j+1];
		}
		F_ave = k_Int(Par, Temp, (a+b)/2.)*sin((a+b)/2.)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;
		a = b;
	}

	return(Answer/pow(2.*M_PI,2));
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements k_Int(long double Par[5], int Temp, long double theta)	//Integrates the k momentum results
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	Elements F_a, F_b, F_ave;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Answer to be returned
	Elements PartialAnswer;	//Answer for sub-interval for determining completeness
	long double x1, x2;	//Abscissa
	long double a = 0, b = 0;//Sub-interval limits of integration
	int Poles;	//Number of poles
	long double zero[3];	//The real part of the signular pole
	long double gamma[3];	//The distance to the singular, maybe
	long double Min_upper;	//The integral has to go at least this far
	long double Width;	//Length of the next sub-interval
	long double Early = 0;	//Early change from one pole to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;//The next width that will be used in the event of an early change of poles
	int i = 0, j, l;	//Counters
	Elements holder;

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);

	Min_upper = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//This the upper bound that the vacuum calls for, Partial/total will promote higher as needed

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0) //j!=0 is because the loop will find other early termination
		Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];

	do
	{
		if(b == 0 && j != 0)	//First pole is closer than zero-64*gamma to lower limit of integration
			Width = zero[i]+Range[j]*gamma[i]-b;
		else if((i < Poles && b+100 < zero[i]-64.*gamma[i]) || a-100 > Min_upper || b+100 < Min_upper || (i > 0 && a-100 > zero[i-1]+64.*gamma[i-1]))	//Middle of nowhere intervals
			Width = 100;
		else if((i < Poles && b+50 < zero[i]-64.*gamma[i]) || a-50 > Min_upper || b+50 < Min_upper || (i > 0 && a-50 > zero[i-1]+64.*gamma[i-1]))
			Width = 50;
		else if((i < Poles && b+10 < zero[i]-64.*gamma[i]) || a-10 > Min_upper || b+10 < Min_upper || (i > 0 && a-10 > zero[i-1]+64.*gamma[i-1]))
			Width = 10;
		else
			Width = 3;

		if(j == 8 && i < Poles) //Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i])) //Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Poles != 0 && j < 8 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && (b != 0 || Width == 0))
		{//Integrating a pole and the width wan't resolved by the first condition after do, 24 lines up
			Width = gamma[i]*(Range[j+1]-Range[j]);
			j++;
		}

		if(NextWidth != 0 && NextWidth == NextWidth)	//Resolving early termination of integrating one pole for the next
		{
			Width = NextWidth;
			NextWidth = 0;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Early != 0 && b+Width > Early)
		{
			Width = Early-b;
			Early = 0;
			j = 8-j;
			NextWidth = Width*gamma[i+1]/gamma[i];
			i++;
		}

		b += Width;	//Set next b keeping in mind the upper limit of integration

		if(b > Min_upper && a != Min_upper)	//If (a,b) incompasses Min_upper, make b equal to Min_upper
			b = Min_upper;

		F_a.null();
		F_b.null();
		for(l = 0; l < 9; l++)
		{
			x1 = (b+a-Disp[l]*(b-a))/2.; //Actual evaluation points
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F_a += Folding(Par, Temp, x1, theta)*pow(x1,2)*w[l+1]; //Evaluate function at x1
			F_b += Folding(Par, Temp, x2, theta)*pow(x2,2)*w[l+1]; //Evaluate function at x2
		}
		F_ave = Folding(Par, Temp, (a+b)/2., theta)*pow((a+b)/2.,2)*w[0]; //Evaluate function at (a+b)/2.
		PartialAnswer = (F_a+F_ave+F_b)*(b-a)/(2.);
		Answer += PartialAnswer;
		a = b;
	}while(b <= Min_upper || abs(PartialAnswer/Answer) >= .0001);

	return(Answer);
}

void Characterize_k_Int(long double Par[5], int Temp, long double theta, long double zero[3], long double gamma[3], int &Poles) //Returns the poles of the k integral's integrands
{
	long double holder;
	long double previous[2];
	int i, j = 0, l;

	if(sqrt(Par[4]) > 2.*Par[2])	//Find the poles and estimate a width (distance to pole)
		zero[2] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
	else
		zero[2] = 0;

	zero[0] = zero[1] = zero[2];	//Start the Newton's search for other k poles
	i = Newtons_Test_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta);
	if(i)
	{
		do
		{
			previous[0] = zero[0];
			previous[1] = zero[1];
			Newtons_k_Int(Par[1], Par[4], Par[3], Par[2], zero, theta);
			j++;
		}while((abs(previous[0]/zero[0]-1.) > .001 || abs(previous[1]/zero[1]-1.) > .001) && ((i >= 2 && j <= 10) || i <= 1));	//Keep going while both poles are not known better than 1MeV

		if(j <= 10)
		{
			zero[0] = abs(zero[0]);	//Through investigation of the results, I found that if both peaks are on one condition,
			zero[1] = abs(zero[1]);	//then the lower one will be negative on the other condition
			Poles = 3;
		}
		else
		{
			zero[0] = zero[2];
			gamma[0] = Folding(Par, Temp, zero[0], theta).Min();
			if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
				gamma[0] = 1e-3;
			Poles = 1;
			return;
		}
	}
	else
	{
		zero[0] = zero[2];
		gamma[0] = Folding(Par, Temp, zero[0], theta).Min();
		if(gamma[0] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[0] = 1e-3;
		Poles = 1;
		return;
	}

	for(i = 2; i >= 0; i--)	//Bubble sort
	{
		for(j = 0; j < i; j++)
		{
			if(zero[j] > zero[j+1])
			{
				holder = zero[j+1];
				zero[j+1] = zero[j];
				zero[j] = holder;
			}
			else if(zero[j] == zero[j+1])	//Remove duplicates
			{
				for(l = j; l < 2; l++)
					zero[l] = zero[l+1];
				zero[2] = sqrt(Par[4]+pow(Par[3],2)); //Need to use the biggest finite value of reason or it will be attempted to be sorted to the bottom when it is invalid
				Poles--;
			}
		}
	}

	for(i = 0; i < Poles; i++)
	{
		gamma[i] = Folding(Par, Temp, zero[i], theta).Min();
		if(gamma[i] < 1e-3)	//If width is smaller than this value, make it this big
			gamma[i] = 1e-3;
	}
	return;
}

void Newtons_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = Energy(M,P/2.,k[0],theta)-.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[0],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[0],4))));
	long double fp1 = pow(pow(2.*k[0],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[0],4)*(4.+sqrt(16.+pow(Lambda/k[0],4)))))/(k[0]*sqrt(8.+pow(Lambda/k[0],4)/2.))+(2.*k[0]+P*cos(theta))/(2.*Energy(M,P/2.,k[0],theta));
	long double f2 = -Energy(M,P/2.,-k[1],theta)-.5*(-sqrt(s+pow(P,2))+pow(pow(2.*k[1],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[1],4))));
	long double fp2 = -pow(pow(2.*k[1],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[1],4)*(4.+sqrt(16.+pow(Lambda/k[1],4)))))/(k[1]*sqrt(8.+pow(Lambda/k[1],4)/2.))-(2.*k[1]-P*cos(theta))/(2.*Energy(M,P/2.,-k[1],theta));

	k[0] -= f1/fp1;
	k[1] -= f2/fp2;
	return;
}

int Newtons_Test_k_Int(long double Lambda, long double s, long double P, long double M, long double k[3], long double theta)
{
	long double f1 = Energy(M,P/2.,k[2],theta)-.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[2],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[2],4))));
	long double fp1 = pow(pow(2.*k[2],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4)))))/(k[2]*sqrt(8.+pow(Lambda/k[2],4)/2.))+(2.*k[2]+P*cos(theta))/(2.*Energy(M,P/2.,k[2],theta));
	long double limit1 = sqrt(pow(M,2)+pow(P,2)/4.)-.5*sqrt(s+pow(P,2))+pow(2,-1.5)*Lambda;
	long double f2 = -Energy(M,P/2.,-k[2],theta)+.5*(sqrt(s+pow(P,2))-pow(pow(2.*k[2],4)+pow(Lambda,4),.25)*sqrt(.5+2./sqrt(16.+pow(Lambda/k[2],4))));
	long double fp2 = -pow(pow(2.*k[2],4)+pow(Lambda,4),-.25)*sqrt((pow(Lambda,4)+4.*pow(k[2],4)*(4.+sqrt(16.+pow(Lambda/k[2],4)))))/(k[2]*sqrt(8.+pow(Lambda/k[2],4)/2.))-(2.*k[2]-P*cos(theta))/(2.*Energy(M,P/2.,k[2],theta));
	long double limit2 = -sqrt(pow(M,2)+pow(P,2)/4.)+.5*sqrt(s+pow(P,2))-pow(2,-1.5)*Lambda;
	long double limitp = P*cos(theta)/sqrt(pow(2.*M,2)+pow(P,2));
	int roots = 0;

	if(f1*limit1 < 0 && f2*limit2 < 0)
		roots = 1;
	else if(fp1*limitp < 0 && limit1*limitp < 0)
		roots = 2;
	else if(fp2*limitp < 0 && limit2*limitp < 0)
		roots = 3;
	else if(f1*limit1 < 0 || f2*limit2 < 0)
		cerr << "Unexepected outcome in Newton_Test_k_Int s=" << s << " P=" << P << " theta=" << theta << endl;

	if(roots == 0)
	{
		k[0] = 0;
		k[1] = 0;
	}

	return(roots);
}

//long double Par[5] = {g, Lambda, M, P, s}
Elements Folding(long double Par[5], int Temp, long double k, long double theta)	//Folding integral, energy integral
{
	if(Temp == 0 && sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta) <= Energy(0,Par[3]/2.,k,theta))	//Let's save some time and just return 0, because it is
		return(Elements(0,0,0));
	else if(Par[4]+pow(Par[3],2) <= 0)
		return(Elements(0,0,0));	//Bad data trap and time saver

	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.15896884339395434764996, 0.1527660420658596667789, 0.142606702173606611776, 0.12875396253933622768, 0.1115666455473339947, 0.0914900216224499995, 0.069044542737641227, 0.0448142267656996003, 0.0194617882297264770}; //Weight of the function at Disp
	long double Range[] = {-64,-8,-1,-.5,0,.5,1,8,64};	//Number of gamma from center
	long double a, b;	//Sub-interval limits of integration
	long double Width;	//Next step size
	long double Max;	//Upper limit of integration
	long double Early = 0;	//Early change from one pole to the next, notes the location of change, 0 means no early change
	long double NextWidth = 0;//The next width that will be used in the event of an early change of poles
	Elements F_a, F_b, F_ave;	//Sum of ordinates*weights
	Elements Answer(0,0,0);	//Results to be returned
	long double x1, x2;	//Abscissa
	long double zero[4];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[4];	//Imaginary part of poles
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles

	Characterize_Folding(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about

	if(Temp != 0)	//Assign limits of integration
	{
		a = b = 0;
		Max = sqrt(Par[4]+pow(Par[3],2));
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta);
		Max = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	i = 0;	//Pole counter
	j = 0;	//Range counter
	while(Poles != 0 && zero[i]+Range[j]*gamma[i] < a) j++;	//Pole doesn't need to be checked as all poles are within the limits of integration

	if(zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1] && i+1 < Poles && j != 0) //j!=0 is because the loop will find other early termination
		Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i];

	do
	{
		if((b == 0 || b == Energy(0,Par[3]/2.,k,theta)) && j != 0)	//First pole is closer than zero-64*gamma to lower limit of integration
			Width = zero[i]+Range[j]*gamma[i]-b;
		else if((i < Poles && b+100 < zero[i]-64.*gamma[i]) || b+100 < Max)	//Middle of nowhere intervals
			Width = 100;
		else if((i < Poles && b+50 < zero[i]-64.*gamma[i]) || b+50 < Max)
			Width = 50;
		else if((i < Poles && b+10 < zero[i]-64.*gamma[i]) || b+10 < Max)
			Width = 10;
		else
			Width = 3;

		if(j == 8 && i < Poles) //Last pole has been integrated, time to reset for the next one
		{
			i++;
			j = 0;
		}

		if(Poles != 0 && (a < zero[i]-64.*gamma[i] && b+Width >= zero[i]-64.*gamma[i])) //Stutter step before the next pole
		{
			Width = zero[i]-64.*gamma[i]-b;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Poles != 0 && j < 8 && (a >= zero[i]-64.*gamma[i] && a <= zero[i]+64.*gamma[i]) && !(b == 0 || b == Energy(0,Par[3]/2.,k,theta)))
		{//Integrating a pole and the width wan't resolved by the first condition after do, 24 lines up
			Width = gamma[i]*(Range[j+1]-Range[j]);
			j++;
		}

		if(NextWidth != 0 && NextWidth == NextWidth)	//Resolving early termination of integrating one pole for the next
		{
			Width = NextWidth;
			NextWidth = 0;
			if(i+1 < Poles && zero[i]+64.*gamma[i] > zero[i+1]-64.*gamma[i+1]) //There exists a next pole and their ranges overlap
				Early = zero[i]+(zero[i+1]-zero[i])/(gamma[i]+gamma[i+1])*gamma[i]; //Sets early terminantion point
		}
		else if(Early != 0 && b+Width > Early)
		{
			Width = Early-b;
			Early = 0;
			j = 8-j;
			NextWidth = Width*gamma[i+1]/gamma[i];
			i++;
		}

		b += Width;	//Set next b keeping in mind the upper limit of integration
		if(b > Max)
			b = Max;

		F_a.null();
		F_b.null();
		for(l = 0; l < 9; l++)	//Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F_a += Elements(1, Potential1(Par,x1,k), Potential2(Par,x1,k))*Folding_Integrand(Par,x1,k,theta,Temp)*w[l+1];
			F_b += Elements(1, Potential1(Par,x2,k), Potential2(Par,x2,k))*Folding_Integrand(Par,x2,k,theta,Temp)*w[l+1];
		}
		F_ave = Elements(1, Potential1(Par,(a+b)/2.,k), Potential2(Par,(a+b)/2.,k))*Folding_Integrand(Par,(a+b)/2.,k,theta,Temp)*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;
		a = b;
	}while(b < Max);

	return(Answer);
}

void Characterize_Folding(long double Par[5], int Temp, long double k, long double theta, long double zero[4], long double gamma[4], int &Poles)
{
	long double Lower, Upper; //Limits of integration in Folding, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(Temp != 0)
	{
		Lower = 0;
		Upper = sqrt(Par[4]+pow(Par[3],2));
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta);
		Upper = sqrt(Par[4]+pow(Par[3],2))-Energy(0,Par[3]/2.,-k,theta);
	}

	zero[0] = .5*(sqrt(Par[4]+pow(Par[3],2))+pow(pow(2.*k,4)+pow(Par[1],4),.25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));	//Potential poles, I know exactly where these are at.
	zero[1] = .5*(sqrt(Par[4]+pow(Par[3],2))-pow(pow(2.*k,4)+pow(Par[1],4),.25)*cos(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[0] = .5*(pow(pow(2.*k,4)+pow(Par[1],4),.25)*sin(.5*atan(pow(Par[1]/(2.*k),2))));
	gamma[1] = .5*(pow(pow(2.*k,4)+pow(Par[1],4),.25)*sin(.5*atan(pow(Par[1]/(2.*k),2))));

	holder = GAMMA;
	zero[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));	//Exact vacuum
	zero[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));	//Start of exact vacuum

	if(Temp != 0)	//media estimate
	{
		holder = Self_Energy(Par[2], zero[2], Energy(0,Par[3]/2.,k,theta), Temp)/sqrt(pow(zero[2],2)-pow(Energy(0,Par[3]/2.,k,theta),2));
		zero[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));
		gamma[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));

		holder = Self_Energy(Par[2], zero[3], Energy(0,Par[3]/2.,-k,theta), Temp)/sqrt(pow(zero[3],2)-pow(Energy(0,Par[3]/2.,k,theta),2));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*cos(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
		gamma[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
	}
	else	//Finish up exact vacuum calculations
	{
		gamma[2] = pow(pow(pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,k,theta),2)-pow(holder,2))));
		gamma[3] = pow(pow(pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2)/2.,2)+.25*(pow(2.*Par[2]*holder,2)-pow(holder,4)),.25)*sin(.5*atan(sqrt(pow(2.*Par[2]*holder,2)-pow(holder,4))/(2.*pow(Energy(Par[2],Par[3]/2.,-k,theta),2)-pow(holder,2))));
		zero[3] = sqrt(Par[4]+pow(Par[3],2))-zero[3];
	}

	for(i = 3; i >= 0; i--)	//Bubble sort
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

	i = j = 0;	//Find the first zero greater than 0
	while(zero[i] < Lower) i++;

	while(zero[i] <= Upper && i < 4)	//Move zeroes up to front of array, count off poles within the limits of integration
	{
		zero[j] = zero[i];
		gamma[j] = gamma[i];
		i++;
		j++;
	}
	Poles = j;

	return;
}

//long double Par[5] = {g, Lambda, M, P, s}
long double Self_Energy(long double M, long double omega, long double k, int Temp)	//Single quark self energy
{
	long double Par[6];	//Momentum dependance parameterization
	long double E_0 = Energy(M,k,0,0); //location of lorentzian
	long double Sigma; //size of energy dependance
	long double b1, b2; //width of lorentzian
	long double Delta; //Exponential parameters

	switch(Temp)
	{
		case 0:
			if(omega>k)
				return(GAMMA*sqrt(pow(omega,2)-pow(k,2)));
			else
				return(0);
			break;
		case 1://235.2MeV
			Sigma = -0.1289680072519721;
			b1 = 3.322825262881092;
			b2 = 2.2878310836782014;
			Delta = 1.228601982782018;
			Par[0] = .7359389831810698;
			Par[1] = 7.487501146014314;
			Par[2] = 1.9490238595657456;
			Par[3] = .700215754;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 2://294MeV
			Sigma = -0.09606152620146369;
			b1 = 3.285053276019642;
			b2 = 1.886285913340202;
			Delta = 1.1858269101609233;
			Par[0] = .7409390219065235;
			Par[1] = 7.450458343071824;
			Par[2] = 1.8620618988580635;
			Par[3] = .860810762;
			Par[4] = 10;
			Par[5] = 3;
			break;
		case 3://362MeV
			Sigma = -0.09933548776506283;
			b1 = 3.2108770392083246;
			b2 = 1.3694064180118886;
			Delta = 1.3043774341616825;
			Par[0] = .7426375963204489;
			Par[1] = 7.698646415632565;
			Par[2] = 1.771465704769189;
			Par[3] = .608717852;
			Par[4] = 10;
			Par[5] = 3;
			break;
		default:
			Sigma = 0;
			b1 = 0;
			b2 = 0;
			Delta = 0;
			Par[0] = 1;
			Par[1] = 1;
			Par[2] = 1;
			Par[3] = 0;
			Par[4] = 0;
			Par[5] = 0;
			break;
	}

	return((Par[0]*exp(-pow(k/Par[1],2))+(1-Par[0])*exp(-pow(k/Par[2],2))+Par[3])*(Sigma*exp(Delta+(b1-b2)*(omega-E_0)*E_0/2.-sqrt(b1*b2*pow((omega-E_0)*E_0,2)+pow(Delta+(b1-b2)*(omega-E_0)*E_0/2.,2))))*M/E_0);

}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	return(sqrt(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta)));
}

long double Fermi(long double omega, int T)	//Fermi factor
{
	long double Temp;

	switch(T)
	{
		case 0:
			return(0);
			break;
		case 1:
			Temp = .2352;
			break;
		case 2:
			Temp = .294;
			break;
		case 3:
			Temp = .392;
			break;
		default:
			Temp = 1;
	}

	return(1./(1.+exp(omega/Temp)));
}

long double Potential_on(long double Par[5])	//On-shell potential for the on-shell T-Matrix
{
	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(Par[4]-4.*pow(Par[2],2),2)),2));
}

long double Potential1(long double Par[5], long double omega, long double k)	//Potiential for the numerator of the boson spectrum
{
	return(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2),2)));
}

long double Potential2(long double Par[5], long double omega, long double k)	//Potiential for the denominator of the T-Matrix and boson spectrum
{
	return(Par[0]*pow(pow(Par[1],4)/(pow(Par[1],4)+pow(pow(2.*omega-sqrt(Par[4]+pow(Par[3],2)),2)-4.*pow(k,2),2)),2));
}

long double Quark_Spectrum(long double omega, long double k, long double M, int Temp)	//Single quark spectral function
{
	return(Self_Energy(M, omega, k, Temp)/(pow(pow(omega,2)-pow(k,2)-pow(M,2),2)+pow(Self_Energy(M, omega, k, Temp),2)));
}

long double Spin_Sum(long double Par[5])	//Spinor sum, depends on spin and other quantum numbers of the boson (scalar, pseudo-scale, vector, axial vector), stricktly scalar for now
{
	return(Par[4]);
}

long double Folding_Integrand(long double Par[5], long double omega, long double k, long double theta, int Temp)	//Integrand of the folding integral
{
	return(-Quark_Spectrum(omega, Energy(0, Par[3]/2., k, theta), Par[2], Temp)*Quark_Spectrum(sqrt(Par[4]+pow(Par[3],2))-omega, Energy(0, Par[3]/2., -k, theta), Par[2], Temp)*Spin_Sum(Par));
}
