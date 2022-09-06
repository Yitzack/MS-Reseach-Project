//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<cfloat>
#include<complex>
#include<queue>
#include"Around.h"
using namespace std;

struct Dev_Pointer
{
	double* omega;
	double* Fermi;
	double* ImSelf;
	double* ReSelf;
	double* q;
	double* Ordinate;
	double* F;
	double* Par;
};

//Integrals that define results
Around Dispersion(Dev_Pointer, double[], int, double, double, double);			//Dispersion relation for turning ImG_12 into ReG_12
Around Dispersion(Dev_Pointer, double[], int, double, double, double, double, double, Around, int, int);	//Dispersion relation for turning ImG_12 into ReG_12
Around k0_Int(Dev_Pointer, double[], int, double, double);					//k0 integral aka energy integral
Around k0_Int(Dev_Pointer, double[], int, double, double, double, double, int, int);	//k0 integral aka energy integral

__global__ void k0_omega_Fermi_16(double*, double*, double*, double*);				//Energy and Fermi function for 16th order
__global__ void k0_Vaccum_ImSelf_16(double*, double*, double*);					//ImSelf for Vacuum and 16th order
__global__ void k0_Vaccum_ReSelf_16(double*, double*, double*);					//ReSelf for Vacuum and 16th order
__global__ void k0_194_ImSelf_16(double*, double*, double*, double*);				//ImSelf for T=194 MeV and 16th order
__global__ void k0_194_ReSelf_16(double*, double*, double*, double*);				//ReSelf for T=194 MeV and 16th order
__global__ void k0_Ordinate_16(double*, double*, double*, double*, double*, double*, double*);	//Ordinate for 16th order
__global__ void k0_Reduce_16(double*, double*);							//Reduce 16th order

__global__ void k0_omega_Fermi_37(double*, double*, double*, double*);				//Energy and Fermi function for 37th order
__global__ void k0_Vaccum_ImSelf_37(double*, double*, double*);					//ImSelf for Vacuum and 37th order
__global__ void k0_Vaccum_ReSelf_37(double*, double*, double*);					//ReSelf for Vacuum and 37th order
__global__ void k0_194_ImSelf_37(double*, double*, double*, double*);				//ImSelf for T=194 MeV and 37th order
__global__ void k0_194_ReSelf_37(double*, double*, double*, double*);				//ReSelf for T=194 MeV and 37th order
__global__ void k0_Ordinate_37(double*, double*, double*, double*, double*, double*, double*);	//Ordinate for 37th order
__global__ void k0_Reduce_37(double*, double*);	//Reduce 37th order

__global__ void k0_omega_Fermi_97(double*, double*, double*, double*);				//Energy and Fermi function for 97th order
__global__ void k0_Vaccum_ImSelf_97(double*, double*, double*);					//ImSelf for Vacuum and 97th order
__global__ void k0_Vaccum_ReSelf_97(double*, double*, double*);					//ReSelf for Vacuum and 97th order
__global__ void k0_194_ImSelf_97(double*, double*, double*, double*);				//ImSelf for T=194 MeV and 97th order
__global__ void k0_194_ReSelf_97(double*, double*, double*, double*);				//ReSelf for T=194 MeV and 97th order
__global__ void k0_Ordinate_97(double*, double*, double*, double*, double*, double*, double*);	//Ordinate for 97th order
__global__ void k0_Reduce_97(double*, double*);	//Reduce 97th order

//Functions for finding points of interest in the k0 integral
void Characterize_k0_Int(double[], int, double, double, double[], double[], int&);	//Returns the poles of the k0 integral's integrands
double Newton_Method_k0(double, double[], double, double, int, double (*)(double[], double, double, double, int));	//Returns the k0 of the on-shell peak using Newton's method on 1/f
double omega_Width(double, double[], double, double, int, double (*)(double[], double, double, double, int));	//Returns the width of on-shell peak using the assumption of a breit-wigner peak 

//Functions for finding points of interest in the dispersion integral
void Characterize_Dispersion(double[], int, double, double, double, double[], double[], int&);
double sp_Width(double[], double, double, double, int, double (*)(double[], double, double, double, int));	//Breit-Wigner width of the peak

//Functions that return physics for the integrand
void ImSelf_Energy(double, double, double[], int, double[]);		//Returns the imaginary single quark self-energies for both quarks, contains an alternate T=194 MeV solution
double ImSelf_Energy(double, double, double, int);			//Returns the imaginary single quark self-energies for one quark, contains an alternate T=194 MeV solution
void ReSelf_Energy(double, double, double[], int, double[]);		//Returns the real single quark self-energies for both quarks, contains an alternate T=194 MeV solution
void Self_Energy(double, double, double[], int, double[], double[]);	//Returns the complex single quark self-energies for both quarks, is a simple Breit-Wigner self-energy and alternate to those above
double Energy(double, double, double, double);			//Single quark energy, also used to return total momentum by setting M=0
double Fermi(double, int);						//Fermi function
double Set_Temp(int);							//Decodes 0-4 into numeric temprature for Fermi factor
double Imk0_Integrand(double[], double, double, double, int);	//Integrand of the k0 integral for positive energy
__device__ __host__ double sq(double);

void George(double M, double omega[], double k[], int Temp, double Results[])	//Single quark self energy
{
	static double Sigma[2];		//Strength
	static double x0[2], x1[2];	//Centrality markers
	static double gamma[2];		//Width
	static double Shift, M_T;
	static double k_old[2];		//Note on validity of k

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
				Sigma[0] = .257498/sqrt(sq(k[0])+pow(1.33201, 2))+.00762638;
				Sigma[1] = .257498/sqrt(sq(k[1])+pow(1.33201, 2))+.00762638;
				x0[0] = sqrt(sq(k[0])+pow(1.54778+Shift, 2))+.276509;
				x0[1] = sqrt(sq(k[1])+pow(1.54778+Shift, 2))+.276509;
				x1[0] = sqrt(sq(k[0])+pow(1.49799+Shift, 2))+.246719;
				x1[1] = sqrt(sq(k[1])+pow(1.49799+Shift, 2))+.246719;
				gamma[0] = .658734/sqrt(sq(k[0])+pow(3.35217, 2))+.0815109;
				gamma[1] = .658734/sqrt(sq(k[1])+pow(3.35217, 2))+.0815109;
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571/sqrt(sq(k[0])+pow(1.17821, 2))+.00762638;
				Sigma[1] = .212571/sqrt(sq(k[1])+pow(1.17821, 2))+.00762638;
				x0[0] = sqrt(sq(k[0])+pow(1.57536+Shift, 2))+.259147;
				x0[1] = sqrt(sq(k[1])+pow(1.57536+Shift, 2))+.259147;
				x1[0] = sqrt(sq(k[0])+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(sq(k[1])+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .336699/sqrt(sq(k[0])+pow(1.87956, 2))+.0651449;
				gamma[1] = .336699/sqrt(sq(k[1])+pow(1.87956, 2))+.0651449;
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972/sqrt(sq(k[0])+pow(1.41483, 2))+.0101423;
				Sigma[1] = .307972/sqrt(sq(k[1])+pow(1.41483, 2))+.0101423;
				x0[0] = sqrt(sq(k[0])+pow(1.56476+Shift, 2))+.251031;
				x0[1] = sqrt(sq(k[1])+pow(1.56476+Shift, 2))+.251031;
				x1[0] = sqrt(sq(k[0])+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(sq(k[1])+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .550628/sqrt(sq(k[0])+pow(2.43968, 2))+.0981269;
				gamma[1] = .550628/sqrt(sq(k[1])+pow(2.43968, 2))+.0981269;
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131/sqrt(sq(k[0])+pow(1.43308, 2))+.0125796;
				Sigma[1] = .339131/sqrt(sq(k[1])+pow(1.43308, 2))+.0125796;
				x0[0] = sqrt(sq(k[0])+pow(1.55034+Shift, 2))+.257788;
				x0[1] = sqrt(sq(k[1])+pow(1.55034+Shift, 2))+.257788;
				x1[0] = sqrt(sq(k[0])+pow(1.46999+Shift, 2))+.231821;
				x1[1] = sqrt(sq(k[1])+pow(1.46999+Shift, 2))+.231821;
				gamma[0] = .615278/sqrt(sq(k[0])+pow(2.22298, 2))+.143376;
				gamma[1] = .615278/sqrt(sq(k[1])+pow(2.22298, 2))+.143376;
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841/sqrt(sq(k[0])+pow(1.42911, 2))+.0157245;
				Sigma[1] = .304841/sqrt(sq(k[1])+pow(1.42911, 2))+.0157245;
				x0[0] = sqrt(sq(k[0])+pow(1.55511+Shift, 2))+.231105;
				x0[1] = sqrt(sq(k[1])+pow(1.55511+Shift, 2))+.231105;
				x1[0] = sqrt(sq(k[0])+pow(1.44714+Shift, 2))+.20956;
				x1[1] = sqrt(sq(k[1])+pow(1.44714+Shift, 2))+.20956;
				gamma[0] = .862629/sqrt(sq(k[0])+pow(2.67193, 2))+.189598;
				gamma[1] = .862629/sqrt(sq(k[1])+pow(2.67193, 2))+.189598;
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

auto Start_Time = chrono::system_clock::now();
int Allotment = 90;

//Merge Sort. There are a number of semi-sorted lists that need sorting. This will probably beat quick sort under similar conditions. Don't worry about the second element. It doesn't need sorting as either they don't matter or they will be assigned after sorting.
void mergeSort(double List[], int a, int b)
{
	int i, j, k;
	double Temp[(a+b)/2-a+1];

	if(b-a > 1)	//Divide...
	{
		mergeSort(List, a, (a+b)/2);
		mergeSort(List, (a+b)/2+1, b);
	}

	for(i = 0; i <= (a+b)/2-a; i++)	//Copy out the lower half array in prep for copy over
	{
		Temp[i] = List[i+a];
	}

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
#define GAMMA -.015	//Width of single quark propagator
#endif

Around Dispersion(Dev_Pointer Pointers, double Par[], int Temp, double k0, double k, double theta)
{
	double a, b;	//Sub-interval limits of integration
	double Min;	//Lower limit of integration
	double Max = 0;	//Upper limit of principal value integration
	Around ImG12 = k0_Int(Pointers, Par, Temp, k, theta);		//Holder of the ImG12 that belongs to the other half of G12 that is calculated here

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	double Boundary_k_k0[] = {0.0491, .421};
	double Range[] = {-Boundary_k_k0[1], -Boundary_k_k0[0], 0, Boundary_k_k0[0], Boundary_k_k0[1]};	//Number of gamma from center
	double ParLoc[5] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//Local copy of the parameters as Par[4] corrisponds to s and ParLoc[4] is s'

	Around Answer(0, 0);	//Results to be returned
	Around Partial;	//Partial results to examine convergance
	Around Holder;

	double zero[5];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	double gamma[5];	//Imaginary part of poles
	int Poles;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities
	int Principal = 1;	//Principal value integral calculating or not

	Characterize_Dispersion(ParLoc, Temp, k0, k, theta, zero, gamma, Poles);
	double Stops[Poles*5+9];		//Extra stops to ensure correctness

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
	Stops[l] = 4.*pow(k0, 2)-sq(Par[3]);	//Both quarks remain energy positive, should be the center of the fermi functions, (1-f-f)=.5
	Stops[l+1] = 4.*sq(k)+4.*pow(k0, 2)+3.*sq(Par[3])+4.*k*Par[3]*cos(theta)-8.*sqrt(pow(k*k0, 2)+pow(k0*Par[3], 2)+k*Par[3]*pow(k0, 2)*cos(theta));	//Light-like quarks
	Stops[l+2] = 4.*sq(k)+4.*pow(k0, 2)+3.*sq(Par[3])-4.*k*Par[3]*cos(theta)+8.*sqrt(pow(k*k0, 2)+pow(k0*Par[3], 2)-k*Par[3]*pow(k0, 2)*cos(theta));
	Stops[l+3] = Par[4];	//Division by zero of dispersion relation
	Stops[l+4] = 4.*(sq(k)+sq(Par[2])+k*Par[3]*cos(theta));
	Stops[l+5] = 4.*(sq(k)+sq(Par[2])-k*Par[3]*cos(theta));
	Stops[l+6] = -sq(Par[3]);
	Stops[l+7] = 4.*sq(k);

	mergeSort(Stops, 0, l+7);
	Stops[l+8] = Stops[l+7]+100;	//Adds the minimum end point to keep the integration going

	Min = a = b = -sq(Par[3]);	//Start from s'=-P^2
	if(abs(Min-Par[4])<.1 && Min-.1 > -sq(Par[3]))
	{
		Min -= .1;
		a -= .1;
		b -= .1;
	}

	for(i = 0; Stops[i]+3 < Par[4]; i++)	//Replace all stops too close to s with s. Will skip over them in the s' traversal.
		if(abs(b/Par[4]-(double)(1.)) < FLT_EPSILON)
			Stops[i] = Par[4];

	for(i = 0; Stops[i] < Min+3 || Holder == Around(0); i++)	//Remove any stops below the minimum of the limit of integration. Faster to illiminate here than by popping
	{
		i++;
		ParLoc[4] = Stops[i];
		Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
		if(!(Holder == Around(0)))
		{
			i--;
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
			while(i < Intervals && (b == Min || abs(b/a-(double)(1.)) < FLT_EPSILON) || b-a == 0)
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
			Partial = Dispersion(Pointers, Par, Temp, k0, k, theta, a, b, ImG12*Principal, 97, 0);
		else
			Partial = Dispersion(Pointers, Par, Temp, k0, k, theta, a, b, ImG12*Principal, 16, 0);
		Answer += Partial;		//Add the Region to the total
		a = b;
	}while((a < 4.*(sq(k)+sq(Par[2]))+50. && i < Intervals) || Partial/Answer > 1e-6);	//Keep going while intervals aren't exhausted and upper limit of integration not excceeded or until convergance

	if(Max == 0)	//Just in case it terminates before getting to s+100
		Max = a;
	if(abs(ImG12) >= 1e-12)
		return((Answer+ImG12*log(abs((Max-Par[4])/(Par[4]-Min))))/M_PI);
	return(Answer/M_PI);
}

Around Dispersion(Dev_Pointer Pointers, double Par[], int Temp, double k0, double k, double theta, double a, double b, Around ImG12, int order, int deep)
{
//9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight
	double Disp9[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	double Disp97[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	double w63[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	double w97[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight

	Around F[2] = {0, 0};	//Sum of ordinates*weights
	Around Answer(0, 0);	//Results to be returned
	Around Holder;

	double ParLoc[5] = {Par[0], Par[1], Par[2], Par[3], Par[4]};	//Local copy of the parameters as Par[4] corrisponds to s and ParLoc[4] is s'

	switch(order)
	{
	case 16:
		for(int l = 0; l < 5; l++) //Count through points away from center
		{
			ParLoc[4] = (b+a-Disp9[l]*(b-a))/2.;
			Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[l+1];
			ParLoc[4] = (b+a+Disp9[l]*(b-a))/2.;
			Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[l+1];
		}
		ParLoc[4] = (b+a)/2.;
		Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
		F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w9[0];
		F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w16[0];
		break;
	case 97:
		for(int l = 0; l < 32; l++) //Count through points away from center
		{
			ParLoc[4] = (b+a-Disp97[l]*(b-a))/2.;
			Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[l+1];
			ParLoc[4] = (b+a+Disp97[l]*(b-a))/2.;
			Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
			F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[l+1];
			F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[l+1];
		}
		ParLoc[4] = (b+a)/2.;
		Holder = k0_Int(Pointers, ParLoc, Temp, k, theta);
		F[0] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w63[0];
		F[1] += (Holder-ImG12)/(ParLoc[4]-Par[4])*w97[0];
		break;
	}

	Answer = Around(F[1], abs(F[0]-F[1]))*(b-a)/2.;
	if(!(Answer.RelErr() < 1e-6 || Answer.Error() < 1e-7) && deep < 4 && abs(b/a-(double)(1.)) > FLT_EPSILON)
		Answer = Dispersion(Pointers, Par, Temp, k0, k, theta, a, (a+b)/2., ImG12, order, deep+1) + Dispersion(Pointers, Par, Temp, k0, k, theta, (a+b)/2., b, ImG12, order, deep+1);

	return(Answer);
}

Around k0_Int(Dev_Pointer Pointers, double Par[], int Temp, double k, double theta)
{
	double a, b;	//Sub-interval limits of integration
	double Max;	//Upper limit of integration

	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	double Boundary_k_k0[] = {.421, 4.85};
	double Range[] = {-Boundary_k_k0[1], -Boundary_k_k0[0], 0, Boundary_k_k0[0], Boundary_k_k0[1]};	//Number of gamma from center

	Around Answer(0, 0);	//Results to be returned
	Around Partial;		//Partial sum to determine continuation

	double zero[12];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	double gamma[12];	//Imaginary part of poles
	double Caution = (Par[4]-pow(2.*Par[2], 2)>0)?(Par[3]*cos(theta)/2.*sqrt((Par[4]-pow(2.*Par[2], 2))/(Par[4]-pow(Par[3]*sin(theta), 2)))):((Energy(Par[2], Par[3]/2., k, theta)-Energy(0, Par[3]/2., -k, theta))/2.);	//Double on-shell in the k0 direction
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities

	Characterize_k0_Int(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	double Stops[Poles*5+6];					//Intervals that are required by integrating near poles

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
	Stops[l] = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+sq(Par[3]))/2.;	//Lower light-like edge
	Stops[l+1] = sqrt(Par[4]+sq(Par[3]))/2.-Energy(0, Par[3]/2., -k, theta);	//Upper light-like edge
	Stops[l+2] = Energy(0, Par[3]/2., k, theta)+sqrt(Par[4]+sq(Par[3]))/2.;	//Pretty sure this is the negative energy solution of the lower light-like edge
	Stops[l+3] = sqrt(Par[4]+sq(Par[3]))/2.+Energy(0, Par[3]/2., -k, theta);	//Pretty sure this is the negative energy solution of the upper light-like edge
	Stops[l+4] = sqrt(Par[4]+sq(Par[3]))/2.;					//Upper energy boundary (E/2)
	Stops[l+5] = -sqrt(Par[4]+sq(Par[3]))/2.;					//Lower energy boundary (-E/2)

	if(Temp != 0)
	{
		a = b = -sqrt(Par[4]+sq(Par[3]))/2.;	//Lower edge for non-vacuum
		Max = sqrt(Par[4]+sq(Par[3]))/2.;		//Upper edge for non-vacuum
	}
	else
	{
		a = b = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+sq(Par[3]))/2.;	//Lower edge for vacuum
		Max = sqrt(Par[4]+sq(Par[3]))/2.-Energy(0, Par[3]/2., -k, theta);	//Upper edge for vacuum
		if(a>Max)
		{
			a = b = sqrt(Par[4]+sq(Par[3]))/2.-Energy(0, Par[3]/2., -k, theta);	//Lower edge for vacuum
			Max = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+sq(Par[3]))/2.;	//Upper edge for vacuum
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
	while(Stops[j] == Stops[j]+1.)	//Remove Regions that duplicates or below the lower edge
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

		if((abs(a-Caution) < 1 || abs(b-Caution) < 1) && Par[4] > pow(2.*Par[2], 2) && abs(k-.5*sqrt((Par[4]-pow(2.*Par[2], 2))*(Par[4]+sq(Par[3]))/(Par[4]+pow(Par[3]*sin(theta), 2)))) < 1)
			Partial = k0_Int(Pointers, Par, Temp, k, theta, a, b, 97, 0);
		else if((abs(a-Caution) < 1 || abs(b-Caution) < 1))
			Partial = k0_Int(Pointers, Par, Temp, k, theta, a, b, 37, 0);
		else
			Partial = k0_Int(Pointers, Par, Temp, k, theta, a, b, 16, 0);

		Answer += Partial;		//Add the Region to the total
		a = b;
	}while((i < Intervals || abs(Partial/Answer) >= .0001) && a < Max);	//Keep going while intervals aren't exhausted and upper limit of ntegration not excceeded

	return(Answer/M_PI);
}

Around k0_Int(Dev_Pointer Pointers, double Par[], int Temp, double k, double theta, double a, double b, int order, int deep)
{
#ifdef DEBUG
double Disp9[] = {-0.9840853600948424644961729, -0.906179845938663992797627, -0.7541667265708492204408172, -0.538469310105683091036314, -0.2796304131617831934134665, 0, 0.2796304131617831934134665, 0.538469310105683091036314, 0.7541667265708492204408172, 0.906179845938663992797627, 0.9840853600948424644961729};	//Displacement from center
	double Debug[16][130];
#endif
//double w9[] = {0., 0.236926885056189087514264, 0., 0.478628670499366468041292, 0., 128./225., 0., 0.478628670499366468041292, 0., 0.236926885056189087514264, 0.};	//9th order Gauss-Legendre weights
//double w16[]= {0.042582036751081832864509451, 0.11523331662247339402462685, 0.18680079655649265746780003, 0.2410403392286475866999426, 0.27284980191255892234099326, 0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
	double F[2];
	double Par_loc[] = {Par[0], Par[1], Par[2], Par[3], Par[4], k, theta, a, b, Set_Temp(Temp)};
	

	cudaMemcpy((void*)Pointers.Par, (void*)Par_loc, 10*sizeof(double), cudaMemcpyHostToDevice);

	switch(order)
	{
	case 97:
		k0_omega_Fermi_97<<<1,65>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi);	//Energy and Fermi function for 97th order
		switch(Temp)
		{
		case 0:
			k0_Vaccum_ImSelf_97<<<1,130>>>(Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for Vacuum and 97th order
			k0_Vaccum_ReSelf_97<<<1,130>>>(Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for Vacuum and 97th order
			break;
		case 1:
			k0_194_ImSelf_97<<<1,130>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for T=194 MeV and 97th order
			k0_194_ReSelf_97<<<1,130>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for T=194 MeV and 97th order
			break;
		}
		k0_Ordinate_97<<<1,65>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi, Pointers.ImSelf, Pointers.ReSelf, Pointers.Ordinate);	//Ordinate for 97th order
		k0_Reduce_97<<<1,65>>>(Pointers.Ordinate, Pointers.F);								//Reduce 97th order
		break;
	case 37:
		k0_omega_Fermi_37<<<1,25>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi);	//Energy and Fermi function for 37th order
		switch(Temp)
		{
		case 0:
			k0_Vaccum_ImSelf_37<<<1,50>>>(Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for Vacuum and 37th order
			k0_Vaccum_ReSelf_37<<<1,50>>>(Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for Vacuum and 37th order
			break;
		case 1:
			k0_194_ImSelf_37<<<1,50>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for T=194 MeV and 37th order
			k0_194_ReSelf_37<<<1,50>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for T=194 MeV and 37th order
			break;
		}
		k0_Ordinate_37<<<1,25>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi, Pointers.ImSelf, Pointers.ReSelf, Pointers.Ordinate);	//Ordinate for 37th order
		k0_Reduce_37<<<1,25>>>(Pointers.Ordinate, Pointers.F);								//Reduce 37th order
		break;
	case 16:
		k0_omega_Fermi_16<<<1,11>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi);	//Energy and Fermi function for 16th order
#ifdef DEBUG
		for(int i = 0; i < 11; i++)
		{
			double x = (Par_loc[8]+Par_loc[7]+Disp9[i]*(Par_loc[8]-Par_loc[7]))/2.;
			Debug[0][i] = sqrt(Par_loc[4]+sq(Par_loc[3]))/2.+x;
			Debug[0][i+11] = sqrt(Par_loc[4]+sq(Par_loc[3]))/2.-x;
			Debug[1][i] = sqrt(sq(Par_loc[3])/4.+sq(Par_loc[5])+Par_loc[0]*Par_loc[5]*cos(Par_loc[6]));
			Debug[1][i+11] = sqrt(sq(Par_loc[3])/4.+sq(Par_loc[5])-Par_loc[0]*Par_loc[5]*cos(Par_loc[6]));
			Debug[2][i] = Fermi(Debug[0][i],Temp);
			Debug[2][i+11] = Fermi(Debug[0][i+11],Temp);
		}
		cudaMemcpy((void*)Debug[12], (void*)Pointers.Par, 10*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)Debug[3], (void*)Pointers.omega, 22*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)Debug[4], (void*)Pointers.q, 22*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy((void*)Debug[5], (void*)Pointers.Fermi, 22*sizeof(double), cudaMemcpyDeviceToHost);
#endif
		switch(Temp)
		{
		case 0:
			k0_Vaccum_ImSelf_16<<<1,22>>>(Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for Vacuum and 16th order
			k0_Vaccum_ReSelf_16<<<1,22>>>(Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for Vacuum and 16th order
#ifdef DEBUG
			for(int i = 0; i < 11; i++)
			{
				double Store0[2];
				double Store1[2];
				double Store2[2];
				Debug[6][i] = ImSelf_Energy(Par_loc[2], Debug[0][i], Debug[1][i], Temp);
				Debug[6][i+11] = ImSelf_Energy(Par_loc[2], Debug[0][i+11], Debug[1][i+11], Temp);
				Store0[0] = Debug[0][i];
				Store0[1] = Debug[0][i+11];
				Store1[0] = Debug[1][i];
				Store1[1] = Debug[1][i+11];
				George(Par_loc[2], Store0, Store1, Temp, Store2);
				Debug[7][i] = Store2[0];
				Debug[7][i+11] = Store2[1];
			}
			cudaMemcpy((void*)Debug[8], (void*)Pointers.ImSelf, 22*sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy((void*)Debug[9], (void*)Pointers.ReSelf, 22*sizeof(double), cudaMemcpyDeviceToHost);
#endif
			break;
		case 1:
			k0_194_ImSelf_16<<<1,22>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ImSelf);	//ImSelf for T=194 MeV and 16th order
			k0_194_ReSelf_16<<<1,22>>>(Pointers.Par, Pointers.q, Pointers.omega, Pointers.ReSelf);	//ReSelf for T=194 MeV and 16th order
			break;
		}
		k0_Ordinate_16<<<1,11>>>(Pointers.Par, Pointers.omega, Pointers.q, Pointers.Fermi, Pointers.ImSelf, Pointers.ReSelf, Pointers.Ordinate);	//Ordinate for 16th order
#ifdef DEBUG
		for(int i = 0; i < 11; i++)
		{
			Debug[10][i] = -((4.*Debug[6][i]*Debug[6][i+11]*sq(Par_loc[2])*(1.-Debug[2][i]-Debug[2][i+11]))/((sq(sq(Debug[0][i])-sq(Debug[1][i])-sq(Par_loc[2])-2.*Par_loc[2]*Debug[7][i])+sq(Debug[6][i]))*(sq(sq(Debug[0][i+11])-sq(Debug[1][i+11])-sq(Par_loc[2])-2.*Par_loc[2]*Debug[7][i+11])+sq(Debug[6][i+1]))));
			Debug[13][i] = 4.*Debug[6][i]*Debug[6][i+11]*sq(Par_loc[2])*(1.-Debug[2][i]-Debug[2][i+11]);
			Debug[14][i] = (sq(sq(Debug[0][i])-sq(Debug[1][i])-sq(Par_loc[2])-2.*Par_loc[2]*Debug[7][i])+sq(Debug[6][i]));
			Debug[15][i] = (sq(sq(Debug[0][i+11])-sq(Debug[1][i+11])-sq(Par_loc[2])-2.*Par_loc[2]*Debug[7][i+11])+sq(Debug[6][i+1]));
		}
		cudaMemcpy((void*)Debug[11], (void*)Pointers.Ordinate, 11*sizeof(double), cudaMemcpyDeviceToHost);
#endif
		k0_Reduce_16<<<1,11>>>(Pointers.Ordinate, Pointers.F);								//Reduce 16th order
		break;
	}

	cudaMemcpy((void*)F, (void*)Pointers.F, 2*sizeof(double), cudaMemcpyDeviceToHost);
	Around Answer = Around(F[1], abs(F[0]-F[1]))*(b-a)/2.;//Around(F[0])*(b-a)/2.;//
	/*if(Answer.RelErr() > 1e-8 && deep < 4 && abs(b/a-(double)(1.)) > FLT_EPSILON)
		Answer = k0_Int(Par, Temp, k, theta, a, (a+b)/2., order, deep+1) + k0_Int(Par, Temp, k, theta, (a+b)/2., b, order, deep+1);//*/

	return(Answer);

}

//9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight, 11 points
__constant__ double Disp9[] = {-0.9840853600948424644961729, -0.906179845938663992797627, -0.7541667265708492204408172, -0.538469310105683091036314, -0.2796304131617831934134665, 0, 0.2796304131617831934134665, 0.538469310105683091036314, 0.7541667265708492204408172, 0.906179845938663992797627, 0.9840853600948424644961729};	//Displacement from center
__constant__ double w9[] = {0., 0.236926885056189087514264, 0., 0.478628670499366468041292, 0., 128./225., 0., 0.478628670499366468041292, 0., 0.236926885056189087514264, 0.};	//9th order Gauss-Legendre weights
__constant__ double w16[]= {0.042582036751081832864509451, 0.11523331662247339402462685, 0.18680079655649265746780003, 0.2410403392286475866999426, 0.27284980191255892234099326, 0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights
__global__ void k0_omega_Fermi_16(double* Par_globe, double* omega, double* q, double* Fermi)	//Energy and Fermi function for 16th order
{
	__shared__ double x[11];
	__shared__ double Par[10];
	if(threadIdx.x < 10)
		Par[threadIdx.x] = Par_globe[threadIdx.x];

	x[threadIdx.x] = (Par[8]+Par[7]+Disp9[threadIdx.x]*(Par[8]-Par[7]))/2.;	//abscisca
	omega[threadIdx.x] = sqrt(Par[4]+sq(Par[3]))/2.+x[threadIdx.x];
	omega[threadIdx.x+11] = sqrt(Par[4]+sq(Par[3]))/2.-x[threadIdx.x];
	q[threadIdx.x] = sqrt(sq(Par[3])/4.+sq(Par[5])+Par[3]*Par[5]*cos(Par[6]));
	q[threadIdx.x+11] = sqrt(sq(Par[3])/4.+sq(Par[5])-Par[3]*Par[5]*cos(Par[6]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x] = 0;
		else
			Fermi[threadIdx.x] = 1;
	}
	else
		Fermi[threadIdx.x] = 1./(1.+exp(omega[threadIdx.x]/Par[9]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x+11] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x+11] = 0;
		else
			Fermi[threadIdx.x+11] = 1;
	}
	else
		Fermi[threadIdx.x+11] = 1./(1.+exp(omega[threadIdx.x+11]/Par[9]));
}

__global__ void k0_Vaccum_ImSelf_16(double* q, double* omega_globe, double* ImSelf)			//ImSelf for Vacuum and 16th order
{
	__shared__ double omega[22];
	__shared__ double k[22];
	omega[threadIdx.x] = omega_globe[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA;
	else
		ImSelf[threadIdx.x] = 0;
}

__global__ void k0_Vaccum_ReSelf_16(double* q, double* omega_globe, double* ReSelf)			//ReSelf for Vacuum and 16th order
{
	ReSelf[threadIdx.x] = 0;
}

__global__ void k0_194_ImSelf_16(double* Par, double* q, double* omega_globe, double* ImSelf)	//ImSelf for T=194 MeV and 16th order
{
	__shared__ double omega[22];
	__shared__ double k[22];
	__shared__ double ImSigma[22];
	__shared__ static double M;
	__shared__ static double omega0[2];		//Location of central peak
	__shared__ static double Sigma[2];		//Amplitude of energy dependance
	__shared__ static double a[2], b[2];	//Slope of exponential decrease to left and right
	__shared__ static double knee[2];		//Interval to change from left to right side of peak
	__shared__ static double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	__shared__ static double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	omega[threadIdx.x] = omega_globe[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[11] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 11))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/11];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/11] = .569969/sqrt(sq(k[threadIdx.x/11])+sq(1.75236))+.0187484;
		a[threadIdx.x/11] = 12.5349/(sq(k[threadIdx.x/11])+sq(1.63711))+5.026;
		b[threadIdx.x/11] = -291.579/(sq(k[threadIdx.x/11]+15.2519)+sq(.0614821))+3.36681;
		omega0[threadIdx.x/11] = sqrt(sq(1.51443+Shift)+sq(k[threadIdx.x/11]))+.232841;
		knee[threadIdx.x/11] = 3.78956*pow(k[threadIdx.x/11]+1., (double)-.530289)+.305*(tanh((k[threadIdx.x/11]-48.4)/11.1111)+1);
	}

	if((omega[threadIdx.x]-omega0[threadIdx.x/11]+knee[threadIdx.x/11]*(b[threadIdx.x/11]-a[threadIdx.x/11])/(sqrt(a[threadIdx.x/11]*b[threadIdx.x/11])*(a[threadIdx.x/11]+b[threadIdx.x/11])))/knee[threadIdx.x/11] < -4.)
		ImSigma[threadIdx.x] = a[threadIdx.x/11]*(omega[threadIdx.x]-omega0[threadIdx.x/11]+knee[threadIdx.x/11]/sqrt(a[threadIdx.x/11]*b[threadIdx.x/11]));
	else if((omega[threadIdx.x]-omega0[threadIdx.x/11]+knee[threadIdx.x/11]*(b[threadIdx.x/11]-a[threadIdx.x/11])/(sqrt(a[threadIdx.x/11]*b[threadIdx.x/11])*(a[threadIdx.x/11]+b[threadIdx.x/11])))/knee[threadIdx.x/11] > 4.)
		ImSigma[threadIdx.x] = b[threadIdx.x/11]*(omega0[threadIdx.x/11]-omega[threadIdx.x]+knee[threadIdx.x/11]/sqrt(a[threadIdx.x/11]*b[threadIdx.x/11]));
	else	//Lost of precision having been circumvented, the actual value
		ImSigma[threadIdx.x] = -.5*((a[threadIdx.x/11]-b[threadIdx.x/11])*omega0[threadIdx.x/11]-((a[threadIdx.x/11]+b[threadIdx.x/11])*knee[threadIdx.x/11])/sqrt(a[threadIdx.x/11]*b[threadIdx.x/11]))+(a[threadIdx.x/11]-b[threadIdx.x/11])*omega[threadIdx.x]/2-sqrt(pow(((a[threadIdx.x/11]+b[threadIdx.x/11])/2.)*(omega[threadIdx.x]-omega0[threadIdx.x/11]+((a[threadIdx.x/11]-b[threadIdx.x/11])*knee[threadIdx.x/11])/(sqrt(a[threadIdx.x/11]*b[threadIdx.x/11])*(a[threadIdx.x/11]+b[threadIdx.x/11]))), 2)+pow(knee[threadIdx.x/11], 2));

#ifdef HALF
	ImSigma[threadIdx.x] = -M*Sigma[threadIdx.x/11]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#else
	ImSigma[threadIdx.x] = -2.*M*Sigma[threadIdx.x/11]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#endif

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA+ImSigma[threadIdx.x];
	else
		ImSelf[threadIdx.x] = ImSigma[threadIdx.x];
}

__global__ void k0_194_ReSelf_16(double* Par, double* q, double* omega_globe, double* ReSelf)	//ReSelf for T=194 MeV and 16th order
{
	__shared__ double omega[22];
	__shared__ double k[22];
	__shared__ static double M;
	__shared__ static double Sigma[2];		//Strength
	__shared__ static double x0[2], x1[2];	//Centrality markers
	__shared__ static double gamma[2];		//Width
	__shared__ static double Shift, M_T;
	__shared__ static double k_old[2];		//Note on validity of k

	omega[threadIdx.x] = omega[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[11] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 11))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/11];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/11] = .212571/sqrt(sq(k[threadIdx.x])+sq(1.17821))+.00762638;
		x0[threadIdx.x/11] = sqrt(sq(k[threadIdx.x])+sq(1.57536+Shift))+.259147;
		x1[threadIdx.x/11] = sqrt(sq(k[threadIdx.x])+sq(1.50194+Shift))+.222526;
		gamma[threadIdx.x/11] = .336699/sqrt(sq(k[threadIdx.x])+sq(1.87956))+.0651449;
	}

#ifdef HALF
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/11]*(omega[threadIdx.x]-x0[threadIdx.x/11])/(pow(omega[threadIdx.x]-x1[threadIdx.x/11], 2)+gamma[threadIdx.x/11])/2.;
#else
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/11]*(omega[threadIdx.x]-x0[threadIdx.x/11])/(pow(omega[threadIdx.x]-x1[threadIdx.x/11], 2)+gamma[threadIdx.x/11]);
#endif
}

__global__ void k0_Ordinate_16(double* Par, double* omega, double* q, double* fermi, double* ImSelf_globe, double* ReSelf, double* Ordinate)	//Ordinate for 16th order
{
	__shared__ double ImSelf[11][2];
	__shared__ double M[11];
	
	M[threadIdx.x] = Par[2];
	ImSelf[threadIdx.x][0] = ImSelf_globe[threadIdx.x];
	ImSelf[threadIdx.x][1] = ImSelf_globe[threadIdx.x+11];

	Ordinate[threadIdx.x] = -((4.*ImSelf[threadIdx.x][0]*ImSelf[threadIdx.x][1]*sq(M[threadIdx.x])*(1.-fermi[threadIdx.x]-fermi[threadIdx.x+11]))/((sq(sq(omega[threadIdx.x])-sq(q[threadIdx.x])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x])+sq(ImSelf[threadIdx.x][0]))*(sq(sq(omega[threadIdx.x+11])-sq(q[threadIdx.x+11])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x+11])+sq(ImSelf[threadIdx.x][1]))));
}

__global__ void k0_Reduce_16(double* Ordinate_globe, double* Answer)	//Reduce 16th order
{
	__shared__ double lower_order[11];
	__shared__ double higher_order[11];
	__shared__ double Ordinate[11];

	Ordinate[threadIdx.x] = Ordinate_globe[threadIdx.x];
	lower_order[threadIdx.x] = Ordinate[threadIdx.x]*w9[threadIdx.x];	//store the lower order value
	higher_order[threadIdx.x] = Ordinate[threadIdx.x]*w16[threadIdx.x];	//store the higher order value

	for(int i = 0; i < 4; i++)	//reduce
	{
		if(threadIdx.x+pow(2,i) < 11)
		{
			lower_order[threadIdx.x] += lower_order[int(threadIdx.x+pow(2,i))];
			higher_order[threadIdx.x] += higher_order[int(threadIdx.x+pow(2,i))];
		}
	}

	if(threadIdx.x == 0)
	{
		Answer[0] = lower_order[0];
		Answer[1] = higher_order[0];
	}
}

//23th order Gauss-Legendre/37th order Gauss-Kronrod integration, 25 points
__constant__ double Disp37[] = {-0.99693392252959542691235, -0.98156063424671925069055, -0.95053779594312129654906, -0.90411725637047485667847, -0.84355812416115324479214, -0.76990267419430468703689, -0.68405989547005589394493, -0.58731795428661744729670, -0.48133945047815709293594, -0.36783149899818019375269, -0.24850574832046927626779, -0.12523340851146891547244, 0, 0.12523340851146891547244, 0.24850574832046927626779, 0.36783149899818019375269, 0.48133945047815709293594, 0.58731795428661744729670, 0.68405989547005589394493, 0.76990267419430468703689, 0.84355812416115324479214, 0.90411725637047485667847, 0.95053779594312129654906, 0.98156063424671925069055, 0.99693392252959542691235};	//Displacement from center
__constant__ double w23[] = {0, 0.047175336386511827194616, 0, 0.106939325995318430960255, 0, 0.16007832854334622633465, 0, 0.20316742672306592174906, 0, 0.23349253653835480876085, 0, 0.24914704581340278500056, 0, 0.24914704581340278500056, 0, 0.23349253653835480876085, 0, 0.20316742672306592174906, 0, 0.16007832854334622633465, 0, 0.106939325995318430960255, 0, 0.047175336386511827194616, 0};	//23rd order Gauss-Legendre weight
__constant__ double w37[] =  {0.00825771143316839575769392, 0.023036084038982232591085, 0.0389152304692994771150896, 0.053697017607756251228889, 0.0672509070508399303049409, 0.079920275333601701493393, 0.0915494682950492105281719, 0.101649732279060277715689, 0.110022604977644072635907, 0.11671205350175682629358, 0.121626303523948383246100, 0.12458416453615607343731, 0.125556893905474335304296, 0.12458416453615607343731, 0.121626303523948383246100, 0.11671205350175682629358, 0.110022604977644072635907, 0.101649732279060277715689, 0.0915494682950492105281719, 0.079920275333601701493393, 0.0672509070508399303049409, 0.053697017607756251228889, 0.0389152304692994771150896, 0.023036084038982232591085, 0.00825771143316839575769392};	//37th order Gauss-Kronrod weight
__global__ void k0_omega_Fermi_37(double* Par_globe, double* omega, double* q, double* Fermi)	//Energy and Fermi function for 37th order
{
	__shared__ double x[25];
	__shared__ double Par[10];
	if(threadIdx.x < 10)
		Par[threadIdx.x] = Par_globe[threadIdx.x];

	x[threadIdx.x] = (Par[8]+Par[7]+Disp37[threadIdx.x]*(Par[8]-Par[7]))/2.;	//abscisca
	omega[threadIdx.x] = sqrt(Par[4]+sq(Par[3]))/2.+x[threadIdx.x];
	omega[threadIdx.x+25] = sqrt(Par[4]+sq(Par[3]))/2.-x[threadIdx.x];
	q[threadIdx.x] = sqrt(sq(Par[3])/4.+sq(Par[5])+Par[3]*Par[5]*cos(Par[6]));
	q[threadIdx.x+25] = sqrt(sq(Par[3])/4.+sq(Par[5])-Par[3]*Par[5]*cos(Par[6]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x] = 0;
		else
			Fermi[threadIdx.x] = 1;
	}
	else
		Fermi[threadIdx.x] = 1./(1.+exp(omega[threadIdx.x]/Par[9]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x+25] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x+25] = 0;
		else
			Fermi[threadIdx.x+25] = 1;
	}
	else
		Fermi[threadIdx.x+25] = 1./(1.+exp(omega[threadIdx.x+25]/Par[9]));
}

__global__ void k0_Vaccum_ImSelf_37(double* q, double* omega_globe, double* ImSelf)			//ImSelf for Vacuum and 37th order
{
	__shared__ double omega[50];
	__shared__ double k[50];
	omega[threadIdx.x] = omega_globe[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA;
	else
		ImSelf[threadIdx.x] = 0;
}

__global__ void k0_Vaccum_ReSelf_37(double* q, double* omega_globe, double* ReSelf)			//ReSelf for Vacuum and 37th order
{
	ReSelf[threadIdx.x] = 0;
}

__global__ void k0_194_ImSelf_37(double* Par, double* q, double* omega_globe, double* ImSelf)	//ImSelf for T=194 MeV and 37th order
{
	__shared__ double omega[50];
	__shared__ double k[50];
	__shared__ double ImSigma[5022];
	__shared__ static double M;
	__shared__ static double omega0[2];		//Location of central peak
	__shared__ static double Sigma[2];		//Amplitude of energy dependance
	__shared__ static double a[2], b[2];	//Slope of exponential decrease to left and right
	__shared__ static double knee[2];		//Interval to change from left to right side of peak
	__shared__ static double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	__shared__ static double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	omega[threadIdx.x] = omega[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[25] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 25))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/25];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/25] = .569969/sqrt(sq(k[threadIdx.x/25])+sq(1.75236))+.0187484;
		a[threadIdx.x/25] = 12.5349/(sq(k[threadIdx.x/25])+sq(1.63711))+5.026;
		b[threadIdx.x/25] = -291.579/(sq(k[threadIdx.x/25]+15.2519)+sq(.0614821))+3.36681;
		omega0[threadIdx.x/25] = sqrt(sq(1.51443+Shift)+sq(k[threadIdx.x/25]))+.232841;
		knee[threadIdx.x/25] = 3.78956*pow(k[threadIdx.x/25]+1., (double)-.530289)+.305*(tanh((k[threadIdx.x/25]-48.4)/11.1111)+1);
	}

	if((omega[threadIdx.x]-omega0[threadIdx.x/25]+knee[threadIdx.x/25]*(b[threadIdx.x/25]-a[threadIdx.x/25])/(sqrt(a[threadIdx.x/25]*b[threadIdx.x/25])*(a[threadIdx.x/25]+b[threadIdx.x/25])))/knee[threadIdx.x/25] < -4.)
		ImSigma[threadIdx.x] = a[threadIdx.x/25]*(omega[threadIdx.x]-omega0[threadIdx.x/25]+knee[threadIdx.x/25]/sqrt(a[threadIdx.x/25]*b[threadIdx.x/25]));
	else if((omega[threadIdx.x]-omega0[threadIdx.x/25]+knee[threadIdx.x/25]*(b[threadIdx.x/25]-a[threadIdx.x/25])/(sqrt(a[threadIdx.x/25]*b[threadIdx.x/25])*(a[threadIdx.x/25]+b[threadIdx.x/25])))/knee[threadIdx.x/25] > 4.)
		ImSigma[threadIdx.x] = b[threadIdx.x/25]*(omega0[threadIdx.x/25]-omega[threadIdx.x]+knee[threadIdx.x/25]/sqrt(a[threadIdx.x/25]*b[threadIdx.x/25]));
	else	//Lost of precision having been circumvented, the actual value
		ImSigma[threadIdx.x] = -.5*((a[threadIdx.x/25]-b[threadIdx.x/25])*omega0[threadIdx.x/25]-((a[threadIdx.x/25]+b[threadIdx.x/25])*knee[threadIdx.x/25])/sqrt(a[threadIdx.x/25]*b[threadIdx.x/25]))+(a[threadIdx.x/25]-b[threadIdx.x/25])*omega[threadIdx.x]/2-sqrt(pow(((a[threadIdx.x/25]+b[threadIdx.x/25])/2.)*(omega[threadIdx.x]-omega0[threadIdx.x/25]+((a[threadIdx.x/25]-b[threadIdx.x/25])*knee[threadIdx.x/25])/(sqrt(a[threadIdx.x/25]*b[threadIdx.x/25])*(a[threadIdx.x/25]+b[threadIdx.x/25]))), 2)+pow(knee[threadIdx.x/25], 2));

#ifdef HALF
	ImSigma[threadIdx.x] = -M*Sigma[threadIdx.x/25]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#else
	ImSigma[threadIdx.x] = -2.*M*Sigma[threadIdx.x/25]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#endif

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA+ImSigma[threadIdx.x];
	else
		ImSelf[threadIdx.x] = ImSigma[threadIdx.x];
}

__global__ void k0_194_ReSelf_37(double* Par, double* q, double* omega_globe, double* ReSelf)	//ReSelf for T=194 MeV and 37th order
{
	__shared__ double omega[50];
	__shared__ double k[50];
	__shared__ static double M;
	__shared__ static double Sigma[2];		//Strength
	__shared__ static double x0[2], x1[2];	//Centrality markers
	__shared__ static double gamma[2];		//Width
	__shared__ static double Shift, M_T;
	__shared__ static double k_old[2];		//Note on validity of k

	omega[threadIdx.x] = omega[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[25] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 25))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/25];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/25] = .212571/sqrt(sq(k[threadIdx.x])+sq(1.17821))+.00762638;
		x0[threadIdx.x/25] = sqrt(sq(k[threadIdx.x])+sq(1.57536+Shift))+.259147;
		x1[threadIdx.x/25] = sqrt(sq(k[threadIdx.x])+sq(1.50194+Shift))+.222526;
		gamma[threadIdx.x/25] = .336699/sqrt(sq(k[threadIdx.x])+sq(1.87956))+.0651449;
	}

#ifdef HALF
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/25]*(omega[threadIdx.x]-x0[threadIdx.x/25])/(pow(omega[threadIdx.x]-x1[threadIdx.x/25], 2)+gamma[threadIdx.x/25])/2.;
#else
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/25]*(omega[threadIdx.x]-x0[threadIdx.x/25])/(pow(omega[threadIdx.x]-x1[threadIdx.x/25], 2)+gamma[threadIdx.x/25]);
#endif
}

__global__ void k0_Ordinate_37(double* Par, double* omega, double* q, double* fermi, double* ImSelf_globe, double* ReSelf, double* Ordinate)	//Ordinate for 37th order
{
	__shared__ double ImSelf[25][2];
	__shared__ double M[25];
	
	M[threadIdx.x] = Par[2];
	ImSelf[threadIdx.x][0] = ImSelf_globe[threadIdx.x];
	ImSelf[threadIdx.x][1] = ImSelf_globe[threadIdx.x+25];

	Ordinate[threadIdx.x] = -((4.*ImSelf[threadIdx.x][0]*ImSelf[threadIdx.x][1]*sq(M[threadIdx.x])*(1.-fermi[threadIdx.x]-fermi[threadIdx.x+25]))/((sq(sq(omega[threadIdx.x])-sq(q[threadIdx.x])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x])+sq(ImSelf[threadIdx.x][0]))*(sq(sq(omega[threadIdx.x+25])-sq(q[threadIdx.x+25])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x+25])+sq(ImSelf[threadIdx.x][1]))));
}

__global__ void k0_Reduce_37(double* Ordinate_globe, double* Answer)	//Reduce 37th order
{
	__shared__ double lower_order[25];
	__shared__ double higher_order[25];
	__shared__ double Ordinate[25];

	Ordinate[threadIdx.x] = Ordinate_globe[threadIdx.x];
	lower_order[threadIdx.x] = Ordinate[threadIdx.x]*w9[threadIdx.x];	//store the lower order value
	higher_order[threadIdx.x] = Ordinate[threadIdx.x]*w16[threadIdx.x];	//store the higher order value

	for(int i = 0; i < 5; i++)	//reduce
	{
		if(threadIdx.x+pow(2,i) < 25)
		{
			lower_order[threadIdx.x] += lower_order[int(threadIdx.x+pow(2,i))];
			higher_order[threadIdx.x] += higher_order[int(threadIdx.x+pow(2,i))];
		}
	}

	if(threadIdx.x == 0)
	{
		Answer[0] = lower_order[0];
		Answer[1] = higher_order[0];
	}
}

//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration, 65 points
__constant__ double Disp97[] = {-0.99954590212436447863561, -0.99726386184948156354498, -0.99262803526297191268579, -0.98561151154526833540018, -0.97631028361466380719767, -0.96476225558750643077381, -0.95095468484866118538988, -0.93490607593773968917092, -0.91667726665136432427535, -0.89632115576605212396531, -0.87386976894531060612966, -0.84936761373256997013369, -0.82288295013605132164827, -0.79448379596794240696310, -0.76422825199780370415066, -0.73218211874028968038743, -0.69842655779521049288477, -0.66304426693021520097512, -0.62611293770182399782024, -0.58771575724076232904075, -0.54794631419915247868094, -0.50689990893222939002375, -0.46466930848199221775618, -0.42135127613063534536412, -0.37704942115412110544534, -0.33186860228212764977992, -0.28591245858945975941661, -0.23928736225213707454460, -0.19210360898314249727164, -0.14447196158279649348519, -0.09650269687689436580083, -0.04830766568773831623481, 0, 0.04830766568773831623481, 0.09650269687689436580083, 0.14447196158279649348519, 0.19210360898314249727164, 0.23928736225213707454460, 0.28591245858945975941661, 0.33186860228212764977992, 0.37704942115412110544534, 0.42135127613063534536412, 0.46466930848199221775618, 0.50689990893222939002375, 0.54794631419915247868094, 0.58771575724076232904075, 0.62611293770182399782024, 0.66304426693021520097512, 0.69842655779521049288477, 0.73218211874028968038743, 0.76422825199780370415066, 0.79448379596794240696310, 0.82288295013605132164827, 0.84936761373256997013369, 0.87386976894531060612966, 0.89632115576605212396531, 0.91667726665136432427535, 0.93490607593773968917092, 0.95095468484866118538988, 0.96476225558750643077381, 0.97631028361466380719767, 0.98561151154526833540018, 0.99262803526297191268579, 0.99726386184948156354498, 0.99954590212436447863561};	//Displacement from center
__constant__ double w63[] = {0, 0.0070186100094700966004071, 0, 0.0162743947309056706051706, 0, 0.025392065309262059455753, 0, 0.034273862913021433102688, 0, 0.042835898022226680656879, 0, 0.050998059262376176196163, 0, 0.058684093478535547145284, 0, 0.065822222776361846837650, 0, 0.072345794108848506225399, 0, 0.078193895787070306471741, 0, 0.083311924226946755222199, 0, 0.087652093004403811142771, 0, 0.091173878695763884712869, 0, 0.093844399080804565639180, 0, 0.09563872007927485941908, 0, 0.09654008851472780056676, 0, 0.09654008851472780056676, 0, 0.09563872007927485941908, 0, 0.093844399080804565639180, 0, 0.091173878695763884712869, 0, 0.087652093004403811142771, 0, 0.083311924226946755222199, 0, 0.078193895787070306471741, 0, 0.072345794108848506225399, 0, 0.065822222776361846837650, 0, 0.058684093478535547145284, 0, 0.050998059262376176196163, 0, 0.042835898022226680656879, 0, 0.034273862913021433102688, 0, 0.025392065309262059455753, 0, 0.0162743947309056706051706, 0, 0.0070186100094700966004071, 0};	//63rd order Gauss-Legendre weight
__constant__ double w97[] = {0.00122336081795147180029304, 0.0034268187757723709355746, 0.00584173707916669330394798, 0.0081725040385316684143438, 0.0104239873988068188280343, 0.012676054806654402859369, 0.0149361036060860273850968, 0.017149805209784253256086, 0.0192987714303268112944037, 0.021408913184821915955778, 0.0234866596721633245920879, 0.025505695480894652814529, 0.0274520984222104037831477, 0.029336956689620661368616, 0.0311633255619737371711558, 0.032915077643903600263296, 0.0345821227447330341307264, 0.036169769475642299860958, 0.0376791306456133985148960, 0.039099420133306611207482, 0.0404234923703730966723493, 0.041654019985643051398296, 0.0427911155964467469336549, 0.043827544030139749046816, 0.0447586387497669372951992, 0.045585826564547070280575, 0.0463087567380257132403813, 0.046922968281703611103481, 0.0474260618738823823628799, 0.047818908736988472212264, 0.0481009691854577469278465, 0.04827019307577738559871, 0.0483263839865677583754454, 0.04827019307577738559871, 0.0481009691854577469278465, 0.047818908736988472212264, 0.0474260618738823823628799, 0.046922968281703611103481, 0.0463087567380257132403813, 0.045585826564547070280575, 0.0447586387497669372951992, 0.043827544030139749046816, 0.0427911155964467469336549, 0.041654019985643051398296, 0.0404234923703730966723493, 0.039099420133306611207482, 0.0376791306456133985148960, 0.036169769475642299860958, 0.0345821227447330341307264, 0.032915077643903600263296, 0.0311633255619737371711558, 0.029336956689620661368616, 0.0274520984222104037831477, 0.025505695480894652814529, 0.0234866596721633245920879, 0.021408913184821915955778, 0.0192987714303268112944037, 0.017149805209784253256086, 0.0149361036060860273850968, 0.012676054806654402859369, 0.0104239873988068188280343, 0.0081725040385316684143438, 0.00584173707916669330394798, 0.0034268187757723709355746, 0.00122336081795147180029304};	//97th order Gauss-Kronrod weight
__global__ void k0_omega_Fermi_97(double* Par_globe, double* omega, double* q, double* Fermi)	//Energy and Fermi function for 97th order
{
	__shared__ double x[65];
	__shared__ double Par[10];
	if(threadIdx.x < 10)
		Par[threadIdx.x] = Par_globe[threadIdx.x];

	__syncthreads();

	x[threadIdx.x] = (Par[8]+Par[7]+Disp97[threadIdx.x]*(Par[8]-Par[7]))/2.;	//abscisca
	omega[threadIdx.x] = sqrt(Par[4]+sq(Par[3]))/2.+x[threadIdx.x];
	omega[threadIdx.x+65] = sqrt(Par[4]+sq(Par[3]))/2.-x[threadIdx.x];
	q[threadIdx.x] = sqrt(sq(Par[3])/4.+sq(Par[5])+Par[3]*Par[5]*cos(Par[6]));
	q[threadIdx.x+65] = sqrt(sq(Par[3])/4.+sq(Par[5])-Par[3]*Par[5]*cos(Par[6]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x] = 0;
		else
			Fermi[threadIdx.x] = 1;
	}
	else
		Fermi[threadIdx.x] = 1./(1.+exp(omega[threadIdx.x]/Par[9]));

	if(Par[9] == 0)
	{
		if(omega[threadIdx.x+65] >= 0)	//Fermi factor for vacuum
			Fermi[threadIdx.x+65] = 0;
		else
			Fermi[threadIdx.x+65] = 1;
	}
	else
		Fermi[threadIdx.x+65] = 1./(1.+exp(omega[threadIdx.x+65]/Par[9]));
}

__global__ void k0_Vaccum_ImSelf_97(double* q, double* omega_globe, double* ImSelf)			//ImSelf for Vacuum and 97th order
{
	__shared__ double omega[130];
	__shared__ double k[130];
	omega[threadIdx.x] = omega_globe[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA;
	else
		ImSelf[threadIdx.x] = 0;
}

__global__ void k0_Vaccum_ReSelf_97(double* q, double* omega_globe, double* ReSelf)			//ReSelf for Vacuum and 97th order
{
	ReSelf[threadIdx.x] = 0;
}

__global__ void k0_194_ImSelf_97(double* Par, double* q, double* omega_globe, double* ImSelf)	//ImSelf for T=194 MeV and 97th order
{
	__shared__ double omega[130];
	__shared__ double k[130];
	__shared__ double ImSigma[130];
	__shared__ static double M;
	__shared__ static double omega0[2];		//Location of central peak
	__shared__ static double Sigma[2];		//Amplitude of energy dependance
	__shared__ static double a[2], b[2];	//Slope of exponential decrease to left and right
	__shared__ static double knee[2];		//Interval to change from left to right side of peak
	__shared__ static double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	__shared__ static double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	omega[threadIdx.x] = omega[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[65] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 65))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/65];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/65] = .569969/sqrt(sq(k[threadIdx.x/65])+sq(1.75236))+.0187484;
		a[threadIdx.x/65] = 12.5349/(sq(k[threadIdx.x/65])+sq(1.63711))+5.026;
		b[threadIdx.x/65] = -291.579/(sq(k[threadIdx.x/65]+15.2519)+sq(.0614821))+3.36681;
		omega0[threadIdx.x/65] = sqrt(sq(1.51443+Shift)+sq(k[threadIdx.x/65]))+.232841;
		knee[threadIdx.x/65] = 3.78956*pow(k[threadIdx.x/65]+1., (double)-.530289)+.305*(tanh((k[threadIdx.x/65]-48.4)/11.1111)+1);
	}
	__syncthreads();

	if((omega[threadIdx.x]-omega0[threadIdx.x/65]+knee[threadIdx.x/65]*(b[threadIdx.x/65]-a[threadIdx.x/65])/(sqrt(a[threadIdx.x/65]*b[threadIdx.x/65])*(a[threadIdx.x/65]+b[threadIdx.x/65])))/knee[threadIdx.x/65] < -4.)
		ImSigma[threadIdx.x] = a[threadIdx.x/65]*(omega[threadIdx.x]-omega0[threadIdx.x/65]+knee[threadIdx.x/65]/sqrt(a[threadIdx.x/65]*b[threadIdx.x/65]));
	else if((omega[threadIdx.x]-omega0[threadIdx.x/65]+knee[threadIdx.x/65]*(b[threadIdx.x/65]-a[threadIdx.x/65])/(sqrt(a[threadIdx.x/65]*b[threadIdx.x/65])*(a[threadIdx.x/65]+b[threadIdx.x/65])))/knee[threadIdx.x/65] > 4.)
		ImSigma[threadIdx.x] = b[threadIdx.x/65]*(omega0[threadIdx.x/65]-omega[threadIdx.x]+knee[threadIdx.x/65]/sqrt(a[threadIdx.x/65]*b[threadIdx.x/65]));
	else	//Lost of precision having been circumvented, the actual value
		ImSigma[threadIdx.x] = -.5*((a[threadIdx.x/65]-b[threadIdx.x/65])*omega0[threadIdx.x/65]-((a[threadIdx.x/65]+b[threadIdx.x/65])*knee[threadIdx.x/65])/sqrt(a[threadIdx.x/65]*b[threadIdx.x/65]))+(a[threadIdx.x/65]-b[threadIdx.x/65])*omega[threadIdx.x]/2-sqrt(pow(((a[threadIdx.x/65]+b[threadIdx.x/65])/2.)*(omega[threadIdx.x]-omega0[threadIdx.x/65]+((a[threadIdx.x/65]-b[threadIdx.x/65])*knee[threadIdx.x/65])/(sqrt(a[threadIdx.x/65]*b[threadIdx.x/65])*(a[threadIdx.x/65]+b[threadIdx.x/65]))), 2)+pow(knee[threadIdx.x/65], 2));

#ifdef HALF
	ImSigma[threadIdx.x] = -M*Sigma[threadIdx.x/65]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#else
	ImSigma[threadIdx.x] = -2.*M*Sigma[threadIdx.x/65]*exp(ImSigma[threadIdx.x]);	//ImSigma from the in-medium
#endif

	if(sq(omega[threadIdx.x])>=sq(k[threadIdx.x]) && omega[threadIdx.x] >= 0)	//Vacuum width
		ImSelf[threadIdx.x] = sqrt(sq(omega[threadIdx.x])-sq(k[threadIdx.x]))*GAMMA+ImSigma[threadIdx.x];
	else
		ImSelf[threadIdx.x] = ImSigma[threadIdx.x];
}

__global__ void k0_194_ReSelf_97(double* Par, double* q, double* omega_globe, double* ReSelf)	//ReSelf for T=194 MeV and 97th order
{
	__shared__ double omega[130];
	__shared__ double k[130];
	__shared__ static double M;
	__shared__ static double Sigma[2];		//Strength
	__shared__ static double x0[2], x1[2];	//Centrality markers
	__shared__ static double gamma[2];		//Width
	__shared__ static double Shift, M_T;
	__shared__ static double k_old[2];		//Note on validity of k

	omega[threadIdx.x] = omega[threadIdx.x];
	k[threadIdx.x] = q[threadIdx.x];

	if((k[0] != k_old[0] || k[65] != k_old[1]) && (threadIdx.x == 0 || threadIdx.x == 65))
	{
		k_old[threadIdx.x%2] = k[threadIdx.x/65];
		if(threadIdx.x == 0)
		{
			M_T = 1.84184;
			M = Par[2];
			Shift = M-M_T;
		}
		Sigma[threadIdx.x/65] = .212571/sqrt(sq(k[threadIdx.x])+sq(1.17821))+.00762638;
		x0[threadIdx.x/65] = sqrt(sq(k[threadIdx.x])+sq(1.57536+Shift))+.259147;
		x1[threadIdx.x/65] = sqrt(sq(k[threadIdx.x])+sq(1.50194+Shift))+.222526;
		gamma[threadIdx.x/65] = .336699/sqrt(sq(k[threadIdx.x])+sq(1.87956))+.0651449;
	}
	__syncthreads();

#ifdef HALF
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/65]*(omega[threadIdx.x]-x0[threadIdx.x/65])/(pow(omega[threadIdx.x]-x1[threadIdx.x/65], 2)+gamma[threadIdx.x/65])/2.;
#else
	ReSelf[threadIdx.x] = Sigma[threadIdx.x/65]*(omega[threadIdx.x]-x0[threadIdx.x/65])/(pow(omega[threadIdx.x]-x1[threadIdx.x/65], 2)+gamma[threadIdx.x/65]);
#endif
}

__global__ void k0_Ordinate_97(double* Par, double* omega, double* q, double* fermi, double* ImSelf_globe, double* ReSelf, double* Ordinate)	//Ordinate for 97th order
{
	__shared__ double ImSelf[65][2];
	__shared__ double M[65];
	
	M[threadIdx.x] = Par[2];
	ImSelf[threadIdx.x][0] = ImSelf_globe[threadIdx.x];
	ImSelf[threadIdx.x][1] = ImSelf_globe[threadIdx.x+65];

	Ordinate[threadIdx.x] = -((4.*ImSelf[threadIdx.x][0]*ImSelf[threadIdx.x][1]*sq(M[threadIdx.x])*(1.-fermi[threadIdx.x]-fermi[threadIdx.x+65]))/((sq(sq(omega[threadIdx.x])-sq(q[threadIdx.x])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x])+sq(ImSelf[threadIdx.x][0]))*(sq(sq(omega[threadIdx.x+65])-sq(q[threadIdx.x+65])-sq(M[threadIdx.x])-2.*M[threadIdx.x]*ReSelf[threadIdx.x+65])+sq(ImSelf[threadIdx.x][1]))));
}

__global__ void k0_Reduce_97(double* Ordinate_globe, double* Answer)	//Reduce 97th order
{
	__shared__ double lower_order[65];
	__shared__ double higher_order[65];
	__shared__ double Ordinate[65];

	Ordinate[threadIdx.x] = Ordinate_globe[threadIdx.x];
	lower_order[threadIdx.x] = Ordinate[threadIdx.x]*w9[threadIdx.x];	//store the lower order value
	higher_order[threadIdx.x] = Ordinate[threadIdx.x]*w16[threadIdx.x];	//store the higher order value

	for(int i = 0; i < 7; i++)	//reduce
	{
		if(threadIdx.x+pow(2,i) < 65)
		{
			lower_order[threadIdx.x] += lower_order[int(threadIdx.x+pow(2,i))];
			higher_order[threadIdx.x] += higher_order[int(threadIdx.x+pow(2,i))];
		}
	}

	if(threadIdx.x == 0)
	{
		Answer[0] = lower_order[0];
		Answer[1] = higher_order[0];
	}
}

void Characterize_k0_Int(double Par[], int Temp, double k, double theta, double zero[10], double gamma[10], int &Poles)
{
	double Lower;	//Limits of integration in k0_Int, vacuum limits are much smaller
	double holder;
	int i, j;

	if(true)//Temp != 0)
		Lower = -sqrt(Par[4]+sq(Par[3]))/2.;	//Integrate from -E/2 to E/2
	else
		Lower = Energy(0, Par[3]/2., k, theta)-sqrt(Par[4]+sq(Par[3]))/2.;

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
	gamma[0] = k-sqrt(abs(sq(k)-pow(Par[1], 2)/2));
	gamma[1] = k-sqrt(abs(sq(k)-pow(Par[1], 2)/2));
#elif VERSION == 42
	zero[0] = .5*sqrt(complex<double>(4.*sq(k), pow(Par[1], 2))).real();	//Potential poles, I know exactly where these are at.
	zero[1] = -.5*sqrt(complex<double>(4.*sq(k), pow(Par[1], 2))).real();
	gamma[0] = abs(.5*sqrt(complex<double>(4.*sq(k), pow(Par[1], 2))).imag());
	gamma[1] = abs(-.5*sqrt(complex<double>(4.*sq(k), pow(Par[1], 2))).imag());
#endif

	zero[2] = .5*(sqrt(Par[4]+sq(Par[3]))-real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])-k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));	//Exact vacuum poles
	zero[3] = .5*(sqrt(Par[4]+sq(Par[3]))+real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])-k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	zero[4] = .5*(-sqrt(Par[4]+sq(Par[3]))-real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])+k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	zero[5] = .5*(-sqrt(Par[4]+sq(Par[3]))+real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])+k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
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
		gamma[2] = gamma[3] = abs(.5*imag(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])-k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
		gamma[4] = gamma[5] = abs(.5*imag(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])+k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))));
	}

	gamma[6] = ImSelf_Energy(Par[2], sqrt(Par[4]+sq(Par[3]))/2.+zero[6], Energy(0, Par[3], k, theta), Temp)+ImSelf_Energy(Par[2], sqrt(Par[4]+sq(Par[3]))/2.-zero[6], Energy(0, Par[3], -k, theta), Temp)-GAMMA;
	if(!isnan(zero[7]))
		gamma[7] = ImSelf_Energy(Par[2], sqrt(Par[4]+sq(Par[3]))/2.+zero[7], Energy(0, Par[3], k, theta), Temp)+ImSelf_Energy(Par[2], sqrt(Par[4]+sq(Par[3]))/2.-zero[7], Energy(0, Par[3], -k, theta), Temp);

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

double Newton_Method_k0(double k0, double Par[], double k, double theta, int Temp, double (*Folding)(double[], double, double, double, int))	//Newton's method for finding poles of f by looking for zeros of 1/f, much more stable to the point of absolute confidence
{
	double new_k0;
	double h;	//Finite difference
	double Exit;
	int i = 0;
	double Danger[] = {.5*(sqrt(Par[4]+sq(Par[3]))+real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])-k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4)))))), .5*(-sqrt(Par[4]+sq(Par[3]))-real(sqrt(complex<double>(4.*(sq(k)+sq(Par[2])+k*Par[3]*cos(theta))+sq(Par[3])-2.*pow(GAMMA, 2), 2.*sqrt(4.*pow(Par[2]*GAMMA, 2)-pow(GAMMA, 4))))))};

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

double omega_Width(double zero, double Par[], double k, double theta, int Temp, double (*Folding)(double[], double, double, double, int))	//Breit-Wigner width of the peak
{
	return(sqrt(abs(2e-10*Folding(Par, zero, k, theta, Temp)/(Folding(Par, zero-1e-5, k, theta, Temp)-2.*Folding(Par, zero, k, theta, Temp)+Folding(Par, zero+1e-5, k, theta, Temp)))));
}

void Characterize_Dispersion(double Par[], int Temp, double k0, double k, double theta, double zero[], double gamma[], int &Poles)
{
	zero[0] = 4.*(sq(k)+pow(k0, 2)+sq(Par[2])+k*Par[3]*cos(theta)-sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)+4.*k*pow(k0, 2)*Par[3]*cos(theta)));
	zero[1] = 4.*(sq(k)+pow(k0, 2)+sq(Par[2])+k*Par[3]*cos(theta)+sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)+4.*k*pow(k0, 2)*Par[3]*cos(theta))); //Both of the possible on-shell s using positive k^mu
	zero[2] = 4.*(sq(k)+pow(k0, 2)+sq(Par[2])-k*Par[3]*cos(theta)-sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)-4.*k*pow(k0, 2)*Par[3]*cos(theta)));
	zero[3] = 4.*(sq(k)+pow(k0, 2)+sq(Par[2])-k*Par[3]*cos(theta)+sqrt(pow(2.*k*k0, 2)+pow(2.*k0*Par[2], 2)+pow(k0*Par[3], 2)-4.*k*pow(k0, 2)*Par[3]*cos(theta))); //Both of the possible on-shell s using negative k^mu
	zero[4] = .5*(pow(2.*k, 2)+pow(2.*Par[2], 2)-sq(Par[3])+sqrt(pow(sq(Par[3])-pow(2.*k, 2)-pow(2.*Par[2], 2), 2)+pow(2.*Par[3], 2)*(pow(2.*Par[2], 2)+pow(2.*k*sin(theta), 2))));	//Intersection of two on-shells in terms of k, M, P, and theta

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

double sp_Width(double Par[], double k0, double k, double theta, int Temp, double (*Folding)(double[], double, double, double, int))	//Breit-Wigner width of the peak
{
	double y[3];
	y[1] = Folding(Par, k0, k, theta, Temp);
	Par[4] -= 1e-5;
	y[0] = Folding(Par, k0, k, theta, Temp);
	Par[4] += 2e-5;
	y[2] = Folding(Par, k0, k, theta, Temp);
	return(sqrt(abs(2e-10*y[1]/(y[0]-2.*y[1]+y[2]))));
}

__device__ __host__ double sq(double a)
{
	return(a*a);
}

void ImSelf_Energy(double M, double omega[], double k[], int Temp, double Results[])	//Single quark self energy for both quarks
{
	static double omega0[2];		//Location of central peak
	static double Sigma[2];		//Amplitude of energy dependance
	static double a[2], b[2];	//Slope of exponential decrease to left and right
	static double knee[2];		//Interval to change from left to right side of peak
	static double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	static double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	if(sq(omega[0])>=sq(k[0]) && omega[0] >= 0)	//Vacuum width
		Results[0] = sqrt(sq(omega[0])-sq(k[0]))*GAMMA;
	else
		Results[0] = 0;
	if(sq(omega[1])>=sq(k[1]) && omega[1] >= 0)
		Results[1] = sqrt(sq(omega[1])-sq(k[1]))*GAMMA;
	else
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
				Sigma[0] = .569969/sqrt(sq(k[0])+pow(1.75236, 2))+.0187484;
				Sigma[1] = .569969/sqrt(sq(k[1])+pow(1.75236, 2))+.0187484;
				a[0] = 4.689/(sq(k[0])+pow(1.18, 2))+4.59495;
				a[1] = 4.689/(sq(k[1])+pow(1.18, 2))+4.59495;
				b[0] = -70400/(pow(k[0]+20, 2)+pow(130, 2))+6.24;
				b[1] = -70400/(pow(k[1]+20, 2)+pow(130, 2))+6.24;
				omega0[0] = sqrt(pow(1.51443+Shift, 2)+sq(k[0]))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift, 2)+sq(k[1]))+.232841;
				knee[0] = 3.78956*pow(k[0]+1., (double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1., (double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(sq(k[0])+pow(1.75236, 2))+.0187484;
				Sigma[1] = .569969/sqrt(sq(k[1])+pow(1.75236, 2))+.0187484;
				a[0] = 12.5349/(sq(k[0])+pow(1.63711, 2))+5.026;
				a[1] = 12.5349/(sq(k[1])+pow(1.63711, 2))+5.026;
				b[0] = -291.579/(pow(k[0]+15.2519, 2)+pow(.0614821, 2))+3.36681;
				b[1] = -291.579/(pow(k[1]+15.2519, 2)+pow(.0614821, 2))+3.36681;
				omega0[0] = sqrt(pow(1.51443+Shift, 2)+sq(k[0]))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift, 2)+sq(k[1]))+.232841;
				knee[0] = 3.78956*pow(k[0]+1., (double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1., (double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;
			case 2://285MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .625855/sqrt(sq(k[0])+pow(1.8429, 2))+.0249334;
				Sigma[1] = .625855/sqrt(sq(k[1])+pow(1.8429, 2))+.0249334;
				a[0] = 3.3971/(sq(k[0])+pow(1.01744, 2))+3.99561;
				a[1] = 3.3971/(sq(k[1])+pow(1.01744, 2))+3.99561;
				b[0] = -65187.5/(pow(k[0]+3.11711, 2)+pow(101.697, 2))+8.15532;
				b[1] = -65187.5/(pow(k[1]+3.11711, 2)+pow(101.697, 2))+8.15532;
				omega0[0] = sqrt(pow(1.5065+Shift, 2)+sq(k[0]))+.209135;
				omega0[1] = sqrt(pow(1.5065+Shift, 2)+sq(k[1]))+.209135;
				knee[0] = 3.1568*pow(k[0]+1., (double)-.624827)+.197004*(tanh((k[0]-27.1743)/10.0192)+1);
				knee[1] = 3.1568*pow(k[1]+1., (double)-.624827)+.197004*(tanh((k[1]-27.1743)/10.0192)+1);
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .587509/sqrt(sq(k[0])+pow(1.84447, 2))+.0309251;
				Sigma[1] = .587509/sqrt(sq(k[1])+pow(1.84447, 2))+.0309251;
				a[0] = 2.44943/(sq(k[0])+pow(.887313, 2))+3.32859;
				a[1] = 2.44943/(sq(k[1])+pow(.887313, 2))+3.32859;
				b[0] = -4439.38/(pow(k[0]-7.23198, 2)+pow(38.9387, 2))+4.55531;
				b[1] = -4439.38/(pow(k[1]-7.23198, 2)+pow(38.9387, 2))+4.55531;
				omega0[0] = sqrt(pow(1.47725+Shift, 2)+sq(k[0]))+.219181;
				omega0[1] = sqrt(pow(1.47725+Shift, 2)+sq(k[1]))+.219181;
				knee[0] = 3.28564*pow(k[0]+1., (double)-.721321)+.330483*(tanh((k[0]-22.9096)/10.7139)+1);
				knee[1] = 3.28564*pow(k[1]+1., (double)-.721321)+.330483*(tanh((k[1]-22.9096)/10.7139)+1);
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .459303/sqrt(sq(k[0])+pow(1.84321, 2))+.0386564;
				Sigma[1] = .459303/sqrt(sq(k[1])+pow(1.84321, 2))+.0386564;
				a[0] = 1.79149/(sq(k[0])+pow(.764836, 2))+2.66209;
				a[1] = 1.79149/(sq(k[1])+pow(.764836, 2))+2.66209;
				b[0] = -1856.16/(pow(k[0]-8.69519, 2)+pow(26.3551, 2))+3.94631;
				b[1] = -1856.16/(pow(k[1]-8.69519, 2)+pow(26.3551, 2))+3.94631;
				omega0[0] = sqrt(pow(1.45428+Shift, 2)+sq(k[0]))+.197493;
				omega0[1] = sqrt(pow(1.45428+Shift, 2)+sq(k[1]))+.197493;
				knee[0] = 3.06296*pow(k[0]+1., (double)-.917081)+.394833*(tanh((k[0]-19.5932)/12.0494)+1);
				knee[1] = 3.06296*pow(k[1]+1., (double)-.917081)+.394833*(tanh((k[1]-19.5932)/12.0494)+1);
				break;
			default:
				omega0[0]=omega0[1]=1.74727;
				Sigma[0]=Sigma[1]=.344006;
				a[0]=a[1]=9.70298;
				b[0]=b[1]=2.11338;
				knee[0]=knee[1]=3.78966;
		}
	}

	double ImSigma[2];	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
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

double ImSelf_Energy(double M, double omega, double k, int Temp)	//Single quark self energy
{

	double omega0;	//location of central peak
	double Sigma;	//size of energy dependance
	double a, b;	//slope of exponential decrease to left and right
	double knee;	//space to change from left to right side of peak
	double M_T, Shift=0;
	double answer;

	if(pow(omega, 2)>=sq(k) && omega >= 0)
		answer = sqrt(pow(omega, 2)-sq(k))*GAMMA;
	else
		answer = 0;

	if(Temp == 0)
		return(answer);

	switch(Temp)
	{
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(sq(k)+pow(1.75236, 2))+.0187484;
			a = 4.689/(sq(k)+pow(1.18, 2))+4.59495;
			b = -70400/(pow(k+20, 2)+pow(130, 2))+6.24;
			omega0 = sqrt(pow(1.51443+Shift, 2)+sq(k))+.232841;
			knee = 3.78956*pow(k+1., (double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;*/
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(sq(k)+pow(1.75236, 2))+.0187484;
			a = 12.5349/(sq(k)+pow(1.63711, 2))+5.026;
			b = -291.579/(pow(k+15.2519, 2)+pow(.0614821, 2))+3.36681;
			omega0 = sqrt(pow(1.51443+Shift, 2)+sq(k))+.232841;
			knee = 3.78956*pow(k+1., (double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .625855/sqrt(sq(k)+pow(1.8429, 2))+.0249334;
			a = 3.3971/(sq(k)+pow(1.01744, 2))+3.99561;
			b = -65187.5/(pow(k+3.11711, 2)+pow(101.697, 2))+8.15532;
			omega0 = sqrt(pow(1.5065+Shift, 2)+sq(k))+.209135;
			knee = 3.1568*pow(k+1., (double)-.624827)+.197004*(tanh((k-27.1743)/10.0192)+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .587509/sqrt(sq(k)+pow(1.84447, 2))+.0309251;
			a = 2.44943/(sq(k)+pow(.887313, 2))+3.32859;
			b = -4439.38/(pow(k-7.23198, 2)+pow(38.9387, 2))+4.55531;
			omega0 = sqrt(pow(1.47725+Shift, 2)+sq(k))+.219181;
			knee = 3.28564*pow(k+1., (double)-.721321)+.330483*(tanh((k-22.9096)/10.7139)+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .459303/sqrt(sq(k)+pow(1.84321, 2))+.0386564;
			a = 1.79149/(sq(k)+pow(.764836, 2))+2.66209;
			b = -1856.16/(pow(k-8.69519, 2)+pow(26.3551, 2))+3.94631;
			omega0 = sqrt(pow(1.45428+Shift, 2)+sq(k))+.197493;
			knee = 3.06296*pow(k+1., (double)-.917081)+.394833*(tanh((k-19.5932)/12.0494)+1);
			break;
		case 5://40MeV
			M_T = 1.8;
			Shift = M-M_T;
			Sigma = .00386564;
			a = 6.2;
			b = 2.8;
			omega0 = sqrt(pow(1.53+Shift, 2)+sq(k));
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

	double ImSigma;	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
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

void ReSelf_Energy(double M, double omega[], double k[], int Temp, double Results[])	//Single quark self energy
{
	static double Sigma[2];		//Strength
	static double x0[2], x1[2];	//Centrality markers
	static double gamma[2];		//Width
	static double Shift, M_T;
	static double k_old[2];		//Note on validity of k

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
				Sigma[0] = .257498/sqrt(sq(k[0])+pow(1.33201, 2))+.00762638;
				Sigma[1] = .257498/sqrt(sq(k[1])+pow(1.33201, 2))+.00762638;
				x0[0] = sqrt(sq(k[0])+pow(1.54778+Shift, 2))+.276509;
				x0[1] = sqrt(sq(k[1])+pow(1.54778+Shift, 2))+.276509;
				x1[0] = sqrt(sq(k[0])+pow(1.49799+Shift, 2))+.246719;
				x1[1] = sqrt(sq(k[1])+pow(1.49799+Shift, 2))+.246719;
				gamma[0] = .658734/sqrt(sq(k[0])+pow(3.35217, 2))+.0815109;
				gamma[1] = .658734/sqrt(sq(k[1])+pow(3.35217, 2))+.0815109;
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571/sqrt(sq(k[0])+pow(1.17821, 2))+.00762638;
				Sigma[1] = .212571/sqrt(sq(k[1])+pow(1.17821, 2))+.00762638;
				x0[0] = sqrt(sq(k[0])+pow(1.57536+Shift, 2))+.259147;
				x0[1] = sqrt(sq(k[1])+pow(1.57536+Shift, 2))+.259147;
				x1[0] = sqrt(sq(k[0])+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(sq(k[1])+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .336699/sqrt(sq(k[0])+pow(1.87956, 2))+.0651449;
				gamma[1] = .336699/sqrt(sq(k[1])+pow(1.87956, 2))+.0651449;
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972/sqrt(sq(k[0])+pow(1.41483, 2))+.0101423;
				Sigma[1] = .307972/sqrt(sq(k[1])+pow(1.41483, 2))+.0101423;
				x0[0] = sqrt(sq(k[0])+pow(1.56476+Shift, 2))+.251031;
				x0[1] = sqrt(sq(k[1])+pow(1.56476+Shift, 2))+.251031;
				x1[0] = sqrt(sq(k[0])+pow(1.50194+Shift, 2))+.222526;
				x1[1] = sqrt(sq(k[1])+pow(1.50194+Shift, 2))+.222526;
				gamma[0] = .550628/sqrt(sq(k[0])+pow(2.43968, 2))+.0981269;
				gamma[1] = .550628/sqrt(sq(k[1])+pow(2.43968, 2))+.0981269;
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131/sqrt(sq(k[0])+pow(1.43308, 2))+.0125796;
				Sigma[1] = .339131/sqrt(sq(k[1])+pow(1.43308, 2))+.0125796;
				x0[0] = sqrt(sq(k[0])+pow(1.55034+Shift, 2))+.257788;
				x0[1] = sqrt(sq(k[1])+pow(1.55034+Shift, 2))+.257788;
				x1[0] = sqrt(sq(k[0])+pow(1.46999+Shift, 2))+.231821;
				x1[1] = sqrt(sq(k[1])+pow(1.46999+Shift, 2))+.231821;
				gamma[0] = .615278/sqrt(sq(k[0])+pow(2.22298, 2))+.143376;
				gamma[1] = .615278/sqrt(sq(k[1])+pow(2.22298, 2))+.143376;
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841/sqrt(sq(k[0])+pow(1.42911, 2))+.0157245;
				Sigma[1] = .304841/sqrt(sq(k[1])+pow(1.42911, 2))+.0157245;
				x0[0] = sqrt(sq(k[0])+pow(1.55511+Shift, 2))+.231105;
				x0[1] = sqrt(sq(k[1])+pow(1.55511+Shift, 2))+.231105;
				x1[0] = sqrt(sq(k[0])+pow(1.44714+Shift, 2))+.20956;
				x1[1] = sqrt(sq(k[1])+pow(1.44714+Shift, 2))+.20956;
				gamma[0] = .862629/sqrt(sq(k[0])+pow(2.67193, 2))+.189598;
				gamma[1] = .862629/sqrt(sq(k[1])+pow(2.67193, 2))+.189598;
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

void Self_Energy(double M, double omega[], double k[], int Temp, double ImSelf[], double ReSelf[])	//Single quark self energy for both quarks. This one has both imaginary and real parts. It is a simple Breit-Wigner peak and simplier than the other provisioned version
{
	static double omega0[2];	//location of central peak
	static double Sigma[2];	//size of energy dependance
	static double gamma[2];	//space to change from left to right side of peak
	static double k_old[2];

	if(sq(omega[0])>=sq(k[0]))
		ImSelf[0] = sqrt(sq(omega[0])-sq(k[0]))*GAMMA;
	else
		ImSelf[0] = 0;
	if(sq(omega[1])>=sq(k[1]))
		ImSelf[1] = sqrt(sq(omega[1])-sq(k[1]))*GAMMA;
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
				Sigma[0] = .840172/sqrt(sq(k[0])+pow(1.45603, 2))+.021257;
				Sigma[1] = .840172/sqrt(sq(k[1])+pow(1.45603, 2))+.021257;
				//omega0[0] = sqrt(sq(M)+sq(k[0]));
				//omega0[1] = sqrt(sq(M)+sq(k[1]));
				omega0[0] = sqrt(pow(1.99829, 2)+sq(k[0]));
				omega0[1] = sqrt(pow(1.99829, 2)+sq(k[1]));
				gamma[0] = 1.05035*pow(k[0]+1.3891, (double)-1.3891)+.01;
				gamma[1] = 1.05035*pow(k[1]+1.3891, (double)-1.3891)+.01;
				break;
			case 2://285MeV
				Sigma[0] = 1.05337/sqrt(sq(k[0])+pow(1.50861, 2))+.0282696;
				Sigma[1] = 1.05337/sqrt(sq(k[1])+pow(1.50861, 2))+.0282696;
				//omega0[0] = sqrt(sq(M)+sq(k[0]));
				//omega0[1] = sqrt(sq(M)+sq(k[1]));
				omega0[0] = sqrt(pow(1.97732, 2)+sq(k[0]));
				omega0[1] = sqrt(pow(1.97732, 2)+sq(k[1]));
				gamma[0] = 1.4624*pow(k[0]+2.64, (double)-1.41048)+.01;
				gamma[1] = 1.4624*pow(k[1]+2.64, (double)-1.41048)+.01;
				break;
			case 3://320MeV
				Sigma[0] = 1.14064/sqrt(sq(k[0])+pow(1.54999, 2))+.0350631;
				Sigma[1] = 1.14064/sqrt(sq(k[1])+pow(1.54999, 2))+.0350631;
				//omega0[0] = sqrt(sq(M)+sq(k[0]));
				//omega0[1] = sqrt(sq(M)+sq(k[1]));
				omega0[0] = sqrt(pow(1.96823, 2)+sq(k[0]));
				omega0[1] = sqrt(pow(1.96823, 2)+sq(k[1]));
				gamma[0] = 2.07102*pow(k[0]+3.037, (double)-1.46076)+.01;
				gamma[1] = 2.07102*pow(k[1]+3.037, (double)-1.46076)+.01;
				break;
			case 4://400MeV
				Sigma[0] = 1.06073/sqrt(sq(k[0])+pow(1.64912, 2))+.0438288;
				Sigma[1] = 1.06073/sqrt(sq(k[1])+pow(1.64912, 2))+.0438288;
				//omega0[0] = sqrt(sq(M)+sq(k[0]));
				//omega0[1] = sqrt(sq(M)+sq(k[1]));
				omega0[0] = sqrt(pow(1.93309, 2)+sq(k[0]));
				omega0[1] = sqrt(pow(1.93309, 2)+sq(k[1]));
				gamma[0] = 3.42222*pow(k[0]+3.663, (double)-1.56165)+.01;
				gamma[1] = 3.42222*pow(k[1]+3.663, (double)-1.56165)+.01;
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

double Energy(double M, double P, double k, double theta)	//Single quark energy, can return momentum if M=0
{
	if(sq(M)+sq(P)+sq(k)+2.*P*k*cos(theta) < 0)
		return(0.);
	else
		return(sqrt(sq(M)+sq(P)+sq(k)+2.*P*k*cos(theta)));
}

double Set_Temp(int T)
{
	const double Temps[] = {0, .194, .258, .32, .4, .04, .04};
	return(Temps[T]);
}

double Fermi(double omega, int T)	//Fermi factor
{
	double Temp = Set_Temp(T);

	if(Temp == 0)
	{
		if(omega >= 0)	//Fermi factor for vacuum
			return(0);
		else
			return(1);
	}
	return(1./(1.+exp(omega/Temp)));
}

double Imk0_Integrand(double Par[], double k0, double k, double theta, int Temp)	//Integrand of the folding integral for positive energy
{
	static double q[2];
	static double k_old = -1;
	double omega[2] = {sqrt(Par[4]+sq(Par[3]))/2.+k0, sqrt(Par[4]+sq(Par[3]))/2.-k0};
	double fermi[2] = {Fermi(omega[0], Temp), Fermi(omega[1], Temp)};
	double ImSelf[2];
	double ReSelf[2];

	if(k_old != k)
	{
		k_old = k;
		q[0] = Energy(0, Par[3]/2., k, theta);
		q[1] = Energy(0, Par[3]/2., -k, theta);
	}

	//Self_Energy(Par[2], omega, q, Temp, ImSelf, ReSelf);
	ImSelf_Energy(Par[2], omega, q, Temp, ImSelf);
	ReSelf_Energy(Par[2], omega, q, Temp, ReSelf);

	return(-((4.*ImSelf[0]*ImSelf[1]*sq(Par[2])*(1.-fermi[0]-fermi[1]))/((sq(sq(omega[0])-sq(q[0])-sq(Par[2])-2.*Par[2]*ReSelf[0])+sq(ImSelf[0]))*(sq(sq(omega[1])-sq(q[1])-sq(Par[2])-2.*Par[2]*ReSelf[1])+sq(ImSelf[1])))));
}
