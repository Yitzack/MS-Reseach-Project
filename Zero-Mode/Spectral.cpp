//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<cfloat>
#include<complex>
#include<chrono>
#include"Around.h"
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double Set_Mq(long double, long double, long double);			//Momentum dependence for the quark mass, <number> 0 causes it to be constant
void mergeSort(long double[], int, int);

Around theta_Int(long double[], int, long double(*Rho_Int)(long double, long double, long double, long double, int));
Around theta_Int(long double[], int, long double, long double, int, long double(*Rho_Int)(long double, long double, long double, long double, int));

long double RhoInt1(long double, long double, long double, long double, int);
long double RhoInt2(long double, long double, long double, long double, int);
long double RhoInt3(long double, long double, long double, long double, int);

int main(int argc, char* argv[])
{
	char File[70] = "data/Zero_Mode.";  //Name of the file

	char* Process = argv[1];
	strcat(File, argv[3]);	//Appends the temprature to the file name
	strcat(File, ".");
	strcat(File, Process);	//Appends the process number to the file name

	bool Restart = Restart_Check(File, argv[4], argv[5], argv[6], argv[9], argv[10]);	//True if scrapping the file contents and restarting, only if all args are already in file header

	ofstream TPlot;
	if(Restart)	//If starting from the beginning, overwrite
	{
		TPlot.open(File);
		TPlot << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[9] << " " << argv[10] << endl;
	}
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);// */
	//cout << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[9] << " " << argv[10] << endl;

	int i,j;					//Counters
	int Start, Finish;
	Start = atoi(argv[7]);				//Initial starting point
	Finish = atoi(argv[8]);			//Finish at the point given
	Start = Start_Point(Start, File);		//Go find last point written to file and maybe start from there
	if(Finish < Start && Finish >= atoi(argv[7]))	//Missing point directive issued, go back and get the missed points
		Start = atoi(argv[7]);

	const int iProcess = atoi(argv[1]) % atoi(argv[2]);	//Assigned column(s)
	const int Total = atoi(argv[2]);			//Number of concurent threads
	const int Temp = atoi(argv[3]);			//Temprature enumeration
	long double Table[616][5];				//Table of calculated values
	long double Par[5];					//Parameters to be used in calculation {Coupling constant, potential cutoff, quark mass, P, s}
	Around holder[2];					//Calculated value before distribution to Table
	time_t Start_Time, End_Time;				//Time at the start and end of calculation

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	//cout << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = Start; i <= Finish; i++)
	{
		for(j = iProcess; j < 151; j+=Total)	//Does the subset of j that has been assigned to this process
		{
			if(i <= 208)
			{
				Par[3] = i/10.+j/10.;
				Par[4] = -i*j/50.-j*j/100.;
			}
			else
			{
				Par[3] = i+j/10.-187.2;
				Par[4] = -j/5.*(i-187.2)-j*j/100.;
			}

			Par[2] = atof(argv[6]);

			auto Start_Time = chrono::system_clock::now();
			holder[0] = theta_Int(Par, Temp, RhoInt1);
			holder[1] = theta_Int(Par, Temp, RhoInt2);
			holder[2] = theta_Int(Par, Temp, RhoInt3);
			auto End_Time = chrono::system_clock::now();
			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << holder[0] << " " << holder[1] << " " << holder[2] << " " << chrono::duration_cast<chrono::nanoseconds>(End_Time-Start_Time).count()/1000000000. << endl;
		}
		TPlot << endl;
	}

	return(0);
}

int Start_Point(int Start, char File[70])	//Go through and find largest starting point in file and return it, causes it to repeat last line
{
	ifstream TPlot(File);
	char Line[700];
	int Test;

	TPlot.getline(Line, 700);
	if(!TPlot.is_open())
		return(Start);

	TPlot.getline(Line, 700);
	while(!TPlot.eof())
	{
		Test = atoi(Line);
		if(Test > Start)
			Start = Test;
		TPlot.getline(Line,700);
	}

	return(Start);
}

bool Restart_Check(char File[70], char* g, char* Lambda, char* Mq, char* P0, char* fraction)	//Looks for same input values in fist line of file
{
	ifstream InFile(File);

	if(InFile.is_open() == false)
		return(true);

	double g_File;
	double Lambda_File;
	double Mq_File;
	double P0_File;
	double fraction_File;
	InFile >> g_File;
	InFile >> Lambda_File;
	InFile >> Mq_File;
	InFile >> P0_File;
	InFile >> fraction_File;
	InFile.close();

	if(abs(g_File/atof(g)-1.) < .0001 &&
	   abs(Mq_File/atof(Mq)-1.) < .0001 &&
	   (abs(Lambda_File/atof(Lambda)-1.) < .0001 || Lambda_File-atof(Lambda) < .0001) &&
	   abs(P0_File/atof(P0)-1.) < .0001 &&
	   (abs(fraction_File/atof(fraction)-1.) < .0001 || fraction_File-atof(fraction) < .0001))
		return(false);

	InFile.close();
	return(true);
}

long double Set_Mq(long double Mq0, long double P, long double P0)
{
#ifndef BB
	long double Mqf = 1.85;
#else
	long double Mqf = 5.25;
#endif

	return((Mq0*pow(P0,2)+Mqf*pow(P,2))/(pow(P0,2)+pow(P,2)));
}
//\[Rho]1[s_,P_,M_,T_]:=If[s<2 M^2-Sqrt[2] Sqrt[2 M^4+2 M^2 P^2],NIntegrate[\[Rho]Int1[s,P,M,\[Theta],T],{\[Theta],0,\[Pi]}],NIntegrate[\[Rho]Int1[s,P,M,\[Theta],T],{\[Theta],0,1/2 ArcCos[(2 M^2 P^2+4 M^2 s-s^2)/(2 M^2 P^2)]}]+NIntegrate[\[Rho]Int1[s,P,M,\[Theta],T],{\[Theta],\[Pi]-1/2 ArcCos[(2 M^2 P^2+4 M^2 s-s^2)/(2 M^2 P^2)],\[Pi]}]]
//ArcSin[Sqrt[-s]/P]
//1/2 ArcCos[(2 M^2 P^2+4 M^2 s-s^2)/(2 M^2 P^2)]=arccos(1+(2 s)/P^2-s^2/(2 M^2 P^2))/2
//\[Pi]-1/2 ArcCos[(2 M^2 P^2+4 M^2 s-s^2)/(2 M^2 P^2)],
//\[Pi]-ArcSin[Sqrt[-s]/P]}]]

Around theta_Int(long double Par[], int Temp, long double(*Rho_Int)(long double, long double, long double, long double, int))
{
	if(Par[4] >= 0 || Par[4] <= -pow(Par[3],2)+.001)	//Short cut for s>=0 || s=-P^2 answer is zero
		return(Around(0,0));

	long double x1;
	long double a = 0, b;					//Sub-interval limits of integration
	Around Answer(0,0);				//Answer to be returned
	int i, j;						//Counters

	//List of boundaries between subintervals
	long double Range[27] = {0., .05*M_PI, .1*M_PI, .15*M_PI, .2*M_PI, .25*M_PI, .3*M_PI, .35*M_PI, .4*M_PI, .45*M_PI, .5*M_PI, .55*M_PI, .6*M_PI, .65*M_PI, .7*M_PI, .75*M_PI, .8*M_PI, .85*M_PI, .9*M_PI, .95*M_PI, M_PI, asin(sqrt(-Par[4])/Par[3]), acos(1.+2.*Par[4]/pow(Par[3],2)-pow(Par[4]/(Par[2]*Par[3]),2)/2.)/2., M_PI-acos(1.+2.*Par[4]/pow(Par[3],2)-pow(Par[4]/(Par[2]*Par[3]),2)/2.)/2., M_PI-asin(sqrt(-Par[4])/Par[3]),acos(-(sqrt(Par[4]+pow(Par[3],2))/Par[3])), acos((sqrt(Par[4]+pow(Par[3],2))/Par[3]))};

	for(i = 21; i < 27; i++)	//Replace NaN with pi and send to end of range
		if(isnan(Range[i]))
			Range[i] = M_PI;

	//Put in asending order
	mergeSort(Range, 0, 26);

	for(i = 0; i < 27 && Range[i] < M_PI; i++)	//Count through pre-determined intervals
	{
		b = Range[i];	//Upper edge

		Answer += theta_Int(Par, Temp, a, b, 0, Rho_Int);	//Add the subinterval to total of the integral
		a = b;	//Upper edge becomes lower edge
	}

	return(Answer);
}

Around theta_Int(long double Par[], int Temp, long double a, long double b, int deep, long double(*Rho_Int)(long double, long double, long double, long double, int))
{
	if(abs(a/b-1.) < DBL_EPSILON)	//Shortcut zero width interval and prevents returning nan
		return(Around(0,0));

/*9th order Gauss-Legendre integration/16th order Gauss-Kronrod weight
	long double Disp9[] = {0.2796304131617831934134665, sqrt(5.-2.*sqrt(10./7.))/3., 0.7541667265708492204408172, sqrt(5.+2.*sqrt(10./7.))/3., 0.9840853600948424644961729};	//Displacement from center
	long double w9[] = {128./225., 0., (322.+13.*sqrt(70.))/900., 0., (322.-13.*sqrt(70.))/900., 0.};	//9th order Gauss-Legendre weights
	long double w16[]= {0.2829874178574912132042556, 0.27284980191255892234099326, 0.2410403392286475866999426, 0.18680079655649265746780003, 0.11523331662247339402462685, 0.042582036751081832864509451}; //16th order Gauss-Kronrod weights*/
//1st/2nd order Newton-Coates integration
	long double Disp9[] = {1};	//Displacement from center
	long double w9[] = {0,1.};	//1st order Newton-Coates weights
	long double w16[]= {4./3.,1./3.}; //2nd order Newton-Coates weights*/
#if ORDER == 37	//37th order Gauss-Legendre integration
//23th order Gauss-Legendre/37th order Gauss-Kronrod integration
	long double Disp[] = {0.1252334085114689154724414, 0.2485057483204692762677910, 0.3678314989981801937526915, 0.4813394504781570929359436, 0.5873179542866174472967024, 0.6840598954700558939449291, 0.7699026741943046870368938, 0.8435581241611532447921419, 0.9041172563704748566784659, 0.9505377959431212965490602, 0.9815606342467192506905491, 0.9969339225295954269123502};	//Displacement from center
	long double wl[] = {0, 0.2491470458134027850005624, 0, 0.2334925365383548087608499, 0, 0.2031674267230659217490645, 0, 0.1600783285433462263346525, 0, 0.10693932599531843096025472, 0, 0.04717533638651182719461596, 0};	//23rd order Gauss-Legendre weight
	long double wh[] =  {0.12555689390547433530429613, 0.1245841645361560734373125, 0.12162630352394838324609976, 0.1167120535017568262935807, 0.11002260497764407263590740, 0.10164973227906027771568877, 0.091549468295049210528171940, 0.07992027533360170149339261, 0.067250907050839930304940940, 0.05369701760775625122888916, 0.038915230469299477115089632, 0.02303608403898223259108458, 0.0082577114331683957576939224};	//37th order Gauss-Kronrod weight
#elif ORDER == 97
//63rd order Gauss-Legendre/97th order Gauss-Kronrod integration
	long double Disp[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};	//Displacement from center
	long double wl[] = {0, 0.0965400885147278005667648, 0, 0.0956387200792748594190820, 0, 0.09384439908080456563918024, 0, 0.09117387869576388471286858, 0, 0.08765209300440381114277146, 0, 0.08331192422694675522219907, 0, 0.07819389578707030647174092, 0, 0.07234579410884850622539936, 0, 0.06582222277636184683765006, 0, 0.05868409347853554714528364, 0, 0.05099805926237617619616324, 0, 0.04283589802222668065687865, 0, 0.03427386291302143310268773, 0, 0.02539206530926205945575259, 0, 0.016274394730905670605170562, 0, 0.007018610009470096600407064, 0};	//63rd order Gauss-Legendre weight
	long double wh[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};	//97th order Gauss-Kronrod weight
#endif
	long double x1, x2;		//Abscissa
	long double F[2] = {0,0};	//Sum of ordinate*weights
	Around Answer;			//Answer to be returned
	long double Holder;
	int i, j = 0;			//Counters

	//for(j = 0; j < 5; j++)	//Count through points away from center
	{
		x1 = (b+a-Disp[j]*(b-a))/2.;
		x2 = (b+a+Disp[j]*(b-a))/2.;

		Holder = Rho_Int(Par[4], Par[3], Par[2], x1, Temp);
		F[0] += Holder*w9[j+1];
		F[1] += Holder*w16[j+1];

		Holder = Rho_Int(Par[4], Par[3], Par[2], x2, Temp);
		F[0] += Holder*w9[j+1];
		F[1] += Holder*w16[j+1];
	}
	Holder = Rho_Int(Par[4], Par[3], Par[2], (a+b)/2., Temp);
	F[0] += Holder*w9[0];
	F[1] += Holder*w16[0];

	//Answer = Estimation(F[0], F[1])*Around((b-a)/2.);	//Add the subinterval to total of the integral

	if(abs(F[0]-F[1])*Around(2.)/abs(F[0]+F[1]) > 1 && abs(b-a) > FLT_EPSILON && deep > 4)
		Answer = theta_Int(Par, Temp, a, (a+b)/2., deep+1, Rho_Int) + theta_Int(Par, Temp, (a+b)/2., b, deep+1, Rho_Int);
	else
	{
		F[0] = 0;	//Zero out F for a new round
		F[1] = 0;	//Zero out F for a new round

#if ORDER == 37	//Count through points away from center
		for(j = 0; j < 12; j++)
#elif ORDER == 97
		for(j = 0; j < 32; j++)
#endif
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			Holder = Rho_Int(Par[4], Par[3], Par[2], x1, Temp);
			F[0] += Holder*wl[j+1];
			F[1] += Holder*wh[j+1];

			Holder = Rho_Int(Par[4], Par[3], Par[2], x2, Temp);
			F[0] += Holder*wl[j+1];
			F[1] += Holder*wh[j+1];
		}
		Holder = Rho_Int(Par[4], Par[3], Par[2], (a+b)/2., Temp);
		F[0] += Holder*wl[0];
		F[1] += Holder*wh[0];

		Answer = Around(F[1], abs(F[0]-F[1]))*Around((b-a)/2.);	//Add the subinterval to total of the integral
	}

	return(Answer);
}

long double RhoInt1(long double s, long double P, long double M, long double theta, int Temp)
{
	long double T;
	switch(Temp)
	{
	case 0:
		return(0);
	case 1:
		T = .194;
		break;
	case 2:
		T = .258;
		break;
	case 3:
		T = .320;
		break;
	case 4:
		T = .400;
		break;
	}

	if(2.*pow(M,2)-s-2.*sqrt(pow(M,4)+pow(M*P*sin(theta),2)) < 0.)
		return(0);
	if(-(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4*pow(M,2)*s+pow(s,2)-4*pow(M,2)*pow(P,2)*pow(sin(theta),2))))/(2.*(s+pow(P,2)*pow(sin(theta),2))) < 0)
		return(0);

	return((6.*(1./(1.+exp(sqrt(4.*pow(M,2)+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/(2.*T)))-1./(1.+exp((sqrt(pow(P,2)+s)+sqrt(4.*pow(M,2)+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/2.)/T)))*sin(theta)*(-((P*cos(theta)*(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta)))+(sqrt(pow(P,2)+s)*sqrt(4.*pow(M,2)+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/2.))/(M_PI*abs(((sqrt(2)*(-2.*pow(P,3)*cos(theta)-P*s*cos(theta)+2.*pow(P,3)*pow(cos(theta),3)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/sqrt((-(pow(M,2)*pow(P,4))+3.*pow(P,6)-4.*pow(M*P,2)*s+7.*pow(P,4)*s+10.*pow(P*s,2)+2.*pow(s,3)-6.*pow(P*s*cos(theta),2)-4.*pow(P,2)*(pow(P,4)+pow(M,2)*s+2.*pow(P,2)*s)*cos(2.*theta)+pow(M,2)*pow(P,4)*cos(4.*theta)+pow(P,6)*cos(4.*theta)+pow(P,4)*s*cos(4.*theta)+4.*P*cos(theta)*(-pow(P,2)-s+pow(P,2)*cos(2.*theta))*sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/pow(s+pow(P*sin(theta),2),2))+(-(P*s*cos(theta))-sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/sqrt(4.*pow(M,2)+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/(s+pow(P*sin(theta),2)))*sqrt(pow(M,2)+pow(P,2)+(P*cos(theta)*(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta))+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/(4.*pow(s+pow(P*sin(theta),2),2)))*sqrt(4.*pow(M,2)+pow(P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))));
}

long double RhoInt2(long double s, long double P, long double M, long double theta, int Temp)
{
	long double T;
	switch(Temp)
	{
	case 0:
		return(0);
	case 1:
		T = .194;
		break;
	case 2:
		T = .258;
		break;
	case 3:
		T = .320;
		break;
	case 4:
		T = .400;
		break;
	}

	if(s+pow(P*sin(theta),2) < 0 || 2.*pow(M,2)-s-2.*sqrt(pow(M,4)+pow(M*P*sin(theta),2)) < 0)
		return(0);
	if((-P*s*cos(theta)+sqrt((pow(P,2)+s)*(-4*pow(M,2)*s+pow(s,2)-4*pow(M,2)*pow(P,2)*pow(sin(theta),2))))/(2.*(s+pow(P,2)*pow(sin(theta),2))) < 0)
		return(0);

	return((6.*(1./(1.+exp(sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/(2.*T)))-1./(1.+exp((sqrt(pow(P,2)+s)+sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/2.)/T)))*sin(theta)*(-((P*cos(theta)*(P*s*cos(theta)-sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta)))+(sqrt(pow(P,2)+s)*sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/2.))/(M_PI*abs((-((sqrt(2)*(2.*pow(P,3)*cos(theta)+P*s*cos(theta)-2.*pow(P,3)*pow(cos(theta),3)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/sqrt((-(pow(M,2)*pow(P,4))+3.*pow(P,6)-4.*pow(M*P,2)*s+7.*pow(P,4)*s+10.*pow(P*s,2)+2.*pow(s,3)-6.*pow(P*s*cos(theta),2)-4.*pow(P,2)*(pow(P,4)+pow(M,2)*s+2.*pow(P,2)*s)*cos(2.*theta)+pow(M,2)*pow(P,4)*cos(4.*theta)+pow(P,6)*cos(4.*theta)+pow(P,4)*s*cos(4.*theta)-4.*P*cos(theta)*(-pow(P,2)-s+pow(P,2)*cos(2.*theta))*sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/pow(s+pow(P*sin(theta),2),2)))+(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/(s+pow(P*sin(theta),2)))*sqrt(pow(M,2)+pow(P,2)+(P*cos(theta)*(P*s*cos(theta)-sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta))+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/(4.*pow(s+pow(P*sin(theta),2),2)))*sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))));
}

long double RhoInt3(long double s, long double P, long double M, long double theta, int Temp)
{
	long double T;
	switch(Temp)
	{
	case 0:
		return(0);
	case 1:
		T = .194;
		break;
	case 2:
		T = .258;
		break;
	case 3:
		T = .320;
		break;
	case 4:
		T = .400;
		break;
	}

	if(s+pow(P*sin(theta),2) >= 0 || theta > M_PI/2.)
		return(0);

	return((12.*(-(1./(1.+exp(sqrt(4.*pow(M,2)+pow(P,2)+2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))/(2.*T))))+1./(1.+exp((-2.*sqrt(pow(P,2)+s)+sqrt(4.*pow(M,2)+pow(P,2)+2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))/(2.*T))))*(pow(P,2)+P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))-sqrt(pow(P,2)+s)*sqrt(4.*pow(M,2)+pow(P,2)+2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))))/(M_PI*abs((-4.*P*cos(theta)+4.*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))/(2.*sqrt(4.*pow(M,2)+pow(P,2)-2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))-(2.*(P*cos(theta)+sqrt(-(((4.*pow(M,2)-s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))))/sqrt(4.*pow(M,2)+pow(P,2)+2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2))))*sqrt(4.*pow(M,2)+pow(P,2)-2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))*sqrt(4.*pow(M,2)+pow(P,2)+2.*P*cos(theta)*sqrt(((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))+((-4.*pow(M,2)+s)*(pow(P,2)+s))/(s+pow(P*sin(theta),2)))));
}

void mergeSort(long double List[], int a, int b)
{
	int i, j, k;
	long double Temp[(a+b)/2-a+1];

	if(b-a > 1)	//Divide...
	{
		mergeSort(List, a, (a+b)/2);
		mergeSort(List, (a+b)/2+1, b);
	}

	for(i = 0; i <= (a+b)/2-a; i++)	//Copy out the lower half array in prep for copy over
		Temp[i] = List[i+a];

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
