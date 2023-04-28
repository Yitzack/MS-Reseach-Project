//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<complex>
#include<chrono>
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double Set_Mq(long double, long double, long double);			//Momentum dependence for the quark mass, <number> 0 causes it to be constant

long double RhoInt1(long double, long double, long double, long double, int);
long double RhoInt2(long double, long double, long double, long double, int);

int main(int argc, char* argv[])
{
	char File[70] = "data/Zero_Mode";  //Name of the file

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
	Elements<Around> holder;					//Calculated value before distribution to Table
	time_t Start_Time, End_Time;				//Time at the start and end of calculation

	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	//cout << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = Start; i <= Finish; i++)
	{
		for(j = iProcess+151; j < 625; j+=Total)	//Does the subset of j that has been assigned to this process
		{
			if(j <= 150)
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
			}
			else
			{
				Par[3] = i*.8;
				if(i < 0)
					Par[3] = ((long double)(i%7)/8.-.125-floor((long double)(i)/7.))*.8;
#ifndef BB
				if(j <= 161)
					Par[4] = pow((j-151.)/100.,2);
				else if(j <= 190)
					Par[4] = pow((j-161.)/10.+.1,2);
				else if(j <= 390)
					Par[4] = pow((j-190.)/100.+3.,2);
				else if(j <= 575)
					Par[4] = pow((j-390.)/10.+5.,2);
				else
					Par[4] = 552.25+GaussLa[j-576];
#else
				if(j <= 251)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 451)
					Par[4] = pow((j-251.)/100.+10.,2);
				else if(j <= 566)
					Par[4] = pow((j-451.)/10.+12.,2);
				else
					Par[4] = 552.25+GaussLa[j-567];
#endif
			}

			Par[1] = Set_Lambda(atof(argv[5]), Par[3], atof(argv[9]), atof(argv[10]), Temp);
			Par[0] = -Set_C(atof(argv[4]), Par[3], atof(argv[9]), Par[1], atof(argv[10]));
			Par[2] = atof(argv[6]);

			auto Start_Time = chrono::system_clock::now();
			holder = theta_Int(Par, Temp);
			auto End_Time = chrono::system_clock::now();
			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << *holder[0] << " " << *holder[1] << " " << *holder[2] << " " << *holder[3] << " " << *holder[4] << " " << *holder[5] << " " << chrono::duration_cast<chrono::nanoseconds>(End_Time-Start_Time).count()/1000000000. << endl;
			//cout << i << " " << j << " " << Par[3] << " " << Par[4] << " " << *holder[0] << " " << *holder[1] << " " << *holder[2] << " " << *holder[3] << " " << *holder[4] << " " << *holder[5] << " " << chrono::duration_cast<chrono::nanoseconds>(End_Time-Start_Time).count()/1000000000. << endl;
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
	long double Mqf = 1.8;
#else
	long double Mqf = 5.25;
#endif

	return((Mq0*pow(P0,2)+Mqf*pow(P,2))/(pow(P0,2)+pow(P,2)));
}

Elements<Around> theta_Int(long double Par[], int Temp)
{
	if(Par[3] == 0)	//Short cut for P=0, theta integral is analytic
		return(k_Int(Par, Temp, M_PI/2.)*Around(2./pow(2.*M_PI,2)));
	else if(Par[4]+pow(Par[3],2)<1e-9)
	{
		Around Answer[6] = {0,0,0,0,0,0};
		return(Elements<Around>(Answer,6));
	}

	long double x1;
	long double a = 0, b;					//Sub-interval limits of integration
	long double Boundary_theta[] = {1./17., 0.3, 0.08};	//Extra boundary values
	Elements<Around> Answer(6);				//Answer to be returned
	int i, j;						//Counters
	Answer.null();

	if(Par[4] > 0 && Par[3] > sqrt(Par[4]/2.)) //Where the maximum of the theta integral ought to land. It might only be correct for BbS reduction, but is close enough for all other cases. Only valid for s>0 and P>sqrt(s/2)
		x1 = asin(sqrt(Par[4]/2.)/Par[3]);
	else
		x1 = M_PI/10.;	//If it isn't valid, value is needed anyways to split up the integral

	//Don't get too close to the pole or details might get lost
	if(x1>M_PI/10.)
		x1 = M_PI/10.;

	//List of boundaries between subintervals
	long double Range[] = {x1*Boundary_theta[0], x1*Boundary_theta[1], x1, x1*(2.-Boundary_theta[1]), x1*(2.-Boundary_theta[1])*(1.-Boundary_theta[2])+M_PI/2.*Boundary_theta[2], M_PI/2., asin(sqrt(-Par[4])/Par[3]),0,0};

	//Some kind of intersection, probably between the simultanous on-shell and potential peak, don't rightly remember
	Range[7] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[7] = acos((pow(Range[7],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(Range[7]*Par[3]));
	Range[8] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range[8] = acos((pow(Range[8],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(-Range[8]*Par[3]));

	//Bad data trap for NaN and negative boundaries. These are only ones that can NaN or return negative numbers
	if(isnan(Range[6]) || Range[6] < 0) Range[6] = M_PI;
	if(isnan(Range[7])) Range[7] = M_PI;
	if(isnan(Range[8])) Range[8] = M_PI;

	//Put in asending order
	mergeSort(Range, 0, 8);

	for(i = 0; i < 8 && a < M_PI/2.; i++)	//Count through pre-determined intervals
	{
		b = Range[i];	//Upper edge

		Answer += theta_Int(Par, Temp, a, b, 0);	//Add the subinterval to total of the integral

		a = b;	//Upper edge becomes lower edge
	}

	return(Answer*Around(2./pow(2.*M_PI,2)));
}

Elements<Around> theta_Int(long double Par[], int Temp, long double a, long double b, int deep)
{
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
	Elements<Around> F[2] = {Elements<Around>(6),Elements<Around>(6)};//Sum of ordinate*weights
	Elements<Around> Answer(6);	//Answer to be returned
	Answer.null();
	Elements<Around> Holder;
	int i, j = 0;			//Counters

	F[0].null();	//Zero out F for a new round
	F[1].null();	//Zero out F for a new round

	//for(j = 0; j < 5; j++)	//Count through points away from center
	{
		x1 = (b+a-Disp[j]*(b-a))/2.;
		x2 = (b+a+Disp[j]*(b-a))/2.;

		Holder = k_Int(Par, Temp, x1)*Around(sin(x1));
		F[0] += Holder*Around(w9[j+1]);
		F[1] += Holder*Around(w16[j+1]);

		Holder = k_Int(Par, Temp, x2)*Around(sin(x2));
		F[0] += Holder*Around(w9[j+1]);
		F[1] += Holder*Around(w16[j+1]);
	}
	Holder = k_Int(Par, Temp, (a+b)/2.)*Around(sin((a+b)/2.));
	F[0] += Holder*Around(w9[0]);
	F[1] += Holder*Around(w16[0]);

	//Answer = Estimation(F[0], F[1])*Around((b-a)/2.);	//Add the subinterval to total of the integral

	if(abs(F[0]-F[1])*Around(2.)/abs(F[0]+F[1]) > 1 && abs(b-a) > FLT_EPSILON && deep > 4)
		Answer = theta_Int(Par, Temp, a, (a+b)/2., deep+1) + theta_Int(Par, Temp, (a+b)/2., b, deep+1);
	else
	{
		F[0].null();	//Zero out F for a new round
		F[1].null();	//Zero out F for a new round

#if ORDER == 37	//Count through points away from center
		for(j = 0; j < 12; j++)
#elif ORDER == 97
		for(j = 0; j < 32; j++)
#endif
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;
			x2 = (b+a+Disp[j]*(b-a))/2.;

			Holder = k_Int(Par, Temp, x1)*Around(sin(x1));
			F[0] += Holder*Around(wl[j+1]);
			F[1] += Holder*Around(wh[j+1]);

			Holder = k_Int(Par, Temp, x2)*Around(sin(x2));
			F[0] += Holder*Around(wl[j+1]);
			F[1] += Holder*Around(wh[j+1]);
		}
		Holder = k_Int(Par, Temp, (a+b)/2.)*Around(sin((a+b)/2.));
		F[0] += Holder*Around(wl[0]);
		F[1] += Holder*Around(wh[0]);

		Answer = Estimation(F[0], F[1])*Around((b-a)/2.);	//Add the subinterval to total of the integral
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

	return((6.*(1./(1.+exp(sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/(2.*T)))-1./(1.+exp((sqrt(pow(P,2)+s)+sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))/2.)/T)))*sin(theta)*(-((P*cos(theta)*(P*s*cos(theta)-sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta)))+(sqrt(pow(P,2)+s)*sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/2.))/(M_PI*abs((-((sqrt(2)*(2.*pow(P,3)*cos(theta)+P*s*cos(theta)-2.*pow(P,3)*pow(cos(theta),3)+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/sqrt((-(pow(M,2)*pow(P,4))+3.*pow(P,6)-4.*pow(M*P,2)*s+7.*pow(P,4)*s+10.*pow(P*s,2)+2.*pow(s,3)-6.*pow(P*s*cos(theta),2)-4.*pow(P,2)*(pow(P,4)+pow(M,2)*s+2.*pow(P,2)*s)*cos(2.*theta)+pow(M,2)*pow(P,4)*cos(4.*theta)+pow(P,6)*cos(4.*theta)+pow(P,4)*s*cos(4.*theta)-4.*P*cos(theta)*(-pow(P,2)-s+pow(P,2)*cos(2.*theta))*sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/pow(s+pow(P*sin(theta),2),2)))+(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))))/sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2)))/(s+pow(P*sin(theta),2)))*sqrt(pow(M,2)+pow(P,2)+(P*cos(theta)*(P*s*cos(theta)-sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2)))))/(-pow(P,2)-2.*s+pow(P,2)*cos(2.*theta))+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/(4.*pow(s+pow(P*sin(theta),2),2)))*sqrt(4.*pow(M,2)+pow(-(P*s*cos(theta))+sqrt((pow(P,2)+s)*(-4.*pow(M,2)*s+pow(s,2)-4.*pow(M*P*sin(theta),2))),2)/pow(s+pow(P*sin(theta),2),2))));
}

