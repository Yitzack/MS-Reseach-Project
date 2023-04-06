//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<complex>
#include<chrono>
#include"Spectral.h"
using namespace std;

int Start_Point(int, char[70]);						//Find highest line calculated and returns it, as written causes last line to be recalculated
bool Restart_Check(char[70], char*, char*, char*, char*, char*);		//Checks to see if file header matches input parameters and clears it if not
long double Set_Mq(long double, long double, long double);			//Momentum dependence for the quark mass, <number> 0 causes it to be constant
long double Set_Lambda(long double, long double, long double, long double, int);//Momentum dependence for the potential cutoff, <number> 0 causes it to be constant
long double Set_C(long double, long double, long double, long double, long double);//Momentum dependence for the coupling constant, <number> 0 causes it to be constant

int main(int argc, char* argv[])
{
#ifdef BB	//use option -D BB= to activate bottomium macro
	char File[70] = "data/Spectralbb";  //Name of the file
#endif
#ifdef CC	//use option -D CC= to activate charmonium macro
	char File[70] = "data/Spectralcc";
#endif

#if VERSION == EXP	//use option -D VERSION={Exp,22,24,42} to select one of the potentials
	strcat(File,"Exp.");
#elif VERSION == 22
	strcat(File,"22.");
#elif VERSION == 24
	strcat(File,"24.");
#elif VERSION == 42
	strcat(File,"42.");
#endif

#if ORDER == 37	//use option -D ORDER={37,97} to select 37th order or 97th order Gauss-Legendre integration
	strcat(File, "37.");
#elif ORDER == 97
	strcat(File, "97.");
#endif

#ifdef HALF	//use option -D HALF= to divide self-energy in half
	strcat(File, "Half.");
#elif defined QUARTER	//use option -D HALF= to divide self-energy in half
	strcat(File, "Quarter.");
#endif
//	strcat(File, "Lambda.");
//	strcat(File, "500M");

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
		for(j = iProcess; j < 625; j+=Total)	//Does the subset of j that has been assigned to this process
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

long double Set_Lambda(long double G0, long double P, long double P0, long double fraction, int T)
{
	long double G = G0*(pow(P0,2)+(1-fraction)*pow(P,2))/(pow(P0,2)+pow(P,2));
	//long double G = G0;
	long double TempList[] = {0,.194,.258,.32,.4};
	long double Temp = TempList[T];

	/*long double Average_Lorentz_List[124][2] = {{0.,1},{0.1,1.0001803146917454},{0.2,1.0007209081856414},{0.30000000000000004,1.0016207327802138},{0.4,1.002878055663696},{0.5,1.0044904785443696},{0.6000000000000001,1.0064549643442586},{0.7000000000000001,1.0087678702182143},{0.8,1.0114249860172264},{0.9,1.014421577214468},{1.,1.0177524312569883},{1.1,1.0214119062939866},{1.2000000000000002,1.0253939812604214},{1.3,1.029692306356399},{1.4000000000000001,1.0343002530513208},{1.5,1.0392109628492994},{1.6,1.0444173941713806},{1.7000000000000002,1.049912366833576},{1.8,1.05568860372163},{1.9000000000000001,1.0617387693790339},{2.,1.0680555053304766},{2.1,1.0746314620563953},{2.2,1.0814593276142979},{2.3000000000000003,1.0885318529687331},{2.4000000000000004,1.0958418741445708},{2.5,1.1033823313585187},{2.6,1.111146285312763},{2.7,1.119126930853717},{2.8000000000000003,1.127317608209585},{2.9000000000000004,1.1357118120274672},{3.,1.1443031984068421},{3.1,1.153085590185483},{3.2,1.1620529806184545},{3.3000000000000003,1.1711995356830758},{3.4000000000000004,1.1805195951693013},{3.5,1.1900076727209794},{3.6,1.199658454976776},{3.7,1.2094667999457132},{3.8000000000000003,1.219427734738812},{3.9000000000000004,1.2295364527654935},{4.,1.2397883104912935},{4.8,1.3264497284933872},{5.6,1.4200654943600433},{6.4,1.519160875287701},{7.2,1.6226310974541502},{8.,1.7296410269866578},{8.8,1.8395516568015953},{9.600000000000001,1.951867608176621},{10.4,2.0661998229491414},{11.2,2.1822389236298974},{12.,2.2997360103264057},{12.8,2.4184886625357302},{13.600000000000001,2.5383306068682545},{14.400000000000002,2.6591239991060958},{15.2,2.7807535909240655},{16.,2.903122269882575},{16.8,3.026147613155822},{17.6,3.1497591974848635},{18.400000000000002,3.2738964794203413},{19.2,3.3985071100431115},{20.,3.523545583903495},{20.8,3.6489721474087435},{21.6,3.7747519103540004},{22.400000000000002,3.9008541178049714},{23.200000000000003,4.027251549525046},{24.000000000000004,4.153920021585342},{24.8,4.280837970397768},{25.6,4.407986103660698},{26.400000000000002,4.535347105956564},{27.200000000000003,4.6629053892445445},{28.000000000000004,4.790646881179819},{28.8,4.91855884085933},{29.6,5.046629706320803},{30.400000000000002,5.174848958075492},{31.200000000000003,5.303207003267246},{32.,5.431695074701866},{32.8,5.56030514288752},{33.6,5.6890298391629495},{34.4,5.8178623883093135},{35.2,5.946796549302059},{36.,6.07582656307305},{36.8,6.204947106329509},{37.6,6.334153250622065},{38.4,6.4634404259754135},{39.199999999999996,6.592804388496229},{40.,6.722241191457397},{40.8,6.851747159428774},{41.6,6.981318865084651},{42.4,7.110953108368514},{43.2,7.240646897738869},{44.,7.37039743325629},{44.8,7.500202091303125},{45.6,7.630058410753878},{46.4,7.759964080437278},{47.2,7.889916927750691},{48.,8.01991490830447},{48.8,8.14995609648863},{49.6,8.280038676866853},{50.4,8.410160936313938},{51.2,8.540321256822423},{52.,8.670518108912518},{52.8,8.800750045586758},{53.6,8.931015696777315},{54.4,9.061313764239484},{55.2,9.191643016849806},{56.,9.322002286271749},{56.8,9.452390462955641},{57.6,9.582806492442996},{58.4,9.713249371948384},{59.2,9.843718147194732},{60.,9.974211909480223},{60.8,10.104729792957222},{61.6,10.235270972105338},{62.4,10.365834659382706},{63.2,10.49642010304079},{64.,10.627026585089604},{64.8,10.757653419401212},{65.60000000000001,10.888299949940759},{66.4,11.01896554911495},{67.2,11.149649616229004},{68.,11.28035157604377},{68.8,11.411070877425496},{69.6,11.54180699208131},{70.4,11.67255941633873}};
	long double Average_Lorentz;	//Precalulated int_0^(pi/2)sqrt(1+(p*cos(theta)/m)^2)sin(theta), which should be the average Lorentz Factor over the sphere.
	int i = 0;
	while(abs(Average_Lorentz_List[i][0]-P) > 1e-5)
		i++;
	Average_Lorentz = Average_Lorentz_List[i][1];	//Store the average value for the given P*/

#if VERSION == 22
	return(sqrt(pow(1.23792404139016,2)+pow(G*Temp,2)));
#elif VERSION == 24
	return(sqrt(pow(2.024086607315311,2)+pow(G*Temp,2)/2));
#elif VERSION == 42
	return(pow(pow(2.433406283602878*1.25,4)+pow(G*Temp,4),.25));
	//return(pow(pow(5.0,4)+pow(G*Temp,4),.25));
	//return(pow(pow(.5,4)+pow(G*Temp,4),.25));
#elif VERSION == Exp
	return(sqrt(pow(2.6871146143427187`,2)+pow(G*Temp,2)));
#endif
}

long double Set_C(long double f0, long double P, long double P0, long double Lambda, long double fraction)
{
	long double f = (f0*pow(P0,2)+(fraction*(1-f0)+f0)*pow(P,2))/(pow(P0,2)+pow(P,2));

	/*long double Average_Lorentz_List[124][2] = {{0.,1},{0.1,1.0001803146917454},{0.2,1.0007209081856414},{0.30000000000000004,1.0016207327802138},{0.4,1.002878055663696},{0.5,1.0044904785443696},{0.6000000000000001,1.0064549643442586},{0.7000000000000001,1.0087678702182143},{0.8,1.0114249860172264},{0.9,1.014421577214468},{1.,1.0177524312569883},{1.1,1.0214119062939866},{1.2000000000000002,1.0253939812604214},{1.3,1.029692306356399},{1.4000000000000001,1.0343002530513208},{1.5,1.0392109628492994},{1.6,1.0444173941713806},{1.7000000000000002,1.049912366833576},{1.8,1.05568860372163},{1.9000000000000001,1.0617387693790339},{2.,1.0680555053304766},{2.1,1.0746314620563953},{2.2,1.0814593276142979},{2.3000000000000003,1.0885318529687331},{2.4000000000000004,1.0958418741445708},{2.5,1.1033823313585187},{2.6,1.111146285312763},{2.7,1.119126930853717},{2.8000000000000003,1.127317608209585},{2.9000000000000004,1.1357118120274672},{3.,1.1443031984068421},{3.1,1.153085590185483},{3.2,1.1620529806184545},{3.3000000000000003,1.1711995356830758},{3.4000000000000004,1.1805195951693013},{3.5,1.1900076727209794},{3.6,1.199658454976776},{3.7,1.2094667999457132},{3.8000000000000003,1.219427734738812},{3.9000000000000004,1.2295364527654935},{4.,1.2397883104912935},{4.8,1.3264497284933872},{5.6,1.4200654943600433},{6.4,1.519160875287701},{7.2,1.6226310974541502},{8.,1.7296410269866578},{8.8,1.8395516568015953},{9.600000000000001,1.951867608176621},{10.4,2.0661998229491414},{11.2,2.1822389236298974},{12.,2.2997360103264057},{12.8,2.4184886625357302},{13.600000000000001,2.5383306068682545},{14.400000000000002,2.6591239991060958},{15.2,2.7807535909240655},{16.,2.903122269882575},{16.8,3.026147613155822},{17.6,3.1497591974848635},{18.400000000000002,3.2738964794203413},{19.2,3.3985071100431115},{20.,3.523545583903495},{20.8,3.6489721474087435},{21.6,3.7747519103540004},{22.400000000000002,3.9008541178049714},{23.200000000000003,4.027251549525046},{24.000000000000004,4.153920021585342},{24.8,4.280837970397768},{25.6,4.407986103660698},{26.400000000000002,4.535347105956564},{27.200000000000003,4.6629053892445445},{28.000000000000004,4.790646881179819},{28.8,4.91855884085933},{29.6,5.046629706320803},{30.400000000000002,5.174848958075492},{31.200000000000003,5.303207003267246},{32.,5.431695074701866},{32.8,5.56030514288752},{33.6,5.6890298391629495},{34.4,5.8178623883093135},{35.2,5.946796549302059},{36.,6.07582656307305},{36.8,6.204947106329509},{37.6,6.334153250622065},{38.4,6.4634404259754135},{39.199999999999996,6.592804388496229},{40.,6.722241191457397},{40.8,6.851747159428774},{41.6,6.981318865084651},{42.4,7.110953108368514},{43.2,7.240646897738869},{44.,7.37039743325629},{44.8,7.500202091303125},{45.6,7.630058410753878},{46.4,7.759964080437278},{47.2,7.889916927750691},{48.,8.01991490830447},{48.8,8.14995609648863},{49.6,8.280038676866853},{50.4,8.410160936313938},{51.2,8.540321256822423},{52.,8.670518108912518},{52.8,8.800750045586758},{53.6,8.931015696777315},{54.4,9.061313764239484},{55.2,9.191643016849806},{56.,9.322002286271749},{56.8,9.452390462955641},{57.6,9.582806492442996},{58.4,9.713249371948384},{59.2,9.843718147194732},{60.,9.974211909480223},{60.8,10.104729792957222},{61.6,10.235270972105338},{62.4,10.365834659382706},{63.2,10.49642010304079},{64.,10.627026585089604},{64.8,10.757653419401212},{65.60000000000001,10.888299949940759},{66.4,11.01896554911495},{67.2,11.149649616229004},{68.,11.28035157604377},{68.8,11.411070877425496},{69.6,11.54180699208131},{70.4,11.67255941633873}};
	long double Average_Lorentz;	//Precalulated int_0^(pi/2)sqrt(1+(p*cos(theta)/m)^2)sin(theta), which should be the average Lorentz Factor over the sphere.
	int i = 0;
	while(abs(Average_Lorentz_List[i][0]-P) > 1e-5)
		i++;
	Average_Lorentz = Average_Lorentz_List[i][1];	//Store the average value for the given P*/

#if VERSION == 22
	return(181.0395540963551*f*pow(1.23792404139016/Lambda,4));
#elif VERSION == 24
	return(91.39940525625616*f*pow(2.024086607315311/Lambda,8));
#elif VERSION == 42
//	return(f*pow(.5/Lambda,8));
//	return(50.63740814101998*f*pow(2.433406283602878/Lambda,8));
	return(33.06143111723497*f*pow(2.433406283602878*1.25/Lambda,8));
//	return(14.7835935366527*f*pow(5.0/Lambda,8));
#elif VERSION == Exp
	return(66.3202722101555*f);
#endif
}
