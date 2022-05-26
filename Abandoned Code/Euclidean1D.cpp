#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cstring>
using namespace std;

long double Interpolation(long double, long double[][2]);
long double Dispersion(long double, long double[][2]);
long double EuclideanNon(long double, long double, long double[][2]);
long double EuclideanExt(long double, long double, long double);
long double EuclideanInter(long double, long double, long double, long double[][2], long double[][2], long double[][2], long double[][2]);
long double Euclidean(long double, long double, long double);
long double Spectral(long double, long double, long double[][2], long double[][2], long double[][2], long double[][2]);

long double C0; //vacuum coupling constant

int main(int argc, char* argv[])
{
	char FileName[100] = "data/Tables/Spectralcc";
	int length;
	long double ImGT[465][2], ImGvCT[465][2], ImGvLT[465][2], ImGvQT[465][2], ImGVT[465][2];
	long double ReGvCT[465][2], ReGvLT[465][2], ReGvQT[465][2], ReGVT[465][2];
	long double EuclideanNonList[4][11], EuclideanExtList[4][11], EuclideanInterCList[4][11][20], EuclideanInterLList[4][11][20], EuclideanInterQList[4][11][20];
	int Start, End;
	int i, j;
	long double Temp, tau, mod, MQ;
	long double Parameters[4];
	char Garbage[200];
	char* MQString;
	long double Line[9];

	//Open File
	strcat(FileName,argv[1]);
	strcat(FileName,".");
	strcat(FileName,argv[2]);
	strcat(FileName,".");
	strcat(FileName,argv[3]);
	strcat(FileName,".csv");
	ifstream iFile(FileName);
	if(iFile.good() == false)
	{
		cout << FileName << " is an invalid file." << endl;
		return(0);
	}

	//Get Parameters
	iFile.getline(Garbage,70);
	Parameters[0] = atof(strtok(Garbage,","));
	Parameters[1] = atof(strtok(NULL,","));
	MQ = Parameters[2] = atof(strtok(NULL,","));
	Parameters[3] = atof(strtok(NULL,","));

	//Ditch useless data
	for(i = 0; i < 151; i++)
		iFile.getline(Garbage,200);

	//Copy useful data into interpolation data tables (5 ways with 's' varible)
	for(i = 0; i < 465; i++)
	{
		iFile.getline(Garbage,200);
		Line[0] = atof(strtok(Garbage,","));
		for(j = 1; j < 9; j++)
		{
			Line[j] = atof(strtok(NULL,","));
		}
		ImGT[i][0] = ImGvCT[i][0] = ImGvLT[i][0] = ImGvQT[i][0] = ImGVT[i][0] = Line[3];
		ReGvCT[i][0] = ReGvLT[i][0] = ReGvQT[i][0] = ReGVT[i][0] = Line[3];
		ImGT[i][1] = Line[4];
		ImGvCT[i][1] = Line[5];
		ImGvLT[i][1] = Line[6];
		ImGvQT[i][1] = Line[7];
		ImGVT[i][1] = Line[8];
	}
	iFile.close();

	//Calculate 4 real part tables in interpolation form
	for(i = 0; i <= 416; i++)
	{
		ReGvCT[i][1] = Dispersion(ReGvCT[i][0], ImGvCT);
		ReGvLT[i][1] = Dispersion(ReGvLT[i][0], ImGvLT);
		ReGvQT[i][1] = Dispersion(ReGvQT[i][0], ImGvQT);
		ReGVT[i][1] = Dispersion(ReGVT[i][0], ImGVT);
	}

	//Set Temp
	switch(atoi(argv[2]))
	{
		case 0:
		case 1:
			Start = 1;
			End = 4;
			break;
		case 2:
		case 3:
		case 4:
			Start = End = atoi(argv[2]);
			break;
	}
	//Set C0
	if(strncmp(argv[1],"22",2)==0)
		C0 = 116.253;
	else if(strncmp(argv[1],"24",2)==0)
		C0 = 65.8549;
	else if(strncmp(argv[1],"42",2)==0)
		C0 = 38.4541;
	else if(strncmp(argv[1],"Ex",2)==0)
		C0 = 50.3627;

	//Euclidean correlation functions
	for(int T = Start; T <= End; T++)
	{
		switch(T)
		{
			case 1:
				Temp = .194;
				break;
			case 2:
				Temp = .258;
				break;
			case 3:
				Temp = .32;
				break;
			case 4:
				Temp = .4;
				break;
		}

		for(i = 0; i <= 10; i++)
		{
			tau = i/(20.*Temp);
			EuclideanNonList[T-1][i] = EuclideanNon(tau, Temp, ImGT);
			EuclideanExtList[T-1][i] = EuclideanExt(tau, Temp, MQ);
			for(j = 0; j < 20; j++)
			{
				mod = .8+j*.01;
				EuclideanInterCList[T-1][i][j] = EuclideanInter(tau, Temp, mod, ReGvCT, ImGvCT, ReGVT, ImGVT)/4.;
				EuclideanInterLList[T-1][i][j] = EuclideanInter(tau, Temp, mod, ReGvLT, ImGvLT, ReGVT, ImGVT);
				EuclideanInterQList[T-1][i][j] = EuclideanInter(tau, Temp, mod, ReGvQT, ImGvQT, ReGVT, ImGVT);
			}
		}
	}

	//Open File for data export to Mathematica
	length = strlen(FileName);	//Find the end of the string (the location of the null character)
	FileName[length-3] = char(0);	//Backup 3 spaces and insert a new null character removing the "csv" from the string
	strcat(FileName,"m");		//Replace the \0sv with m\0 making it a Mathematica file output
	ofstream oFile(FileName);

	//Export interplation table, {s,ImG,ImGvC,ImGvL,ImGvQ,ImGV,ReGvC,ReGvL,ReGvQ,ReGV}
	oFile << "{{" << Parameters[0] << "," << Parameters[1] << "," << Parameters[2] << "," << Parameters[3] << "}," << endl << "{";
	oFile << setprecision(18);
	for(i = 0; i < 416; i++)
		oFile << "{" << ImGT[i][0] << "," << ImGT[i][1] << "," << ImGvCT[i][1] << "," << ImGvLT[i][1] << "," << ImGvQT[i][1] << "," << ImGVT[i][1] << "," << ReGvCT[i][1] << "," << ReGvLT[i][1] << "," << ReGvQT[i][1] << "," << ReGVT[i][1] << "}," << endl;
	for(; i < 465; i++)
	{
		oFile << "{" << ImGT[i][0] << "," << ImGT[i][1] << "," << ImGvCT[i][1] << "," << ImGvLT[i][1] << "," << ImGvQT[i][1] << "," << ImGVT[i][1] << ",Null,Null,Null,Null}" << flush;
		if(i != 464)
			oFile << "," << endl;
	}
	oFile << "},\n{" << flush;

	//Export EuclideanNon correlation table
	for(int T = Start; T <= End; T++)
	{
		if(Start != End)
			oFile << "{" << flush;
		for(i = 0; i < 11; i++)
		{
			oFile << EuclideanNonList[T-1][i];
			if(i != 10)
				oFile << ",";
		}
		if(Start != End)
			oFile << "}" << flush;
		if(T != End)
			oFile << "," << endl;
	}
	oFile << "},\n{" << flush;

	//Export EuclideanInterC correlation table
	for(int T = Start; T <= End; T++)
	{
		if(Start != End)
			oFile << "{" << flush;
		for(i = 0; i < 11; i++)
		{
			oFile << "{" << flush;
			for(j = 0; j < 20; j++)
			{
				oFile << EuclideanInterCList[T-1][i][j];
				if(j != 19)
					oFile << ",";
			}
			oFile << "}";
			if(i != 10)
				oFile << "," << endl;
		}
		if(Start != End)
			oFile << "}" << flush;
		if(T != End)
			oFile << "," << endl;
	}
	oFile << "},\n{" << flush;

	//Export EuclideanInterL correlation table
	for(int T = Start; T <= End; T++)
	{
		if(Start != End)
			oFile << "{" << flush;
		for(i = 0; i < 11; i++)
		{
			oFile << "{" << flush;
			for(j = 0; j < 20; j++)
			{
				oFile << EuclideanInterLList[T-1][i][j];
				if(j != 19)
					oFile << ",";
			}
			oFile << "}";
			if(i != 10)
				oFile << "," << endl;
		}
		if(Start != End)
			oFile << "}" << flush;
		if(T != End)
			oFile << "," << endl;
	}
	oFile << "},\n{" << flush;

	//Export EuclideanInterQ correlation table
	for(int T = Start; T <= End; T++)
	{
		if(Start != End)
			oFile << "{" << flush;
		for(i = 0; i < 11; i++)
		{
			oFile << "{" << flush;
			for(j = 0; j < 20; j++)
			{
				oFile << EuclideanInterQList[T-1][i][j];
				if(j != 19)
					oFile << ",";
			}
			oFile << "}";
			if(i != 10)
				oFile << "," << endl;
		}
		if(Start != End)
			oFile << "}" << flush;
		if(T != End)
			oFile << "," << endl;
	}
	oFile << "},\n{" << flush;

	//Export EuclideanExt correlation table
	for(int T = Start; T <= End; T++)
	{
		if(Start != End)
			oFile << "{" << flush;
		for(i = 0; i < 11; i++)
		{
			oFile << EuclideanExtList[T-1][i];
			if(i != 10)
				oFile << ",";
		}
		if(Start != End)
			oFile << "}" << flush;
		if(T != End)
			oFile << "," << endl;
	}
	oFile << "}}" << endl;

	cout << "sed -i 's/e/*^/g' " << FileName << endl;
	oFile.close();
	return(0);
}

long double EuclideanExt(long double tau, long double Temp, long double MQ) //given quark mass, tau, and temp, return Euclidean correlation function for non-interacting spectral function covering sqrt(s) > 23.5 GeV
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double Intervals[] = {23.5,100.,200.,300.};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = Intervals[0];
	long double b;
	int i, j;
	long double holder;

	for(i = 1; i < 4; i++)
	{
		b = pow(Intervals[i],2);
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += Euclidean(x1, tau, Temp)*(3.*x1)/(pow(4.*M_PI,2)*sqrt(x1))*sqrt(1.-pow(2.*MQ,2)/x1)*w[j+1];	//Evaluate k integral at x1
			F_b += Euclidean(x2, tau, Temp)*(3.*x2)/(pow(4.*M_PI,2)*sqrt(x2))*sqrt(1.-pow(2.*MQ,2)/x2)*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Euclidean(a/2.+b/2., tau, Temp)*(3.*(a/2.+b/2.))/(pow(4.*M_PI,2)*sqrt(a/2.+b/2.))*sqrt(1.-pow(2.*MQ,2)/(a/2.+b/2.))*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	return(Answer);
}

long double EuclideanInter(long double tau, long double Temp, long double mod, long double ReGv[][2], long double ImGv[][2], long double ReGV[][2], long double ImGV[][2]) //given tables (ReGv, ImGv, ReGV, ImGV), tau, coupling constant modifier, and temp, return Euclidean correlation function for interacting spectral function
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double Intervals[] = {0.,2.99,3.02,3.025,3.03,3.0375,3.04,3.0425,3.05,3.055,3.06,3.09,3.1,3.2,3.3,3.4,3.5,3.6,4.,10.,23.5};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = Intervals[0];
	long double b;
	int i, j;
	long double holder;

	for(i = 1; i < 21; i++)
	{
		b = pow(Intervals[i],2);
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += Euclidean(x1, tau, Temp)*Spectral(x1, mod, ReGv, ImGv, ReGV, ImGV)/(2.*sqrt(x1))*w[j+1];	//Evaluate k integral at x1
			F_b += Euclidean(x2, tau, Temp)*Spectral(x2, mod, ReGv, ImGv, ReGV, ImGV)/(2.*sqrt(x2))*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = Euclidean(a/2.+b/2., tau, Temp)*Spectral(a/2.+b/2., mod, ReGv, ImGv, ReGV, ImGV)/(2.*sqrt(a/2.+b/2.))*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}
	return(Answer);
}

long double EuclideanNon(long double tau, long double Temp, long double ImG[][2]) //given table (ImG), tau, and temp, return Euclidean correlation function for non-interacting spectral function
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double Intervals[] = {0.,3.,3.1,3.2,3.3,3.4,3.5,3.6,4.,10.,23.5};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = Intervals[0];
	long double b;
	int i, j;
	long double holder;

	for(i = 1; i < 11; i++)
	{
		b = pow(Intervals[i],2);
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += -3./M_PI*Euclidean(x1, tau, Temp)*Interpolation(x1, ImG)/(2.*sqrt(x1))*w[j+1];	//Evaluate k integral at x1
			F_b += -3./M_PI*Euclidean(x2, tau, Temp)*Interpolation(x2, ImG)/(2.*sqrt(x2))*w[j+1];	//Evaluate k integral at x3
		}
		F_ave = -3./M_PI*Euclidean(a/2.+b/2., tau, Temp)*Interpolation(a/2.+b/2., ImG)/(2.*sqrt(a/2.+b/2.))*w[0];
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	return(Answer);
}

long double Euclidean(long double s, long double tau, long double Temp)
{
	return(cosh(sqrt(s)*(tau-1./(2.*Temp)))/sinh(sqrt(s)/(2.*Temp)));
}

long double Dispersion(long double s, long double table[][2]) //given table (imagainary function), calculate the dispersion relation for the given s
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204}; //Weight of the function at Disp
	long double DispLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre quadrature
	long double wLa[] = {0.07496328305102102808055, 0.1745735743605928864303, 0.2745074833881225250022, 0.3747323102655645620060, 0.4753412526072084401161, 0.5764380939967183636147, 0.6781307242364945406823, 0.7805307978511547593175, 0.8837542316062452388883, 0.9879219194279636096671, 1.0931605619330277996916, 1.1996035979670979427973, 1.3073922479469277349326, 1.416676687469297701993, 1.5276173754408796787012, 1.640386566702889623924, 1.7551700457872174635214, 1.8721691266543402861779, 1.9916029736088098866132, 2.1137113117669909276048, 2.2387576123844772725684, 2.3670328602831611098048, 2.4988600392644108123394, 2.6345995091430390709, 2.7746554982525006307172, 2.9194840027576204632431, 3.0696024758091833914472, 3.2256018156600758204608, 3.3881613374746331979827, 3.5580676615951707296054, 3.7362388067183244743069, 3.9237552950635210172968, 4.1219008467729629867363, 4.3322164077399479741288, 4.5565730632309056055423, 4.7972722621195591678357, 5.057186469320242487569, 5.3399612774797865633198, 5.6503138450512931300331, 5.9944877492232503537552, 6.3809726096501927329094, 6.8216946862388774056326, 7.3340972531892936469048, 7.9450326451948326187906, 8.6987143462393085933469, 9.6750102652900375180015, 11.039313738067347840094, 13.220456867750092021034, 17.982575250664959108273};	//Weights for 95th order Gauss-Laguerre quadrature
	long double Range[] = {0, 3.24, 12.96, 25, 100, 225, 552.25};
	long double Answer;
	long double F_a, F_b, F_ave;
	long double x1, x2;
	long double a = 0, b;
	int i, j;

	Answer = Interpolation(s, table)*log(abs((552.25-s)/(s)));

	if(abs(s) < pow(.001,2)*.1)	//Take linear limits to subvert issues with the endpoints of log((552.25-s)/(s-a)) causing issues
		return(2.*Dispersion(s+.001, table)-Dispersion(s+.002, table));
	else if(abs(s-552.25) < .00001)
		return(-Dispersion(pow(23.498,2), table)+2.*Dispersion(pow(23.499,2), table));

	b = 0;
	for(i = 0; i < 7; i++)
	{
		b = Range[i];

		F_a = 0;
		F_b = 0;
		for(j = 0; j < 9; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.; //Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			F_a += (Interpolation(x1, table)-Interpolation(s, table))/(x1-s)*w[j+1]; //Evaluate function at x1
			F_b += (Interpolation(x2, table)-Interpolation(s, table))/(x2-s)*w[j+1]; //Evaluate function at x2
		}
		F_ave = (Interpolation((a+b)/2., table)-Interpolation(s, table))/((a+b)/2.-s)*w[0]; //Evaluate function at (a+b)/2.
		Answer += (F_a+F_ave+F_b)*(b-a)/2.;
		a = b;
	}

	F_a = 0;
	for(j = 0; j < 49; j++)
	{
		x1 = 552.25+DispLa[j];
		F_a += table[415+j][1]/(x1-s)*wLa[j]; //Evaluate function at x1
	}
	Answer += F_a;

	return(Answer/M_PI);
}

long double Spectral(long double s, long double mod, long double ReGvT[][2], long double ImGvT[][2], long double ReGVT[][2], long double ImGVT[][2])//given tables (ReGv, ImGv, ReGV, ImGV), coupling constant modifier, and s, return interacting spectral function
{
	long double ReGv = Interpolation(s, ReGvT);
	long double ImGv = Interpolation(s, ImGvT);
	long double ReGV = Interpolation(s, ReGVT);
	long double ImGV = Interpolation(s, ImGVT);

	return((3./M_PI*C0*mod*(ImGV*mod*(pow(ReGv,2)-pow(ImGv,2))-2.*ImGv*ReGv*(ReGV*mod-1.)))/(pow(ImGV*mod,2)+pow(ReGV*mod-1.,2)));
}

long double Interpolation(long double s, long double table[][2])
{
	int i;
	if(s <= 9)
		i = 10.*sqrt(s);
	else if(s <= 25)
		i = 100.*sqrt(s)-270.;
	else if(s <= 552.25)
		i = 180.+10.*sqrt(s);
	else
	{
		i = 416;
		while(table[i][0] <= s)
			i++;
	}

	long double xi[4];
	long double yi[4];
	long double denominator;
	long double coefficents[4];

	if(i >= 1 && i <= 464)
		i--;
	else if(i > 464)
		i -= 2;

	xi[0] = table[i][0]-table[i+1][0];
	xi[1] = 0;
	xi[2] = table[i+2][0]-table[i+1][0];
	xi[3] = table[i+3][0]-table[i+1][0];
	yi[0] = table[i][1];
	yi[1] = table[i+1][1];
	yi[2] = table[i+2][1];
	yi[3] = table[i+3][1];
	s -= table[i+1][0];

	denominator = xi[0]*xi[2]*xi[3]*(xi[0]-xi[2])*(xi[0]-xi[3])*(xi[2]-xi[3]);
	coefficents[0] = xi[0]*(xi[0]-xi[3])*xi[3]*(yi[1]-yi[2])+pow(xi[2],2)*(xi[3]*(yi[0]-yi[1])+xi[0]*(yi[1]-yi[3]))+xi[2]*(pow(xi[3],2)*(-yi[0]+yi[1])+pow(xi[0],2)*(-yi[1]+yi[3]));
	coefficents[1] = -(xi[0]*xi[3]*(pow(xi[0],2)-pow(xi[3],2))*(yi[1]-yi[2]))+xi[2]*(pow(xi[3],3)*(yi[0]-yi[1])+pow(xi[0],3)*(yi[1]-yi[3]))+pow(xi[2],3)*(xi[3]*(-yi[0]+yi[1])+xi[0]*(-yi[1]+yi[3]));
	coefficents[2] = pow(xi[0],2)*(xi[0]-xi[3])*pow(xi[3],2)*(yi[1]-yi[2])+pow(xi[2],3)*(pow(xi[3],2)*(yi[0]-yi[1])+pow(xi[0],2)*(yi[1]-yi[3]))+pow(xi[2],2)*(pow(xi[3],3)*(-yi[0]+yi[1])+pow(xi[0],3)*(-yi[1]+yi[3]));
	coefficents[3] = xi[0]*(xi[0]-xi[2])*xi[2]*(xi[0]-xi[3])*(xi[2]-xi[3])*xi[3]*yi[1];
	coefficents[0] /= denominator;
	coefficents[1] /= denominator;
	coefficents[2] /= denominator;
	coefficents[3] /= denominator;

	return(coefficents[0]*pow(s,3)+coefficents[1]*pow(s,2)+coefficents[2]*s+coefficents[3]);
}
