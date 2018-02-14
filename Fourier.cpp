#include<fstream>
#include<iostream>
#include<cmath>
#include<complex>
#include<cstring>
#include<cstdlib>
#include<omp.h>
#include<iomanip>
using namespace std;

bool ReadIn(long double***[], int[], int[], char*);	//Read in from file the contents of the matrix. It includes 2 offsets for multiple file read in and for padding due to negative momentum not being in the file
void Init(long double***[], int[], int[]);	//Initialize the table to the correct size
void Validate(long double***[], int[], int[]);	//Checks the points determine if they are valid or if direct calls to the spectral function are needed
long double Spectral(long double***[], long double, long double, long double, int);	//The tabulated points and 2 inputs and returns the bicubic interpolation of that input
long double Bicubic(long double***[], long double, long double, long double, int);	//Bicubic interpolation
long double Correlator(long double(*)(long double***[], long double, long double, int), long double***[], long double, int);	//Evaluates the spatial correlator
long double SpatialNeg(long double***[], long double, long double, int);	//1D integral for spatial integrator with fixed j, negative s
long double SpatialDebug(long double***[], long double, long double, int);	//1D integral for spatial integrator with fixed j, negative s
long double SpatialPos(long double***[], long double, long double, int);	//1D integral for spatial integrator with fixed E, positive s
long double SpatialVac(long double***[], long double, long double, int);	//1D analytic integral for vacuum on assumption of Lorentz invariance
long double Euclidean(long double***[], long double, long double, int);		//Kernal for eucledian-time correlator
long double Spectral_Analytic(long double, long double);	//Analytic spectral function fit
long double ImG(long double[4], long double);	//The 5 inputs to get the spectral function back
long double ReGV1(long double[4], long double);
long double ImGV1(long double[4], long double);
long double ReGV2(long double[4], long double);
long double ImGV2(long double[4], long double);
complex<long double> atanh(complex<long double>);
long double Tail_Factor(long double);

char* Process;
long double List[16];

int main(int argc, char* argv[])	//Process, # of Process, Output file, Input file, Temp
{
	char* File = new char[150];	//Name of the file
	strcpy(File, argv[3]);
	Process = argv[1];
	strcat(File, ".");
	strcat(File, argv[5]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot(File);
	TPlot << setprecision(18);
	long double*** Table[15];	//The table of values computed by Spectral and Dispersion
	long double z, tau;	//The position value of the spactial correlator and tau of the euclidean-time correlator
	long double holder[5];
	int N[15] = {151,151,465,151,151,465,151,151,465,151,151,465,151,151,465};	//The size of the table
	int M[15] = {209,581,752,209,581,752,209,581,752,209,581,752,209,581,752};	//The size of the table
	const int iProcess = atoi(argv[1]);
	const int Total = atoi(argv[2]);
	const int Temp = atoi(argv[5]);

	for(int i = 0; i < 9; i++)
		List[i] = atof(argv[i+6]);

	/*Init(Table, N, M);
	if(!ReadIn(Table, N, M, argv[4]))
		return(0);*/
	//Validate(Table, N, M);

//Debugging code, used for checking the contents of the Table to ensure that it was filled correctly
	/*const long double EList[] = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.540308, 2.545308, 2.550308, 2.555308, 2.560308, 2.565308, 2.570308, 2.575308, 2.580308, 2.585308, 2.590308, 2.595308, 2.600308, 2.605308, 2.610308, 2.615308, 2.620308, 2.625308, 2.630308, 2.635308, 2.640308, 2.645308, 2.650308, 2.655308, 2.660308, 2.665308, 2.670308, 2.675308, 2.680308, 2.685308, 2.690308, 2.695308, 2.700308, 2.705308, 2.710308, 2.715308, 2.720308, 2.725308, 2.730308, 2.735308, 2.740308, 2.745308, 2.750308, 2.755308, 2.760308, 2.765308, 2.770308, 2.775308, 2.780308, 2.785308, 2.790308, 2.795308, 2.800308, 2.805308, 2.810308, 2.815308, 2.820308, 2.825308, 2.830308, 2.835308, 2.840308, 2.845308, 2.850308, 2.855308, 2.860308, 2.865308, 2.870308, 2.875308, 2.880308, 2.885308, 2.890308, 2.895308, 2.900308, 2.905308, 2.910308, 2.915308, 2.920308, 2.925308, 2.930308, 2.935308, 2.940308, 2.945308, 2.950308, 2.955308, 2.960308, 2.965308, 2.970308, 2.975308, 2.980308, 2.985308, 2.990308, 2.995308, 3.000308, 3.005308, 3.010308, 3.015308, 3.020308, 3.025308, 3.030308, 3.035308, 3.040308, 3.045308, 3.050308, 3.055308, 3.060308, 3.065308, 3.070308, 3.075308, 3.080308, 3.085308, 3.090308, 3.095308, 3.100308, 3.105308, 3.110308, 3.115308, 3.120308, 3.125308, 3.130308, 3.135308, 3.140308, 3.145308, 3.150308, 3.155308, 3.160308, 3.165308, 3.170308, 3.175308, 3.180308, 3.185308, 3.190308, 3.195308, 3.200308, 3.205308, 3.210308, 3.215308, 3.220308, 3.225308, 3.230308, 3.235308, 3.240308, 3.245308, 3.250308, 3.255308, 3.260308, 3.265308, 3.270308, 3.275308, 3.280308, 3.285308, 3.290308, 3.295308, 3.300308, 3.305308, 3.310308, 3.315308, 3.320308, 3.325308, 3.330308, 3.335308, 3.340308, 3.345308, 3.350308, 3.355308, 3.360308, 3.365308, 3.370308, 3.375308, 3.380308, 3.385308, 3.390308, 3.395308, 3.400308, 3.405308, 3.410308, 3.415308, 3.420308, 3.425308, 3.430308, 3.435308, 3.440308, 3.445308, 3.450308, 3.455308, 3.460308, 3.465308, 3.470308, 3.475308, 3.480308, 3.485308, 3.490308, 3.495308, 3.500308, 3.505308, 3.510308, 3.515308, 3.520308, 3.525308, 3.530308, 3.535308, 3.540308, 3.55, 3.56375, 3.5775, 3.59125, 3.605, 3.61875, 3.6325, 3.64625, 3.66, 3.67375, 3.6875,  3.70125, 3.715, 3.72875, 3.7425, 3.75625, 3.77, 3.78375, 3.7975, 3.81125, 3.825, 3.83875, 3.8525, 3.86625, 3.88, 3.89375, 3.9075, 3.92125, 3.935, 3.94875, 3.9625, 3.97625, 3.99, 4.00375, 4.0175, 4.03125, 4.045, 4.05875, 4.0725, 4.08625, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21., 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22., 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23., 23.1, 23.2, 23.3, 23.4, 23.5};
	for(int i = 1; i < 376; i++)
	{
		for(int j = 0; j < 462; j++)
			//cout << (i-1)*.8 << " " << EList[j] << " " << Table[i][j][0] << " " << Table[i][j][1] << " " << Table[i][j][2] << " " << Table[i][j][3] << " " << Table[i][j][4] << endl;
		//cout << endl;
	}//*/

//Debugging code, used to ensure that the interpolations or othe approximations are working correctly
	/*//cout << setprecision(18);
	for(long double j = 0; j <= 0; j+=.4)
	{
		for(long double i = 0; i <= 23.5; i+=.0002)
			//cout << j << " " << i << " " << Spectral(Table, i, j, -1, 1) << endl;
		//cout << endl;
	}//*/
	/*for(long double j = 0; j <= 25; j+=.4)
	{
		for(long double i = 0; i < 400; i+=.25)
		{
			if(j < 13)
				//cout << j << " " << i << " " << Spectral(Table, i, j, -1, 0) << endl;
			if(j >= 13)
				//cout << j << " " << i << " " << Spectral(Table, i, j, -1, 1) << endl;
		}
		//cout << endl;
	}//*/

//Debugging code, examines the difference between interpolations and a finite-width approximation
	/*for(long double i = 3.02; i <= 3.06; i+=.00005)
	{
		for(long double j = 70; j <= 75; j+=.02)
			//cout << i << " " << j << " " << Spectral(Table,i,j) << " " << Analytic(i,Epsilon) << endl;
		//cout << endl;
	}//*/

//The actual program
	/*TPlot << "#" << flush;
	for(int i = 0; i < argc; i++)
		TPlot << argv[i] << " " << flush;
	TPlot << endl;

	if(atoi(argv[1]) == 0)
	{
		TPlot << endl;
		for(long double P = 0; P <= 100.1; P += .8)
		{
			for(long double roots = 0; roots <= 23.5; roots += .01)
				TPlot << roots << " " << P << " " << Spectral_Analytic(roots*roots, P) << endl;
			TPlot << endl;
		}
	}*/

	#pragma omp parallel for private(tau, z, holder)
	for(int i = 290*iProcess/Total; i <= 290*(iProcess+1)/Total; i+=1)
	{
		z = .3+i*.02;
		tau = i*.008;
		holder[0] = Correlator(SpatialDebug, Table, z, Temp);
		holder[1] = 0;//Correlator(SpatialNeg, Table, z, Temp);
		holder[2] = Correlator(SpatialPos, Table, z, Temp);
		holder[3] = Correlator(SpatialVac, Table, z, Temp);
		holder[4] = Correlator(Euclidean, Table, tau, Temp);
		#pragma omp critical
		{
			TPlot << z << " " << holder[0] << " " << holder[1] << " " << holder[2] << " " << holder[3] << " " << tau << " " << holder[4] << endl;
		}
	}//*/

	return(0);
}

long double Correlator(long double(*Kernal)(long double***[], long double, long double, int), long double*** Table[], long double z, int Temp)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double distance[] = {5e-2, 2e-2, 1.5e-2, 1e-2, 2.5e-3, 1e-4, 1e-5, 1e-6};	//Stride of the integral
	long double Answer = 0;
	long double F_a, F_b, F_ave;
	long double a = 0;
	long double b = 0;
	int i, j;
	long double holder;
	//ofstream oTable("sTable", ios::app);

	if(Kernal != SpatialNeg)
	{
		for(i = 0; i < 8; i++)
		{
			b = 3.040308 - distance[i];
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		b = 3.040308;
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			holder = Kernal(Table, x1, z, Temp);
			F_a += holder*w[j+1];	//Evaluate k integral at x1
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x1 << " " << holder << endl;
			holder = Kernal(Table, x2, z, Temp);
			F_b += holder*w[j+1];	//Evaluate k integral at x3
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x2 << " " << holder << endl;
		}
		F_ave = Kernal(Table, a/2.+b/2., z, Temp);
		if(Kernal == SpatialDebug)
			cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;

		for(i = 7; i >= 0; i--)
		{
			b = 3.040308 + distance[i];
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		for(b = 3.1; b <= 3.4; b += .05)
		{
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		for(b = 3.4; b <= 3.45; b += .01)
		{
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		for(i = 7; i >= 0; i--)
		{
			b = 3.6 + distance[i];
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
				if(Kernal == SpatialDebug)
					cout << "1 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}

		b = 4;
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			holder = Kernal(Table, x1, z, Temp);
			F_a += holder*w[j+1];	//Evaluate k integral at x1
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x1 << " " << holder << endl;
			holder = Kernal(Table, x2, z, Temp);
			F_b += holder*w[j+1];	//Evaluate k integral at x3
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x2 << " " << holder << endl;
		}
		F_ave = Kernal(Table, a/2.+b/2., z, Temp);
		if(Kernal == SpatialDebug)
			cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;

		b = 10;
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			holder = Kernal(Table, x1, z, Temp);
			F_a += holder*w[j+1];	//Evaluate k integral at x1
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x1 << " " << holder << endl;
			holder = Kernal(Table, x2, z, Temp);
			F_b += holder*w[j+1];	//Evaluate k integral at x3
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x2 << " " << holder << endl;
		}
		F_ave = Kernal(Table, a/2.+b/2., z, Temp);
		if(Kernal == SpatialDebug)
			cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;

		b = 22.5;
		F_a = F_b = 0;	//Start integration at 0
		for(j = 0; j < 24; j++)
		{
			x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[j]*(b-a))/2.;

			holder = Kernal(Table, x1, z, Temp);
			F_a += holder*w[j+1];	//Evaluate k integral at x1
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x1 << " " << holder << endl;
			holder = Kernal(Table, x2, z, Temp);
			F_b += holder*w[j+1];	//Evaluate k integral at x3
			if(Kernal == SpatialDebug)
				cout << "1 " << z << " " << x2 << " " << holder << endl;
		}
		F_ave = Kernal(Table, a/2.+b/2., z, Temp);
		if(Kernal == SpatialDebug)
			cout << "1 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
	}
	else
	{
		long double stride = 1.;
		do
		{
			b += stride;
			F_a = F_b = 0;	//Start integration at 0
			for(j = 0; j < 24; j++)
			{
				x1 = (b+a-Disp[j]*(b-a))/2.;	//Actual evaluation points
				x2 = (b+a+Disp[j]*(b-a))/2.;

				holder = Kernal(Table, x1, z, Temp);
				F_a += holder*w[j+1];	//Evaluate k integral at x1
					cout << "2 " << z << " " << x1 << " " << holder << endl;
				holder = Kernal(Table, x2, z, Temp);
				F_b += holder*w[j+1];	//Evaluate k integral at x3
					cout << "2 " << z << " " << x2 << " " << holder << endl;
			}
			F_ave = Kernal(Table, a/2.+b/2., z, Temp);
				cout << "2 " << z << " " << a/2.+b/2. << " " << F_ave << endl;
			Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
			a = b;
		}while(b < 15.);
	}

	return(Answer);
}

long double SpatialNeg(long double*** Table[], long double x, long double z, int Temp)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 95th order Gaussian quadrature for this interval
	long double x2;
	long double Answer = 0;
	long double stride = 2.*M_PI/z;	//Stride of the integral
	long double F_a, F_b, F_ave;
	long double a = x;
	long double b = x;	//The location of where the current stride should have ended
	int i;

	do
	{
		b += stride;

		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[i]*(b-a))/2.;

			F_a += 2.*Spectral(Table, pow(x,2)-2.*x*x1, x1, z, 1)*cos(x1*z)/(x1-x)*w[i+1];	//Evaluate k integral at x1
			F_b += 2.*Spectral(Table, pow(x,2)-2.*x*x2, x2, z, 1)*cos(x2*z)/(x2-x)*w[i+1];	//Evaluate k integral at x3
		}
		F_ave = 2.*Spectral(Table, pow(x,2)-2.*x*(a+b)/2., (a+b)/2., z, 1)*cos((a+b)*z/2.)/(a+b-x);
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/2.;//*/
		a = b;
	}while(b < 2001.*M_PI/z);

	b = 2001.*M_PI/z;

	F_a = F_b = 0;	//Start integration at 0
	for(i = 0; i < 24; i++)
	{
		x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
		x2 = (b+a+Disp[i]*(b-a))/2.;

		F_a += 2.*Spectral(Table, pow(x,2)-2.*x*x1, x1, z, 1)*cos(x1*z)/(x1-x)*w[i+1];	//Evaluate k integral at x1
		F_b += 2.*Spectral(Table, pow(x,2)-2.*x*x2, x2, z, 1)*cos(x2*z)/(x2-x)*w[i+1];	//Evaluate k integral at x3
	}
	F_ave = 2.*Spectral(Table, pow(x,2)-2.*x*(a+b)/2., (a+b)/2., z, 1)*cos((a+b)*z/2.)/(a+b-x);
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/2.;//*/

	return(Answer);	//return the best estimate of the integral on the interval*/
}

long double SpatialDebug(long double*** Table[], long double roots, long double z, int Temp)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double Answer = 0;
	long double stride = 2*M_PI/z;	//Stride of the integral
	long double F_a, F_b, F_ave;
	long double a = 0;
	long double b = 0;
	int i;

//This code is for integating out to a z dependant boundary
	for(int j = 0; j < 1000; j++)	//need to start the count off from where it left off in the previous integration block
	{
		b += stride;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[i]*(b-a))/2.;

			F_a += (Spectral(Table, roots, x1, z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(x1,2)+pow(roots,2))*cos(z*x1)*w[i+1];	//Evaluate k integral at x1
			F_b += (Spectral(Table, roots, x2, z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(x2,2)+pow(roots,2))*cos(z*x2)*w[i+1];	//Evaluate k integral at x3
		}
		F_ave = (Spectral(Table, roots, a/2.+b/2., z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(a/2.+b/2.,2)+pow(roots,2))*cos(z*(a/2.+b/2.));
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}	//For the bulk of the integral where either the result is well approximated by either the finite or zero width analytic result

	b += stride/2.;	//Evaluate an extra half stride, should result in a result between the extrems of a quarter and three quarter strides
	F_a = F_b = 0;	//Start integration at 0
	for(i = 0; i < 24; i++)
	{
		x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
		x2 = (b+a+Disp[i]*(b-a))/2.;

		F_a += (Spectral(Table, roots, x1, z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(x1,2)+pow(roots,2))*cos(z*x1)*w[i+1];	//Evaluate k integral at x1
		F_b += (Spectral(Table, roots, x2, z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(x2,2)+pow(roots,2))*cos(z*x2)*w[i+1];	//Evaluate k integral at x3
	}
	F_ave = (Spectral(Table, roots, a/2.+b/2., z, 2)-Spectral(Table, roots, 0., -1, 2))*4.*roots/(pow(a/2.+b/2.,2)+pow(roots,2))*cos(z*(a/2.+b/2.));
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);//*/

	return(Answer);	//return the best estimate of the integral on the interval*/
}
long double SpatialPos(long double*** Table[], long double roots, long double z, int Temp)
{
	long double Disp[] = {0.06342068498268678602883,  0.1265859972696720510680, 0.1892415924618135864853,  0.2511351786125772735072, 0.3120175321197487622079,  0.3716435012622848888637, 0.4297729933415765246586,  0.4861719414524920421770, 0.5406132469917260665582,  0.5928776941089007124559, 0.6427548324192376640569,  0.6900438244251321135048, 0.7345542542374026962137,  0.7761068943454466350181, 0.8145344273598554315395,  0.8496821198441657010349, 0.8814084455730089100370,  0.9095856558280732852130, 0.9341002947558101490590,  0.9548536586741372335552, 0.9717622009015553801400,  0.9847578959142130043593, 0.9937886619441677907601,  0.9988201506066353793618};	//Dispacement from center
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of data point
	long double x1;	//These are the two other points required for 97th order Gaussian quadrature for this interval
	long double x2;
	long double Answer = 0;
	long double stride = 2*M_PI/z;	//Stride of the integral
	long double F_a, F_b, F_ave;
	long double a = 0;
	long double b = 0;
	int i;

//This code is for integating out to a z dependant boundary
	for(int j = 0; j < 1000; j++)	//need to start the count off from where it left off in the previous integration block
	{
		b += stride;
		F_a = F_b = 0;	//Start integration at 0
		for(i = 0; i < 24; i++)
		{
			x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
			x2 = (b+a+Disp[i]*(b-a))/2.;

			F_a += Spectral(Table, roots, x1, z, 2)*4.*roots/(pow(x1,2)+pow(roots,2))*cos(z*x1)*w[i+1];	//Evaluate k integral at x1
			F_b += Spectral(Table, roots, x2, z, 2)*4.*roots/(pow(x2,2)+pow(roots,2))*cos(z*x2)*w[i+1];	//Evaluate k integral at x3
		}
		F_ave = Spectral(Table, roots, a/2.+b/2., z, 2)*4.*roots/(pow(a/2.+b/2.,2)+pow(roots,2))*cos(z*(a/2.+b/2.));
		Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);
		a = b;
	}	//For the bulk of the integral where either the result is well approximated by either the finite or zero width analytic result

	b += stride/2.;	//Evaluate an extra half stride, should result in a result between the extrems of a quarter and three quarter strides
	F_a = F_b = 0;	//Start integration at 0
	for(i = 0; i < 24; i++)
	{
		x1 = (b+a-Disp[i]*(b-a))/2.;	//Actual evaluation points
		x2 = (b+a+Disp[i]*(b-a))/2.;

		F_a += Spectral(Table, roots, x1, z, 2)*4.*roots/(pow(x1,2)+pow(roots,2))*cos(z*x1)*w[i+1];	//Evaluate k integral at x1
		F_b += Spectral(Table, roots, x2, z, 2)*4.*roots/(pow(x2,2)+pow(roots,2))*cos(z*x2)*w[i+1];	//Evaluate k integral at x3
	}
	F_ave = Spectral(Table, roots, a/2.+b/2., z, 2)*4.*roots/(pow(a/2.+b/2.,2)+pow(roots,2))*cos(z*(a/2.+b/2.));
	Answer += (F_a+w[0]*F_ave+F_b)*(b-a)/(2.);//*/

	return(Answer);	//return the best estimate of the integral on the interval*/
}

long double SpatialVac(long double*** Table[], long double roots, long double z, int Temp)
{
	return(2.*M_PI*exp(-roots*z)*Spectral(Table, roots, 0, z, 2));	//return the integral for vacuum from 0 to infinity
}

long double Euclidean(long double*** Table[], long double roots, long double tau, int Temp)
{
	long double T; //T_c = .196GeV = 196MeV

	switch(Temp)
	{
		case 0:
			T = 0;
			break;
		case 1:
			T = .196*1.2;
			break;
		case 2:
			T = .196*1.5;
			break;
		case 3:
			T = .196*2.;
			break;
	}

	if(tau > 1./(2.*T))
		return(0);

	return(cosh(roots*(tau-1./(2.*T)))/sinh(roots/(2.*T))*Spectral(Table, roots, 0, 0, 2)/(2.*M_PI));	//return the integral for vacuum from 0 to infinity
}

complex<long double> atanh(complex<long double> x)
{
	return(complex<long double>(.5,0)*log((complex<long double>(1.,0)+x)/(complex<long double>(1.,0)-x)));
}

long double Spectral_Analytic(long double s, long double p)
{
	long double a = 0.1+(1.09682*exp(-0.0854832*pow(p,2))+0.553167*exp(-0.00242106*pow(p,2)))*List[0];
	long double c = -2.6+(-2.09722*exp(-0.0873686*p)-2.08312*exp(-0.0873686*p))*List[1];
	long double s1 = 3.7249;
	long double s2 = pow(2.85+(0.163437*exp(-0.00115935*pow(p,2))+0.425939*exp(-0.00115935*pow(p,2)))*List[2],2);
	long double b = (c-1.)/(s1-s2);
	long double knee1 = 0.0295;
	long double knee2 = 0.01466;

	long double MJPsi = 3.04031+(0.316713*exp(-0.00197248*pow(p,2))+0.0740995*exp(-0.000182019*pow(p,2)))*List[3];
	long double f1 = a*pow(MJPsi,2)-a*s1+c;
	long double f2 = b*pow(MJPsi,2)-b*s1+c;
	long double f3 = 1.;
	long double f4 = (f1+f2)/2.-(f1-f2)/2.*tanh((pow(MJPsi,2)-s1)/knee1);
	long double f5_A1 = exp((f4+f3)/2.-(f4-f3)/2.*tanh((pow(MJPsi,2)-s2)/knee2));

	f1 = a*s-a*s1+c;
	f2 = b*s-b*s1+c;
	f4 = (f1+f2)/2.-(f1-f2)/2.*tanh((s-s1)/knee1);
	long double f5 = exp((f4+f3)/2.-(f4-f3)/2.*tanh((s-s2)/knee2));

	long double Mq = 1.8;
	long double SigmaRe = 23.744345908864457-1.4154246741829764e6/(2890.3345511938423+pow(244.93882585792915+p,2)*List[4]);
	long double A1 = (9.17+(-1459.29/(369.206+pow(-33.6271+p,2))-218.311/(56.0751+pow(-3.01972+p,2)))*List[5])*f5_A1*MJPsi;
	long double G1 = (0.0026+(0.0140119*exp(-38.3527*pow(p,2))+0.0563881*exp(-0.00556731*pow(p,2)))*List[6])*f5;
	long double A2 = 0.449881+(303.749/(2.67973e-10+pow(81.4093+p,2))-131.725/(2450.75+pow(18.1846+p,2)))*List[7];
	long double G2 = (0.00325+(0.0145782*exp(-0.119134*pow(p,2))+0.00567183*exp(-0.00744749*pow(p,2)))*List[8])*f5;

	/*long double a = 0.1;
	long double c = -2.6;
	long double s1 = 3.7249;
	long double s2 = 8.1225;
	long double b = .818628343;
	long double knee1 = 0.0295;
	long double knee2 = 0.01466;

	long double MJPsi = 3.04031+0.316713*exp(-0.00197248*pow(p,2));
	long double f1 = a*pow(MJPsi,2)-a*s1+c;
	long double f2 = b*pow(MJPsi,2)-b*s1+c;
	long double f3 = 1.;
	long double f4 = (f1+f2)/2.-(f1-f2)/2.*tanh((pow(MJPsi,2)-s1)/knee1);
	long double f5_A1 = exp((f4+f3)/2.-(f4-f3)/2.*tanh((pow(MJPsi,2)-s2)/knee2));

	f1 = a*s-a*s1+c;
	f2 = b*s-b*s1+c;
	f4 = (f1+f2)/2.-(f1-f2)/2.*tanh((s-s1)/knee1);
	long double f5 = exp((f4+f3)/2.-(f4-f3)/2.*tanh((s-s2)/knee2));

	long double Mq = 1.8;
	long double SigmaRe = 23.744345908864457;
	long double A1 = 9.17*f5_A1*MJPsi;
	long double G1 = 0.0026*f5;
	long double A2 = 0.449881;
	long double G2 = 0.00325*f5;*/

	long double sigmaBound = (A1*sqrt(pow(MJPsi,2)/(pow(MJPsi,2)+pow(p,2))*(s+pow(p,2)))*G1)/(M_PI*(pow(s-pow(MJPsi,2),2)+pow(MJPsi,2)/(pow(MJPsi,2)+pow(p,2))*(s+pow(p,2))*pow(G1,2)));
	long double sigmaNonFit = A2*(SigmaRe+s)/2.*imag(sqrt(complex<long double>(s-pow(2.*Mq,2),4*sqrt(pow(2.*Mq,2)/(pow(2.*Mq,2)+pow(p,2))*(s+pow(p,2)))*G2)/complex<long double>(s,4*sqrt(pow(2.*Mq,2)/(pow(2.*Mq,2)+pow(p,2))*(s+pow(p,2)))*G2))*atanh(sqrt(complex<long double>(s,4*sqrt(pow(2.*Mq,2)/(pow(2.*Mq,2)+pow(p,2))*(s+pow(p,2)))*G2)/complex<long double>(s-pow(2.*Mq,2),4*sqrt(pow(2.*Mq,2)/(pow(2.*Mq,2)+pow(p,2))*(s+pow(p,2)))*G2))));

	/*long double a = 2.16;
	long double c = -4.84;
	long double s1 = 3.7249;
	long double s2 = 8.1225;
	long double b = (c-1.)/(s1-s2);
	long double knee1 = 0.0295;
	long double knee2 = 0.01466;

	long double f1 = a*s-a*s1+c;
	long double f2 = b*s-b*s1+c;
	long double f3 = 1.;
	long double f4 = (f1+f2)/2.-(f1-f2)/2.*tanh((s-s1)/knee1);
	long double f5 = exp((f4+f3)/2.-(f4-f3)/2.*tanh((s-s2)/knee2));

	long double Mq = 1.8;
	long double MJPsi = 3.04031;
	long double SigmaRe = 0;
	long double A1 = 67;
	long double G1 = .0026*f5;
	long double A2 = .469224;
	long double G2 = .00325*f5;

	long double sigmaBound = (A1*sqrt(s)*G1)/(M_PI*(pow(s-pow(MJPsi,2),2)+s*pow(G1,2)));
	long double sigmaNonFit = A2*(SigmaRe+s)/2.*imag(sqrt(complex<long double>(s-pow(2.*Mq,2),4*sqrt(s)*G2)/complex<long double>(s,4*sqrt(s)*G2))*atanh(sqrt(complex<long double>(s,4*sqrt(s)*G2)/complex<long double>(s-pow(2.*Mq,2),4*sqrt(s)*G2))));*/

	return(sigmaBound-sigmaNonFit);
}

long double Spectral(long double*** Table[], long double roots, long double p, long double z, int Specify)
{
	if(Specify == 2)
		return(Spectral_Analytic(pow(roots,2),p));
	else
		return(Spectral_Analytic(roots,p));
	/*if(Specify == 2)
		roots *= roots;
	return(-18./M_PI*(Bicubic(Table,roots,p,z,Specify)-26.796184939764153*imag(pow(complex<long double>(Bicubic(Table,roots,p,z,Specify+9),Bicubic(Table,roots,p,z,Specify+3)),2)/complex<long double>(1.-Bicubic(Table,roots,p,z,Specify+12),-Bicubic(Table,roots,p,z,Specify+6)))));*/
}

long double Bicubic(long double*** Table[], long double s, long double p, long double z, int Specify)
{
	long double Interpolation;
	long double t, u;
	int i, j;

	if(z != 0 && p > ((long double)(int(550*z/(2.*M_PI)))-.25)*2.*M_PI/z && Specify%3==2)
		p = ((long double)(int(550*z/(2.*M_PI)))-.25)*2.*M_PI/z;
	else if(z != 0 && p > ((long double)(int((787.2-j/10.)*z/(2.*M_PI)))-.25)*2.*M_PI/z && Specify%3!=2)
		p = ((long double)(int((787.2-j/10.)*z/(2.*M_PI)))-.25)*2.*M_PI/z;

	t = p/.8;	//returns the p index with the fractional part
	i = t;		//p index in the Table
	t -= i;		//Removes the index leaving the fractional part
	if(p >= 600.8)	//Resolves P=600.8, which is at the upper edge of the P Table
	{
		i--;
		t++;
	}	//This only works if s >= 0, this will be chucked if it isn't

	if(s <= 552.25 && s >= 0)
	{
		if(sqrt(s) < 3)	//These are to give the fractional distance from one sqrt(s)-point to the next+the index
			u = sqrt(s)/.1;
		else if(sqrt(s) < 5)
			u = 30+(sqrt(s)-3)/.01;
		else
			u = 230+(sqrt(s)-5)/.1;
		j = u;	//sqrt(s) index in the Table
		u -= j;
	}
	else if(s > 552.25)
	{
		long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration, something cleaver for the Gauss-Laugare region
		s = s-552.25;
		j = 0;
		u = 1;	//When above 23.5, I should only ask for points that I have. I only need those points in one end or the other. It is most likely they'll be in the upper end the way I've written it.
		while(j < 49 && !((abs(s-GaussLa[j+1]) < .00001) || (abs(s-GaussLa[j]) < .00001)))
			j++;
		if(abs(s-GaussLa[j]) < .00001) //Written for case where the point is in the lower end
			u = 0;
		j += 416;
	}
	else if(s<432.64-pow(p,2))	//Resolve t, u, i, and j for s < 0;
	{
		Specify -= 2;
		t = 10.*sqrt(abs(s+pow(p,2)));
		u = 10.*(p-sqrt(abs(s+pow(p,2))));
		i = t;
		t -= i;
		j = u;
		u -= j;
	}
	else
	{
		Specify -= 1;
		t = 187.2+sqrt(s+pow(p,2));
		u = 10.*(p-sqrt(s+pow(p,2)));
		i = t;
		if(i < 208)
			i = 208;
		t -= i;
		i -= 208;	//Remove offset from larger table that doesn't care so much about derivatives
		j = u;
		u -= j;
	}

	/*if(Table[Specify][j][i][4] == 0 && Table[Specify][j][i+1][4] == 0 && Table[Specify][j+1][i][4] == 0 && Table[Specify][j+1][i+1][4] == 0)	//If any of the required points have been invalidated, calculate points from integrals.
	{
		Interpolation = Table[Specify][j][i][0]*(1.-u)*(1.-t)+Table[Specify][j+1][i][0]*u*(1.-t)+Table[Specify][j][i+1][0]*(1.-u)*t+Table[Specify][j+1][i+1][0]*u*t;	//A bilinear interpolation
	}
	else
	{*/
		long double f[16] = {Table[Specify][j][i][0], Table[Specify][j][i+1][0], Table[Specify][j+1][i][0], Table[Specify][j+1][i+1][0], Table[Specify][j][i][1], Table[Specify][j][i+1][1], Table[Specify][j+1][i][1], Table[Specify][j+1][i+1][1], Table[Specify][j][i][2], Table[Specify][j][i+1][2], Table[Specify][j+1][i][2], Table[Specify][j+1][i+1][2], Table[Specify][j][i][3], Table[Specify][j][i+1][3], Table[Specify][j+1][i][3], Table[Specify][j+1][i+1][3]};	//fetch the data points and store them
	if(sqrt(s) >= 3 && sqrt(s) <= 5)	//resolve the width of the sqrt(s) to keep things more correct over trying to resolve it with the table load
		for(int k = 8; k < 16; k++)
			f[k] *= .01;
	else if(s > 0);
		for(int k = 8; k < 16; k++)
			f[k] *= .1;

		long double a[16] = {	f[0],	//Calculate the coeffecients of the function
					f[4],
					-3.*f[0]+3.*f[1]-2.*f[4]-f[5],
					2.*f[0]-2.*f[1]+f[4]+f[5],
					f[8],
					f[12],
					-3.*f[8]+3.*f[9]-2.*f[12]-f[13],
					2.*f[8]-2.*f[9]+f[12]+f[13],
					-3.*f[0]+3.*f[2]-2.*f[8]-f[10],
					-3.*f[4]+3.*f[6]-2.*f[12]-f[14],
					9.*f[0]-9.*f[1]-9.*f[2]+9.*f[3]+6.*f[4]+3.*f[5]-6.*f[6]-3.*f[7]+6.*f[8]-6.*f[9]+3.*f[10]-3.*f[11]+4.*f[12]+2.*f[13]+2.*f[14]+f[15],
					-6.*f[0]+6.*f[1]+6.*f[2]-6.*f[3]-3.*f[4]-3.*f[5]+3.*f[6]+3.*f[7]-4.*f[8]+4.*f[9]-2.*f[10]+2.*f[11]-2.*f[12]-2.*f[13]-f[14]-f[15],
					2.*f[0]-2.*f[2]+f[8]+f[10],
					2.*f[4]-2.*f[6]+f[12]+f[14],
					-6.*f[0]+6.*f[1]+6.*f[2]-6.*f[3]-4.*f[4]-2.*f[5]+4.*f[6]+2.*f[7]-3.*f[8]+3.*f[9]-3.*f[10]+3.*f[11]-2.*f[12]-f[13]-2.*f[14]-f[15],
					4.*f[0]-4.*f[1]-4.*f[2]+4.*f[3]+2.*f[4]+2.*f[5]-2.*f[6]-2.*f[7]+2.*f[8]-2.*f[9]+2.*f[10]-2.*f[11]+f[12]+f[13]+f[14]+f[15]};

		Interpolation = a[0]+a[1]*t+a[2]*pow(t,2)+a[3]*pow(t,3)+a[4]*u+a[5]*t*u+a[6]*pow(t,2)*u+a[7]*pow(t,3)*u+a[8]*pow(u,2)+a[9]*t*pow(u,2)+a[10]*pow(t*u,2)+a[11]*pow(t,3)*pow(u,2)+a[12]*pow(u,3)+a[13]*t*pow(u,3)+a[14]*pow(t,2)*pow(u,3)+a[15]*pow(t*u,3);	//Output the function, bicubic interpolation
	//}

	return(Interpolation);
}

void Validate(long double*** Table[], int M[], int N[])
{
	const long double E[] = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.540308, 2.545308, 2.550308, 2.555308, 2.560308, 2.565308, 2.570308, 2.575308, 2.580308, 2.585308, 2.590308, 2.595308, 2.600308, 2.605308, 2.610308, 2.615308, 2.620308, 2.625308, 2.630308, 2.635308, 2.640308, 2.645308, 2.650308, 2.655308, 2.660308, 2.665308, 2.670308, 2.675308, 2.680308, 2.685308, 2.690308, 2.695308, 2.700308, 2.705308, 2.710308, 2.715308, 2.720308, 2.725308, 2.730308, 2.735308, 2.740308, 2.745308, 2.750308, 2.755308, 2.760308, 2.765308, 2.770308, 2.775308, 2.780308, 2.785308, 2.790308, 2.795308, 2.800308, 2.805308, 2.810308, 2.815308, 2.820308, 2.825308, 2.830308, 2.835308, 2.840308, 2.845308, 2.850308, 2.855308, 2.860308, 2.865308, 2.870308, 2.875308, 2.880308, 2.885308, 2.890308, 2.895308, 2.900308, 2.905308, 2.910308, 2.915308, 2.920308, 2.925308, 2.930308, 2.935308, 2.940308, 2.945308, 2.950308, 2.955308, 2.960308, 2.965308, 2.970308, 2.975308, 2.980308, 2.985308, 2.990308, 2.995308, 3.000308, 3.005308, 3.010308, 3.015308, 3.020308, 3.025308, 3.030308, 3.035308, 3.040308, 3.045308, 3.050308, 3.055308, 3.060308, 3.065308, 3.070308, 3.075308, 3.080308, 3.085308, 3.090308, 3.095308, 3.100308, 3.105308, 3.110308, 3.115308, 3.120308, 3.125308, 3.130308, 3.135308, 3.140308, 3.145308, 3.150308, 3.155308, 3.160308, 3.165308, 3.170308, 3.175308, 3.180308, 3.185308, 3.190308, 3.195308, 3.200308, 3.205308, 3.210308, 3.215308, 3.220308, 3.225308, 3.230308, 3.235308, 3.240308, 3.245308, 3.250308, 3.255308, 3.260308, 3.265308, 3.270308, 3.275308, 3.280308, 3.285308, 3.290308, 3.295308, 3.300308, 3.305308, 3.310308, 3.315308, 3.320308, 3.325308, 3.330308, 3.335308, 3.340308, 3.345308, 3.350308, 3.355308, 3.360308, 3.365308, 3.370308, 3.375308, 3.380308, 3.385308, 3.390308, 3.395308, 3.400308, 3.405308, 3.410308, 3.415308, 3.420308, 3.425308, 3.430308, 3.435308, 3.440308, 3.445308, 3.450308, 3.455308, 3.460308, 3.465308, 3.470308, 3.475308, 3.480308, 3.485308, 3.490308, 3.495308, 3.500308, 3.505308, 3.510308, 3.515308, 3.520308, 3.525308, 3.530308, 3.535308, 3.540308, 3.55, 3.56375, 3.5775, 3.59125, 3.605, 3.61875, 3.6325, 3.64625, 3.66, 3.67375, 3.6875,  3.70125, 3.715, 3.72875, 3.7425, 3.75625, 3.77, 3.78375, 3.7975, 3.81125, 3.825, 3.83875, 3.8525, 3.86625, 3.88, 3.89375, 3.9075, 3.92125, 3.935, 3.94875, 3.9625, 3.97625, 3.99, 4.00375, 4.0175, 4.03125, 4.045, 4.05875, 4.0725, 4.08625, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21., 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22., 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23., 23.1, 23.2, 23.3, 23.4, 23.5};
	int i,j,k;
	long double Test;
	//long double Average;

	for(k = 0; k < 1; k++)
	{
		for(j = 0; j < N[k]-1; j++)
		{
			for(i = 0; i < M[k]-1; i++)
			{
				Test = Spectral(Table, j+.5, i, -1., k);	//Checks the point on the mid-point E
				//Average = (Spectral(Table, j, i, -1., k)+Spectral(Table, j+1, i, -1., k))/2.;
				if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
				{
					Table[k][i][j][4] = 0;	//Invalidates the data around the point
					Table[k][i+1][j][4] = 0;
					Table[k][i][j+1][4] = 0;
					Table[k][i+1][j+1][4] = 0;
				}

				Test = Spectral(Table, j, i+.5, -1., k);	//Checks the point on the mid-point i
				//Average = (Spectral(Table, j, i, -1., k)+Spectral(Table, j, i+1., -1., k))/2.;
				if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
				{
					Table[k][i][j][4] = 0;	//Invalidates the data around the point
					Table[k][i+1][j][4] = 0;
					Table[k][i][j+1][4] = 0;
					Table[k][i+1][j+1][4] = 0;
				}

				Test = Spectral(Table, j+.5, i+.5, -1., k);	//Checks the point on the mid-point E,i
				//Average = (Spectral(Table, j, i, -1., k)+Spectral(Table, j+1, i, -1., k)+Spectral(Table, j, i+1., -1., k)+Spectral(Table, j+1, i+1., -1., k))/4.;
				if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
				{
					Table[k][i][j][4] = 0;	//Invalidates the data around the point
					Table[k][i+1][j][4] = 0;
					Table[k][i][j+1][4] = 0;
					Table[k][i+1][j+1][4] = 0;
				}
			}
		}
	}
	for(j = 0; j < M[1]-2; j++)
	{
		for(long double P = 0; P < .8*N[1]-1.6; P += .8)
		{
			i = P/.8;	//returns the p index without the fractional part

			Test = Spectral(Table, (E[j]+E[j+1])/2., P, -1., 1);	//Checks the point on the mid-point E
			//Average = (Spectral(Table, E[j], P, -1., 2)+Spectral(Table, E[j+1], P, -1., 2))/2.;
			if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
			{
				Table[1][j][i][4] = 0;	//Invalidates the data around the point
				Table[1][j][i+1][4] = 0;
				Table[1][j+1][i][4] = 0;
				Table[1][j+1][i+1][4] = 0;
			}

			Test = Spectral(Table, E[j], P+.4, -1., 1);	//Checks the point on the mid-point P
			//Average = (Spectral(Table, E[j], P, -1., 2)+Spectral(Table, E[j], P+.8, -1., 2))/2.;
			if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
			{
				Table[1][j][i][4] = 0;	//Invalidates the data around the point
				Table[1][j][i+1][4] = 0;
				Table[1][j+1][i][4] = 0;
				Table[1][j+1][i+1][4] = 0;
			}

			Test = Spectral(Table, (E[j]+E[j+1])/2., P+.4, -1., 1);	//Checks the point on the mid-point E,P
			//Average = (Spectral(Table, E[j], P, -1., 2)+Spectral(Table, E[j+1], P, -1., 2)+Spectral(Table, E[j], P+.8, -1., 2)+Spectral(Table, E[j+1], P+.8, -1., 2))/4.;
			if(Test < 0)// || abs(Average-Test)/Average < 1.15)	//The Spectral function must be positive and resonably close to the linear interpolation
			{
				Table[1][j][i][4] = 0;	//Invalidates the data around the point
				Table[1][j][i+1][4] = 0;
				Table[1][j+1][i][4] = 0;
				Table[1][j+1][i+1][4] = 0;
			}
		}
	}
}

bool ReadIn(long double*** Table[], int N[], int M[], char* FileReadIn)
{
	char* FileName1 = new char[150];
	char* FileName2 = new char[150];
	char* FileName3 = new char[150];
	char* FileName4 = new char[150];
	char* FileName5 = new char[150];
	strcpy(FileName1,FileReadIn);
	strcpy(FileName2,FileReadIn);
	strcpy(FileName3,FileReadIn);
	strcpy(FileName4,FileReadIn);
	strcpy(FileName5,FileReadIn);
	ifstream File1(strcat(FileName1,"ImG"));
	ifstream File2(strcat(FileName2,"ImGV1"));
	ifstream File3(strcat(FileName3,"ImGV2"));
	ifstream File4(strcat(FileName4,"ReGV1"));
	ifstream File5(strcat(FileName5,"ReGV2"));
	int i,j,k,m;	//Counters
	long double holder;

	if(File1.is_open() == false || File2.is_open() == false || File3.is_open() == false || File4.is_open() == false || File5.is_open() == false)
		return(false);

	for(m = 0; m < 3; m++)	//Table count
	{
		for(k = 0; k < M[m]; k++)	//i/P count
		{
			for(j = 0; j < N[m]; j++)	//j/sqrt(s) count
			{
				for(i = 0; i < 4; i++)	//Value type count
				{
					File1 >> Table[m][j][k][i];
					File2 >> Table[m+3][j][k][i];
					File3 >> Table[m+6][j][k][i];
					if(j < N[m+9])
					{
						File4 >> Table[m+9][j][k][i];
						File5 >> Table[m+12][j][k][i];
					}
					if(m == 2)
					{
						if(i == 1 || i == 3)
						{
							Table[m][j][k][i] *= .8;	//*=dP/di
							Table[m+3][j][k][i] *= .8;
							Table[m+6][j][k][i] *= .8;
							if(j < N[m+9])
							{
								Table[m+9][j][k][i] *= .8;
								Table[m+12][j][k][i] *= .8;
							}
						}
					}
				}
				Table[m][j][k][4] = 1;
				Table[m+3][j][k][4] = 1;
				Table[m+6][j][k][4] = 1;
				if(j < N[m+9])
				{
					Table[m+9][j][k][4] = 1;
					Table[m+12][j][k][4] = 1;
				}
			}
		}
	}

	return(true);
}

void Init(long double*** Table[], int N[], int M[])
{
	int i, j, k;

	for(k = 0; k < 15; k++)
	{
		Table[k] = new long double**[N[k]];
		for(i = 0; i < N[k]; i++)
		{
			Table[k][i] = new long double*[M[k]];
			for(j = 0; j < M[k]; j++)
				Table[k][i][j] = new long double[5];
		}
	}

	return;
}
