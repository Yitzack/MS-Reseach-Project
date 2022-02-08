#include<fstream>
#ifndef SPECTRAL_NON
#define SPECTRAL_NON
using namespace std;

class Spectral_Non{
	public:
		long double Spectral(long double s, long double P);
		long double Spatial(long double z);
		long double Euclidean(long double tau, long double P);
		long double Euclidean(long double tau, long double P, long double T);
		void Print(ostream&); //Print the parameters
		void Add(long double, int);		//Add to a parameter
		void Add(long double, int, int);
		void Replace(long double, int);	//Replace a parameter outright
		void Replace(long double, int, int);
		long double Read(int);			//Read a parameter
		long double Read(int, int);
		void Normal(int, long double[2], long double[2]);

		Spectral_Non(long double[2][3], long double, bool);
		Spectral_Non(long double[6], long double, bool);
		Spectral_Non(long double[2][3], int, bool);
		Spectral_Non(long double[6], int, bool);
	private:
		long double Parameters[2][3];
		long double Temp;
		bool Vacuum;

		long double Spatial_sInt(long double z, long double P);
		long double Spatial_P0Int(long double z, long double P0);

		long double SpatialGeneralKernel(long double s, long double P, long double z);
		long double SpatialCutoffKernel(long double s, long double P0, long double z);
		long double EuclideanKernel(long double s, long double P, long double tau, long double T);

		long double Q(long double, long double, long double, long double);
		long double Uniform();
};

Spectral_Non::Spectral_Non(long double Parm[2][3], long double T, bool Vac)
{
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			Parameters[i][j] = Parm[i][j];
	
	Temp = T;
	Vacuum = Vac;
}

Spectral_Non::Spectral_Non(long double Parm[6], long double T, bool Vac)
{
	for(int i = 0; i < 6; i++)
		Parameters[i/3][i%3] = Parm[i];
	
	Temp = T;
	Vacuum = Vac;
}

Spectral_Non::Spectral_Non(long double Parm[2][3], int T, bool Vac)
{
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			Parameters[i][j] = Parm[i][j];
	
	switch(T)
	{
		case 0:
			Temp = 0;
			break;
		case 1:
			Temp = .194;
			break;
		case 2:
			Temp = .258;
			break;
		case 3:
			Temp = .320;
			break;
		case 4:
			Temp = .400;
			break;
	}

	Vacuum = Vac;
}

Spectral_Non::Spectral_Non(long double Parm[6], int T, bool Vac)
{
	for(int i = 0; i < 6; i++)
		Parameters[i/3][i%3] = Parm[i];
	
	switch(T)
	{
		case 0:
			Temp = 0;
			break;
		case 1:
			Temp = .194;
			break;
		case 2:
			Temp = .258;
			break;
		case 3:
			Temp = .320;
			break;
		case 4:
			Temp = .400;
			break;
	}

	Vacuum = Vac;
}

void Spectral_Non::Normal(int i, long double Range0[2], long double Range1[2])
{
	long double uniform[2] = {Uniform(),Uniform()};
	long double normal[2] = {sqrt(-2.*log(uniform[0]))*cos(2.*M_PI*uniform[1]),sqrt(-2.*log(uniform[0]))*sin(2.*M_PI*uniform[1])};
	long double mu[2] = {Parameters[i][1], Parameters[i][2]};
	long double test[2] = {mu[0]+.1*normal[0],mu[1]+.1*normal[1]};

	while(test[0] < Range0[0] || test[0] > Range0[1] || test[1] < Range1[0] || test[1] > Range1[1])
	{
		uniform[0] = Uniform();
		uniform[1] = Uniform();
		normal[0] = sqrt(-2.*log(uniform[0]))*cos(2.*M_PI*uniform[1]);
		normal[1] = sqrt(-2.*log(uniform[0]))*sin(2.*M_PI*uniform[1]);
		test[0] = mu[0]+.1*normal[0];
		test[1] = mu[1]+.1*normal[1];
	}

	Parameters[i][1] = test[0];
	Parameters[i][2] = test[1];
}

void Spectral_Non::Print(ostream& Stream)
{
	for(int i = 0; i < 2; i++)
		for(int j = 1; j < 3; j++)
		{
			Stream << Parameters[i][j];
			if(!(i == 1 && j == 2))
				Stream << ",";
		}
	Stream << flush;
}

void Spectral_Non::Add(long double Value, int i)		//Add to a parameter
{
	Parameters[i/3][i%3] += Value;
}

void Spectral_Non::Add(long double Value, int i, int j)
{
	Parameters[i][j] += Value;
}

void Spectral_Non::Replace(long double Value, int i)	//Replace a parameter outright
{
	Parameters[i/3][i%3] = Value;
}

void Spectral_Non::Replace(long double Value, int i, int j)
{
	Parameters[i][j] = Value;
}

long double Spectral_Non::Read(int i)			//Read a parameter
{
	return(Parameters[i/3][i%3]);
}

long double Spectral_Non::Read(int i, int j)
{
	return(Parameters[i][j]);
}

long double Spectral_Non::Spatial(long double z)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	/*long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double Disp[] = {0.27963041316178305, 0.538469310105683, 0.7541667265708494, 0.9061798459386639, 0.9840853600948425};	//Displacement from center for unknown order Gauss-Kronrod integration from Mathematica
	long double w[] = {0.2829874178574912, 0.2728498019125589, 0.24104033922864776, 0.1868007965564926, 0.11523331662247445, 0.0425820367510818};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = 187.*M_PI/(2.*z);	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	long double Intervals = 2.*M_PI/z;
	int i, j, l;		//Counting varibles

	a = 0;
	b = M_PI/(2.*z);

	do
	{
		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 5; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;
			F += w[l+1]*Spatial_sInt(z, x1);
			F += w[l+1]*Spatial_sInt(z, x2);
		}
		F += w[0]*Spatial_sInt(z, a/2.+b/2.);

		Answer += F*(b-a)/2.;
		a = b;
		b += Intervals;
	}while(a < Max);
	Answer += Spatial_P0Int(z, Max);

	return(Answer);
}

long double Spectral_Non::Spatial_sInt(long double z, long double P)
{
	long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	/*long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = 552.25;	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int Intervals;
	int i = 1, j, l;		//Counting varibles
	long double Stops[15] = {0, 0, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 552.25};

	Stops[0] = pow(2.*Q(P, Parameters[0][0], Parameters[0][1], Parameters[0][2]),2);
	Stops[1] = Stops[0]/2.+Stops[2]/2.;
	Intervals = 14;

	a = Stops[0];
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 24; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*SpatialGeneralKernel(x1, P, z)*Spectral(x1, P);
			F += w[l+1]*SpatialGeneralKernel(x2, P, z)*Spectral(x2, P);
		}
		F += w[0]*SpatialGeneralKernel(a/2.+b/2., P, z)*Spectral(a/2.+b/2., P);

		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Spectral_Non::Spatial_P0Int(long double z, long double P0)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = 552.25;	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int Intervals;
	int i = 1, j, l;		//Counting varibles
	long double Stops[15] = {0, 0, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 552.25};

	Stops[0] = pow(2.*Q(P0, Parameters[0][0], Parameters[0][1], Parameters[0][2]),2);
	Stops[1] = Stops[0]/2.+Stops[2]/2.;
	Intervals = 14;

	a = Stops[0];
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 9; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*SpatialCutoffKernel(x1, P0, z)*Spectral(x1, P0);
			F += w[l+1]*SpatialCutoffKernel(x2, P0, z)*Spectral(x2, P0);
		}
		F += w[0]*SpatialCutoffKernel(a/2.+b/2., P0, z)*Spectral(a/2.+b/2., P0);
		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Spectral_Non::Euclidean(long double tau, long double P, long double T)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = pow(400.,2);	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int i = 1, j, l;		//Counting varibles
	long double Stops[44] = {0, 0, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 604, 704, 804, 904, 1004, 2004, 3004, 4004, 5004, 6004, 7004, 8004, 9004, 10004, 10004, 20004, 30004, 40004, 50004, 60004, 70004, 80004, 90004, 100004, 110004, 120004, 130004, 140004, 150004, 160000};

	Stops[0] = pow(2.*Q(P, Parameters[0][0], Parameters[0][1], Parameters[0][2]),2);
	Stops[1] = Stops[0]/2.+Stops[2]/2.;

	tau /= T;

	a = Stops[0];
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 9; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*EuclideanKernel(x1, P, tau, T)*Spectral(x1, P);
			F += w[l+1]*EuclideanKernel(x2, P, tau, T)*Spectral(x2, P);
		}
		F += w[0]*EuclideanKernel(a/2.+b/2., P, tau, T)*Spectral(a/2.+b/2., P);

		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Spectral_Non::Euclidean(long double tau, long double P)
{
	/*long double Disp[] = {0.06342068498268678602883, 0.1265859972696720510680, 0.1892415924618135864853, 0.2511351786125772735072, 0.3120175321197487622079, 0.3716435012622848888637, 0.4297729933415765246586, 0.4861719414524920421770, 0.5406132469917260665582, 0.5928776941089007124559, 0.6427548324192376640569, 0.6900438244251321135048, 0.7345542542374026962137, 0.7761068943454466350181, 0.8145344273598554315395, 0.8496821198441657010349, 0.8814084455730089100370, 0.9095856558280732852130, 0.9341002947558101490590, 0.9548536586741372335552, 0.9717622009015553801400, 0.9847578959142130043593, 0.9937886619441677907601, 0.9988201506066353793618};	//Displacement from center for 97th order Gauss-Legendre integration
	long double w[] = {0.06346328140479059771825, 0.06333550929649174859084, 0.06295270746519569947440, 0.06231641732005726740108, 0.06142920097919293629683, 0.06029463095315201730311, 0.05891727576002726602453, 0.05730268153018747548516, 0.05545734967480358869043, 0.05338871070825896852794, 0.05110509433014459067462, 0.04861569588782824027765, 0.04593053935559585354250, 0.04306043698125959798835, 0.04001694576637302136861, 0.03681232096300068981947, 0.03345946679162217434249, 0.02997188462058382535069, 0.02636361892706601696095, 0.02264920158744667649877, 0.01884359585308945844445, 0.01496214493562465102958, 0.01102055103159358049751, 0.007035099590086451473451, 0.003027278988922905077481};	//Weight of the function at Disp*/
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177};	//Displacement from center for 37th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204};	//Weight of the function at Disp*/
	long double a, b;	//Sub-interval limits of integration
	long double Max = pow(400.,2);	//Upper limit of integration
	long double F;	//Sum of ordinates*weights
	long double Answer = 0;	//Results to be returned
	long double x1, x2;	//Abscissa
	long double holder;
	pair<long double, long double> zero;
	int i = 1, j, l;		//Counting varibles
	long double Stops[44] = {0, 0, 15, 18, 21, 24, 34, 44, 54, 104, 204, 304, 404, 504, 604, 704, 804, 904, 1004, 2004, 3004, 4004, 5004, 6004, 7004, 8004, 9004, 10004, 10004, 20004, 30004, 40004, 50004, 60004, 70004, 80004, 90004, 100004, 110004, 120004, 130004, 140004, 150004, 160000};

	Stops[0] = pow(2.*Q(P, Parameters[0][0], Parameters[0][1], Parameters[0][2]),2);
	Stops[1] = Stops[0]/2.+Stops[2]/2.;

	tau /= Temp;

	a = Stops[0];
	do
	{
		b = Stops[i];
		i++;

		if(b > Max)
			b = Max;

		F = 0;

		for(l = 0; l < 9; l++) //Integrate the sub-interval
		{
			x1 = (b+a-Disp[l]*(b-a))/2.;
			x2 = (b+a+Disp[l]*(b-a))/2.;

			F += w[l+1]*EuclideanKernel(x1, P, tau, Temp)*Spectral(x1, P);
			F += w[l+1]*EuclideanKernel(x2, P, tau, Temp)*Spectral(x2, P);
		}
		F += w[0]*EuclideanKernel(a/2.+b/2., P, tau, Temp)*Spectral(a/2.+b/2., P);

		Answer += F*(b-a)/2.;
		a = b;
	}while(a < Max);

	return(Answer);
}

long double Spectral_Non::SpatialGeneralKernel(long double s, long double P, long double z)
{
	return(cos(P*z)/(s+pow(P,2)));
}

long double Spectral_Non::SpatialCutoffKernel(long double s, long double P0, long double z)
{
	return(-((z*(-2.*P0*z*(12852.-5868.*pow(sqrt(s)*z,2)+1737.*pow(sqrt(s)*z,4)-78.*pow(sqrt(s)*z,6)+pow(P0*z,8)+pow(sqrt(s)*z,8)+3.*pow(P0*z,4)*(579.+26.*pow(sqrt(s)*z,2)+2.*pow(sqrt(s)*z,4))+pow(P0,6)*(78.*pow(z,6)+4.*s*pow(z,8))+pow(P0,2)*(5868.*pow(z,2)+2418.*s*pow(z,4)-78.*pow(s,2)*pow(z,6)+4.*pow(s,3)*pow(z,8)))*cos(P0*z)+(-5400.+33948.*pow(sqrt(s)*z,2)-16074.*pow(sqrt(s)*z,4)+2301.*pow(sqrt(s)*z,6)-88.*pow(sqrt(s)*z,8)+pow(P0*z,10)+pow(sqrt(s)*z,10)+pow(P0*z,8)*(84.+5.*pow(sqrt(s)*z,2))+pow(P0*z,6)*(2037.+164.*pow(sqrt(s)*z,2)+10.*pow(sqrt(s)*z,4))+pow(P0*z,4)*(10782.+3975.*pow(sqrt(s)*z,2)-12.*pow(sqrt(s)*z,4)+10.*pow(sqrt(s)*z,6))+pow(P0*z,2)*(27756.-4716.*pow(sqrt(s)*z,2)+4239.*pow(sqrt(s)*z,4)-180.*pow(sqrt(s)*z,6)+5.*pow(sqrt(s)*z,8)))*sin(P0*z)))/(pow(P0*z,12)+6.*pow(P0*z,10)*(15.+pow(sqrt(s)*z,2))+3.*pow(P0*z,8)*(819.+90.*pow(sqrt(s)*z,2)+5.*pow(sqrt(s)*z,4))+pow(-36.+216.*pow(sqrt(s)*z,2)-45.*pow(sqrt(s)*z,4)+pow(sqrt(s)*z,6),2)+4.*pow(P0*z,6)*(4878.+1593.*pow(sqrt(s)*z,2)+45.*pow(sqrt(s)*z,4)+5.*pow(sqrt(s)*z,6))+3.*pow(P0*z,4)*(16632.+6120.*pow(sqrt(s)*z,2)+2610.*pow(sqrt(s)*z,4)-60.*pow(sqrt(s)*z,6)+5.*pow(sqrt(s)*z,8))+6.*pow(P0*z,2)*(2592.+12312.*pow(sqrt(s)*z,2)-3060.*pow(sqrt(s)*z,4)+1062.*pow(sqrt(s)*z,6)-45.*pow(sqrt(s)*z,8)+pow(sqrt(s)*z,10)))));
}

long double Spectral_Non::EuclideanKernel(long double s, long double P, long double tau, long double T)
{
	return(cosh(sqrt(s+pow(P,2))*(tau-1./(2.*T)))/(2.*sqrt(s+pow(P,2))*sinh(sqrt(s+pow(P,2))/(2.*T))));
}

long double Spectral_Non::Spectral(long double s, long double P)
{
	if(Vacuum)
		P = 0;

	static long double old_P = P;
	static long double M = Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
	static long double n = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);

	if(P != old_P)
	{
		old_P = P;
		M = Q(P,Parameters[0][0],Parameters[0][1],Parameters[0][2]);
		n = Q(P,Parameters[1][0],Parameters[1][1],Parameters[1][2]);
	}

	if(s < pow(2.*M,2))
		return(0);
	return((3.*s)/(8.*pow(M_PI,2))*sqrt(1.-pow(pow(2.*M,2)/s,n)));
}

long double Spectral_Non::Q(long double P, long double Q0, long double QV, long double P0)
{
	return((Q0*pow(P0,2)+QV*pow(P,2))/(pow(P0,2)+pow(P,2)));
}

long double Spectral_Non::Uniform()
{
	return((long double)(rand())/(long double)(RAND_MAX));
}

#endif
