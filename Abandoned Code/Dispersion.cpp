//This program is to calculate the dispersion relation, yielding 2 real parts.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cstring>
#include<cstdlib>
using namespace std;

bool ReadIn(long double***[], int[], int[], char*);	//Read in from file the contents of the matrix. It includes 2 offsets for multiple file read in and for padding due to negative momentum not being in the file
void Init(long double***[], int[], int[]);	//Initialize the table to the correct size
long double Imaginary(long double***[], long double, long double);	//The tabulated points and 2 inputs and returns the bicubic interpolation of that input
long double Real(long double***[], long double, long double);	//Returns the real part from the imaginary part by dispersion relation

char* Process;

int main(int argc, char* argv[])	//Process#, # of Process, output file name, Input file name, starting point
{
	char* File = new char[30];	//Name of the file
	strcpy(File, argv[3]);
	Process = argv[1];
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot;
	if(atoi(argv[5]) == 0)	//If starting from the beginning, overwrite
		TPlot.open(File);
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);
	TPlot << setprecision(18);
	long double*** Table[3];	//The table of values computed by Spectral
	long double holder[616];
	int N[3] = {151,151,465}, M[3] = {209,581,752};	//The size of the table
	const int iProcess = atoi(argv[1]);
	const int Total = atoi(argv[2]);
	int i, j;
	long double s, P;

	Init(Table, N, M);
	if(!ReadIn(Table, N, M, argv[4]))
		return(0);

	for(i = atoi(argv[5]); i <= 788; i++)	//Argv[5] allows to restart where ever
	{
		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					P = i/10.+j/10.;
					s = -i*j/50.-pow(j/10.,2);
				}
				else
				{
					P = i+j/10.-187.2;
					s = -j/5.*(i-187.2)-pow(j/10.,2);
				}
			}
			else
			{
				P = i*.8;
#ifndef BB
				if(j <= 181)
					s = pow((j-151.)/10.,2);
				else if(j <= 381)
					s = pow((j-181.)/100.+3.,2);
				else
					s = pow((j-381.)/10.+5.,2);

#else
				if(j <= 251)
					s = pow((j-151.)/10.,2);
				else if(j <= 451)
					s = pow((j-251.)/100.+10.,2);
				else
					s = pow((j-451.)/10.+12.,2);
#endif
			}


			holder[j] = Real(Table, s, P);
		}

		for(j = iProcess; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					P = i/10.+j/10.;
					s = -i*j/50.-pow(j/10.,2);
				}
				else
				{
					P = i+j/10.-187.2;
					s = -j/5.*(i-187.2)-pow(j/10.,2);
				}
			}
			else
			{
				P = i*.8;
#ifndef BB
				if(j <= 181)
					s = pow((j-151.)/10.,2);
				else if(j <= 381)
					s = pow((j-181.)/100.+3.,2);
				else
					s = pow((j-381.)/10.+5.,2);
#else
				if(j <= 251)
					s = pow((j-151.)/10.,2);
				else if(j <= 451)
					s = pow((j-251.)/100.+10.,2);
				else
					s = pow((j-451.)/10.+12.,2);
#endif
			}
			TPlot << i << " " << j << " " << P << " " << s << " " << holder[j] << endl;
		}
	}

	TPlot.close();//*/

	return(0);
}

long double Real(long double*** Table[], long double s, long double p)
{
	long double Disp[] = {0.1603586456402253758680961, 0.3165640999636298319901173, 0.4645707413759609457172671, 0.6005453046616810234696382, 0.7209661773352293786170959, 0.8227146565371428249789225, 0.9031559036148179016426609, 0.9602081521348300308527788, 0.9924068438435844031890177}; //Displacement from center for 35th order Gauss-Legendre integration
	long double w[] = {8589934592./53335593025., 0.1589688433939543476499564, 0.1527660420658596667788554, 0.1426067021736066117757461, 0.1287539625393362276755158, 0.1115666455473339947160239, 0.09149002162244999946446209, 0.06904454273764122658070826, 0.04481422676569960033283816, 0.01946178822972647703631204}; //Weight of the function at Disp
	long double DispLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre quadrature
	long double wLa[] = {0.07496328305102102808055, 0.1745735743605928864303, 0.2745074833881225250022, 0.3747323102655645620060, 0.4753412526072084401161, 0.5764380939967183636147, 0.6781307242364945406823, 0.7805307978511547593175, 0.8837542316062452388883, 0.9879219194279636096671, 1.0931605619330277996916, 1.1996035979670979427973, 1.3073922479469277349326, 1.416676687469297701993, 1.5276173754408796787012, 1.640386566702889623924, 1.7551700457872174635214, 1.8721691266543402861779, 1.9916029736088098866132, 2.1137113117669909276048, 2.2387576123844772725684, 2.3670328602831611098048, 2.4988600392644108123394, 2.6345995091430390709, 2.7746554982525006307172, 2.9194840027576204632431, 3.0696024758091833914472, 3.2256018156600758204608, 3.3881613374746331979827, 3.5580676615951707296054, 3.7362388067183244743069, 3.9237552950635210172968, 4.1219008467729629867363, 4.3322164077399479741288, 4.5565730632309056055423, 4.7972722621195591678357, 5.057186469320242487569, 5.3399612774797865633198, 5.6503138450512931300331, 5.9944877492232503537552, 6.3809726096501927329094, 6.8216946862388774056326, 7.3340972531892936469048, 7.9450326451948326187906, 8.6987143462393085933469, 9.6750102652900375180015, 11.039313738067347840094, 13.220456867750092021034, 17.982575250664959108273};	//Weights for 95th order Gauss-Laguerre quadrature
	long double Range[] = {0, 3.24, 12.96, 25, 100, 225, 552.25};
	long double Answer;
	long double F_a, F_b, F_ave;
	long double x1, x2;
	long double a, b;
	int i, j;

	if(p > 600.8 || s > 552.25)	//No meaningful dispersion can be calculuated for P>600.8 as only s < 0 is avalible in this region
		return(0);

	if(p > 15)
	{
		Answer = Imaginary(Table, s, p)*log(abs((552.25-s)/(s+30.*p-225.)));
		a = -15.*(-15.+2.*p);
	}
	else
	{
		Answer = Imaginary(Table, s, p)*log(abs((552.25-s)/(s+pow(p,2))));
		a = -pow(p,2);
	}

	if(abs(s+pow(p,2)) < pow(.001,2)*.1 && p <= 150)	//Take linear limits to subvert issues with the endpoints of log((552.25-s)/(s-a)) causing issues
		return(2.*Real(Table, s+.001, p)-Real(Table, s+.002, p));
	else if(abs(s+30.*p-225.) < pow(.001,2)*.1 && p > 150)
		return(2.*Real(Table, s+.001, p)-Real(Table, s+.002, p));
	else if(abs(s-552.25) < .00001)
		return(-Real(Table, pow(23.498,2), p)+2.*Real(Table, pow(23.499,2), p));

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

			F_a += (Imaginary(Table, x1, p)-Imaginary(Table, s, p))/(x1-s)*w[j+1]; //Evaluate function at x1
			F_b += (Imaginary(Table, x2, p)-Imaginary(Table, s, p))/(x2-s)*w[j+1]; //Evaluate function at x2
		}
		F_ave = (Imaginary(Table, (a+b)/2., p)-Imaginary(Table, s, p))/((a+b)/2.-s)*w[0]; //Evaluate function at (a+b)/2.
		Answer += (F_a+F_ave+F_b)*(b-a)/(2.);
		a = b;
	}

	F_a = 0;
	for(j = 0; j < 49; j++)
	{
		x1 = 552.25+DispLa[j];
		F_a += Imaginary(Table, x1, p)/(x1-s)*wLa[j]; //Evaluate function at x1
	}
	Answer += F_a;

	return(Answer/M_PI);
}

long double Imaginary(long double*** Table[], long double s, long double p)
{
	long double t, u;
	int i, j;
	int Specify = 2;	//The specify varible was to specify between 3 different table. At the moment, I only have one. But I don't know if I'll get a second, so I'll leave as much structure as possible for the possiblity of getting it.

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
		Specify = 0;
		t = 10.*sqrt(abs(s+pow(p,2)));
		u = 10.*(p-sqrt(abs(s+pow(p,2))));
		i = t;
		t -= i;
		j = u;
		u -= j;
	}
	else
	{
		Specify = 1;
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

	if(Specify != 2 && j >= 150)	//Probably correct for such a large j
		return(0);

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

	return(a[0]+a[1]*t+a[2]*pow(t,2)+a[3]*pow(t,3)+a[4]*u+a[5]*t*u+a[6]*pow(t,2)*u+a[7]*pow(t,3)*u+a[8]*pow(u,2)+a[9]*t*pow(u,2)+a[10]*pow(t*u,2)+a[11]*pow(t,3)*pow(u,2)+a[12]*pow(u,3)+a[13]*t*pow(u,3)+a[14]*pow(t,2)*pow(u,3)+a[15]*pow(t*u,3));	//Output the function, bicubic interpolation
}

bool ReadIn(long double*** Table[], int N[], int M[], char* FileReadIn)
{
	ifstream File(FileReadIn);
	int i,j,k,m;	//Counters
	long double Holder;

	if(File.good() == false)
		return(false);

	for(m = 0; m < 3; m++)	//Table count
	{
		for(k = 0; k < M[m]; k++)	//i/P count
		{
			for(j = 0; j < N[m]; j++)	//j/s count
			{
				for(i = 0; i < 4; i++)	//Value type count
				{
					File >> Table[m][j][k][i];
					if(m == 2)
					{
						if(i == 1 || i == 3)
							Table[m][j][k][i] *= .8;	//*=dP/di
					}
				}
			}
		}
	}

	return(true);
}

void Init(long double*** Table[], int N[], int M[])
{
	int i, j, k;

	for(k = 0; k < 3; k++)
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
