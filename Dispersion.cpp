//This program is to calculate the dispersion relation, yielding 2 real parts.
#include<iostream>
#include<fstream>
#inclde<cmath>
using namespace std;

bool ReadIn(long double***[], int[], int[], char*);	//Read in from file the contents of the matrix. It includes 2 offsets for multiple file read in and for padding due to negative momentum not being in the file
void Init(long double***[], int[], int[]);	//Initialize the table to the correct size
long double Imaginary(long double***[], long double, long double, long double, int);	//The tabulated points and 2 inputs and returns the bicubic interpolation of that input

int main(int argc, char* argv[])
{
	char* File = new char[25];	//Name of the file
	strcpy(File, argv[3]);
	Process = argv[1];
	strcat(File, ".");
	strcat(File, argv[5]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot(File);
	TPlot << setprecision(18);
	long double*** Table[2];	//The table of values computed by Spectral
	long double z, tau;	//The position value of the spactial correlator and tau of the euclidean-time correlator
	long double holder[5];
	int N[2] = {789,752}, M[2] = {151,463};	//The size of the table
	const int iProcess = atoi(argv[1]);
	const int Total = atoi(argv[2]);

	Init(Table, N, M);
	if(!ReadIn(Table, N, M, argv[4]);)
		return(0);
}

long double Imaginary(long double*** Table[], long double E, long double p, long double z, int Specify)
{
	long double Interpolation;
	long double t, u;
	int i, j;

	if(Specify == 1)
	{
		t = p/.8;	//returns the p index with the fractional part
		if(E < 2.5)	//These are to give the fractional distance from one E-point to the next+the index
			u = E/.1;
		else if(E < 2.540308)
			u = 25+(E-2.5)/.040308;
		else if(E < 3.540308)
			u = 26+(E-2.540308)/.005;
		else if(E < 3.55)
			u = 226+(E-3.540308)/.009692;
		else if(E < 4.1)
			u = 227+(E-3.55)/.01375;
		else
			u = 267+(E-4.1)/.1;
		i = t;	//p index in the Table
		j = u;	//E index in the Table
		t -= i;	//Removes the index leaving the fractional part
		u -= j;

		if(z > 0 && p > ((long double)(int(600*z/(2.*M_PI)))-.25)*2.*M_PI/z)
	                return(Spectral(Table, E, ((long double)(int(600*z/(2.*M_PI)))-.25)*2.*M_PI/z, z, Specify));
	}
	else
	{
		if(z > 0 && p > 788.-fmod(acos(cos((600.8+p)*z))+M_PI/2.,M_PI/2.)*2.*M_PI/z)
	                return(Spectral(Table, E, 788.-fmod(acos(cos((600.8+p)*z))+M_PI/2.,M_PI/2.)*2.*M_PI/z, z, Specify));
		i = p;
		j = E;
		t = p - i;
		u = E - j;
	}

	if(Table[Specify][i][j][4] == 0 || Table[Specify][i+1][j][4] == 0 || Table[Specify][i][j+1][4] == 0 || Table[Specify][i+1][j+1][4] == 0)	//If any of the required points have been invalidated, calculate points from integrals.
	{
		Interpolation = Table[Specify][i][j][0]*(1.-u)*(1.-t)+Table[Specify][i][j+1][0]*u*(1.-t)+Table[Specify][i+1][j][0]*(1.-u)*t+Table[Specify][i+1][j+1][0]*u*t;	//A bilinear interpolation
	}
	else
	{
		long double f[16] = {Table[Specify][i][j][0], Table[Specify][i+1][j][0], Table[Specify][i][j+1][0], Table[Specify][i+1][j+1][0], Table[Specify][i][j][1], Table[Specify][i+1][j][1], Table[Specify][i][j+1][1], Table[Specify][i+1][j+1][1], Table[Specify][i][j][2], Table[Specify][i+1][j][2], Table[Specify][i][j+1][2], Table[Specify][i+1][j+1][2], Table[Specify][i][j][3], Table[Specify][i+1][j][3], Table[Specify][i][j+1][3], Table[Specify][i+1][j+1][3]};	//fetch the data points and store them
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
	}

	return(Interpolation);
}

bool ReadIn(long double*** Table[], int N[], int M[], char* FileReadIn)
{
	ifstream File(FileReadIn);
	int i,j,k,m;	//Counters
	long double Holder;

	File.open(strcat(FileReadIn, ".xml"));
	if(File.fail())
	{
		cout << "Make sure the *.xml file has the same name as the Spectral function file dingus." << endl;
		return(false);
	}

	for(m = 0; m < 2; m++)	//Table count
	{
		for(i = 0; i < 4; i++)	//Value type count
		{
			for(j = 0; j < N[i]; j++)	//j/sqrt(s) count
			{
				for(k = 0; k < M[i]; k++)	//i/P count
				{
					File >> Table[m][j][k][i];
					if(m == 1)
					{
						if(i == 1 || i == 3)
							Table[m][j][k][i] *= .8;	//*=dP/di
						if(i == 2 || i == 3)
						{
							if(k < 25)
								Table[m][j][k][i] *= .1;
							else if(k < 26)
								Table[m][j][k][i] *= .040308;
							else if(k < 226)
								Table[m][j][k][i] *= .005;
							else if(k < 227)
								Table[m][j][k][i] *= .009692;
							else if(k < 267)
								Table[m][j][k][i] *= .01375;
							else
								Table[m][j][k][i] *= .1;
						}
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

	for(k = 0; k < 2; k++)
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
