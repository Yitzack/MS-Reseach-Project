//This program is written to integrate a 2-D function which is a form factor and propagator. The function is a function of theta and k and will not hold the correct properties if the average angle approximation is applied.
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstring>
#include<omp.h>
#include"Spectral.h"
using namespace std;

void Plot(long double, int);	//Plots the functions from table
long double RandFloat(long double, long double);
bool Poll(long double[2][84]);

char* Process;

int main(int argc, char* argv[])
{
/*#ifndef BB	//use option -D BB= to activate BB macro
	char File[30] = "Spectralcc.";  //Name of the file
#else
     	char File[30] = "Spectralbb.";  //Name of the file
#endif
	Process = argv[1];
	strcat(File, argv[3]);
	strcat(File, ".");
	strcat(File, Process);			//Appends the process number to the file name
	ofstream TPlot("Table");
	if(atoi(argv[6]) == 0)	//If starting from the beginning, overwrite
		TPlot.open(File);
	else	//If not starting from the beginning, append
		TPlot.open(File, ios::app);
	const int Temp = atoi(argv[3]);
	int i,j;	//counters
	const int iProcess = atoi(argv[1]) % atoi(argv[2]);
	const int Total = atoi(argv[2]);
	long double Table[616][3];//*/
	long double Par[5] = {-158.90117114622294, 2.643945190802571, 1.8, 0, 0};
	Elements holder[28];
	long double GaussLa[] = {0.0292089494940390418, 0.1539325380822080769, 0.3784519114339929046, 0.703043968841429832, 1.12804449030959115901, 1.65388906539884363591, 2.28111923347644653209, 3.01038628120128830529, 3.84245522739668292116, 4.77820943138205453677, 5.81865597642423461728, 6.96493193346708690195, 8.2183116110416122313, 9.58021491185883249065, 11.0522169380215279328, 12.63605901385725832108, 14.33366132857440339499, 16.14713744153402449126, 18.07881094274913343943, 20.13123462273780157763, 22.3072125823387678126, 24.60982580889231094881, 27.04246186610561423232, 29.60884949880154539486, 32.31309915127963456172, 35.15975065392247902555, 38.15382966748456817771, 41.3009149171740471975, 44.60721884062876818128, 48.0796850753673570501, 51.72610731101421216486, 55.55527556274067844963, 59.5771580886221159235, 63.80313029304261238365, 68.24626653908353044698, 72.92171766800947991981, 77.84720759844820215182, 83.04369909859864667464, 88.53630611197943572002, 94.35557619641319288989, 100.53934816696116679177, 107.13554136224855814149, 114.20653122712858723725, 121.83639878660318539969, 130.14381522449526055617, 139.30719756334274304328, 149.62081975792771442406, 161.64877015704720903095, 176.84630940701588372409};	//Displacement from 0 for Gauss-Laguerre integration

	/*TPlot << setprecision(18);	//18 digits is the "Number of decimal digits that can be rounded into a floating-point and back without change in the number of decimal digits" for long double.
	for(i = atoi(argv[6]); i <= 788; i++)	//Argv[6] allows to restart where ever
	{
		#pragma omp parallel for
		for(j = iProcess+151; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			Par[1] = 2.643945190802571;
			Par[2] = 1.8;
			if(argc == 4)
				switch(Temp)
				{
					case 1:
						Par[1] *= exp(-4./35.);
						Par[2] = 1.848;//1.902;
						break;
					case 2:
						Par[1] *= exp(-2./7.);
						Par[2] = 1.719;//1.777;
						break;
					case 3:
						Par[1] *= exp(-4./7.);
						Par[2] = 1.563;//1.652;
						break;
				}
			else
			{
				Par[1] *= atof(argv[4]);
				Par[2] = atof(argv[5]);
			}

			if(j <= 150)
			{
				cout << "If you come down this way, you have to alter the curvilinar system used to sample. You will want to use lines of constant s, or sqrt(E^2+P^2). You must do this before you do anything in s<0GeV^2 region" << endl;
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-j*j;
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-j*j;
				}
			}
			else
			{
				Par[3] = i*.8;
#ifndef BB
				if(j <= 181)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 381)
					Par[4] = pow((j-181.)/100.+3.,2);
				else if(j <= 566)
					Par[4] = pow((j-381.)/10.+5.,2);
#else
				if(j <= 251)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 451)
					Par[4] = pow((j-251.)/100.+10.,2);
				else if(j <= 566)
					Par[4] = pow((j-451.)/10.+12.,2);
#endif
				else
					Par[4] = 552.25+GaussLa[j-567];
			}

			if(j > 150 && i > 751)
			{
				Table[j][0] = Table[j][1] = Table[j][2] = 0;
			}
			else
			{
				holder = theta_Int(Par, Temp);
				Table[j][0] = holder.store(0);
				Table[j][1] = holder.store(1);
				Table[j][2] = holder.store(2);
			}
		}

		for(j = iProcess+151; j < 616; j+=Total)	//Does the subset of E that has been assigned to this process
		{
			if(j <= 150)
			{
				if(i <= 208)
				{
					Par[3] = i/10.+j;
					Par[4] = -i*j/5.-j*j;
				}
				else
				{
					Par[3] = i+j-187.2;
					Par[4] = -2.*j*(i-187.2)-j*j;
				}
			}
			else
			{
				Par[3] = i*.8;
#ifndef BB
				if(j <= 181)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 381)
					Par[4] = pow((j-181.)/100.+3.,2);
				else if(j <= 566)
					Par[4] = pow((j-381.)/10.+5.,2);
#else
				if(j <= 251)
					Par[4] = pow((j-151.)/10.,2);
				else if(j <= 451)
					Par[4] = pow((j-251.)/100.+10.,2);
				else if(j <= 566)
					Par[4] = pow((j-451.)/10.+12.,2);
#endif
				else
					Par[4] = 552.25+GaussLa[j-567];
			}

			TPlot << i << " " << j << " " << Par[3] << " " << Par[4] << " " << Table[j][0] << " " << Table[j][1] << " " << Table[j][2] << endl;
		}
		TPlot << endl;
	}

	//TPlot << "#Potiential Cutoff = " << Par[1] << " Mass = " << Par[2] << endl;
	TPlot.close();//*/

	cerr << setprecision(18);
	long double error[2][84];
	//long double Previous[] = {1.39461036156073218, 6.332062117923765, 7.01846419642699572, 17.3474665827736626, 20.0872380683069175, 29.8258531041717639, 46.704784395474316, 50.713420893829978, 51.278938722070447, 0.043158410019226763, 2.57317145318693376, 12.3076691166712644, 14.4865250760850598, 23.5028086983717032, 23.5103874011593978, 24.8348072100687947, 29.1168316508529893, 0.00594803112485079734, 0.00739439903373680452, 0.380668112496927022, 0.458735834559289807, 0.826951279575519325, 1.56823398390205848};	//well with P=20 and 600 at sqrt(s)=2
	//long double Previous[] = {0.123580343549897199, 0.541741934373112939, 7.69782055247488887, 8.72622414180869948, 11.6402506402220173, 26.0054177880140755, 33.3572779676262976, 40.5098268871267691, 42.4422552553474719, 1.44398460776065712, 1.91040368553818341, 2.5609217895117776, 3.00573801292983285, 9.12381233459490903, 10.1514320355772167, 12.4660249393208086, 15.3692872963081938, 0.0600075833235828519, 0.0786135669260155975, 0.669714802844075505, 0.794854226355800526, 1.2758869640745808, 1.56895224717512517};	//exceptional with P=200 at sqrt(s)=2
	long double Previous[] = {0.0428680552208614973, 0.351360127706610333, 0.400000000000000022, 0.41560939420741491, 2.4889659741563509, 27.5736341741087481, 30.3611216126246189, 64.1888719042480827, 65.7977831341683622, 0.0173088638452181461, 0.836829860598354491, 1.19162592788339509, 1.91451047357362453, 3.14427392232442597, 4.44058339437223329, 11.2562239844174467, 39.9320959431304673, 0.00853563204339899112, 0.0103170045937223034, 0.010382813423813732, 0.540182131169703438, 0.739622753171926845, 1.5};	//Pretty shitty at 23.5GeV
	long double slist[] = {.01, 3.24, 4., 12.96, 25.};//, 100., 552.25};
	long double Plist[] = {0, 20, 200, 600};
	int count;
	int i, j = 0;

	for(i = 0; i < 23; i++)
		Boundary[i] = Previous[i];

	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 5; j++)
		{
			Par[3] = Plist[i];
			Par[4] = slist[j];
			holder[i+4*j] = theta_Int(Par, 0);
			//cout << Par[3] << " " << sqrt(Par[4]) << " " << holder[i+4*j].store(0) << " " << holder[i+4*j].store(1) << " " << holder[i+4*j].store(2) << endl;
		}
	}

	cout << setprecision(18);
	for(i = 0; i < 20; i++)
	{
		error[0][3*i] = abs(holder[i].store(0)/holder[int(floor(i/4.))*4].store(0)-1.);
		error[0][3*i+1] = abs(holder[i].store(1)/holder[int(floor(i/4.))*4].store(1)-1.);
		error[0][3*i+2] = abs(holder[i].store(2)/holder[int(floor(i/4.))*4].store(2)-1.);
		cout << error[0][3*i] << " " << error[0][3*i+1] << " " << error[0][3*i+2] << " " << flush;
	}
	cout << setprecision(18);
	for(i = 0; i < 23; i++)
		cout << Previous[i] << " " << flush;
	cout << endl;
	srand(time(NULL)+atoi(argv[1])*30+100*atoi(argv[2]));

	do
	{
		count = rand()%23;
		switch(count)
		{
			case 0:
				Boundary[count] = RandFloat(0,Boundary[1]);
				break;
			case 8:
				Boundary[count] = RandFloat(Boundary[7],Boundary[8]+1.);
				break;
			case 9:
				Boundary[count] = RandFloat(0,Boundary[10]);
				break;
			case 16:
				Boundary[count] = RandFloat(Boundary[15],Boundary[16]+1.);
				break;
			case 17:
				Boundary[count] = RandFloat(0,Boundary[18]);
				break;
			case 22:
				Boundary[count] = RandFloat(Boundary[21],M_PI/2.);
				break;
			default:
				Boundary[count] = RandFloat(Boundary[count-1],Boundary[count+1]);
				break;
		}

		for(i = 0; i < 4; i++)
		{
			for(j = 0; j < 5; j++)
			{
				Par[3] = Plist[i];
				Par[4] = slist[j];
				holder[i+4*j] = theta_Int(Par, 0);
			}
		}

		cout << setprecision(18);
		for(i = 0; i < 20; i++)
		{
			error[1][3*i] = abs(holder[i].store(0)/holder[int(floor(i/4.))*4].store(0)-1.);
			error[1][3*i+1] = abs(holder[i].store(1)/holder[int(floor(i/4.))*4].store(1)-1.);
			error[1][3*i+2] = abs(holder[i].store(2)/holder[int(floor(i/4.))*4].store(2)-1.);
			cout << error[1][3*i] << " " << error[1][3*i+1] << " " << error[1][3*i+2] << " " << flush;
		}
		cout << setprecision(18);
		for(i = 0; i < 23; i++)
			cout << Boundary[i] << " " << flush;
		cout << endl;

		if(Poll(error))
		{	//reject
			cout << "attempt rejected" << endl;
			Boundary[count] = Previous[count];
		}
		else	//keep
		{
			cout << "attempt accepted" << endl;
			Previous[count] = Boundary[count];
			for(i = 0; i < 28; i++)
				error[0][i] = error[1][i];
		}

		j++;
	}while(j < 2000);//*/

	return(0);
}

bool Poll(long double error[2][84])
{
	long double count = 0;

	for(int i = 0; i < 60; i++)
	{
		switch(int(floor(i/12.)))	//{.1, 1.8, 2, 3.6, 5, 10, 23.5}
		{
			case 0:	//sqrt(s)==.1 GeV has 4 votes
				if(i%12 <= 2);	//Nan, no useful information
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 4.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 2.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 1:	//sqrt(s)==1.8 GeV has 4 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 4.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 2.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 2:	//sqrt(s)==2 GeV has 4 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					return(true);	//this one has veto rights
				else
					count -= 4.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 3:	//sqrt(s)==3.6 GeV has 4 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 4.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 2.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 4:	//sqrt(s)==5 GeV has 2 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 2.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 2.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 5:	//sqrt(s)==10 GeV has 2 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 2.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 1.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
			case 6:	//sqrt(s)==23.5 GeV has 2 votes
				if(i%12 <= 2);
				else if(abs(error[0][i]/error[1][i]-1.) < 1e-2);	//Don't care (vote present)
				else if(error[0][i] < error[1][i])	//Count reject conditions
					count += 2.*abs(error[0][i]/error[1][i]-1.);
				else
					count -= 1.*abs(error[0][i]/error[1][i]-1.);	//Count accept conditions
				break;
		}
	}

	long double rand = RandFloat(0,3);
	cout << rand << " " << count << endl;
	if(rand < count)
		return(true);
	return(false);
}

long double RandFloat(long double a, long double b)
{
	return(a+(b-a)*((long double)rand())/((long double)RAND_MAX));
}
