//This program is written to integrate a 3-D function which is a form factor and propagator. The function is a function of theta, k, and k_0 and will not hold the correct properties if the average angle approximation is applied.
#include<cmath>
#include<cstdlib>
#include<cfloat>
#include<complex>
#include<queue>
#include"Elements.h"
#include"Region.h"
using namespace std;

//Integrals that define results
Elements<Around> theta_Int(long double[], int);		//Theta integral
long double* k_Int(long double[], int, long double, int &);	//k integral
long double* k0_Int(long double[], int, long double, long double, int &);	//k0 integral aka energy integral
void Eval_Integral(long double[], Region&, int);

//Functions for finding points of interest in the k integral
void Characterize_k_Int(long double[], int, long double, long double[], long double[], int&);	//Returns the poles of the k integral's integrands
bool Newton_Method_k(long double&, long double, long double, long double, long double, long double, long double(*)(long double, long double), long double(*)(long double, long double, long double, long double, long double));			//Returns the k-intesection of two features of interest: potential peak, on-shell peak, light-like boundaries
long double V_Plus(long double, long double);			//Potiential peaks, matches a potential not in use
long double V_Minus(long double, long double);
long double Emm(long double, long double, long double, long double, long double);		//on-shell peaks
long double Epm(long double, long double, long double, long double, long double);
long double mEmp(long double, long double, long double, long double, long double);
long double Emp(long double, long double, long double, long double, long double);
long double mEpp(long double, long double, long double, long double, long double);
long double Epp(long double, long double, long double, long double, long double);
long double Upper_Bound(long double, long double, long double, long double, long double);	//light-like boundaries, where the imaginary self-energy has a known cusp
long double Lower_Bound(long double, long double, long double, long double, long double);

//Functions for finding points of interest in the k0 integral
void Characterize_k0_Int(long double[], int, long double, long double, long double[], long double[], int&);	//Returns the poles of the k0 integral's integrands
long double Newton_Method_k0(long double, long double[], long double, long double, int, int, long double (*)(long double[], long double, long double, long double, int, int));	//Returns the k0 of the on-shell peak using Newton's method on 1/f
long double omega_Width(long double, long double[], long double, long double, int, int, long double (*)(long double[], long double, long double, long double, int, int));	//Returns the width of on-shell peak using the assumption of a breit-wigner peak 

//Functions that return physics for the integrand
void ImSelf_Energy(long double, long double, long double[], long double[],int, long double[]);			//Returns the imaginary single quark self-energies for both quarks, contains an alternate T=194 MeV solution
long double ImSelf_Energy(long double, long double, long double, long double[], int);				//Returns the imaginary single quark self-energies for one quark, contains an alternate T=194 MeV solution
void ReSelf_Energy(long double, long double[], long double[], int, long double[]);					//Returns the real single quark self-energies for both quarks, contains an alternate T=194 MeV solution
void Self_Energy(long double, long double, long double[], long double[],int, long double[], long double[]);	//Returns the complex single quark self-energies for both quarks, is a simple Breit-Wigner self-energy and alternate to those above
long double Energy(long double, long double, long double, long double);						//Single quark energy, also used to return total momentum by setting M=0
long double Fermi(long double, int);											//Fermi function
long double Set_Temp(int);												//Decodes 0-4 into numeric temprature for Fermi factor
long double Potential1(long double[], long double, long double);							//One vertex of the potiential
long double Potential2(long double[], long double, long double);							//Two vertices of the potiential
long double Non_Interacting_Trace(long double[], long double, long double, long double);				//Non-interacting trace, depends on quantum numbers of boson (scalar, pseudo-scalar, vector, axial vector)
long double Interacting_Linear_Trace(long double[], long double, long double, long double);				//Linear (linear in sqrt(s)) contribution to the interacting trace
long double Interacting_Quad_Trace(long double[], long double, long double, long double);				//Quadratic contribution to the interacting trace
long double Imk0_Integrand(long double[], long double, long double, long double, int, int);				//Integrand of the k0 integral for positive energy

//Merge Sort. There are a number of semi-sorted lists that need sorting. This will probably beat quick sort under similar conditions.
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

#ifndef GAMMA	//use option -D GAMMA=<number> to alter single particle vacuum width, default value is 15MeV
#define GAMMA -.015	//Width of single quark propagator
#endif

Elements<Around> Int(long double Par[], int Temp)
{
	long double x1;					//Abscissa
	long double Boundary_theta[] = {1./17., 0.3, 0.08};	//Extra boundary values
	long double a, b, c, d, e, f;				//Boundaries
	Elements<Around> Answer(0,0,0,0,0,0);			//Answer to be returned
	Elements<long double> Total(0,0,0,0,0,0);			//Running Total
	Elements<long double> Error(0,0,0,0,0,0);			//Running Error
	int i, j, l;						//Counters

	if(Par[4] > 0 && Par[3] > sqrt(Par[4]/2.)) //Where the maximum of the theta integral ought to land. It might only be correct for BbS reduction, but is close enough for all other cases. Only valid for s>0 and P>sqrt(s/2)
		x1 = asin(sqrt(Par[4]/2.)/Par[3]);
	else
		x1 = M_PI/10.;	//If it isn't valid, value is needed anyways to split up the integral

	//Don't get too close to the pole or details might get lost
	if(x1>M_PI/10.)
		x1 = M_PI/10.;

	//List of boundaries between subintervals
	long double Range_theta[] = {x1*Boundary_theta[0], x1*Boundary_theta[1], x1, x1*(2.-Boundary_theta[1]), x1*(2.-Boundary_theta[1])*(1.-Boundary_theta[2])+M_PI/2.*Boundary_theta[2], M_PI/2., asin(sqrt(-Par[4])/Par[3]),0,0};
	long double* Range_k;
	long double* Range_k0;
	int Range_k_Elements;
	int Range_k0_Elements;
	priority_queue<long double, vector<long double>, greater<long double>> k_Stops;
	priority_queue<long double, vector<long double>, greater<long double>> k0_Stops;
	queue<Region> Initial_Regions;
	priority_queue<Region> Evaluated_Regions;
	queue<Region> Indivisible_Regions;
	Region Consideration[9];

	//Some kind of intersection, probably between the simultanous on-shell and potential peak, don't rightly remember
	Range_theta[7] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range_theta[7] = acos((pow(Range_theta[7],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(Range_theta[7]*Par[3]));
	Range_theta[8] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);
	Range_theta[8] = acos((pow(Range_theta[8],2)+pow(Par[2],2)-Par[4]-(long double).75*pow(Par[3],2))/(-Range_theta[8]*Par[3]));

	//Bad data trap for NaN and negative boundaries. These are only ones that can NaN or return negative numbers
	if(isnan(Range_theta[6]) || Range_theta[6] < 0) Range_theta[6] = M_PI;
	if(isnan(Range_theta[7])) Range_theta[7] = M_PI;
	if(isnan(Range_theta[8])) Range_theta[8] = M_PI;

	//Put in asending order
	mergeSort(Range_theta, 0, 8);

	for(i = 0; i < 9 && Range_theta[i] <= M_PI/2.; i++)	//Count through pre-determined intervals
	{
		Range_k = k_Int(Par, Temp, Range_theta[i], Range_k_Elements);	//Get the list of points of interest for each theta
		for(j = 0; j < Range_k_Elements; j++)					//And move them into k_Stops priority queue
		{
			while(Range_k[j] < 0 || isnan(Range_k[j])) j++;
			k_Stops.push(Range_k[j]);
		}
		delete Range_k;
	}

	Range_k_Elements = k_Stops.size();		//Copy k_Stops into Range_k so that it can be cycled over several times
	Range_k = new long double[Range_k_Elements];
	for(i = 0; i < Range_k_Elements; i++)
	{
		Range_k[i] = k_Stops.top();
		k_Stops.pop();
	}

	for(i = 0; i < 9 && Range_theta[i] <= M_PI/2.; i++)	//Count through pre-determined intervals
	{
		for(j = 0; j < Range_k_Elements; j++)
		{
			Range_k0 = k0_Int(Par, Temp, Range_theta[i], Range_k[j], Range_k0_Elements);	//Get the list of points of interest for each theta
			for(l = 0; l < Range_k0_Elements; l++)					//And move them into k0_Stops priority queue
			{
				while(Range_k0[j] < -sqrt(Par[4]+pow(Par[3],2)) || Range_k0[j] > sqrt(Par[4]+pow(Par[3],2)) || isnan(Range_k0[j])) j++;
				k0_Stops.push(Range_k0[j]);
			}
			delete Range_k0;
		}
	}

cout << i << " " << Range_k_Elements << " " << k0_Stops.size() << endl;
	e = 0;
	do
	{
		c = 0;
		while((abs(k0_Stops.top()/c-1.) < FLT_EPSILON || k0_Stops.top()==e) && !k0_Stops.empty())	//Work through the list k0_Stops
			k0_Stops.pop();
		if(k0_Stops.empty())	//Being sure not to pop and empty queue (that eats RAM and is unhelpful)
			break;
		f = k0_Stops.top();
		k0_Stops.pop();
		for(j = 0; j < Range_k_Elements; j++)
		{
			a = 0;
			d = Range_k[j];
			for(i = 1; i < 9 && Range_theta[i] <= M_PI/2.; i++)	//Work through the list of theta
			{
				b = Range_theta[i];
				if((abs(.5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(a),2)))-c) < 1.5 || 
					abs(.5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(b),2)))-c) < 1.5 || 
					abs(.5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(a),2)))-d) < 1.5 || 
					abs(.5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(b),2)))-d) < 1.5) &&
					(abs(.5*Par[3]*cos(a)*sqrt((Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(a),2)))-e) < 1.5 || 
					abs(.5*Par[3]*cos(b)*sqrt((Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(b),2)))-e) < 1.5 || 
					abs(.5*Par[3]*cos(a)*sqrt((Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(a),2)))-f) < 1.5 || 
					abs(.5*Par[3]*cos(b)*sqrt((Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(b),2)))-f) < 1.5))	//If any corner of a region is within 1.5 GeV of the simultaous on-shell, it and all of its descendants shall use the 97th order integral
				{
					Initial_Regions.emplace(a, b, c, d, e, f, 97);	//Add the region to the list of initial regions
				}
				else{Initial_Regions.emplace(a, b, c, d, e, f);}	//Default to 37th order for all descendants
				a = b;
			}
			c = d;
		}
		e = f;
	}while(!k0_Stops.empty());

	while(!Initial_Regions.empty())	//Go through the list of regions and give them an intial evaluation before storing them to the evaluated region list
	{
		Consideration[0] = Initial_Regions.front();
		Initial_Regions.pop();
		Eval_Integral(Par, Consideration[0], Temp);
		Evaluated_Regions.push(Consideration[0]);
		Total += Consideration[0].Int;	//Establish a running total for the total and error
		Error += Consideration[0].Err;
	}

cout << Par[3] << " " << Par[4] << " " << Total[0] << " " << Total[1] << " " << Total[2] << " " << Total[3] << " " << Total[4] << " " << Total[5] << " " << (Error/Total)[0] << " " << (Error/Total)[1] << " " << (Error/Total)[2] << " " << (Error/Total)[3] << " " << (Error/Total)[4] << " " << (Error/Total)[5] << " " << Evaluated_Regions.top().Area() << " " << (abs(Evaluated_Regions.top().Err/Error))[0] << " " << (abs(Evaluated_Regions.top().Err/Error))[1] << " " << (abs(Evaluated_Regions.top().Err/Error))[2] << " " << (abs(Evaluated_Regions.top().Err/Error))[3] << " " << (abs(Evaluated_Regions.top().Err/Error))[4] << " " << (abs(Evaluated_Regions.top().Err/Error))[5] << " " << Evaluated_Regions.size()+Indivisible_Regions.size() << endl;
	while((abs(Error/Total) >= 1e-6) && !Evaluated_Regions.empty() && abs(Evaluated_Regions.top().Err/Error)/Evaluated_Regions.top().Area() >= 1e-6)
	{//While the evaluated regions aren't empty and accuracy goals aren't reached (absolute error on left and all regions relative error condition on right)
		Consideration[0] = Evaluated_Regions.top();	//Consider the top region
		Evaluated_Regions.pop();
		Total -= Consideration[0].Int;	//Running total update
		Error -= Consideration[0].Err;
		if(((abs(Consideration[0].xErr/Consideration[0].yErr-1.) < 1. && abs(Consideration[0].zErr/Consideration[0].yErr-1.) < 1.) || Consideration[0].Err == 0) && Consideration[0].yDeep < 4 && Consideration[0].xDeep < 4 && Consideration[0].zDeep < 4)	//Divide both dimensions, they're roughly equally bad. Consideration[0].Err == 0 subdivides unevaluated regions only if they aren't 4 levels deep
		{
			Consideration[1] = Region(Consideration[0].x1, (Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].y1, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].z1, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].order);
			Consideration[2] = Region((Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].x2, Consideration[0].y1, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].z1, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].order);
			Consideration[3] = Region(Consideration[0].x1, (Consideration[0].x1+Consideration[0].x2)/2., (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].y2, Consideration[0].z1, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].order);
			Consideration[4] = Region((Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].x2, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].y2, Consideration[0].z1, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].order);
			Consideration[5] = Region(Consideration[0].x1, (Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].y1, (Consideration[0].y1+Consideration[0].y2)/2., (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].z2, Consideration[0].order);
			Consideration[6] = Region((Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].x2, Consideration[0].y1, (Consideration[0].y1+Consideration[0].y2)/2., (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].z2, Consideration[0].order);
			Consideration[7] = Region(Consideration[0].x1, (Consideration[0].x1+Consideration[0].x2)/2., (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].y2, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].z2, Consideration[0].order);
			Consideration[8] = Region((Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].x2, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].y2, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].z2, Consideration[0].order);
			for(i = 1; i <= 8; i++)
			{
				Eval_Integral(Par, Consideration[i], Temp);			//Evaluate the subdivided regions
				Consideration[i].xDeep = Consideration[0].xDeep+1;	//Increament the depth
				Consideration[i].yDeep = Consideration[0].yDeep+1;
				Consideration[i].zDeep = Consideration[0].zDeep+1;
				Total += Consideration[i].Int;			//Update the running total
				Error += Consideration[i].Err;
				Evaluated_Regions.push(Consideration[i]);		//Push the subdivided regions to the queue
			}
		}
		else if(Consideration[0].xErr > Consideration[0].yErr && Consideration[0].xErr > Consideration[0].zErr && Consideration[0].xDeep < 4)	//x is worst, divide x
		{
			Consideration[1] = Region(Consideration[0].x1, (Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].y1, Consideration[0].y2, Consideration[0].z1, Consideration[0].z2, Consideration[0].order);
			Consideration[2] = Region((Consideration[0].x1+Consideration[0].x2)/2., Consideration[0].x2, Consideration[0].y1, Consideration[0].y2, Consideration[0].z1, Consideration[0].z2, Consideration[0].order);
			for(i = 1; i <= 2; i++)
			{
				Eval_Integral(Par, Consideration[i], Temp);
				Consideration[i].xDeep = Consideration[0].xDeep+1;
				Total += Consideration[i].Int;
				Error += Consideration[i].Err;
				Evaluated_Regions.push(Consideration[i]);
			}
		}
		else if(Consideration[0].yErr > Consideration[0].xErr && Consideration[0].yErr > Consideration[0].zErr && Consideration[0].yDeep < 4)	//y is worst, divide y
		{
			Consideration[1] = Region(Consideration[0].x1, Consideration[0].x2, Consideration[0].y1, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].z1, Consideration[0].z2, Consideration[0].order);
			Consideration[2] = Region(Consideration[0].x1, Consideration[0].x2, (Consideration[0].y1+Consideration[0].y2)/2., Consideration[0].y2, Consideration[0].z1, Consideration[0].z2, Consideration[0].order);
			for(i = 1; i <= 2; i++)
			{
				Eval_Integral(Par, Consideration[i], Temp);
				Consideration[i].yDeep = Consideration[0].yDeep+1;
				Total += Consideration[i].Int;
				Error += Consideration[i].Err;
				Evaluated_Regions.push(Consideration[i]);
			}
		}
		else if(Consideration[0].zErr > Consideration[0].xErr && Consideration[0].zErr > Consideration[0].yErr && Consideration[0].zDeep < 4)	//z is worst, divide z
		{
			Consideration[1] = Region(Consideration[0].x1, Consideration[0].x2, Consideration[0].y1, Consideration[0].y2, Consideration[0].z1, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].order);
			Consideration[2] = Region(Consideration[0].x1, Consideration[0].x2, Consideration[0].y1, Consideration[0].y2, (Consideration[0].z1+Consideration[0].z2)/2., Consideration[0].z2, Consideration[0].order);
			for(i = 1; i <= 2; i++)
			{
				Eval_Integral(Par, Consideration[i], Temp);
				Consideration[i].yDeep = Consideration[0].zDeep+1;
				Total += Consideration[i].Int;
				Error += Consideration[i].Err;
				Evaluated_Regions.push(Consideration[i]);
			}
		}
		else	//if none of the above happen, the region is propbably indivisible, store it that queue so it doesn't come back to the top of the priority queue
		{
			Indivisible_Regions.push(Consideration[0]);
			Total += Consideration[0].Int;
			Error += Consideration[0].Err;
		}
cout << Par[3] << " " << Par[4] << " " << Total[0] << " " << Total[1] << " " << Total[2] << " " << Total[3] << " " << Total[4] << " " << Total[5] << " " << (Error/Total)[0] << " " << (Error/Total)[1] << " " << (Error/Total)[2] << " " << (Error/Total)[3] << " " << (Error/Total)[4] << " " << (Error/Total)[5] << " " << Evaluated_Regions.top().Area() << " " << (abs(Evaluated_Regions.top().Err/Error))[0] << " " << (abs(Evaluated_Regions.top().Err/Error))[1] << " " << (abs(Evaluated_Regions.top().Err/Error))[2] << " " << (abs(Evaluated_Regions.top().Err/Error))[3] << " " << (abs(Evaluated_Regions.top().Err/Error))[4] << " " << (abs(Evaluated_Regions.top().Err/Error))[5] << " " << Evaluated_Regions.size()+Indivisible_Regions.size() << endl;
	}

	Total = Elements<long double>(0,0,0,0,0,0);	//Total's error estimate has probably been messed up by adding and subtracting regions of the running total
	Error = Elements<long double>(0,0,0,0,0,0);
	while(!Indivisible_Regions.empty())	//Resum the running total so that the error estimates aren't inflated
	{
		Total += Indivisible_Regions.front().Int;
		Error += Indivisible_Regions.front().Err;
		Indivisible_Regions.pop();
	}
	while(!Evaluated_Regions.empty())
	{
		Total += Evaluated_Regions.top().Int;
		Error += Evaluated_Regions.top().Err;
		Evaluated_Regions.pop();
	}
	Answer = Elements<Around>(Around(Total[0], Error[0]), Around(Total[1], Error[1]), Around(Total[2], Error[2]), Around(Total[3], Error[3]), Around(Total[4], Error[4]), Around(Total[5], Error[5]));//Assemble the answer

	return(Answer/pow(2.*M_PI,2)*2.);
}

void Eval_Integral(long double Par[], Region& Stuff, int Temp)
{
//Displacements, weights, and error weights for 37th order Guass-Kronrod
	const long double Disp37[] = {0.1252334085114689154724414, 0.2485057483204692762677910, 0.3678314989981801937526915, 0.4813394504781570929359436, 0.5873179542866174472967024, 0.6840598954700558939449291, 0.7699026741943046870368938, 0.8435581241611532447921419, 0.9041172563704748566784659, 0.9505377959431212965490602, 0.9815606342467192506905491, 0.9969339225295954269123502};
	const long double w37[] = {0.12555689390547433530429613, 0.1245841645361560734373125, 0.12162630352394838324609976, 0.1167120535017568262935807, 0.11002260497764407263590740, 0.10164973227906027771568877, 0.091549468295049210528171940, 0.07992027533360170149339261, 0.067250907050839930304940940, 0.05369701760775625122888916, 0.038915230469299477115089632, 0.02303608403898223259108458, 0.0082577114331683957576939224};
	const long double errw37[] = {0.12555689390547433530429613, -0.1245628812772467115632500, 0.12162630352394838324609976, -0.11678048303659798246726915, 0.11002260497764407263590740, -0.10151769444400564403337569, 0.091549468295049210528171940, -0.08015805320974452484125992, 0.067250907050839930304940940, -0.05324230838756217973136555, 0.038915230469299477115089632, -0.024139252347529594603531381, 0.0082577114331683957576939224};
//Displacements, weights, and error weights for 97th order Guass-Kronrod
	const long double Disp97[] = {0.0483076656877383162348126, 0.0965026968768943658008313, 0.1444719615827964934851864, 0.1921036089831424972716416, 0.2392873622521370745446032, 0.2859124585894597594166071, 0.3318686022821276497799168, 0.3770494211541211054453355, 0.4213512761306353453641194, 0.4646693084819922177561782, 0.5068999089322293900237475, 0.5479463141991524786809395, 0.5877157572407623290407455, 0.6261129377018239978202384, 0.6630442669302152009751152, 0.6984265577952104928847701, 0.7321821187402896803874267, 0.7642282519978037041506601, 0.7944837959679424069630973, 0.8228829501360513216482688, 0.8493676137325699701336930, 0.8738697689453106061296618, 0.8963211557660521239653072, 0.9166772666513643242753457, 0.9349060759377396891709191, 0.9509546848486611853898828, 0.9647622555875064307738119, 0.9763102836146638071976696, 0.9856115115452683354001750, 0.9926280352629719126857912, 0.9972638618494815635449811, 0.9995459021243644786356103};
	const long double w97[] = {0.048326383986567758375445434, 0.0482701930757773855987121, 0.048100969185457746927846544, 0.04781890873698847221226358, 0.047426061873882382362879950, 0.04692296828170361110348071, 0.046308756738025713240381298, 0.04558582656454707028057546, 0.044758638749766937295199192, 0.04382754403013974904681615, 0.042791115596446746933654925, 0.04165401998564305139829641, 0.040423492370373096672349269, 0.03909942013330661120748213, 0.037679130645613398514895974, 0.03616976947564229986095839, 0.034582122744733034130726383, 0.03291507764390360026329648, 0.031163325561973737171155849, 0.02933695668962066136861561, 0.027452098422210403783147707, 0.02550569548089465281452890, 0.023486659672163324592087913, 0.02140891318482191595577752, 0.019298771430326811294403740, 0.01714980520978425325608583, 0.014936103606086027385096751, 0.01267605480665440285936888, 0.010423987398806818828034251, 0.008172504038531668414343805, 0.0058417370791666933039479766, 0.003426818775772370935574576, 0.0012233608179514718002930372};
	const long double errw97[] = {0.048326383986567758375445434, -0.0482698954389504149680528, 0.048100969185457746927846544, -0.04781981134228638720681842, 0.047426061873882382362879950, -0.04692143079910095453569952, 0.046308756738025713240381298, -0.04558805213121681443229312, 0.044758638749766937295199192, -0.04382454897426406209595531, 0.042791115596446746933654925, -0.04165790424130370382390266, 0.040423492370373096672349269, -0.03909447565376369526425879, 0.037679130645613398514895974, -0.03617602463320620636444097, 0.034582122744733034130726383, -0.03290714513245824657435359, 0.031163325561973737171155849, -0.02934713678891488577666802, 0.027452098422210403783147707, -0.02549236378148152338163434, 0.023486659672163324592087913, -0.021426984837404764701101129, 0.019298771430326811294403740, -0.017124057703237179846601906, 0.014936103606086027385096751, -0.012716010502607656596383706, 0.010423987398806818828034251, -0.008101890692374002190826757, 0.0058417370791666933039479766, -0.003591791233697725664832488, 0.0012233608179514718002930372};

	long double x1, x2, y1, y2, z1, z2;	//Abscissa
	long double k011, k012, k021, k022;	//On-shell energy transfer at abscissa

	long double a = Stuff.x1;	//Limits of integration
	long double b = Stuff.x2;
	long double c = Stuff.y1;
	long double d = Stuff.y2;
	long double e = Stuff.z1;
	long double f = Stuff.z2;

	Elements<long double>Holder[8];	//8 holder varibles

	Stuff.Int = Elements<long double>(0,0,0,0,0,0);	//Intialize integration varibles
	Stuff.xErr = Elements<long double>(0,0,0,0,0,0);
	Stuff.yErr = Elements<long double>(0,0,0,0,0,0);
	Stuff.zErr = Elements<long double>(0,0,0,0,0,0);
	Stuff.Err = Elements<long double>(0,0,0,0,0,0);

	if(Stuff.order == 37)
	{
		for(int i = 0; i < 12; i++)//for(int l = 0; l < 12; l+=2)// //Count through points away from center
		{
			x1 = (b+a-Disp37[i]*(b-a))/2.;	//theta
			x2 = (b+a+Disp37[i]*(b-a))/2.;
			for(int j = 0; j < 12; j++)
			{
				y1 = (d+c-Disp37[j]*(d-c))/2.;	//k
				y2 = (d+c+Disp37[j]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;
				k022 = (Energy(Par[2], Par[3]/2., y2, x2)-Energy(Par[2], Par[3]/2., -y2, x2))/2.;

				for(int k = 0; k < 12; k++)
				{
					z1 = (f+e-Disp37[k]*(f-e))/2.;	//k0
					z2 = (f+e+Disp37[k]*(f-e))/2.;

					Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
					Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
					Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
					Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z1, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);
					Holder[4] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
					Holder[5] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z2, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
					Holder[6] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z2, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
					Holder[7] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z2, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);

					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[0];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[1];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[2];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[3];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[4];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[5];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[6];
					Stuff.Int += w37[i+1]*w37[j+1]*w37[k+1]*Holder[7];

					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[0];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[1];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[2];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[3];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[4];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[5];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[6];
					Stuff.xErr += errw37[i+1]*w37[j+1]*w37[k+1]*Holder[7];

					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[0];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[1];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[2];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[3];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[4];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[5];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[6];
					Stuff.yErr += errw37[j+1]*w37[i+1]*w37[k+1]*Holder[7];

					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[0];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[1];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[2];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[3];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[4];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[5];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[6];
					Stuff.zErr += errw37[k+1]*w37[j+1]*w37[i+1]*Holder[7];

					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[0];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[1];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[2];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[3];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[4];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[5];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[6];
					Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[k+1]*Holder[7];
				}
			}
		}
		for(int i = 0; i < 12; i++)
		{
			for(int j = 0; j < 12; j++)
			{
				x1 = (b+a-Disp37[i]*(b-a))/2.;	//theta
				x2 = (b+a+Disp37[i]*(b-a))/2.;

				y1 = (d+c-Disp37[j]*(d-c))/2.;	//k
				y2 = (d+c+Disp37[j]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;
				k022 = (Energy(Par[2], Par[3]/2., y2, x2)-Energy(Par[2], Par[3]/2., -y2, x2))/2.;

				z1 = (f+e)/2.;	//k0

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z1, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);

				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[0];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[1];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[2];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[3];

				Stuff.xErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[0];
				Stuff.xErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[1];
				Stuff.xErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[2];
				Stuff.xErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[3];

				Stuff.yErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[0];
				Stuff.yErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[1];
				Stuff.yErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[2];
				Stuff.yErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[3];

				Stuff.zErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[0];
				Stuff.zErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[1];
				Stuff.zErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[2];
				Stuff.zErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[3];

				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[0];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[1];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[2];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[3];

				x1 = (b+a-Disp37[i]*(b-a))/2.;	//theta
				x2 = (b+a+Disp37[i]*(b-a))/2.;

				y1 = (d+c)/2.;	//k

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;

				z1 = (f+e-Disp37[j]*(f-e))/2.;	//k0
				z2 = (f+e+Disp37[j]*(f-e))/2.;

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z2, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);

				Stuff.Int += w37[i+1]*w37[0]*w37[j+1]*Holder[0];
				Stuff.Int += w37[i+1]*w37[0]*w37[j+1]*Holder[1];
				Stuff.Int += w37[i+1]*w37[0]*w37[j+1]*Holder[2];
				Stuff.Int += w37[i+1]*w37[0]*w37[j+1]*Holder[3];

				Stuff.xErr += errw37[i+1]*w37[0]*w37[j+1]*Holder[0];
				Stuff.xErr += errw37[i+1]*w37[0]*w37[j+1]*Holder[1];
				Stuff.xErr += errw37[i+1]*w37[0]*w37[j+1]*Holder[2];
				Stuff.xErr += errw37[i+1]*w37[0]*w37[j+1]*Holder[3];

				Stuff.yErr += errw37[0]*w37[i+1]*w37[j+1]*Holder[0];
				Stuff.yErr += errw37[0]*w37[i+1]*w37[j+1]*Holder[1];
				Stuff.yErr += errw37[0]*w37[i+1]*w37[j+1]*Holder[2];
				Stuff.yErr += errw37[0]*w37[i+1]*w37[j+1]*Holder[3];

				Stuff.zErr += errw37[j+1]*w37[0]*w37[i+1]*Holder[0];
				Stuff.zErr += errw37[j+1]*w37[0]*w37[i+1]*Holder[1];
				Stuff.zErr += errw37[j+1]*w37[0]*w37[i+1]*Holder[2];
				Stuff.zErr += errw37[j+1]*w37[0]*w37[i+1]*Holder[3];

				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[0];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[1];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[2];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[3];

				x1 = (b+a)/2.;	//theta

				y1 = (d+c-Disp37[i]*(d-c))/2.;	//k
				y2 = (d+c+Disp37[i]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;

				z1 = (f+e-Disp37[j]*(f-e))/2.;	//k0
				z2 = (f+e+Disp37[j]*(f-e))/2.;

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z2, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);

				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[0];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[1];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[2];
				Stuff.Int += w37[i+1]*w37[j+1]*w37[0]*Holder[3];

				Stuff.xErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[0];
				Stuff.xErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[1];
				Stuff.xErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[2];
				Stuff.xErr += errw37[0]*w37[j+1]*w37[i+1]*Holder[3];

				Stuff.yErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[0];
				Stuff.yErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[1];
				Stuff.yErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[2];
				Stuff.yErr += errw37[i+1]*w37[j+1]*w37[0]*Holder[3];

				Stuff.zErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[0];
				Stuff.zErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[1];
				Stuff.zErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[2];
				Stuff.zErr += errw37[j+1]*w37[i+1]*w37[0]*Holder[3];

				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[0];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[1];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[2];
				Stuff.Err += errw37[i+1]*errw37[j+1]*errw37[0]*Holder[3];
			}
		}
		for(int i = 0; i < 12; i++)
		{
			x1 = (b+a)/2.;	//theta
			y1 = (d+c)/2.;	//k
			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;

			z1 = (f+e-Disp37[i]*(f-e))/2.;	//k0
			z2 = (f+e+Disp37[i]*(f-e))/2.;

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);

			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[0];
			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[1];

			Stuff.yErr += errw37[0]*w37[i+1]*w37[0]*Holder[0];
			Stuff.yErr += errw37[0]*w37[i+1]*w37[0]*Holder[1];

			Stuff.zErr += errw37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.zErr += errw37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[0];
			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[1];

			x1 = (b+a)/2.;	//theta

			y1 = (d+c-Disp37[i]*(d-c))/2.;	//k
			y2 = (d+c+Disp37[i]*(d-c))/2.;

			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
			k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;

			z1 = (f+e)/2.;	//k0

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);

			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[0];
			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[1];

			Stuff.yErr += errw37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.yErr += errw37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.zErr += errw37[0]*w37[0]*w37[i+1]*Holder[0];
			Stuff.zErr += errw37[0]*w37[0]*w37[i+1]*Holder[1];

			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[0];
			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[1];

			x1 = (b+a-Disp37[i]*(b-a))/2.;	//theta
			x2 = (b+a+Disp37[i]*(b-a))/2.;

			y1 = (d+c)/2.;	//k

			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
			k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;

			z1 = (f+e)/2.;	//k0

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);

			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.Int += w37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[0];
			Stuff.xErr += errw37[0]*w37[i+1]*w37[0]*Holder[1];

			Stuff.yErr += errw37[0]*w37[i+1]*w37[0]*Holder[0];
			Stuff.yErr += errw37[0]*w37[i+1]*w37[0]*Holder[1];

			Stuff.zErr += errw37[i+1]*w37[0]*w37[0]*Holder[0];
			Stuff.zErr += errw37[i+1]*w37[0]*w37[0]*Holder[1];

			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[0];
			Stuff.Err += errw37[i+1]*errw37[0]*errw37[0]*Holder[1];
		}

		x1 = (b+a)/2.;	//theta
		y1 = (d+c)/2.;	//k
		k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
		z1 = (f+e)/2.;	//k0

		Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);

		Stuff.Int += w37[0]*w37[0]*w37[0]*Holder[0];
		Stuff.xErr += errw37[0]*w37[0]*w37[0]*Holder[0];
		Stuff.yErr += w37[0]*errw37[0]*w37[0]*Holder[0];
		Stuff.zErr += w37[0]*errw37[0]*w37[0]*Holder[0];
		Stuff.Err += errw37[0]*errw37[0]*errw37[0]*Holder[0];
	}
	else if(Stuff.order == 97)
	{
		for(int i = 0; i < 32; i++)//for(int l = 0; l < 32; l+=2)// //Count through points away from center
		{
			x1 = (b+a-Disp97[i]*(b-a))/2.;	//theta
			x2 = (b+a+Disp97[i]*(b-a))/2.;
			for(int j = 0; j < 32; j++)
			{
				y1 = (d+c-Disp97[j]*(d-c))/2.;	//k
				y2 = (d+c+Disp97[j]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;
				k022 = (Energy(Par[2], Par[3]/2., y2, x2)-Energy(Par[2], Par[3]/2., -y2, x2))/2.;

				for(int k = 0; k < 32; k++)
				{
					z1 = (f+e-Disp97[k]*(f-e))/2.;	//k0
					z2 = (f+e+Disp97[k]*(f-e))/2.;

					Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
					Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
					Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
					Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z1, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);
					Holder[4] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
					Holder[5] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z2, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
					Holder[6] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z2, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
					Holder[7] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z2, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);

					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[0];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[1];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[2];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[3];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[4];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[5];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[6];
					Stuff.Int += w97[i+1]*w97[j+1]*w97[k+1]*Holder[7];

					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[0];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[1];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[2];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[3];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[4];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[5];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[6];
					Stuff.xErr += errw97[i+1]*w97[j+1]*w97[k+1]*Holder[7];

					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[0];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[1];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[2];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[3];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[4];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[5];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[6];
					Stuff.yErr += errw97[j+1]*w97[i+1]*w97[k+1]*Holder[7];

					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[0];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[1];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[2];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[3];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[4];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[5];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[6];
					Stuff.zErr += errw97[k+1]*w97[j+1]*w97[i+1]*Holder[7];

					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[0];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[1];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[2];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[3];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[4];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[5];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[6];
					Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[k+1]*Holder[7];
				}
			}
		}
		for(int i = 0; i < 32; i++)
		{
			for(int j = 0; j < 32; j++)
			{
				x1 = (b+a-Disp97[i]*(b-a))/2.;	//theta
				x2 = (b+a+Disp97[i]*(b-a))/2.;

				y1 = (d+c-Disp97[j]*(d-c))/2.;	//k
				y2 = (d+c+Disp97[j]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;
				k022 = (Energy(Par[2], Par[3]/2., y2, x2)-Energy(Par[2], Par[3]/2., -y2, x2))/2.;

				z1 = (f+e)/2.;	//k0

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k022, y2, x2), Potential1(Par,k022, y2), Interacting_Linear_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Interacting_Quad_Trace(Par, k022, y2, x2)*Potential1(Par,k022, y2), Potential2(Par,k022, x2)))*Imk0_Integrand(Par, z1, y2, x2, Temp, 0)*pow(y2,2)*sin(x2);

				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[0];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[1];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[2];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[3];

				Stuff.xErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[0];
				Stuff.xErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[1];
				Stuff.xErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[2];
				Stuff.xErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[3];

				Stuff.yErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[0];
				Stuff.yErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[1];
				Stuff.yErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[2];
				Stuff.yErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[3];

				Stuff.zErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[0];
				Stuff.zErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[1];
				Stuff.zErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[2];
				Stuff.zErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[3];

				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[0];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[1];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[2];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[3];

				x1 = (b+a-Disp97[i]*(b-a))/2.;	//theta
				x2 = (b+a+Disp97[i]*(b-a))/2.;

				y1 = (d+c)/2.;	//k

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;

				z1 = (f+e-Disp97[j]*(f-e))/2.;	//k0
				z2 = (f+e+Disp97[j]*(f-e))/2.;

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z2, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);

				Stuff.Int += w97[i+1]*w97[0]*w97[j+1]*Holder[0];
				Stuff.Int += w97[i+1]*w97[0]*w97[j+1]*Holder[1];
				Stuff.Int += w97[i+1]*w97[0]*w97[j+1]*Holder[2];
				Stuff.Int += w97[i+1]*w97[0]*w97[j+1]*Holder[3];

				Stuff.xErr += errw97[i+1]*w97[0]*w97[j+1]*Holder[0];
				Stuff.xErr += errw97[i+1]*w97[0]*w97[j+1]*Holder[1];
				Stuff.xErr += errw97[i+1]*w97[0]*w97[j+1]*Holder[2];
				Stuff.xErr += errw97[i+1]*w97[0]*w97[j+1]*Holder[3];

				Stuff.yErr += errw97[0]*w97[i+1]*w97[j+1]*Holder[0];
				Stuff.yErr += errw97[0]*w97[i+1]*w97[j+1]*Holder[1];
				Stuff.yErr += errw97[0]*w97[i+1]*w97[j+1]*Holder[2];
				Stuff.yErr += errw97[0]*w97[i+1]*w97[j+1]*Holder[3];

				Stuff.zErr += errw97[j+1]*w97[0]*w97[i+1]*Holder[0];
				Stuff.zErr += errw97[j+1]*w97[0]*w97[i+1]*Holder[1];
				Stuff.zErr += errw97[j+1]*w97[0]*w97[i+1]*Holder[2];
				Stuff.zErr += errw97[j+1]*w97[0]*w97[i+1]*Holder[3];

				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[0];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[1];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[2];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[3];

				x1 = (b+a)/2.;	//theta

				y1 = (d+c-Disp97[i]*(d-c))/2.;	//k
				y2 = (d+c+Disp97[i]*(d-c))/2.;

				k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
				k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;

				z1 = (f+e-Disp97[j]*(f-e))/2.;	//k0
				z2 = (f+e+Disp97[j]*(f-e))/2.;

				Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);
				Holder[2] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
				Holder[3] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z2, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);

				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[0];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[1];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[2];
				Stuff.Int += w97[i+1]*w97[j+1]*w97[0]*Holder[3];

				Stuff.xErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[0];
				Stuff.xErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[1];
				Stuff.xErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[2];
				Stuff.xErr += errw97[0]*w97[j+1]*w97[i+1]*Holder[3];

				Stuff.yErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[0];
				Stuff.yErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[1];
				Stuff.yErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[2];
				Stuff.yErr += errw97[i+1]*w97[j+1]*w97[0]*Holder[3];

				Stuff.zErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[0];
				Stuff.zErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[1];
				Stuff.zErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[2];
				Stuff.zErr += errw97[j+1]*w97[i+1]*w97[0]*Holder[3];

				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[0];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[1];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[2];
				Stuff.Err += errw97[i+1]*errw97[j+1]*errw97[0]*Holder[3];
			}
		}
		for(int i = 0; i < 32; i++)
		{
			x1 = (b+a)/2.;	//theta
			y1 = (d+c)/2.;	//k
			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;

			z1 = (f+e-Disp97[i]*(f-e))/2.;	//k0
			z2 = (f+e+Disp97[i]*(f-e))/2.;

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z2, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);

			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[0];
			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[1];

			Stuff.yErr += errw97[0]*w97[i+1]*w97[0]*Holder[0];
			Stuff.yErr += errw97[0]*w97[i+1]*w97[0]*Holder[1];

			Stuff.zErr += errw97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.zErr += errw97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[0];
			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[1];

			x1 = (b+a)/2.;	//theta

			y1 = (d+c-Disp97[i]*(d-c))/2.;	//k
			y2 = (d+c+Disp97[i]*(d-c))/2.;

			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
			k012 = (Energy(Par[2], Par[3]/2., y2, x1)-Energy(Par[2], Par[3]/2., -y2, x1))/2.;

			z1 = (f+e)/2.;	//k0

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k012, y2, x1), Potential1(Par,k012, y2), Interacting_Linear_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Interacting_Quad_Trace(Par, k012, y2, x1)*Potential1(Par,k012, y2), Potential2(Par,k012, x1)))*Imk0_Integrand(Par, z1, y2, x1, Temp, 0)*pow(y2,2)*sin(x1);

			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[0];
			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[1];

			Stuff.yErr += errw97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.yErr += errw97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.zErr += errw97[0]*w97[0]*w97[i+1]*Holder[0];
			Stuff.zErr += errw97[0]*w97[0]*w97[i+1]*Holder[1];

			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[0];
			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[1];

			x1 = (b+a-Disp97[i]*(b-a))/2.;	//theta
			x2 = (b+a+Disp97[i]*(b-a))/2.;

			y1 = (d+c)/2.;	//k

			k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
			k021 = (Energy(Par[2], Par[3]/2., y1, x2)-Energy(Par[2], Par[3]/2., -y1, x2))/2.;

			z1 = (f+e)/2.;	//k0

			Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);
			Holder[1] = (Elements<long double>(2., Non_Interacting_Trace(Par, k021, y1, x2), Potential1(Par,k021, y1), Interacting_Linear_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Interacting_Quad_Trace(Par, k021, y1, x2)*Potential1(Par,k021, y1), Potential2(Par,k021, x2)))*Imk0_Integrand(Par, z1, y1, x2, Temp, 0)*pow(y1,2)*sin(x2);

			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.Int += w97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[0];
			Stuff.xErr += errw97[0]*w97[i+1]*w97[0]*Holder[1];

			Stuff.yErr += errw97[0]*w97[i+1]*w97[0]*Holder[0];
			Stuff.yErr += errw97[0]*w97[i+1]*w97[0]*Holder[1];

			Stuff.zErr += errw97[i+1]*w97[0]*w97[0]*Holder[0];
			Stuff.zErr += errw97[i+1]*w97[0]*w97[0]*Holder[1];

			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[0];
			Stuff.Err += errw97[i+1]*errw97[0]*errw97[0]*Holder[1];
		}

		x1 = (b+a)/2.;	//theta
		y1 = (d+c)/2.;	//k
		k011 = (Energy(Par[2], Par[3]/2., y1, x1)-Energy(Par[2], Par[3]/2., -y1, x1))/2.;
		z1 = (f+e)/2.;	//k0

		Holder[0] = (Elements<long double>(2., Non_Interacting_Trace(Par, k011, y1, x1), Potential1(Par,k011, y1), Interacting_Linear_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Interacting_Quad_Trace(Par, k011, y1, x1)*Potential1(Par,k011, y1), Potential2(Par,k011, x1)))*Imk0_Integrand(Par, z1, y1, x1, Temp, 0)*pow(y1,2)*sin(x1);

		Stuff.Int += w97[0]*w97[0]*w97[0]*Holder[0];
		Stuff.xErr += errw97[0]*w97[0]*w97[0]*Holder[0];
		Stuff.yErr += w97[0]*errw97[0]*w97[0]*Holder[0];
		Stuff.zErr += w97[0]*errw97[0]*w97[0]*Holder[0];
		Stuff.Err += errw97[0]*errw97[0]*errw97[0]*Holder[0];
	}

	Stuff.Int = Stuff.Int*Stuff.Area()/8.;		//final weights for size and abs() for error estimates
	Stuff.xErr = abs(Stuff.xErr)*Stuff.Area()/8.;
	Stuff.yErr = abs(Stuff.yErr)*Stuff.Area()/8.;
	Stuff.zErr = abs(Stuff.zErr)*Stuff.Area()/8.;
	Stuff.Err = abs(Stuff.Err)*Stuff.Area()/8.;

	return;
}

long double* k_Int(long double Par[], int Temp, long double theta, int & Stop_Num)
{
	//Extra boundaries that insert extra intervals around peaks. Used a machine learn algorithm of sorts to minimize error to pick these values.
	long double Range[] = {-.421, 0, .421};	//Number of gamma from center

	int Poles;		//Number of poles
	long double zero[26];	//The real part of the signular pole
	long double gamma[26];	//The distance to the singular, maybe
	long double Max;
	int i, j, l;		//Counters, would use 'k', but 'k' is occupied by relative 3-momenta in other parts of program
	int Intervals;		//Number of intervals recorded in Stops

	Characterize_k_Int(Par, Temp, theta, zero, gamma, Poles);	//Find the location of the complex poles
	long double* Stops = new long double[Poles*17+12];		//List of pre-determined subintervals

	l = 0;
	for(i = 0; i < Poles; i++)	//Counting through the poles
	{
		if(!isnan(gamma[i]))
			for(j = 0; j < 3; j++)
			{
				Stops[l] = zero[i]+gamma[i]*Range[j];	//Add all intervals for simultanous on-shell
				l++;
			}
		else	//If a bad pole, at least get the central value for it
		{
			Stops[l] = zero[i];
			l++;
		}
	}

	//More intervals from features not already considered
	Max = Stops[l] = .5*sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//k for which quarks are simultanous light-like, highest k needed for vacuum
	if(isnan(Stops[l]))	//If meson is space-like, keep absolute value of it anyways even though it probably does nothing
		Stops[l] = .5*sqrt(-Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));
	Stops[l+1] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]-pow(2.*Par[2],2)+pow(Par[3]*cos(theta),2)));	//On-shells leaving the positive energy range
	Stops[l+2] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]-pow(2.*Par[2],2)+pow(Par[3]*cos(theta),2)));
	Stops[l+3] = sqrt(4.*pow(Par[3],4)+8.*pow(Par[3],2)*Par[4]+4.*pow(Par[4],2)-pow(Par[1],4))/pow(256.*pow(Par[3],4)+512.*pow(Par[3],2)*Par[4]+256.*pow(Par[4],2),(long double).25);	//Potiential leaving the positive energy range
	Stops[l+4] = abs((pow(Par[2],2)*Par[3]*cos(theta)+sqrt((Par[4]+pow(Par[3],2))*(pow(Par[2],4)+(Par[4]+pow(Par[3]*sin(theta),2))*(Par[4]-2.*pow(Par[2],2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2))));	//On-shell leaving the time-like range
	Stops[l+5] = abs((pow(Par[2],2)*Par[3]*cos(theta)-sqrt((Par[4]+pow(Par[3],2))*(pow(Par[2],4)+(Par[4]+pow(Par[3]*sin(theta),2))*(Par[4]-2.*pow(Par[2],2)))))/(2.*(Par[4]+pow(Par[3]*sin(theta),2))));
	Stops[l+6] = .5*abs(Par[3]*cos(theta)+sqrt(Par[4]+pow(Par[3]*cos(theta),2)));	//Photon point leaving positive energy range. Not sure what photon point
	Stops[l+7] = .5*abs(Par[3]*cos(theta)-sqrt(Par[4]+pow(Par[3]*cos(theta),2)));
	Stops[l+8] = .5*abs(Par[3]*cos(theta)+sqrt(3.*pow(Par[3],2)+4.*Par[4]+pow(Par[3]*cos(theta),2)));
	Stops[l+9] = .5*abs(Par[3]*cos(theta)-sqrt(3.*pow(Par[3],2)+4.*Par[4]+pow(Par[3]*cos(theta),2)));

	for(i = 0; i < l+10; i++)	//Removes stops that are NaN or bigger than necessary
	{
		if(isnan(Stops[i]))
			Stops[i] = -1;
		else if(isinf(Stops[i]) || Stops[i] > 100)
			Stops[i] = 100;
	}

	mergeSort(Stops, 0, l+9);	//Sort the list of sub-intervals

	i = 0;
	j = 0;
	while(Stops[j] <= 0)	//Skip past negative sub-intervals and form NaN
		j++;
	for(; j < l+10; j++)
	{
		if((i > 0 && Stops[i-1] != Stops[j]) || i == 0)	//Removes duplicates, faster to remove duplicates than to evaluate zero width interval
		{
			Stops[i] = Stops[j];
			i++;
		}
		else if(Stops[j] != Stops[j])
			break;
	}
	Stop_Num = i;

	return(Stops);
}

long double* k0_Int(long double Par[], int Temp, long double k, long double theta, int & Stop_Num)
{
	long double zero[12];	//Real part of poles, up to 2 come from potential and up to 2 come from single quark spectrum
	long double gamma[12];	//Imaginary part of poles
	long double a, b;
	long double Max;
	int Poles = 0;		//Number of poles with real parts between 0 and E
	int i, j, l;		//Counting varibles
	int Intervals;		//Number of intervals required by poles and discontinuities

	Characterize_k0_Int(Par, Temp, k, theta, zero, gamma, Poles);	//Get the poles that I have to be concerned about
	long double* Stops = new long double[Poles*17+6];			//Intervals that are required by integrating near poles

	l = 0;
	for(i = 0; i < Poles; i++)
	{
		if(!isnan(zero[i]))	//At lease insert the central point of the pole if the width isn't properly measured
		{
			Stops[l] = zero[i];
			l++;
		}
	}
	Stops[l] = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;	//Lower light-like edge
	Stops[l+1] = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);	//Upper light-like edge
	Stops[l+2] = Energy(0,Par[3]/2.,k,theta)+sqrt(Par[4]+pow(Par[3],2))/2.;	//Pretty sure this is the negative energy solution of the lower light-like edge
	Stops[l+3] = sqrt(Par[4]+pow(Par[3],2))/2.+Energy(0,Par[3]/2.,-k,theta);	//Pretty sure this is the negative energy solution of the upper light-like edge
	Stops[l+4] = sqrt(Par[4]+pow(Par[3],2))/2.;					//Upper energy boundary (E/2)
	Stops[l+5] = -sqrt(Par[4]+pow(Par[3],2))/2.;					//Lower energy boundary (-E/2)

	if(Temp != 0)
	{
		a = b = -sqrt(Par[4]+pow(Par[3],2))/2.;	//Lower edge for non-vacuum
		Max = sqrt(Par[4]+pow(Par[3],2))/2.;		//Upper edge for non-vacuum
	}
	else
	{
		a = b = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;	//Lower edge for vacuum
		Max = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);	//Upper edge for vacuum
	}

	for(i = 0; i < l+6; i++)
	{
		if(Stops[i] < a)
			Stops[i] = -Stops[i];
	}

	mergeSort(Stops, 0, l+5);	//Sort the subintervals

	i = 0;
	j = 0;
	while(Stops[j] == Stops[j]+1.)	//Remove subintervals that duplicates or below the lower edge
		j++;
	for(; j < l+6; j++)
	{
		if(((i > 0 && Stops[i-1] != Stops[j]) || i == 0))	//Remove dublicates and intervals above the upper edge
		{
			Stops[i] = Stops[j];
			i++;
		}
	}
	Stop_Num = i;	//Record the number of intervals

	return(Stops);
}

void Characterize_k_Int(long double Par[], int Temp, long double theta, long double zero[], long double gamma[], int &Poles)
{
	long double holder;	//Holder for bubble sort at the end
	int i, j, l;		//Counter

	Poles = 2;	//Two hard coded poles
	zero[0] = .5*Par[3]*abs(cos(theta));	//Supposedly this is near intersection of 2 on-shells. I don't recongize it, lines 430-435 does that
	gamma[0] = .05;
	zero[1] = Par[2];
	gamma[1] = Par[1];

	if(Par[4]-pow(2.*Par[2],2) > 0.)
	{
		zero[Poles] = .5*sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)));	//relative 3-momentum for which both quarks are on-shell
		gamma[Poles] = 2.*ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp);	//Twice the self-energy of the quarks. Should be sum of self-energy of the two quarks as I don't think the self-energy of both quarks are equal to each other for P!=0.
		Poles++;
	}

	zero[Poles] = Par[2];	//Use Newtons's method to find intersection of features of interest. This one is the positive pole of the potential and one of the peaks of the propagator, I think its the negative pole of the anti-aligned quark.
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;	//Try again with a different seed in case the first missed
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Plus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epm))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEmp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Emp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, mEpp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	zero[Poles] = Par[2];
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}
	zero[Poles] = 10;
	if(Newton_Method_k(zero[Poles], Par[4], Par[3], theta, Par[2], Par[1], V_Minus, Epp))
	{
		gamma[Poles] = ImSelf_Energy(Par[2], Energy(Par[2], Par[3]/2., zero[Poles], theta), zero[Poles], Par, Temp)+sqrt(complex<long double>(pow(2.*zero[Poles],2),pow(Par[2],2))).imag();
		Poles++;
	}

	for(i = 0; i < Poles; i++)	//Move negative results to positive result. Should probably drop it, but might have ben missed or otherwise caught a feature that is interesting, just with the wrong sign
		zero[i] = abs(zero[i]);

	for(i = Poles-1; i >= 0; i--)	//Bubble sort of a self-written pair object
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

	i = 0;
	for(j = 0; j < Poles; j++)
	{
		if(zero[j] > 1000)	//Don't bother with points beyond 1 TeV. All integrals should end well before then
			break;
		if(((i > 0 && zero[i-1] != zero[j]) || i == 0) && !isnan(zero[j]))	//Remove duplicates and NaN
		{
			zero[i] = zero[j];
			gamma[i] = gamma[j];
			i++;
		}
	}
	Poles = i;

	return;
}

bool Newton_Method_k(long double& k, long double s, long double P, long double theta, long double M, long double Lambda, long double(*V)(long double, long double), long double(*k0)(long double, long double, long double, long double, long double))	//Returns the k-intesection of a potiential and on-shell peak
{
	long double newk;
	const long double h = 1e-4;	//Size of finite difference
	bool Success = true;
	int i = 0;

	newk = k - 2.*h*(V(k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(k-h, Lambda)+V(k+h, Lambda)-k0(s, P, k+h, theta, M));	//First iteration of Newton's method with finite differences for derivatives

	while(abs(1.-newk/k) > 1e-5 && i <= 10)	//Allow up to 12 iteration steps, but stop early if last step was small
	{
		k = newk;
		newk = k - 2.*h*(V(k, Lambda)-k0(s, P, k, theta, M))/(k0(s, P, k-h, theta, M)-V(k-h, Lambda)+V(k+h, Lambda)-k0(s, P, k+h, theta, M));
		i++;
	}

	if(abs(1.-newk/k) > 1e-2 || newk < 0 || newk != newk)	//Soundness of value (positive and not NaN) and degree of improvement on last interation (last iteration should have been small)
		Success = false;
	if(abs(1.-V(newk, Lambda)/k0(s, P, newk, theta, M)) > 1e-2)	//Closeness to solution
		Success = false;

	k = newk;

	return(Success);	//Note success or failure of solution find
}

long double V_Plus(long double k, long double Lambda)
{
#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	return(k);
#elif VERSION == 22
	return(k);
#elif VERSION == 24
	return(k);
#elif VERSION == 42
	return(.5*sqrt(complex<long double>(4.*(pow(k, 2)), pow(Lambda, 2))).real());
#endif
}

long double V_Minus(long double k, long double Lambda)
{
#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	return(-k);
#elif VERSION == 22
	return(-k);
#elif VERSION == 24
	return(-k);
#elif VERSION == 42
	return(-.5*sqrt(complex<long double>(4.*(pow(k, 2)), pow(Lambda, 2))).real());
#endif
}

long double Emm(long double s, long double P, long double k, long double theta, long double M)	//peak of the vacuum on-shells
{
	return(sqrt(s+pow(P,2))/2.-Energy(M,P/2.,-k,theta));
}

long double Epm(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.+Energy(M,P/2.,-k,theta));
}

long double mEmp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P,2))/2.-Energy(M,P/2.,k,theta));
}

long double Emp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.-Energy(M,P/2.,k,theta));
}

long double mEpp(long double s, long double P, long double k, long double theta, long double M)
{
	return(-sqrt(s+pow(P,2))/2.+Energy(M,P/2.,k,theta));
}

long double Epp(long double s, long double P, long double k, long double theta, long double M)
{
	return(sqrt(s+pow(P,2))/2.+Energy(M,P/2.,k,theta));
}

long double Upper_Bound(long double s, long double P, long double k, long double theta, long double M)	//Vacuum boundaries
{
	return(sqrt(s+pow(P,2))/2.-Energy(0,P/2.,-k,theta));
}

long double Lower_Bound(long double s, long double P, long double k, long double theta, long double M)
{
	return(Energy(0,P/2.,k,theta)-sqrt(s+pow(P,2))/2.);
}

void Characterize_k0_Int(long double Par[], int Temp, long double k, long double theta, long double zero[10], long double gamma[10], int &Poles)
{
	long double Lower, Upper;	//Limits of integration in k0_Int, vacuum limits are much smaller
	long double holder;
	int i, j;

	if(Temp != 0)
	{
		Lower = -sqrt(Par[4]+pow(Par[3],2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3],2))/2.;	//Integrate from -E/2 to E/2
	}
	else
	{
		Lower = Energy(0,Par[3]/2.,k,theta)-sqrt(Par[4]+pow(Par[3],2))/2.;
		Upper = sqrt(Par[4]+pow(Par[3],2))/2.-Energy(0,Par[3]/2.,-k,theta);
	}

#if VERSION == EXP	//use option -D VERSION={Exp, 22, 24, 42} to select one of the potentials
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
#elif VERSION == 22
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
#elif VERSION == 24
	zero[0] = k;	//Potential poles, I know exactly where these are at.
	zero[1] = -k;
#elif VERSION == 42
	zero[0] = .5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).real();	//Potential poles, I know exactly where these are at.
	zero[1] = -.5*sqrt(complex<long double>(4.*pow(k, 2), pow(Par[1], 2))).real();
#endif

	zero[2] = .5*(sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));	//Exact vacuum poles
	zero[3] = .5*(sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[4] = .5*(-sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));
	zero[5] = .5*(-sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))));

	if(Temp != 0)	//media estimate
	{
		zero[2] = Newton_Method_k0(zero[2], Par, k, theta, Temp, 2, Imk0_Integrand);
		zero[3] = Newton_Method_k0(zero[3], Par, k, theta, Temp, 2, Imk0_Integrand);
		zero[4] = Newton_Method_k0(zero[4], Par, k, theta, Temp, 1, Imk0_Integrand);
		zero[5] = Newton_Method_k0(zero[5], Par, k, theta, Temp, 1, Imk0_Integrand);
	}

	for(i = 0; i < 6; i++)
	{
		if(zero[i] < Lower)
			zero[i] = -zero[i];
	}

	for(i = 5; i >= 0; i--)	//Bubble sort
	{
		for(j = 0; j < i; j++)
		{
			if(zero[j] > zero[j+1])
			{
				holder = zero[j+1];
				zero[j+1] = zero[j];
				zero[j] = holder;
			}
		}
	}

	Poles = 6;

	return;
}

long double Newton_Method_k0(long double k0, long double Par[], long double k, long double theta, int Temp, int j, long double (*Folding)(long double[], long double, long double, long double, int, int))	//Newton's method for finding poles of f by looking for zeros of 1/f, much more stable to the point of absolute confidence
{
	long double new_k0;
	long double h;	//Finite difference
	long double Exit;
	int i = 0;
	long double Danger[] = {.5*(sqrt(Par[4]+pow(Par[3],2))+real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)-k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4)))))),.5*(-sqrt(Par[4]+pow(Par[3],2))-real(sqrt(complex<long double>(4.*(pow(k,2)+pow(Par[2],2)+k*Par[3]*cos(theta))+pow(Par[3],2)-2.*pow(GAMMA,2),2.*sqrt(4.*pow(Par[2]*GAMMA,2)-pow(GAMMA,4))))))};

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

	new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp, j)-1./Folding(Par, k0-h, k, theta, Temp, j))/((1./Folding(Par, k0-h, k, theta, Temp, j)-2./Folding(Par, k0, k, theta, Temp, j)+1./Folding(Par, k0+h, k, theta, Temp, j)));	//First iteration of Netwon's method using finite differences

	while(/*abs(1.-new_k0/k0) > Exit && */i < 100)	//Allow up to 12 iterations to find the pole
	{
		k0 = new_k0;
		new_k0 = k0 - .5*h*(1./Folding(Par, k0+h, k, theta, Temp, j)-1./Folding(Par, k0-h, k, theta, Temp, j))/((1./Folding(Par, k0-h, k, theta, Temp, j)-2./Folding(Par, k0, k, theta, Temp, j)+1./Folding(Par, k0+h, k, theta, Temp, j)));
		i++;
	}

	return(k0);
}

long double omega_Width(long double zero, long double Par[], long double k, long double theta, int Temp, int i, long double (*Folding)(long double[], long double, long double, long double, int, int))	//Breit-Wigner width of the peak
{
	const long double h = 1e-13;
	return(sqrt(abs(24.*h*h*Folding(Par, zero, k, theta, Temp, i)/(Folding(Par, zero-2.*h, k, theta, Temp, i)+16.*Folding(Par, zero-h, k, theta, Temp, i)-30.*Folding(Par, zero, k, theta, Temp, i)+16.*Folding(Par, zero+h, k, theta, Temp, i)-1.*Folding(Par, zero+2.*h, k, theta, Temp, i)))));
}

void ImSelf_Energy(long double M, long double omega[], long double k[], long double Par[], int Temp, long double Results[])	//Single quark self energy for both quarks
{
	static long double omega0[2];		//Location of central peak
	static long double Sigma[2];		//Amplitude of energy dependance
	static long double a[2], b[2];	//Slope of exponential decrease to left and right
	static long double knee[2];		//Interval to change from left to right side of peak
	static long double M_T, Shift;	//Default quark mass, shfift from default quark mass to given quark mass
	static long double k_old[2];		//Previous value of k to know if the parmeters need to recalculated

	if(pow(omega[0],2)>=pow(k[0],2) && omega[0] >= 0)	//Vacuum width
		Results[0] = sqrt(pow(omega[0],2)-pow(k[0],2))*GAMMA;
	else
		Results[0] = 0;
	if(pow(omega[1],2)>=pow(k[1],2) && omega[1] >= 0)
		Results[1] = sqrt(pow(omega[1],2)-pow(k[1],2))*GAMMA;
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
				Sigma[0] = .569969/sqrt(pow(k[0],2)+pow(1.75236,2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1],2)+pow(1.75236,2))+.0187484;
				a[0] = 4.689/(pow(k[0],2)+pow(1.18,2))+4.59495;
				a[1] = 4.689/(pow(k[1],2)+pow(1.18,2))+4.59495;
				b[0] = -70400/(pow(k[0]+20,2)+pow(130,2))+6.24;
				b[1] = -70400/(pow(k[1]+20,2)+pow(130,2))+6.24;
				omega0[0] = sqrt(pow(1.51443+Shift,2)+pow(k[0],2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift,2)+pow(k[1],2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1.,(long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1.,(long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .569969/sqrt(pow(k[0],2)+pow(1.75236,2))+.0187484;
				Sigma[1] = .569969/sqrt(pow(k[1],2)+pow(1.75236,2))+.0187484;
				a[0] = 12.5349/(pow(k[0],2)+pow(1.63711,2))+5.026;
				a[1] = 12.5349/(pow(k[1],2)+pow(1.63711,2))+5.026;
				b[0] = -291.579/(pow(k[0]+15.2519,2)+pow(.0614821,2))+3.36681;
				b[1] = -291.579/(pow(k[1]+15.2519,2)+pow(.0614821,2))+3.36681;
				omega0[0] = sqrt(pow(1.51443+Shift,2)+pow(k[0],2))+.232841;
				omega0[1] = sqrt(pow(1.51443+Shift,2)+pow(k[1],2))+.232841;
				knee[0] = 3.78956*pow(k[0]+1.,(long double)-.530289)+.305*(tanh((k[0]-48.4)/11.1111)+1);
				knee[1] = 3.78956*pow(k[1]+1.,(long double)-.530289)+.305*(tanh((k[1]-48.4)/11.1111)+1);
				break;
			case 2://285MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .625855/sqrt(pow(k[0],2)+pow(1.8429,2))+.0249334;
				Sigma[1] = .625855/sqrt(pow(k[1],2)+pow(1.8429,2))+.0249334;
				a[0] = 3.3971/(pow(k[0],2)+pow(1.01744,2))+3.99561;
				a[1] = 3.3971/(pow(k[1],2)+pow(1.01744,2))+3.99561;
				b[0] = -65187.5/(pow(k[0]+3.11711,2)+pow(101.697,2))+8.15532;
				b[1] = -65187.5/(pow(k[1]+3.11711,2)+pow(101.697,2))+8.15532;
				omega0[0] = sqrt(pow(1.5065+Shift,2)+pow(k[0],2))+.209135;
				omega0[1] = sqrt(pow(1.5065+Shift,2)+pow(k[1],2))+.209135;
				knee[0] = 3.1568*pow(k[0]+1.,(long double)-.624827)+.197004*(tanh((k[0]-27.1743)/10.0192)+1);
				knee[1] = 3.1568*pow(k[1]+1.,(long double)-.624827)+.197004*(tanh((k[1]-27.1743)/10.0192)+1);
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .587509/sqrt(pow(k[0],2)+pow(1.84447,2))+.0309251;
				Sigma[1] = .587509/sqrt(pow(k[1],2)+pow(1.84447,2))+.0309251;
				a[0] = 2.44943/(pow(k[0],2)+pow(.887313,2))+3.32859;
				a[1] = 2.44943/(pow(k[1],2)+pow(.887313,2))+3.32859;
				b[0] = -4439.38/(pow(k[0]-7.23198,2)+pow(38.9387,2))+4.55531;
				b[1] = -4439.38/(pow(k[1]-7.23198,2)+pow(38.9387,2))+4.55531;
				omega0[0] = sqrt(pow(1.47725+Shift,2)+pow(k[0],2))+.219181;
				omega0[1] = sqrt(pow(1.47725+Shift,2)+pow(k[1],2))+.219181;
				knee[0] = 3.28564*pow(k[0]+1.,(long double)-.721321)+.330483*(tanh((k[0]-22.9096)/10.7139)+1);
				knee[1] = 3.28564*pow(k[1]+1.,(long double)-.721321)+.330483*(tanh((k[1]-22.9096)/10.7139)+1);
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .459303/sqrt(pow(k[0],2)+pow(1.84321,2))+.0386564;
				Sigma[1] = .459303/sqrt(pow(k[1],2)+pow(1.84321,2))+.0386564;
				a[0] = 1.79149/(pow(k[0],2)+pow(.764836,2))+2.66209;
				a[1] = 1.79149/(pow(k[1],2)+pow(.764836,2))+2.66209;
				b[0] = -1856.16/(pow(k[0]-8.69519,2)+pow(26.3551,2))+3.94631;
				b[1] = -1856.16/(pow(k[1]-8.69519,2)+pow(26.3551,2))+3.94631;
				omega0[0] = sqrt(pow(1.45428+Shift,2)+pow(k[0],2))+.197493;
				omega0[1] = sqrt(pow(1.45428+Shift,2)+pow(k[1],2))+.197493;
				knee[0] = 3.06296*pow(k[0]+1.,(long double)-.917081)+.394833*(tanh((k[0]-19.5932)/12.0494)+1);
				knee[1] = 3.06296*pow(k[1]+1.,(long double)-.917081)+.394833*(tanh((k[1]-19.5932)/12.0494)+1);
				break;
		}
	}

	long double ImSigma[2];	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
	if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] < -4.)
		ImSigma[0] = a[0]*(omega[0]-omega0[0]+knee[0]/sqrt(a[0]*b[0]));
	else if((omega[0]-omega0[0]+knee[0]*(b[0]-a[0])/(sqrt(a[0]*b[0])*(a[0]+b[0])))/knee[0] > 4.)
		ImSigma[0] = b[0]*(omega0[0]-omega[0]+knee[0]/sqrt(a[0]*b[0]));
	else	//Lost of precision having been circumvented, the actual value
		ImSigma[0] = -.5*((a[0]-b[0])*omega0[0]-((a[0]+b[0])*knee[0])/sqrt(a[0]*b[0]))+(a[0]-b[0])*omega[0]/2-sqrt(pow(((a[0]+b[0])/2.)*(omega[0]-omega0[0]+((a[0]-b[0])*knee[0])/(sqrt(a[0]*b[0])*(a[0]+b[0]))),2)+pow(knee[0],2));

	if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] < -4.)
		ImSigma[1] = a[1]*(omega[1]-omega0[1]+knee[1]/sqrt(a[1]*b[1]));
	else if((omega[1]-omega0[1]+knee[1]*(b[1]-a[1])/(sqrt(a[1]*b[1])*(a[1]+b[1])))/knee[1] > 4.)
		ImSigma[1] = b[1]*(omega0[1]-omega[1]+knee[1]/sqrt(a[1]*b[1]));
	else
		ImSigma[1] = -.5*((a[1]-b[1])*omega0[1]-((a[1]+b[1])*knee[1])/sqrt(a[1]*b[1]))+(a[1]-b[1])*omega[1]/2-sqrt(pow(((a[1]+b[1])/2.)*(omega[1]-omega0[1]+((a[1]-b[1])*knee[1])/(sqrt(a[1]*b[1])*(a[1]+b[1]))),2)+pow(knee[1],2));

#ifdef HALF
	Results[0] += -M*Sigma[0]*exp(ImSigma[0]);	//ImSigma from the in-medium
	Results[1] += -M*Sigma[1]*exp(ImSigma[1]);
#else
	Results[0] += -2.*M*Sigma[0]*exp(ImSigma[0]);
	Results[1] += -2.*M*Sigma[1]*exp(ImSigma[1]);
#endif
	return;
}

long double ImSelf_Energy(long double M, long double omega, long double k, long double Par[], int Temp)	//Single quark self energy
{

	long double omega0;	//location of central peak
	long double Sigma;	//size of energy dependance
	long double a, b;	//slope of exponential decrease to left and right
	long double knee;	//space to change from left to right side of peak
	long double M_T, Shift;
	long double answer;

	if(pow(omega,2)>=pow(k,2) && omega >= 0)
		answer = sqrt(pow(omega,2)-pow(k,2))*GAMMA;
	else
		answer = 0;

	if(Temp == 0)
		return(answer);

	switch(Temp)
	{
		/*case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k,2)+pow(1.75236,2))+.0187484;
			a = 4.689/(pow(k,2)+pow(1.18,2))+4.59495;
			b = -70400/(pow(k+20,2)+pow(130,2))+6.24;
			omega0 = sqrt(pow(1.51443+Shift,2)+pow(k,2))+.232841;
			knee = 3.78956*pow(k+1.,(long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;*/
		case 1://194MeV
			M_T = 1.84184;
			Shift = M-M_T;
			Sigma = .569969/sqrt(pow(k,2)+pow(1.75236,2))+.0187484;
			a = 12.5349/(pow(k,2)+pow(1.63711,2))+5.026;
			b = -291.579/(pow(k+15.2519,2)+pow(.0614821,2))+3.36681;
			omega0 = sqrt(pow(1.51443+Shift,2)+pow(k,2))+.232841;
			knee = 3.78956*pow(k+1.,(long double)-.530289)+.305*(tanh((k-48.4)/11.1111)+1);
			break;
		case 2://285MeV
			M_T = 1.69584;
			Shift = M-M_T;
			Sigma = .625855/sqrt(pow(k,2)+pow(1.8429,2))+.0249334;
			a = 3.3971/(pow(k,2)+pow(1.01744,2))+3.99561;
			b = -65187.5/(pow(k+3.11711,2)+pow(101.697,2))+8.15532;
			omega0 = sqrt(pow(1.5065+Shift,2)+pow(k,2))+.209135;
			knee = 3.1568*pow(k+1.,(long double)-.624827)+.197004*(tanh((k-27.1743)/10.0192)+1);
			break;
		case 3://320MeV
			M_T = 1.59439;
			Shift = M-M_T;
			Sigma = .587509/sqrt(pow(k,2)+pow(1.84447,2))+.0309251;
			a = 2.44943/(pow(k,2)+pow(.887313,2))+3.32859;
			b = -4439.38/(pow(k-7.23198,2)+pow(38.9387,2))+4.55531;
			omega0 = sqrt(pow(1.47725+Shift,2)+pow(k,2))+.219181;
			knee = 3.28564*pow(k+1.,(long double)-.721321)+.330483*(tanh((k-22.9096)/10.7139)+1);
			break;
		case 4://400MeV
			M_T = 1.48038;
			Shift = M-M_T;
			Sigma = .459303/sqrt(pow(k,2)+pow(1.84321,2))+.0386564;
			a = 1.79149/(pow(k,2)+pow(.764836,2))+2.66209;
			b = -1856.16/(pow(k-8.69519,2)+pow(26.3551,2))+3.94631;
			omega0 = sqrt(pow(1.45428+Shift,2)+pow(k,2))+.197493;
			knee = 3.06296*pow(k+1.,(long double)-.917081)+.394833*(tanh((k-19.5932)/12.0494)+1);
			break;
		case 5://40MeV
			M_T = 1.8;
			Shift = M-M_T;
			Sigma = .00386564;
			a = 6.2;
			b = 2.8;
			omega0 = sqrt(pow(1.53+Shift,2)+pow(k,2));
			knee = .56;
			break;
	}

	long double ImSigma;	//Calculation of the argument to the exponential, these first 2 are approximations to avoid catastrophic loss of precision
	if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee < -4.)
		ImSigma = a*(omega-omega0+knee/sqrt(a*b));
	else if((omega-omega0+knee*(b-a)/(sqrt(a*b)*(a+b)))/knee > 4.)
		ImSigma = b*(omega0-omega+knee/sqrt(a*b));
	else
		ImSigma = -.5*((a-b)*omega0-((a+b)*knee)/sqrt(a*b))+(a-b)*omega/2-sqrt(pow(((a+b)/2.)*(omega-omega0+((a-b)*knee)/(sqrt(a*b)*(a+b))),2)+pow(knee,2));

#ifdef HALF
	answer += -M*Sigma*exp(ImSigma);
#else
	answer += -2.*M*Sigma*exp(ImSigma);
#endif

	return(answer);
}

void ReSelf_Energy(long double M, long double omega[], long double k[], int Temp, long double Results[])	//Single quark self energy
{
	static long double Sigma[2];		//Strength
	static long double x0[2], x1[2];	//Centrality markers
	static long double gamma[2];		//Width
	static long double Shift, M_T;
	static long double k_old[2];		//Note on validity of k

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
				Sigma[0] = .257498/sqrt(pow(k[0],2)+pow(1.33201,2))+.00762638;
				Sigma[1] = .257498/sqrt(pow(k[1],2)+pow(1.33201,2))+.00762638;
				x0[0] = sqrt(pow(k[0],2)+pow(1.54778+Shift,2))+.276509;
				x0[1] = sqrt(pow(k[1],2)+pow(1.54778+Shift,2))+.276509;
				x1[0] = sqrt(pow(k[0],2)+pow(1.49799+Shift,2))+.246719;
				x1[1] = sqrt(pow(k[1],2)+pow(1.49799+Shift,2))+.246719;
				gamma[0] = .658734/sqrt(pow(k[0],2)+pow(3.35217,2))+.0815109;
				gamma[1] = .658734/sqrt(pow(k[1],2)+pow(3.35217,2))+.0815109;
				break;*/
			case 1://194MeV
				M_T = 1.84184;
				Shift = M-M_T;
				Sigma[0] = .212571/sqrt(pow(k[0],2)+pow(1.17821,2))+.00762638;
				Sigma[1] = .212571/sqrt(pow(k[1],2)+pow(1.17821,2))+.00762638;
				x0[0] = sqrt(pow(k[0],2)+pow(1.57536+Shift,2))+.259147;
				x0[1] = sqrt(pow(k[1],2)+pow(1.57536+Shift,2))+.259147;
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194+Shift,2))+.222526;
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194+Shift,2))+.222526;
				gamma[0] = .336699/sqrt(pow(k[0],2)+pow(1.87956,2))+.0651449;
				gamma[1] = .336699/sqrt(pow(k[1],2)+pow(1.87956,2))+.0651449;
				break;
			case 2://258MeV
				M_T = 1.69584;
				Shift = M-M_T;
				Sigma[0] = .307972/sqrt(pow(k[0],2)+pow(1.41483,2))+.0101423;
				Sigma[1] = .307972/sqrt(pow(k[1],2)+pow(1.41483,2))+.0101423;
				x0[0] = sqrt(pow(k[0],2)+pow(1.56476+Shift,2))+.251031;
				x0[1] = sqrt(pow(k[1],2)+pow(1.56476+Shift,2))+.251031;
				x1[0] = sqrt(pow(k[0],2)+pow(1.50194+Shift,2))+.222526;
				x1[1] = sqrt(pow(k[1],2)+pow(1.50194+Shift,2))+.222526;
				gamma[0] = .550628/sqrt(pow(k[0],2)+pow(2.43968,2))+.0981269;
				gamma[1] = .550628/sqrt(pow(k[1],2)+pow(2.43968,2))+.0981269;
				break;
			case 3://320MeV
				M_T = 1.59439;
				Shift = M-M_T;
				Sigma[0] = .339131/sqrt(pow(k[0],2)+pow(1.43308,2))+.0125796;
				Sigma[1] = .339131/sqrt(pow(k[1],2)+pow(1.43308,2))+.0125796;
				x0[0] = sqrt(pow(k[0],2)+pow(1.55034+Shift,2))+.257788;
				x0[1] = sqrt(pow(k[1],2)+pow(1.55034+Shift,2))+.257788;
				x1[0] = sqrt(pow(k[0],2)+pow(1.46999+Shift,2))+.231821;
				x1[1] = sqrt(pow(k[1],2)+pow(1.46999+Shift,2))+.231821;
				gamma[0] = .615278/sqrt(pow(k[0],2)+pow(2.22298,2))+.143376;
				gamma[1] = .615278/sqrt(pow(k[1],2)+pow(2.22298,2))+.143376;
				break;
			case 4://400MeV
				M_T = 1.48038;
				Shift = M-M_T;
				Sigma[0] = .304841/sqrt(pow(k[0],2)+pow(1.42911,2))+.0157245;
				Sigma[1] = .304841/sqrt(pow(k[1],2)+pow(1.42911,2))+.0157245;
				x0[0] = sqrt(pow(k[0],2)+pow(1.55511+Shift,2))+.231105;
				x0[1] = sqrt(pow(k[1],2)+pow(1.55511+Shift,2))+.231105;
				x1[0] = sqrt(pow(k[0],2)+pow(1.44714+Shift,2))+.20956;
				x1[1] = sqrt(pow(k[1],2)+pow(1.44714+Shift,2))+.20956;
				gamma[0] = .862629/sqrt(pow(k[0],2)+pow(2.67193,2))+.189598;
				gamma[1] = .862629/sqrt(pow(k[1],2)+pow(2.67193,2))+.189598;
				break;
		}
	}

#ifdef HALF
	Results[0] = Sigma[0]*(omega[0]-x0[0])/(pow(omega[0]-x1[0],2)+gamma[0])/2.;
	Results[1] = Sigma[1]*(omega[1]-x0[1])/(pow(omega[1]-x1[1],2)+gamma[1])/2.;
#else
	Results[0] = Sigma[0]*(omega[0]-x0[0])/(pow(omega[0]-x1[0],2)+gamma[0]);
	Results[1] = Sigma[1]*(omega[1]-x0[1])/(pow(omega[1]-x1[1],2)+gamma[1]);
#endif
	return;
}

void Self_Energy(long double M, long double omega[], long double k[], long double Par[], int Temp, long double ImSelf[], long double ReSelf[])	//Single quark self energy for both quarks. This one has both imaginary and real parts. It is a simple Breit-Wigner peak and simplier than the other provisioned version
{
	static long double omega0[2];	//location of central peak
	static long double Sigma[2];	//size of energy dependance
	static long double gamma[2];	//space to change from left to right side of peak
	static long double k_old[2];

	if(pow(omega[0],2)>=pow(k[0],2))
		ImSelf[0] = sqrt(pow(omega[0],2)-pow(k[0],2))*GAMMA;
	else
		ImSelf[0] = 0;
	if(pow(omega[1],2)>=pow(k[1],2))
		ImSelf[1] = sqrt(pow(omega[1],2)-pow(k[1],2))*GAMMA;
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
				Sigma[0] = .840172/sqrt(pow(k[0],2)+pow(1.45603,2))+.021257;
				Sigma[1] = .840172/sqrt(pow(k[1],2)+pow(1.45603,2))+.021257;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.99829,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.99829,2)+pow(k[1],2));
				gamma[0] = 1.05035*pow(k[0]+1.3891,(long double)-1.3891)+.01;
				gamma[1] = 1.05035*pow(k[1]+1.3891,(long double)-1.3891)+.01;
				break;
			case 2://285MeV
				Sigma[0] = 1.05337/sqrt(pow(k[0],2)+pow(1.50861,2))+.0282696;
				Sigma[1] = 1.05337/sqrt(pow(k[1],2)+pow(1.50861,2))+.0282696;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.97732,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.97732,2)+pow(k[1],2));
				gamma[0] = 1.4624*pow(k[0]+2.64,(long double)-1.41048)+.01;
				gamma[1] = 1.4624*pow(k[1]+2.64,(long double)-1.41048)+.01;
				break;
			case 3://320MeV
				Sigma[0] = 1.14064/sqrt(pow(k[0],2)+pow(1.54999,2))+.0350631;
				Sigma[1] = 1.14064/sqrt(pow(k[1],2)+pow(1.54999,2))+.0350631;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.96823,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.96823,2)+pow(k[1],2));
				gamma[0] = 2.07102*pow(k[0]+3.037,(long double)-1.46076)+.01;
				gamma[1] = 2.07102*pow(k[1]+3.037,(long double)-1.46076)+.01;
				break;
			case 4://400MeV
				Sigma[0] = 1.06073/sqrt(pow(k[0],2)+pow(1.64912,2))+.0438288;
				Sigma[1] = 1.06073/sqrt(pow(k[1],2)+pow(1.64912,2))+.0438288;
				//omega0[0] = sqrt(pow(M,2)+pow(k[0],2));
				//omega0[1] = sqrt(pow(M,2)+pow(k[1],2));
				omega0[0] = sqrt(pow(1.93309,2)+pow(k[0],2));
				omega0[1] = sqrt(pow(1.93309,2)+pow(k[1],2));
				gamma[0] = 3.42222*pow(k[0]+3.663,(long double)-1.56165)+.01;
				gamma[1] = 3.42222*pow(k[1]+3.663,(long double)-1.56165)+.01;
				break;
		}
	}

#ifdef HALF
	ImSelf[0] += -M*Sigma[0]*omega[0]*gamma[0]/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)));
	ImSelf[1] += -M*Sigma[1]*omega[1]*gamma[1]/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)));
	ReSelf[0] += Sigma[0]*(omega[0]-omega0[0])/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)))/2.;
	ReSelf[1] += Sigma[1]*(omega[1]-omega0[1])/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)))/2.;
#else
	ImSelf[0] += -2.*M*Sigma[0]*omega[0]*gamma[0]/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)));
	ImSelf[1] += -2.*M*Sigma[1]*omega[1]*gamma[1]/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)));
	ReSelf[0] += Sigma[0]*(omega[0]-omega0[0])/(M_PI*(pow(omega[0]-omega0[0],2)+pow(omega[0]*gamma[0],2)));
	ReSelf[1] += Sigma[1]*(omega[1]-omega0[1])/(M_PI*(pow(omega[1]-omega0[1],2)+pow(omega[1]*gamma[1],2)));
#endif

	return;
}

long double Energy(long double M, long double P, long double k, long double theta)	//Single quark energy, can return momentum if M=0
{
	if(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta) < 0)
		return(0.);
	else
		return(sqrt(pow(M,2)+pow(P,2)+pow(k,2)+2.*P*k*cos(theta)));
}

long double Set_Temp(int T)
{
	const long double Temps[] = {0,.194,.258,.32,.4,.04,.04};
	return(Temps[T]);
}

long double Fermi(long double omega, int T)	//Fermi factor
{
	static long double Temp = Set_Temp(T);

	if(Temp == 0)
	{
		if(omega >= 0)	//Fermi factor for vacuum
			return(0);
		else
			return(1);
	}
	return(1./(1.+exp(omega/Temp)));
}

long double Potential1(long double Par[], long double k0, long double k)	//Single vertex of potiential without coupling constant
{
#if VERSION == Exp
	return(exp(-abs(-4.*pow(k0,2)+4.*pow(k,2))/pow(Par[1],2)));
#elif VERSION == 22
	return(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))));
#elif VERSION == 42
	return(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)));
#elif VERSION == 24
	return(pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#endif
}

long double Potential2(long double Par[], long double k0, long double k)	//Two vertices of potential with coupling constant
{
#if VERSION == Exp
	return(Par[0]*exp(-2.*abs(-4.*pow(k0,2)+4.*pow(k,2))/pow(Par[1],2)));
#elif VERSION == 22
	return(Par[0]*pow(pow(Par[1],2.)/(pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),2));
#elif VERSION == 42
	return(Par[0]*pow(pow(Par[1],4.)/(pow(Par[1],4.)+pow(-4.*pow(k0,2)+4.*pow(k,2),2)),2));
#elif VERSION == 24
	return(Par[0]*pow(2.*pow(Par[1],2.)/(2.*pow(Par[1],2.)+abs(-4.*pow(k0,2)+4.*pow(k,2))),4));
#endif
}

long double Non_Interacting_Trace(long double Par[], long double k0, long double k , long double theta)
{
	return((Par[4]/4.+pow(k,2)-pow(k0,2)+pow(Par[2],2))/(pow(Par[2],2)));
}

long double Interacting_Linear_Trace(long double Par[], long double k0, long double k , long double theta)
{
	return(sqrt(3.*Par[4]/(8.*pow(Par[2],2))));
}

long double Interacting_Quad_Trace(long double Par[], long double k0, long double k , long double theta)
{
	return(Par[4]/4.-pow(k0,2)+pow(k,2))/(2.*pow(Par[2],2));
}

long double Imk0_Integrand(long double Par[], long double k0, long double k, long double theta, int Temp, int Factor)	//Integrand of the folding integral for positive energy
{
	static long double q[2] = {Energy(0, Par[3]/2., k, theta),Energy(0, Par[3]/2., -k, theta)};
	static long double k_old = k;
	long double omega[2] = {sqrt(Par[4]+pow(Par[3],2))/2.+k0,sqrt(Par[4]+pow(Par[3],2))/2.-k0};
	long double fermi[2] = {Fermi(omega[0], Temp),Fermi(omega[1], Temp)};
	long double ImSelf[2];
	long double ReSelf[2];
	complex<long double> M(Par[2],0);
	complex<long double> gamma(GAMMA,0);

	if(k_old != k)
	{
		k_old = k;
		q[0] = Energy(0, Par[3]/2., k, theta);
		q[1] = Energy(0, Par[3]/2., -k, theta);
	}

	//Self_Energy(Par[2], omega, q, Par, Temp, ImSelf, ReSelf);
	ImSelf_Energy(Par[2], omega, q, Par, Temp, ImSelf);
	ReSelf_Energy(Par[2], omega, q, Temp, ReSelf);

	switch(Factor)
	{
	default:
	case 0:
		return(-((4.*ImSelf[0]*ImSelf[1]*pow(Par[2],2)*(1.-fermi[0]-fermi[1]))/((pow(pow(omega[0],2)-pow(q[0],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[0],2)+pow(ImSelf[0],2))*(pow(pow(omega[1],2)-pow(q[1],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[1],2)+pow(ImSelf[1],2)))));
		break;
	case 1:
		return(ImSelf[0]/(pow(pow(omega[0],2)-pow(q[0],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[0],2)+pow(ImSelf[0],2)));
		break;
	case 2:
		return(ImSelf[1]/(pow(pow(omega[1],2)-pow(q[1],2)-pow(Par[2],2)-2.*Par[2]*ReSelf[1],2)+pow(ImSelf[1],2)));
		break;//*/
	}
}

/*long double Imk0_Integrand(long double Par[], long double k0, long double k, long double theta, int Temp, int Factor)	//Integrand of the folding integral for positive energy
{
	static long double q[2] = {Energy(0, Par[3]/2., k, theta),Energy(0, Par[3]/2., -k, theta)};
	static long double k_old = k;
	long double omegap[2] = {sqrt(Par[4]+pow(Par[3],2))/2.+k0,sqrt(Par[4]+pow(Par[3],2))/2.-k0};
	long double omegan[2] = {-omegap[0],-omegap[1]};
	long double fermi[2] = {Fermi(omegap[0], Temp),Fermi(omegap[1], Temp)};
	long double ImSelfp[2];
	long double ReSelfp[2];
	long double ImSelfn[2];
	long double ReSelfn[2];

	if(k_old != k)
	{
		k_old = k;
		q[0] = Energy(0, Par[3]/2., k, theta);
		q[1] = Energy(0, Par[3]/2., -k, theta);
	}

	//Self_Energy(Par[2], omega, q, Par, Temp, ImSelf, ReSelf);
	ImSelf_Energy(Par[2], omegap, q, Par, Temp, ImSelfp);
	ReSelf_Energy(Par[2], omegap, q, Temp, ReSelfp);
	ImSelf_Energy(Par[2], omegan, q, Par, Temp, ImSelfn);
	ReSelf_Energy(Par[2], omegan, q, Temp, ReSelfn);

	switch(Factor)
	{
	default:
	case 0:
		return(-((4.*(ImSelfp[0]+ImSelfn[0])*(ImSelfp[1]+ImSelfn[1])*pow(Par[2],2)*(1.-fermi[0]-fermi[1]))/((pow(pow(omegap[0],2)-pow(q[0],2)-pow(Par[2],2)-2.*Par[2]*(ReSelfp[0]-ReSelfn[0]),2)+pow(ImSelfn[0]+ImSelfp[0],2))*(pow(pow(omegap[1],2)-pow(q[1],2)-pow(Par[2],2)-2.*Par[2]*(ReSelfp[1]-ReSelfn[1]),2)+pow(ImSelfp[1]+ImSelfn[1],2)))));
		break;
	case 1:
		return((ImSelfp[0]+ImSelfn[0])/(pow(pow(omegap[0],2)-pow(q[0],2)-pow(Par[2],2)-2.*Par[2]*(ReSelfp[0]-ReSelfn[0]),2)+pow(ImSelfp[0]+ImSelfn[0],2)));
		break;
	case 2:
		return((ImSelfp[1]+ImSelfn[1])/(pow(pow(omegap[1],2)-pow(q[1],2)-pow(Par[2],2)-2.*Par[2]*(ReSelfp[1]-ReSelfn[1]),2)+pow(ImSelfp[1]+ImSelfn[1],2)));
		break;
	}
}*/
