#ifndef INTERPOLATION
#define INTERPOLATION
long double i_k(long double, long double, long double, long double);			//Conversion from k to i counter that may be needed with the interpolation.
long double i_k(long double, long double, long double, long double, long double);	//Interpolation doesn't need it to be an interpolation, therefor their not methods of the class.
long double i_k_wrap(long double, long double[], long double);			//Wrapper for i_k() to select the correct one.

template <class T>
class Interpolation
{
	public:
		~Interpolation();				//Deconstructor
		Interpolation();				//Default Constructor
		Interpolation(T**, int xSize, int ySize);	//Constructor with long double array
		T operator()(long double x, long double y);	//Interpolation evaluation
		void operator=(Interpolation<T>);		//Assignment operator
		bool is_ready();
	private:
		bool ready;
		T** control_points;
		int*** offset;
		long double xRange, yRange;
		long double Basis0(long double);
		long double Basis1(long double);
		long double Basis2(long double);
		long double Basis3(long double);
		long double Basisn(long double);
		long double Basis_Wrapper(long double, int);
};

template <class T>
Interpolation<T>::~Interpolation()
{
	if(ready)
	{
		for(int i = 0; i <= xRange; i++)
		{
			for(int j = 0; j < yRange; j++)
			{
				delete offset[i][j];
			}
			delete offset[i];
			delete control_points[i];
		}
		delete offset;
		delete control_points;
	}
}

template <class T>
Interpolation<T>::Interpolation()
{
	ready = false;
}

template <class T>
Interpolation<T>::Interpolation(T** Control, int xSize, int ySize)	//I really wanted to derive the control points myself, but it is a lot easier to scrape them from Mathematica than figure out the matrix coefficents and inversion.
{
	xRange = xSize-1;
	yRange = ySize-1;

	control_points = new T*[xSize];
	offset = new int**[xSize];
	for(int i = 0; i < xSize; i++)
	{
		control_points[i] = new T[ySize];
		offset[i] = new int*[ySize];
		for(int j = 0; j < ySize; j++)
		{
			control_points[i][j] = Control[i][j];
			offset[i][j] = new int[2];
			offset[i][j][0] = i-5;
			offset[i][j][1] = j-5;
	}	}

	ready = true;
}

template <class T>
void Interpolation<T>::operator=(Interpolation<T> A)
{
	xRange = A.xRange;
	yRange = A.yRange;

	control_points = new T*[int(xRange)+1];
	offset = new int**[int(xRange)+1];
	for(int i = 0; i <= xRange; i++)
	{
		control_points[i] = new T[int(yRange)+1];
		offset[i] = new int*[int(yRange)+1];
		for(int j = 0; j <= yRange; j++)
		{
			control_points[i][j] = A.control_points[i][j];
			offset[i][j] = new int[2];
			offset[i][j][0] = i-5;
			offset[i][j][1] = j-5;
	}	}

	ready = true;
}

template <class T>
T Interpolation<T>::operator()(long double x, long double y)
{
	int Basisi[4] = {4,4,4,4};		//Basis Functions in the i direction
	int Basisj[4] = {4,4,4,4};		//Basis Functions in the j direction
	long double zx[4][4], zy[4][4];	//z from the x direction and y direction
	int offset_i = -1, offset_j = -1;	//Index offsets in calling up control points, normally -1, but can be other values on the ends
	T Answer = T(0);

	if(x < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisi[0] = 0;
		Basisi[1] = 1;
		Basisi[2] = 2;
		Basisi[3] = 3;
		if(x < 1)
			offset_i = 0;
	}
	else if(x < 3)
	{
		Basisi[0] = 1;
		Basisi[1] = 2;
		Basisi[2] = 3;
	}
	else if(x < 4)
	{
		Basisi[0] = 2;
		Basisi[1] = 3;
	}
	else if(x < 5)
	{
		Basisi[0] = 3;
	}
	else if(xRange-2 <= x)
	{
		Basisi[0] = 3;
		Basisi[1] = 2;
		Basisi[2] = 1;
		Basisi[3] = 0;
	}
	else if(xRange-3 <= x)
	{
		Basisi[1] = 3;
		Basisi[2] = 2;
		Basisi[3] = 1;
	}
	else if(xRange-4 <= x)
	{
		Basisi[2] = 3;
		Basisi[3] = 2;
	}
	else if(xRange-5 <= x)
	{
		Basisi[3] = 3;
	}

	if(y < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisj[0] = 0;
		Basisj[1] = 1;
		Basisj[2] = 2;
		Basisj[3] = 3;
		if(y < 1)
			offset_j = 0;
	}
	else if(y < 3)
	{
		Basisj[0] = 1;
		Basisj[1] = 2;
		Basisj[2] = 3;
	}
	else if(y < 4)
	{
		Basisj[0] = 2;
		Basisj[1] = 3;
	}
	else if(y < 5)
	{
		Basisj[0] = 3;
	}
	else if(yRange-2 <= y)
	{
		Basisj[0] = 3;
		Basisj[1] = 2;
		Basisj[2] = 1;
		Basisj[3] = 0;
	}
	else if(yRange-3 <= y)
	{
		Basisj[1] = 3;
		Basisj[2] = 2;
		Basisj[3] = 1;
	}
	else if(yRange-4 <= y)
	{
		Basisj[2] = 3;
		Basisj[3] = 2;
	}
	else if(yRange-5 <= y)
	{
		Basisj[3] = 3;
	}

	if(int(x)+3+offset_i >= xRange+1)
		offset_i = xRange-4-int(x)-offset_i;
	if(int(y)+3+offset_j >= yRange+1)
		offset_j = yRange-4-int(y)-offset_j;

	for(int i_count = 3; i_count >= 0; i_count--)	//Evaluate the Basis Functions
		for(int j_count = 3; j_count >= 0; j_count--)
		{
			if(Basisj[j_count] != 4 && y < 5)
				zy[i_count][j_count] = Basis_Wrapper(y, Basisj[j_count]);
			else if(Basisj[j_count] != 4 && yRange-5 <= y)
				zy[i_count][j_count] = Basis_Wrapper(yRange-y, Basisj[j_count]);
			else
				zy[i_count][j_count] = Basis_Wrapper(y-offset[int(x)][int(y)][1]-j_count, Basisj[j_count]);	//Some how the fractional part of the arguments were trading 

			if(Basisi[i_count] != 4 && x < 5)
				zx[i_count][j_count] = Basis_Wrapper(x, Basisi[i_count]);
			else if(Basisi[i_count] != 4 && xRange-5 <= x)
				zx[i_count][j_count] = Basis_Wrapper(xRange-x, Basisi[i_count]);
			else
				zx[i_count][j_count] = Basis_Wrapper(x-offset[int(x)][int(y)][0]-i_count, Basisi[i_count]);

			Answer += zx[i_count][j_count]*zy[i_count][j_count]*control_points[int(x)+i_count+offset_i][int(y)+j_count+offset_j];
		}

	return(Answer);
}

template <class T>
long double Interpolation<T>::Basis_Wrapper(long double x, int Basis)
{
	switch(Basis)
	{
	case 0:
		return(Basis0(x));
		break;
	case 1:
		return(Basis1(x));
		break;
	case 2:
		return(Basis2(x));
		break;
	case 3:
		return(Basis3(x));
		break;
	case 4:
	default:
		return(Basisn(x));
		break;
	}
}

template <class T>
long double Interpolation<T>::Basis0(long double x)
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3));
	return(0);
}

template <class T>
long double Interpolation<T>::Basis1(long double x)
{
	if(0 <= x && x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(0 <= x && x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

template <class T>
long double Interpolation<T>::Basis2(long double x)
{
	if(0 <= x && x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(0 <= x && x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(0 <= x && x <= 4)
		return(-pow((x-4.)/2.,3));
	return(0);
}

template <class T>
long double Interpolation<T>::Basis3(long double x)
{
	if(0 <= x && x < 2)
		return(pow(x,3)/24.);
	else if(0 <= x && x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(0 <= x && x < 4)
		return(11.*pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(0 <= x && x <= 5)
		return(-pow((x-5.),3)/6.);
	return(0);
}

template <class T>
long double Interpolation<T>::Basisn(long double x)
{
	if(2 <= x && x < 3)
		return(pow(x-2.,3)/6.);
	else if(2 <= x && x < 4)
		return(-pow(x,3)/2.+5*pow(x,2)-16.*x+50./3.);
	else if(2 <= x && x < 5)
		return(pow(x,3)/2.-7*pow(x,2)+32.*x-142./3.);
	else if(2 <= x && x <= 6)
		return(-pow((x-6.),3)/6.);
	return(0);
}

long double i_k_wrap(long double k, long double Par[], long double theta)
{
	if(Par[4] > pow(Par[2]*2.,2) && sqrt((Par[4]-pow(2.*Par[2],2))*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)))/2. >= .5)
		return(i_k(k, Par[4], Par[3], theta, Par[2]));
	else if(Par[4] > 0 && sqrt(Par[4]*(Par[4]+pow(Par[3],2))/(Par[4]+pow(Par[3]*sin(theta),2)))/2. >= .5)
		return(i_k(k, Par[4], Par[3], theta));
	else
		return(k*10.);
}

long double i_k(long double k, long double s, long double P, long double theta)
{
	return((24000.*k)/(40.*k+100.*sqrt((s*(pow(P,2)+s))/(s+pow(P,2)*pow(sin(theta),2)))-2.*k*sqrt((s*(pow(P,2)+s))/(s+pow(P,2)*pow(sin(theta),2)))+(s*(pow(P,2)+s))/(s+pow(P,2)*pow(sin(theta),2))));
}

long double i_k(long double k, long double s, long double P, long double theta, long double M)
{
	long double on_shell = sqrt((s-pow(2.*M,2))*(s+pow(P,2))/(s+pow(P*sin(theta),2)));
	long double photon = sqrt(s*(s+pow(P,2))/(s+pow(P*sin(theta),2)));
	return(1200.*k*(100.-on_shell+photon)/(200.*k+2.*k*photon+500.*on_shell-12.*k*on_shell+5.*on_shell*photon));
}
#endif
