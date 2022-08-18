long double Basis0(long double x)
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3));
	return(0);
}

long double Basis1(long double x)
{
	if(0 <= x && x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(0 <= x && x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

long double Basis2(long double x)
{
	if(0 <= x && x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(0 <= x && x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(0 <= x && x <= 4)
		return(-pow((x-4.)/2.,3));
	return(0);
}

long double Basis3(long double x)
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

long double Basisn(long double x)
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

template <class T>
class Interpolation
{
	public:
		Interpolation(T**, int xSize, int ySize);	//Constructor with long double array
		T operator()(long double x, long double y);
	private:
		T** control_points;
		int*** offset;
		long double xRange, yRange;
};


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
}

template <class T>
T Interpolation<T>::operator()(long double x, long double y)
{
	long double (*Basisi[4])(long double) = {Basisn, Basisn, Basisn, Basisn}; //Basis Functions in the i direction
	long double (*Basisj[4])(long double) = {Basisn, Basisn, Basisn, Basisn}; //Basis Functions in the i direction
	long double zx[4][4], zy[4][4];					 //z from the x direction and y direction
	int offset_i = -1, offset_j = -1;					 //Index offsets in calling up control points, normally -1, but can be 0 or -2 on the ends
	T Answer = T(0);

	if(x < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisi[0] = Basis0;
		Basisi[1] = Basis1;
		Basisi[2] = Basis2;
		Basisi[3] = Basis3;
		if(x < 1)
			offset_i = 0;
	}
	else if(x < 3)
	{
		Basisi[0] = Basis1;
		Basisi[1] = Basis2;
		Basisi[2] = Basis3;
	}
	else if(x < 4)
	{
		Basisi[0] = Basis2;
		Basisi[1] = Basis3;
	}
	else if(x < 5)
	{
		Basisi[0] = Basis3;
	}
	else if(xRange-2 <= x)
	{
		Basisi[0] = Basis3;
		Basisi[1] = Basis2;
		Basisi[2] = Basis1;
		Basisi[3] = Basis0;
	}
	else if(xRange-3 <= x)
	{
		Basisi[1] = Basis3;
		Basisi[2] = Basis2;
		Basisi[3] = Basis1;
	}
	else if(xRange-4 <= x)
	{
		Basisi[2] = Basis3;
		Basisi[3] = Basis2;
	}
	else if(xRange-5 <= x)
	{
		Basisi[3] = Basis3;
	}

	if(y < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisj[0] = Basis0;
		Basisj[1] = Basis1;
		Basisj[2] = Basis2;
		Basisj[3] = Basis3;
		if(y < 1)
			offset_j = 0;
	}
	else if(y < 3)
	{
		Basisj[0] = Basis1;
		Basisj[1] = Basis2;
		Basisj[2] = Basis3;
	}
	else if(y < 4)
	{
		Basisj[0] = Basis2;
		Basisj[1] = Basis3;
	}
	else if(y < 5)
	{
		Basisj[0] = Basis3;
	}
	else if(yRange-2 <= y)
	{
		Basisj[0] = Basis3;
		Basisj[1] = Basis2;
		Basisj[2] = Basis1;
		Basisj[3] = Basis0;
	}
	else if(yRange-3 <= y)
	{
		Basisj[1] = Basis3;
		Basisj[2] = Basis2;
		Basisj[3] = Basis1;
	}
	else if(yRange-4 <= y)
	{
		Basisj[2] = Basis3;
		Basisj[3] = Basis2;
	}
	else if(yRange-5 <= y)
	{
		Basisj[3] = Basis3;
	}

	if(int(x)+3+offset_i >= xRange+1)
		offset_i = xRange-4-int(x)-offset_i;
	if(int(y)+3+offset_j >= yRange+1)
		offset_j = yRange-4-int(y)-offset_j;

	for(int i_count = 3; i_count >= 0; i_count--)	//Evaluate the Basis Functions
		for(int j_count = 3; j_count >= 0; j_count--)
		{
			if(Basisj[j_count] != Basisn && y < 5)
				zy[i_count][j_count] = Basisj[j_count](y);
			else if(Basisj[j_count] != Basisn && yRange-5 <= y)
				zy[i_count][j_count] = Basisj[j_count](yRange-y);
			else
				zy[i_count][j_count] = Basisj[j_count](y-offset[int(x)][int(y)][1]-j_count);	//Some how the fractional part of the arguments were trading 

			if(Basisi[i_count] != Basisn && x < 5)
				zx[i_count][j_count] = Basisi[i_count](x);
			else if(Basisi[i_count] != Basisn && xRange-5 <= x)
				zx[i_count][j_count] = Basisi[i_count](xRange-x);
			else
				zx[i_count][j_count] = Basisi[i_count](x-offset[int(x)][int(y)][0]-i_count);

			Answer += zx[i_count][j_count]*zy[i_count][j_count]*control_points[int(x)+i_count+offset_i][int(y)+j_count+offset_j];
		}

	return(Answer);
}

