#ifndef REGION
#define REGION

class Region	//region of integral
{
	public:
	long double x1, x2;				//x coordinates of the region
	long double y1, y2;				//y coordinates of the region
	long double z1, z2;				//z coordinates of the region
	Elements<long double> Int;			//Integral of the region
	Elements<long double> xErr, yErr, zErr;	//Error estimate in x and y direction. I'm not sure if they actually means anything.
	Elements<long double> Err;			//Error estimate of the integral
	int order;					//Order of integation to be used in the region
	int xDeep, yDeep, zDeep;			//How deep into the recurssion the region is used
	bool operator>(const Region&) const;		//Greater than operator, both based off of the error
	bool operator<(const Region&) const;		//Less than operator
	long double Area() const;			//Returns the area of the rectangle
	Region();					//Default constructor
	Region(long double, long double, long double, long double, long double, long double, int);	//Constructor given (x1,y1), (x2,y2) of rectangle and default order
	Region(const Region&);				//Copy constructor
	void operator=(const Region&);		//Assignment operator
};

Region::Region()
{
	xDeep = 0;
	yDeep = 0;
	zDeep = 0;
	order = 37;
}

Region::Region(long double a, long double b, long double c, long double d,  long double e, long double f, int precision = 37)
{
	xDeep = 0;
	yDeep = 0;
	zDeep = 0;
	x1 = a;
	x2 = b;
	y1 = c;
	y2 = d;
	z1 = e;
	z2 = f;
	order = precision;
}

Region::Region(const Region& A)
{
	xDeep = A.xDeep;
	yDeep = A.yDeep;
	zDeep = A.zDeep;
	x1 = A.x1;
	x2 = A.x2;
	y1 = A.y1;
	y2 = A.y2;
	z1 = A.z1;
	z2 = A.z2;
	order = A.order;
	Int = A.Int;
	xErr = A.xErr;
	yErr = A.yErr;
	zErr = A.zErr;
	Err = A.Err;
}

void Region::operator=(const Region& x)
{
	x1 = x.x1;
	x2 = x.x2;
	y1 = x.y1;
	y2 = x.y2;
	z1 = x.z1;
	z2 = x.z2;
	Int = x.Int;
	xErr = x.xErr;
	yErr = x.yErr;
	zErr = x.zErr;
	Err = x.Err;
	order = x.order;
	xDeep = x.xDeep;
	yDeep = x.yDeep;
	zDeep = x.zDeep;
}

bool Region::operator>(const Region& x) const
{
	//return((Err/Int)/Area()>(x.Err/x.Int)/x.Area());
	//return(Err/Int>x.Err/x.Int);
	//return(Err>x.Err);
	return(Err/Area()>x.Err/x.Area());
}

bool Region::operator<(const Region& x) const
{
	//return((Err/Int)/Area()<(x.Err/x.Int)/x.Area());
	//return(Err/Int<x.Err/x.Int);
	//return(Err<x.Err);
	return(Err/Area()<x.Err/x.Area());
}

long double Region::Area() const
{
	return((x2-x1)*(y2-y1)*(z2-z1));
}

#endif
