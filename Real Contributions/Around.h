#include<cmath>
//#include"Elements.h"
using namespace std;

#ifndef AROUND
#define AROUND

class Around
{
	public:
		Around();				//Default constructor
		Around(long double);			//Constructor with x+/-0
		Around(long double, long double);	//Constructor with x+/-y
		Around(const Around&, const Around&);	//Constructor with two Arounds x+/-sqrt(deltax^2+y^2), to accomadate error of error integration accumalation
		Around(const Around&);			//Copy constructor
		Around& operator=(Around);		//Assignment
		Around& operator+=(Around);		//Add and assign
		Around operator+(Around);		//Summation
		Around operator+(long double);
		Around operator-(Around);		//Difference
		Around operator-(long double);
		Around operator*(Around);		//Product
		Around operator*(long double);
		Around operator/(Around);		//Quotient with uncertain divisor
		Around operator/(long double);	//Quotient with exact divisor
		Around abs(Around&);			//Absolute value
		Around abs();
		bool operator==(Around);		//Equality test
		bool operator>=(Around);		//Greater than or equal test
		bool operator>(Around);		//Greater than test
		bool operator==(long double);		//Equality to exact value. I wanted to be clever with if the value is within a standard deviation, call it equal. Caused problems in the algorithm cutting off early because small value appeared to be zero when more distance required
		bool operator>=(long double);
		bool operator>(long double);
		bool isnan();				//Is either part of Around() nan?
		long double Value();			//Return value
		long double Error();			//Return error estimate
		friend ostream& operator<<(ostream&, const Around&);	//Write Around to stream conformal to Mathematica standard
	private:
		long double value;
		long double error;
};

Around::Around()
{
	value = 0;
	error = 0;
}

Around::Around(long double V)
{
	value = V;
	error = 0;
}

Around::Around(long double V, long double E)
{
	value = V;
	error = E;
}

Around::Around(const Around& V, const Around& E)
{
	value = V.value;
	error = sqrt(pow(V.error,2)+pow(E.value,2));
}

Around::Around(const Around& A)
{
	value = A.value;
	error = A.error;
}

Around& Around::operator=(Around A)
{
	error = A.error;
	value = A.value;
	return(*this);
}

Around& Around::operator+=(Around A)
{
	error = sqrt(pow(error,2)+pow(A.error,2));
	value += A.value;
	return(*this);
}

Around Around::operator+(Around A)
{
	Around B;
	B.error = sqrt(pow(error,2)+pow(A.error,2));
	B.value = value + A.value;
	return(B);
}

Around Around::operator+(long double A)
{
	Around B;
	B.error = error;
	B.value = value + A;
	return(B);
}

Around operator+(long double A, Around B)	//exact calling sum, turn it around have the uncertain call the sum
{
	return(Around(A,0)+B);
}

Around Around::operator-(Around A)
{
	Around B;
	B.error = sqrt(pow(error,2)+pow(A.error,2));
	B.value = value - A.value;
	return(B);
}

Around Around::operator-(long double A)
{
	Around B;
	B.error = error;
	B.value = value - A;
	return(B);
}

Around operator-(long double A, Around B)
{
	return(Around(A,0)-B);
}

Around Around::operator*(Around A)
{
	Around B;
	B.error = sqrt(pow(A.value*error,2)+pow(value*A.error,2));
	B.value = value * A.value;
	return(B);
}

Around Around::operator*(long double A)
{
	Around B;
	B.error = error * A;
	B.value = value * A;
	return(B);
}

Around operator*(long double A, Around B)
{
	return(B*A);
}

Around Around::operator/(Around A)
{
	Around B;
	B.error = sqrt(pow(error/A.value,2)+pow(value*A.error/pow(A.value,2),2));
	B.value = value / A.value;
	return(B);
}

Around Around::operator/(long double A)
{
	Around B;
	B.error = error / A;
	B.value = value / A;
	return(B);
}

Around operator/(long double A, Around B)	//Quotent with exact dividend.
{
	return(Around(A,0)/B);
}

Around abs(Around A)
{
	return(A.abs());
}

Around Around::abs(Around& A)
{
	return(A.abs());
}

Around Around::abs()
{
	return(Around(std::abs(value), error));
}

bool Around::operator==(Around A)
{
	if(std::abs(A.value-this->value) < sqrt(pow(this->error,2)+pow(A.error,2)))
		return(true);
	return(false);
}

bool Around::operator>(Around A)
{
	if(value > A.value && std::abs(A.value-this->value) > sqrt(pow(this->error,2)+pow(A.error,2)))
		return(true);
	return(false);
}

bool Around::operator>=(Around A)
{
	if(*this > A || *this == A)
		return(true);
	return(false);
}

bool Around::operator==(long double A)
{
	if(value == A)
		return(true);
	return(false);
}

bool Around::operator>(long double A)
{
	if(value > A)
		return(true);
	return(false);
}

bool Around::operator>=(long double A)
{
	if(*this > A || *this == A)
		return(true);
	return(false);
}

bool isnan(Around A)
{
	return(A.isnan());
}

bool Around::isnan()
{
	return(std::isnan(value) || std::isnan(error));
}

long double Around::Value()
{
	return(value);
}

long double Around::Error()
{
	return(error);
}

ostream& operator<<(ostream& os, const Around& A)
{
	os << "Around[" << A.value << "," << A.error << "]" << flush;
	return(os);
}

#endif
