#include<cstdlib>
using namespace std;

class Elements
{
	public:
		Elements();				//Default constructor
		Elements(long double, long double, long double, long double);	//Constructor from 5 elements
		Elements(long double[4]);		//Array constructor
		Elements(const Elements&);		//Copy constructor
		void operator=(const Elements&);	//Assignment
		void operator+=(const Elements &);	//Accumalate and assign
		bool operator==(long double);		//This is looking for all components == 0, not any scalar. So, its actually looking for the 0 vector
		bool operator>=(long double);		//This is about accuracy, so all components must pass
		bool operator>(long double);		//This is about accuracy, so all components must pass
		bool isnan();				//Equavalent to isnan(long double) but for vector but called as A.isnan instead of isnan(A)
		Elements operator+(Elements);		//Sum of vectors
		Elements operator-(Elements);		//Difference of vectors
		Elements operator/(Elements);		//This is about accuracy, not the correct definition of division, so it is an element-wise division
		Elements operator*(long double);	//Scalar multiply
		Elements operator/(long double);	//Scalar divide
		Elements abs();			//Absolute value of all components
		void null();				//Make the element the zero vector
		long double operator[](int);		//Returns the element at int. This is not the correct way to do this as it should return a pointer to the component so it can be altered. I can get away with it as I'm only printing the the contents to an output stream.
	private:
		long double Array[4];			//The vector itself
};

Elements::Elements()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
	Array[3] = 0;
}

Elements::Elements(long double A, long double B, long double C, long double D)
{
	Array[0] = A;
	Array[1] = B;
	Array[2] = C;
	Array[3] = D;
}

Elements::Elements(long double A[4])
{
	Array[0] = A[0];
	Array[1] = A[1];
	Array[2] = A[2];
	Array[3] = A[3];
}

Elements::Elements(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
}

void Elements::operator=(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
}

void Elements::operator+=(const Elements &A)
{
	Array[0] += A.Array[0];
	Array[1] += A.Array[1];
	Array[2] += A.Array[2];
	Array[3] += A.Array[3];
}

bool Elements::operator==(long double A)
{
	return(Array[0] == A &&
		Array[1] == A &&
		Array[2] == A &&
		Array[3] == A);
}

bool Elements::operator>=(long double A)
{
	return(std::abs(Array[0]) >= A ||
		std::abs(Array[1]) >= A ||
		std::abs(Array[2]) >= A ||
		std::abs(Array[3]) >= A);
}

bool Elements::operator>(long double A)
{
	return(std::abs(Array[0]) > A ||
		std::abs(Array[1]) > A ||
		std::abs(Array[2]) > A ||
		std::abs(Array[3]) > A);
}

bool Elements::isnan()
{
	return(std::isnan(Array[0]) ||
		std::isnan(Array[1]) ||
		std::isnan(Array[2]) ||
		std::isnan(Array[3]));
}

Elements Elements::operator+(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] + A.Array[0];
	B.Array[1] = Array[1] + A.Array[1];
	B.Array[2] = Array[2] + A.Array[2];
	B.Array[3] = Array[3] + A.Array[3];
	return(B);
}

Elements Elements::operator-(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] - A.Array[0];
	B.Array[1] = Array[1] - A.Array[1];
	B.Array[2] = Array[2] - A.Array[2];
	B.Array[3] = Array[3] - A.Array[3];
	return(B);
}

Elements Elements::operator/(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] / A.Array[0];
	B.Array[1] = Array[1] / A.Array[1];
	B.Array[2] = Array[2] / A.Array[2];
	B.Array[3] = Array[3] / A.Array[3];
	return(B);
}

Elements Elements::operator*(long double A)
{
	Elements B;
	B.Array[0] = Array[0] * A;
	B.Array[1] = Array[1] * A;
	B.Array[2] = Array[2] * A;
	B.Array[3] = Array[3] * A;
	return(B);
}

Elements Elements::operator/(long double A)
{
	Elements B;
	B.Array[0] = Array[0] / A;
	B.Array[1] = Array[1] / A;
	B.Array[2] = Array[2] / A;
	B.Array[3] = Array[3] / A;
	return(B);
}

Elements abs(Elements A)
{
	return(A.abs());
}

Elements Elements::abs()
{
	Elements B;
	B.Array[0] = std::abs(Array[0]);
	B.Array[1] = std::abs(Array[1]);
	B.Array[2] = std::abs(Array[2]);
	B.Array[3] = std::abs(Array[3]);
	return(B);
}

void Elements::null()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
	Array[3] = 0;
}

long double Elements::operator[](int i)
{
	return(Array[i]);
}
