#include<complex>
#include<cstdlib>
using namespace std;

class Elements
{
	public:
		Elements operator+(Elements);
		void operator+=(const Elements &);
		bool operator>=(long double);
		bool operator==(long double);
		bool operator>(long double);
		Elements operator-(Elements);
		Elements operator/(Elements);	//This is about accuracy, not the correct definition of division
		Elements operator*(long double);
		Elements operator*(complex<long double>);
		Elements operator/(long double);
		Elements operator+(long double);
		Elements operator-(long double);
		Elements abs();
		//long double abs(long double&);
		void operator=(const Elements&);
		Elements();
		Elements(long double, long double, long double);
		Elements(complex<long double>, complex<long double>, complex<long double>);
		Elements(long double[3]);
		Elements(complex<long double>[3]);
		Elements(const Elements&);
		void null();
		complex<long double> store(int);
	private:
		complex<long double> Array[3];
};

/*long double Elements::abs(long double& A)
{
	return(A<0?-A:A);
}*/

complex<long double> Elements::store(int i)
{
	return(Array[i]);
}

void Elements::null()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
}

Elements::Elements(long double A, long double B, long double C)
{
	Array[0] = A;
	Array[1] = B;
	Array[2] = C;
}

Elements::Elements(complex<long double> A, complex<long double> B, complex<long double> C)
{
	Array[0] = A;
	Array[1] = B;
	Array[2] = C;
}

Elements::Elements(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
}

Elements::Elements(long double A[3])
{
	Array[0] = A[0];
	Array[1] = A[1];
	Array[2] = A[2];
}

Elements::Elements(complex<long double> A[3])
{
	Array[0] = A[0];
	Array[1] = A[1];
	Array[2] = A[2];
}

Elements::Elements()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
}

Elements Elements::abs()
{
	Elements B;
	B.Array[0] = std::abs(Array[0]);
	B.Array[1] = std::abs(Array[1]);
	B.Array[2] = std::abs(Array[2]);
	return(B);
}

bool Elements::operator==(long double A)
{
	return(Array[0] == A &&
		Array[1] == A &&
		Array[2] == A);
}

bool Elements::operator>=(long double A)
{
	return(std::abs(Array[0]) >= A ||
		std::abs(Array[1]) >= A ||
		std::abs(Array[2]) >= A);
}

bool Elements::operator>(long double A)
{
	return(std::abs(Array[0]) > A ||
		std::abs(Array[1]) > A ||
		std::abs(Array[2]) > A);
}

void Elements::operator=(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
}

void Elements::operator+=(const Elements &A)
{
	#pragma omp atomic
	Array[0] += A.Array[0];
	#pragma omp atomic
	Array[1] += A.Array[1];
	#pragma omp atomic
	Array[2] += A.Array[2];
}

Elements Elements::operator+(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] + A.Array[0];
	B.Array[1] = Array[1] + A.Array[1];
	B.Array[2] = Array[2] + A.Array[2];
	return(B);
}

Elements Elements::operator/(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] / A.Array[0];
	B.Array[1] = Array[1] / A.Array[1];
	B.Array[2] = Array[2] / A.Array[2];
	return(B);
}

Elements Elements::operator-(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] - A.Array[0];
	B.Array[1] = Array[1] - A.Array[1];
	B.Array[2] = Array[2] - A.Array[2];
	return(B);
}

Elements Elements::operator+(long double A)
{
	Elements B;
	B.Array[0] = Array[0] + A;
	B.Array[1] = Array[1] + A;
	B.Array[2] = Array[2] + A;
	return(B);
}

Elements Elements::operator-(long double A)
{
	Elements B;
	B.Array[0] = Array[0] - A;
	B.Array[1] = Array[1] - A;
	B.Array[2] = Array[2] - A;
	return(B);
}

Elements Elements::operator*(long double A)
{
	Elements B;
	B.Array[0] = Array[0] * A;
	B.Array[1] = Array[1] * A;
	B.Array[2] = Array[2] * A;
	return(B);
}

Elements Elements::operator*(complex<long double> A)
{
	Elements B;
	B.Array[0] = Array[0] * A;
	B.Array[1] = Array[1] * A;
	B.Array[2] = Array[2] * A;
	return(B);
}

Elements Elements::operator/(long double A)
{
	Elements B;
	B.Array[0] = Array[0] / A;
	B.Array[1] = Array[1] / A;
	B.Array[2] = Array[2] / A;
	return(B);
}

Elements abs(Elements A)
{
	return(A.abs());
}
