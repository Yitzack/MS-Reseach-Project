#include<cmath>
using namespace std;

class Elements
{
	public:
		Elements operator+(Elements);
		void operator+=(const Elements &);
		bool operator>=(long double);
		bool operator>(long double);
		Elements operator-(Elements);
		Elements operator/(Elements);
		Elements operator*(long double);
		Elements operator/(long double);
		Elements operator+(long double);
		Elements operator-(long double);
		Elements abs();
		void operator=(const Elements&);
		Elements();
		Elements(long double[5]);
		Elements(const Elements&);
		//long double abs(long double A){return(abs(A));}
		void null();
		long double* store();
	private:
		long double Array[5];
};

long double* Elements::store()
{
	return(Array);
}

void Elements::null()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
	Array[3] = 0;
	Array[4] = 0;
}

Elements::Elements(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
	Array[4] = A.Array[4];
}

Elements::Elements(long double A[5])
{
	Array[0] = A[0];
	Array[1] = A[1];
	Array[2] = A[2];
	Array[3] = A[3];
	Array[4] = A[4];
}

Elements::Elements()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
	Array[3] = 0;
	Array[4] = 0;
}

Elements Elements::abs()
{
	Elements B;
	B.Array[0] = std::abs(Array[0]);
	B.Array[1] = std::abs(Array[1]);
	B.Array[2] = std::abs(Array[2]);
	B.Array[3] = std::abs(Array[3]);
	B.Array[4] = std::abs(Array[4]);
	return(B);
}

bool Elements::operator>=(long double A)
{
	return(Array[0] >= A &&
		Array[1] >= A &&
		Array[2] >= A &&
		Array[3] >= A &&
		Array[4] >= A);
}

bool Elements::operator>(long double A)
{
	return(Array[0] > A &&
		Array[1] > A &&
		Array[2] > A &&
		Array[3] > A &&
		Array[4] > A);
}

void Elements::operator=(const Elements &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
	Array[4] = A.Array[4];
}

void Elements::operator+=(const Elements &A)
{
	Array[0] += A.Array[0];
	Array[1] += A.Array[1];
	Array[2] += A.Array[2];
	Array[3] += A.Array[3];
	Array[4] += A.Array[4];
}

Elements Elements::operator+(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] + A.Array[0];
	B.Array[1] = Array[1] + A.Array[1];
	B.Array[2] = Array[2] + A.Array[2];
	B.Array[3] = Array[3] + A.Array[3];
	B.Array[4] = Array[4] + A.Array[4];
	return(B);
}

Elements Elements::operator/(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] / A.Array[0];
	B.Array[1] = Array[1] / A.Array[1];
	B.Array[2] = Array[2] / A.Array[2];
	B.Array[3] = Array[3] / A.Array[3];
	B.Array[4] = Array[4] / A.Array[4];
	return(B);
}

Elements Elements::operator-(Elements A)
{
	Elements B;
	B.Array[0] = Array[0] - A.Array[0];
	B.Array[1] = Array[1] - A.Array[1];
	B.Array[2] = Array[2] - A.Array[2];
	B.Array[3] = Array[3] - A.Array[3];
	B.Array[4] = Array[4] - A.Array[4];
	return(B);
}

Elements Elements::operator+(long double A)
{
	Elements B;
	B.Array[0] = Array[0] + A;
	B.Array[1] = Array[1] + A;
	B.Array[2] = Array[2] + A;
	B.Array[3] = Array[3] + A;
	B.Array[4] = Array[4] + A;
	return(B);
}

Elements Elements::operator-(long double A)
{
	Elements B;
	B.Array[0] = Array[0] - A;
	B.Array[1] = Array[1] - A;
	B.Array[2] = Array[2] - A;
	B.Array[3] = Array[3] - A;
	B.Array[4] = Array[4] - A;
	return(B);
}

Elements Elements::operator*(long double A)
{
	Elements B;
	B.Array[0] = Array[0] * A;
	B.Array[1] = Array[1] * A;
	B.Array[2] = Array[2] * A;
	B.Array[3] = Array[3] * A;
	B.Array[4] = Array[4] * A;
	return(B);
}

Elements Elements::operator/(long double A)
{
	Elements B;
	B.Array[0] = Array[0] / A;
	B.Array[1] = Array[1] / A;
	B.Array[2] = Array[2] / A;
	B.Array[3] = Array[3] / A;
	B.Array[4] = Array[4] / A;
	return(B);
}
