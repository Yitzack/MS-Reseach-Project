#include<cstdlib>
#include"Around.h"
using namespace std;

#ifndef ELEMENTS
#define ELEMENTS

template<class T>
class Elements
{
	public:
		Elements();				//Default constructor
		Elements(T, T, T, T);			//Constructor from 4 elements
		Elements(T[4]);			//Array constructor
		Elements(const Elements&);		//Copy constructor
		void operator=(const Elements&);	//Assignment
		void operator+=(const Elements &);	//Accumalate and assign
		void operator-=(const Elements &);	//Deaccumalate and assign
		bool operator==(T);			//This is looking for all components == 0, not any scalar. So, its actually looking for the 0 vector
		bool operator>=(T);			//This is about accuracy, so all components must pass
		bool operator>(const Elements<T>) const;//This is about accuracy, so all components must pass
		bool operator<(const Elements<T>) const;//This is about accuracy, so all components must pass
		bool operator>(T);			//This is about accuracy, so all components must pass
		bool operator<(T);			//This is about accuracy, so all components must pass
		bool isnan();				//Equavalent to isnan(T) but for vector but called as A.isnan instead of isnan(A)
		Elements<T> operator+(Elements);	//Sum of vectors
		Elements<T> operator+(T);		//Add a number to elements of vector
		Elements<T> operator-(Elements);	//Difference of vectors
		Elements<T> operator-(T);		//Subtract a number from elements of vector
		Elements<T> operator/(Elements<T>) const;//This is about accuracy, not the correct definition of division, so it is an element-wise division
		Elements<T> operator/(T);		//Scalar divide
		Elements<T> operator*(T);		//Scalar multiply
		Elements<T> operator*(Elements);	//Vector multiply (not to be confused with cross product, but element by element multiply
		Elements<T> abs(const Elements<T>&);	//Absolute value of all elements
		Elements<T> abs() const;		//Absolute value of all elements
		void null();				//Make the element the zero vector
		T operator[](int);			//Returns the element at int. This is not the correct way to do this as it should return a pointer to the component so it can be altered. I can get away with it as I'm only printing the the contents to an output stream.
	private:
		T Array[4];			//The vector itself
};

template <class T>
Elements<T>::Elements()
{
	Array[0] = 0;
	Array[1] = 0;
	Array[2] = 0;
	Array[3] = 0;
}

template <class T>
Elements<T>::Elements(T A, T B, T C, T D)
{
	Array[0] = A;
	Array[1] = B;
	Array[2] = C;
	Array[3] = D;
}

template <class T>
Elements<T>::Elements(T A[4])
{
	Array[0] = A[0];
	Array[1] = A[1];
	Array[2] = A[2];
	Array[3] = A[3];
}

template <class T>
Elements<T>::Elements(const Elements<T> &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
}

template <class T>
void Elements<T>::operator=(const Elements<T> &A)
{
	Array[0] = A.Array[0];
	Array[1] = A.Array[1];
	Array[2] = A.Array[2];
	Array[3] = A.Array[3];
}

template <class T>
void Elements<T>::operator+=(const Elements<T> &A)
{
	Array[0] += A.Array[0];
	Array[1] += A.Array[1];
	Array[2] += A.Array[2];
	Array[3] += A.Array[3];
}

template <class T>
void Elements<T>::operator-=(const Elements<T> &A)
{
	Array[0] -= A.Array[0];
	Array[1] -= A.Array[1];
	Array[2] -= A.Array[2];
	Array[3] -= A.Array[3];
}

template <class T>
bool Elements<T>::operator==(T A)
{
	return(Array[0] == A &&
		Array[1] == A &&
		Array[2] == A &&
		Array[3] == A);
}

template <class T>
bool Elements<T>::operator>=(T A)
{
	using std::abs;
	return(abs(Array[0]) >= A ||
		abs(Array[1]) >= A ||
		abs(Array[2]) >= A ||
		abs(Array[3]) >= A);
}

template <class T>
bool Elements<T>::operator>(const Elements<T> A) const
{
	T lhs, rhs;	//The error of the elements remain roughly propitional to each other throught out the execution of the algorithm. The sum of the elements does a better job of ordering the objects than an or operation on the elements individually.
	lhs = Array[0]+Array[1]+Array[2]+Array[3];
	rhs = A.Array[0]+A.Array[1]+A.Array[2]+A.Array[3];
	return(lhs > rhs);
}

template <class T>
bool Elements<T>::operator<(const Elements<T> A) const
{
	T lhs, rhs;	//The error of the elements remain roughly propitional to each other throught out the execution of the algorithm. The sum of the elements does a better job of ordering the objects than an or operation on the elements individually.
	lhs = Array[0]+Array[1]+Array[2]+Array[3];
	rhs = A.Array[0]+A.Array[1]+A.Array[2]+A.Array[3];
	return(lhs < rhs);
}

template <class T>
bool Elements<T>::operator>(T A)
{
	using std::abs;
	return(abs(Array[0]) > A ||
		abs(Array[1]) > A ||
		abs(Array[2]) > A ||
		abs(Array[3]) > A);
}

template <class T>
bool Elements<T>::operator<(T A)
{
	using std::abs;
	return(abs(Array[0]) < A ||
		abs(Array[1]) < A ||
		abs(Array[2]) < A ||
		abs(Array[3]) < A);
}

template <class T>
bool Elements<T>::isnan()
{
	using std::isnan;
	return(isnan(Array[0]) ||
		isnan(Array[1]) ||
		isnan(Array[2]) ||
		isnan(Array[3]));
}

template <class T>
Elements<T> Elements<T>::operator+(Elements<T> A)
{
	Elements<T> B;
	B.Array[0] = Array[0] + A.Array[0];
	B.Array[1] = Array[1] + A.Array[1];
	B.Array[2] = Array[2] + A.Array[2];
	B.Array[3] = Array[3] + A.Array[3];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator+(T A)
{
	Elements<T> B;
	B.Array[0] = Array[0] + A;
	B.Array[1] = Array[1] + A;
	B.Array[2] = Array[2] + A;
	B.Array[3] = Array[3] + A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator-(Elements<T> A)
{
	Elements<T> B;
	B.Array[0] = Array[0] - A.Array[0];
	B.Array[1] = Array[1] - A.Array[1];
	B.Array[2] = Array[2] - A.Array[2];
	B.Array[3] = Array[3] - A.Array[3];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator-(T A)
{
	Elements<T> B;
	B.Array[0] = Array[0] - A;
	B.Array[1] = Array[1] - A;
	B.Array[2] = Array[2] - A;
	B.Array[3] = Array[3] - A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator/(Elements<T> A) const
{
	Elements<T> B;
	B.Array[0] = Array[0] / A.Array[0];
	B.Array[1] = Array[1] / A.Array[1];
	B.Array[2] = Array[2] / A.Array[2];
	B.Array[3] = Array[3] / A.Array[3];
	return(B);
}

template <class T>
Elements<T> operator*(long double A, Elements<T> B)
{
	return(B*A);
}

template <class T>
Elements<T> Elements<T>::operator*(T A)
{
	Elements<T> B;
	B.Array[0] = Array[0] * A;
	B.Array[1] = Array[1] * A;
	B.Array[2] = Array[2] * A;
	B.Array[3] = Array[3] * A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator*(Elements<T> A)
{
	Elements<T> B;
	B.Array[0] = Array[0] * A.Array[0];
	B.Array[1] = Array[1] * A.Array[1];
	B.Array[2] = Array[2] * A.Array[2];
	B.Array[3] = Array[3] * A.Array[3];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator/(T A)
{
	Elements<T> B;
	B.Array[0] = Array[0] / A;
	B.Array[1] = Array[1] / A;
	B.Array[2] = Array[2] / A;
	B.Array[3] = Array[3] / A;
	return(B);
}

template <typename T>
Elements<T> Elements<T>::abs(const Elements<T>& A)
{
	using std::abs;
	return(Elements<T>(abs(A[0]),abs(A[1]),abs(A[2]),abs(A[3])));
}

template <typename T>
Elements<T> Elements<T>::abs() const
{
	using std::abs;
	return(Elements<T>(abs(Array[0]),abs(Array[1]),abs(Array[2]),abs(Array[3])));
}

template <typename T>
Elements<T> abs(const Elements<T>& A)
{
	return(A.abs());
}

template <class T>
void Elements<T>::null()
{
	Array[0] = T(0);
	Array[1] = T(0);
	Array[2] = T(0);
	Array[3] = T(0);
}

template <class T>
T Elements<T>::operator[](int i)
{
	return(Array[i]);
}

#endif
