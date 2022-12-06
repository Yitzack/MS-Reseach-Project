#include<cstdlib>
#include"Around.h"
using namespace std;

#ifndef ELEMENTS
#define ELEMENTS

template<class T>
class Elements
{
	public:
		~Elements();					//Deconstructor
		Elements();					//Default constructor
		Elements(int);					//Null constructor with a size of array
		Elements(T[], int);				//Array constructor
		Elements(const Elements&);			//Copy constructor
		Elements<T> operator=(const Elements&);	//Assignment
		Elements<T> operator+=(const Elements &);	//Accumalate and assign
		Elements<T> operator-=(const Elements &);	//Deaccumalate and assign
		bool operator==(T);				//This is looking for all components == 0, not any scalar. So, its actually looking for the 0 vector
		bool operator>=(T);				//This is about accuracy, so all components must pass
		bool operator>(const Elements<T>) const;	//This is about accuracy, so all components must pass
		bool operator<(const Elements<T>) const;	//This is about accuracy, so all components must pass
		bool operator>(T);				//This is about accuracy, so all components must pass
		bool operator<(T);				//This is about accuracy, so all components must pass
		bool isnan();					//Equavalent to isnan(T) but for vector but called as A.isnan instead of isnan(A)
		Elements<T> concat(const Elements<T>) const;	//Concats two elements together and then returns the result
		Elements<T> operator+(Elements);		//Sum of vectors
		Elements<T> operator+(T);			//Add a number to elements of vector
		Elements<T> operator-(Elements);		//Difference of vectors
		Elements<T> operator-(T);			//Subtract a number from elements of vector
		Elements<T> operator/(Elements<T>) const;	//This is about accuracy, not the correct definition of division, so it is an element-wise division
		Elements<T> operator/(T) const;		//Scalar divide
		Elements<T> operator*(T);			//Scalar multiply
		Elements<T> operator*(Elements);		//Vector multiply (not to be confused with cross product, but element by element multiply
		Elements<T> abs(const Elements<T>&) const;	//Absolute value of all elements
		Elements<T> abs() const;			//Absolute value of all elements
		void null();					//Make Array the zero vector
		int size(){return(Size);};			//returns the size of Array
		T* operator[](int);				//Returns the element at int. This is not the correct way to do this as it should return a pointer to the component so it can be altered. I can get away with it as I'm only printing the the contents to an output stream.
	private:
		T* Array;	//The vector itself
		int Size;	//Size of the vector
};

template <class T>
Elements<T>::~Elements()
{
	Size = 0;
	delete Array;
}

template <class T>
Elements<T>::Elements()
{
	Size = 0;
}

template <class T>
Elements<T>::Elements(int N)
{
	Size = N;
	Array = new T[Size];
}

template <class T>
Elements<T>::Elements(T A[], int N)
{
	Size = N;
	Array = new T[Size];
	for(int i = 0; i < Size; i++)
		Array[i] = A[i];
}

template <class T>
Elements<T>::Elements(const Elements<T> &A)
{
	if(&A != this)
	{
		Size = A.Size;
		Array = new T[Size];
		for(int i = 0; i < Size; i++)
			Array[i] = A.Array[i];
	}
}

template <class T>
Elements<T> Elements<T>::operator=(const Elements<T> &A)
{
	if(&A != this)
	{
		Size = A.Size;
		delete Array;
		Array = new T[Size];
		for(int i = 0; i < Size; i++)
			Array[i] = A.Array[i];
	}
	return(*this);
}

template <class T>
Elements<T> Elements<T>::operator+=(const Elements<T> &A)
{
	for(int i = 0; i < Size; i++)
		Array[i] += A.Array[i];
	return(*this);
}

template <class T>
Elements<T> Elements<T>::operator-=(const Elements<T> &A)
{
	for(int i = 0; i < Size; i++)
		Array[i] -= A.Array[i];
	return(*this);
}

template <class T>
bool Elements<T>::operator==(T A)
{
	int i = 0;
	while(i < Size)
	{
		if(Array[i] != A)
			return(false);
	}
	return(true);
}

template <class T>
bool Elements<T>::operator>=(T A)
{
	using std::abs;
	int i = 0;
	while(i < Size)
	{
		if(abs(Array[i]) >= A)
			return(true);
	}
	return(false);
}

template <class T>
bool Elements<T>::operator>(const Elements<T> A) const
{
	T lhs = T(0), rhs = T(0);	//The error of the elements remain roughly propitional to each other through out the execution of the algorithm. The sum of the elements does a better job of ordering the objects than an || operation on the elements individually.
	for(int i = 0; i < Size; i++)
	{
		lhs += Array[i];
		rhs += A.Array[i];
	}
	return(lhs > rhs);
}

template <class T>
bool Elements<T>::operator<(const Elements<T> A) const
{
	T lhs = T(0), rhs = T(0);	//The error of the elements remain roughly propitional to each other through out the execution of the algorithm. The sum of the elements does a better job of ordering the objects than an || operation on the elements individually.
	for(int i = 0; i < Size; i++)
	{
		lhs += Array[i];
		rhs += A.Array[i];
	}
	return(lhs < rhs);
}

template <class T>
bool Elements<T>::operator>(T A)
{
	using std::abs;
	int i = 0;
	while(i < Size)
	{
		if(abs(Array[i]) > A)
			return(true);
	}
	return(false);
}

template <class T>
bool Elements<T>::operator<(T A)
{
	using std::abs;
	int i = 0;
	while(i < Size)
	{
		if(abs(Array[i]) < A)
			return(true);
	}
	return(false);
}

template <class T>
bool Elements<T>::isnan()
{
	using std::isnan;
	int i = 0;
	while(i < Size)
	{
		if(isnan(Array[i]))
			return(true);
	}
	return(false);
}

template <class T>
Elements<T> Elements<T>::concat(Elements<T> A) const
{
	Elements<T> B(Size+A.Size);
	int i, j;

	for(i = 0; i < Size; i++)
		B.Array[i] = Array[i];

	j = i;
	for(i = 0; i < Size; i++)
	{
		B.Array[j] = A.Array[i];
		j++;
	}

	return(B);
}

template <class T>
Elements<T> Elements<T>::operator+(Elements<T> A)
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] + A.Array[i];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator+(T A)
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] + A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator-(Elements<T> A)
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] - A.Array[i];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator-(T A)
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] - A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator/(Elements<T> A) const
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] / A.Array[i];
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
	Elements<T> B(Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i]*A;
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator*(Elements<T> A)
{
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] * A.Array[i];
	return(B);
}

template <class T>
Elements<T> Elements<T>::operator/(T A) const
{
	Elements<T> B(Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = Array[i] / A;
	return(B);
}

template <typename T>
Elements<T> Elements<T>::abs(const Elements<T>& A) const
{
	using std::abs;
	Elements<T> B(A.Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = abs(Array[i]);
	return(B);
}

template <typename T>
Elements<T> Elements<T>::abs() const
{
	using std::abs;
	Elements<T> B(Size);
	for(int i = 0; i < Size; i++)
		B.Array[i] = abs(Array[i]);
	return(B);
}

template <typename T>
Elements<T> abs(const Elements<T>& A)
{
	return(A.abs());
}

template <class T>
void Elements<T>::null()
{
	for(int i = 0; i < Size; i++)
		Array[i] = 0;
}

template <class T>
T* Elements<T>::operator[](int i)
{
	return(&Array[i]);
}

#endif
