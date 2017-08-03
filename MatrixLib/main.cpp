#include <iostream>
#include "matrix.hpp"

using namespace std;

int main() {
	double mData[] = {1,2,3,0,4,5,1,0,6};// { 1, 3, 3, 1, 4, 3, 1, 3, 4 };
	Matrix L, U, z, x, id;
	Matrix m(3, 3, mData);
//	m.printMatrix();
//
//	cout << "Det: " << m.det() << endl;
//
//	id.identity(3, 3);
//
//	m.luDecomp(L, U);
//	cout << endl;
//	L.printMatrix();
//	cout << endl;
//	U.printMatrix();
//
//	cout << endl;
//
//	double fData[] = {-4,-6,-15};
//	Matrix f(1, 3, fData);
//
//	z = m.forwardSub(L, id.getSubMatrix(0, 0, 0, 2));
//	x = m.backwardSub(U, z);
//
//	cout << endl;
//	z.printMatrix();
//	cout << endl;
//	x.printMatrix();

//	cout << endl;
//	m.inverse().printMatrix();
//	cout << endl;
//	m.cofactor().printMatrix();

	z = m^(-1);
	z.printMatrix();

	cout << endl;
	m.inverse().printMatrix();
}
