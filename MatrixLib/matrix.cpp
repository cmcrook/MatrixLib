#include "matrix.hpp"
#include <iostream>
#include <cmath>

Matrix operator*(double c, Matrix m) {
	for (int i = 0; i < m.getHeight(); i++) {
		for (int j = 0; j < m.getWidth(); j++) {
			m.setElement(j, i, c * m.getElement(j, i));
		}
	}
	return m;
}

Matrix operator*(Matrix m, double c) {
	return c * m;
}

Matrix::Matrix(int width, int height, double elements[]) {
	this->width = width;
	this->height = height;
	matrix.insert(matrix.begin(), elements, elements + width * height);
}

Matrix::Matrix(int width, int height, std::vector<double> elements) {
	this->width = width;
	this->height = height;
	matrix = elements;
}

//Generates a matrix of specified width and height with each element as val
void Matrix::fill(int width, int height, double val) {
	this->width = width;
	this->height = height;
	matrix.assign(width * height, val);
}

//Generates a matrix of zeros of specified width and height
void Matrix::zeros(int width, int height) {
	fill(width, height, 0);
}

//Generates an identity matrix of specified width and height, if matrix is not square then enters largest possible diagonal with 1s
void Matrix::identity(int width, int height) {
	//Creates blank matrix full of zeros
	zeros(width, height);

	//Adds 1s on largest diagonal
	for (int i = 0; i < ((height <= width) ? height : width); i++) {
		matrix.insert(matrix.begin() + width * i + i, 1);
		matrix.erase(matrix.begin() + width * i + i + 1);
	}
}

//Generates a matrix of ones of specified width and height
void Matrix::ones(int width, int height) {
	fill(width, height, 1);
}

//Returns the width of the matrix
int Matrix::getWidth() {
	return width;
}

//Returns the height of the matrix
int Matrix::getHeight() {
	return height;
}

//Returns the number of elements in the matrix
int Matrix::getSize() {
	return matrix.size();
}

//Returns the value of the element at indices i and j
double Matrix::getElement(int i, int j) {
	int index = width * j + i;
	return matrix.at(index);
}

//Returns a vector of the elements in the matrix
std::vector<double> Matrix::getElements() {
	return matrix;
}

//Returns a matrix
Matrix Matrix::getSubMatrix(int iW, int fW, int iH, int fH) {
	int newWidth = fW - iW + 1;
	int newHeight = fH - iH + 1;
	std::vector<double> elements;
	for (int j = 0; j < newHeight; j++) {
		for (int i = 0; i < newWidth; i++) {
			elements.push_back(getElement(i + iW, j + iH));
		}
	}
	Matrix sub(newWidth, newHeight, elements);
	return sub;
}

//Sets the value of the element at indices i and j to val
void Matrix::setElement(int i, int j, double val) {
	if (i > width - 1) {
		for (int m = 0; m < height; m++) {
			for (int n = 0; n < i - width + 1; n++) {
				int index = width * (m + 1) + (i - width + 1) * m + n;
				matrix.insert(matrix.begin() + index, 0);
			}
		}
		width = i + 1;
	}

	if (j > height - 1) {
		for (int m = 0; m < j - height + 1; m++) {
			for (int n = 0; n < width; n++) {
				int index = width * height + m * width + n;
				matrix.insert(matrix.begin() + index, 0);
			}
		}
		height = j + 1;
	}

	matrix.insert(matrix.begin() + width * j + i, val);
	matrix.erase(matrix.begin() + width * j + i + 1);
}

//
void Matrix::setElements(std::vector<double> elements) {
	matrix = elements;
}

//Inserts a matrix m into the the instance at a upper left index denoted by x and y
//If the matrix does not fit within the size of the instance at the indices x and y then
//An error is printed to the console and nothing is performed
void Matrix::insertSubMatrix(int x, int y, Matrix m) {
	int widthLowerBound = x;
	int widthUpperBound = x + m.width;
	int heightLowerBound = y;
	int heightUpperBound = y + m.height;
	if (widthLowerBound >= 0 && widthUpperBound <= width
			&& heightLowerBound >= 0 && heightUpperBound <= height) {
		for (int i = 0; i < m.getWidth(); i++) {
			for (int j = 0; j < m.getHeight(); j++) {
				setElement(i + x, j + y, m.getElement(i, j));
			}
		}
	} else {
		std::cout
				<< "INSERT_SUB_MATRIX::Index out of bounds or sub matrix is larger than matrix!"
				<< std::endl;
	}
}

//Decomposes matrix into two triangular matrices denoted a L for Lower and U for upper
void Matrix::luDecomp(Matrix& L, Matrix& U) {
	L.identity(width, height);
	U = *this;
	if (width == height) {
		int rWidth = width - 1, rHeight = height - 1;
		for (int i = 1; i < height; i++) {
			double cL = U.getElement(0, i) / U.getElement(0, 0);
			for (int j = 0; j < width; j++) {
				U.setElement(j, i,
						U.getElement(j, i) - cL * U.getElement(j, 0));
			}
			L.setElement(0, i, cL);
		}

		if (rWidth >= 2 && rHeight >= 2) {
			Matrix rL, rU;
			rL.zeros(rWidth, rHeight);
			rU.zeros(rWidth, rHeight);
			U.getSubMatrix(1, rWidth, 1, rHeight).luDecomp(rL, rU);
			L.insertSubMatrix(1, 1, rL);
			U.insertSubMatrix(1, 1, rU);
		}
	}
}

//Decomposes matrix into a upper triangular matrix R and an orthogonal matrix Q
void Matrix::qrDecomp(Matrix& Q, Matrix& R) {
	Matrix U, col;
	Q.zeros(width, height);
	U.zeros(width, height);
	col.zeros(1, height);

	for (int i = 0; i < width; i++) {
		col = getSubMatrix(i, i, 0, height - 1);
		for (int j = 0; j < i; j++) {
			col = col - U.getSubMatrix(j, j, 0, U.getHeight() - 1).proj(
					getSubMatrix(i, i, 0, height - 1));
		}
		Q.insertSubMatrix(i, 0, col.unitize());
		U.insertSubMatrix(i, 0, col);
	}
	R = Q.transpose() * (*this);
}

//Projects the vector m onto this vector
//Only valid if both are same size and one dimensional on a side
Matrix Matrix::proj(Matrix m) {
	Matrix projection;
	Matrix u = *this;

	double p = (u.transpose() * m).getElement(0, 0)
			/ (u.transpose() * u).getElement(0, 0);
	projection = p * u;

	return projection;
}

//Computes the norm of the vector
//Only valid for one dimensional matrices
double Matrix::norm() {
	if ((height == 1) ^ (width == 1)) {
		double sum = 0;
		for (int i = 0; i < ((height > width) ? height : width); i++) {
			sum += pow(matrix.at(i), 2);
		}
		return sqrt(sum);
	} else {
		return -1;
	}
}

//Returns a unit vector of this instance
Matrix Matrix::unitize() {
	if ((height == 1) ^ (width == 1)) {
		double normal = norm();
		Matrix n;
		n.zeros(width, height);
		for (int i = 0; i < ((height > width) ? height : width); i++) {
			n.setElement((width > height) ? i : 0, (height > width) ? i : 0,
					matrix.at(i) / normal);
		}
		return n;
	} else {
		return *this;
	}
}

Matrix Matrix::luSolve(Matrix L, Matrix U, Matrix m) {
	if (m.width == 1) {
		return backwardSub(U, forwardSub(L, m));
	} else {
		std::cout << "System cannot be solved!" << std::endl;
		return *(new Matrix());
	}
}

//Performs forward substitution to solve for in L*z = m where L is a lower triangular matrix
Matrix Matrix::forwardSub(Matrix L, Matrix m) {
	Matrix z;
	z.zeros(1, L.getHeight());

	double sum = 0;

	for (int i = 0; i < L.height; i++) {
		for (int j = 0; j < i; j++) {
			sum += L.getElement(j, i) * z.getElement(0, j);
		}
		z.setElement(0, i, (m.getElement(0, i) - sum) / L.getElement(i, i));
		sum = 0;
	}

	return z;
}

//Performs backwards substitution to solve for x in U*x = m where U is an upper triangular matrix
Matrix Matrix::backwardSub(Matrix U, Matrix m) {
	Matrix x;
	x.zeros(1, U.getHeight());

	double sum = 0;

	for (int i = U.height - 1; i >= 0; i--) {
		for (int j = U.height - 1; j >= i; j--) {
			sum += U.getElement(j, i) * x.getElement(0, j);
		}
		x.setElement(0, i, (m.getElement(0, i) - sum) / U.getElement(i, i));
		sum = 0;
	}

	return x;
}

Matrix Matrix::operator+(Matrix m) {
	if (width == m.getWidth() && height == m.getHeight()) {
		Matrix sum;
		sum.zeros(width, height);

		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				sum.setElement(i, j, getElement(i, j) + m.getElement(i, j));
			}
		}

		return sum;
	} else {
		std::cout << "OPERATOR+::Matrix dimensions do not match!" << std::endl;
		return *(new Matrix());
	}
}

Matrix Matrix::operator-(Matrix m) {
	if (width == m.getWidth() && height == m.getHeight()) {
		Matrix sum;
		sum.zeros(width, height);

		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				sum.setElement(i, j, getElement(i, j) - m.getElement(i, j));
			}
		}

		return sum;
	} else {
		std::cout << "OPERATOR-::Matrix dimensions do not match!" << std::endl;
		return *(new Matrix());
	}
}

Matrix Matrix::operator*(Matrix k) {
	if (width != k.height) {
		std::cout << "OPERATOR*::Matrix dimensions do not agree!" << std::endl;
	}
	std::vector<double> productVec;
	double elementSum = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < k.width; j++) {
			for (int m = 0; m < width; m++) {
				elementSum += this->getElement(m, i) * k.getElement(j, m);
			}
			productVec.insert(productVec.begin() + i * k.width + j, elementSum);
			elementSum = 0;
		}
	}
	Matrix product(k.width, height, productVec);
	return product;
}

Matrix Matrix::operator^(int i) {
	Matrix m = *this;

	if (i == 0) {
		m.identity(width, height);
		return m;
	}

	if (i < 0) {
		m = m.inverse();
		i *= -1;
	}

	Matrix n = m;
	for (int j = 0; j < i - 1; j++) {
		n = n * m;
	}

	return n;
}

void Matrix::operator=(Matrix k) {
	this->width = k.getWidth();
	this->height = k.getHeight();
	setElements(k.getElements());
}

bool Matrix::operator==(Matrix m) {
	if (width == m.getWidth() && height == m.getHeight()) {
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				if (getElement(i, j) != m.getElement(i, j)) {
					return false;
				}
			}
		}
		return true;
	} else {
		return false;
	}
}

Matrix Matrix::transpose() {
	std::vector<double> elements;
	elements.assign(width * height, 0);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			elements.insert(elements.begin() + height * i + j,
					getElement(i, j));
		}
	}
	Matrix product(height, width, elements);
	return product;
}

double Matrix::det() {
	if (width == height) {
		Matrix L, U;
		luDecomp(L, U);

		double productU = 1;

		for (int i = 0; i < L.getWidth(); i++) {
			productU *= U.getElement(i, i);
		}

		return productU;
	} else {
		std::cout << "DET::Matrix sides are not equal!" << std::endl;
		return -1;
	}
}

Matrix Matrix::inverse() {
	double determinant = det();
	if (determinant != 0) {
		Matrix L, U, id, result;
		result.zeros(width, height);
		id.identity(width, height);

		luDecomp(L, U);

		for (int i = 0; i < width; i++) {
			result.insertSubMatrix(i, 0,
					luSolve(L, U, id.getSubMatrix(i, i, 0, height - 1)));
		}

		return result;
	} else {
		std::cout << "INVERSE::Matrix is singular!" << std::endl;
		return *(new Matrix());
	}
}

Matrix Matrix::adj() {
	return (inverse() * det());
}

Matrix Matrix::cofactor() {
	return adj().transpose();
}

//Prints the matrix to the console
void Matrix::printMatrix() {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			std::cout << matrix.at(width * i + j);
			if (j != width - 1) {
				std::cout << " ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
