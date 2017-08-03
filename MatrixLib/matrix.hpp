#pragma once
#include <vector>

class Matrix {
	int width, height;
	std::vector<double> matrix;

	public:
		Matrix(){
			height = 0;
			width = 0;
		};
		Matrix(int, int, double[]);
		Matrix(int, int, std::vector<double>);

		~Matrix(){}

		void fill(int, int, double);
		void zeros(int, int);
		void identity(int, int);
		void ones(int, int);

		int getWidth();
		int getHeight();
		int getSize();
		double getElement(int, int);
		std::vector<double> getElements();
		Matrix getSubMatrix(int, int , int ,int);

		void setElement(int, int, double);
		void setElements(std::vector<double>);
		void insertSubMatrix(int, int, Matrix);

		void luDecomp(Matrix& L, Matrix& U);
		Matrix luSolve(Matrix, Matrix, Matrix);
		Matrix forwardSub(Matrix, Matrix);
		Matrix backwardSub(Matrix, Matrix);

		//Overridden operators: NEED TO FIGURE OUT ORDER INDEPENDENCE
		Matrix operator*(Matrix);
		Matrix operator^(int);
		void operator=(Matrix);

		Matrix transpose();
		double det();
		Matrix inverse();
		Matrix adj();
		Matrix cofactor();

		void printMatrix();
};
