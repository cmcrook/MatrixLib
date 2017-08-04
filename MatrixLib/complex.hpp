#pragma once
#include <iostream>
#include <vector>

class Complex{
	double a, b;

	public:
	Complex(double, double);
	~Complex(){}

	Complex conjug();
	double real();
	double imag();
	void polar(double&, double&);
	double mag();
};


Complex operator+(Complex, Complex);
Complex operator-(Complex, Complex);
Complex operator*(Complex, Complex);
Complex operator/(Complex, Complex);
Complex operator+(Complex, double);
Complex operator-(Complex, double);
Complex operator*(Complex, double);
Complex operator/(Complex, double);
Complex operator+(double, Complex);
Complex operator-(double, Complex);
Complex operator*(double, Complex);
Complex operator/(double, Complex);
