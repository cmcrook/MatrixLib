#include "complex.hpp"
#include <iostream>
#include <cmath>

Complex operator*(Complex a, Complex b) {
	return Complex(a.real() * b.real() - a.imag() * b.imag(),
			a.real() * b.imag() + a.imag() * b.real());
}

Complex operator+(Complex a, Complex b) {
	return Complex(a.real() + b.real(), a.imag() + b.imag());
}

Complex operator-(Complex a, Complex b) {
	return Complex(a.real() - b.real(), a.imag() - b.imag());
}

Complex operator/(Complex a, Complex b) {
	double denominator = pow(b.real(), 2) - pow(b.imag(), 2);
	Complex numerator = b.conjug()*a;
	return numerator/denominator;
}

Complex operator+(Complex a, double b) {
	return Complex(a.real() + b, a.imag());
}

Complex operator-(Complex a, double b) {
	return Complex(a.real() - b, a.imag());
}

Complex operator*(Complex a, double b) {
	return Complex(a.real() * b, a.imag() * b);
}

Complex operator/(Complex a, double b) {
	return Complex(a.real() / b, a.imag() / b);
}

Complex operator+(double a, Complex b) {
	return b+a;
}

Complex operator-(double a, Complex b) {
	return b-a;
}

Complex operator*(double a, Complex b) {
	return b*a;
}

Complex operator/(double a, Complex b) {
	return b/a;
}

Complex::Complex(double a, double b) {
	this->a = a;
	this->b = b;
}

Complex Complex::conjug() {
	b *= -1;
	return *this;
}

double Complex::real() {
	return a;
}

double Complex::imag() {
	return b;
}

double Complex::mag(){
	return sqrt(pow(a,2) + pow(b,2));
}

void Complex::polar(double& r, double& theta) {
	r = mag();
	theta = atan(b/a);
}
