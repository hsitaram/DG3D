#include<iostream>
#include<cmath>
#define PI 3.14159265

#ifndef LEGENDRE_H
#define LEGENDRE_H
class Legendre
{
	private:
		int n;
		double *coeffs;
		double *der_coeffs;
		double *der_coeffs_1;
		
		double evaluate(double *coeffs,int n,double x);
		double findzero(double *coeffs,int n,double guess); 

	public:
		Legendre(int polyorder);
		void gauss_lobatto_points(double *points,double *weights,int n);

};
#endif
