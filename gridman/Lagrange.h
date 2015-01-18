#include<iostream>
#include<cmath>
#define DIM 1

#ifndef LAGRANGE_H
#define LAGRANGE_H
class Lagrange
{

	private:

		int porder;
		double *points;
		int support;

	public:

		void assign_lagrange_points(double *p,int mypoint,int n);
		double find_value_at_x(double x);
		double find_der_value_at_x(double x);

};
#endif
