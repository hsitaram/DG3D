#include"Legendre.h"

Legendre::Legendre(int polyorder)
{
	n = polyorder;

	double *coeffs1,*coeffs2;

	
	coeffs     = new double[n+1];
	coeffs1    = new double[n+1];
	coeffs2    = new double[n+1];
	
	der_coeffs = new double[n];
	der_coeffs_1 = new double[n-1];

	for(int i=0;i<n+1;i++)
	{
		coeffs[i]     = 0.0;
		coeffs1[i]    = 0.0;
		coeffs2[i]    = 0.0;
	}
	
	for(int i=0;i<n;i++)
	{	
		der_coeffs[i] = 0.0;
	}
	for(int i=0;i<n-1;i++)
	{
		der_coeffs_1[i]=0.0;
	}

	coeffs2[0]=1.0;
	coeffs1[1]=1.0;

	for(int i=2;i<(n+1);i++)
	{

		for(int j=0;j<(i+1);j++)
		{
			coeffs[j]=0.0;
		}

		for(int j=0;j<i;j++)
		{
			coeffs[j+1] += coeffs1[j]*(2*double(i)-1)/double(i);
		}

		for(int j=0;j<(i-1);j++)
		{
			coeffs[j] -=  (double(i)-1)/double(i)*coeffs2[j];
		}

		//exchange
		for(int j=0;j<i;j++)
		{
			coeffs2[j] = coeffs1[j];
		}
		for(int j=0;j<(i+1);j++)
		{
			coeffs1[j] = coeffs[j];
		}
	}

	for(int i=0;i<n;i++)
	{
		der_coeffs[i]=(i+1)*coeffs[i+1];
	}


	for(int i=0;i<n-1;i++)
	{
		der_coeffs_1[i]=(i+1)*coeffs2[i+1];
	}
}
//==========================================================================================
double Legendre::evaluate(double *coeffs,int n,double x)
{
	double sum=0.0;

	for(int i=0;i<(n+1);i++)
	{
		sum += coeffs[i]*pow(x,double(i));
	}

	return(sum);
}
//==========================================================================================
double Legendre::findzero(double *coeffs,int n,double guess)
{
	double *coeffs_der;
	int itmax=50;
	double errmin=1e-6;
	double err;
	int it;
	double soln;
	double dsoln;

	double val,derval;

	coeffs_der = new double[n];

	for(int i=0;i<n;i++)
	{
		coeffs_der[i]=(i+1)*coeffs[i+1];
	}

	//newton method
	it   = 0;
	soln = guess;
	err = fabs(evaluate(coeffs,n,soln));
	while( (it < itmax) && (err > errmin) )
	{
		val    = evaluate(coeffs,n,soln);
		derval = evaluate(coeffs_der,n-1,soln);
		dsoln  = -val/derval;	

		soln   += dsoln;

		err = fabs(val);
	}

	return(soln);
}
//==========================================================================================
void Legendre::gauss_lobatto_points(double *points,double *weights,int n)
{
	//n is the number of integration points
	double guess;
	points[0]=-1.0; points[n-1]=1.0;

	//linear case
	if(n==2)
	{
		weights[0]=1.0;
		weights[n-1]=1.0;
	}
	else
	{
		for(int i=1;i<=(n-2);i++)
		{
			guess=cos(PI-double(i)*PI/double(n-1));
			points[i]=findzero(der_coeffs,n-2,guess);
		}

		for(int i=0;i<n;i++)
		{
			weights[i]=2.0/double(n)/double(n-1);
			weights[i]=weights[i]/pow(evaluate(coeffs,n-1,points[i]),2.0);
		}
	}

}
//==========================================================================================
