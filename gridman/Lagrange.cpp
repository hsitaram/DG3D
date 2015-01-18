#include"Lagrange.h"

//================================================================
void Lagrange::assign_lagrange_points(double *p,int mypoint,int n)
{
	porder = n-1;
	points = new double[porder+1];
        support = mypoint;

        for(int i=0;i<(porder+1);i++)	
	{
		points[i]=p[i];
	}
}
//================================================================
double Lagrange::find_value_at_x(double x)
{
	double val;

	val=1;

	for(int i=0;i<(porder+1);i++)
	{
		if(i != support)
		{
			val=val*(x-points[i])/(points[support]-points[i]);
		}
	}

	return(val);	
}
//================================================================
double Lagrange::find_der_value_at_x(double x)
{
	double val,sum,dtr;
	
	sum=0.0;
	dtr=1.0;
	for(int i=0;i<(porder+1);i++)
	{
		if(i != support)
		{
			val=1.0;
			for(int j=0;j<(porder+1);j++)
			{
				
				if( (j != i) && (j != support))
				{
					val=val*(x-points[j]);
				}


			}		

			dtr=dtr*(points[support]-points[i]);
			sum=sum+val;
		}
	}

	return(sum/dtr);

}
//================================================================
