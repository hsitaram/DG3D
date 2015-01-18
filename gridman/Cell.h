#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>
#include<vector>
#include"Face.h"
#include"Lagrange.h"
#include"Legendre.h"
#define NCPTS 8
#define NCFACES 6
#ifndef CELL_CLASS
#define CELL_CLASS
#define XDIR 0
#define YDIR 1
#define ZDIR 2
class Cell
{
	private:
		double m_volume;
		int m_porder;
		Lagrange *m_N_x;
		Lagrange *m_N_y;
		Lagrange *m_N_z;

		int m_nqpts;
		
		double *m_qpts;
		double *m_qwts;

		void m_master_to_physical(double a,double b,double c,double &x,double &y,double &z);
		void m_getcomponentindices(int index,int &a,int &b,int &c);

		double m_master_linearbasis(double x,int num)
		{
			double result;
			result = (num==0)?0.5*(1-x):0.5*(1+x);
			return(result);
		}
		
		double m_master_linearbasisderiv(int num)
		{
			double result;
			result = (num==0)?-0.5:0.5;
			return(result);
		}
		double m_deriv(int xyznum,int abcnum,double a,double b,double c);
		double m_Jacobian(double Jac[NDIM][NDIM],double a,double b,double c);

		double m_dxda(double a,double b,double c)
		{
			return(m_dphysical_dmaster(0,a,b,c));
		}
		double m_dydb(double a,double b,double c)
		{
			return(m_dphysical_dmaster(1,a,b,c));
		}

		double m_dzdc(double a,double b,double c)
		{
			return(m_dphysical_dmaster(2,a,b,c));
		}

		double m_dphysical_dmaster(int dir,double a,double b,double c);
		
	
	public:
		void m_physical_to_master(double x,double y,double z,double &a,double &b,double &c);
		void getpreviouslayerofpoints(int *faceindices,int *prvslayer);
		double *eltnodes_x;
		double *eltnodes_y;
		double *eltnodes_z;
		double coord[NCPTS][NDIM];
		int nodeids[NCPTS];
		int faceids[NCFACES];
		double centroid[NDIM];

		void setvolume(double v)
		{
			m_volume=v;
		}
		double getvolume(){return m_volume;};
		void setbasis(int order);
		void printcellvtu();
		void computelocalmassmat(int ncomp,double *localmassmat);
		void computelocalconvtermAX(int ncomp,double *localsoln,
				double *localflux,double *localrhs);
		void getfaceflux(int findex,double n[NDIM],double area,double *localsoln,
		double *localvel,double *localfaceflux);
		void getinterfaceflux(int ncomp,int findex,double area,int *faceindices,
		double *localf,double *localfaceflux);
		std::vector<int> getDofsforDirichlet(int findex);

		void findcentroid()
		{
			centroid[0]=0.0; centroid[1]=0.0; centroid[2]=0.0;
			for(int i=0;i<NCPTS;i++)
			{
				centroid[0] += coord[i][0];
				centroid[1] += coord[i][1];
				centroid[2] += coord[i][2];
			}

			centroid[0] *=0.125;
			centroid[1] *=0.125;
			centroid[2] *=0.125;
		}

		double getsolnatxyz(std::vector<double> solnvec,double x,
				double y,double z);
		void invertmat(double M[NDIM][NDIM],double Minv[NDIM][NDIM]);
		void returnlocalfacepoints(int findex,int *indices,double *points);
		void getindicesfromxyz(int n,double *coord,int *indices);

};
#endif
