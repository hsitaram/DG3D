#include"Fieldvar.h"
#include"Gridman.h"

#define NCVAR 5
#define MASS 0
#define XMOM 1
#define YMOM 2
#define ZMOM 3
#define ENRG 4
#define G_AIR 1.4

#ifndef EULER_MANAGER_H
#define EULER_MANAGER_H
class Euler_manager
{
	private:
		Gridman *m_mesh;
		int m_porder;
		double m_P0,m_rho0;
		double m_dt;
		
		void m_computefaceterms(); 
		void m_computeconvtermAX();
		void m_computeboundaryterms();
		void m_computeinviscidflux(int n,double *localsoln,
				double *localflux);
		void m_getwavespeed(int n,double normal[3],double *localsoln_l,
		double *localsoln_r,double *localwavespeed,double *localcourantnumber);
		void m_computeLaxFriedrichflux(int n,double normal[3],
				double *wavespeeds,double *localsoln_l,
				double *localsoln_r,double *localfaceflux);
		double m_computepressure(double u[NCVAR]);
		void m_getavgsolnvalue(int cellid,double u[NCVAR]);
		void m_computetvdflux(int n,double normal[3],double *limiter
		,double *localsoln_l,double *localsoln_r,double *localfaceflux);
		double m_vanleerlimiter(double r,double mod_c);

	public:
		Fieldvar *uvar;
		Fieldvar *rhsvar;
		Fieldvar *massmat;
		Fieldvar *uvar_prvs;
		void init(std::string meshfilename,int porder,double dt);
		void computeRHS();
		void computeMassMat();

		Gridman * getmeshptr(){return m_mesh;}

};
#endif
