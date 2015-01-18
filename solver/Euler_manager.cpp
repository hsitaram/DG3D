#include"Euler_manager.h"

//====================================================================
void Euler_manager::init(std::string meshfilename,int porder,double dt)
{
	int ncells,ndof,nlocaldof,nqpts;
	int offs,ind;
	double x,y,z;
	double u0,x0,y0,z0,a0,b0,c0;
	double pL,rhoL,pR,rhoR;
	double expval;
	double delta;

	m_porder = porder;
	m_dt     = dt;

	m_mesh = new Gridman;
	m_mesh->readgridfile(meshfilename);
	m_mesh->printmeshvtu();

	ncells = m_mesh->getncells();

	for(int i=0;i<ncells;i++)
	{
		m_mesh->Cvec[i].setbasis(porder);
	}

	nqpts     = (porder+1);
	nlocaldof = nqpts*nqpts*nqpts;
	ndof      = nlocaldof*ncells;

	uvar      = new Fieldvar;
	uvar_prvs = new Fieldvar;
	rhsvar    = new Fieldvar;
	massmat   = new Fieldvar;

	uvar->setvarprops("U",ndof,NCVAR);
	uvar_prvs->setvarprops("U_prvs",ndof,NCVAR);
	rhsvar->setvarprops("RHS",ndof,NCVAR);
	massmat->setvarprops("Mass_mat",ndof,NCVAR);

	//Problem specific initialization
	for(int i=0;i<ndof;i++)
	{
		for(int k=0;k<NCVAR;k++)
		{
			massmat->realarray[NCVAR*i+k]  = 0.0;
			rhsvar->realarray [NCVAR*i+k]  = 0.0;
		}
	}

	offs = 0;
	u0 = 0.1;
	x0 = 0.5; y0   = 0.5;    z0 = 0.5;
	a0 = 0.25; b0  = 0.25;   c0 = 0.25;
	m_P0 = 1.0; m_rho0 = 1.0;
	for(int i=0;i<ncells;i++)
	{
		for(int c=0;c<nqpts;c++)
		{
			for(int b=0;b<nqpts;b++)
			{
				for(int a=0;a<nqpts;a++)
				{
					ind = c*nqpts*nqpts + b*nqpts + a;
					
					x = m_mesh->Cvec[i].eltnodes_x[ind];
					y = m_mesh->Cvec[i].eltnodes_y[ind];
					z = m_mesh->Cvec[i].eltnodes_z[ind];
					expval =  (x-x0)*(x-x0)/a0/a0;
					//expval += (y-y0)*(y-y0)/b0/b0;
					
					uvar->realarray[(offs+ind)*NCVAR+MASS] = u0*exp(-5.0*expval); 
					uvar->realarray[(offs+ind)*NCVAR+XMOM] = u0*exp(-5.0*expval); 
					uvar->realarray[(offs+ind)*NCVAR+YMOM] = u0*exp(-5.0*expval); 
					uvar->realarray[(offs+ind)*NCVAR+ZMOM] = u0*exp(-5.0*expval);
					uvar->realarray[(offs+ind)*NCVAR+ENRG] = (G_AIR*m_P0/m_rho0)
										*u0*exp(-5.0*expval);

				}
			}
		}

		offs += nlocaldof;
	}

}
//====================================================================
void Euler_manager::computeRHS()
{
	int nlocaldof,ncells;

	nlocaldof = (m_porder+1)*(m_porder+1)*(m_porder+1);
	ncells = m_mesh->getncells();

	m_computeconvtermAX();
	m_computefaceterms();
	m_computeboundaryterms();

	/*for(int i=0;i<ncells;i++)
	{
		std::cout<<"cell "<<i<<"\n**********\n";
		for(int k=0;k<nlocaldof;k++)
		{
			std::cout<<"dofnum:"<<k<<"\n\n";
			for(int c=0;c<NCVAR;c++)
			{
				std::cout<<"comp:"<<c<<"\trhs:"
					<<rhsvar->realarray[NCVAR*(i*nlocaldof+k)+c]<<"\n";
			}
		}
	}*/
	
}
//====================================================================
void Euler_manager::computeMassMat()
{
	int ncells,nlocaldof,offset;
	double *localmassmat;

	ncells=m_mesh->getncells();
	nlocaldof = (m_porder+1)*(m_porder+1)*(m_porder+1);

	localmassmat = new double[nlocaldof*NCVAR];

	for(int i=0;i<ncells;i++)
	{
		offset=i*nlocaldof;
		m_mesh->Cvec[i].computelocalmassmat(NCVAR,localmassmat);

		for(int k=0;k<nlocaldof;k++)
		{
			for(int c=0;c<NCVAR;c++) //c for component
			{
				massmat->realarray[NCVAR*(offset+k)+c] 
					= localmassmat[NCVAR*k+c];
			}
		}
	}

	for(int i=0;i<ncells;i++)
	{
		//std::cout<<"cell :"<<i<<"\n***********\n";

		for(int k=0;k<nlocaldof;k++)
		{
			for(int c=0;c<NCVAR;c++)
			{
				//std::cout<<"k:"<<k<<"\tmassmat:"<<massmat->realarray[NCVAR*(i*nlocaldof+k)+c]<<"\n";
			}
		}
	}
}
//====================================================================
void Euler_manager::m_computeconvtermAX()
{
	int ncells,nlocaldof,offset;
	double *localsoln,*localrhs,*localflux;

	ncells = m_mesh->getncells();
	nlocaldof = (m_porder+1)*(m_porder+1)*(m_porder+1);

	localsoln  = new double[nlocaldof*NCVAR];
	localrhs   = new double[nlocaldof*NCVAR];
	localflux  = new double[nlocaldof*NCVAR*NDIM];

	//std::cout<<"computing conv term\n";
	for(int i=0;i<ncells;i++)
	{
		offset=i*nlocaldof;
		for(int k=0;k<nlocaldof;k++)
		{
			for(int c=0;c<NCVAR;c++)
			{
				localsoln[k*NCVAR+c] = uvar->realarray[NCVAR*(offset+k)+c];
			}
		}

		m_computeinviscidflux(nlocaldof,localsoln,localflux);

		m_mesh->Cvec[i].computelocalconvtermAX(NCVAR,localsoln,localflux,localrhs);
		
		for(int k=0;k<nlocaldof;k++)
		{
			for(int c=0;c<NCVAR;c++)
			{
				rhsvar->realarray[NCVAR*(offset+k)+c] = localrhs[NCVAR*k+c];
			}
		}
	
	}
	for(int i=0;i<ncells;i++)
	{
		//std::cout<<"cell :"<<i<<"\n***********\n";

		for(int k=0;k<nlocaldof;k++)
		{

			for(int c=0;c<NCVAR;c++)
			{
				//std::cout<<"k: "<<k<<"\tconvAX :"<<rhsvar->realarray[NCVAR*(i*nlocaldof+k)+c]<<"\n";
				//std::cout<<"k: "<<k<<"\tlocalflux :"<<localflux[NCVAR*NDIM*k+NDIM*c+0]<<"\t";
				//std::cout<<localflux[NCVAR*NDIM*k+NDIM*c+1]<<"\t"<<localflux[NCVAR*NDIM*k+NDIM*c+2]<<"\n";
			}
		}
	}
	//std::cout<<"computed conv term\n";
}
//====================================================================
void Euler_manager::m_computeinviscidflux(int n,double *localsoln,double *localflux)
{
	double rho,u,v,w,p,e;

	double Fx[NCVAR],Fy[NCVAR],Fz[NCVAR];

	for(int i=0;i<n;i++)
	{
		rho = localsoln[i*NCVAR+MASS];
		u   = localsoln[i*NCVAR+XMOM];
		v   = localsoln[i*NCVAR+YMOM];
		w   = localsoln[i*NCVAR+ZMOM];
		p   = localsoln[i*NCVAR+ENRG];


		/*Fx[MASS]=m_rho0*u;	Fy[MASS]=m_rho0*v;	Fz[MASS]=m_rho0*w;
		Fx[XMOM]=p/m_rho0;	Fy[XMOM]=0.0;		Fz[XMOM]=0.0;
		Fx[YMOM]=0.0;		Fy[YMOM]=p/m_rho0;	Fz[YMOM]=0.0;
		Fx[ZMOM]=0.0;		Fy[ZMOM]=0.0;		Fz[ZMOM]=p/m_rho0;
		Fx[ENRG]=G_AIR*m_P0*u;	Fy[ENRG]=G_AIR*m_P0*v;	Fz[ENRG]=G_AIR*m_P0*w;*/
		
		/*Fx[MASS]=2.0*rho+u;	Fy[MASS]=2.0*rho;		Fz[MASS]=0.0;
		Fx[XMOM]=3.0*rho+2.0*u;	Fy[XMOM]=2.0*u;			Fz[XMOM]=0.0;
		Fx[YMOM]=0.0;		Fy[YMOM]=1.0*v;			Fz[YMOM]=0.0;
		Fx[ZMOM]=0.0;		Fy[ZMOM]=3.0*w;			Fz[ZMOM]=0.0;
		Fx[ENRG]=3.0*p;		Fy[ENRG]=3.0*p;			Fz[ENRG]=0.0;*/
		
		/*Fx[MASS]=2.0*rho+u;	Fy[MASS]=0.0;			Fz[MASS]=0.0;
		Fx[XMOM]=3.0*rho+2.0*u;	Fy[XMOM]=0.0;			Fz[XMOM]=0.0;
		Fx[YMOM]=0.0;		Fy[YMOM]=0.0;			Fz[YMOM]=0.0;
		Fx[ZMOM]=0.0;		Fy[ZMOM]=0.0;			Fz[ZMOM]=0.0;
		Fx[ENRG]=3.0*p;		Fy[ENRG]=0.0;			Fz[ENRG]=0.0;*/
		
		Fx[MASS]=0.5*rho*rho;	Fy[MASS]=0.0;		        Fz[MASS]=0.0;
		Fx[XMOM]=0.0;		Fy[XMOM]=0.0;			Fz[XMOM]=0.0;
		Fx[YMOM]=0.0;		Fy[YMOM]=0.0;			Fz[YMOM]=0.0;
		Fx[ZMOM]=0.0;		Fy[ZMOM]=0.0;			Fz[ZMOM]=0.0;
		Fx[ENRG]=0.0;		Fy[ENRG]=0.0;			Fz[ENRG]=0.0;

		for(int j=0;j<NCVAR;j++)
		{
			localflux[i*NCVAR*NDIM+j*NDIM+0]=Fx[j];
			localflux[i*NCVAR*NDIM+j*NDIM+1]=Fy[j];
			localflux[i*NCVAR*NDIM+j*NDIM+2]=Fz[j];
			
			//std::cout<<"j:"<<j<<"\t"<<"flux:"<<localflux[i*NCVAR*NDIM+j*NDIM+0]<<"\n";
		}

	}
}
//====================================================================
double Euler_manager::m_computepressure(double u[NCVAR])
{
	double rho,p,rhoe,one_by_rho;

	rho = u[MASS];
	one_by_rho = 1.0/rho;
	rhoe = u[ENRG];
	
	p = (G_AIR-1)*(rhoe - 0.5*one_by_rho*
			(u[XMOM]*u[XMOM]+u[YMOM]*u[YMOM]+u[ZMOM]*u[ZMOM]));

	return(p);

}
//====================================================================
void Euler_manager::m_computefaceterms()
{
	int nintfaces,findex,offset_l,offset_r;
	int lcell,rcell;
	double *localsoln_l,*localsoln_r;
	double *delta_ql,*delta_qr,*limiter;
	double *localfaceflux_l;
	double *localfaceflux_r;
	double *localf,*localwavespeed,*localcourantnumber;
	int nlocaldof,nfacedof;
	int *faceindices_l,*faceindices_lprv;
	int *faceindices_r,*faceindices_rprv;
	double *facepoints_l;
	int sindexl,sindexr,sindexl_prv,sindexr_prv;
	double uL[NCVAR],uR[NCVAR];
	double eps,mod_c;
	double r;

	eps = 1e-6;

	nlocaldof = (m_porder+1)*(m_porder+1)*(m_porder+1);
	nfacedof  = (m_porder+1)*(m_porder+1);

	faceindices_l = new int[nfacedof];
	faceindices_r = new int[nfacedof];
	
	faceindices_lprv = new int[nfacedof];
	faceindices_rprv = new int[nfacedof];
	facepoints_l  = new double[NDIM*nfacedof];

	localsoln_l = new double[NCVAR*nfacedof];
	localsoln_r = new double[NCVAR*nfacedof];
	
	delta_ql    = new double[NCVAR*nfacedof];
	delta_qr    = new double[NCVAR*nfacedof];
	limiter     = new double[NCVAR*nfacedof];
	localf      = new double[NCVAR*nfacedof];

	localfaceflux_l = new double[NCVAR*nfacedof];
	localfaceflux_r = new double[NCVAR*nfacedof];

	localwavespeed = new double[nfacedof];
	localcourantnumber = new double[nfacedof];

	nintfaces = m_mesh->interiorfaceids.size();

	for(int i=0;i<nintfaces;i++)
	{
		
		findex = m_mesh->interiorfaceids[i];

		lcell = m_mesh->Fvec[findex].lcellid;
		rcell = m_mesh->Fvec[findex].rcellid;

		offset_l = lcell*nlocaldof;
		offset_r = rcell*nlocaldof;

		m_mesh->Cvec[lcell].returnlocalfacepoints(findex,
				faceindices_l,facepoints_l);

		m_mesh->Cvec[rcell].getindicesfromxyz(nfacedof,
				facepoints_l,faceindices_r);


		m_mesh->Cvec[lcell].getpreviouslayerofpoints(faceindices_l,faceindices_lprv);
		m_mesh->Cvec[rcell].getpreviouslayerofpoints(faceindices_r,faceindices_rprv);

		for(int k=0;k<nfacedof;k++)
		{
			sindexl = faceindices_l[k];
			sindexr = faceindices_r[k];
			
			sindexl_prv = faceindices_lprv[k];
			sindexr_prv = faceindices_rprv[k];

			//std::cout<<"sindexl:"<<sindexl<<"\t"<<"sindexl_prv:"<<sindexl_prv<<"\n";
			//std::cout<<"sindexr:"<<sindexr<<"\t"<<"sindexl_prv:"<<sindexr_prv<<"\n";

			for(int c=0;c<NCVAR;c++)
			{
				localsoln_l[k*NCVAR+c] = uvar->realarray[NCVAR*(offset_l+sindexl)+c];
				localsoln_r[k*NCVAR+c] = uvar->realarray[NCVAR*(offset_r+sindexr)+c];
				
				
				//localsoln_l[k*NCVAR+c] = 0.5*(uvar->realarray[NCVAR*(offset_l+sindexl)+c] +
				//				uvar->realarray[NCVAR*(offset_l+sindexl_prv)+c]);

				//localsoln_r[k*NCVAR+c] = 0.5*(uvar->realarray[NCVAR*(offset_r+sindexr)+c] +
				//				uvar->realarray[NCVAR*(offset_r+sindexr_prv)+c]);


				//localsoln_l[k*NCVAR+c] = uL[c];

				delta_ql[k*NCVAR+c] = localsoln_r[k*NCVAR+c]
					- uvar->realarray[NCVAR*(offset_l+sindexl_prv)+c];

				delta_qr[k*NCVAR+c] = uvar->realarray[NCVAR*(offset_r+sindexr_prv)+c] - 
					localsoln_l[k*NCVAR+c];
				//localsoln_r[k*NCVAR+c] = uR[c];

				//std::cout<<"localsoln_l and uL:"<<localsoln_l[k*NCVAR+c]<<"\t"<<uL[c]<<"\n";
			}

		}
		

		m_getwavespeed(nfacedof,m_mesh->Fvec[findex].normal,
				localsoln_l,localsoln_r,localwavespeed,localcourantnumber);
		
		/*std::cout<<"face "<<i<<"\n";
		for(int k=0;k<nfacedof;k++)
		{
			std::cout<<"localcourantnumber of "<<k<<":"<<localcourantnumber[k]<<"\n";
			//std::cout<<"localwavespeed of "<<k<<":"<<localwavespeed[k]<<"\n";
		}*/


		for(int k=0;k<nfacedof;k++)
		{
			for(int c=0;c<NCVAR;c++)
			{
				if(delta_qr[k*NCVAR+c] <= eps)
				{
					if(delta_ql[k*NCVAR+c] <= eps)
					{
						limiter[k*NCVAR+c] = 0.0;
					}
					else
					{
						limiter[k*NCVAR+c] = 2.0; 
					}
				}
				else
				{
					r = delta_ql[k*NCVAR+c]/delta_qr[k*NCVAR+c];
					limiter[k*NCVAR+c] = m_vanleerlimiter(r,
							localcourantnumber[k]);

					if(r<0)
					{
						std::cout<<"limiter:"<<limiter[k*NCVAR+c]<<
							"\t"<<delta_ql[k*NCVAR+c]<<"\t"<<delta_qr[k*NCVAR+c]<<"\n";
					}

				}


			}

			
		}

		m_computetvdflux(nfacedof,m_mesh->Fvec[findex].normal,limiter,
				localsoln_l,localsoln_r,localf);

		//m_computeLaxFriedrichflux(nfacedof,m_mesh->Fvec[findex].normal,localwavespeed,
		//		localsoln_l,localsoln_r,localf);
		  //std::cout<<"localwavespeed:\n";
		  /*for(int j=0;j<nfacedof;j++)
		  std::cout<<localwavespeed[j]<<"\t";
		  std::cout<<"\n";*/

		//m_mesh->Cvec[lcell].getfaceflux(findex,m_mesh->Fvec[findex].normal,
		//		m_mesh->Fvec[findex].area,localsoln_l,localvel_l,localfaceflux_l);

		//m_mesh->Cvec[rcell].getfaceflux(findex,m_mesh->Fvec[findex].normal,
		//	m_mesh->Fvec[findex].area,localsoln_r,localvel_r,localfaceflux_r);

		m_mesh->Cvec[lcell].getinterfaceflux(NCVAR,findex,
				m_mesh->Fvec[findex].area,faceindices_l,localf,localfaceflux_l);

		m_mesh->Cvec[rcell].getinterfaceflux(NCVAR,findex,
				m_mesh->Fvec[findex].area,faceindices_r,localf,localfaceflux_r);

		/*std::cout<<"localfaceflux_l:\n";
		  for(int j=0;j<nfacedof;j++)
		  std::cout<<localfaceflux_l[j]<<"\t";
		  std::cout<<"\n";

		  std::cout<<"localfaceflux_r:\n";
		  for(int j=0;j<nfacedof;j++)
		  std::cout<<localfaceflux_r[j]<<"\t";
		  std::cout<<"\n";*/

		for(int k=0;k<nfacedof;k++)
		{
			sindexl = faceindices_l[k];
			sindexr = faceindices_r[k];

			for(int c=0;c<NCVAR;c++)
			{
				rhsvar->realarray[NCVAR*(offset_l+sindexl)+c] -= localfaceflux_l[NCVAR*k+c];
				rhsvar->realarray[NCVAR*(offset_r+sindexr)+c] += localfaceflux_r[NCVAR*k+c];
			}

			//std::cout<<"k:"<<k<<"\n";
		}

	}

	//std::cout<<"finished face flux calc\n";
	delete(faceindices_l);
	delete(faceindices_r);
	
	delete(faceindices_lprv);
	delete(faceindices_rprv);
	delete(facepoints_l);

	delete(localsoln_l);
	delete(localsoln_r);
	
	delete(delta_ql);
	delete(delta_qr);
	delete(limiter);
	delete(localf);  

	delete(localfaceflux_l);
	delete(localfaceflux_r);

	delete(localwavespeed);
	delete(localcourantnumber);

}
//================================================================================
double Euler_manager::m_vanleerlimiter(double r,double mod_c)
{
	double val;

	if(r <= 0 || r>1)
	{
		val = 2.0;
	}
	else
	{
		val= 1.0 - (1.0-mod_c)*2.0*r/(1.0 + r); 
	}


	return(val);

}
//================================================================================
void Euler_manager::m_getavgsolnvalue(int cellid,double u[NCVAR])
{
	int nlocaldof,offset;

	nlocaldof = (m_porder+1)*(m_porder+1)*(m_porder+1);

	for(int c=0;c<NCVAR;c++)
	{
		u[c]=0.0;
	}

	offset = cellid*nlocaldof;

	for(int i=0;i<nlocaldof;i++)
	{
		for(int c=0;c<NCVAR;c++)
		{
			u[c] += uvar->realarray[NCVAR*(offset+i)+c];
		}
	}

	for(int c=0;c<NCVAR;c++)
	{
		u[c]=u[c]/double(nlocaldof);
	}

	//std::cout<<"cellid:"<<cellid<<"\t"<<u[0]<<"\t"<<u[1]<<"\t"<<u[2]<<"\t"<<u[3]<<"\t"<<u[4]<<"\n";
}
//================================================================================
void Euler_manager::m_getwavespeed(int n,double normal[3],double *localsoln_l,
		double *localsoln_r,double *localwavespeed,double *localcourantnumber)
{
	double rho,p,rho_u,rho_v,rho_w,rho_e;
	double al,ar,vl_n,vr_n;

	for(int i=0;i<n;i++)
	{
		/*rho   = localsoln_l[i*NCVAR+0];
		rho_u = localsoln_l[i*NCVAR+1];
		rho_v = localsoln_l[i*NCVAR+2];
		rho_w = localsoln_l[i*NCVAR+3];
		rho_e = localsoln_l[i*NCVAR+4];
		p     = (G_AIR-1)*(rho_e
				-0.5*(rho_u*rho_u+rho_v*rho_v+rho_w*rho_w)/rho);
		al = sqrt(G_AIR*p/rho);

		vl_n   = (rho_u*normal[0]+rho_v*normal[1]+rho_w*normal[2])/rho;

		rho   = localsoln_r[i*NCVAR+0];
		rho_u = localsoln_r[i*NCVAR+1];
		rho_v = localsoln_r[i*NCVAR+2];
		rho_w = localsoln_r[i*NCVAR+3];
		rho_e = localsoln_r[i*NCVAR+4];
		p     = (G_AIR-1)*(rho_e
				-0.5*(rho_u*rho_u+rho_v*rho_v+rho_w*rho_w)/rho);
		ar = sqrt(G_AIR*p/rho);

		vr_n   = (rho_u*normal[0]+rho_v*normal[1]+rho_w*normal[2])/rho;

		localwavespeed[i]=0.5*((al+ar)+(fabs(vl_n)+fabs(vr_n)));*/

		//std::cout<<"localsoln_l:"<<localsoln_l[i*NCVAR+0]<<"\t"<<localsoln_l[i*NCVAR+1]<<
		//	"\t"<<localsoln_l[i*NCVAR+2]<<"\t"<<localsoln_l[i*NCVAR+3]<<"\t"<<localsoln_l[i*NCVAR+4]<<"\n";
		//std::cout<<"localwavespeed:"<<(al+ar)<<"\t"<<fabs(vl_n)+fabs(vr_n)<<"\n";
		//localwavespeed[i]=sqrt(G_AIR*m_P0/m_rho0);
		localwavespeed[i]=0.5*(localsoln_l[i*NCVAR]+localsoln_r[i*NCVAR]);
		localcourantnumber[i] = localwavespeed[i]*m_dt/0.005;

	}

}
//================================================================================
void Euler_manager::m_computeLaxFriedrichflux(int n,double normal[3],double *wavespeeds
		,double *localsoln_l,double *localsoln_r,double *localfaceflux)
{
	double *localfaceflux_l;
	double *localfaceflux_r;

	double lflux_n,rflux_n;

	localfaceflux_l = new double[n*NCVAR*NDIM];
	localfaceflux_r = new double[n*NCVAR*NDIM];

	m_computeinviscidflux(n,localsoln_l,localfaceflux_l);
	m_computeinviscidflux(n,localsoln_r,localfaceflux_r);
	
//	std::cout<<"normal:"<<normal[0]<<"\t"<<normal[1]<<"\t"<<normal[2]<<"\n";

	for(int i=0;i<n;i++)
	{
//		std::cout<<"wavespeed of "<<i<<"\t"<<wavespeeds[i]<<"\n";

		for(int j=0;j<NCVAR;j++)
		{
			lflux_n =   localfaceflux_l[i*NCVAR*NDIM+j*NDIM+0]*normal[0] +
				    localfaceflux_l[i*NCVAR*NDIM+j*NDIM+1]*normal[1] +
				    localfaceflux_l[i*NCVAR*NDIM+j*NDIM+2]*normal[2];

			rflux_n =   localfaceflux_r[i*NCVAR*NDIM+j*NDIM+0]*normal[0] +
				    localfaceflux_r[i*NCVAR*NDIM+j*NDIM+1]*normal[1] +
				    localfaceflux_r[i*NCVAR*NDIM+j*NDIM+2]*normal[2];

			localfaceflux[i*NCVAR+j]  = 0.5*(lflux_n + rflux_n);
			localfaceflux[i*NCVAR+j] -= 0.5*wavespeeds[i]*(localsoln_r[i*NCVAR+j] - localsoln_l[i*NCVAR+j]);


		}
	}

}
//================================================================================
void Euler_manager::m_computetvdflux(int n,double normal[3],double *limiter
		,double *localsoln_l,double *localsoln_r,double *localfaceflux)
{
	double *localfaceflux_l;
	double *localfaceflux_r;

	double lflux_n,rflux_n;

	localfaceflux_l = new double[n*NCVAR*NDIM];
	localfaceflux_r = new double[n*NCVAR*NDIM];

	m_computeinviscidflux(n,localsoln_l,localfaceflux_l);
	m_computeinviscidflux(n,localsoln_r,localfaceflux_r);

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<NCVAR;j++)
		{
			lflux_n =   localfaceflux_l[i*NCVAR*NDIM+j*NDIM+0]*normal[0] +
				    localfaceflux_l[i*NCVAR*NDIM+j*NDIM+1]*normal[1] +
				    localfaceflux_l[i*NCVAR*NDIM+j*NDIM+2]*normal[2];

			rflux_n =   localfaceflux_r[i*NCVAR*NDIM+j*NDIM+0]*normal[0] +
				    localfaceflux_r[i*NCVAR*NDIM+j*NDIM+1]*normal[1] +
				    localfaceflux_r[i*NCVAR*NDIM+j*NDIM+2]*normal[2];

			localfaceflux[i*NCVAR+j]  = 0.5*(lflux_n + rflux_n);
			localfaceflux[i*NCVAR+j] -= 0.5*limiter[i*NCVAR+j]*(rflux_n-lflux_n);

		}
	}

}
//================================================================================
void Euler_manager::m_computeboundaryterms()
{
	double *localsoln,*localvel,*localfaceflux,*localf;
	std::vector<int> indices;
	int nboundaries,findex,cellid,nfacedof,nlocaldof;
	int offset;
	double sign;
	double u[NCVAR];
	double p;
	int *faceindices;

	nfacedof = (m_porder+1)*(m_porder+1);
	nlocaldof = nfacedof*(m_porder+1);
	nboundaries = m_mesh->boundaryvec.size();

	localsoln     = new double[NCVAR*nfacedof];
	localf        = new double[NCVAR*nfacedof];
	localfaceflux = new double[NCVAR*nfacedof];
	localvel      = new double[NDIM*nfacedof];
	faceindices   = new int[nfacedof];

	//for time being	
	m_mesh->bconditions.resize(nboundaries);

	m_mesh->bconditions[0]=BC_SYMM;
	m_mesh->bconditions[1]=BC_SYMM;
	//m_mesh->bconditions[2]=BC_INFLOW_SUBSONIC;	
	m_mesh->bconditions[2]=BC_DIRICHLET;	
	m_mesh->bconditions[3]=BC_DIRICHLET;	
	m_mesh->bconditions[4]=BC_SYMM;	
	m_mesh->bconditions[5]=BC_SYMM;	

	for(int i=0;i<nboundaries;i++)
	{
		//std::cout<<"boundary num:"<<i<<"\n";
		if(m_mesh->bconditions[i] == BC_WALL)
		{
			for(unsigned int j=0;j<m_mesh->boundaryvec[i].Fvecids.size();j++)
			{
				findex=m_mesh->boundaryvec[i].Fvecids[j];

				//std::cout<<"findex:"<<findex<<"\n";

				cellid=(m_mesh->Fvec[findex].lcellid==-1)?
					m_mesh->Fvec[findex].rcellid:m_mesh->Fvec[findex].lcellid;

				//std::cout<<"cellid:"<<cellid<<"\n";
				
				sign=(m_mesh->Fvec[findex].lcellid==cellid)?1.0:-1.0;
			
				indices=m_mesh->Cvec[cellid].getDofsforDirichlet(findex);

				for(int k=0;k<nfacedof;k++)
				{
					faceindices[k]=indices[k];
					//std::cout<<"indices at "<<k<<":"<<faceindices[k]<<"\n";
				}
				
				offset = cellid*nlocaldof;
				for(unsigned int k=0;k<nfacedof;k++)
				{
					for(int c=0;c<NCVAR;c++)
					{
						u[c] = uvar->realarray[NCVAR*(offset+indices[k])+c];
					}
					p = m_computepressure(u);

					localf[k*NCVAR+MASS]=0.0;
					localf[k*NCVAR+XMOM]=p*m_mesh->Fvec[findex].normal[0];
					localf[k*NCVAR+YMOM]=p*m_mesh->Fvec[findex].normal[1];
					localf[k*NCVAR+ZMOM]=p*m_mesh->Fvec[findex].normal[2];
					localf[k*NCVAR+ENRG]=0.0;

					//std::cout<<"localf at "<<k<<":"<<localf[0]<<"\t"<<localf[1]<<
					//	"\t"<<localf[2]<<"\t"<<localf[3]<<"\t"<<localf[4]<<"\n";
				}

				m_mesh->Cvec[cellid].getinterfaceflux(NCVAR,findex,m_mesh->Fvec[findex].area,faceindices,
						localf,localfaceflux);
				for(int k=0;k<nfacedof;k++)
				{
					for(int c=0;c<NCVAR;c++)
					{
						rhsvar->realarray[NCVAR*(offset+indices[k])+c] -= 
							sign*localfaceflux[NCVAR*k+c];
					}
				}
			}

		}
		else
			if(m_mesh->bconditions[i] == BC_DIRICHLET)
			{
				for(unsigned int j=0;j<m_mesh->boundaryvec[i].Fvecids.size();j++)
				{
					findex=m_mesh->boundaryvec[i].Fvecids[j];

					cellid=(m_mesh->Fvec[findex].lcellid==-1)?
						m_mesh->Fvec[findex].rcellid:m_mesh->Fvec[findex].lcellid;

					indices=m_mesh->Cvec[cellid].getDofsforDirichlet(findex);

					offset = cellid*nlocaldof;
					for(unsigned int k=0;k<indices.size();k++)
					{
						for(int c=0;c<NCVAR;c++)
						{
							rhsvar->realarray[NCVAR*(offset+indices[k])+c]=0.0;
						}
					}
				}

			}
			else
				if(m_mesh->bconditions[i] == BC_FLUX)
				{

				}
				else
					if(m_mesh->bconditions[i] == BC_SYMM)
					{
						//zero flux
						//no need to do anything
					}
					else
					{
						std::cout<<"Boundary condition code "<<
							m_mesh->bconditions[i]<<" not implemented yet\n";
						exit(0);
					}	
	}

	//std::cout<<"finished calling boundary terms\n";
}
//================================================================================
