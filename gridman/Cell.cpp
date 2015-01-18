#include"Cell.h"

//==================================================================================
void Cell::setbasis(int order)
{
	int index;
	double x,y,z;

	m_porder = order;
	m_nqpts = m_porder+1;

	m_N_x = new Lagrange[m_porder+1];
	m_N_y = new Lagrange[m_porder+1];
	m_N_z = new Lagrange[m_porder+1];

	m_qpts = new double[m_porder+1];
	m_qwts = new double[m_porder+1];

	eltnodes_x = new double[m_nqpts*m_nqpts*m_nqpts];
	eltnodes_y = new double[m_nqpts*m_nqpts*m_nqpts];
	eltnodes_z = new double[m_nqpts*m_nqpts*m_nqpts];

	Legendre leg_obj(m_nqpts-1);
	leg_obj.gauss_lobatto_points(m_qpts,m_qwts,m_nqpts);

	for(int k=0;k<m_nqpts;k++)
	{
		for(int j=0;j<m_nqpts;j++)
		{
			for(int i=0;i<m_nqpts;i++)
			{
				index = k*m_nqpts*m_nqpts+j*m_nqpts+i;

				m_master_to_physical(m_qpts[i],m_qpts[j],m_qpts[k],
						x,y,z);

				eltnodes_x[index] = x;
				eltnodes_y[index] = y;
				eltnodes_z[index] = z;
			}
		}
	}

	for(int i=0;i<m_nqpts;i++)
	{
		m_N_x[i].assign_lagrange_points(m_qpts,i,m_porder+1);
		m_N_y[i].assign_lagrange_points(m_qpts,i,m_porder+1);
		m_N_z[i].assign_lagrange_points(m_qpts,i,m_porder+1);
	}


}
//==================================================================================
void Cell::m_master_to_physical(double a,double b,double c,double &x,double &y,double &z)
{
	int linearbasisnum=2;
	double Na,Nb,Nc;
	int index;

	x=0.0; y=0.0; z=0.0;

	for(int k=0;k<linearbasisnum;k++)
	{
		Nc = m_master_linearbasis(c,k);

		for(int j=0;j<linearbasisnum;j++)
		{
			Nb = m_master_linearbasis(b,j);

			for(int i=0;i<linearbasisnum;i++)
			{
				index = linearbasisnum*(k*linearbasisnum + j) + i;
				Na = m_master_linearbasis(a,i);

				x += coord[index][0]*Na*Nb*Nc;
				y += coord[index][1]*Na*Nb*Nc;
				z += coord[index][2]*Na*Nb*Nc;

			}
		}
	}
}
//==================================================================================
double Cell::m_Jacobian(double Jac[NDIM][NDIM],double a,double b,double c)
{
	double det;

	for(int i=0;i<NDIM;i++)
	{
		for(int j=0;j<NDIM;j++)
		{
			Jac[i][j]=m_deriv(i,j,a,b,c);	
		}
	}

	det =   Jac[0][0]*(Jac[1][1]*Jac[2][2] - Jac[1][2]*Jac[2][1]);
	det += -Jac[0][1]*(Jac[1][0]*Jac[2][2] - Jac[1][2]*Jac[2][0]);
	det +=  Jac[0][2]*(Jac[1][0]*Jac[2][1] - Jac[1][1]*Jac[2][0]);

	return(fabs(det));


}
//==================================================================================
double Cell::m_deriv(int xyznum,int abcnum,double a,double b,double c)
{
	int linearbasisnum=2;
	int index;
	double Na,Nb,Nc,val;

	val = 0.0;
	for(int k=0;k<linearbasisnum;k++)
	{
		Nc = (abcnum==2)?m_master_linearbasisderiv(k):m_master_linearbasis(c,k);

		for(int j=0;j<linearbasisnum;j++)
		{
			Nb = (abcnum==1)?m_master_linearbasisderiv(j):m_master_linearbasis(b,j);

			for(int i=0;i<linearbasisnum;i++)
			{
				Na = (abcnum==0)?m_master_linearbasisderiv(i)
					:m_master_linearbasis(b,i);
			
				index = linearbasisnum*(k*linearbasisnum+j)+i;

				val += coord[index][xyznum]*Na*Nb*Nc;
			}
		}
	}

	return(val);
}
//==================================================================================
double Cell::m_dphysical_dmaster(int dir,double a,double b,double c)
{
	double sum,half;

	half = 0.5;

	if(dir == 0)
	{
		sum  = coord[0][dir] * half*(1-b) * half*(1-c) * (-half);
		sum += coord[1][dir] * half*(1-b) * half*(1-c) *  (half);
		sum += coord[2][dir] * half*(1+b) * half*(1-c) * (-half);
		sum += coord[3][dir] * half*(1+b) * half*(1-c) *  (half);
		
		sum  = coord[4][dir] * half*(1-b) * half*(1+c) * (-half);
		sum += coord[5][dir] * half*(1-b) * half*(1+c) *  (half);
		sum += coord[6][dir] * half*(1+b) * half*(1+c) * (-half);
		sum += coord[7][dir] * half*(1+b) * half*(1+c) *  (half);
	}
	else
	if(dir == 1)
	{
		sum  = coord[0][dir] * half*(1-a) * half*(1-c) * (-half);
		sum += coord[1][dir] * half*(1+a) * half*(1-c) * (-half);
		sum += coord[2][dir] * half*(1-a) * half*(1-c) *  (half);
		sum += coord[3][dir] * half*(1+a) * half*(1-c) *  (half);
		
		sum  = coord[4][dir] * half*(1-a) * half*(1+c) * (-half);
		sum += coord[5][dir] * half*(1+a) * half*(1+c) * (-half);
		sum += coord[6][dir] * half*(1-a) * half*(1+c) *  (half);
		sum += coord[7][dir] * half*(1+a) * half*(1+c) *  (half);
		
	}
	else
	if(dir == 2)
	{
		sum  = coord[0][dir] * half*(1-a) * half*(1-b) * (-half);
		sum += coord[1][dir] * half*(1+a) * half*(1-b) * (-half);
		sum += coord[2][dir] * half*(1-a) * half*(1+b) * (-half);
		sum += coord[3][dir] * half*(1+a) * half*(1+b) * (-half);
		
		sum  = coord[4][dir] * half*(1-a) * half*(1-b) * (half);
		sum += coord[5][dir] * half*(1+a) * half*(1-b) * (half);
		sum += coord[6][dir] * half*(1-a) * half*(1+b) * (half);
		sum += coord[7][dir] * half*(1+a) * half*(1+b) * (half);
	}
	else
	{
		std::cerr<<"dir "<<dir<<" not valid "<<
			"in function m_dphysical_dmaster\n";
	}

	return(sum);

}
//==================================================================================
void Cell::printcellvtu()
{
	char underscore;
	int offset;
	int intval,index;
	double x,y,z;
	double val;

	int nodenum,cellnum;

	nodenum = m_nqpts*m_nqpts*m_nqpts;
	cellnum = m_porder*m_porder*m_porder;

	std::ofstream outfile;
	outfile.open("meshcell.vtu");

	offset=0;
	underscore='_';

	outfile<<"<?xml version=\"1.0\"?>\n";	
	outfile<<"<VTKFile type=\"UnstructuredGrid\" ";
	outfile<<"version=\"0.1\" byte_order=\"LittleEndian\">\n";	
	outfile<<"<UnstructuredGrid>\n";
	outfile<<"<Piece NumberOfPoints=\""<<nodenum<<"\" NumberOfCells=\""
		<<cellnum<<"\">\n";

	outfile<<"<PointData>\n";
	outfile<<"<DataArray type=\"Float64\" Name=\"circumrad\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</PointData>\n";
	offset = offset + sizeof(int) + sizeof(double)*nodenum;
	
	outfile<<"<Points>\n";
	outfile<<"<DataArray name=\"Position\" type=\"Float64\" ";
	outfile<<"NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</Points>\n";

	offset = offset + sizeof(int) + NDIM*sizeof(double)*nodenum;

	outfile<<"<Cells>\n";
	outfile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""
		<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + NCPTS*sizeof(int)*cellnum;
	outfile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""
		<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + sizeof(int)*cellnum;	
	outfile<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""
		<<offset<<"\"/>\n";

	outfile<<"</Cells>\n";
	outfile<<"</Piece>\n";
	outfile<<"</UnstructuredGrid>\n";
	outfile<<"<AppendedData encoding=\"raw\">\n";

	outfile.write(&underscore,sizeof(char));

	intval=sizeof(double)*nodenum;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<nodenum;i++)
	{
		val=i;
		outfile.write((char *)&val,sizeof(double));	
	}

	//writing node data======================================
	intval=NDIM*sizeof(double)*nodenum;
	outfile.write((char *)&intval,sizeof(int));

	std::cout<<"Points\n";
	for(int k=0;k<m_nqpts;k++)
	{
		for(int j=0;j<m_nqpts;j++)
		{ 
			for(int i=0;i<m_nqpts;i++)
			{
				index = k*m_nqpts*m_nqpts+j*m_nqpts+i;

				x = eltnodes_x[index]; 
				y = eltnodes_y[index];
				z = eltnodes_z[index];

				outfile.write((char *)&x,sizeof(double));
				outfile.write((char *)&y,sizeof(double));
				outfile.write((char *)&z,sizeof(double));

				std::cout<<eltnodes_x[index]<<"\t"<<eltnodes_y[index]<<"\t"
					<<eltnodes_z[index]<<"\n";
			}
		}
	}
	//=======================================================

	std::cout<<"Cells\n";
	//writing cell nodes=====================================
	intval=NCPTS*sizeof(int)*cellnum;
	outfile.write((char *)&intval,sizeof(int));

	for(int k=0;k<m_nqpts-1;k++)
	{
		for(int j=0;j<m_nqpts-1;j++)
		{
			for(int i=0;i<m_nqpts-1;i++)
			{
				index=k*m_nqpts*m_nqpts+j*m_nqpts+i;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=k*m_nqpts*m_nqpts+j*m_nqpts+i+1;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=k*m_nqpts*m_nqpts+(j+1)*m_nqpts+i+1;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=k*m_nqpts*m_nqpts+(j+1)*m_nqpts+i;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				
				
				
				index=(k+1)*m_nqpts*m_nqpts + j*m_nqpts + i;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=(k+1)*m_nqpts*m_nqpts+j*m_nqpts+i+1;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=(k+1)*m_nqpts*m_nqpts+(j+1)*m_nqpts+i+1;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";
				
				index=(k+1)*m_nqpts*m_nqpts+(j+1)*m_nqpts+i;
				outfile.write((char *)&index,sizeof(int));
				std::cout<<eltnodes_x[index]<<"\t"
					<<eltnodes_y[index]<<"\t"<<eltnodes_z[index]<<"\n";


				std::cout<<"\n\n";

				

			}
		}
	}
	//=======================================================

	//writing cell offsets===================================
	intval=sizeof(int)*cellnum;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<cellnum;i++)
	{
		intval=(i+1)*NCPTS;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	//writing cell types=====================================
	intval=sizeof(int)*cellnum;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<cellnum;i++)
	{
		intval=12;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	outfile.close();

}
//==================================================================================
void Cell::computelocalmassmat(int ncomp,double *localmassmat)
{
	int index;
	double basisval;
	double Jac_det;
	double Jac[NDIM][NDIM];
	double value;

	Jac_det = m_Jacobian(Jac,1.0,1.0,1.0);

	for(int c=0;c<m_nqpts;c++)
	{
		for(int b=0;b<m_nqpts;b++)
		{
			for(int a=0;a<m_nqpts;a++)
			{
				index = c*m_nqpts*m_nqpts+b*m_nqpts+a;

				Jac_det = m_Jacobian(Jac,m_qpts[a],m_qpts[b],
						m_qpts[c]);
				basisval = 1.0;
				value = Jac_det*
					basisval*m_qwts[a]*m_qwts[b]*m_qwts[c];

				for(int i=0;i<ncomp;i++)
				{
					localmassmat[index*ncomp+i]=value;
				}


			}
		}
	}
}
//==================================================================================
void Cell::computelocalconvtermAX(int ncomp,double *localsoln,
		double *localflux,double *localrhs)
{
	double deriv_val_a,deriv_val_b,deriv_val_c;
	int l,m,n; //test function indices, n is the quickest
	int index;
	int nlocaldof;
	double Jac[NDIM][NDIM];
	double Jacinv[NDIM][NDIM];
	double Jac_det;
	double dadx,dady,dadz;
	double dbdx,dbdy,dbdz;
	double dcdx,dcdy,dcdz;
	double value;

	nlocaldof = m_nqpts*m_nqpts*m_nqpts;

	for(int i=0;i<nlocaldof;i++) //test function loop
	{
		//test function indices
		n = i/(m_nqpts*m_nqpts);
		m = (i%(m_nqpts*m_nqpts))/m_nqpts;
		l = (i%(m_nqpts*m_nqpts))%m_nqpts; //quickest index

		for(int k=0;k<ncomp;k++)
		{
			localrhs[i*ncomp+k] = 0.0;
		}
		for(int j=0;j<m_nqpts;j++)  //trial function loop
		{
			index = n*m_nqpts*m_nqpts + m*m_nqpts + j;
			Jac_det = m_Jacobian(Jac,m_qpts[j],m_qpts[m],m_qpts[n]);
			invertmat(Jac,Jacinv);
			
			deriv_val_a = m_N_x[l].find_der_value_at_x(m_qpts[j]);
			dadx = Jacinv[0][0];
			dady = Jacinv[0][1];
			dadz = Jacinv[0][2];
			
			for(int k=0;k<ncomp;k++)
			{
				value = 0.0;
				value += localflux[ncomp*NDIM*index+k*NDIM+XDIR]*dadx;
				value += localflux[ncomp*NDIM*index+k*NDIM+YDIR]*dady;
				value += localflux[ncomp*NDIM*index+k*NDIM+ZDIR]*dadz;

				localrhs[i*ncomp+k] += value*
					deriv_val_a*Jac_det*
					m_qwts[j]*m_qwts[m]*m_qwts[n];
			}

		}
		
		for(int j=0;j<m_nqpts;j++)  //trial function loop
		{
			index = n*m_nqpts*m_nqpts + j*m_nqpts + l;
			Jac_det = m_Jacobian(Jac,m_qpts[l],m_qpts[j],m_qpts[n]);
			invertmat(Jac,Jacinv);
			
			deriv_val_b = m_N_y[m].find_der_value_at_x(m_qpts[j]);
			dbdx = Jacinv[1][0];
			dbdy = Jacinv[1][1];
			dbdz = Jacinv[1][2];
			
			for(int k=0;k<ncomp;k++)
			{
				value  = 0.0;
				value += localflux[ncomp*NDIM*index+k*NDIM+XDIR]*dbdx;
				value += localflux[ncomp*NDIM*index+k*NDIM+YDIR]*dbdy;
				value += localflux[ncomp*NDIM*index+k*NDIM+ZDIR]*dbdz;
				
				localrhs[i*ncomp+k] += value*
					deriv_val_b*Jac_det*
					m_qwts[l]*m_qwts[j]*m_qwts[n];
			}
		}
		
		for(int j=0;j<m_nqpts;j++)  //trial function loop
		{
			index = j*m_nqpts*m_nqpts + m*m_nqpts + l;
			Jac_det = m_Jacobian(Jac,m_qpts[l],m_qpts[m],m_qpts[j]);
			invertmat(Jac,Jacinv);
			
			deriv_val_c = m_N_z[n].find_der_value_at_x(m_qpts[j]);
			dcdx = Jacinv[2][0];
			dcdy = Jacinv[2][1];
			dcdz = Jacinv[2][2];
			
			for(int k=0;k<ncomp;k++)
			{
				value  = 0.0;
				value += localflux[ncomp*NDIM*index+k*NDIM+XDIR]*dcdx;
				value += localflux[ncomp*NDIM*index+k*NDIM+YDIR]*dcdy;
				value += localflux[ncomp*NDIM*index+k*NDIM+ZDIR]*dcdz;
				
				localrhs[i*ncomp+k] += value*
					deriv_val_c*Jac_det*
					m_qwts[l]*m_qwts[m]*m_qwts[j];
			}
		}
	}	

}
//==================================================================================
void Cell::returnlocalfacepoints(int findex,int *indices,double *points)
{
	int index;
	double x,y,z;
	int k,f;

	//find which face we are talking about
	for(f=0;f<NCFACES;f++)
	{
		if(faceids[f]==findex)
		{
			break;
		}
	}

	switch(f)
	{
		case 0:
			k=0;
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = b*m_nqpts + a;
					indices[k]=index;
					m_master_to_physical(m_qpts[a],m_qpts[b],-1.0,
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;

		case 1:
			k=0;
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = (m_nqpts-1)*m_nqpts*m_nqpts+ b*m_nqpts + a;
					indices[k]=index;
					m_master_to_physical(m_qpts[a],m_qpts[b],1.0,
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;
		case 2:
			k=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts;
					indices[k]=index;
					m_master_to_physical(-1.0,m_qpts[b],m_qpts[c],
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;
		case 3:
			k=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts + (m_nqpts-1);
					indices[k]=index;
					m_master_to_physical(1.0,m_qpts[b],m_qpts[c],
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;
		case 4:
			k=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + a;
					indices[k]=index;
					m_master_to_physical(m_qpts[a],-1.0,m_qpts[c],
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;
		case 5:
			k=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + (m_nqpts-1)*m_nqpts + a;
					indices[k]=index;
					m_master_to_physical(m_qpts[a],1.0,m_qpts[c],
							x,y,z);
					points[NDIM*k+0] = x; 
					points[NDIM*k+1] = y;
					points[NDIM*k+2] = z;

					k++;
				}
			}
			break;

		default:
			std::cout<<"Face not found\n";
			exit(0);
			break;
	}

}
//==================================================================================
void Cell::getpreviouslayerofpoints(int *faceindices,int *prvslayer)
{
	int a,b,c;
	int sum_a,sum_b,sum_c;
	int max_sum;

	//check first three points to know which face
	sum_a=0; sum_b=0; sum_c=0;
	max_sum=m_nqpts*m_nqpts*(m_nqpts-1);
	
	for(int i=0;i<m_nqpts*m_nqpts;i++)
	{
		m_getcomponentindices(faceindices[i],a,b,c);

		sum_a += a;
		sum_b += b;
		sum_c += c;
	}

	//std::cout<<"sum_a:"<<sum_a<<"\tsum_b:"<<sum_b<<"\tsum_c:"<<sum_c<<"\n";

	if(sum_a == 0)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]+1;
		}
	}

	if(sum_a == max_sum)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]-1;
		}
	}

	if(sum_b == 0)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]+m_nqpts;
		}
	}
	if(sum_b == max_sum)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]-m_nqpts;
		}

	}
	if(sum_c == 0)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]+m_nqpts*m_nqpts;
		}

	}
	if(sum_c == max_sum)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			prvslayer[i]=faceindices[i]-m_nqpts*m_nqpts;
		}
	}
	
}
//==================================================================================
void Cell::getfaceflux(int findex,double n[NDIM],double area,double *localsoln,
		double *localvel,double *localfaceflux)
{
	double Jac_det;
	int index,f;

	//this is an approximation. Only exact if physical element is
	//also cartesian.
	Jac_det = 0.25*area; 

	//find which face we are talking about
	for(f=0;f<NCFACES;f++)
	{
		if(faceids[f]==findex)
		{
			break;
		}
	}
	 
	//initialize all faceflux terms to 0
	for(int i=0;i<m_nqpts*m_nqpts*m_nqpts;i++)
	{
		localfaceflux[i]=0.0;
	}

	switch(f)
	{
		case 0:
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = b*m_nqpts + a;
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[a]*m_qwts[b];	
				}
			}
			break;

		case 1:
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = (m_nqpts-1)*m_nqpts*m_nqpts+ b*m_nqpts + a;
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[a]*m_qwts[b];	
				}
			}
			break;
		case 2:
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts;
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[b]*m_qwts[c];	
				}
			}
			break;
		case 3:
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts + (m_nqpts-1);
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[b]*m_qwts[c];	
				}
			}
			break;
		case 4:
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + a;
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[c]*m_qwts[a];	
				}
			}
			break;
		case 5:
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + (m_nqpts-1)*m_nqpts + a;
					localfaceflux[index] = (localvel[NDIM*index+0]*n[0]+
						localvel[NDIM*index+1]*n[1]+
					localvel[NDIM*index+2]*n[2])*Jac_det*localsoln[index]
						*m_qwts[c]*m_qwts[a];	
				}
			}
			break;

		default:
			std::cout<<"Face not found\n";
			exit(0);
			break;


	}

}
//==================================================================================
void Cell::m_getcomponentindices(int index,int &a,int &b,int &c)
{
	c = index/(m_nqpts*m_nqpts);
	b = (index%(m_nqpts*m_nqpts))/m_nqpts;
	a = (index%(m_nqpts*m_nqpts))%m_nqpts;
}
//==================================================================================
void Cell::getinterfaceflux(int ncomp,int findex,double area,int *faceindices,
		double *localf,double *localfaceflux)
{
	double Jac_det;
	int f;
	int a,b,c;

	//this is an approximation. Only exact if physical element is
	//also cartesian.
	Jac_det = 0.25*area; 

	//find which face we are talking about
	for(f=0;f<NCFACES;f++)
	{
		if(faceids[f]==findex)
		{
			break;
		}
	}

	if(f==0 || f==1)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			m_getcomponentindices(faceindices[i],a,b,c);
			
			for(int k=0;k<ncomp;k++)
			{
				localfaceflux[i*ncomp+k] = Jac_det
					*localf[i*ncomp+k]
				*m_qwts[a]*m_qwts[b];	
			}
		}
	}
	if(f==2 || f==3)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			m_getcomponentindices(faceindices[i],a,b,c);

			for(int k=0;k<ncomp;k++)
			{
				localfaceflux[i*ncomp+k] = Jac_det*localf[i*ncomp+k]
				*m_qwts[b]*m_qwts[c];	
			}
		}

	}
	if(f==4 || f==5)
	{
		for(int i=0;i<m_nqpts*m_nqpts;i++)
		{
			m_getcomponentindices(faceindices[i],a,b,c);
			
			for(int k=0;k<ncomp;k++)
			{
				localfaceflux[i*ncomp+k] = Jac_det*localf[i*ncomp+k]
				*m_qwts[c]*m_qwts[a];
			}	
		}

	}

}
//==================================================================================
std::vector<int> Cell::getDofsforDirichlet(int findex)
{
	int f,pos,index;
	std::vector<int> indices;

	indices.resize(m_nqpts*m_nqpts);

	for(f=0;f<NCFACES;f++)
	{
		if(faceids[f]==findex)
		{
			break;
		}
	}

	switch(f)
	{
		case 0:
			pos=0;
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = b*m_nqpts + a;
					indices[pos++]=index;
				}
			}
			break;

		case 1:
			pos=0;
			for(int b=0;b<m_nqpts;b++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = (m_nqpts-1)*m_nqpts*m_nqpts+ b*m_nqpts + a;
					indices[pos++]=index;
				}
			}
			break;
		case 2:
			pos=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts;
					indices[pos++]=index;
				}
			}
			break;
		case 3:
			pos=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int b=0;b<m_nqpts;b++)
				{
					index = c*m_nqpts*m_nqpts + b*m_nqpts + (m_nqpts-1);
					indices[pos++]=index;
				}
			}
			break;
		case 4:
			pos=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + a;
					indices[pos++]=index;
				}
			}
			break;
		case 5:
			pos=0;
			for(int c=0;c<m_nqpts;c++)
			{
				for(int a=0;a<m_nqpts;a++)
				{
					index = c*m_nqpts*m_nqpts + (m_nqpts-1)*m_nqpts + a;
					indices[pos++]=index;
				}
			}
			break;

		default:
			std::cout<<"Face not found\n";
			exit(0);
			break;


	}

	return(indices);
	 

}
//==================================================================================
void Cell::getindicesfromxyz(int np,double *coord,int *indices)
 {
	 double x,y,z;
	 double a,b,c;
	 bool lflag,mflag,nflag,found;
	 double eps;

	 int l,m,n;

	 eps = 1e-4;

	 for(int i=0;i<np;i++)
	 {
		x = coord[NDIM*i+0];
		y = coord[NDIM*i+1];
		z = coord[NDIM*i+2];

		m_physical_to_master(x,y,z,a,b,c);

		//std::cout<<"a,b,c:"<<a<<"\t"<<b<<"\t"<<c<<"\n";
		lflag = false; mflag = false; nflag = false;
		found = false;
		for(int q=0;q<m_nqpts;q++)
		{

			if(fabs(m_qpts[q]-a) < eps)
			{
				l = q;
				lflag = true;
			}
			if(fabs(m_qpts[q]-b) < eps)
			{
				m = q;
				mflag = true;
			}
			if(fabs(m_qpts[q]-c) < eps)
			{
				n = q;
				nflag = true;
			}

			found = lflag & mflag & nflag;
			if(found)
			{
				break;
			}

		}
		if(found)
		{
			//std::cout<<"lmn:"<<l<<"\t"<<m<<"\t"<<n<<"\n";
			indices[i] = n*m_nqpts*m_nqpts + m*m_nqpts + l;
		}
		else
		{
			//std::cout<<"face point not found\n";
			indices[i] = -1;
		}
	 }
 }	
//==================================================================================
void Cell::m_physical_to_master(double x,double y,double z,double &a,double &b,double &c)
{
	double mindist;
	double dist;
	double a_g,b_g,c_g;
	int lmin,mmin,nmin;
	int index;
	double det,jac[NDIM][NDIM],jacinv[NDIM][NDIM];
	double err[NDIM],x_g,y_g,z_g;
	double change[NDIM];
	double errnorm,tol;
	int MAXIT,it;

	MAXIT = 10;
	tol = 1e-4;
	errnorm = 3e8;

	mindist=3e8;
	lmin = 0; mmin=0; nmin=0;
	for(int n=0;n<m_nqpts;n++)
	{
		for(int m=0;m<m_nqpts;m++)
		{
			for(int l=0;l<m_nqpts;l++)
			{
				index = n*m_nqpts*m_nqpts+m*m_nqpts+l;

				dist = (x-eltnodes_x[index])*(x-eltnodes_x[index]) +
					(y-eltnodes_y[index])*(y-eltnodes_y[index]) +
					(z-eltnodes_z[index])*(z-eltnodes_z[index]);

				if(dist < mindist)
				{	
					mindist = dist;
					lmin = l; mmin=m; nmin=n;
				}

			}
		}
	}

	a_g = m_qpts[lmin];
	b_g = m_qpts[mmin];
	c_g = m_qpts[nmin];

	it = 0;
	while(it < MAXIT && errnorm > tol)
	{
		//find jacobiapn
		det = m_Jacobian(jac,a_g,b_g,c_g);
		invertmat(jac,jacinv);

		m_master_to_physical(a_g,b_g,c_g,x_g,y_g,z_g);

		err[0] = x - x_g;
		err[1] = y - y_g;
		err[2] = z - z_g;

		errnorm = err[0]*err[0]+err[1]*err[1]+err[2]*err[2];

		for(int i=0;i<NDIM;i++)
		{
			change[i]=0.0;
			for(int j=0;j<NDIM;j++)
			{
				change[i] += jacinv[i][j]*err[j];
			}
		}

		a_g += change[0];
		b_g += change[1];
		c_g += change[2];

		it++;
	}

	a=a_g; b=b_g; c=c_g;


}
//==================================================================================
void Cell::invertmat(double M[NDIM][NDIM],double Minv[NDIM][NDIM])
{
	double a1,a2,a3;
	double b1,b2,b3;
	double c1,c2,c3;
	double det,detinv;
	
	a1 = M[0][0]; a2=M[0][1]; a3=M[0][2];
	b1 = M[1][0]; b2=M[1][1]; b3=M[1][2];
	c1 = M[2][0]; c2=M[2][1]; c3=M[2][2];

	det = a1*(b2*c3-b3*c2) + a2*(b3*c1-b1*c3) + a3*(b1*c2-b2*c1);
	if(det != 0.0)
	{
		detinv = 1.0/det;

		Minv[0][0] = (b2*c3-b3*c2)*detinv;
		Minv[0][1] = (c2*a3-c3*a2)*detinv;
		Minv[0][2] = (a2*b3-a3*b2)*detinv;

		Minv[1][0] = (b3*c1-b1*c3)*detinv;
		Minv[1][1] = (c3*a1-c1*a3)*detinv;
		Minv[1][2] = (a3*b1-a1*b3)*detinv;

		Minv[2][0] = (b1*c2-b2*c1)*detinv;
		Minv[2][1] = (c1*a2-c2*a1)*detinv;
		Minv[2][2] = (a1*b2-a2*b1)*detinv;
	}
	else
	{
		std::cerr<<"Determinant is 0 while inverting 3 x 3 matrix\n";
	}

}
//==================================================================================
