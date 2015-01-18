#include"Euler_timestepper.h"

//===================================================================
void Euler_timestepper::init(std::string inpfilename)
{
	std::string meshfilename;
	std::string temp;
	std::ifstream infile;
	std::vector<Fieldvar *> fieldvarptrs;
	
	infile.open(inpfilename.c_str());

	infile>>temp>>meshfilename;
	infile>>temp>>m_porder;
	infile>>temp>>m_maxsteps;
	infile>>temp>>m_dt;
	infile>>temp>>m_fintime;
	infile>>temp>>m_timealgo;

	fieldvarptrs.resize(0);

	m_eulerptr = new Euler_manager;
	m_outputptr   = new Euler_outputs;

	m_eulerptr->init(meshfilename,m_porder,m_dt);

	m_uvar      = m_eulerptr->uvar;
	m_uvar_prvs = m_eulerptr->uvar_prvs;
	m_rhsvar    = m_eulerptr->rhsvar;
	m_massmat   = m_eulerptr->massmat;

	fieldvarptrs.push_back(m_uvar);

	m_outputptr->meshptr = m_eulerptr->getmeshptr(); 	
	m_outputptr->porder  = m_porder;
	m_outputptr->fieldvarptrs = fieldvarptrs;
	m_outputptr->writevtufile("file_00000.vtu");
	m_outputptr->writeonedfile("initial.dat");

	infile.close();


}
//===============================================================
void Euler_timestepper::ExplicitEE()
{
	std::string numberstring;
	char num[10];
	std::string fname;
	int output_it;
	int it=0;
	double t=0.0;
	int ndim,np;
	output_it=m_maxsteps/10;

	np = m_uvar->getsize();
	ndim = m_uvar->getdim();

	if(output_it == 0)
	{
		output_it=1;
	}

	while(it < m_maxsteps && t < m_fintime)
	{
		m_eulerptr->computeRHS();
		m_eulerptr->computeMassMat();
		
		for(int i=0;i<np;i++)
		{
			//watch out for rhs sign
			for(int k=0;k<ndim;k++)
			{
				m_uvar->realarray[i*ndim+k] 
					+= m_dt*
					m_rhsvar->realarray[i*ndim+k]
					/m_massmat->realarray[i*ndim+k];
			}
		}

		/*for(int i=0;i<m_uvar->getsize();i++)
		{
			std::cout<<"i:"<<i<<"\n";
			for(int k=0;k<m_uvar->getdim();k++)
			{
				std::cout<<"k:"<<k<<"\t"
					<<m_uvar->realarray[i*m_uvar->getdim()+k]<<"\n";
			}
		}*/

		it++;
		t=t+m_dt;

		if(it%output_it == 0)
		{
			sprintf(num,"%5.5d",it);
			numberstring = num;

			fname = "file_"+numberstring+".vtu";
			m_outputptr->writevtufile(fname);
		
		}

		std::cout<<"it:"<<it<<"\ttime:"<<t<<"\n";
	}
	sprintf(num,"%5.5d",it+1);
	numberstring = num;
	fname = "file_"+numberstring+".vtu";
	m_outputptr->writevtufile(fname);
	m_outputptr->writeonedfile("final.dat");

}
//==============================================================
void Euler_timestepper::ExplicitRK2()
{
	std::string numberstring;
	char num[10];
	std::string fname;
	int output_it;
	int it=0;
	double t=0.0;
	output_it=m_maxsteps/10;
	int ndim,np;

	if(output_it == 0)
	{
		output_it=1;
	}

	np = m_uvar->getsize();
	ndim = m_uvar->getdim();


	while(it < m_maxsteps && t < m_fintime)
	{
		m_uvar_prvs->copy(m_uvar);

		//stage 1
		m_eulerptr->computeRHS();
		m_eulerptr->computeMassMat();
		for(int i=0;i<np;i++)
		{
			for(int k=0;k<ndim;k++)
			{
				m_uvar->realarray[i*ndim+k] 
					+= m_dt*
					m_rhsvar->realarray[i*ndim+k]
					/m_massmat->realarray[i*ndim+k];
			}


		}
		std::cout<<"finished stage 1\n";

		//stage 2
		m_eulerptr->computeRHS();
		m_eulerptr->computeMassMat();
		for(int i=0;i<np;i++)
		{
			for(int k=0;k<ndim;k++)
			{
				m_uvar->realarray[i*ndim+k] = 0.5*(m_uvar->realarray[i*ndim+k]+
						m_uvar_prvs->realarray[i*ndim+k]);

				m_uvar->realarray[i*m_uvar->getdim()+k] 
					+= 0.5*m_dt*
					m_rhsvar->realarray[i*m_uvar->getdim()+k]
					/m_massmat->realarray[i*m_uvar->getdim()+k];


			}
		}
		std::cout<<"finished stage 2\n";

		/*for(int i=0;i<m_uvar->getsize();i++)
		  {
		  std::cout<<"i:"<<i<<"\n";
		  for(int k=0;k<m_uvar->getdim();k++)
		  {
		  std::cout<<"k:"<<k<<"\t"
		  <<m_uvar->realarray[i*m_uvar->getdim()+k]<<"\n";
		  }
		  }*/

		it++;
		t=t+m_dt;

		if(it%output_it == 0)
		{
			sprintf(num,"%5.5d",it);
			numberstring = num;

			fname = "file_"+numberstring+".vtu";
			m_outputptr->writevtufile(fname);

		}

		std::cout<<"it:"<<it<<"\ttime:"<<t<<"\n";
	}
	sprintf(num,"%5.5d",it+1);
	numberstring = num;
	fname = "file_"+numberstring+".vtu";
	m_outputptr->writevtufile(fname);
	m_outputptr->writeonedfile("final.dat");

}
//==============================================================
