#include<iostream>
#include<cmath>
#include<vector>

#ifndef FIELDVAR_H
#define FIELDVAR_H
class Fieldvar
{
	private:
		std::string m_variablename;
		int m_dim;
		int m_num;

	public:
		double *realarray;
		void setvarprops(std::string name,int num,int dim=1)
		{
			m_variablename = name;
			m_num=num;
			m_dim=dim;

			realarray = new double[m_num*m_dim];
		}

		std::string getname()
		{
			return(m_variablename);
		}
		int getsize()
		{
			return(m_num);
		}
		int getdim()
		{
			return(m_dim);
		}

		/*Fieldvar operator=(const Fieldvar &rhs)
		{
			this->m_dim = rhs.m_dim;
			this->m_num = rhs.m_num;

			for(int i=0;i<rhs.m_num*rhs.m_dim;i++)
			{
				this->realarray[i]=rhs.realarray[i];
			}

			return(*this);
		}*/

		void copy(Fieldvar *rhs)
		{
			for(int i=0;i<m_num*m_dim;i++)
			{
				this->realarray[i]=rhs->realarray[i];
			}
		}

};
#endif
