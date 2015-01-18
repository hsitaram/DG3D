#include"Euler_manager.h"
#include"Euler_outputs.h"

#ifndef EULER_TIMESTEPPER_H
#define EULER_TIMESTEPPER_H
class Euler_timestepper
{
	private:
		Fieldvar *m_uvar;
		Fieldvar *m_uvar_prvs;
		Fieldvar *m_rhsvar;
		Fieldvar *m_massmat;
		Fieldvar *m_velvar;

		int m_porder;
		int m_maxsteps;
		double m_dt;
		double m_fintime;
		std::string m_timealgo;
		Euler_manager *m_eulerptr;
		Euler_outputs *m_outputptr;

	public:

		void init(std::string inpfilename);
		void ExplicitEE();
		void ExplicitRK2();
		void Explicit_timestepping()
		{
			if(m_timealgo == "EE")
			{
				ExplicitEE();
			}
			else
			if(m_timealgo == "RK2")
			{
				ExplicitRK2();
			}
			else
			{
				std::cout<<"Time stepping algorithm "<<m_timealgo
					<<" not implemented yet\n";
			}
		}


};
#endif
