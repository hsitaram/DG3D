#include"Face.h"
#ifndef BOUNDARY_CLASS
#define BOUNDARY_CLASS
class Boundary
{
	private:
		std::string m_boundaryname;
		bool procboundaryflag;

	public:
		std::vector<int> Fvecids;
		void setboundaryname(std::string bname)
		{
			//todo : set procboundaryflag only if 
			//bname is of the
			//form proc_a_b
			m_boundaryname=bname;
			procboundaryflag=false;
		        Fvecids.resize(0);	
		}
		std::string getboundaryname()
		{
			return(m_boundaryname);
		}
		bool isprocboundary()
		{
			return procboundaryflag;
		}

};
#endif
