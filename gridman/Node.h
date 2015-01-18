#include<iostream>
#include<cmath>
#ifndef NODE_CLASS
#define NODE_CLASS
class Node
{
	private:
		double m_x,m_y,m_z;
	public:
		void setcoord(double x,double y,double z)
		{
			m_x=x; m_y=y; m_z=z;
		}
		void getcoord(double &x,double &y,double &z)
		{
			x=m_x; y=m_y; z=m_z;
		}
};
#endif
