#include"Gridman.h"

int main()
{
	double Jac[NDIM][NDIM],det;
	double a,b,c;
	int indices[4];
	double coord[12];

	Gridman obj;
	obj.readgridfile("meshgen/boxmesh.dat");
	obj.writelogfile();

        obj.Cvec[0].setbasis(2);
        //obj.Cvec[1].setbasis(2);
	obj.Cvec[0].printcellvtu();


	obj.Cvec[0].m_physical_to_master(0.35,0.06,0.08,a,b,c);
	std::cout<<"a,b,c:"<<a<<"\t"<<b<<"\t"<<c<<"\n";


	coord[0] = 0.5; coord[1] = 0.0;	coord[2] = 0.0;
	coord[3] = 0.5; coord[4] = 0.05; coord[5] = 0.0;
	
	coord[6] = 0.5; coord[7] = 0.0;	coord[8] = 0.05;
	coord[9] = 0.5; coord[10] = 0.05; coord[11] = 0.15;
	
	//obj.Cvec[1].getindicesfromxyz(4,coord,indices);

	//std::cout<<"indices:"<<indices[0]<<"\t"<<indices[1]<<"\t"<<indices[2]
	//	<<"\t"<<indices[3]<<"\n";


	return(0);
}
