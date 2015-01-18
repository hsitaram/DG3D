#include"Node.h"
#include"Cell.h"
#include"Face.h"
#include<string>
#include<vector>
#include<fstream>
#include"Boundary.h"
#define BC_DIRICHLET 0
#define BC_FLUX 1
#define BC_SYMM 2
#define BC_WALL 3
#define BC_INFLOW_SUBSONIC 4
#define BC_OUTFLOW_SUPERSONIC 5
#ifndef GRIDMAN_CLASS
#define GRIDMAN_CLASS
class Gridman
{

	private:
		std::vector<int> m_pointsperbc;
		std::vector<int *> m_bcpoints;
		int m_nnodes,m_ncells,m_nfaces;
		int m_nboundaries;

		std::vector<std::string> m_bcnames;
		void m_findfaces();
		bool m_isfaceadded(int fnodes[NFPTS],int nfaces,int &loc);
		int m_addface(int fnodes[NFPTS],int &top,int cellid);
		void m_findcellvolumes();
		double m_findpyramidvolume(double base_area,double base_norma[NDIM],
				double side[NDIM]);
		void m_findboundaryfaces();
		void m_orientfacenormals();

		std::vector< std::vector<int> > m_nodeidfaceidsmap;

	public:
	
		Cell *Cvec;
		Face *Fvec;
		Node *Nvec;
		std::vector<Boundary> boundaryvec;
		std::vector<int> interiorfaceids;
		std::vector<int> boundaryfaceids;
		std::vector<int> bconditions;


		void readgridfile(std::string gridfilename);
		void writelogfile();
		static void binarysearch(int arr[],int val,int start,int end,
				int size,int &pos,bool &found);
		void printmeshvtu();

		int getncells(){return m_ncells;}
		int getnfaces(){return m_nfaces;}

};
#endif
