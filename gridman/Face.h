#include<iostream>
#include<cmath>
#define NFPTS 4
#define NDIM 3
#ifndef FACE_CLASS
#define FACE_CLASS
class Face
{
	private:
		bool m_normalfound;
		double m_nodecoords[NFPTS][NDIM];
		void m_findareaandnormal();
		double m_findareaoftriangle(double x1,double y1,double z1,
				double x2,double y2,double z2,double x3,
				double y3,double z3);

	public:
		int lcellid,rcellid;
		int nodes[NFPTS];
		double normal[NDIM];
		double centroid[NDIM];
		double area;
		
		static void st_intsort(int arr[],int n);
		void setnodes(int arr[NFPTS],double coords[NFPTS][NDIM]);
		bool isitsame(int n[NFPTS]);

};
#endif
