#include"Face.h"

//==========================================================
void Face::setnodes(int arr[NFPTS],double coords[NFPTS][NDIM])
{
	nodes[0]=arr[0]; nodes[1]=arr[1];
	nodes[2]=arr[2]; nodes[3]=arr[3];

	centroid[0]=0.0; centroid[1]=0.0; centroid[2]=0.0;
	for(int i=0;i<NFPTS;i++)
	{
		for(int k=0;k<NDIM;k++)
		{
			m_nodecoords[i][k]=coords[i][k];
			centroid[k] += coords[i][k];
		}
	}

	centroid[0] *= 0.25;
	centroid[1] *= 0.25;
	centroid[2] *= 0.25;

	lcellid=-1;
	rcellid=-1;
	m_normalfound=false;
	m_findareaandnormal();
}
//==========================================================
bool Face::isitsame(int n[NFPTS])
{
	bool flag=true;

	int ncopy[NFPTS],nodescopy[NFPTS];

	for(int i=0;i<NFPTS;i++)
	{
		ncopy[i]=n[i];
		nodescopy[i]=nodes[i];
	}

	st_intsort(ncopy,NFPTS);
	st_intsort(nodescopy,NFPTS);

	for(int i=0;i<NFPTS;i++)	
	{
		if(ncopy[i] != nodescopy[i])
		{
			flag=false;
			break;
		}
	}
	return(flag);
}
//==========================================================
void Face::m_findareaandnormal()
{
	double x1,y1,z1;	
	double x2,y2,z2;	
	double x3,y3,z3;	
	double x4,y4,z4;

	x1=m_nodecoords[0][0]; y1=m_nodecoords[0][1]; z1=m_nodecoords[0][2];
	x2=m_nodecoords[1][0]; y2=m_nodecoords[1][1]; z2=m_nodecoords[1][2];
	x3=m_nodecoords[3][0]; y3=m_nodecoords[3][1]; z3=m_nodecoords[3][2];
	x4=m_nodecoords[2][0]; y4=m_nodecoords[2][1]; z4=m_nodecoords[2][2];

	area=m_findareaoftriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)+
		m_findareaoftriangle(x1,y1,z1,x4,y4,z4,x3,y3,z3);

}
//==========================================================
void Face::st_intsort(int arr[],int n)
{
	int temp;
	//bubble sort
	for(int i=0;i<n-1;i++)
	{
		for(int j=0;j<n-i-1;j++)
		{
			if(arr[j] > arr[j+1])
			{
				temp     = arr[j];
				arr[j]   = arr[j+1];
				arr[j+1] = temp;	
			}
		}
	}

}
//===============================================================
double Face::m_findareaoftriangle(double x1,double y1,double z1,
				double x2,double y2,double z2,
				double x3,double y3,double z3)
{
	double a[NDIM],b[NDIM];
	double mod_a,mod_b,angle;
	double crossp[NDIM],mod_crossp;

	a[0]=(x2-x1); a[1]=(y2-y1); a[2]=(z2-z1);
	b[0]=(x3-x2); b[1]=(y3-y2); b[2]=(z3-z2);

	mod_a = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	mod_b = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);

	angle = acos((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/mod_a/mod_b);

	if(!m_normalfound)
	{
		crossp[0] = a[1]*b[2]-a[2]*b[1];
		crossp[1] = a[2]*b[0]-a[0]*b[2];
		crossp[2] = a[0]*b[1]-a[1]*b[0];

		mod_crossp = sqrt(crossp[0]*crossp[0]+crossp[1]*crossp[1]+crossp[2]*crossp[2]);

		for(int i=0;i<NDIM;i++)
		{
			normal[i]=crossp[i]/mod_crossp;
		}
		m_normalfound=true;
	}

	return(0.5*mod_a*mod_b*sin(angle));

}
//==========================================================
