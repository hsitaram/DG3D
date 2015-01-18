#include<fstream>
#include<iostream>
#include<vector>

int main()
{
	int l,m,n;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double dx,dy,dz,x,y,z;
	int count,origin_node,node;

	xmin=0.0; xmax=1.0;
	ymin=0.0; ymax=1.0;
	zmin=0.0; zmax=0.1;

	l=10; m=10; n=2;

	dx = (xmax-xmin)/double(l-1);
	dy = (ymax-ymin)/double(m-1);
	dz = (zmax-zmin)/double(n-1);

	std::ofstream outfile("boxmesh.dat");

	count=0;
	outfile<<"Nodes\t"<<l*m*n<<"\n";
	for(int k=0;k<n;k++)
	{
		z=k*dz;

		for(int j=0;j<m;j++)
		{
			y=j*dy;

			for(int i=0;i<l;i++)
			{
				x=i*dx;
				outfile<<x<<"\t"<<y<<"\t"<<z<<"\n";
			}
		}
	}

	outfile<<"Cells\t"<<(l-1)*(m-1)*(n-1)<<"\n";

	for(int k=0;k<n-1;k++)
	{
		for(int j=0;j<m-1;j++)
		{
			for(int i=0;i<l-1;i++)
			{
				origin_node = k*m*l +j*l + i;

				outfile<<origin_node<<"\t"<<origin_node+1<<"\t";
				outfile<<origin_node+l<<"\t"<<origin_node+l+1<<"\t";
				
				outfile<<origin_node+m*l<<"\t"<<origin_node+m*l+1<<"\t";
				outfile<<origin_node+l+m*l<<"\t"<<origin_node+l+m*l+1<<"\n";
			}
		}
	}

	outfile<<"Boundaries\t"<<6<<"\n";

	outfile<<"Bottom\t"<<l*m<<"\n";
	for(int j=0;j<m;j++)
	{
		for(int i=0;i<l;i++)
		{
			node=0*m*l+j*l+i;
			outfile<<node<<"\n";
		}
	}
	
	outfile<<"Top\t"<<l*m<<"\n";
	for(int j=0;j<m;j++)
	{
		for(int i=0;i<l;i++)
		{
			node=(n-1)*m*l+j*l+i;
			outfile<<node<<"\n";
		}
	}
	
	outfile<<"Left\t"<<n*m<<"\n";
	for(int k=0;k<n;k++)
	{
		for(int j=0;j<m;j++)
		{
			node=k*m*l+j*l+0;
			outfile<<node<<"\n";
		}
	}
	
	outfile<<"Right\t"<<n*m<<"\n";
	for(int k=0;k<n;k++)
	{
		for(int j=0;j<m;j++)
		{
			node=k*m*l+j*l+l-1;
			outfile<<node<<"\n";
		}
	}
	
	outfile<<"Front\t"<<n*l<<"\n";
	for(int k=0;k<n;k++)
	{
		for(int i=0;i<l;i++)
		{
			node=k*m*l+0*l+i;
			outfile<<node<<"\n";
		}
	}
	
	outfile<<"Back\t"<<n*l<<"\n";
	for(int k=0;k<n;k++)
	{
		for(int i=0;i<l;i++)
		{
			node=k*m*l+(m-1)*l+i;
			outfile<<node<<"\n";
		}
	}



	outfile.close();

	return(0);
}
