#include"Euler_outputs.h"

//======================================================================================
void Euler_outputs::writevtufile(std::string fname)
{
	char underscore;
	int offset,offs;
	int intval;
	double x,y,z;
	int numnodes,numcells;
	int nlocaldof,nlocaldofcells;
	int nqpts;
	int nvals;

	nqpts = porder+1;

	std::ofstream outfile;
	outfile.open(fname.c_str());

	offset=0;
	underscore='_';

	nlocaldof = nqpts*nqpts*nqpts;
	nlocaldofcells = porder*porder*porder;
	numnodes = meshptr->getncells()*nlocaldof;

	numcells = meshptr->getncells();

	outfile<<"<?xml version=\"1.0\"?>\n";	
	outfile<<"<VTKFile type=\"UnstructuredGrid\" ";
	outfile<<"version=\"0.1\" byte_order=\"LittleEndian\">\n";	
	outfile<<"<UnstructuredGrid>\n";
	outfile<<"<Piece NumberOfPoints=\""<<numnodes<<"\" NumberOfCells=\""
		<<numcells*nlocaldofcells<<"\">\n";

	outfile<<"<PointData>\n";
	for(unsigned int i=0;i<fieldvarptrs.size();i++)
	{
		outfile<<"<DataArray type=\"Float64\" Name=";
		outfile<<"\""<<fieldvarptrs[i]->getname()<<"\" ";
	        outfile<<"NumberOfComponents="<<"\""<<fieldvarptrs[i]->getdim()<<"\" ";	
		outfile<<"format=\"appended\" offset=\""<<offset<<"\"/>\n";
		offset = offset + sizeof(int) + 
			sizeof(double)*fieldvarptrs[i]->getdim()*numnodes;
	}
	outfile<<"</PointData>\n";
	
	outfile<<"<Points>\n";
	outfile<<"<DataArray name=\"Position\" type=\"Float64\" ";
	outfile<<"NumberOfComponents=\"3\" format=\"appended\" offset=\""
		<<offset<<"\"/>\n";
	outfile<<"</Points>\n";
	offset = offset + sizeof(int) + NDIM*sizeof(double)*numnodes;

	outfile<<"<Cells>\n";
	outfile<<"<DataArray type=\"Int32\" Name=\"connectivity\" ";
	outfile<<"format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + sizeof(int)*NCPTS*nlocaldofcells*numcells;
	outfile<<"<DataArray type=\"Int32\" Name=\"offsets\" "; 
	outfile<<"format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + sizeof(int)*nlocaldofcells*numcells;	
	outfile<<"<DataArray type=\"Int32\" Name=\"types\" "; 
	outfile<<"format=\"appended\" offset=\""<<offset<<"\"/>\n";

	outfile<<"</Cells>\n";
	outfile<<"</Piece>\n";
	outfile<<"</UnstructuredGrid>\n";
	outfile<<"<AppendedData encoding=\"raw\">\n";

	outfile.write(&underscore,sizeof(char));

	std::cout<<"Finished writing text section\n";

	//write fieldvars
	for(unsigned int i=0;i<fieldvarptrs.size();i++)
	{
		nvals = fieldvarptrs[i]->getdim()*numnodes;
		intval=sizeof(double)*nvals;
		outfile.write((char *)&intval,sizeof(int));

		std::cout<<"writing "<<fieldvarptrs[i]->getname()<<"\n";

		for(int j=0;j<nvals;j++) 
		{
			outfile.write((char *)&fieldvarptrs[i]->realarray[j],
					sizeof(double));
		}

	}

	std::cout<<"finished writing fieldvars\n";

	//writing node data======================================
	intval=NDIM*sizeof(double)*numnodes;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<numcells;i++)
	{
		for(int j=0;j<nlocaldof;j++)
		{
			x = meshptr->Cvec[i].eltnodes_x[j];
			y = meshptr->Cvec[i].eltnodes_y[j];
			z = meshptr->Cvec[i].eltnodes_z[j];

			outfile.write((char *)&x,sizeof(double));
			outfile.write((char *)&y,sizeof(double));
			outfile.write((char *)&z,sizeof(double));
		}
	}
	//=======================================================

	std::cout<<"Finished writing node data\n";

	//writing cell nodes=====================================
	intval=sizeof(int)*NCPTS*nlocaldofcells*numcells;
	outfile.write((char *)&intval,sizeof(int));

	offs=0;
	for(int i=0;i<numcells;i++)
	{
		for(int c=0;c<porder;c++)
		{
		for(int b=0;b<porder;b++)
		{
		for(int a=0;a<porder;a++)
		{
			intval = offs+c*nqpts*nqpts+b*nqpts+a;
			outfile.write((char *)&intval,sizeof(int));

			intval = offs+c*nqpts*nqpts+b*nqpts+a+1;
			outfile.write((char *)&intval,sizeof(int));
			
			intval = offs+c*nqpts*nqpts+(b+1)*nqpts+a+1;
			outfile.write((char *)&intval,sizeof(int));
			
			intval = offs+c*nqpts*nqpts+(b+1)*nqpts+a;
			outfile.write((char *)&intval,sizeof(int));
			
			
			intval = offs+(c+1)*nqpts*nqpts+b*nqpts+a;
			outfile.write((char *)&intval,sizeof(int));

			intval = offs+(c+1)*nqpts*nqpts+b*nqpts+a+1;
			outfile.write((char *)&intval,sizeof(int));
			
			intval = offs+(c+1)*nqpts*nqpts+(b+1)*nqpts+a+1;
			outfile.write((char *)&intval,sizeof(int));
			
			intval = offs+(c+1)*nqpts*nqpts+(b+1)*nqpts+a;
			outfile.write((char *)&intval,sizeof(int));
		
		}
		}
		}

		offs += nlocaldof;
	}
	//=======================================================

	std::cout<<"Finished with connectivity\n";

	//writing cell offsets===================================
	intval=sizeof(int)*numcells*nlocaldofcells;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<numcells*nlocaldofcells;i++)
	{
		intval=(i+1)*NCPTS;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	//writing cell types=====================================
	intval=sizeof(int)*numcells*nlocaldofcells;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<numcells*nlocaldofcells;i++)
	{
		intval=12;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	std::cout<<"Done writing "<<fname<<"\n";

	outfile.close();
	
}
//======================================================================================
void Euler_outputs::writeonedfile(std::string fname)
{
	double x;
	int numcells,nlocaldof,index,nqpts;
	int ndim;
	std::ofstream outfile;
	outfile.open(fname.c_str());

	numcells = meshptr->getncells();
	nlocaldof = (porder+1)*(porder+1)*(porder+1);
	nqpts = porder+1;

	ndim=fieldvarptrs[0]->getdim();

	std::cout<<"writing one d file\n";
	
	for(int i=0;i<numcells;i++)
	{
		for(int j=0;j<(porder+1);j++)
		{
			index = nqpts*nqpts*(0.5*nqpts)+nqpts*(0.5*nqpts)+j;
			x = meshptr->Cvec[i].eltnodes_x[index];
			outfile<<x<<"\t";

			for(int c=0;c<ndim;c++)
			{
				outfile<<fieldvarptrs[0]->realarray[ndim*(i*nlocaldof+index)+c]<<"\t";
			}

			outfile<<"\n";
		}
	}

	std::cout<<"finished one d file\n";

	outfile.close();
}
//======================================================================================
