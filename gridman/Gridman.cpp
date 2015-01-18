#include"Gridman.h"

//===============================================================
void Gridman::readgridfile(std::string gridfilename)
{
	std::string tempstr;
	int tempint;
	double x,y,z;
	double nbfaces_est,nbpoints;

	std::ifstream infile(gridfilename.c_str());

	infile>>tempstr;
	infile>>m_nnodes;

	Nvec = new Node[m_nnodes];

	for(int i=0;i<m_nnodes;i++)
	{
		infile>>x>>y>>z;
		Nvec[i].setcoord(x,y,z);	
	}

	infile>>tempstr;
	infile>>m_ncells;
	Cvec = new Cell[m_ncells];

	for(int i=0;i<m_ncells;i++)
	{
		for(int k=0;k<NCPTS;k++)
		{
			infile>>tempint;
			Cvec[i].nodeids[k] = tempint;
			Nvec[tempint].getcoord(Cvec[i].coord[k][0],Cvec[i].coord[k][1],
					Cvec[i].coord[k][2]);
		}

		Cvec[i].findcentroid();	
	}

	infile>>tempstr;
	infile>>m_nboundaries;

	m_pointsperbc.resize(m_nboundaries);
	m_bcpoints.resize(m_nboundaries);
	m_bcnames.resize(m_nboundaries);
	boundaryvec.resize(m_nboundaries);

	for(int i=0;i<m_nboundaries;i++)
	{
		infile>>m_bcnames[i]>>m_pointsperbc[i];
		boundaryvec[i].setboundaryname(m_bcnames[i]);

		m_bcpoints[i] = new int[m_pointsperbc[i]];

		for(int k=0;k<m_pointsperbc[i];k++)
		{
			infile>>tempint;
			m_bcpoints[i][k]=tempint;
		}	
	}

	nbpoints = 0;
	for(int i=0;i<m_nboundaries;i++)
	{
		nbpoints += m_pointsperbc[i];
	}

	//estimate number of faces
	nbfaces_est = nbpoints;

	m_nfaces = 0.5*(NCFACES*m_ncells - nbfaces_est) + nbfaces_est;

	m_nodeidfaceidsmap.resize(m_nnodes);
	for(int i=0;i<m_nnodes;i++)
	{
		m_nodeidfaceidsmap[i].resize(0);
	}

	Fvec = new Face[m_nfaces];
	infile.close();
	
	std::cout<<"finished mesh file reading\n";
	m_findfaces();
	std::cout<<"finished finding faces\n";
	m_findcellvolumes();
	std::cout<<"finished finding cell volumes\n";
	m_findboundaryfaces();
	std::cout<<"finished finding boundary faces\n";
	m_orientfacenormals();
	std::cout<<"oriented face normals\n";
	printmeshvtu();
	std::cout<<"finished printing mesh\n";
	


}
//===============================================================
void Gridman::m_findfaces()
{
	int fnodes[NFPTS];
	int top,loc;
	bool boundaryfaceflag;

	top=0;
	for(int i=0;i<m_ncells;i++)
	{
		//6 faces
		
		//face 0
		fnodes[0]=Cvec[i].nodeids[0];		
		fnodes[1]=Cvec[i].nodeids[1];		
		fnodes[2]=Cvec[i].nodeids[2];		
		fnodes[3]=Cvec[i].nodeids[3];
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[0]=loc;

		//face 1
		fnodes[0]=Cvec[i].nodeids[4];		
		fnodes[1]=Cvec[i].nodeids[5];		
		fnodes[2]=Cvec[i].nodeids[6];		
		fnodes[3]=Cvec[i].nodeids[7];	
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[1]=loc;
	

		//face 2
		fnodes[0]=Cvec[i].nodeids[0];		
		fnodes[1]=Cvec[i].nodeids[4];		
		fnodes[2]=Cvec[i].nodeids[2];		
		fnodes[3]=Cvec[i].nodeids[6];	
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[2]=loc;
	

		//face 3
		fnodes[0]=Cvec[i].nodeids[1];		
		fnodes[1]=Cvec[i].nodeids[5];		
		fnodes[2]=Cvec[i].nodeids[3];		
		fnodes[3]=Cvec[i].nodeids[7];	
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[3]=loc;

		
		//face 4
		fnodes[0]=Cvec[i].nodeids[0];		
		fnodes[1]=Cvec[i].nodeids[1];		
		fnodes[2]=Cvec[i].nodeids[4];		
		fnodes[3]=Cvec[i].nodeids[5];	
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[4]=loc;
	
	
		//face 5
		fnodes[0]=Cvec[i].nodeids[2];		
		fnodes[1]=Cvec[i].nodeids[3];		
		fnodes[2]=Cvec[i].nodeids[6];		
		fnodes[3]=Cvec[i].nodeids[7];	
		loc=m_addface(fnodes,top,i);
		Cvec[i].faceids[5]=loc;

			
	}	
	m_nfaces=top;

	boundaryfaceids.resize(0);
	interiorfaceids.resize(0);

	for(int i=0;i<m_nfaces;i++)
	{
		boundaryfaceflag=(Fvec[i].lcellid == -1 
				|| Fvec[i].rcellid == -1)?true:false;

		if(boundaryfaceflag)
		{
			boundaryfaceids.push_back(i);
		}
		else
		{
			interiorfaceids.push_back(i);
		}
	}
}
//===============================================================
bool Gridman::m_isfaceadded(int fnodes[NFPTS],int nfaces,int &loc)
{
	int fnodescopy[NFPTS];
	int min,fid;
	bool flag=false;

	for(int i=0;i<NFPTS;i++)
	{
		fnodescopy[i]=fnodes[i];
	}

	Face::st_intsort(fnodescopy,NFPTS);
	min = fnodescopy[0];

	loc=-1;
	for(unsigned int i=0;i<m_nodeidfaceidsmap[min].size();i++)
	{
		fid = m_nodeidfaceidsmap[min][i];

		if(Fvec[fid].isitsame(fnodes))
		{
			flag=true;
			loc=fid;
			break;
		}
	}
	if(!flag)
	{
		m_nodeidfaceidsmap[min].push_back(nfaces);
	}

	return(flag);
}
//===============================================================
int Gridman::m_addface(int fnodes[NFPTS],int &top,int cellid)
{
		int loc;
		double coords[NFPTS][NDIM];

		for(int j=0;j<NFPTS;j++)
		{
			Nvec[fnodes[j]].getcoord(coords[j][0],coords[j][1],
					coords[j][2]);
		}	

		if(!m_isfaceadded(fnodes,top,loc))
		{
			Fvec[top].setnodes(fnodes,coords);
			Fvec[top].lcellid = cellid;
			loc=top;
			top++;
		}
		else
		{
			Fvec[loc].rcellid = cellid;	
		}

		return(loc);

}
//===============================================================
void Gridman::m_findcellvolumes()
{
	double side[NDIM];
	double vol;

	for(int i=0;i<m_ncells;i++)
	{
		side[0]=Cvec[i].coord[4][0]-Cvec[i].coord[0][0];
		side[1]=Cvec[i].coord[4][1]-Cvec[i].coord[0][1];
		side[2]=Cvec[i].coord[4][2]-Cvec[i].coord[0][2];

		vol = m_findpyramidvolume(Fvec[Cvec[i].faceids[0]].area,
				Fvec[Cvec[i].faceids[0]].normal,side);

		side[0]=Cvec[i].coord[4][0]-Cvec[i].coord[1][0];
		side[1]=Cvec[i].coord[4][1]-Cvec[i].coord[1][1];
		side[2]=Cvec[i].coord[4][2]-Cvec[i].coord[1][2];

		vol += m_findpyramidvolume(Fvec[Cvec[i].faceids[3]].area,
				Fvec[Cvec[i].faceids[3]].normal,side);
		
		side[0]=Cvec[i].coord[4][0]-Cvec[i].coord[2][0];
		side[1]=Cvec[i].coord[4][1]-Cvec[i].coord[2][1];
		side[2]=Cvec[i].coord[4][2]-Cvec[i].coord[2][2];

		vol += m_findpyramidvolume(Fvec[Cvec[i].faceids[5]].area,
				Fvec[Cvec[i].faceids[5]].normal,side);

		Cvec[i].setvolume(vol);

	}	
}
//===============================================================
double Gridman::m_findpyramidvolume(double base_area,double base_normal[NDIM],double side[NDIM])
{
	double height;

	height = fabs(base_normal[0]*side[0]+base_normal[1]*side[1]+base_normal[2]*side[2]);
	return(0.33333*height*base_area);

}
//===============================================================
void Gridman::m_findboundaryfaces()
{

	int pos,i;
	bool found;
	bool allthreefound;

	for(int i=0;i<m_nboundaries;i++)
	{
		Face::st_intsort(m_bcpoints[i],m_pointsperbc[i]);
	}

	for(int b=0;b<m_nboundaries;b++)
	{
		for(unsigned int ind=0;ind<boundaryfaceids.size();ind++)
		{
			i = boundaryfaceids[ind];
			allthreefound=false;

			if(Fvec[i].lcellid == -1 || Fvec[i].rcellid == -1)
			{
				allthreefound=true;
				for(int k=0;k<NFPTS-1;k++)
				{
					binarysearch(m_bcpoints[b],Fvec[i].nodes[k],0,
							m_pointsperbc[b]-1,m_pointsperbc[b],pos,found);

					if(!found)
					{
						allthreefound=false;
						break;
					}

				}
			}

			if(allthreefound)
			{
				boundaryvec[b].Fvecids.push_back(i);	
			}
		}
	}

}
//===============================================================
void Gridman::writelogfile()
{
	double x,y,z;
	std::ofstream outfile("mesh.log");

	outfile<<"No: of nodes:"<<m_nnodes<<"\n";
	outfile<<"No: of faces:"<<m_nfaces<<"\n";
	outfile<<"No: of cells:"<<m_ncells<<"\n";

	outfile<<"Nodes\n";
	for(int i=0;i<m_nnodes;i++)
	{
		Nvec[i].getcoord(x,y,z);
		outfile<<i<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\n";
	}

	outfile<<"Faces\n";
	for(int i=0;i<m_nfaces;i++)
	{
		outfile<<i<<"\t"<<Fvec[i].nodes[0]<<"\t"<<Fvec[i].nodes[1]
			<<"\t"<<Fvec[i].nodes[2]<<"\t"<<Fvec[i].nodes[3]<<"\t"
			<<Fvec[i].lcellid<<"\t"<<Fvec[i].rcellid<<"\n";	
	}

	outfile<<"Cells\n";
	for(int i=0;i<m_ncells;i++)
	{
		outfile<<i<<"\t";
		for(int k=0;k<NCPTS;k++)
		{
			outfile<<Cvec[i].nodeids[k]<<"\t";
		}
		outfile<<"\n";
		outfile<<"\t";
		for(int k=0;k<NCFACES;k++)
		{
			outfile<<Cvec[i].faceids[k]<<"\t";
		}
		outfile<<"\n";
	}

	outfile<<"Face areas and normals\n";
	for(int i=0;i<m_nfaces;i++)
	{
		outfile<<i<<"\t"<<Fvec[i].area<<"\t"<<Fvec[i].normal[0]<<
			"\t"<<Fvec[i].normal[1]<<"\t"<<Fvec[i].normal[2]<<"\n";
	}

	outfile<<"Cell volumes\n";
	for(int i=0;i<m_ncells;i++)
	{
		outfile<<i<<"\t"<<Cvec[i].getvolume()<<"\n";
	}

	outfile<<"Boundaries\n";
	for(int i=0;i<m_nboundaries;i++)
	{
		outfile<<m_bcnames[i]<<" faces\n";
		for(unsigned int k=0;k<boundaryvec[i].Fvecids.size();k++)
		{
			outfile<<boundaryvec[i].Fvecids[k]<<"\n";
		}
	}


	outfile.close();
}
//===============================================================
void Gridman::binarysearch(int arr[],int val,int start,int end,int size,int &pos,bool &found)
{

	int mid,n;

	found = false;
	pos   = -1;

	mid = (start+end)/2;

	n = (end-start+1);

	if(val < arr[start] || val > arr[end])
	{
		found=false;
	}
	else
	{
		if(n==1)
		{
			if(val == arr[mid])
			{
				found=true;
				pos=mid;
			}
		}

		else
		{
			if(val > arr[mid])
			{
				binarysearch(arr,val,mid+1,end,size,pos,found);
			}
			else
				if(val < arr[mid])
				{
					binarysearch(arr,val,start,mid-1,size,pos,found);
				}
				else if(val == arr[mid])
				{
					found=true;
					pos=mid;
				}
		}
	}

}
//===============================================================
void Gridman::m_orientfacenormals()
{
	int findex;
	int lcell,rcell,cid;
	double vec[NDIM];
	double sign,sign1;

	//interior faces
	for(unsigned int i=0;i<interiorfaceids.size();i++)
	{
		findex=interiorfaceids[i];

		lcell = Fvec[findex].lcellid;
		rcell = Fvec[findex].rcellid;

		vec[0] = Cvec[rcell].centroid[0] - Cvec[lcell].centroid[0];
		vec[1] = Cvec[rcell].centroid[1] - Cvec[lcell].centroid[1];
		vec[2] = Cvec[rcell].centroid[2] - Cvec[lcell].centroid[2];

		sign = (vec[0]*Fvec[findex].normal[0] + vec[1]*Fvec[findex].normal[1] +
				vec[2]*Fvec[findex].normal[2] > 0)?1.0:-1.0;

		
		Fvec[findex].normal[0] *= sign;
		Fvec[findex].normal[1] *= sign;
		Fvec[findex].normal[2] *= sign;
	}

	//boundary faces
	for(unsigned int i=0;i<boundaryfaceids.size();i++)
	{
		findex = boundaryfaceids[i];

		lcell = Fvec[findex].lcellid;
		rcell = Fvec[findex].rcellid;

		cid   = (rcell == -1)?lcell:rcell;
		sign =  (rcell == -1)?1.0:-1.0;

		vec[0] = Fvec[findex].centroid[0] - Cvec[cid].centroid[0];
		vec[1] = Fvec[findex].centroid[1] - Cvec[cid].centroid[1];
		vec[2] = Fvec[findex].centroid[2] - Cvec[cid].centroid[2];

		sign1 = (vec[0]*Fvec[findex].normal[0] + vec[1]*Fvec[findex].normal[1] +
				vec[2]*Fvec[findex].normal[2] > 0)?1.0:-1.0;

		Fvec[findex].normal[0] *= sign1*sign;
		Fvec[findex].normal[1] *= sign1*sign;
		Fvec[findex].normal[2] *= sign1*sign;
	}

}
//===============================================================
void Gridman::printmeshvtu()
{
	char underscore;
	int offset;
	int intval;
	double x,y,z;

	std::ofstream outfile;
	outfile.open("mesh.vtu");

	offset=0;
	underscore='_';

	outfile<<"<?xml version=\"1.0\"?>\n";	
	outfile<<"<VTKFile type=\"UnstructuredGrid\" ";
	outfile<<"version=\"0.1\" byte_order=\"LittleEndian\">\n";	
	outfile<<"<UnstructuredGrid>\n";
	outfile<<"<Piece NumberOfPoints=\""<<m_nnodes<<"\" NumberOfCells=\""<<m_ncells<<"\">\n";

	/*outfile<<"<CellData>\n";
	outfile<<"<DataArray type=\"Float64\" Name=\"circumrad\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</CellData>\n";
	offset = offset + sizeof(int) + sizeof(double)*celldata;*/
	
	outfile<<"<Points>\n";
	outfile<<"<DataArray name=\"Position\" type=\"Float64\" ";
	outfile<<"NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
	outfile<<"</Points>\n";

	offset = offset + sizeof(int) + NDIM*sizeof(double)*m_nnodes;

	outfile<<"<Cells>\n";
	outfile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + NCPTS*sizeof(int)*m_ncells;
	outfile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	offset = offset + sizeof(int) + sizeof(int)*m_ncells;	
	outfile<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""<<offset<<"\"/>\n";

	outfile<<"</Cells>\n";
	outfile<<"</Piece>\n";
	outfile<<"</UnstructuredGrid>\n";
	outfile<<"<AppendedData encoding=\"raw\">\n";

	outfile.write(&underscore,sizeof(char));
	/*intval=sizeof(double)*m_ncells;
	outfile.write((char *)&intval,sizeof(int));
	for(unsigned int i=0;i<trianglelist.size();i++)
	{
		outfile.write((char *)&trianglelist[i].circrad,sizeof(double));	
	}*/

	//writing node data======================================
	intval=NDIM*sizeof(double)*m_nnodes;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<m_nnodes;i++)
	{
		Nvec[i].getcoord(x,y,z);
		
		outfile.write((char *)&x,sizeof(double));
		outfile.write((char *)&y,sizeof(double));
		outfile.write((char *)&z,sizeof(double));
	}
	//=======================================================

	//writing cell nodes=====================================
	intval=NCPTS*sizeof(int)*m_ncells;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<m_ncells;i++)
	{
		outfile.write((char *)&Cvec[i].nodeids[0],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[1],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[3],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[2],sizeof(int));
		
		outfile.write((char *)&Cvec[i].nodeids[4],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[5],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[7],sizeof(int));
		outfile.write((char *)&Cvec[i].nodeids[6],sizeof(int));
	}
	//=======================================================

	//writing cell offsets===================================
	intval=sizeof(int)*m_ncells;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<m_ncells;i++)
	{
		intval=(i+1)*NCPTS;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	//writing cell types=====================================
	intval=sizeof(int)*m_ncells;
	outfile.write((char *)&intval,sizeof(int));
	for(int i=0;i<m_ncells;i++)
	{
		intval=12;
		outfile.write((char *)&intval,sizeof(int));
	}	
	//=======================================================

	outfile.close();

}
//=======================================================================================
