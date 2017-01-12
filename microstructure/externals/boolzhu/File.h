//#####################################################################
// File
// Bo Zhu, MIT, 06/15
//#####################################################################
#ifndef __File_h__
#define __File_h__
#include "Common.h"
#include "Lattice.h"
#include "Field.h"
#include "FaceField.h"
#include "Mesh.h"
#include "Particles.h"
#include "AleGrid.h"
#include <ostream>
#include <istream>
#include <fstream>
#include <iostream>
#include <cerrno>
#include <climits>
#include <cstdio>
#ifdef WIN32
#include <windows.h>
#endif
#ifdef USE_TINY_OBJ_LOADER
#include "tiny_obj_loader.h"
#endif

const static bool file_verbose=false;
////binary io
template<class T_VAL> inline void Write_Binary(std::ostream& output,const T_VAL& data)
{output.write(reinterpret_cast<const char*>(&data),sizeof(T_VAL));}

template<class T_VAL> inline void Read_Binary(std::istream& input,T_VAL& data)
{input.read(reinterpret_cast<char*>(&data),sizeof(T_VAL));}

template<class T_VAL> inline void Write_Binary_Array(std::ostream& output,const T_VAL* array,const int n)
{if(n>0)output.write(reinterpret_cast<const char*>(array),n*sizeof(T_VAL));}

template<class T_VAL> inline void Read_Binary_Array(std::istream& input,T_VAL* array,const int n)
{if(n>0)input.read(reinterpret_cast<char*>(array),n*sizeof(T_VAL));}

template<class T_VAL> inline bool Write_Binary_To_File(const std::string& file_name,const T_VAL& data)
{std::ofstream output(file_name,std::ios::binary);if(!output)return false;Write_Binary(output,data);return true;}

template<class T_VAL> inline bool Read_Binary_From_File(const std::string& file_name,T_VAL& data)
{std::ifstream input(file_name,std::ios::binary);
if(!input)return false;
Read_Binary(input,data);
input.clear();input.close();
return true;}

////text io
template<class T_VAL> inline void Write_Text(std::ostream& output,const T_VAL& data){output<<data;}

template<class T_VAL> inline void Read_Text(std::istream& input,T_VAL& data){input>>data;}

template<class T_ARRAY> inline void Write_Text_Array(std::ostream& output,const T_ARRAY& array,const int n,char separator='\n')
{for(int i=0;i<n;i++)output<<array[i]<<separator;}

template<class T_ARRAY> inline void Read_Text_Array(std::istream& input,T_ARRAY& array,const int n,char separator='\n')
{for(int i=0;i<n;i++)input>>array[i];}

template<class T_VAL> inline bool Write_Text_To_File(const std::string& file_name,const T_VAL& data)
{std::ofstream output(file_name);if(!output)return false;Write_Text(output,data);return true;}

template<class T_VAL> inline bool Read_Text_From_File(const std::string& file_name,T_VAL& data)
{std::ifstream input(file_name);if(!input)return false;Read_Text(input,data);return true;}

template<class T_VAL> inline bool Write_Text_Array_To_File(const std::string& file_name,const T_VAL& array,const int n,char separator='\n')
{std::ofstream output(file_name);if(!output)return false;Write_Text_Array(output,array,n,separator);return true;}

template<class T_ARRAY> inline bool Read_Text_Array_From_File(const std::string& file_name,T_ARRAY& array,const int n)
{std::ifstream input(file_name);if(!input)return false;Read_Text_Array(input,array,n);return true;}

template<class T_VAL> inline bool Append_Text_To_File(const std::string& file_name,const T_VAL& data)
{std::ofstream output(file_name,std::ios_base::app);if(!output)return false;Write_Text(output,data);return true;}

//template<class T_VAL> inline void Write_To_Text_File(const std::string& filename,const T_VAL& data)
//{std::ostream output(filename);output<<d1;delete output;}

////file operations
#ifdef WIN32
inline bool Directory_Exists(const char* dirname)
{DWORD attr=GetFileAttributes(dirname);return((attr!=-1)&&(attr&FILE_ATTRIBUTE_DIRECTORY));}

inline bool Create_Directory(const std::string& dirname)
{
    if(!Directory_Exists(dirname.c_str())){size_t pos=0;
        do{pos=dirname.find_first_of("\\/",pos+1);
        if(!Directory_Exists(dirname.substr(0,pos).c_str())){
            if(CreateDirectory(dirname.substr(0,pos).c_str(),NULL)==0 && ERROR_ALREADY_EXISTS!=GetLastError()){
                std::cerr<<"Failed!"<<std::endl;return false;}}}while(pos!=std::string::npos);}
    return true;
}
#endif

inline bool File_Exists(const std::string& file_name)
{return !std::ifstream(file_name.c_str()).fail();}

inline std::string File_Extension_Name(const std::string& file_name)
{size_t p=file_name.rfind('.');if(p!=std::string::npos)return file_name.substr(p+1);else return "";}

////Lattice io
template<int d> inline void Write_Binary(std::ostream& output,const Lattice<d>& lattice)
{Write_Binary(output,lattice.cell_counts);Write_Binary(output,lattice.dx);Write_Binary(output,lattice.domain_min);}

template<int d> inline void Read_Binary(std::istream& input,Lattice<d>& lattice)
{Read_Binary(input,lattice.cell_counts);Read_Binary(input,lattice.dx);Read_Binary(input,lattice.domain_min);
lattice.Initialize(lattice.cell_counts,lattice.dx,lattice.domain_min);}

////Field io
template<class T_VAL,int d> inline void Write_Binary(std::ostream& output,const Field<T_VAL,d>& field)
{Write_Binary(output,field.counts);Write_Binary_Array(output,&field.array[0],(int)field.array.size());}

template<class T_VAL,int d> inline void Read_Binary(std::istream& input,Field<T_VAL,d>& field)
{Read_Binary(input,field.counts);field.Resize(field.counts);Read_Binary_Array(input,&field.array[0],(int)field.array.size());}

////FaceField io
template<class T_VAL,int d> inline void Write_Binary(std::ostream& output,const FaceField<T_VAL,d>& face_field)
{Write_Binary(output,face_field.mac_grid.grid.cell_counts);for(int i=0;i<d;i++)Write_Binary(output,face_field.face_fields[i]);}

template<class T_VAL,int d> inline void Read_Binary(std::istream& input,FaceField<T_VAL,d>& face_field)
{Type_Define_VectorDi(d);VectorDi cell_counts;Read_Binary(input,cell_counts);face_field.Resize(cell_counts);
for(int i=0;i<d;i++)Read_Binary(input,face_field.face_fields[i]);}

////Mesh io
template<int d,int e_d> inline void Write_Binary_Simplicial_Mesh(std::ostream& output,const SimplicialMesh<d,e_d>& mesh)
{int vtx_n=(int)mesh.vertices.size();Write_Binary(output,vtx_n);Write_Binary_Array(output,&mesh.vertices[0],vtx_n);
int e_n=(int)mesh.elements.size();Write_Binary(output,e_n);if(e_n>0)Write_Binary_Array(output,&mesh.elements[0],e_n);}

template<int d,int e_d> inline void Read_Binary_Simplicial_Mesh(std::istream& input,SimplicialMesh<d,e_d>& mesh)
{int vtx_n=0;Read_Binary(input,vtx_n);mesh.vertices.resize(vtx_n);Read_Binary_Array(input,&mesh.vertices[0],vtx_n);
int e_n=0;Read_Binary(input,e_n);if(e_n>0){mesh.elements.resize(e_n);Read_Binary_Array(input,&mesh.elements[0],e_n);}}

#define Read_Write_Mesh(MeshType) \
template<int d> inline void Write_Binary(std::ostream& output,const MeshType& mesh) \
{Write_Binary_Simplicial_Mesh(output,mesh);} \
template<int d> inline void Read_Binary(std::istream& input,MeshType& mesh) \
{Read_Binary_Simplicial_Mesh(input,mesh);} \
template<int d> inline void Write_Text(std::ostream& output,const MeshType& mesh) \
{Write_Text_Simplicial_Mesh(output,mesh);} \
template<int d> inline void Read_Text(std::istream& input,MeshType& mesh) \
{Read_Text_Simplicial_Mesh(input,mesh);}

Read_Write_Mesh(VolumetricMesh<d>)
Read_Write_Mesh(SurfaceMesh<d>)
Read_Write_Mesh(TetrahedronMesh<d>)
Read_Write_Mesh(TriangleMesh<d>)
Read_Write_Mesh(SegmentMesh<d>)
Read_Write_Mesh(QuadMesh<d>)

////mesh io, text file
template<int d,int e_d> inline void Write_Text_Simplicial_Mesh(std::ostream& output,const SimplicialMesh<d,e_d>& mesh)
{int vtx_n=(int)mesh.vertices.size();Write_Text(output,vtx_n);Write_Text(output,'\n');
if(vtx_n>0){for(int i=0;i<vtx_n;i++){Write_Text_Array(output,mesh.vertices[i],d,' ');Write_Text(output,'\n');}}
int e_n=(int)mesh.elements.size();Write_Text(output,'\n');Write_Text(output,e_n);Write_Text(output,'\n');
if(e_n>0){for(int i=0;i<e_n;i++){Write_Text_Array(output,mesh.elements[i],e_d,' ');Write_Text(output,'\n');}}}

template<int d,int e_d> inline void Read_Text_Simplicial_Mesh(std::istream& input,SimplicialMesh<d,e_d>& mesh)
{int vtx_n=0;Read_Text(input,vtx_n);
if(vtx_n>0){mesh.vertices.resize(vtx_n);for(int i=0;i<vtx_n;i++)Read_Text_Array(input,mesh.vertices[i],d);}
int e_n=0;Read_Text(input,e_n);
if(e_n>0){mesh.elements.resize(e_n);for(int i=0;i<e_n;i++)Read_Text_Array(input,mesh.elements[i],e_d);}}

////Points io
template<int d> inline void Write_Binary(std::ostream& output,const Points<d>& points)
{	
	int n=points.Size();Write_Binary(output,n);
	if(n>0){
		Write_Binary_Array(output,&(*points.X())[0],n);}
}
template<int d> inline void Read_Binary(std::istream& input,Points<d>& points)
{
	int n=0;Read_Binary(input,n);
	points.Resize(n);
	if(n>0){Read_Binary_Array(input,&(*points.X())[0],n);}
	points.Initialize_Attributes();
}

////Particles io
template<int d> inline void Write_Binary(std::ostream& output,const Particles<d>& particles)
{	
	int n=particles.Size();Write_Binary(output,n);
	if(n>0){
		Write_Binary_Array(output,&(*particles.X())[0],n);
		Write_Binary_Array(output,&(*particles.V())[0],n);
		Write_Binary_Array(output,&(*particles.F())[0],n);
		Write_Binary_Array(output,&(*particles.Mass())[0],n);
		Write_Binary_Array(output,&(*particles.Color())[0],n);
	}
}
template<int d> inline void Read_Binary(std::istream& input,Particles<d>& particles)
{
	int n=0;Read_Binary(input,n);
	particles.Resize(n);
	if(n>0){Read_Binary_Array(input,&(*particles.X())[0],n);
		Read_Binary_Array(input,&(*particles.V())[0],n);
		Read_Binary_Array(input,&(*particles.F())[0],n);
		Read_Binary_Array(input,&(*particles.Mass())[0],n);
		Read_Binary_Array(input,&(*particles.Color())[0],n);}
	particles.Initialize_Attributes();
}

////AleGrid io
template<int d> inline void Write_Binary(std::ostream& output,const AleGrid<d>& ale_grid)
{
	Write_Binary(output,ale_grid.grid);
	Write_Binary(output,ale_grid.x);
	Write_Binary(output,ale_grid.use_singularity);
	if(ale_grid.use_singularity){
		Write_Binary(output,ale_grid.cell_singularity);
		Write_Binary_Array(output,&ale_grid.head[0],(int)ale_grid.head.size());
		Write_Binary_Array(output,&ale_grid.next[0],(int)ale_grid.next.size());}
}

template<int d> inline void Read_Binary(std::istream& input,AleGrid<d>& ale_grid)
{
	Read_Binary(input,ale_grid.grid);
	std::cout<<"read ale_grid: "<<ale_grid.grid.cell_counts.transpose()<<std::endl;
	Read_Binary(input,ale_grid.x);
	Read_Binary(input,ale_grid.use_singularity);
	if(ale_grid.use_singularity){
		ale_grid.Initialize_Singularity();
		Read_Binary(input,ale_grid.cell_singularity);
		Read_Binary_Array(input,&ale_grid.head[0],(int)ale_grid.head.size());
		Read_Binary_Array(input,&ale_grid.next[0],(int)ale_grid.next.size());}
}

////Read obj
#ifdef USE_TINY_OBJ_LOADER
inline void Read_From_Obj_File(const std::string& file_name,/*rst*/Array<SurfaceMesh<2>*>& meshes)
{std::cerr<<"Read_From_Obj_File<2> not implemented"<<std::endl;}

inline void Read_From_Obj_File(const std::string& file_name,/*rst*/Array<SurfaceMesh<3>*>& meshes)
{
    Array<tinyobj::shape_t> shapes;
    Array<tinyobj::material_t> materials;
	std::string base;size_t l;
	if ((l=file_name.find_last_of('/')) != std::string::npos)
		base=file_name.substr(0, l+1);
	else if((l=file_name.find_last_of('\\'))!=std::string::npos)
		base=file_name.substr(0, l+1);
	std::string err=tinyobj::LoadObj(shapes,materials,file_name.c_str(),base.c_str());
    
    std::cout<<"#shapes="<<shapes.size()<<" #materials="<<materials.size()<<std::endl;
    
	meshes.resize((int)shapes.size());
    for(size_type i=0;i<shapes.size();i++){
		meshes[i]=new SurfaceMesh<3>();
        tinyobj::mesh_t& mesh=shapes[i].mesh;
        for(size_type j=0;j<mesh.positions.size()/3;j++){
            meshes[i]->vertices.push_back(Vector3((T)mesh.positions[j*3],
			(T)mesh.positions[j*3+1],(T)mesh.positions[j*3+2]));}
        for(size_type j=0;j<mesh.indices.size()/3;j++){
			meshes[i]->elements.push_back(Vector3i(mesh.indices[j*3],
			mesh.indices[j*3+1],mesh.indices[j*3+2]));}}
}
#endif

////write command line to txt file
inline void Write_To_Text_File(const std::string& file_name,int argc,char* argv[])
{std::ofstream output(file_name,std::ios::out);if(!output)return;for(int i=0;i<argc;i++){output<<argv[i]<<" ";}output.close();}

////output functions for opengl viewer (basically convert everything to 3D in order to be visualized in the viewer)
template<int d> void Write_Lattice(const std::string& file_name,const Lattice<d>& lattice)
{
	Lattice<3> lattice_to_write;lattice.Convert_To(lattice_to_write);
    Write_Binary_To_File(file_name,lattice_to_write);
    if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Read_Lattice(const std::string& file_name,Lattice<d>& lattice)
{
	Lattice<3> lattice_to_read;Read_Binary_From_File(file_name,lattice_to_read);
	lattice_to_read.Convert_To(lattice);
    if(file_verbose)std::cout<<"read from file "<<file_name<<std::endl;
}

template<class ArrayT,class TV_INT> void Write_Indirect_Vector_Field(const std::string& file_name,const ArrayT& u,const Array<int>& map,const TV_INT& counts)
{
	const int d=(int)counts.size();
    Vector3i counts_3d;Convert_Vector(counts,counts_3d,1);
    Field<Vector3,3> displacements(counts_3d);displacements.Fill(Vector3::Zero());
    Copy_Indirect(u,d,map,displacements.array);
    Write_Binary_To_File(file_name,displacements);
    if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<class ArrayT> void Write_VectorN_To_Field1D(const std::string& file_name,const ArrayT& u,const int d)
{
	const int n=(int)u.size()/d;Field<Vector3,1> displacements(Vec1i(n));displacements.Fill(Vector3::Zero());
	Copy(u,d,displacements.array);Write_Binary_To_File(file_name,displacements);
}

template<class T_VAL> void Write_Array_To_Field1D(const std::string& file_name,const Array<T_VAL>& array)
{
	int n=(int)array.size();Field<T_VAL,1> field(Vec1i(n));field.array=array;
	Write_Binary_To_File(file_name,field);
}

template<class T_VAL,class TV_INT> void Write_Field(const std::string& file_name,const Array<T_VAL>& array,const TV_INT& cell_counts)
{
    Vector3i counts_3d;Convert_Vector(cell_counts,counts_3d,1);
    Field<T_VAL,3> field_3d(counts_3d);
    std::copy(array.begin(),array.end(),field_3d.array.begin());
    Write_Binary_To_File(file_name,field_3d);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<class T_VAL,int d> void Write_Vector_Field(const std::string& file_name,const Array<Vector<T_VAL,d> >& array,const Vector<int,d>& cell_counts)
{
	if(d==3)Write_Field(file_name,array,cell_counts);
	else if(d==2){	////d=2
		Vector3i counts_3d;Convert_Vector(cell_counts,counts_3d,1);
		Field<Vector<T_VAL,3>,3> vector_field_3d(counts_3d);vector_field_3d.Fill(Vector<T_VAL,3>::Zero());	
		for(size_type i=0;i<array.size();i++)for(int j=0;j<d;j++)vector_field_3d.array[i][j]=array[i][j];
		Write_Binary_To_File(file_name,vector_field_3d);
		if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
	}
	else{std::cout<<"Write_Vector_Field supports 2d and 3d only!"<<std::endl;}////d=1
}

template<class T_MATRIX,int d> void Write_Tensor_Field(const std::string& file_name,const Array<T_MATRIX>& array,const Vector<int,d>& cell_counts)
{
	if(d==3)Write_Field(file_name,array,cell_counts);
	else if(d==2){	////d=2
		Vector3i counts_3d;Convert_Vector(cell_counts,counts_3d,1);
		Field<Matrix3,3> tensor_field_3d(counts_3d);tensor_field_3d.Fill(Matrix3::Zero());	
		for(size_type i=0;i<array.size();i++)for(int j=0;j<d;j++)for(int k=0;k<d;k++)tensor_field_3d.array[i](j,k)=(T)array[i](j,k);
		Write_Binary_To_File(file_name,tensor_field_3d);
		if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
	}
	else{std::cout<<"Write_Tensor_Field supports 2d and 3d only!"<<std::endl;}	////d=1
}

template<class T_VAL,class TV_INT> void Read_Field(const std::string& file_name,Array<T_VAL>& array,TV_INT& cell_counts)
{
    Field<T_VAL,3> field_3d;Read_Binary_From_File(file_name,field_3d);
	Convert_Vector(field_3d.counts,cell_counts,1);array.resize(field_3d.array.size());
	std::copy(field_3d.array.begin(),field_3d.array.end(),array.begin());
	if(file_verbose)std::cout<<"read from file "<<file_name<<std::endl;
}

template<class T_VAL,int d> void Read_Vector_Field(const std::string& file_name,Array<Vector<T_VAL,d> >& array,Vector<int,d>& cell_counts)
{
	if(d==3)Read_Field(file_name,array,cell_counts);
	else if(d==2){	////d=2
		Vector3i counts_3d;Field<Vector3,3> vector_field_3d;Read_Field(file_name,vector_field_3d.array,counts_3d);
		Convert_Vector(counts_3d,cell_counts,1);
		array.resize(vector_field_3d.array.size());
		for(size_type i=0;i<array.size();i++)for(int j=0;j<d;j++)array[i][j]=(T)vector_field_3d.array[i][j];
		if(file_verbose)std::cout<<"read from file "<<file_name<<std::endl;}
	else{std::cout<<"Read_Vector_Field 1D not implemented"<<std::endl;}
}

template<class T_MATRIX,int d> void Read_Tensor_Field(const std::string& file_name,Array<T_MATRIX>& array,Vector<int,d>& cell_counts)
{
	if(d==3)Read_Field(file_name,array,cell_counts);
	else if(d==2){	////d=2
		Vector3i counts_3d;Field<Matrix3,3> tensor_field_3d;Read_Field(file_name,tensor_field_3d.array,counts_3d);
		Convert_Vector(counts_3d,cell_counts,1);
		array.resize(tensor_field_3d.array.size());
		for(size_type i=0;i<array.size();i++)for(int j=0;j<d;j++)for(int k=0;k<d;k++)array[i](j,k)=(T)tensor_field_3d.array[i](j,k);
		if(file_verbose)std::cout<<"read from file "<<file_name<<std::endl;}
	else{std::cout<<"Read_Tensor_Field 1D not implemented"<<std::endl;}
}

template<int d> void Write_Face_Field_To_Node_Field(const std::string& file_name,const FaceField<T,d>& velocity,const MacGrid<d>& mac_grid)
{Type_Define_VectorD(d);
	Vector3i counts_3d;Convert_Vector(mac_grid.grid.node_counts,counts_3d,1);
	Field<Vector3,3> v3(counts_3d);v3.Fill(Vector3::Zero());
	Interpolation<d> intp(mac_grid);Field<VectorD,d> node_velocity(mac_grid.grid.node_counts);
	intp.Interpolate_Faces_To_Nodes(velocity,node_velocity);
	Convert_Vector_Array(node_velocity.array,v3.array,(T)0);
	Write_Binary_To_File(file_name,v3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Face_Field_To_Cell_Field(const std::string& file_name,const FaceField<T,d>& velocity,const MacGrid<d>& mac_grid)
{Type_Define_VectorD(d);
	Vector3i counts_3d;Convert_Vector(mac_grid.grid.cell_counts,counts_3d,1);
	Interpolation<d> intp(mac_grid);Field<VectorD,d> cell_velocity(mac_grid.grid.cell_counts);
	intp.Interpolate_Faces_To_Cells(velocity,cell_velocity);
	Field<Vector3,3> v3(counts_3d);Convert_Vector_Array(cell_velocity.array,v3.array,(T)0);
	Write_Binary_To_File(file_name,v3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<class TV,int d> void Write_Cell_Vector_Field_To_Node_Field(const std::string& file_name,const Field<TV,d>& velocity,const Lattice<d>& grid)
{Type_Define_VectorD(d);
	Vector3i counts_3d;Convert_Vector(grid.node_counts,counts_3d,1);
	Field<Vector3,3> v3(counts_3d);v3.Fill(Vector3::Zero());
	Interpolation<d> intp(grid);Field<VectorD,d> node_velocity(grid.node_counts);
	intp.Interpolate_Cells_To_Nodes(velocity,node_velocity);
	Convert_Vector_Array(node_velocity.array,v3.array,(T)0);
	Write_Binary_To_File(file_name,v3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<class T_VAL,int d> void Write_Face_Field(const std::string& file_name,const FaceField<T_VAL,d>& field,const MacGrid<d>& mac_grid)
{Type_Define_VectorD(d);
	Vector3i counts_3d;Convert_Vector(mac_grid.grid.cell_counts,counts_3d,1);
	FaceField<T_VAL,3> field_3(counts_3d);field_3.Fill((T_VAL)0);
	for(int i=0;i<d;i++)std::copy(field.face_fields[i].array.begin(),field.face_fields[i].array.end(),field_3.face_fields[i].array.begin());
	Write_Binary_To_File(file_name,field_3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Points(const std::string& file_name,const Points<d>& points)
{Type_Define_VectorD(d);
	Points<3> points_3d;points_3d.Resize((int)points.Size());
	Convert_Vector_Array(*points.X(),*points_3d.X(),(T)0);
	Write_Binary_To_File(file_name,points_3d);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Particles(const std::string& file_name,const Particles<d>& particles)
{Type_Define_VectorD(d);
	Particles<3> particles_3d;particles_3d.Resize((int)particles.Size());
	Convert_Vector_Array(*particles.X(),*particles_3d.X(),(T)0);
	Convert_Vector_Array(*particles.V(),*particles_3d.V(),(T)0);
	Convert_Vector_Array(*particles.F(),*particles_3d.F(),(T)0);
	std::copy(particles.Mass()->begin(),particles.Mass()->end(),particles_3d.Mass()->begin());
	std::copy(particles.Color()->begin(),particles.Color()->end(),particles_3d.Color()->begin());
	Write_Binary_To_File(file_name,particles_3d);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Mesh(const std::string& file_name,const SurfaceMesh<d>& mesh)
{typedef typename IF<d==2,SegmentMesh<3>,TriangleMesh<3> >::TYPE Mesh3Type;
	Mesh3Type mesh3;Convert_Mesh(mesh,mesh3);
	Write_Binary_To_File(file_name,mesh3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Mesh(const std::string& file_name,const VolumetricMesh<d>& mesh)
{typedef typename IF<d==2,TriangleMesh<3>,VolumetricMesh<3> >::TYPE Mesh3Type;
	Mesh3Type mesh3;Convert_Mesh(mesh,mesh3);
	Write_Binary_To_File(file_name,mesh3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

template<int d> void Write_Mesh(const std::string& file_name,const QuadMesh<d>& mesh)
{typedef QuadMesh<3> Mesh3Type;
	Mesh3Type mesh3;Convert_Mesh(mesh,mesh3);
	Write_Binary_To_File(file_name,mesh3);
	if(file_verbose)std::cout<<"write to file "<<file_name<<std::endl;
}

#endif