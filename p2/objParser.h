#pragma once
//remember to delete these out of main
#include <string>
#include <vector>

struct Face
{
	int v;
	int vt;
	int vn;
	Face(int _v, int _vt, int _vn) 
	:
	v(_v),
	vt(_vt),
	vn(_vn)
	{}
};

struct Vert
{
	double x;
	double y;
	double z;
	Vert(double _x, double _y, double _z)
	:
	x(_x),
	y(_y),
	z(_z)
	{}

};

class ObjParser
{
public:
	//singleton
	inline static ObjParser& openFile(std::string&& fileName) 
	{ 
		static ObjParser toBeReturned(fileName.c_str()); 
		return toBeReturned; 
	}

	inline const std::vector<Vert> getVerts() const { return this->verts; }
	inline const std::vector<Face>   getFaces() const { return this->faces; }

	//deleted funcs
	ObjParser() = delete;
	ObjParser(const ObjParser&) = delete;
	ObjParser(ObjParser&&) = delete;
	ObjParser& operator=(const ObjParser&) = delete;
	ObjParser& operator=(ObjParser&&) = delete;

private:
	//read in file, fill verts and faces
	ObjParser(const char* fileName);


	std::vector<Vert> verts;
	std::vector<Face> faces;
};
