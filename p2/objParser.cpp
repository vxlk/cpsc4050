#include "objParser.h"
#include <fstream>
#include <iostream>

constexpr auto VERTEX = "v";
constexpr auto FACE = "f";

ObjParser::ObjParser(const char* fileName)
{
	//open file
	//read verts -> emplace into array
	//read faces -> emplace into array
	//skip smoothing groups for now i guess idk
	std::ifstream objFile(fileName);
	if (!objFile) std::cout << "bad file name\n";

	//for now .. ignore everything but the vertices
	std::string buffer = "";
	while (objFile >> buffer)
	{
		if (buffer == VERTEX)
		{
			//for 3 iterations .. read into the verts
			//x
			objFile >> buffer;
			double x = std::stod(buffer);

			//y
			objFile >> buffer;
			double y = std::stod(buffer);

			//z
			objFile >> buffer;
			double z = std::stod(buffer);

			verts.emplace_back(Vert(x, y, z));
		}

		if (buffer == FACE)
		{
			//for v, vt, && vn
			for (int index = 0; index < 3; ++index)
			{
				//setup variables
				std::string readInNumber = "";
				std::vector<std::string> container;
				container.reserve(3);

				//read into buffer
				objFile >> buffer;

				//parse out / -> get integers out
				int i = 0;
				while (buffer[i] != '\0')
				{
					while (buffer[i] != '/' && buffer[i] != '\0')
					{
						readInNumber += buffer[i];
						++i;
					}
					container.emplace_back(readInNumber);
					readInNumber.clear();
					if (buffer[i] != '\0') ++i;
				}

				//transfer to face container
				faces.emplace_back(Face(stoi(container[0]), stoi(container[1]), stoi(container[2])));
			}
		}
	}

	
	//<debug>-----------------------------------------------------------------------
	for (int i = 1; i <= faces.size(); ++i)
	{
		std::cout << faces[i-1].v << "/" << faces[i-1].vt << "/" << faces[i-1].vn;
		if (i % 3 == 0) std::cout << "\n";
		else std::cout << " ";
	}

	for (int i = 0; i < verts.size(); ++i)
		std::cout << verts[i].x << " " << verts[i].y << " " << verts[i].z << "\n";
	//</debug>----------------------------------------------------------------------
	
}