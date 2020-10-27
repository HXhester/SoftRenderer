#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char* filename) :verts(), texCoords(), normals(), faces_vert(), faces_texCoord(), faces_normal(), diffuseMap() {

	std::ifstream input;
	input.open(filename, std::ifstream::in);
	if (input.fail()) {
		printf("fail to read.\n");
		return;
	}

	std::string line;
	while (!input.eof())
	{
		std::getline(input, line);
		std::istringstream iss(line.c_str());

		char trash;
		if (!line.compare(0, 2, "v ")) {
			Vec3f v;
			iss >> trash;
			iss >> v.x >> v.y >> v.z;

			verts.push_back(v);
		}
		else if (!line.compare(0, 3, "vt ")) {
			Vec2f tex;
			iss >> trash >> trash;
			iss >> tex.x >> tex.y;

			texCoords.push_back(tex);
		}
		else if (!line.compare(0, 3, "vn ")) {
			Vec3f norm;
			iss >> trash >> trash;
			iss >> norm.x >> norm.y >> norm.z;

			normals.push_back(norm);
		}
		else if (!line.compare(0, 2, "f ")) {
			int ivert, itex, inormal;
			std::vector<int> face_vert;
			std::vector<int> face_tex;
			std::vector<int> face_norm;
			iss >> trash;
			while (iss >> ivert >> trash >> itex >> trash >> inormal) {
				ivert--;
				itex--;
				inormal--;
				face_vert.push_back(ivert);
				face_tex.push_back(itex);
				face_norm.push_back(inormal);
			}

			faces_vert.push_back(face_vert);
			faces_texCoord.push_back(face_tex);
			faces_normal.push_back(face_norm);
		}
	}

	std::cerr << "#v #" << verts.size() << " #f #" << faces_vert.size() << std::endl;

	load_texture(filename, "_diffuse.tga", diffuseMap);
}

Model::~Model() {
}

int Model::nfaces() {
	return faces_vert.size();
}

int Model::nverts() {
	return verts.size();
}

std::vector<int> Model::face_vert(int i) {
	return faces_vert[i];
}

Vec3f Model::vert(int i) {
	return verts[i];
}

Vec2f Model::texCoord(int iface, int ivert) {
	return texCoords[faces_texCoord[iface][ivert]];
}

Vec3f Model::normal(int iface, int ivert) {
	return normals[faces_normal[iface][ivert]];
}

void Model::load_texture(std::string filename, std::string suffix, TGAImage& image) {
	size_t dot = filename.find_last_of(".");
	if (dot == std::string::npos) return;
	std::string texfile = filename.substr(0, dot) + suffix;
	std::cerr << "texture file " << texfile << " loading " << (image.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
	image.flip_vertically();
}

// µã²ÉÑù
TGAColor Model::diffuse(Vec2f uv) {
	return diffuseMap.get(uv.u * diffuseMap.get_width(), uv.v * diffuseMap.get_height());
}