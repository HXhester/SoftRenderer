#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

// Vertex, texcoord, normal, faces
class Model {
private:
    std::vector<Vec3f> verts;
    std::vector<Vec2f> texCoords;
    std::vector<Vec3f> normals;
    std::vector<std::vector<int>> faces_vert;
    std::vector<std::vector<int>> faces_texCoord;
    std::vector<std::vector<int>> faces_normal;

    TGAImage diffuseMap;

    void load_texture(const std::string filename, const std::string suffix, TGAImage& image);

public:
    Model(const char* filename);
    ~Model();

    int nverts();
    int nfaces();
    Vec3f vert(int i);
    std::vector<int> face_vert(int i);

    Vec2f texCoord(int iface, int ivert);
    Vec3f normal(int iface, int ivert);

    TGAColor diffuse(Vec2f texCoord);

    //Vec3f normal(int i);
};


#endif
