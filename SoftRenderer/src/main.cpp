#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

#define PI 3.14159265

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const int width = 800;
const int height = 800;
const int depth = 255;

int* zBuffer = NULL;
Model* model = NULL;
const char* filename = "resources/african_head.obj";
Matrix MVP = Matrix::identity(4);

float fov = 60.f;
float near = 0.1;
float far = 100.f;

Vec3f lightDir(1, -1, 1);
Vec3f camera(0, 0, 4);
Vec3f origin(0, 0, 0);

Vec3f barycentric(Vec3f* pts, Vec2i P) {
	Vec3f u = Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]) ^ Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]);
	/* `pts` and `P` has integer value as coordinates
	   so `abs(u[2])` < 1 means `u[2]` is 0, that means
	   triangle is degenerate, in this case return something with negative coordinates */
	if (std::abs(u[2]) < 1) return Vec3f(-1, 1, 1);
	return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	bool steep = false;

	if (abs(y1 - y0) > abs(x1 - x0)) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}

	if (x0 > x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}


	int dx = x1 - x0;
	int dy = y1 - y0;

	float error = 0;
	float diffErr = abs(dy) * 2;

	int y = y0;
	for (int x = x0; x <= x1; x++) {
		error += diffErr;
		if (abs(error) > dx) {
			y += (dy > 0) ? 1 : -1;
			error -= dx * 2;
		}

		if (steep) {
			image.set(y, x, color);
		}
		else
		{
			image.set(x, y, color);
		}
	}

}

// TODO: improve z buffer mapping
void triangle(Vec3f* pts, int* zBuffer, TGAImage& image,Vec2f* texCoord, float* intensity) {

	Vec2f bboxmin(image.get_width() - 1, image.get_height() - 1);
	Vec2f bboxmax(0, 0);
	Vec2f bboxclamp(image.get_width() - 1, image.get_height() - 1);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::max(0.f, std::min(pts[i][j], bboxmin[j]));
			bboxmax[j] = std::min(bboxclamp[j], std::max(pts[i][j], bboxmax[j]));
		}
	}

	Vec2i P;
	for (P.x = bboxmin.x; P.x < bboxmax.x + 1; P.x++) {
		for (P.y = bboxmin.y; P.y < bboxmax.y + 1; P.y++) { 
			Vec3f bary = barycentric(pts, P);
			if (bary.x >= 0 && bary.y >= 0 && bary.z >= 0) {
				float z = bary[0] * pts[0].z + bary[1] * pts[1].z + bary[2] * pts[2].z;

				Vec2f uv = texCoord[0]*bary[0] + texCoord[1] * bary[1] + texCoord[2] * bary[2];
				float intens= intensity[0]*bary[0] + intensity[1] * bary[1] + intensity[2] * bary[2];
				if (zBuffer[P.x + P.y * width] < z) {
					zBuffer[P.x + P.y * width] = z;
					
					TGAColor c = model->diffuse(uv);
					image.set(P.x, P.y, TGAColor(255,255,255) * intens);
				}			
			}
		}
	}
}

Vec3f world2Screen(Vec3f p) {
	return Vec3f((p.x + 1.) * width / 2, (p.y + 1.) * height / 2, p.z);
}

//TODO: check why this is performance heavy
void mulMatrixVector(Matrix& m, const Vec3f& v, Vec3f& out) {
	out.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3];
	out.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3];
	out.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3];
	float w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3];

	if (w != 1) {
		out.x /= w;
		out.y /= w;
		out.z /= w;
	}
}

void drawWireframe(TGAImage & image) {
	for (int i = 0; i < model->nfaces(); i++) {
		std::vector<int> vertices = model->face_vert(i);
		for (int j = 2; j >= 0; j--) {
			Vec3f v0 = model->vert(vertices[j]);
			Vec3f v1 = model->vert(vertices[(j + 1) % 3]);

			int x0 = (v0.x + 1.) * width / 2.;
			int y0 = (v0.y + 1.) * height / 2.;

			int x1 = (v1.x + 1.) * width / 2.;
			int y1 = (v1.y + 1.) * height / 2.;

			line(x0, y0, x1, y1, image, white);
		}
	}
}

void drawSolid(TGAImage& image, int* zBuffer) {

	for (int i = 0; i < model->nfaces(); i++) {
		std::vector<int> vertices = model->face_vert(i);
		Vec3f worldCoords[3];
		Vec3f screenCoords[3];
		Vec2f texCoords[3];
		float intensity[3];
		for (int j = 2; j >= 0; j--) {
			Vec3f v = model->vert(vertices[j]);
			worldCoords[j] = v;
			screenCoords[j] = Vec3f(MVP * v);
			//screenCoords[j] = world2Screen(Vec3f(v.x / (1 - v.z / 5.f), v.y / (1 - v.z / 5.f), v.z / (1 - v.z / 5.f)));
			texCoords[j] = model->texCoord(i, j);
			intensity[j] = model->normal(i, j) * lightDir;
		}

		triangle(screenCoords, zBuffer, image, texCoords, intensity);
	}
}

void setProjectionMatrix(const float& fov, const float& near, const float& far, Matrix& M) {
	float scale = 1.f / tan(fov * 0.5 * PI / 180);
	M[0][0] = scale;
	M[1][1] = scale;
	M[2][2] = -far / (far - near); // used to remap z to [0,1] 
	M[3][2] = -far * near / (far - near); // used to remap z [0,1] 
	M[2][3] = -1; // set w = -z 
	M[3][3] = 0;
}

// WIP transform
Matrix modelM(Vec3f pos, Vec3f rot, Vec3f scale) {
	Matrix T = Matrix::identity(4);
	Matrix Rx = Matrix::identity(4);
	Matrix Ry = Matrix::identity(4);
	Matrix Rz = Matrix::identity(4);
	Matrix S = Matrix::identity(4);

	T[0][3] = pos.x;
	T[1][3] = pos.y;
	T[2][3] = pos.z;

	S[0][0] = scale.x;
	S[1][1] = scale.y;
	S[2][2] = scale.z;

	Rx[1][1] = std::cos((double)rot.x * PI / 180);
	Rx[1][2] = -std::sin((double)rot.x * PI / 180);
	Rx[2][1] = std::sin((double)rot.x * PI / 180);
	Rx[2][2] = std::cos((double)rot.x * PI / 180);

	Ry[0][0] = std::cos((double)rot.y * PI / 180);
	Ry[2][0] = -std::sin((double)rot.y * PI / 180);
	Ry[0][2] = std::sin((double)rot.y * PI / 180);
	Ry[2][2] = std::cos((double)rot.y * PI / 180);

	Rz[0][0] = std::cos((double)rot.z * PI / 180);
	Rz[1][0] = std::sin((double)rot.z * PI / 180);
	Rz[0][1] = -std::sin((double)rot.z * PI / 180);
	Rz[1][1] = std::cos((double)rot.z * PI / 180);

	Matrix result(4, 4);
	result = Rx * Ry * Rz * S * T;
	return result;
}

Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
	Matrix result = Matrix::identity(4);
	Matrix obj2world = Matrix::identity(4);
	Matrix world2cam = Matrix::identity(4);

	Vec3f z = (eye - center).normalize();
	Vec3f x = (up ^ z).normalize();
	Vec3f y = (z ^ x).normalize();

	for (int i = 0; i < 3; i++) {
		world2cam[0][i] = x[i];
		world2cam[1][i] = y[i];
		world2cam[2][i] = z[i];
		obj2world[i][3] = -center[i];
	}
	result = world2cam * obj2world;
	return result;
}

Matrix viewport(int x, int y, int w, int h) {
	Matrix result = Matrix::identity(4);
	result[0][3] = x + w / 2.f;
	result[1][3] = y + h / 2.f;
	result[2][3] = depth / 2.f;

	result[0][0] = w / 2.f;
	result[1][1] = h / 2.f;
	result[2][2] = depth / 2.f;
	return result;
}


int main(int argc, char** argv) {

	model = new Model(filename);
	TGAImage image(width, height, TGAImage::RGB);
	zBuffer = new int[width * height];

	Matrix modelViewM = lookat(camera, origin, Vec3f(0, 1.f, 0));
	Matrix viewportM = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
	Matrix projectionM = Matrix::identity(4);
	projectionM[3][2] = -1.f / (camera - origin).magnitude();

	// Be cautious about the order!
	MVP = viewportM * projectionM * modelViewM;
	std::cerr << modelViewM << std::endl;
	std::cerr << projectionM << std::endl;
	std::cerr << viewportM << std::endl;
	std::cerr << MVP << std::endl;

	drawSolid(image, zBuffer);

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("out/output.tga");
	

	// dump z-buffer (debugging purposes only)
	TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			zbimage.set(i, j, TGAColor(zBuffer[i + j * width]));
		}
	}
	zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	zbimage.write_tga_file("out/zbuffer.tga");


	delete model;
	delete[] zBuffer;
	return 0;
}