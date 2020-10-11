#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "include/stb_image_write.h"

#include "../src/converter.hpp"

#ifndef CHANNEL_NUM
#define CHANNEL_NUM 4
#endif

int main(int argc, char **argv)
{
	Converter::Image img;
	Converter::Image faces[Converter::FACE_NUM];

	int w, h, bpp;
	for (int i = 0; i < Converter::FACE_NUM; ++i)
	{
		faces[i].img = stbi_load(argv[i + 1], &w, &h, &bpp, CHANNEL_NUM);
		faces[i].h = h, faces[i].w = w;
	}

	Converter::Face face = Converter::Face(faces);

	img = face.getFace(Converter::FRONT);
	stbi_write_png("out/front.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

	Converter::Cube cube = face.toCube();
	img = cube.getCubeMap();
	stbi_write_png("out/cubemap.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

	Converter::Equi equi = cube.toEqui();
	img = equi.getEqui();
	stbi_write_png("out/equi.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

	img = equi.toStereo(Converter::TOP).getStereo();
	stbi_write_png("out/stereo.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

	return 0;
}