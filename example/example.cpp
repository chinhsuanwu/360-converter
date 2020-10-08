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

    for (int i = 0; i < Converter::FACE_NUM; ++i)
    {
        int w, h, bpp;
        faces[i].img = stbi_load(argv[i + 1], &w, &h, &bpp, CHANNEL_NUM);
        faces[i].h = h;
        faces[i].w = w;
    }

    Converter::Cube cube = Converter::Cube();
    cube.setCube(faces);

    img = cube.getFace(Converter::FRONT);
    stbi_write_png("out/front.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

    img = cube.getCubeMap();
    stbi_write_png("out/cubemap.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

    Converter::Equi equi = cube.toEqui();

    img = equi.getEqui();
    stbi_write_png("out/equi.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

    Converter::Cube cube_ = equi.toCube();

    img = cube_.getFace(Converter::LEFT);
    stbi_write_png("out/left.png", img.w, img.h, CHANNEL_NUM, img.img, img.w * CHANNEL_NUM);

    return 0;
}