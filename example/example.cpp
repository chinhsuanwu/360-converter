#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "include/stb_image_write.h"

#define CONVERTER_HPP
#include "../src/converter.hpp"

#ifndef CHANNEL_NUM
#define CHANNEL_NUM 4
#endif


int main(int argc, char** argv) {
    Converter::Image faces[Converter::FACE_NUM];

    for (int i = 0; i < Converter::FACE_NUM; ++i) {
        int w, h, bpp;
        faces[i].img = stbi_load(argv[i+1], &w, &h, &bpp, CHANNEL_NUM);
        faces[i].h = h;
        faces[i].w = w;
    }

    Converter::Cube cube = Converter::Cube();
    cube.setCube(faces);

    Converter::Image face = cube.getFace(5);
    stbi_write_png("out/face.png", face.w, face.h, CHANNEL_NUM, face.img, face.w*CHANNEL_NUM);
    
    Converter::Image cubemap = cube.getCubeMap();
    stbi_write_png("out/cubemap.png", cubemap.w, cubemap.h, CHANNEL_NUM, cubemap.img, cubemap.w*CHANNEL_NUM);

    Converter::Equi equi = cube.toEqui();

    Converter::Image equi_ = equi.getEqui();
    stbi_write_png("out/equi.png", equi_.w, equi_.h, CHANNEL_NUM, equi_.img, equi_.w*CHANNEL_NUM);

    return 0;
}