# 360-converter

## Overview
Input one of the following type 360 image and convert it to another one.
| Faces | Cubemap | Equirectangular | Stereographic |
| :---: | :---: | :---:| :---:|
| ![](./assets/faces.png) | ![](./assets/cubemap.png) | ![](./assets/equi.png) | ![](./assets/stereo.png)|

## Prerequisites
This project has been tested on:
- Ubuntu 18.04
- G++ 7.5.0
- Clang 6.0.0

The only thing you need is a C/C++ compiler, no other package do you need to install.

## Run example
```
git clone https://github.com/chinhsuanwu/360-converter.git
cd 360-converter
make
make run
```

## Usage
Just directly include the [`converter.hpp`](https://github.com/chinhsuanwu/360-converter/blob/master/src/converter.hpp) in folder `src`
```c++
#include "src/converter.hpp"
```
For example, the conversion from `Face` to `Cube`, `Equi` and `Stereo` can be shown as below, other conversions can be deduced:
![](./assets/usage.png)

- To convert from `Face` to `Equi`
    ```c++
    Converter::Equi equi = face.toEqui();
    ```
- To get image
    ```c++
    Converter::Image img = equi.getEqui();
    ```
- For `Stereo` images, you can assign `TOP` or `DOWN`, while the default will be `DOWN`
    ```c++
    Converter::Image img = cube.toStereo(Converter::TOP).getStereo();
    ```
    | `TOP` | `DOWN` |
    | :---: | :---: |
    | ![](assets/stereo_top.png) | ![](assets/stereo.png) |
- To write image
    ```c++
    stbi_write_png("out/equi.png", img.w, img.h, CHANNEL_NUM, img.img, img.w*CHANNEL_NUM);
    ```

Find out more at [example.cpp](https://github.com/chinhsuanwu/360-converter/blob/master/example/example.cpp)

## Acknowledgement
This project is using [stb](https://github.com/nothings/stb) library for image I/O, great thanks to their excellent work.