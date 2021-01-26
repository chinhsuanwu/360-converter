# 360-converter

## Overview
Input one of the following type 360 image and convert it to another one.
| Faces | Cubemap | Equirectangular | Stereographic |
| :---: | :---: | :---:| :---:|
| ![](https://user-images.githubusercontent.com/67839539/105892418-cdb67580-604c-11eb-95b9-877cafa2dfca.png) | ![](https://user-images.githubusercontent.com/67839539/105892517-efaff800-604c-11eb-8963-0a4ea77752ac.png) | ![](https://user-images.githubusercontent.com/67839539/105892651-1c640f80-604d-11eb-9d1d-d01e9ff5b02b.png) | ![](https://user-images.githubusercontent.com/67839539/105892704-2f76df80-604d-11eb-99f2-0f9dd3a20a22.png) |

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
Just directly include the [`converter.hpp`](https://github.com/chinhsuanwu/360-converter/blob/master/src/converter.hpp) in folder `src`.
```c++
#include "src/converter.hpp"
```
For example, the conversion from `Face` to `Cube`, `Equi` and `Stereo` can be shown as below, other conversions can be deduced:
 <p align="center"><img src="https://user-images.githubusercontent.com/67839539/105900200-ad8bb400-6056-11eb-834b-cadf902cab89.png" height="70%" width="70%"></p>

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
    | <p align="center"><img src="https://user-images.githubusercontent.com/67839539/105893800-7f09db00-604e-11eb-9787-308b4376b133.png" height="70%" width="70%"></p> | <p align="center"><img src="https://user-images.githubusercontent.com/67839539/105892704-2f76df80-604d-11eb-99f2-0f9dd3a20a22.png" height="70%" width="70%"></p> |

- To load image
    ```c++
    img.img = stbi_load("assets/equi.png", &w, &h, &bpp, CHANNEL_NUM);
	img.w = w, img.h = h;
    ```
- To write image
    ```c++
    stbi_write_png("out/equi.png", img.w, img.h, CHANNEL_NUM, img.img, img.w*CHANNEL_NUM);
    ```

Find out more at [example.cpp](https://github.com/chinhsuanwu/360-converter/blob/master/example/example.cpp).

## Acknowledgement
This project is using [stb](https://github.com/nothings/stb) library for image I/O in [example.cpp](https://github.com/chinhsuanwu/360-converter/blob/master/example/example.cpp), great thanks to their excellent work.