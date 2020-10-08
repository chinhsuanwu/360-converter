#ifndef CONVERTER_HPP
#define CONVERTER_HPP

#include <assert.h>
#include <stdint.h>
#include <math.h>
#endif

#ifndef CHANNEL_NUM
#define CHANNEL_NUM 4
#endif

namespace Converter
{
    enum Face
    {
        FRONT = 0,
        RIGHT,
        BACK,
        LEFT,
        TOP,
        DOWN,
        FACE_NUM
    };

    struct Coord
    {
        Face face; // the face of the cube
        double x;  // the x Coordinate
        double y;  // the y Coordinate
    };

    struct Image
    {
        unsigned int h;
        unsigned int w;
        uint8_t *img;

        /**
         * Suppose to put in GPU using CUDA
         */
        void copyTo(Image dst, unsigned int start_h, unsigned int start_w)
        {
            for (unsigned int i = start_h; i < start_h + h; ++i)
                for (unsigned int j = start_w; j < start_w + w; ++j)
                    for (int k = 0; k < CHANNEL_NUM; ++k)
                        dst.img[CHANNEL_NUM * (i * dst.w + j) + k] = img[CHANNEL_NUM * ((i - start_h) * w + (j - start_w)) + k];
        }
    };

    class Cube;
    class Equi;
    class Stereo;
    class Cube
    {
    public:
        Cube(){};

        void setCube(Image t_faces[FACE_NUM])
        {
            for (int i = 0; i < FACE_NUM; ++i)
            {
                assert(t_faces[i].h == t_faces[i].w);
                faces[i] = t_faces[i];
            }
            cube_h = faces[0].h;
            cube_w = faces[0].w;
        }

        Image getFace(int i)
        {
            return faces[i];
        }

        Image getCubeMap()
        {
            cubemap = {cube_h * 3, cube_w * 4};
            cubemap.img = new uint8_t[cubemap.w * cubemap.h * CHANNEL_NUM];
            faces[FRONT].copyTo(cubemap, cube_h, cube_w);     // front
            faces[RIGHT].copyTo(cubemap, cube_h, 2 * cube_w); // right
            faces[BACK].copyTo(cubemap, cube_h, 3 * cube_w);  // left
            faces[LEFT].copyTo(cubemap, cube_h, 0);           // back
            faces[TOP].copyTo(cubemap, 0, cube_w);            // top
            faces[DOWN].copyTo(cubemap, 2 * cube_h, cube_w);  // down

            return cubemap;
        }

        inline bool isInRange(const double l, const double theta, const double u)
        {
            return (theta >= l && theta < u);
        }

        const Face getFaceID(const double theta, const double phi)
        {
            Face faceID;
            double n_theta; // M_PI_4 / (double) cube_h

            if (isInRange(-M_PI_4, theta, M_PI_4))
            {
                faceID = FRONT;
                n_theta = theta;
            }
            else if (isInRange(M_PI_4, theta, M_PI_2 + M_PI_4))
            {
                faceID = RIGHT;
                n_theta = theta - M_PI_2;
            }
            else if (isInRange(-(M_PI_2 + M_PI_4), theta, -M_PI_4))
            {
                faceID = LEFT;
                n_theta = theta + M_PI_2;
            }
            else
            {
                faceID = BACK;
                if (theta > 0.0)
                {
                    n_theta = theta - M_PI;
                }
                else
                {
                    n_theta = theta + M_PI;
                }
            }

            double phiThr = atan2(1.0, 1.0 / cos(n_theta));
            if (phi > phiThr)
            {
                faceID = DOWN;
            }
            else if (phi < -phiThr)
            {
                faceID = TOP;
            }
            return faceID;
        }

        void transform(const double axis, const double px, const double py, const double radian, double *coord_x, double *coord_y)
        {
            double radius = (double)cube_h / 2.0;
            double ratio = radius / axis;

            *coord_x = ratio * px;
            *coord_y = ratio * py;

            // Rotation
            double tmp = *coord_x;
            *coord_x = *coord_x * cos(radian) - *coord_y * sin(radian);
            *coord_y = tmp * sin(radian) + *coord_y * cos(radian);

            // Translation
            *coord_x += radius;
            *coord_y += radius;
        }

        void setCEMap()
        {
            equi_h = cube_h * 2;
            equi_w = cube_w * 4;

            CE_map = new Coord[equi_h * equi_w];

            unsigned int pix = 0;
            for (int i = 0; i < equi_h; ++i)
            {
                for (int j = 0; j < equi_w; ++j)
                {
                    double u = ((2.0 * j) / equi_w - 1.0);
                    double v = ((2.0 * i) / equi_h - 1.0);
                    double theta = u * M_PI, phi = v * M_PI_2;
                    double x = cos(phi) * cos(theta), y = sin(phi), z = cos(phi) * sin(theta);

                    Face faceID = getFaceID(theta, phi);

                    double coord_x, coord_y;
                    switch (faceID)
                    {
                    case FRONT:
                        transform(x, z, y, 0.0, &coord_x, &coord_y);
                        break;
                    case RIGHT:
                        transform(z, y, x, M_PI_2, &coord_x, &coord_y);
                        break;
                    case BACK:
                        transform(x, y, z, -M_PI_2, &coord_x, &coord_y);
                        break;
                    case LEFT:
                        transform(z, x, y, M_PI, &coord_x, &coord_y);
                        break;
                    case TOP:
                        transform(y, z, x, M_PI, &coord_x, &coord_y);
                        break;
                    case DOWN:
                        transform(y, x, z, -M_PI_2, &coord_x, &coord_y);
                        break;
                    default:
                        break;
                    }

                    CE_map[pix].face = faceID;
                    CE_map[pix].x = coord_x;
                    CE_map[pix++].y = coord_y;
                }
            }
        }

        Coord getCECoord(int i, int j)
        {
            return CE_map[i * equi_w + j];
        }

        Coord *getCEMap()
        {
            return CE_map;
        }

        Equi toEqui();

    private:
        // input size
        unsigned int cube_w, cube_h;
        unsigned int equi_w, equi_h;

        // assets
        Image faces[FACE_NUM];
        Image cubemap;
        Coord *CE_map;
    };

    class Equi
    {
    public:
        Equi(){};

        void setEqui(Image t_equi)
        {
            equi = t_equi;
            equi_h = equi.h;
            equi_w = equi.w;
        }

        void setCEMap(Coord *t_CE_map)
        {
            CE_map = t_CE_map;
        }

        Image getEqui()
        {
            return equi;
        }

        Cube toCube();

    private:
        // input size
        unsigned int equi_w, equi_h;

        // assets
        Image equi;
        Coord *CE_map;
    };

    class Stereo
    {
    public:
        Stereo(){};

    private:
        // input size
        unsigned int stereo_w, stereo_h;

        // assets
        Image stereo;
    };

    Equi Cube::toEqui()
    {
        Cube::setCEMap();
        Image equi_;
        equi_.h = equi_h;
        equi_.w = equi_w;
        equi_.img = new uint8_t[equi_w * equi_h * CHANNEL_NUM];

        for (int i = 0; i < equi_h; ++i)
        {
            for (int j = 0; j < equi_w; ++j)
            {
                Coord coord = getCECoord(i, j);
                for (int k = 0; k < CHANNEL_NUM; ++k)
                {
                    equi_.img[CHANNEL_NUM * (i * equi_.w + j) + k] =
                        faces[coord.face].img[CHANNEL_NUM * ((unsigned long)coord.y * cube_w + (unsigned long)coord.x) + k];
                }
            }
        }
        Equi equi = Equi();
        equi.setEqui(equi_);
        equi.setCEMap(Cube::getCEMap());
        return equi;
    }

    Cube Equi::toCube()
    {
        Cube cube = Cube();
        return cube;
    }
} // namespace Converter