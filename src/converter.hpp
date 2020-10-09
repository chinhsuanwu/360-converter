#ifndef CONVERTER_HPP
#define CONVERTER_HPP

#include <assert.h>
#include <stdint.h>
#include <math.h>

#ifndef CHANNEL_NUM
#define CHANNEL_NUM 4
#endif

namespace Converter
{
    enum FaceID
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
        FaceID face; // the face of the face
        double x;    // the x Coordinate
        double y;    // the y Coordinate
    };

    struct Image
    {
        unsigned int h;
        unsigned int w;
        uint8_t *img;
    };

    class Converter
    {
    public:
        Converter()
        {
            FE_map = FS_map = NULL;
        }

        inline bool isInRange(const double l, const double theta, const double u)
        {
            return (theta >= l && theta < u);
        }

        const FaceID getFaceID(const double theta, const double phi)
        {
            FaceID faceID;
            double n_theta;

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
                n_theta = (theta > 0.0 ? theta - M_PI : theta + M_PI);
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
            double radius = (double)face_h / 2.0;
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

        void setFEMap(Coord *t_FE_map = NULL)
        {
            if (t_FE_map)
            {
                FE_map = t_FE_map;
                return;
            }

            assert(face_h == face_w);
            equi_h = face_h * 2;
            equi_w = face_w * 4;

            FE_map = new Coord[equi_h * equi_w];

            unsigned int pix = 0;
            for (int i = 0; i < equi_h; ++i)
            {
                for (int j = 0; j < equi_w; ++j)
                {
                    double u = ((2.0 * j) / equi_w - 1.0);
                    double v = ((2.0 * i) / equi_h - 1.0);
                    double theta = u * M_PI, phi = v * M_PI_2;
                    double x = cos(phi) * cos(theta), y = sin(phi), z = cos(phi) * sin(theta);

                    FaceID faceID = getFaceID(theta, phi);

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

                    FE_map[pix].face = faceID;
                    FE_map[pix].x = coord_x;
                    FE_map[pix++].y = coord_y;
                }
            }
        }

        Coord getFECoord(int i, int j)
        {
            return FE_map[i * equi_w + j];
        }

        Coord *getFEMap()
        {
            return FE_map;
        }

    protected:
        unsigned int face_h, face_w;
        unsigned int cube_h, cube_w;
        unsigned int equi_h, equi_w;
        unsigned int stereo_h, stereo_w;
        Coord *FE_map, *FS_map;
    };

    class Face;
    class Cube;
    class Equi;
    class Stereo;

    class Face : public Converter
    {
    public:
        Face(Image t_faces[FACE_NUM])
        {
            for (int i = 0; i < FACE_NUM; ++i)
            {
                assert(t_faces[i].h == t_faces[i].w);
                faces[i] = t_faces[i];
            }

            face_h = faces[0].h;
            face_w = faces[0].w;
        };

        Image getFace(int faceID)
        {
            return faces[faceID];
        }

        void copy(Image src, Image dst, unsigned int start_h, unsigned int start_w)
        {
            for (unsigned int i = start_h; i < start_h + src.h; ++i)
                for (unsigned int j = start_w; j < start_w + src.w; ++j)
                    for (int k = 0; k < CHANNEL_NUM; ++k)
                        dst.img[CHANNEL_NUM * (i * dst.w + j) + k] = src.img[CHANNEL_NUM * ((i - start_h) * src.w + (j - start_w)) + k];
        }

        Cube toCube();
        Equi toEqui();

    private:
        Image faces[FACE_NUM];
    };

    class Cube : public Converter
    {
    public:
        Cube(Image t_cubemap) : cubemap(t_cubemap)
        {
            cube_h = cubemap.h;
            cube_w = cubemap.w;
        };

        Image getCubeMap()
        {
            return cubemap;
        }

        void crop(Image src, Image dst, unsigned int start_h, unsigned int start_w)
        {
            for (unsigned int i = start_h; i < start_h + dst.h; ++i)
                for (unsigned int j = start_w; j < start_w + dst.w; ++j)
                    for (int k = 0; k < CHANNEL_NUM; ++k)
                        dst.img[CHANNEL_NUM * ((i - start_h) * dst.w + (j - start_w)) + k] = src.img[CHANNEL_NUM * (i * src.w + j) + k];
        }

        Face toFace();
        Equi toEqui();

    private:
        Image cubemap;
    };

    class Equi : public Converter
    {
    public:
        Equi(Image t_equi) : equi(t_equi)
        {
            equi_h = equi.h;
            equi_w = equi.w;
        };

        Image getEqui()
        {
            return equi;
        }

        Face toFace();
        Cube toCube();

    private:
        Image equi;
    };

    class Stereo : public Converter
    {
    public:
        Stereo(){};

    private:
        Image stereo;
    };

    Cube Face::toCube()
    {
        Image img = {face_h * 3, face_w * 4, new uint8_t[img.w * img.h * CHANNEL_NUM]};
        copy(faces[FRONT], img, face_h, face_w);
        copy(faces[RIGHT], img, face_h, 2 * face_w);
        copy(faces[BACK], img, face_h, 3 * face_w);
        copy(faces[LEFT], img, face_h, 0);
        copy(faces[TOP], img, 0, face_w);
        copy(faces[DOWN], img, 2 * face_h, face_w);

        return Cube(img);
    }

    Equi Face::toEqui()
    {
        if (!FE_map)
            setFEMap();

        equi_h = face_h * 2;
        equi_w = face_w * 4;

        Image img = {equi_h, equi_w, new uint8_t[equi_w * equi_h * CHANNEL_NUM]};

        for (int i = 0; i < equi_h; ++i)
        {
            for (int j = 0; j < equi_w; ++j)
            {
                Coord coord = getFECoord(i, j);
                for (int k = 0; k < CHANNEL_NUM; ++k)
                {
                    img.img[CHANNEL_NUM * (i * img.w + j) + k] =
                        faces[coord.face].img[CHANNEL_NUM * ((unsigned long)coord.y * face_w + (unsigned long)coord.x) + k];
                }
            }
        }

        Equi equi = Equi(img);
        equi.setFEMap(FE_map);

        return equi;
    }

    Face Cube::toFace()
    {
        face_h = cubemap.h / 3;
        face_w = cubemap.w / 4;

        Image faces[FACE_NUM];
        for (int i = 0; i < FACE_NUM; ++i) 
            faces[i] = (Image) {face_h, face_w, new uint8_t[face_w * face_h * CHANNEL_NUM]};

        crop(cubemap, faces[FRONT], face_h, face_w);
        crop(cubemap, faces[RIGHT], face_h, 2 * face_w);
        crop(cubemap, faces[BACK], face_h, 3 * face_w);
        crop(cubemap, faces[LEFT], face_h, 0);
        crop(cubemap, faces[TOP], 0, face_w);
        crop(cubemap, faces[DOWN], 2 * face_h, face_w);

        return Face(faces);
    }

    Equi Cube::toEqui()
    {
        return toFace().toEqui();
    }

    Face Equi::toFace()
    {
        if (!FE_map)
            setFEMap();

        face_h = equi.h / 2;
        face_w = equi.w / 4;

        Image faces[FACE_NUM];
        for (int i = 0; i < FACE_NUM; ++i)
            faces[i] = (Image) {face_h, face_w, new uint8_t[face_w * face_h * CHANNEL_NUM]};

        for (int i = 0; i < equi.h; ++i)
        {
            for (int j = 0; j < equi.w; ++j)
            {
                Coord coord = getFECoord(i, j);
                for (int k = 0; k < CHANNEL_NUM; ++k)
                {
                    faces[coord.face].img[CHANNEL_NUM * ((unsigned long)coord.y * face_w + (unsigned long)coord.x) + k] =
                        equi.img[CHANNEL_NUM * (i * equi.w + j) + k];
                }
            }
        }

        Face face = Face(faces);
        face.setFEMap(FE_map);

        return face;
    }

    Cube Equi::toCube()
    {
        return toFace().toCube();
    }

} // namespace Converter
#endif