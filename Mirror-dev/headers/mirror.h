#pragma once
#include "Point3D.h"
#include "Triangulation.h"
#include <vector>

namespace Shapes3D
{
    class mirror
    {
    public:
        mirror();
        ~mirror();

        std::vector<Point3D> readFile();
        std::vector<Point3D> getPlane(std::vector<Point3D> readvector);
        std::vector<Point3D> reflectPoint(std::vector<Point3D> point_vector);
        std::vector<Point3D> writePlane(std::vector<Point3D> point_vector);
        void plotPlane(std::vector<Point3D> reflect);
        void plot(std::vector<Point3D> reflect);
    };
}
