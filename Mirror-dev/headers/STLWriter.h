// STLWriter.h
#pragma once
#include <string>
#include "Triangulation.h"

namespace Shapes3D
{
    class STLWriter
    {
    public:
        STLWriter();
        ~STLWriter();

        vector<Point3D> writeSTL(const std::string &filePath, Triangulation &triangulationObj);
    };
}
