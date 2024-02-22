
#pragma once

#include <string>

#include <vector>

#include "Triangulation.h"

namespace Shapes3D

{

    // Class for reading triangulation data from STL files
    class STLReader
    {
    public:
        STLReader();
        ~STLReader();
         Triangulation readSTL(const std::string &filePath);
    };

}