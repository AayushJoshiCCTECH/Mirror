#include <iostream>
#include "./headers/Triangulation.h"
#include "./headers/Triangle.h"
#include "./headers/STLReader.h"
#include "./headers/STLWriter.h"
#include "./headers/mirror.h"

using namespace Shapes3D;
// Define a typedef for a vector of Triangle objects
typedef vector<Triangle> Triangles;

int main()
{
    std::string inputFilePath = "stlFiles/cube1.stl";
    std::string outputFilePath = "textFiles/original.txt";  

    Shapes3D::STLReader stlReader;
    Shapes3D::STLWriter stlWriter;
    Shapes3D::mirror m;


    Triangulation triangulation = stlReader.readSTL(inputFilePath);
    vector<Point3D> one=stlWriter.writeSTL(outputFilePath, triangulation);
    
    vector<Point3D> second = m.getPlane(one);
    m.writePlane(second);
    vector<Point3D> third = m.reflectPoint(second);
    m.plot(third);
    return 0;
}
