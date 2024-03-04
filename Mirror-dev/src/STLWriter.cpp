// STLWriter.cpp
#include <iostream>
#include <fstream>
#include "../headers/STLWriter.h"

using namespace std;
using namespace Shapes3D;

STLWriter::STLWriter() {}
STLWriter::~STLWriter() {}

vector<Point3D> STLWriter::writeSTL(const std::string &filePath, Triangulation &triangulationObj)
{
    ofstream outFile(filePath);
    vector<Point3D> mypoints;
    if (!outFile.is_open())
    {
        cout << "Error while opening the file for writing." << endl;
    }

    // Get references to the unique points in the triangulation object
    vector<Point3D> &points = triangulationObj.uniquePoints();

    // Write the triangulation data to the file
    for (const Point3D &point : points)
    {
        //point.x(),point.y()point.z()
        mypoints.push_back(point);
        outFile << point.x() << " " << point.y() << " " << point.z() << endl;
    }
    return mypoints;
}
