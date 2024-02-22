# 3D Mirror Image Generator

This project provides a tool to generate a mirrored version of a 3D object represented in an STL file, with respect to a user-defined plane.

## Overview

The program reads an input STL file containing the geometry of a 3D object. It prompts the user to define a plane by specifying three points. Using the provided plane, the program computes the mirror image of the object and saves it to a new STL file.

## Features

- **STL File Support**: Reads geometry data from STL files.
- **Interactive Interface**: User-friendly interface for defining reflection planes.
- **Mirror Image Generation**: Computes mirror image of the object with respect to the defined plane.
- **Output File**: Saves the mirrored geometry to a new STL file.

## Usage

1. **Compile**: Compile the project using your preferred C++ compiler.
2. **Run**: Execute the compiled binary.
3. **Input**: Provide the path to the input STL file.
4. **Define Plane**: Define a plane by specifying three points.
5. **Generate Mirror Image**: The program computes the mirror image of the object and saves it to a new STL file.
6. **View Results**: Open the generated STL file in a 3D viewer to inspect the mirrored object.

## Dependencies

- C++ Compiler
- STL file reader (included)
- Other dependencies if any (e.g., GUI libraries for interactive interface)

## Example Usage

```bash
./mirror_generator input.stl
```
## Mirror Class

```c++
//Mirror.h
#pragma once


class Point3D
{
public:
    Point3D();
    Point3D(double inX, double inY, double inZ);
    ~Point3D();

    double x() const;
    double y() const;
    double z() const;

    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    bool operator<(const Point3D &other) const;
    
private:
    double mX;
    double mY;
    double mZ;
};
```
## Mirror.cpp

```c++
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include "..\headers\Point3D.h"
#include "..\headers\mirror.h"
using namespace std;
using namespace Shapes3D;
 mirror::mirror()
    {
    }

    mirror::~mirror()
    {
    }

std::vector<Point3D> mirror::getPlane(std::vector<Point3D> readvector)
{
    double vx=1, vy=1, vz=1, px=1, py=1, pz=1;
    Point3D Plane_Point3Ds(px, py, pz);
    Point3D Vector_Point3Ds(vx, vy, vz);
    readvector.insert(readvector.begin(), Vector_Point3Ds);
    readvector.insert(readvector.begin() + 1, Plane_Point3Ds);
    return readvector;
}

std::vector<Point3D> mirror::reflectPoint(std::vector<Point3D> Point3D_vector)
{
    const double a = Point3D_vector[0].x();
    const double b = Point3D_vector[0].y();
    const double c = Point3D_vector[0].z();
    const double div = sqrt(pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0));
    Point3D normalised_vector(a / div, b / div, c / div);
    double x1 = Point3D_vector[1].x(), y1 = Point3D_vector[1].y(), z1 = Point3D_vector[1].z();
    std::vector<Point3D> reflected_shape;
    for (int i = 2; i < Point3D_vector.size(); i++)
    {
        double tx = Point3D_vector[i].x();

        double ty = Point3D_vector[i].y();

        double tz = Point3D_vector[i].z();

        double dist = abs(((a * tx) + (b * ty) + (c * tz) - (a * x1) - (b * y1) - (c * z1)) / (div));

        reflected_shape.push_back(Point3D(Point3D_vector[i].x() + (2 * dist * normalised_vector.x()), Point3D_vector[i].y() + (2 * dist * normalised_vector.y()), Point3D_vector[i].z() + (2 * dist * normalised_vector.z())));

    }

    return reflected_shape;

}

void mirror::plotPlane(std::vector<Point3D> reflect)
{
    std::ofstream out("textFiles/plane.txt");
    for (auto i : reflect)
    {
        out << i.x() << " " << i.y() << " " << i.z() << std::endl;
    }
    out.close();
}
std::vector<Point3D> mirror::writePlane(std::vector<Point3D> Point3D_vector)
{
    const double a = Point3D_vector[0].x();
    const double b = Point3D_vector[0].y();
    const double c = Point3D_vector[0].z();
    const double div = sqrt(pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0));
    Point3D normalised_vector(a / div, b / div, c / div);
    double x1 = Point3D_vector[1].x(), y1 = Point3D_vector[1].y(), z1 = Point3D_vector[1].z();
    std::vector<Point3D> planar_shape;
    for (int i = 2; i < Point3D_vector.size(); i++)
    {
        double tx = Point3D_vector[i].x();
        double ty = Point3D_vector[i].y();
        double tz = Point3D_vector[i].z();
        double dist = abs(((a * tx) + (b * ty) + (c * tz) - (a * x1) - (b * y1) - (c * z1)) / (div));
        planar_shape.push_back(Point3D(Point3D_vector[i].x() + (dist * normalised_vector.x()), Point3D_vector[i].y() + (dist * normalised_vector.y()), Point3D_vector[i].z() + (dist * normalised_vector.z())));
    }
    plotPlane(planar_shape);
    return planar_shape;
}

void mirror::plot(std::vector<Point3D> reflect)
{
    std::ofstream out("textFiles/reflection.txt");
    for (auto i : reflect)
    {
        out << i.x() << " " << i.y() << " " << i.z() << std::endl;
    }
    out.close();
}
```



