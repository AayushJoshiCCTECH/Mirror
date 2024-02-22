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