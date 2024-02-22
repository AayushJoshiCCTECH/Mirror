#include "../headers/Point3D.h"

Point3D::Point3D() : mX(0), mY(0), mZ(0) {}
Point3D::Point3D(double inX, double inY, double inZ) : mX(inX), mY(inY), mZ(inZ) {}
Point3D::~Point3D() {}

double Point3D::x() const
{
    return mX;
}

double Point3D::y() const
{
    return mY;
}

double Point3D::z() const
{
    return mZ;
}

void Point3D::setX(double newX)
{
    mX = newX;
}

void Point3D::setY(double newY)
{
    mY = newY;
}

void Point3D::setZ(double newZ)
{
    mZ = newZ;
}


// Overloaded less-than operator for comparing Point3D objects

bool Point3D::operator<(const Point3D &other) const
{
    if (mX < other.mX)
        return true;
    if (mX > other.mX)
        return false;
    if (mY < other.mY)
        return true;
    if (mY > other.mY)
        return false;
    return mZ < other.mZ;
}
