#include "..\headers\Triangulation.h"

using namespace Shapes3D;

Triangulation::Triangulation()
{
}

Triangulation::~Triangulation()
{
}

vector<Point3D> &Triangulation::uniquePoints()
{
    // Return reference to the vector of unique points
    return mUniquePoints;
}

vector<Triangle> &Triangulation::triangles()
{
    // Return reference to the vector of triangles
    return mTriangles;
}
