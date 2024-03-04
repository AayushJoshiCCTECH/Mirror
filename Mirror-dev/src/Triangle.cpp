#include "../headers/Triangle.h"
using namespace Shapes3D;

Triangle::Triangle() {}
Triangle::Triangle(const int inIndexFirst, const int inIndexSecond, const int inIndexThird) : mIndexOne(inIndexFirst), mIndexTwo(inIndexSecond), mIndexThree(inIndexThird) {}
Triangle::~Triangle() {}

int Triangle::indexOne() const
{
    return mIndexOne;
}

int Triangle::indexTwo() const
{
    return mIndexTwo;
}

int Triangle::indexThree() const
{
    return mIndexThree;
}

void Triangle::setindexOne(int index1)
{
    mIndexOne = index1;
}

void Triangle::setindexTwo(int index2)
{
    mIndexTwo = index2;
}

void Triangle::setindexThree(int index3)
{
    mIndexThree = index3;
}