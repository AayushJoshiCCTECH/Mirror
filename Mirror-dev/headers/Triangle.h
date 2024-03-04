#pragma once
#include "Point3D.h"
namespace Shapes3D
{
    // Class representing a triangle with three indices
    class Triangle
    {
    public:
        Triangle();
        Triangle(const int inIndexFirst, const int inIndexsecond, const int inIndexthird);
        ~Triangle();

        int indexOne() const;
        int indexTwo() const;
        int indexThree() const;

        void setindexOne(int indexFirst);
        void setindexTwo(int indexSecond);
        void setindexThree(int indexThird);

    private:
        int mIndexOne;
        int mIndexTwo;
        int mIndexThree;
    };
}
