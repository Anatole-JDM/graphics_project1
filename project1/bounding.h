#pragma once
#include "vector.h"

class BoundingBox {
public:
    Vector Bmin;
    Vector Bmax;

    BoundingBox();
    void next(const Vector& point);
    bool intersect(const Ray& ray, double& tmin, double& tmax) const;
};