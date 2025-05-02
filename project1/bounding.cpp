#include "bounding.h"
#include "trianglemesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "vector.h"


BoundingBox::BoundingBox() : Bmin(Vector(1e30, 1e30, 1e30)), Bmax(Vector(-1e30, -1e30, -1e30)) {}

void BoundingBox::next(const Vector& point) {
    for (int i = 0; i < 3; i++) {
        Bmin[i] = std::min(Bmin[i], point[i]);
        Bmax[i] = std::max(Bmax[i], point[i]);
    }
}

bool BoundingBox::intersect(const Ray& ray, double& tmin, double& tmax) const {
    tmin = -1e30;
    tmax = 1e30;

    for (int i = 0; i < 3; i++) {
        double t0 = (Bmin[i] - ray.origin[i]) / ray.direction[i];
        double t1 = (Bmax[i] - ray.origin[i]) / ray.direction[i];
        if (t0 > t1) {
            std::swap(t0, t1);
        }
        tmin = std::max(tmin, t0);
        tmax = std::min(tmax, t1);

        if (tmin > tmax) {
            return false;
        }
    }

    return true;
}

void TriangleMesh::computeBoundingBox() {
    box = BoundingBox();
    for (const auto& vertex : vertices) {
        box.next(vertex);
    }
}
