#pragma once
#include "vector.h"
#include "trianglemesh.h"
#include <vector>

class Scene {
    public:
        Scene(const Vector& lightPos, double intensity);

        Vector getColor(const Ray& ray, int depth) const;
        bool findClosestIntersection(const Ray& ray, const Geometry*& hit_object, Vector& P, Vector& N, double& t_hit) const;
        double visibility(const Vector& P, const Vector& N, const Vector& S) const;
        Vector refraction(const Vector& I, const Vector& N, double n1, double n2, bool& tir) const;
        void add(Geometry* object);

        std::vector<Geometry*> objects;
        Vector lightPosition;
        double lightIntensity;
};