#include "scene.h"
#include <cmath>
#include <algorithm>
#include <random>


static thread_local std::mt19937 rng(std::random_device{}());
static thread_local std::uniform_real_distribution<double> uni01(0.0, 1.0);

Vector random_cos(const Vector& N) {
    double r1 = uni01(rng);
    double r2 = uni01(rng);
    double z = sqrt(r2);
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    
    Vector T1;
    Vector T2;

    if (fabs(N[0]) < fabs(N[1]) && fabs(N[0]) < fabs(N[2])) { T1 = Vector(0, -N[2], N[1]); } 
    else if (fabs(N[1]) < fabs(N[2]))                       { T1 = Vector(-N[2], 0, N[0]); } 
    else                                                    { T1 = Vector(-N[1], N[0], 0); }

    T1.normalize();
    T2 = cross(N, T1);

    return x*T1 + y*T2 + z*N;
}

Scene::Scene(const Vector& lightPos, double intensity) : lightPosition(lightPos), lightIntensity(intensity) {}

void Scene::add(Geometry* object) { objects.push_back(object); }


bool Scene::findClosestIntersection(const Ray& ray, const Geometry*& target, Vector& P, Vector& N, double& intersection) const {
    intersection = 999999999999;
    target = nullptr;

    for (const auto& object : objects) {
        double t;
        Vector P2, N2;
        if (object->intersect(ray, t, P2, N2) && t < intersection && t > 0) {
            intersection = t;
            target = object;
            P = P2;
            N = N2;}}

    return !(target == nullptr);
}

double Scene::visibility(const Vector& P, const Vector& N, const Vector& S) const {
    Vector P2 = P + 1e-9 * N;
    Vector oi = S - P2;
    double d = oi.norm();
    oi.normalize();
    Ray RAY_OF_DEATH(P2, oi);

    for (const auto& object : objects) {
        double t;
        Vector tempP, tempN;
        if (object->intersect(RAY_OF_DEATH, t, tempP, tempN) && t < d) {
            return 0.0;
        }
    }

    return 1.0;
}

Vector Scene::refraction(const Vector& I, const Vector& N, double n1, double n2, bool& total) const {
    double eta = n1 / n2;
    double cosi = dot(I, N);
    double sin2t = square(eta) * (1.0 - square(cosi));
    if (sin2t > 1.0) { 
        total = true;
        return Vector(0,0,0);}
    else {

        double cost = sqrt(1.0 - sin2t);
        Vector T2 = eta * (I - cosi * N);
        Vector nc = -cost * N;

        Vector T = T2 + nc;
        T.normalize();
        total = false;
        return T;
    }
}


Vector Scene::getColor(const Ray& ray, int depth) const {
    if (depth < 0) {return Vector(0, 0, 0);}

    const Geometry* target;
    Vector P, N;
    double thit;

    if (!findClosestIntersection(ray, target, P, N, thit)) {return Vector(0, 0, 0);}

    Vector L0(0, 0, 0);

    if (target->material == TRANSPARENT) {
        double n1, n2;
        Vector NN;
        if (dot(ray.direction, N) < 0) {
            n1 = target->refractiveIndex;
            n2 = target->insideRefractiveIndex;
            NN = N;
        } else {
            n1 = target->insideRefractiveIndex;
            n2 = target->refractiveIndex;
            NN = -N;
        }

        double k0 = square(n1 - n2) / (square(n1 + n2));
        double R2 = k0 + (1 - k0) * pow(1 - fabs(dot(ray.direction, NN)), 5);

        double u = uni01(rng);
        if (u < R2) {
            Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
            R.normalize();
            return getColor(Ray(P + 1e-4 * R, R), depth - 1);
        } else {
            bool total = false;
            Vector T = refraction(ray.direction, NN, n1, n2, total);
            if (total) {
                Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
                R.normalize();
                return getColor(Ray(P + 1e-4 * R, R), depth - 1);
            }
            return getColor(Ray(P + 1e-4 * T, T), depth - 1);
        }
    } else if (target->material == NORMAL) {
        Vector L = lightPosition - P;
        double dist2 = L.norm2();
        L.normalize();
    
        L0 += (target->albedo / M_PI) * (lightIntensity / (4 * M_PI * dist2)) * visibility(P, N, lightPosition) * std::max(dot(N, L), 0.0);
    
        if (depth > 0) {
            Vector randomDir = random_cos(N);
            Ray randomRay(P + 1e-4 * randomDir, randomDir);
            L0 += target->albedo * getColor(randomRay, depth - 1);
        }
    
        return L0;

    } else if (target->material == MIRROR) {
        Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
        R.normalize();
        return getColor(Ray(P + 1e-4 * R, R), depth - 1);
    }

    return Vector(0, 0, 0);

}