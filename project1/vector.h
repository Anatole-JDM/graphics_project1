#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double norm2() const;
    double norm() const;
    void normalize();
    double operator[](int i) const;
    double& operator[](int i);
    Vector operator-() const;
    double data[3];
    Vector& operator+=(const Vector& other);
    friend Vector operator*(const Vector& a, const Vector& b); 
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const double a, const Vector& b);
Vector operator*(const Vector& a, const double b);
Vector operator/(const Vector& a, const double b);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

class Ray {
public:
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {}
    Vector origin;
    Vector direction;
}; // taken live during the lecture 1

enum MaterialType {NORMAL, MIRROR, TRANSPARENT};

class Geometry {
    public:
        virtual bool intersect(const Ray& ray, double& t, Vector& P, Vector& N) const = 0;
    
        Vector albedo;
        MaterialType material;
        double refractiveIndex;      
        double insideRefractiveIndex;
        bool invertNormals;           

        Geometry(const Vector& albedo, MaterialType material = NORMAL, double refractiveIndex = 1.0, double insideRefractiveIndex = 1.0, bool invertNormals = false)
            : albedo(albedo), material(material), refractiveIndex(refractiveIndex), insideRefractiveIndex(insideRefractiveIndex), invertNormals(invertNormals) {}
    };


class Sphere : public Geometry {
    public:
        Sphere(const Vector& center, double radius, const Vector& albedo, MaterialType material = NORMAL,
                double refractiveIndex = 1.0, double insideRefractiveIndex = 1.0, bool invertNormals = false, double kr = 0.0)
            : Geometry(albedo, material, refractiveIndex, insideRefractiveIndex, invertNormals), 
                center(center), radius(radius), n1(refractiveIndex), n2(insideRefractiveIndex), kr(kr) {}
    
        bool intersect(const Ray& ray, double& t, Vector& P, Vector& N) const override;
    
        const Vector& getCenter() const { return center; }
        double getN1() const { return n1; }
        double getN2() const { return n2; }
        Vector center;
        double radius;
        double n1; 
        double n2;
        double kr;
    };

double square(double x);
double visibility(const Vector& P, const Vector& N, const Vector& S, const std::vector<Sphere>& scene);
bool findClosestIntersection(const Ray& ray, const std::vector<Sphere>& scene, 
                            const Sphere*& hit_sphere, Vector& P, Vector& N, double& t_hit);
Vector refraction(const Vector& incident, const Vector& normal, double n1, double n2);

Vector getColor(const Ray& ray, const std::vector<Sphere>& scene, int depth, 
               const Vector& Light, double LightIntensity);