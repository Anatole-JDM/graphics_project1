#include "vector.h"
#include <vector>
#include <algorithm>
#include <cmath>

// Vector implementation
Vector::Vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

double Vector::norm2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
}

double Vector::norm() const {
    return sqrt(norm2());
}

void Vector::normalize() {
    double n = norm();
    data[0] /= n;
    data[1] /= n;
    data[2] /= n;
}

double Vector::operator[](int i) const { 
    return data[i]; 
}

double& Vector::operator[](int i) { 
    return data[i]; 
}

Vector Vector::operator-() const {
    return Vector(-data[0], -data[1], -data[2]);
}

Vector& Vector::operator+=(const Vector& other) {
    data[0] += other[0];
    data[1] += other[1];
    data[2] += other[2];
    return *this;
}


// Vector operators
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}

Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double square(double x) {
    return x * x;
}


bool Sphere::intersect(const Ray& ray, double& t, Vector& P, Vector& N) const {
    Vector arc = ray.origin - center;
    double b = 2 * dot(ray.direction, arc);
    double delta = square(b) - 4 * (dot(arc, arc) - radius * radius);

    if (delta < 0) {
        return false;
    }

    double t1 = (-b - sqrt(delta)) / 2;
    double t2 = (-b + sqrt(delta)) / 2;

    if (t2 < 0) {
        return false;
    }
    if (t1 > 0) {
        t = t1;
    } else {
        t = t2;
    }

    P = ray.origin + t * ray.direction;
    N = P - center;
    N.normalize();
    return true;
}
