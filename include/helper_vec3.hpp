#ifndef HELPER_VEC3_HPP
#define HELPER_VEC3_HPP

#include <cmath>
#include <iostream>

// =========================================================
// Shared Vec3 utilities for all modules
// =========================================================

struct Vec3 {
    double x, y, z;

    // addition
    Vec3 operator+(const Vec3& b) const {
        return {x + b.x, y + b.y, z + b.z};
    }

    // subtraction
    Vec3 operator-(const Vec3& b) const {
        return {x - b.x, y - b.y, z - b.z};
    }

    // scalar multiply
    Vec3 operator*(double s) const {
        return {x * s, y * s, z * s};
    }

    // scalar divide
    Vec3 operator/(double s) const {
        return {x / s, y / s, z / s};
    }
};

// ---------------------------------------------------------
// Dot product
// ---------------------------------------------------------
inline double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// ---------------------------------------------------------
// Vector norm
// ---------------------------------------------------------
inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

// ---------------------------------------------------------
// Normalize vector
// ---------------------------------------------------------
inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    if (n < 1e-12) {
        std::cerr << "ERROR: normalize() encountered zero-length vector.\n";
        exit(1);
    }
    return v / n;
}

// ---------------------------------------------------------
// Cross product
// ---------------------------------------------------------
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

#endif

