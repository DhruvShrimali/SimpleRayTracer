#pragma once

#include "shape.hpp"


// Checks if a ray intersects a triangle defined by vertices v0, v1, v2
// Returns true if intersection occurs, and sets t to the distance along the ray
bool rayTriangleIntersect(vect& orig, vect& dir, vect& v0, vect& v1, vect& v2, double& t);

class Tetrahedron : public Shape {
public:
    vect top;   // Top vertex of the tetrahedron (apex)
    vect b1;    // Base vertex 1
    vect b2;    // Base vertex 2
    vect b3;    // Base vertex 3
    bool isTriangle = false; // Flag to indicate if this object should be treated as a single triangle

    // Constructor with apex and base vertices, full properties (energy_retention, roughness, color, transparency, refractive index, light flag)
    Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3,
                float energy_retention, float rog, colour col, float transp, float refr, bool lig);

    // Constructor with apex and base, without transparency/refractive index
    Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3,
                float energy_retention, float rog, colour col, bool lig);

    // Constructor with only base vertices, full properties (for flat triangles)
    Tetrahedron(const vect& b1, const vect& b2, const vect& b3,
                float energy_retention, float rog, colour col, float transp, float refr, bool lig);

    // Constructor with only base vertices, without transparency/refractive index (for flat triangles)
    Tetrahedron(const vect& b1, const vect& b2, const vect& b3,
                float energy_retention, float rog, colour col, bool lig);

    // Ray-tetrahedron intersection
    // Returns face number of triangle if intersection occurs, sets t0/t1 to distances along the ray
    int intersect(vect& rayorig, vect& raydir, double& t0, double &t1) override;

    // Returns the surface normal at a given intersection point given face number of triangle
    vect getNormal(const vect& intersection, int latest_hit_face) const override;

    // Calculates the reflected ray direction given an incident ray and surface normal
    vect reflect(vect& raydir, vect& normal) override;

    // Calculates the refracted ray direction given an incident ray, surface normal, and medium entry/exit
    vect refract(vect& raydir, vect& normal, bool enteringMedium) override;
};
