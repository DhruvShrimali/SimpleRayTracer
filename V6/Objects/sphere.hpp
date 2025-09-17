#pragma once

#include <cmath>
#include "shape.hpp"
#include "../Core/vectorImplementation.hpp"


class Sphere : public Shape {
public:
    vect center;      // Center position of the sphere in 3D space
    double radius;    // Radius of the sphere
    double radius2;   // Squared radius, precomputed for faster intersection calculations

    // Constructor with full properties (reflection, roughness, color, transparency, refractive index, light flag)
    Sphere(const vect& c, double r, float energy_retention, float rog, colour col, float transp, float refr, bool lig);

    // Constructor for simpler spheres (without transparency and refractive index)
    Sphere(const vect& c, double r, float energy_retention, float rog, colour col, bool lig);

    // Ray-sphere intersection
    // Computes if and where a ray intersects the sphere
    int intersect(vect& rayorig, vect& raydir, double& t0, double& t1) override;

    // Returns the surface normal at a given intersection point
    vect getNormal(const vect& intersection, int latest_hit_face) const override;

    // Calculates the reflected ray direction given an incident ray and surface normal
    vect reflect(vect& raydir, vect& normal) override;

    // Calculates the refracted ray direction given an incident ray, surface normal, and medium entry/exit
    vect refract(vect& raydir, vect& normal, bool enteringMedium) override;
};
