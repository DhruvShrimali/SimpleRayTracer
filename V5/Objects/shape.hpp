#pragma once

#include "../Core/vectorImplementation.hpp"
#include "../Core/colourImplementation.hpp"

class Shape{
    public:
    float energy_retention; // Energy loss factor
                            // 0 = absorbs all light like a black hole
                            // 1 = reflects/refracts all rays without energy loss
                            // For light sources, this value controls brightness and can be >1.
    float rough;            // Surface roughness
                            // 0 = perfect mirror
                            // Higher values → more scattered/refused reflection, rough surface
    colour clr;             // Base color of the shape
    float transparency;     // Transparency factor (0 = opaque, 1 = fully transparent)
    float refr_index;       // Refractive index (controls bending of light through transparent objects)
    bool light;             // True if the shape is a light source
    bool isSphere;          // Flag to distinguish if the object is a sphere (true) or another shape like tetrahedron

    // Constructor to initialize all properties
    Shape(float energy_retention, float rough, colour col, float transp, float refr, bool lig, bool sphere);

    // Virtual methods for ray interactions
    virtual int intersect(vect& rayorig, vect& raydir, double& t0, double& t1) = 0; // Intersection with a ray
    virtual vect getNormal(const vect& intersection, const int latest_hit_face) const = 0; // Surface normal at intersection
    virtual vect reflect(vect& raydir, vect& normal) = 0; // Reflection ray calculation
    virtual vect refract(vect& raydir, vect& normal, bool enteringMedium) = 0; // Refraction ray calculation

    virtual ~Shape() {} // Virtual destructor
};