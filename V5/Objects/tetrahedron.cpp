#include <cmath>
#include "tetrahedron.hpp"
#include "shape.hpp"

bool rayTriangleIntersect(vect& orig, vect& dir, vect& v0, vect& v1, vect& v2, double& t){
    double EPSILON = 1e-8;
    vect edge1 = v1 - v0;
    vect edge2 = v2 - v0;
    vect h = cross(dir, edge2);
    double a = dot(edge1, h);
    
    if (a > -EPSILON && a < EPSILON)
        return false; // This means the ray is parallel to the triangle.

    double f = 1.0 / a;
    vect s = orig - v0;
    double u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    vect q = cross(s, edge1);
    double v = f * dot(dir, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    t = f * dot(edge2, q);
    return t > EPSILON; // The intersection is valid if t > 0
}

Tetrahedron::Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3, float energy_retention, float rog, colour col, float transp, float refr, bool lig): 
    Shape(energy_retention, rog, col, transp, refr, lig, false), 
    top(t), 
    b1(b1), 
    b2(b2), 
    b3(b3) {}
Tetrahedron::Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3, float energy_retention, float rog, colour col, bool lig): 
    Shape(energy_retention, rog, col, 0.0, 1.0, lig, false), 
    top(t), 
    b1(b1), 
    b2(b2), 
    b3(b3) {}
Tetrahedron::Tetrahedron(const vect& b1, const vect& b2, const vect& b3, float energy_retention, float rog, colour col, float transp, float refr, bool lig): 
    Shape(energy_retention, rog, col, transp, refr, lig, false), 
    top(b3),
    b1(b1), 
    b2(b2), 
    b3(b3) {
        isTriangle = true;
    }
Tetrahedron::Tetrahedron(const vect& b1, const vect& b2, const vect& b3, float energy_retention, float rog, colour col, bool lig): 
    Shape(energy_retention, rog, col, 0.0, 1.0, lig, false), 
    top(b3),
    b1(b1), 
    b2(b2), 
    b3(b3) {
        isTriangle = true;
    }


int Tetrahedron::intersect(vect& rayorig, vect& raydir, double& t0, double &t1){
    bool hit = 0;
    double t;
    
    // Check intersection with each triangle
    if(isTriangle){
        if (rayTriangleIntersect(rayorig, raydir, b1, b2, b3, t)) {
            t0 = t;
            hit = 4;
        }
    }
    else{
        if (rayTriangleIntersect(rayorig, raydir, top, b1, b2, t)) {
            t0 = t;
            hit = 1;
        }
        if (rayTriangleIntersect(rayorig, raydir, top, b2, b3, t)) {
            if (t < t0) t0 = t;
            hit = 2;
        }
        if (rayTriangleIntersect(rayorig, raydir, top, b3, b1, t)) {
            if (t < t0) t0 = t;
            hit = 3;
        }
        if (rayTriangleIntersect(rayorig, raydir, b1, b2, b3, t)) {
            if (t < t0) t0 = t;
            hit = 4;
        }
    }

    return hit;
}

vect Tetrahedron::getNormal(const vect& intersection, int latest_hit_face) const {
    vect normal;

    // Calculate the normal based on the latest hit face
    switch (latest_hit_face) {
        case 1:
            // Triangle formed by top, b1, b2
            normal = cross(b1 - top, b2 - top);
            if (dot(normal, intersection - top) < 0) {
                normal = normal * -1; // Flip normal if it's facing inward
            }
            break;

        case 2:
            // Triangle formed by top, b2, b3
            normal = cross(b2 - top, b3 - top);
            if (dot(normal, intersection - top) < 0) {
                normal = normal * -1; // Flip normal if it's facing inward
            }
            break;

        case 3:
            // Triangle formed by top, b3, b1
            normal = cross(b3 - top, b1 - top);
            if (dot(normal, intersection - top) < 0) {
                normal = normal * -1; // Flip normal if it's facing inward
            }
            break;

        case 4:
            // Triangle formed by b1, b2, b3
            normal = cross(b2 - b1, b3 - b1);
            if (dot(normal, intersection - b1) < 0) {
                normal = normal * -1; // Flip normal if it's facing inward
            }
            break;

        default:
            // Handle error case
            return vect(0, 0, 0);
    }
    normalize(normal);
    return normal;
}

vect Tetrahedron::reflect(vect& raydir, vect& normal){
    double mult = 2 * dot(raydir, normal);
    vect reflectedDir = raydir - (normal * mult);
    reflectedDir = reflectedDir + (randomVector()*rough);
    return reflectedDir;
}

vect Tetrahedron::refract(vect& raydir, vect& normal, bool enteringMedium){
    double n = enteringMedium ? 1.0 / refr_index : refr_index;
    normalize(normal);
    normalize(raydir);
    double cosI = dot(normal, raydir);
    double sinT2 = n * n * (1.0 - cosI * cosI);

    if (sinT2 > 1.0) {
        // Total internal reflection, reverse normal direction and compute reflected direction
        normal = normal*-1; // Reverse the normal direction
        return reflect(raydir, normal); // Compute the reflected direction
    }
    double cosT = sqrt(1.0 - sinT2);

    vect temp = raydir*n;
    vect temp2 = normal * ((n * cosI) - cosT);
    return temp + temp2;
}