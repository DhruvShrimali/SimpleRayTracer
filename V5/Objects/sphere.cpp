#include <cmath>
#include "sphere.hpp"
#include "shape.hpp"
#include "../Core/vectorImplementation.hpp"


Sphere::Sphere(const vect& c, double r, float energy_retention, float rog, colour col, float transp, float refr, bool lig):
    Shape(energy_retention, rog, col, transp, refr, lig, true), 
    center(c), 
    radius(r) {
    radius2 = r * r;
}
Sphere::Sphere(const vect& c, double r, float energy_retention, float rog, colour col,bool lig):
    Shape(energy_retention, rog, col, 0.0, 1.0, lig, true), 
    center(c), 
    radius(r) {
    radius2 = r * r;
}

int Sphere::intersect(vect& rayorig, vect& raydir, double& t0, double& t1){
    vect l = center - rayorig;
    double tca = dot(l, raydir);
    if (tca < 0) return 0;

    double d2 = l.length2() - (tca * tca);
    if (d2 > radius2) return 0;

    double thc = sqrt(radius2 - d2);
    t0 = tca - thc;
    t1 = tca + thc;
    return 1;
}

vect Sphere::getNormal(const vect& intersection, int latest_hit_face) const {
    vect normal = intersection - center;
    normalize(normal);
    return normal;
}

vect Sphere::reflect(vect& raydir, vect& normal){
    double mult = 2 * dot(raydir, normal);
    vect reflectedDir = raydir - (normal * mult);
    reflectedDir = reflectedDir + (randomVector()*rough);
    return reflectedDir;
}

vect Sphere::refract(vect& raydir, vect& normal, bool enteringMedium){
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
