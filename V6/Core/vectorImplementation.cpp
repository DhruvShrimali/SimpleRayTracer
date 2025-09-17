#include <cmath>
#include "../Utils/random.h"
#include "vectorImplementation.hpp"


vect::vect(){
    x=0;
    y=0;
    z=0;
}
vect::vect(double k){
    x=k;
    y=k;
    z=k;
}
vect::vect(double i, double j, double k){
    x=i;
    y=j;
    z=k;
}


double vect::length2(){ 
    return (x * x) + (y * y) + (z * z); 
}
double vect::length(){ 
    return sqrt(length2()); 
}

bool vect::operator==(const vect& other){
    return (x == other.x && y == other.y && z == other.z);
}


vect operator-(const vect& a, const vect& b) {
    return vect(a.x - b.x, a.y - b.y, a.z - b.z);
}

vect operator-(const vect& a) {
    return vect(-a.x, -a.y, -a.z);
}

vect operator+(const vect& a, const vect& b) {
    return vect(a.x + b.x, a.y + b.y, a.z + b.z);
}

void operator*=(vect& v, double c) {
    v.x *= c;
    v.y *= c;
    v.z *= c;
}

vect operator*(const vect& v, double c) {
    return vect(v.x * c, v.y * c, v.z * c);
}

void normalize(vect& a) {
    double nor2 = a.length2();
    if (nor2 > 0) {
        double invNor = 1 / sqrt(nor2);
        a *= invNor;
    }
}

double dot(const vect &a, const vect &b){
    double i = (a.x) * (b.x);
    double j = (a.y) * (b.y);
    double k = (a.z) * (b.z);
    return i+j+k;
}

vect cross(const vect &a, const vect &b){
    double i = (a.y * b.z) - (b.y * a.z);
    double j = (b.x * a.z) - (a.x * b.z);
    double k = (a.x * b.y) - (b.x * a.y);
    return vect(i, j, k);
}

vect randomVector(){
    // normalize(rand);
    return vect(randomFloat(), randomFloat(), randomFloat());;
}