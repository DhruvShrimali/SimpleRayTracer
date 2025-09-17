#pragma once

#include "vectorImplementation.hpp"

class vect{
    public:
    double x, y, z;

    vect();
    vect(double k);
    vect(double i, double j, double k);

    
    double length2();
    double length();

    bool operator==(const vect& other);
    
};

vect operator-(const vect& a, const vect& b) ;

vect operator-(const vect& a) ;

vect operator+(const vect& a, const vect& b) ;

void operator*=(vect& v, double c) ;

vect operator*(const vect& v, double c) ;

void normalize(vect& a) ;

double dot(const vect &a, const vect &b);

vect cross(const vect &a, const vect &b);

vect randomVector();