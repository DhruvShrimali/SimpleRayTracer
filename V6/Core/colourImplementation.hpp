#pragma once

class colour{
    public:
    double r;
    double g;
    double b;

    colour();
    colour(short int k);
    colour(int i, int j, int k);

    void add(colour &clr);
    void free_add(colour &clr);
    void mult(double k);
    void mult(int i, int j, int k);
    void mult(colour &clr);
    void div(double k);
    void trim();
};