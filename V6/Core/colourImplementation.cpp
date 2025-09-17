#include "colourImplementation.hpp"
#include "../config.hpp"

double correct_colour(double original_value){
    if(IMG_BITS == 8){
        return original_value;
    }
    else if(IMG_BITS == 16){
        original_value = original_value*COLOR_MAX/255;
        return original_value;
    }
    else{
        return 0.0;
    }
}  

colour::colour(){
    r=0;
    g=0;
    b=0;
}
colour::colour(short int k){
    r=correct_colour(k);
    g=correct_colour(k);
    b=correct_colour(k);
}
colour::colour(int i, int j, int k){
    r=correct_colour(i);
    g=correct_colour(j);
    b=correct_colour(k);
}

void colour::add(colour &clr){
    r+=(clr.r);
    g+=(clr.g);
    b+=(clr.b);
    if (r>COLOR_MAX)
        r=COLOR_MAX;
    if (g>COLOR_MAX)
        g=COLOR_MAX;
    if (b>COLOR_MAX)
        b=COLOR_MAX;
}

void colour::free_add(colour &clr){
    r+=(clr.r);
    g+=(clr.g);
    b+=(clr.b);
}

void colour::mult(double k){
    double r1, g1, b1;
    r1=r*k;
    g1=g*k;
    b1=b*k;
    r=r1;
    g=g1;
    b=b1;
}
void colour::mult(int i, int j, int k){
    int r1, g1, b1;
    r1=r*i;
    g1=g*j;
    b1=b*k;
    r=r1/(double)COLOR_MAX;
    g=g1/(double)COLOR_MAX;
    b=b1/(double)COLOR_MAX;
}
void colour::mult(colour &clr){
    int r1, g1, b1;
    r1=r*(clr.r);
    g1=g*(clr.g);
    b1=b*(clr.b);
    r=r1/(double)COLOR_MAX;
    g=g1/(double)COLOR_MAX;
    b=b1/(double)COLOR_MAX;
}
void colour::div(double k){
    double r1, g1, b1;
    r1=r/k;
    g1=g/k;
    b1=b/k;
    r=r1;
    g=g1;
    b=b1;
}
void colour::trim(){
    if(r < 0){
        r = 0;
    }
    else if (r > COLOR_MAX){
        r = COLOR_MAX;
    }
    if(g < 0){
        g = 0;
    }
    else if (g > COLOR_MAX){
        g = COLOR_MAX;
    }
    if(b < 0){
        b = 0;
    }
    else if (b > COLOR_MAX){
        b = COLOR_MAX;
    }
}
