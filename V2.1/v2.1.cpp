//room depth = 240
//room height = 350
//room width = 300
//(x, y, z)
// camera bottom left = (170, 110, 0)
// camera top left = (170, 140, 0)
// camera top right = (130, 140, 0)
// camera bottom right = (130, 110, 0)
// camera rays origin = (150, 125, -20)

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <string>
#include <limits>
#include <chrono>

#define MAXBOUNCE 5

class colour{
    public:
    int r;
    int g;
    int b;

    colour(){
        r=0;
        g=0;
        b=0;
    }
    colour(int k){
        r=k;
        g=k;
        b=k;
    }
    colour(int i, int j, int k){
        r=i;
        g=j;
        b=k;
    }

    void add(colour &clr){
        r+=(clr.r);
        g+=(clr.g);
        b+=(clr.b);
        if (r>255)
            r=255;
        if (g>255)
            g=255;
        if (b>255)
            b=255;
    }
    void mult(double k){
        double r1, g1, b1;
        r1=r*k;
        g1=g*k;
        b1=b*k;
        r=(int)r1;
        g=(int)g1;
        b=(int)b1;
    }
    void mult(int i, int j, int k){
        int r1, g1, b1;
        r1=r*i;
        g1=g*j;
        b1=b*k;
        r=(int)floor(r1/255);
        g=(int)floor(g1/255);
        b=(int)floor(b1/255);
    }
    void mult(colour &clr){
        int r1, g1, b1;
        r1=r*(clr.r);
        g1=g*(clr.g);
        b1=b*(clr.b);
        r=(int)floor(r1/255);
        g=(int)floor(g1/255);
        b=(int)floor(b1/255);
    }
};

class vect{
    public:
    double x, y, z;

    vect(){
        x=0;
        y=0;
        z=0;
    }
    vect(double k){
        x=k;
        y=k;
        z=k;
    }
    vect(double i, double j, double k){
        x=i;
        y=j;
        z=k;
    }

    
    double length2(){ 
        return (x * x) + (y * y) + (z * z); 
    }
    double length(){ 
        return sqrt(length2()); 
    }
    
};

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

double dot(vect &a, vect &b){
    double i = (a.x) * (b.x);
    double j = (a.y) * (b.y);
    double k = (a.z) * (b.z);
    return i+j+k;
}

vect cross(vect &a, vect &b){
    double i = (a.y * b.z) - (b.y * a.z);
    double j = (b.x * a.z) - (a.x * b.z);
    double k = (a.x * b.y) - (b.x * a.y);
    return vect(i, j, k);
}

class Sphere{
    public:
    vect center;
    double radius;
    double radius2;
    float refl; //between 0 and 1
    colour clr;
    bool light;

    Sphere(const vect& c, double r, float refl, colour col, bool lig) : 
        center(c), 
        radius(r), 
        refl(refl),
        clr(col), 
        light(lig){
        radius2 = r * r;
    }

    bool intersect(vect& rayorig, vect& raydir, double& t0, double& t1){
        vect l = center - rayorig;
        double tca = dot(l, raydir);
        if (tca < 0) return false;

        double d2 = l.length2() - (tca * tca);
        if (d2 > radius2) return false;

        double thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;
        return true;
    }
    vect reflect(vect& raydir, vect& normal) {
        double mult = 2 * dot(raydir, normal);
        vect reflectedDir = raydir - (normal * mult);
        return reflectedDir;
    }

};


colour recursive(vect& orig, vect& dir, int bounces, std::vector<Sphere> *spheres){
    if(bounces >= MAXBOUNCE){
        return colour(0);
    }

    double tnear = INFINITY;
    Sphere* sphere = nullptr;

    std::vector<Sphere>::iterator itr;
    for (itr = spheres->begin(); itr < spheres->end(); itr++){
        double t0 = INFINITY, t1 = INFINITY;
        if (itr->intersect(orig, dir, t0, t1)) {
            if(itr->light){
                colour temp = itr->clr;
                temp.mult(itr->refl);
                return temp;
            }
            else{
                if (t0 < 0) t0 = t1;
                if (t0 < tnear) {
                    tnear = t0;
                    sphere = &(*itr);
                }
            }
        }
    }
    if (!sphere){
        return colour(0); // Default color if no intersection
    }

    if(sphere->clr.r==0 && sphere->clr.g==0 && sphere->clr.b==0){
        return colour(0);
    }
    // Intersection point
    vect intersection = orig + (dir * tnear);
    // Normal at intersection point
    vect normal = intersection - sphere->center;
    normalize(normal);

    // Reflected ray direction
    vect reflectDir = sphere->reflect(dir, normal);
    normalize(reflectDir);
    
    // Reflected ray origin
    vect reflectOrig = intersection + (normal * 0.0001); // Start slightly outside the sphere

    // Calculate color recursively
    colour reflectionColor = recursive(reflectOrig, reflectDir, bounces + 1, spheres);
    colour finalColor = colour(0);
    reflectionColor.mult(sphere->clr);
    finalColor.add(reflectionColor);

    // Compute final color considering reflection and refraction
    finalColor.mult(sphere->refl);
    return finalColor;
}


//Function that checks k*k light rays around ray and returns average of those colors
colour getAverageColour(vect& orig, vect& dir, std::vector<Sphere> *spheres, int numSamples, int k) {
    colour avgColour(0);
    double half = numSamples/(2 * k);
    dir.x -= half;
    dir.y -= half; 
    k=k*10; //experimental value
    for (int i = 0; i < numSamples; ++i) {
        for (int j = 0; j < numSamples; ++j) {
            double x = dir.x + ((double)j/(double)k);
            double y = dir.y + ((double)i/(double)k);
            vect temp_dir(x, y, dir.z);
            vect dirSample = temp_dir - orig;
            normalize(dirSample); 
            colour clr = recursive(orig, dirSample, 0, spheres);
            avgColour.r += clr.r;
            avgColour.g += clr.g;
            avgColour.b += clr.b;
        }
    }
    int  divider = numSamples * numSamples;
    avgColour.r /= divider;
    avgColour.g /= divider;
    avgColour.b /= divider;

    return avgColour;
}

int main() {
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<Sphere> spheres;

    //             center(x, y, z), radius, intensity(for lights), colour(r, g, b), isLight
    spheres.push_back(Sphere(vect(100, 400, 120), 100, 3, colour(255, 255, 255), true)); //light
    spheres.push_back(Sphere(vect(200, 400, 120), 100, 3, colour(255, 255, 255), true)); //light

    //             center(x, y, z), radius, reflectivity, colour(r, g, b), isLight
    spheres.push_back(Sphere(vect(60, 40, 210), 35, 0.99, colour(250, 255, 255), false)); // three spheres
    spheres.push_back(Sphere(vect(50, 130, 120), 45, 0.99, colour(200, 200, 0), false));
    spheres.push_back(Sphere(vect(200, 70, 200), 70, 0.99, colour(255, 200, 255), false));
    
    spheres.push_back(Sphere(vect(150, -5000, 120), 5000, 0.95, colour(255, 255, 255), false)); //bottom floor
    spheres.push_back(Sphere(vect(150, 5350, 120), 5000, 0.95, colour(255, 255, 255), false)); //ceiling
    spheres.push_back(Sphere(vect(150, 175, 5240), 5000, 0.95, colour(200, 255, 200), false)); //front wall
    spheres.push_back(Sphere(vect(150, 175, -5021), 5000, 0.95, colour(255, 255, 255), false)); //back wall
    spheres.push_back(Sphere(vect(5300, 175, 120), 5000, 0.95, colour(255, 200, 200), false)); //left wall
    spheres.push_back(Sphere(vect(-5000, 175, 120), 5000, 0.95, colour(200, 200, 255), false)); //right wall
    
    std::ofstream MyFile("intermediate.txt");

    int ray_per_pixel = 20; // Adjust the number of samples per pixel
    int neighbour_per_pixel = 5; // Adjust to smoothen image
    
    vect orig(150, 125, -20);
    for (int i = 0; i < 30 * ray_per_pixel; ++i) {
        for (int j = 0; j < 40 * ray_per_pixel; ++j) {
            vect dir(170.0 - ((double)j / (double)ray_per_pixel), 140.0 - ((double)i / (double)ray_per_pixel), 0);
            colour clr = getAverageColour(orig, dir, &spheres, neighbour_per_pixel, ray_per_pixel);
            
            MyFile << clr.r << " " << clr.g << " " << clr.b << ",";
            
        }
        MyFile << "\n";    
    }

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate duration
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - begin);

    // Output execution time
    std::cout << "Execution time: " << duration.count() / 60 << " minutes and " << duration.count() % 60 << " seconds." << std::endl;
    return 0;
}
