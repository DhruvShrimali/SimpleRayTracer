//room depth = 240
//room height = 350
//room width = 300
//(x, y, z)
// camera bottom left = (170, 110, 0)
// camera top left = (170, 140, 0)
// camera top right = (130, 140, 0)
// camera bottom right = (130, 110, 0)
// camera rays origin = (150, 130, -20)

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

#define MAXBOUNCE 15


float randomFloat(){
    return (float)((rand() % 10)/10.0f);
}

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

vect getRandDir(){
    vect rand(randomFloat(), randomFloat(), randomFloat());
    // normalize(rand);
    return rand;
}

class Sphere{
    public:
    vect center;
    double radius;
    double radius2;
    float refl; //between 0 and 1
    float rough;
    colour clr;
    float transparency; // between 0 and 1
    float refr_index;
    bool light;

    Sphere(const vect& c, double r, float refl, float rog, colour col, float transp, float refr, bool lig) : 
        center(c), 
        radius(r), 
        refl(refl),
        rough(rog), 
        clr(col), 
        transparency(transp),
        refr_index(refr),
        light(lig){
        radius2 = r * r;
    }
    Sphere(const vect& c, double r, float refl, float rog, colour col, float transp, bool lig) : 
        center(c), 
        radius(r), 
        refl(refl),
        rough(rog), 
        clr(col), 
        transparency(transp),
        light(lig){
        radius2 = r * r;
        refr_index=1.0;
    }
    Sphere(const vect& c, double r, float refl, float rog, colour col, bool lig) : 
        center(c), 
        radius(r), 
        refl(refl),
        rough(rog), 
        clr(col), 
        light(lig){
        radius2 = r * r;
        transparency=0.0;
        refr_index=1.0;
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
        reflectedDir = reflectedDir + (getRandDir()*rough);
        return reflectedDir;
    }

    vect refract(vect& raydir, vect& normal, bool enteringMedium){
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
    bool enteringMedium = (dot(normal, dir) < 0); // Ray is entering medium if dot product is positive
    // Calculate reflectOrig based on whether the ray is entering or exiting the sphere
    vect reflectOrig;
    if (enteringMedium) // when ray goes from outside to inside
        reflectOrig = intersection + (normal * 0.0001); // Start slightly outside the sphere
    else // when ray goes from inside to outside
        reflectOrig = intersection - (normal * 0.0001); // Start slightly outside the sphere

    // Calculate color recursively
    colour reflectionColor = recursive(reflectOrig, reflectDir, bounces + 1, spheres);
    colour finalColor = colour(0);
    reflectionColor.mult(1.0 - sphere->transparency);
    reflectionColor.mult(sphere->clr);
    finalColor.add(reflectionColor);

    
    if (sphere->transparency > 0.01) {
        vect refractDir;
        vect refractOrig;
        if (enteringMedium) // when ray entering sphere
            refractOrig = intersection - (normal * 0.0001); // Start slightly inside the sphere
        else // when ray exiting sphere
            refractOrig = intersection + (normal * 0.0001); // Start slightly outside the sphere
        refractDir = sphere->refract(dir, normal, enteringMedium);
        normalize(refractDir);

        colour refractionColor = recursive(refractOrig, refractDir, bounces + 1, spheres);
        refractionColor.mult(sphere->transparency);
        refractionColor.mult(sphere->clr);
        finalColor.add(refractionColor);
    }

    // Compute final color considering reflection and refraction
    finalColor.mult(sphere->refl);
    return finalColor;
}


colour getAverageColour(vect& orig, vect& dir, std::vector<Sphere> *spheres, int numSamples, int k) {
    colour avgColour(0);
    double half = numSamples/(2 * k);
    dir.x -= half;
    dir.y -= half; 
    k=k*10;
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

    avgColour.r /= (numSamples * numSamples);
    avgColour.g /= (numSamples * numSamples);
    avgColour.b /= (numSamples * numSamples);

    return avgColour;
}

int main() {
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<Sphere> spheres;
    srand(time(0));
    //             center(x, y, z), radius, illumination, roughness, colour(r, g, b), isLight
    spheres.push_back(Sphere(vect(100, 400, 120), 100, 3, 0, colour(255, 255, 255), true)); //light
    spheres.push_back(Sphere(vect(200, 400, 120), 100, 3, 0, colour(255, 255, 255), true)); //light

    //             center(x, y, z), radius, reflectivity, roughness, colour(r, g, b), isLight
    spheres.push_back(Sphere(vect(60, 40, 210), 35, 0.999, 0, colour(250, 255, 255), 0.9, 1.5, false)); // three spheres
    spheres.push_back(Sphere(vect(50, 130, 120), 45, 0.95, 0, colour(200, 200, 0), false));
    spheres.push_back(Sphere(vect(200, 70, 200), 70, 0.999, 0, colour(255, 200, 255), false));
    
    spheres.push_back(Sphere(vect(150, -5000, 120), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //bottom floor
    spheres.push_back(Sphere(vect(150, 5350, 120), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //ceiling
    spheres.push_back(Sphere(vect(150, 175, 5240), 5000, 0.99, 1.4, colour(200, 255, 200), false)); //front wall
    // spheres.push_back(Sphere(vect(150, 175, -5021), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //back wall
    spheres.push_back(Sphere(vect(5300, 175, 120), 5000, 0.99, 1.4, colour(255, 200, 200), false)); //left wall
    spheres.push_back(Sphere(vect(-5000, 175, 120), 5000, 0.99, 1.4, colour(200, 200, 255), false)); //right wall
    
    std::ofstream MyFile("intermediate.txt");

    int ray_per_pixel = 20; // Adjust the number of samples per pixel
    int neighbour_per_pixel = 7; // Adjust to smoothen image
    unsigned long long int time=1;
    time *= spheres.size(); // O(n)
    time *= ray_per_pixel; // O(n^2)
    time *= ray_per_pixel;
    time *= MAXBOUNCE; // O(n^2)
    time *= MAXBOUNCE;
    time *= neighbour_per_pixel; // O(n^2)
    time *= neighbour_per_pixel;
    time /= 60500; //obtained experimentally
    int time_min = time/60;
    time = time%60;
    int time_hour = time_min/60;
    time_min = time_min%60;

    printf("Please wait for less than %d hr, %d min, %lld sec\n",time_hour, time_min, time);
    std::vector<std::string> IO;
    
    vect orig(150, 125, -20);
    for (int i = 0; i < 30 * ray_per_pixel; ++i) {
        for (int j = 0; j < 40 * ray_per_pixel; ++j) {
            vect dir(170.0 - ((double)j / (double)ray_per_pixel), 140.0 - ((double)i / (double)ray_per_pixel), 0);
            colour clr = getAverageColour(orig, dir, &spheres, neighbour_per_pixel, ray_per_pixel); // Adjust the number of samples
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