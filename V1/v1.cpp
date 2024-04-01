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
    bool clr;
    bool light;

    Sphere(const vect& c, double r, bool col, bool lig) : 
        center(c), 
        radius(r), 
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

};


int main() {
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<Sphere> spheres;
    srand(time(0));

    //             center(x, y, z), radius, colour, isLight
    spheres.push_back(Sphere(vect(60, 40, 210), 35, true, false)); // three spheres
    spheres.push_back(Sphere(vect(50, 130, 120), 45, true, false));
    spheres.push_back(Sphere(vect(200, 70, 200), 70, true, false));
    
    spheres.push_back(Sphere(vect(150, -5000, 120), 5000, false, false)); //bottom floor
    spheres.push_back(Sphere(vect(150, 5350, 120), 5000, false, false)); //ceiling
    spheres.push_back(Sphere(vect(150, 175, 5240), 5000, false, false)); //front wall
    spheres.push_back(Sphere(vect(150, 175, -5021), 5000, false, false)); //back wall
    spheres.push_back(Sphere(vect(5300, 175, 120), 5000, false, false)); //left wall
    spheres.push_back(Sphere(vect(-5000, 175, 120), 5000, false, false)); //right wall
    
    std::ofstream MyFile("intermediate.txt");

    int ray_per_pixel = 2; // Adjust the number of samples per pixel
    
    vect orig(150, 125, -20);
    for (int i = 0; i < 30 * ray_per_pixel; ++i) {
        for (int j = 0; j < 40 * ray_per_pixel; ++j) {
            vect dir(170.0 - (j / ray_per_pixel), 140.0 - (i / ray_per_pixel), 0);
            dir = dir-orig;
            normalize(dir);
            bool flag;
            double tnear = INFINITY;
            Sphere* sphere = nullptr;

            std::vector<Sphere>::iterator itr;
            for (itr = spheres.begin(); itr < spheres.end(); itr++){
                double t0 = INFINITY, t1 = INFINITY;
                if (itr->intersect(orig, dir, t0, t1)) {
                    if (t0 < 0) t0 = t1;
                    if (t0 < tnear) {
                        tnear = t0;
                        sphere = &(*itr);
                    }
                }
            }
            if(!sphere){
                flag=false;
            }
            else{
                flag = sphere->clr;
            }
            if(flag){
                MyFile << "##";
            }
            else{
                MyFile << "  ";
            }
            
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