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
#include <thread>

#define MAXBOUNCE 20 // Defines the maxiumum number of reflection calculated (affects brightness of scene)
#define TIMELIMIT 180 // Define the maximum amount of time you want to wait for a render
#define TIMEMODIFIER 5000000 // Adjust the to get approximate time taken to execute at your machine
#define ray_per_pixel 20 // Adjust the number of samples per pixel
#define neighbour_per_pixel  7 // Adjust to smoothen image

float randomFloat(){
    // Convert to [-0.5, 0.5] range
    return (float)(((rand() % 1000)-500)/1000.0f);
}

class colour{
    public:
    double r;
    double g;
    double b;

    colour(){
        r=0;
        g=0;
        b=0;
    }
    colour(short int k){
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
    void free_add(colour &clr){
        r+=(clr.r);
        g+=(clr.g);
        b+=(clr.b);
    }
    void mult(double k){
        double r1, g1, b1;
        r1=r*k;
        g1=g*k;
        b1=b*k;
        r=r1;
        g=g1;
        b=b1;
    }
    void mult(int i, int j, int k){
        int r1, g1, b1;
        r1=r*i;
        g1=g*j;
        b1=b*k;
        r=floor(r1/255);
        g=floor(g1/255);
        b=floor(b1/255);
    }
    void mult(colour &clr){
        int r1, g1, b1;
        r1=r*(clr.r);
        g1=g*(clr.g);
        b1=b*(clr.b);
        r=floor(r1/255);
        g=floor(g1/255);
        b=floor(b1/255);
    }
    void div(double k){
        double r1, g1, b1;
        r1=r/k;
        g1=g/k;
        b1=b/k;
        r=r1;
        g=g1;
        b=b1;
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

    bool operator==(const vect& other){
        return (x == other.x && y == other.y && z == other.z);
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

vect getRandDir(){
    vect rand(randomFloat(), randomFloat(), randomFloat());
    // normalize(rand);
    return rand;
}

bool rayTriangleIntersect(vect& orig, vect& dir, vect& v0, vect& v1, vect& v2, double& t) {
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

class Object{
    public:
    float refl; //between 0 and 1
    float rough;
    colour clr;
    float transparency; // between 0 and 1
    float refr_index;
    bool light;
    bool isSphere;       // Flag to distinguish if the object is a sphere

    Object(float refl, float rough, colour col, float transp, float refr, bool lig, bool sphere): 
        refl(refl), 
        rough(rough), 
        clr(col), 
        transparency(transp), 
        refr_index(refr), 
        light(lig), 
        isSphere(sphere) {}

    // Virtual methods for intersection, normal, reflection, and refraction
    virtual int intersect(vect& rayorig, vect& raydir, double& t0, double& t1) = 0;
    virtual vect getNormal(const vect& intersection, const int latest_hit_face) const = 0;
    virtual vect reflect(vect& raydir, vect& normal) = 0;
    virtual vect refract(vect& raydir, vect& normal, bool enteringMedium) = 0;

    virtual ~Object() {} // Virtual destructor
};

class Sphere : public Object{
    public:
    vect center;
    double radius;
    double radius2;

    Sphere(const vect& c, double r, float refl, float rog, colour col, float transp, float refr, bool lig):
        Object(refl, rog, col, transp, refr, lig, true), 
        center(c), 
        radius(r) {
        radius2 = r * r;
    }
    Sphere(const vect& c, double r, float refl, float rog, colour col,bool lig):
        Object(refl, rog, col, 0.0, 1.0, lig, true), 
        center(c), 
        radius(r) {
        radius2 = r * r;
    }

    int intersect(vect& rayorig, vect& raydir, double& t0, double& t1) override {
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

    vect getNormal(const vect& intersection, int latest_hit_face) const override {
        vect normal = intersection - center;
        normalize(normal);
        return normal;
    }

    vect reflect(vect& raydir, vect& normal) override {
        double mult = 2 * dot(raydir, normal);
        vect reflectedDir = raydir - (normal * mult);
        reflectedDir = reflectedDir + (getRandDir()*rough);
        return reflectedDir;
    }

    vect refract(vect& raydir, vect& normal, bool enteringMedium) override {
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

class Tetrahedron : public Object{
    public:
    vect top;
    vect b1;
    vect b2;
    vect b3;
    int latest_hit_face = -1;
    bool isTriangle = false;

    Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3, float refl, float rog, colour col, float transp, float refr, bool lig): 
        Object(refl, rog, col, transp, refr, lig, false), 
        top(t), 
        b1(b1), 
        b2(b2), 
        b3(b3) {}
    Tetrahedron(const vect& t, const vect& b1, const vect& b2, const vect& b3, float refl, float rog, colour col, bool lig): 
        Object(refl, rog, col, 0.0, 1.0, lig, false), 
        top(t), 
        b1(b1), 
        b2(b2), 
        b3(b3) {}
    Tetrahedron(const vect& b1, const vect& b2, const vect& b3, float refl, float rog, colour col, float transp, float refr, bool lig): 
        Object(refl, rog, col, transp, refr, lig, false), 
        top(b3),
        b1(b1), 
        b2(b2), 
        b3(b3) {
            isTriangle = true;
        }
    Tetrahedron(const vect& b1, const vect& b2, const vect& b3, float refl, float rog, colour col, bool lig): 
        Object(refl, rog, col, 0.0, 1.0, lig, false), 
        top(b3),
        b1(b1), 
        b2(b2), 
        b3(b3) {
            isTriangle = true;
        }


    int intersect(vect& rayorig, vect& raydir, double& t0, double &t1) override {
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

    vect getNormal(const vect& intersection, int latest_hit_face) const override {
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
                return vect(0, 0, 0); // or any default value
        }
        normalize(normal);
        return normal;
    }

    vect reflect(vect& raydir, vect& normal) override {
        double mult = 2 * dot(raydir, normal);
        vect reflectedDir = raydir - (normal * mult);
        reflectedDir = reflectedDir + (getRandDir()*rough);
        return reflectedDir;
    }

    vect refract(vect& raydir, vect& normal, bool enteringMedium) override {
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

colour recursive(vect& orig, vect& dir, int bounces, std::vector<Object*> *objects){
    if(bounces >= MAXBOUNCE){
        return colour(0);
    }

    double tnear = INFINITY;
    Object* object = nullptr;
    int latest_hit_face = 0;

    for (auto itr = objects->begin(); itr < objects->end(); itr++) {
        double t0 = INFINITY, t1 = INFINITY;
        Object* obj = *itr; // Dereference to get the Object pointer
        

        // Check if the object is a sphere
        if(obj->isSphere) {
            int intersect_check = static_cast<Sphere*>(obj)->intersect(orig, dir, t0, t1);
            if (intersect_check) {
                if (obj->light) {
                    colour temp = obj->clr;
                    temp.mult(obj->refl);
                    return temp; // Return light color if the object is a light source
                } else {
                    if (t0 < 0) t0 = t1; // Use t1 if t0 is behind the ray origin
                    if (t0 < tnear) {
                        tnear = t0; // Update the closest intersection distance
                        object = obj; // Update the closest object hit
                        latest_hit_face = intersect_check;
                    }
                }
            }
        } 
        // The object is a tetrahedron
        else{
            int intersect_check = static_cast<Tetrahedron*>(obj)->intersect(orig, dir, t0, t1);
            if (intersect_check) {
                if (obj->light) {
                    colour temp = obj->clr;
                    temp.mult(obj->refl);
                    return temp; // Return light color if the object is a light source
                } else {
                    if (t0 < tnear) {
                        tnear = t0; // Update the closest intersection distance
                        object = obj; // Update the closest object hit
                        latest_hit_face = intersect_check;
                    }
                }
            }
        }
    }
    if (!object){
        return colour(0); // Default color if no intersection
    }

    if(object->clr.r<=0 && object->clr.g<=0 && object->clr.b<=0){
        return colour(0);
    }

    // Intersection point
    vect intersection = orig + (dir * tnear);
    // Normal at intersection point
    vect normal = object->getNormal(intersection, latest_hit_face);

    // Reflected ray direction
    vect reflectDir = object->reflect(dir, normal);
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
    colour reflectionColor = recursive(reflectOrig, reflectDir, bounces + 1, objects);
    colour finalColor = colour(0);
    reflectionColor.mult(1.0 - object->transparency);
    reflectionColor.mult(object->clr);
    finalColor.add(reflectionColor);

    
    if (object->transparency > 0.01) {
        vect refractDir;
        vect refractOrig;
        if (enteringMedium) // when ray entering sphere
            refractOrig = intersection - (normal * 0.0001); // Start slightly inside the sphere
        else // when ray exiting sphere
            refractOrig = intersection + (normal * 0.0001); // Start slightly outside the sphere
        refractDir = object->refract(dir, normal, enteringMedium);
        normalize(refractDir);

        colour refractionColor = recursive(refractOrig, refractDir, bounces + 1, objects);
        refractionColor.mult(object->transparency);
        refractionColor.mult(object->clr);
        finalColor.add(refractionColor);
    }

    // Compute final color considering reflection and refraction
    finalColor.mult(object->refl);
    return finalColor;
}


colour getAverageColour(vect& orig, vect& dir, std::vector<Object*> *objects, int neighborPerPixel, int rayPerPixel) {
    colour avgColour(0);

    rayPerPixel*=10; // Generate 10 times closer than standard
    double half = neighborPerPixel/(2 * rayPerPixel); // Divide by two to get left/top half of the ray direction distance
    dir.x -= half;
    dir.y -= half; 
    
    int divisor = (neighborPerPixel * neighborPerPixel);
    // int EPOCHS = 4;
    // for(int epoch = 0; epoch < EPOCHS; epoch++){
        for (int i = 0; i < neighborPerPixel; ++i) {
            for (int j = 0; j < neighborPerPixel; ++j) {
                double x = dir.x + ((double)j/(double)rayPerPixel);
                double y = dir.y + ((double)i/(double)rayPerPixel);
                vect temp_dir(x, y, dir.z);
                vect dirSample = temp_dir - orig;
                normalize(dirSample); 
                colour clr = recursive(orig, dirSample, 0, objects);
                clr.r /= divisor;
                clr.g /= divisor;
                clr.b /= divisor;
                avgColour.free_add(clr);
            }
        }
    // }

    // avgColour.r /= EPOCHS;
    // avgColour.g /= EPOCHS;
    // avgColour.b /= EPOCHS;

    return avgColour;
}



// Function to render a portion of the image using thread
void renderImagePart(int startRow, int endRow, int rayPerPixel, int neighborPerPixel, std::vector<Object*> *objects, colour **image, int WIDTH) {
    vect orig(150, 125, -20);
    for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < WIDTH; ++j) {
            vect dir(170.0 - ((double)j / (double)rayPerPixel), 140.0 - ((double)i / (double)rayPerPixel), 0);
            colour clr = getAverageColour(orig, dir, objects, neighborPerPixel, rayPerPixel); // Adjust the number of samples
            image[i][j].free_add(clr);
        }
    }
}


int main(int argc, char* argv[]) {
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<Object*> objects;
    srand(time(0));
    //             center(x, y, z), radius, illumination, roughness, colour(r, g, b), isLight
    objects.push_back(new Sphere(vect(100, 400, 180), 100, 10, 0, colour(255, 255, 255), true)); //light
    objects.push_back(new Sphere(vect(200, 400, 180), 100, 10, 0, colour(255, 255, 255), true)); //light
    


    //                          center(x, y, z), radius, reflectivity, roughness, colour(r, g, b), isLight
    objects.push_back(new Sphere(vect(60, 40, 210), 35, 0.999, 0, colour(250, 255, 255), 0.9, 1.5, false)); // three spheres
    objects.push_back(new Sphere(vect(50, 130, 120), 45, 0.95, 0, colour(200, 200, 0), false));
    objects.push_back(new Sphere(vect(200, 70, 220), 70, 0.999, 0, colour(255, 200, 255), false));

    //                              top(x,y,z),         base1(x,y,z),       base2(x,y,z),       base3(x,y,z), reflectivity, roughness, colour(r, g, b), isLight
    objects.push_back(new Tetrahedron(vect(150, 150, 90), vect(180, 120, 70), vect(120, 120, 70), vect(145, 117, 50), 0.999, 0, colour(255, 255, 255), 0.8, 1.5, false));
    objects.push_back(new Tetrahedron(vect(170, 30, 160), vect(170, 30, 80), vect(130, 30, 160), 0.999, 0, colour(150, 150, 255), false));
    objects.push_back(new Tetrahedron(vect(170, 30, 80), vect(130, 30, 80), vect(130, 30, 160), 0.999, 0, colour(150, 150, 255), false));
    
    objects.push_back(new Sphere(vect(150, -5000, 120), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //bottom floor
    objects.push_back(new Sphere(vect(150, 5350, 120), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //ceiling
    objects.push_back(new Sphere(vect(150, 175, 5240), 5000, 0.99, 1.4, colour(200, 255, 200), false)); //front wall
    objects.push_back(new Sphere(vect(150, 175, -5021), 5000, 0.99, 1.4, colour(255, 255, 255), false)); //back wall
    objects.push_back(new Sphere(vect(5300, 175, 120), 5000, 0.99, 1.4, colour(255, 200, 200), false)); //left wall
    objects.push_back(new Sphere(vect(-5000, 175, 120), 5000, 0.99, 1.4, colour(200, 200, 255), false)); //right wall

    int normal = 0;
    int transp = 0;
    for (auto& object : objects) {
        if(object->transparency > 0.01){
            transp++;
        }
        else{
            normal++;
        }
    }
    
    int rpp = ray_per_pixel;
    int npp = neighbour_per_pixel; // Adjust to smoothen image

    unsigned long long int time=1;
    time *= (normal+(2*transp))*objects.size(); // O(n)
    time *= ray_per_pixel; // O(n^2)
    time *= ray_per_pixel;
    time *= MAXBOUNCE; // O(n^2)
    time *= MAXBOUNCE;
    time *= neighbour_per_pixel; // O(n^2)
    time *= neighbour_per_pixel;
    time /= TIMEMODIFIER; //obtained experimentally
    unsigned long long int expected_time = time;

    int EPOCH_LIMIT = INT32_MAX;
    if(TIMELIMIT<time){
        std::cout<<"Warning: Time limit not sufficient"<<std::endl;
        // exit(EXIT_FAILURE);
    }
    else{
        time = TIMELIMIT;
    }
    int time_remaining = TIMELIMIT;
    int time_min = time/60;
    time = time%60;
    int time_hour = time_min/60;
    time_min = time_min%60;

    
    printf("Please wait for less than %d hr, %d min, %lld sec\n",time_hour, time_min, time);

    int HEIGHT = 30 * ray_per_pixel;
    int WIDTH = 40 * ray_per_pixel;

    // Create array to hold color data
    colour** imageData = new colour*[HEIGHT];
    for (int i = 0; i < HEIGHT; ++i) {
        imageData[i] = new colour[WIDTH];
    }

    
    // Define parameters for parallelization
    int numThreads = std::thread::hardware_concurrency();
    numThreads = 4 * numThreads;
    int rowsPerThread = HEIGHT / numThreads;
    
    int EPOCHS=0;
    auto epoch_time = std::chrono::high_resolution_clock::now();
    auto epoch_duration = std::chrono::duration_cast<std::chrono::seconds>(epoch_time - begin); 
    bool time_calc_correct = true;
    bool checked = false;
    int time_modifier = -1;
    do{
        std::vector<std::thread> threads;
        for (int i = 0; i < numThreads; ++i) {
            int startRow = i * rowsPerThread;
            int endRow = (i == numThreads - 1) ? HEIGHT : startRow + rowsPerThread;
            threads.emplace_back([startRow, endRow, rpp, npp, &objects, imageData, WIDTH]() {
                renderImagePart(startRow, endRow, rpp, npp,  &objects, imageData, WIDTH);
            });
        }

        // Join threads
        for (std::thread& t : threads) {
            t.join();
        }
        ++EPOCHS;
        epoch_time = std::chrono::high_resolution_clock::now();
        epoch_duration = std::chrono::duration_cast<std::chrono::seconds>(epoch_time - begin); 
        if(!checked){
            int diff = abs((int)expected_time - epoch_duration.count());
            float percent_err = diff/(float)expected_time;
            if(percent_err>0.2){
                time_calc_correct = false;
                time_modifier =(int)(TIMEMODIFIER * (float)expected_time/ (float)epoch_duration.count());
                time = epoch_duration.count();
            }
            checked=true;
        }
    } 
    // while(epoch_duration.count() < TIMELIMIT);
    while(epoch_duration.count() < ((int)time_remaining-(int)time));


    // Write image data to file
    std::ofstream outputFile("intermediate.txt");
    for (int i = 0; i < HEIGHT; ++i) {
        for (int j = 0; j < WIDTH; ++j) {
            colour temp_clr;
            temp_clr.free_add(imageData[i][j]);
            temp_clr.r /= EPOCHS;
            temp_clr.g /= EPOCHS;
            temp_clr.b /= EPOCHS;

            if(temp_clr.r < 0){
                temp_clr.r = 0;
            }
            else if(temp_clr.r > 255){
                temp_clr.r = 255;
            }
            if(temp_clr.g < 0){
                temp_clr.g = 0;
            }
            else if(temp_clr.g > 255){
                temp_clr.g = 255;
            }
            if(temp_clr.b < 0){
                temp_clr.b = 0;
            }
            else if(temp_clr.b > 255){
                temp_clr.b = 255;
            }
            outputFile << std::fixed << std::setprecision(0) << temp_clr.r << " " << temp_clr.g << " " << temp_clr.b << ",";
        }
        outputFile << "\n";
    }
    outputFile.close();

    // Clean up
    for (int i = 0; i < HEIGHT; ++i) {
        delete[] imageData[i];
    }
    delete[] imageData;
    

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate duration
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - begin);

    // Output execution time
    std::cout << "Execution time: " << duration.count() / 60 << " minutes and " << duration.count() % 60 << " seconds. Epochs completed: " << EPOCHS << std::endl;
    if(!time_calc_correct){
        std::cout << "WARNING: Please change 'TIME_MODIFIER' value in v4.1.sh to: " << time_modifier << std::endl;
    }
    return 0;
}