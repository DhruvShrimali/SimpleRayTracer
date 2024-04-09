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
#include <random>
#include <sstream>

#define MAXBOUNCE 15

void playSound() {
    #ifdef _WIN32
        // For Windows
        system("PowerShell -c (New-Object Media.SoundPlayer \"C:\\Windows\\Media\\tada.wav\").PlaySync();");
    #elif __APPLE__
        // For macOS
        system("afplay /System/Library/Sounds/Glass.aiff");
    #endif
}


float randomFloat(){
    return (float)(((rand() % 20)-10.0f)/10.0);
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
    colour(int i,int j,int k){
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
    float light_loss; //between 0 and 1
    float rough;
    colour clr;
    float transparency; // between 0 and 1
    float refr_index;
    bool light;

    Sphere(const vect& c, double r, float light_loss, float rog, colour col, float transp, float refr, bool lig) : 
        center(c), 
        radius(r), 
        light_loss(light_loss),
        rough(rog), 
        clr(col), 
        transparency(transp),
        refr_index(refr),
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
        reflectedDir = reflectedDir + (getRandDir()*rough);
        normalize(reflectedDir);
        return reflectedDir;
    }

    vect refract(vect& raydir, vect& normal, bool enteringMedium){
        double n = enteringMedium ? 1.0 / refr_index : refr_index;
        normalize(raydir);
        double cosI = dot(normal, raydir);
        double sinT2 = n * n * (1.0 - cosI * cosI);
        if (sinT2 > 1.0) {
            normal = normal*-1;
            return reflect(raydir, normal);
        }
        double cosT = sqrt(1.0 - sinT2);
        return (raydir*n) + (normal * ((n * cosI) - cosT));
    }
};

colour recursiveColourFinder(vect& orig, vect& dir, int bounces, std::vector<Sphere> *spheres){
    if(bounces >= MAXBOUNCE){
        return colour(0);
    }

    double tnear = INFINITY;
    Sphere* sphere = nullptr;

    std::vector<Sphere>::iterator itr;
    for (itr = spheres->begin(); itr < spheres->end(); itr++){
        double t0 = INFINITY, t1 = INFINITY;
        if (itr->intersect(orig, dir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &(*itr);
            }
        }
    }
    if (!sphere){
        return colour(0); // Default color if no intersection
    }
    if(sphere->light){
        colour temp = sphere->clr;
        temp.mult(sphere->light_loss);
        // if(bounces==0)
        //     printf("Light: %d %d %d\n", temp.r, temp.g, temp.b);
        return temp;
    }

    if(sphere->clr.r<=0 && sphere->clr.g<=0 && sphere->clr.b<=0){
        return colour(0);
    }

    // Intersection point
    vect intersection = orig + (dir * tnear);
    // Normal at intersection point
    vect normal = intersection - sphere->center;
    normalize(normal);

    colour finalColor = colour(0);
    
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
    colour reflectionColor = recursiveColourFinder(reflectOrig, reflectDir, bounces + 1, spheres);
    
    reflectionColor.mult(1.0 - sphere->transparency);
    reflectionColor.mult(sphere->clr);
    finalColor.add(reflectionColor);


    if (sphere->transparency > 0) {
        vect refractDir;
        vect refractOrig;
        if (enteringMedium) // when ray entering sphere
            refractOrig = intersection - (normal * 0.0001); // Start slightly inside the sphere
        else // when ray exiting sphere
            refractOrig = intersection + (normal * 0.0001); // Start slightly outside the sphere
        refractDir = sphere->refract(dir, normal, enteringMedium);
        normalize(refractDir);

        colour refractionColor = recursiveColourFinder(refractOrig, refractDir, bounces + 1, spheres);
        refractionColor.mult(sphere->transparency);
        refractionColor.mult(sphere->clr);
        finalColor.add(refractionColor);
    }

    // Compute final color considering reflection and refraction
    finalColor.mult(sphere->light_loss);
    return finalColor;
}


colour getAverageColour(vect& orig, vect& dir, std::vector<Sphere> *spheres, int numSamples, int k) {
    colour avgColour(0);
    double half = numSamples/(2 * k); 
    int divider = numSamples * numSamples;
    dir.x -= half; //get leftmost ray
    dir.y -= half; //get topmost ray
    k=k*numSamples;
    for (int i = 0; i < numSamples; ++i) { //a
        for (int j = 0; j < numSamples; ++j) { //b
            double x = dir.x + ((double)j/(double)k);
            double y = dir.y + ((double)i/(double)k);
            vect temp_dir(x, y, dir.z);
            vect dirSample = temp_dir - orig;
            normalize(dirSample); 
            colour clr = recursiveColourFinder(orig, dirSample, 0, spheres);
            avgColour.r += clr.r;
            avgColour.g += clr.g;
            avgColour.b += clr.b;
        }
    }

    avgColour.r /= divider;
    avgColour.g /= divider;
    avgColour.b /= divider;

    return avgColour;
}

// Function to render a portion of the image using thread
void renderImagePart(int startRow, int endRow, int rayPerPixel, int neighborPerPixel, std::vector<Sphere> &spheres, colour **image, int WIDTH, int iter) {
    vect orig(150, 125, -20);
    float divider = 1.0/(float)iter;
    for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < WIDTH; ++j) {
            vect dir(170.0 - ((double)j / (double)rayPerPixel), 140.0 - ((double)i / (double)rayPerPixel), 0);
            image[i][j]=colour(0);
            for(int k=0; k<iter; k++){
                colour clr = getAverageColour(orig, dir, &spheres, neighborPerPixel, rayPerPixel); // Adjust the number of samples
                image[i][j].r += clr.r;
                image[i][j].g += clr.g;
                image[i][j].b += clr.b;
            }
            image[i][j].mult(divider);
        }
    }
}


int main() {
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<Sphere> spheres;
    srand(time(0));
    //center(x, y, z), radius, illumination, roughness, colour(r, g, b), transparency, refr_index, isLight
    spheres.push_back(Sphere(vect(100, 400, 180), 100, 10, 0, colour(255, 255, 255), 0, 0, true)); //light
    spheres.push_back(Sphere(vect(200, 400, 180), 100, 10, 0, colour(255, 255, 255), 0, 0, true)); //light

    //center(x, y, z), radius, light_loss, roughness, colour(r, g, b), transparency, refr_index, isLight
    spheres.push_back(Sphere(vect(150, -5000, 120), 5000, 0.99, 1.4, colour(255, 255, 255), 0, 0, false)); //bottom floor
    spheres.push_back(Sphere(vect(150, 5350, 120), 5000, 0.99, 1.4, colour(255, 255, 255), 0, 0, false)); //ceiling
    spheres.push_back(Sphere(vect(150, 175, 5240), 5000, 0.99, 1.4, colour(200, 255, 200), 0, 0, false)); //front wall
    spheres.push_back(Sphere(vect(150, 175, -5031), 5000, 0.99, 1.4, colour(255, 255, 255), 0, 0, false)); //back wall
    spheres.push_back(Sphere(vect(5300, 175, 120), 5000, 0.99, 1.4, colour(200, 200, 255), 0, 0, false)); //left wall
    spheres.push_back(Sphere(vect(-5000, 175, 120), 5000, 0.99, 1.4, colour(255, 200, 200), 0, 0, false)); //right wall

    // Read sphere initialization values from a text file
    std::ifstream inputFile("sphere_initialization.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open sphere initialization file." << std::endl;
        return 1; // Return error code
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        if (tokens.size() != 12) {
            std::cerr << "Invalid number of values in sphere initialization line." << std::endl;
            continue; // Skip this line and proceed to the next
        }
        vect center(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
        double radius = std::stod(tokens[3]);
        float light_loss = std::stof(tokens[4]);
        float rough = std::stof(tokens[5]);
        int r = std::stoi(tokens[6]);
        int g = std::stoi(tokens[7]);
        int b = std::stoi(tokens[8]);
        float transparency = std::stof(tokens[9]);
        float refr_index = std::stof(tokens[10]);
        bool light = std::stoi(tokens[11]);
        spheres.push_back(Sphere(center, radius, light_loss, rough, colour(r, g, b), transparency, refr_index, light));
    }

    inputFile.close();
    std::cout<<"IO complete"<<std::endl;
    

    int ray_per_pixel = 20; // Adjust the number of samples per pixel
    int neighbour_per_pixel = 7; // Adjust to smoothen image
    int iterations = 2; // Number of times to iterate over all pixels

    unsigned long long int time=1;
    time *= iterations; // O(n)
    time *= spheres.size(); // O(n)
    time *= ray_per_pixel; // O(n^2)
    time *= ray_per_pixel;
    time *= MAXBOUNCE; // O(n^2)
    time *= MAXBOUNCE;
    time *= neighbour_per_pixel; // O(n^2)
    time *= neighbour_per_pixel;
    time /= 356125; //obtained experimentally
    int time_min = time/60;
    time = time%60;
    int time_hour = time_min/60;
    time_min = time_min%60;

    printf("Expected wait time: %d hr, %d min, %lld sec\n",time_hour, time_min, time);
    printf("Usually takes a fraction of above time if number of clear spheres are less\n");

    int HEIGHT = 30 * ray_per_pixel;
    int WIDTH = 40 * ray_per_pixel;

    // Create array to hold color data
    colour** imageData = new colour*[HEIGHT];
    for (int i = 0; i < HEIGHT; ++i) {
        imageData[i] = new colour[WIDTH];
    }
    
    // Define parameters for parallelization
    int numThreads = std::thread::hardware_concurrency();
    numThreads = 8 * numThreads;
    int rowsPerThread = HEIGHT / numThreads;

    std::vector<std::thread> threads;

    for (int i = 0; i < numThreads; ++i) {
        int startRow = i * rowsPerThread;
        int endRow = (i == numThreads - 1) ? HEIGHT : startRow + rowsPerThread;
        threads.emplace_back([startRow, endRow, ray_per_pixel, neighbour_per_pixel, &spheres, imageData, WIDTH, iterations]() {
            renderImagePart(startRow, endRow, ray_per_pixel, neighbour_per_pixel, spheres, imageData, WIDTH, iterations);
        });
    }

    // Join threads
    for (std::thread& t : threads) {
        t.join();
    }


    // Write image data to file
    std::ofstream outputFile("intermediate.txt");
    for (int i = 0; i < HEIGHT; ++i) {
        for (int j = 0; j < WIDTH; ++j) {
            outputFile << imageData[i][j].r << " " << imageData[i][j].g << " " << imageData[i][j].b << ",";
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
    std::cout << "Execution time: " << duration.count() / 60 << " minutes and " << duration.count() % 60 << " seconds." << std::endl;
    // system("say Image created");
    playSound();


    return 0;
}
