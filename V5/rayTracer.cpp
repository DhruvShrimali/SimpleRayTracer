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
#include "Utils/random.h"
#include "Utils/updateConfig.hpp"
#include "Core/vectorImplementation.hpp"
#include "Core/colourImplementation.hpp"
#include "Objects/shape.hpp"
#include "Objects/sphere.hpp"
#include "Objects/tetrahedron.hpp"
#include "config.hpp"


double closest_object(vect& orig, vect& dir, std::vector<Shape*> *objects){
    double tnear = INFINITY;

    for (auto itr = objects->begin(); itr < objects->end(); itr++) {
        double t0 = INFINITY, t1 = INFINITY;
        Shape* obj = *itr; // Dereference to get the Object pointer

        // Check if the object is a sphere
        if(obj->isSphere) {
            int intersect_check = static_cast<Sphere*>(obj)->intersect(orig, dir, t0, t1);
            if (intersect_check) {
                if (t0 < 0) t0 = t1; // Use t1 if t0 is behind the ray origin
                if (t0 < tnear) {
                    tnear = t0; // Update the closest intersection distance
                }
            }
        } 
        // The object is a tetrahedron
        else{
            int intersect_check = static_cast<Tetrahedron*>(obj)->intersect(orig, dir, t0, t1);
            if (intersect_check) {
                if (t0 < tnear) {
                    tnear = t0; // Update the closest intersection distance
                }
            }
        }
    }
    return tnear;
}

colour recursive(vect& orig, vect& dir, int bounces, std::vector<Shape*> *objects){
    if(bounces >= MAXBOUNCE){
        return colour(0);
    }

    double tnear = INFINITY;
    Shape* object = nullptr;
    int latest_hit_face = 0;

    for (auto itr = objects->begin(); itr < objects->end(); itr++) {
        double t0 = INFINITY, t1 = INFINITY;
        Shape* obj = *itr; // Dereference to get the Object pointer
        

        // Check if the object is a sphere
        if(obj->isSphere) {
            int intersect_check = static_cast<Sphere*>(obj)->intersect(orig, dir, t0, t1);
            if (intersect_check) {
                if (obj->light) {
                    colour temp = obj->clr;
                    temp.mult(obj->energy_retention);
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
                    temp.mult(obj->energy_retention);
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

    if (object->clr.r <= 0 && object->clr.g <= 0 && object->clr.b <= 0) 
    {
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

    // If refraction is not considered or reflection is more probable

    // Calculate color recursively
    colour reflectionColor = recursive(reflectOrig, reflectDir, bounces + 1, objects);

    reflectionColor.mult(1.0 - object->transparency);
    reflectionColor.mult(object->clr);

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
        reflectionColor.free_add(refractionColor);
    }

    // Compute final color considering energy loss
    reflectionColor.mult(object->energy_retention);
    return reflectionColor;
}


colour getAverageColour(vect& orig, vect& dir, std::vector<Shape*> *objects, int neighborPerPixel, int rayPerPixel) {
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
                colour clr(0);
                if(DOF){
                    double objectDistance = closest_object(orig, dirSample, objects);

                    double blurAmount = std::abs(objectDistance - FOCAL_DISTANCE);

                    double A = FOCAL_LENGTH * FOCAL_LENGTH/ FSTOP;

                    double bottom = abs(FOCAL_LENGTH - FOCAL_DISTANCE) * objectDistance;

                    // Circle of confusion radius
                    double coc = A * blurAmount / bottom;

                    // Use coc to jitter
                    double randomOffsetX = guassianRandomFloat() * coc;
                    double randomOffsetY = guassianRandomFloat() * coc;

                    vect temp_dir_dof(x+randomOffsetX, y+randomOffsetY, dir.z);
                    vect dirSample_dof = temp_dir_dof - orig;
                    normalize(dirSample_dof); 
                    clr = recursive(orig, dirSample_dof, 0, objects);
                }
                else{
                    clr = recursive(orig, dirSample, 0, objects);
                }
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
void renderImagePart(int startRow, int endRow, int rayPerPixel, int neighborPerPixel, std::vector<Shape*> *objects, colour **image, int WIDTH) {
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

    std::vector<Shape*> objects;
    srand(time(0));
    //             center(x, y, z), radius, illumination, roughness, colour(r, g, b), isLight
    objects.push_back(new Sphere(vect(100, 400, 180), 100, 1, 0, colour(255, 255, 255), true)); //light
    objects.push_back(new Sphere(vect(200, 400, 180), 100, 1, 0, colour(255, 255, 255), true)); //light
    

    //                          center(x, y, z), radius, energy retention, roughness, colour(r, g, b), isLight
    objects.push_back(new Sphere(vect(60, 40, 210), 35, 0.999, 0, colour(255, 255, 255), 0.9, 1.5, false)); // three spheres
    objects.push_back(new Sphere(vect(50, 130, 120), 45, 0.95, 0, colour(200, 200, 0), false));
    objects.push_back(new Sphere(vect(200, 70, 220), 70, 0.999, 0, colour(255, 200, 255), false));

    //                              top(x,y,z),         base1(x,y,z),       base2(x,y,z),       base3(x,y,z), energy retention, roughness, colour(r, g, b), isLight
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
    time *= (normal+(2*transp))*objects.size(); // Crude estimation: O(n)
    time *= ray_per_pixel; // O(n^2)
    time *= ray_per_pixel;
    time *= MAXBOUNCE; // O(n^2)
    time *= MAXBOUNCE;
    time *= neighbour_per_pixel; // O(n^2)
    time *= neighbour_per_pixel;
    time /= TIME_MODIFIER; //obtained experimentally
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
    numThreads = 2 * numThreads;
    int rowsPerThread = HEIGHT / numThreads;
    
    int EPOCHS=0;
    auto epoch_time = std::chrono::high_resolution_clock::now();
    auto epoch_duration = std::chrono::duration_cast<std::chrono::seconds>(epoch_time - begin); 
    bool time_calc_correct = true;
    bool checked = false;
    int time_modifier = -1;
    int expected_epochs = 0;
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
            if(percent_err>0.1){
                time_calc_correct = false;
                time_modifier =(int)(TIME_MODIFIER * (float)expected_time/ (float)epoch_duration.count());
                time = epoch_duration.count();
            }
            checked=true;
            time = epoch_duration.count();
            expected_epochs = (int)(time_remaining / time);
        }
        std::cout<<"Epoch "<<EPOCHS<<" / "<<expected_epochs<<" completed\n";
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

            temp_clr.trim();

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
        updateConfig("config.hpp", "TIME_MODIFIER", std::to_string(time_modifier));
    }
    return 0;
}