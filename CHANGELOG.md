# Simple Ray Tracer

All notable changes to this project will be documented here.


## V1
- Implemented basic ray-sphere intersection
- Output limited to boolean colors:<br>
&emsp;- true → white<br>
&emsp;- false → black<br>
- Intersection logic derived from geometric solution (<a href="https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html">see scratchapixel reference</a>)
![Ray-Sphere intersection geometric solution](https://www.scratchapixel.com/images/ray-simple-shapes/raysphereisect1.png?)<br>

How to run:
```bash
g++ -std=c++17 v1.cpp -o a.exe
./a.exe
```

#### V1.1
- Typecasted rays_per_pixel to double
- It creates a huge difference, which is apparent when rays_per_pixel is increased

How to run:
```bash
g++ -std=c++17 v1.1.cpp -o a.exe
./a.exe
```

#### V1.2
- Shows colours instead of plain black and white
- Outputs a PNG file for convenience 
- Especially useful when rays_per_pixel is greater than 3

How to run:
```bash
g++ -std=c++17 v1.2.cpp -o a.exe
./a.exe
python v1.2.py
```

## V2
- Returns light only when it reaches a light source through reflections (realistic behavior)
- Uses a recursive function to reflect light upon hitting the nearest object
- If the object is a light source, returns the light’s colour
- If the recursion depth exceeds MAXBOUNCE, returns black
- Returns black when the ray does not intersect any object
- Reflection logic was implemented with the help of <a href="https://math.stackexchange.com/questions/2334939/reflection-of-line-on-a-sphere">this answer</a>

How to run:
```bash
g++ -std=c++17 v2.cpp -o a.exe
./a.exe
python v2.py
```

#### V2.1
- Improves on V2 by removing image granularity
- For each ray in V2, traces neighbour_per_pixel × neighbour_per_pixel rays
- Averages their colour values and assigns it as the original ray’s colour
- Produces a slight blurring effect, making images appear smoother

How to run:
```bash
g++ -std=c++17 v2.1.cpp -o a.exe
./a.exe
python v2.1.py
```

#### V2.2
- Adds diffusion (roughness) to surfaces, similar to real-life objects
- Achieved by adding a random 3D vector to the perfect reflection vector
- The "rough" value in the object’s constructor controls the intensity of this effect
- Allows adjusting surfaces from perfectly smooth to highly rough/diffused

How to run:
```bash
g++ -std=c++17 v2.2.cpp -o a.exe
./a.exe
python v2.2.py
```

## V3
- Adds support for transparent objects
- Rays split into reflected and refracted rays when hitting a transparent surface
- Refracted-to-reflected ratio can be adjusted in the object’s constructor
- Refractive index parameter influences the refracted ray’s path
- Total internal reflection is handled by the recursion depth limit (MAXBOUNCE)
- Includes a shell script for easier execution<br>

<ins><b>⚠️ Note:</b></ins><br>

- This version is much slower due to nearly double recursion (reflection + refraction).
- Recommended to skip this version and use V3.1, which has the same functionality but is much faster.

How to run:
```bash
sh v3.sh
```

#### V3.1
- Parallelized ray-tracing calculations
- A standard number of threads is obtained, and a set of rows is assigned to each thread
- Threads cannot write directly to the output file (to avoid python parsing issues or misplaced data)
- A 2D array is used to temporarily store RGB values for each pixel
- Threads work concurrently and store results in the array
- After all child threads join, the 2D array is written to the intermediate .txt file
- The Python program then generates the final image from this file
- Performance boost: Reduced execution time from 50 minutes to 2.5 minutes (≈20× faster)

How to run:
```bash
sh v3.1.sh
```

#### V3.2
- Adapted to meet HCI project requirements
- Added a text file input for specifying sphere parameters
- Introduced an iterations parameter → runs the tracer multiple times and averages results for better image quality
- Parameters like MAXBOUNCE, neighbour_per_pixel, and iterations are tuned for balanced performance
- Plays a sound notification once rendering completes<br>

<ins><b>📄 Note:</ins></b><br>

- Refer to documentation.pdf for the correct input file format

How to run:
```bash
sh v3.2.sh
```

## V4
- Added support for tetrahedrons (and thus triangular faces)
- Enables transparent 3D objects using tetrahedrons
- Time modifier auto-corrects after the first run if error > 10%
- Improved time calculation logic to account for ray splitting in transparent objects (approximate but practical)
- Added option to set desired execution time (default = 180 seconds)

How to run:
```bash
sh v4.sh
```

## V5
- Added real camera-like functions (focal length, focal distance, f-stops) for more realistic rendering.
- Implemented autocorrection of time modifier for improved runtime estimation.
- Switched output from 8-bit to 16-bit colour images, increasing visual fidelity.
- Refactored codebase: split into multiple files with headers, added dedicated config file for parameters.
- Introduced a cleaner shell script for easier execution and workflow management.<br>

<ins><b>📄 Note:</ins></b><br>

- Lights’ intensity is now fixed. For scenes with non-transparent objects, you can decrease MAX_DEPTH to achieve almost the same quality in much less time.
- Decreasing MAX_DEPTH reduces illumination slightly; this can be compensated by increasing light intensity.
- Prefer using spheres over tetrahedrons when possible, as sphere intersection logic is much faster.

How to run:
```bash
sh v5.sh
```
