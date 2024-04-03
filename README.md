# Simple Ray Tracer
I decided to make this ray tracer as part of my Human Computer Intercation course. Ray tracing, a foundational method in computer graphics, is employed to create lifelike images by simulating the trajectory of light rays interacting with scene elements.<br>
The original objectives were:
<ol>
<li>Implementing a simple ray-tracing algorithm suitable for a simple scene.</li>
<li>Rendering a static scene with a room, movable shapes (constrained to a couple balls for now), a fixed light source, and a camera.</li>
<li>Generating an image render depicting the scene with basic lighting effects and object interactions.</li>
</ol>

Required python librares for png image generation (V1.2 onwards):
```
pip install pillow
pip install numpy
```

Sample Render:<br>


## V1
This tracer only checks if light rays are intersecting with spheres. It initially only has bool colours true and false representing white and black.<br>
The intersection logic was made using the figure given below:<br><br>
![Ray-Sphere intersection geometric solution](https://www.scratchapixel.com/images/ray-simple-shapes/raysphereisect1.png?)<br>
Credit: <a href="https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html">scratchapixel</a>

How to run:
```
g++ -std=c++17 v1.cpp -o a.exe
./a.exe
```

#### V1.1
This version was created to highligh the importance of typecasting rays_per_pixel. Typecasting created huge difference, as can be seen in the output, wich wil be more apparent if rays_per_pixel is increased.

How to run:
```
g++ -std=c++17 v1.1.cpp -o a.exe
./a.exe
```

#### V1.2
This version improves upon the previous versions by showing colours instead of plain black and white. This also outputs a png file, making it more convenient as rays_per_pixel become larger than 3.

How to run:
```
g++ -std=c++17 v1.2.cpp -o a.exe
./a.exe
python v1.2.py
```

## V2
This is a realistic tracer that returns light only when it reaches a light source through reflections, which is what happens in real life. I uses a recursive function that reflects light on hitting the nearest object. If the object is a light, it returns the light's colour, or if the number of bounces (depth of recursion) exceeds the limit MAXBOUNCE, it returns black colour. It also returns black in case the ray does not intersect any object.<br>

Reflection logic was implemented with the help of this answer: <a href="https://math.stackexchange.com/questions/2334939/reflection-of-line-on-a-sphere">Link to answer</a>

How to run:
```
g++ -std=c++17 v2.cpp -o a.exe
./a.exe
python v2.py
```

#### V2.1
This version improves on V2 by removing granularity of images generated by V2. For each ray in V2, it traces neighbour_per_pixel X neighbour_per_pixel rays, and averages their colour values and returns it as the original ray's colour. This achieves a slight blurring effect, which makes the images appear smoother.

How to run:
```
g++ -std=c++17 v2.1.cpp -o a.exe
./a.exe
python v2.1.py
```

#### V2.2
This version adds diffusion or roughness to surfaces, similar to objects found in real life. It achieves this effect by adding a random 3D vector to the perfect reflection vector. This effect can be increased or reduced by changing the "rough" value in constructor of object.

How to run:
```
g++ -std=c++17 v2.2.cpp -o a.exe
./a.exe
python v2.2.py
```

## V3
This version further improves upon the previous tracer by allowing transparent objects. Transparent objects now lead to the ray splitting into a reflected and refracted ray. The ratio of refracted to reflected can be adjusted in the constructor. The refractive index parameter also affects the refrcated ray's path.
Total internal reflection is handled automatically because there is a limit on the number of bounces. Once the limit is reached, it returns a black colour.
This version also has a shell script, which makes executing it easier.<br>
<ins><b>Note:</b></ins> This version is quite slower than previous versions due to almost double recursion at each step (reflection + refaction). I would recommend skipping this version and run version 3.1. It has the same functionalities, but way faster.

How to run:
```
sh v3.sh
```

#### V3.1
The objective of this version is to parallelize the ray-tracing calculations. I first get some standard number of threads and assign few rows to each thread. Since like previous versions, each thread cannot output directly to the output file (otherwise it would lead to data not conforming to python program's reading syntax or wrong data at wrong location). Thus a 2D array is created to temporarily store RGB values for each pixel. This allows the threads work concurrently. After all child threads join, the 2D array is printed to the intermediate txt file, and the same python program generates the corresponding image. <br>
This method has allowed me to speed up the time taken for image generation from 50 minutes to 2.5 minutes - almost a 20 fold decrease in execution time.

How to run:
```
sh v3.1.sh
```

