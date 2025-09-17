#pragma once

constexpr int MAXBOUNCE = 20;            // Max number of reflections a ray can undergo (affects brightness and realism)
constexpr int TIMELIMIT = 600;           // Maximum render time in seconds

constexpr int ray_per_pixel = 20;        // Number of rays traced per pixel (higher → smoother image, slower render)
constexpr int neighbour_per_pixel = 7;   // Number of rays in neighborhood for anti-aliasing (higher → smoother edges)
constexpr int IMG_BITS = 16;             // Bit depth of output image (8 or 16 bits per channel)

constexpr bool DOF = true;               // Enable Depth of Field (true → objects at focal distance sharp, others blurred)
constexpr double FOCAL_LENGTH = 10.0;    // Camera focal length (controls field of view)
constexpr double FSTOP = 0.4;            // Camera aperture (smaller → shallower DOF, larger → deeper focus)
constexpr double FOCAL_DISTANCE = 160.0; // Distance at which camera is perfectly focused (used with DOF)

//Do not change the things below this line
constexpr int TIME_MODIFIER = 3529980;
constexpr int COLOR_MAX = (1 << IMG_BITS) - 1;
