#include "shape.hpp"

Shape::Shape(float energy_retention, float rough, colour col, float transp, float refr, bool lig, bool sphere):
    energy_retention(energy_retention), 
    rough(rough), 
    clr(col), 
    transparency(transp), 
    refr_index(refr), 
    light(lig), 
    isSphere(sphere) {}