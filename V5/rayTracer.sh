#!/bin/bash

# Create ObjectFiles directory if it doesn't exist
mkdir -p ObjectFiles

# Compiler and flags
CXX=g++
CXXFLAGS="-std=c++17 -Ofast"

# Source files
SRC_FILES=(
    "Utils/random.cpp"
    "Utils/updateConfig.cpp"
    "Core/colourImplementation.cpp"
    "Core/vectorImplementation.cpp"
    "Objects/shape.cpp"
    "Objects/sphere.cpp"
    "Objects/tetrahedron.cpp"
    "rayTracer.cpp"
)

# Compile each source file into ObjectFiles/*.o
for src in "${SRC_FILES[@]}"; do
    obj="ObjectFiles/$(basename "${src%.*}.o")"
    echo "Compiling $src -> $obj"
    $CXX $CXXFLAGS -c "$src" -o "$obj"
done

# Collect all object files
OBJ_FILES=$(ls ObjectFiles/*.o)

# Link all object files into executable
EXE="rayTracer.exe"
echo "Linking .o files -> $EXE"
$CXX $CXXFLAGS $OBJ_FILES -o $EXE

# Run the executable
if [ -f "$EXE" ]; then
    echo "Running $EXE"
    ./$EXE
    python3 rayTracer.py
else
    echo "Executable not created!"
fi

