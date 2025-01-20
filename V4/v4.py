import numpy as np
import os
from PIL import Image

# Your code to generate image_data

# Read RGB values from the text file
with open('intermediate.txt', 'r') as file:
    lines = file.readlines()

# Initialize an empty list to store the image data
image_data = []

# Parse each line to extract RGB values
for line in lines:
    rgb_values = line.strip().split(',')
    row_data = []
    for rgb in rgb_values:
        if not rgb:
            continue  # Skip if the string is empty
        r, g, b = map(int, rgb.split())
        row_data.append([r, g, b])  # Append RGB values to the row data
    image_data.append(row_data)  # Append the row data to the image data

# Convert the list of lists to a NumPy array
image_array = np.array(image_data, dtype=np.uint8)

# Create PIL Image from the numpy array
image_pil = Image.fromarray(image_array)

image_pil.save('output.png')

# Remove the temporary file
os.remove('intermediate.txt')
