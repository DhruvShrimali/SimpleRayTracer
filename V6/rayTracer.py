import os
import re
import shutil
import numpy as np
from PIL import Image
import cv2

# -----------------------------
# Configuration and Constants
# -----------------------------
CONFIG_FILE = "config.hpp"
INTERMEDIATE_FILE = "intermediate.txt"
OUTPUT_FILE = "output.png"
VAR_NAME = "IMG_BITS"


# -----------------------------
# Utility Functions
# -----------------------------
def get_img_bits(config_path: str, var_name: str) -> int:
    """Read config.hpp and return the value of a constexpr int variable."""
    with open(config_path, "r") as f:
        content = f.read()
    
    pattern = rf"constexpr\s+int\s+{var_name}\s*=\s*(\d+)\s*;"
    match = re.search(pattern, content)
    
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"{var_name} not found in {config_path}")


def read_image_data(file_path: str) -> np.ndarray:
    """Read RGB values from a text file and return as a NumPy array."""
    image_data = []
    with open(file_path, 'r') as f:
        for line in f:
            rgb_values = line.strip().split(',')
            row_data = []
            for rgb in rgb_values:
                if not rgb:
                    continue
                r, g, b = map(int, rgb.split())
                row_data.append([r, g, b])
            image_data.append(row_data)
    return np.array(image_data)


def save_image(image_array: np.ndarray, bit_depth: int, output_file: str):
    """Save the image array to a file based on bit depth."""
    if bit_depth == 8:
        img = Image.fromarray(image_array.astype(np.uint8))
        img.save(output_file)
    elif bit_depth == 16:
        # OpenCV expects BGR for saving
        image_bgr = cv2.cvtColor(image_array.astype(np.uint16), cv2.COLOR_RGB2BGR)
        cv2.imwrite(output_file, image_bgr)
    else:
        raise ValueError(f"Unsupported bit depth: {bit_depth}")


def cleanup(files_to_remove: list, dirs_to_remove: list):
    """Remove specified files and directories."""
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
    for directory in dirs_to_remove:
        if os.path.exists(directory):
            shutil.rmtree(directory)


# -----------------------------
# Main Execution
# -----------------------------
def main():
    try:
        img_bits = get_img_bits(CONFIG_FILE, VAR_NAME)
        # print(f"{VAR_NAME} = {img_bits}")
    except ValueError as e:
        print(e)
        return

    image_array = read_image_data(INTERMEDIATE_FILE)
    
    save_image(image_array, img_bits, OUTPUT_FILE)

    # Cleanup temporary files
    cleanup([INTERMEDIATE_FILE], ["ObjectFiles"])


if __name__ == "__main__":
    main()
