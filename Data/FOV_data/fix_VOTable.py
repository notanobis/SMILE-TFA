# This script adds the attribute "arraysize" to the FIELD entry for "values" in all VOTable files in a directory.
# This fixes the format of the initial tables, and each row can now be read as an array of 52 values.

import os
import re


# Directory containing your VOTable files
directory = "./Data/FOV_data/01.05.2026-03.05.2026/latep"

# Regular expression pattern to find the FIELD entry for "values"
field_pattern = re.compile(r'(<FIELD name="values".*?)(?<!arraysize="32")(?= width=")')

# Regular expression pattern to fix leading spaces in azimuth values
azimuth_pattern = re.compile(r'(<PARAM name="azimuth".*?value=")\s+(-?\d+\.\d+)', re.DOTALL)


# Iterate over all XML files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".xml"):
        file_path = os.path.join(directory, filename)

        # Read file contents
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()

        # Check if the field already has arraysize="52", if not, add it
        updated_content = field_pattern.sub(r'\1 arraysize="32"', content)

        # Remove leading spaces in the first azimuth value
        updated_content = azimuth_pattern.sub(r'\1\2', updated_content)
        # If changes were made, overwrite the file
        if updated_content != content:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(updated_content)
            print(f"Updated {filename}")
        else:
            print(f"No changes needed for {filename}")
