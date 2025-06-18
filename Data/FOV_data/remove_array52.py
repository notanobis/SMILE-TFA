import os
import re

# Directory containing your VOTable files
directory = "./Data/FOV_data/01.08.2025-03.08.2025-Res0.5/mhd"

# Regular expression pattern to remove arraysize="52"
remove_arraysize_52_pattern = re.compile(r'\s*arraysize="125"')

# Iterate over all XML files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".xml"):
        file_path = os.path.join(directory, filename)

        # Read file contents
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()

        # Remove arraysize="52" if it exists
        updated_content = remove_arraysize_52_pattern.sub("", content)

        # If changes were made, overwrite the file
        if updated_content != content:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(updated_content)
            print(f"Updated {filename}: Removed arraysize='125'")
        else:
            print(f"No changes needed for {filename}")