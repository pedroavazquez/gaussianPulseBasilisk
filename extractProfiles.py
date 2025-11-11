import numpy as np
from pathlib import Path

# Receive a list of filenames corresponding to the same value of 't'
def process_profiles_t(folder_path, t, filenames):
    """
    Process a list of filenames corresponding to the same value of 't'.
    This function can be customized to perform any operation on the files.
    """
    midpoints = []

    for file in filenames:
        try:
            name = str(folder_path) +str("/") + str(file)
            with open(name, 'r') as f:
                # 1. Read the entire file and split it into "blocks"
                #    A blank line in a file usually means two newline characters ('\n\n').
                #    .strip() removes any leading/trailing whitespace from the whole file content.
                content = f.read().strip()
                blocks = content.split('\n\n')
                
                if (len(blocks) > 1):
                    for i, block in enumerate(blocks):
                            # Split the block into individual lines.
                            lines = block.strip().split('\n')

                            # Ensure the block has exactly two lines.
                            if len(lines) != 2:
                                print(f"Warning: Skipping malformed block #{i+1}. Expected 2 lines, found {len(lines)}.")
                                continue

                            try:
                                # 3. Parse the coordinates from the two lines
                                #    .split() without arguments handles any amount of whitespace.
                                x1, y1 = map(float, lines[0].split())
                                x2, y2 = map(float, lines[1].split())

                                # 4. Calculate the midpoint
                                mid_x = (x1 + x2) / 2
                                mid_y = (y1 + y2) / 2

                                # 5. Add the midpoint as a tuple to our list
                                midpoints.append((mid_x, mid_y))
                                #print(f"  - Processed block #{i+1}: ({x1}, {y1}) and ({x2}, {y2}) -> Midpoint: ({mid_x}, {mid_y})")

                            except (ValueError, IndexError):
                                # This catches errors if a line doesn't contain two valid numbers.
                                print(f"Warning: Skipping block #{i+1} due to invalid coordinate format.")
                                continue

        except FileNotFoundError:
            print(f"Error: The file '{file_path}' was not found.")
            # Exit gracefully if the file doesn't exist
            exit()

    # 6. Sort the list of midpoints
    #    By default, sorted() on a list of tuples will sort by the first element,
    #    which is exactly what is needed.
    return sorted(midpoints)

# --- Configuration ---
# Use a Path object for the folder. '.' means the current directory.
folder_path = Path('gaussianPulse') 
file_pattern = 'profile-*'

# --- Script Logic ---
print(f"Searching for files in '{folder_path.resolve()}' matching pattern '{file_pattern}'...")

# 1. Use the .glob() method of the Path object.
#    This returns a generator of Path objects for each matching file.
matching_paths = folder_path.glob(file_pattern)

# 2. Get the .name attribute from each Path object and sort the list.
#    This is slightly cleaner than using os.path.basename().
sorted_filenames = sorted([p.name for p in matching_paths])

# 3. Convert the sorted list to a list of lists. Each sublist contains the filenames sharing the value of 't'.
from collections import defaultdict

# defaultdict(list) creates a dictionary where if a key is accessed for the
# first time, it will automatically be created with an empty list as its value.
grouped_files = defaultdict(list)

for filename in sorted_filenames:
    try:
        # Split the string by the '-' delimiter.
        # e.g., "profile-1.5-10.dat" -> ['profile', '1.5', '10.dat']
        parts = filename.split('-')
        
        # The float value is the second part (index 1)
        # We use this as the key for our dictionary.
        float_key = float(parts[1])
        
        # Add the full, original filename to the list for this key.
        grouped_files[float_key].append(filename)

    except (IndexError, ValueError):
        # This block will catch filenames that don't match the expected
        # "profile-float-integer" pattern, like "profile-invalid-format.txt".
        print(f"Warning: Skipping malformed filename: {filename}")

# File for Helio
fH = open("profileHelio.dat", "w")

# For each value of 't' (key) get a list of tuples with the coordinates of the profile, computed as the midpoints
# of the pairs of points in each file.
for key in sorted(grouped_files.keys()):
    
    # Compute the coordinates of the midpoints for each file corresponding to the same value of 't'
    midpoints = process_profiles_t(folder_path, key, grouped_files[key])
    midpoints.insert(0,(0.,1.))

    # Interpolate the midpoints to create a smooth profile for each value of 't'
    x_values = []
    y_values = []
    if len(midpoints) > 2:
        fH.write(f"{key:8.6f}")
        fH.write("\t")
        for point in midpoints:
            x_values.append(point[0])
            y_values.append(point[1])
    

        # Use numpy's interpolation function to create a smooth curve
        from scipy.interpolate import interp1d
        interpolator = interp1d(x_values, y_values, kind='linear', fill_value='extrapolate')


        # Generate a dense set of x values for the smooth curve
        #x_dense = np.linspace(min(x_values), max(x_values), num=100)
        x_dense = np.arange(min(x_values), max(x_values), 0.1)
        y_dense = interpolator(x_dense)

        # Print or save the interpolated profile
        if len(x_dense) > 0 and len(y_dense) > 0:
            for x, y in zip(x_dense, y_dense):
                fH.write(f"{y:8.6f}  ")
            #print(f"Interpolated profile for t={key}:")
            #nameOutput = f"profile-{key:<6.4f}.dat"
            #with open(nameOutput, 'w') as f:
            #    for x, y in zip(x_dense, y_dense):
                    #f.write(f"{x:.4e} {y:.4e}\n")
            #print(f"Profile saved to {nameOutput}")
            #for x, y in zip(x_dense, y_dense): print(f"{x:.4f} {y:.4f}")


        fH.write("\n")

fH.close()
