# create_movie_and_save_clean_data.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Configuration ---
INPUT_FILE = "profileHelio.dat"          # Original, potentially messy data
SANITIZED_DATA_FILE = "data_clean.txt"  # File to save the cleaned data to
OUTPUT_MOVIE_FILE = "output.mp4"        # Final movie file
X_STEP = 0.1
FPS = 25

def load_and_sanitize_data(filename):
    """
    Reads a space-separated data file, handling rows with inconsistent lengths.
    It finds the minimum number of columns and truncates all rows to this length.
    """
    print(f"Loading data from '{filename}'...")
    
    raw_data_lines = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            try:
                row = [float(val) for val in line.strip().split()]
                raw_data_lines.append(row)
            except ValueError as e:
                print(f"Warning: Could not parse a line. Skipping. Error: {e}")

    if not raw_data_lines:
        raise ValueError("No valid data found in the file.")

    row_lengths = [len(row) for row in raw_data_lines]
    min_len = min(row_lengths)
    max_len = max(row_lengths)

    if min_len != max_len:
        print("------------------------------------------------------------------")
        print("! WARNING: Inconsistent number of columns found in the data file.")
        print(f"  Shortest row has {min_len} columns. Longest has {max_len} columns.")
        print(f"  Truncating all rows to the shortest length ({min_len} columns).")
        print("------------------------------------------------------------------")
        
        sanitized_data = [row[:min_len] for row in raw_data_lines]
    else:
        print("Data file has consistent row lengths.")
        sanitized_data = raw_data_lines

    return np.array(sanitized_data)

# --- 1. Load and Sanitize the Data ---
try:
    # `data` is now a clean, rectangular NumPy array
    data = load_and_sanitize_data(INPUT_FILE)
except (ValueError, FileNotFoundError) as e:
    print(f"Error: {e}")
    exit()

# --- 2. Save the Sanitized Data to a New File (NEW STEP) ---
try:
    print(f"Saving the cleaned and truncated data to '{SANITIZED_DATA_FILE}'...")
    # Use numpy.savetxt to write the array to a text file.
    # `delimiter=' '` ensures the values are space-separated.
    np.savetxt(SANITIZED_DATA_FILE, data, delimiter=' ')
    print("Clean data file saved successfully.")
except Exception as e:
    print(f"Error: Could not save the sanitized data file. Reason: {e}")


# --- 3. Prepare Data for Plotting ---
# The rest of the script proceeds as before, using the clean `data` array
time_values = data[:, 0]
y_data = data[:, 1:]

num_x_points = y_data.shape[1]
x_values = np.arange(0, num_x_points * X_STEP, X_STEP)

# --- 4. Set up the Plot ---
fig, ax = plt.subplots()
line, = ax.plot(x_values, y_data[0, :], color='g') # Green for go!

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(True)

y_min = np.min(y_data)
y_max = np.max(y_data)
ax.set_ylim(y_min * 1.1, y_max * 1.1)

time_template = "Time = %.2f s"
title = ax.set_title(time_template % time_values[0])

# --- 5. Define the Animation Update Function ---
def update(frame_number):
    line.set_ydata(y_data[frame_number, :])
    title.set_text(time_template % time_values[frame_number])
    return line, title

# --- 6. Create and Save the Animation ---
ani = animation.FuncAnimation(
    fig, update, frames=len(time_values), interval=1000 / FPS, blit=True
)

print(f"Saving animation to '{OUTPUT_MOVIE_FILE}'... (This may take a moment)")
try:
    ani.save(OUTPUT_MOVIE_FILE, writer='ffmpeg', fps=FPS, dpi=150)
    print("Movie created successfully!")
except FileNotFoundError:
    print("\nError: 'ffmpeg' not found.")
    print("Please install ffmpeg and ensure it's in your system's PATH.")