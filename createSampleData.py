import numpy as np

# --- Parameters ---
num_time_steps = 100  # Number of lines in the file (tn)
num_x_points = 200    # Number of y-values per line
x_step = 0.1
output_filename = "data.txt"

# Create the x-axis (this is not saved in the file, but used for generation)
x = np.arange(0, num_x_points * x_step, x_step)

# Create a file to write to
with open(output_filename, 'w') as f:
    # Loop through each time step
    for i in range(num_time_steps):
        t = i * 0.05  # Current time value
        
        # Generate a function, e.g., a traveling, decaying sine wave
        # y(x,t) = exp(-t) * sin(2*pi*(x - t))
        y_values = np.exp(-t) * np.sin(2 * np.pi * (x - t))
        
        # Create the line to write: t y0 y1 y2 ...
        line_data = np.hstack([t, y_values])
        
        # Write the line to the file, with space-separated values
        f.write(" ".join(map(str, line_data)) + "\n")

print(f"Sample file '{output_filename}' created with {num_time_steps} time steps.")