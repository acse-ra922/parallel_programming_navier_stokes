import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Specify the folder path containing the .dat files
folder_path = "../code/out"

# Regular expression pattern to extract the second number in the file name
pattern = re.compile(r"_\d+_(\d+)")

# Initialize dictionaries to store the grouped data for each variable
P_data = {}
v_data = {}
u_data = {}

# Iterate over the files in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith(".dat"):
        # Extract the second number from the file name using the pattern
        match = pattern.search(file_name)
        if match:
            second_number = int(match.group(1))

            # Determine the variable type based on the first character in the file name
            variable_type = file_name[0]

            # Read the data from the file
            data = np.loadtxt(os.path.join(folder_path, file_name))

            # Reshape the data array to have 2 dimensions
            data = np.reshape(data, (1, -1))

            # Append the data to the corresponding variable dictionary
            if variable_type == "P":
                if second_number in P_data:
                    if data.shape[1] > 0 and P_data[second_number].shape[1] > 0:
                        if data.shape[1] == P_data[second_number].shape[1]:
                            P_data[second_number] = np.concatenate((P_data[second_number], data))
                else:
                    P_data[second_number] = data
            elif variable_type == "v":
                if second_number in v_data:
                    if data.shape[1] > 0 and v_data[second_number].shape[1] > 0:
                        if data.shape[1] == v_data[second_number].shape[1]:
                            v_data[second_number] = np.concatenate((v_data[second_number], data))
                else:
                    v_data[second_number] = data
            elif variable_type == "u":
                if second_number in u_data:
                    if data.shape[1] > 0 and u_data[second_number].shape[1] > 0:
                        if data.shape[1] == u_data[second_number].shape[1]:
                            u_data[second_number] = np.concatenate((u_data[second_number], data))
                else:
                    u_data[second_number] = data

# Create an animation for each variable group
for second_number in sorted(P_data.keys()):
    # Combine the data for the current variable group
    P_combined_data = np.concatenate(list(P_data[second_number]))
    v_combined_data = np.concatenate(list(v_data[second_number]))
    u_combined_data = np.concatenate(list(u_data[second_number]))

    # Create a figure and axis for the animation
    fig, ax = plt.subplots()

    # Create an empty line to be updated in the animation
    line, = ax.plot([], [], lw=2)

    # Set the axis limits
    ax.set_xlim(np.min(P_combined_data), np.max(P_combined_data))
    ax.set_ylim(np.min(v_combined_data), np.max(v_combined_data))

    # Animation update function
def update(frame):
        n, m = P_combined_data.shape
        Y, X = np.meshgrid(np.arange(m), np.arange(n))  # Grid point indices
        U = u_combined_data.reshape(P_combined_data.shape)  # Reshape u and v to match the shape of P
        V = v_combined_data.reshape(P_combined_data.shape)
        ax.clear()  # Clear the current axes
        ax.imshow(P_combined_data, cmap='coolwarm')
        Q = ax.quiver(Y, X, U, V,units='xy', color='black')  # Draw the quiver plot
        return Q

    # Create the animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani = animation.FuncAnimation(fig, update_quiver, fargs=(ax,), frames=range(100), blit=False)
ani.save('animation.mp4', writer=writer)
