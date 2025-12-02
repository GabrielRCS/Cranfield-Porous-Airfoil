import numpy as np
import pyvista as pv
import h5py
import matplotlib.pyplot as plt

# --- Load your 2D boolean mask --- and 
mask_2d = np.loadtxt("./geometry_flat.dat")   # array of 0s and 1s

    
# --- Wrap as a VTK image grid ---
#nx, ny = 241, 41
nx, ny = 16*150+1, 8*150+1
mask_2d = mask_2d.reshape((nx, ny))
grid = pv.ImageData()
grid.dimensions = (nx + 1, ny + 1, 1)
grid.spacing = (1, 1, 1)      # pixel size; adjust as needed
grid.origin = (0, 0, 0)


    

# --- Add the mask data ---
grid.cell_data["mask"] = mask_2d.flatten(order="F")



# --- Optional preview ---
grid.plot(scalars="mask", cmap="gray")

