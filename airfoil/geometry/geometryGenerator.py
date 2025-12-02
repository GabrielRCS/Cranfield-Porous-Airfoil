import numpy as np
import math
from scipy.io import loadmat
from scipy.spatial import distance
from matplotlib.path import Path
from shapely.geometry import Point, Polygon

AoA_sweep = np.linspace(0,28,15)

for x in AoA_sweep:
    # Clear workspace (not needed in Python, but included for clarity)
    # %% GENERAL FLOW CONSTANTS - INPUTS
    scale_x = 16  # The length of the domain in the x direction (in NACA length units)
    scale_y = 8   # The length of the domain in the y direction
    lattice_resolution = 150
    meshRefinement = 2

    AoA_deg = x  # Angle of Attack in degrees; WARNING: negative for Palabos
    Pore_AoA_deg = 0  # Angle of Attack of pores in NACA referential
    Pore_relative_positions = []  # [0, 0.02, 0.03, 0.04, 0.05, 0.06]
    pore_width = 0.005

    # %% GENERAL FLOW CONSTANTS - OUTPUTS
    lx = lattice_resolution * scale_x  # number of cells in x-direction
    ly = lattice_resolution * scale_y  # number of cells in y-direction
    obst_x = lx // 4 + 1  # position of the NACA; (exact
    obst_y = ly // 2 + 3  # y-symmetry is avoided)
    obst_r = ly // 10 + 1  # This is cheating
    AoA = math.radians(AoA_deg)
    Pore_AoA_deg = Pore_AoA_deg + AoA_deg
    Pore_AoA = math.radians(Pore_AoA_deg)
    n_pores = len(Pore_relative_positions)

    # Define 2D rotation matrix
    RotAoA = np.array([[math.cos(AoA), -math.sin(AoA)], [math.sin(AoA), math.cos(AoA)]])
    RotPores = np.array([[math.cos(Pore_AoA), -math.sin(Pore_AoA)], [math.sin(Pore_AoA), math.cos(Pore_AoA)]])

    # Load geometry
    geom = np.loadtxt("NACA 0012 dat.txt") # Only change if you want to use a different airfoil

    # %% Define grid
    y, x = np.meshgrid(np.arange(1, ly + 1), np.arange(1, lx + 1))

    # %% Position the NACA appropriately
    center_x = np.mean(geom[:, 0])
    center_y = np.mean(geom[:, 1])
    geom = geom - [center_x, center_y]
    geom = np.dot(geom, RotAoA)
    geom = geom + [center_x, center_y]

    # Shift the NACA accordingly
    geom = geom * [lattice_resolution, lattice_resolution]
    geom = geom + [obst_x, obst_y]
    xgeom = geom[:, 0]
    ygeom = geom[:, 1]

    # Create shapely polygon  for the NACA geometry
    polygon = Polygon(zip(xgeom, ygeom))
    # Check if each point in the grid is inside the NACA polygon
    obst = np.array([polygon.contains(Point(xy)) for xy in np.column_stack((x.ravel(), y.ravel()))]).reshape(x.shape)

    # %% Add the porosity (commented out as in the original MATLAB code)
    # poreMask = np.zeros((lx, ly))
    # for i in range(n_pores):
    #     Pore = np.array([[-1, -pore_width/2], [-1, pore_width/2], [1, pore_width/2], [1, -pore_width/2]])
    #     Pore = np.dot(Pore, RotPores.T)
    #     Pore = Pore + [center_x, center_y + Pore_relative_positions[i]]
    #     Pore = Pore * [scaleFactor_x, scaleFactor_y] + [obst_x, obst_y]
    #     xp = Pore[:, 0]
    #     yp = Pore[:, 1]
    #     pore_path = Path(np.column_stack((xp, yp)))
    #     poreGeom = pore_path.contains_points(np.column_stack((x.ravel(), y.ravel()))).reshape(x.shape)
    #     poreMask += poreGeom
    # poreMask = 1 - poreMask
    # obst = obst * poreMask

    # %% Generate the geometry in a format appropriate for Palabos
    geometry = obst.astype(int)
    print("Original shape:", geometry.shape)
    padded = np.pad(geometry, pad_width=((0,1),(0,1)), mode='constant', constant_values=0)
    print("Padded shape:", padded.shape)
    geometry = padded.flatten()

    np.savetxt(f'geometry_{int(x)}.dat', geometry, fmt='%d', newline=' ')