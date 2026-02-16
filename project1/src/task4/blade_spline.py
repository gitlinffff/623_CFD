import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Load blade coordinates from blade.txt
# Coordinates start and end at the same point (the trailing tip) 
coords = np.loadtxt('../../data/blade.txt')

# Calculate the distance between consecutive points
ds = np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))

# Create the parameter s (cumulative distance) starting at 0
s = np.insert(np.cumsum(ds), 0, 0)

# Create splines for x and y
spline_x = CubicSpline(s, coords[:, 0], bc_type='periodic')
spline_y = CubicSpline(s, coords[:, 1], bc_type='periodic')




if __name__ == "__main__":

    # Generate 500 smooth points along the blade
    s_fine = np.linspace(0, s[-1], 500)
    x_fine = spline_x(s_fine)
    y_fine = spline_y(s_fine)

    # plot the 500 smooth points
    plt.figure(figsize=(8, 6))
    plt.plot(x_fine, y_fine, label='Cubic Spline Interpolation', color='blue',alpha=0.7)
    plt.plot(coords[:, 0], coords[:, 1], 'o', label='Original Points', color='red',
             markersize=1)
    
    for i, (x, y) in enumerate(coords):
        plt.text(x, y, f'{i}', fontsize=6, ha='right', va='bottom')
    
    plt.title('Blade Profile Spline Interpolation')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.axis('equal')
    plt.legend()
    plt.grid()
    plt.show()
