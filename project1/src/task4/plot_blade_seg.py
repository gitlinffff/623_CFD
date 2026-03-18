import numpy as np
import matplotlib.pyplot as plt
from projection import read_blade_segments

blade_seg = read_blade_segments('../../data/bladeupper.txt', '../../data/bladelower.txt')

#x = blade_seg[:,:,0]; y = blade_seg[:,:,1]

plt.figure(figsize=(8, 6))
for seg in blade_seg:
    plt.plot(seg[:,0], seg[:,1], label='Original Blade Segment', color='red', alpha=0.5)


#plt.plot(x, y, label='Cubic Spline Interpolation', color='blue',alpha=0.7)

plt.title('Blade Profile Spline Interpolation')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.axis('equal')
plt.legend()
plt.grid()
plt.show()
