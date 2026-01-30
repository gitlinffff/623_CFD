import numpy as np
from generate_matrices import generate_matrices
import sys

def verify_mesh(fname):
    print(f"Verifying mesh: {fname}")
    mats = generate_matrices(fname)
    I2E = mats['I2E']
    In = mats['In']
    B2E = mats['B2E']
    Bn = mats['Bn']
    In_len = mats['In_len']
    Bn_len = mats['Bn_len']
    Ne = 0
    if len(I2E) > 0:
        Ne = max(Ne, np.max(I2E[:, [0, 2]]))
    if len(B2E) > 0:
        Ne = max(Ne, np.max(B2E[:, 0]))
    print(f"Number of elements: {Ne}")
    element_errors = np.zeros((Ne, 2))
    # Process interior faces
    for i in range(len(I2E)):
        elemL = I2E[i, 0] - 1 # Convert to 0-based
        elemR = I2E[i, 2] - 1 # Convert to 0-based
        normal = In[i]
        length = In_len[i]
        element_errors[elemL] += normal * length
        element_errors[elemR] -= normal * length
    # Process boundary faces
    for i in range(len(B2E)):
        elem = B2E[i, 0] - 1 # Convert to 0-based
        normal = Bn[i]
        length = Bn_len[i]
        element_errors[elem] += normal * length
    # Compute magnitude of error for each element
    error_magnitudes = np.linalg.norm(element_errors, axis=1)
    max_error = np.max(error_magnitudes)
    error_elements = np.sum(error_magnitudes > 1e-12)
    
    print(f"Maximum error magnitude: {max_error}")
    print(f"Number of elements with errors: {error_elements}")
    # Check if errors are close to machine precision
    if max_error < 1e-12:
        print("Verification PASSED: Errors are within machine precision.")
    else:
        print("Verification FAILED: Errors exceed machine precision.")
    print("-" * 30)

if __name__ == "__main__":
    verify_mesh("test.gri")
    verify_mesh("initial_mesh_3.gri")
