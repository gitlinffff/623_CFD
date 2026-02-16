# AE623 Project-2: Euler Finite-Volume Solver

A 2D Euler solver for unstructured triangular meshes with Paraview export. All paths are relative; clone and run anywhere.

---

## Prerequisites

- **CMake** (≥ 3.10)
- **C++ compiler** with C++11 support (g++, clang++)

Check versions:
```bash
cmake --version
g++ --version   # or: clang++ --version
```

---
## Step 1: Download ZIP and Extract

```bash
# 1. Go to the GitHub repository page
# 2. Click "Code" → "Download ZIP"
# 3. Extract the ZIP file
# 4. Enter the extracted folder
cd <extracted-folder-name>
```

You should see:
```
.
├── CMakeLists.txt
├── main.cpp
├── mesh/
│   └── initial_mesh_3.gri
├── data/
├── tests/
└── ...
```



---

## Step 2: Create Build Directory

```bash
mkdir build
cd build
```

---

## Step 3: Configure with CMake

```bash
cmake ..
```

Expected output:
```
-- Configuring done
-- Generating done
-- Build files have been written to: ...
```

---

## Step 4: Compile

```bash
make
```

Or parallel build:
```bash
make -j4
```

Expected: `main` and test executables in `build/`.

---

## Step 5: Run the Solver

**Option A (recommended):**
```bash
make run
```
This runs `main` with working directory set to the project root, so `mesh/` and `data/` paths resolve correctly.

**Option B (manual):**
```bash
cd ..          # go to project root
./build/main   # run from project root
```

> ⚠️ Do **not** run `./main` from inside `build/` — the mesh path `mesh/initial_mesh_3.gri` is relative to the project root.

---

## Step 6: View Results in Paraview

1. Open [Paraview](https://www.paraview.org/).
2. File → Open → select `data/solution.vtu`.
3. Click **Apply** in the Properties panel.
4. (Optional) Filters → Alphabetical → **Cell Data to Point Data** for smoother rendering.
5. Use the color map to visualize `rho`, `u`, `v`, or `p`.

---

## Running Tests

From the `build/` directory:

```bash
make run_test        # testFlux
make run_test_mesh   # test_mesh
make run_test_verify # test_verify
```

Or run directly (from project root):
```bash
./build/tests/testFlux
./build/tests/test_mesh
./build/tests/test_verify
```

---

## Project Structure

| Path | Description |
|------|-------------|
| `mesh/` | `.gri` mesh files (Gambit format) |
| `data/` | Output VTU files (created automatically) |
| `build/` | Build artifacts (do not commit) |
| `tests/` | Unit tests |

---

## Rebuild After Code Changes

```bash
cd build
make
make run
```

Full clean rebuild:
```bash
cd build
rm -rf *
cmake ..
make
```

