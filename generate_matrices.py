import numpy as np
import matplotlib.pyplot as plt
import sys

def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    print("Reading vertices!")
    V = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    print("Finished reading vertices.")
    # read boundaries
    print("Reading boundaries!")
    NB = int(f.readline())
    B = []; Bname = []
    B2 = {}
    for i in range(NB):
        s = f.readline().split(); Nb = int(s[0]); name = s[2]; Bname.append(name)
        for n in range(Nb):
            s = f.readline().strip().split()
            n1 = int(s[0]) - 1
            n2 = int(s[1]) - 1
            Bi = np.array([n1, n2])
            faceID = "%d %d"%(min(n1, n2), max(n1, n2))
            B2[faceID] = name
        B.append(Bi)
    print("Finished reading boundaries.")
    # read elements
    print("Reading elements!")
    Ne0 = 0; E = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0==0) else np.concatenate((E,Ei), axis=0)
        Ne0 += ne
    print("Finished reading elements.")
    # read periodic groups - store doubly as dict for easier access (e.g. if 1 and 4 are pair, map 1 to 4 AND 4 to 1)
    NPG = int(f.readline().split()[0])
    PG = []; PGtype = []
    for i in range(NPG):
        s = f.readline().split(); Npgn = int(s[0]); PGtype.append(s[1])
        nmap = {}
        for n in range(Npgn):
            s = f.readline().split(); n1 = int(s[0]) - 1; n2 = int(s[1]) - 1
            nmap[n1] = n2
            nmap[n2] = n1
        PG.append(nmap)
    f.close()
    Mesh = {'V': V, 'E': E, 'B': B, 'Bname': Bname, 'B2': B2, 'PG': PG, 'PGtype': PGtype}
    return Mesh

def generate_matrices(fname):
    print(f"Generating matrices from file '{fname}'")
    mesh = readgri(fname)
    # mesh.V     = node coordinates
    # mesh.E     = list of elements as triangles
    # mesh.B     = list of lists of boundary edges, separated by type
    # mesh.Bname = list of names of boundary edge types
    # print(mesh)

    # I2E: a mapping from interior faces (edges in 2D) to elements
    I2E = []

    # B2E: a mapping from boundary faces to elements and boundary groups
    B2E = []

    # In: normal vectors for interior faces
    In = []

    # Bn: normal vectors for boundary faces
    Bn = []

    # Area: the area of each element
    Area = []

    # keep track of visited edges through a hash table
    visited = {}

    # append a copy of the first column to the triangles matrix - this will
    # simplify some of the logic in the loop
    tris = np.array(mesh['E'])
    tris = np.append(tris, np.atleast_2d(tris[:,0]).T, axis=1)

    # NOTE: this logic is currently explicitly written to handle triangular
    # elements only. It is possible it could successfully be abstracted to
    # elements with arbitrary numbers of faces by replacing all instances of '3'
    # with a call to tris.shape[1], but I have not tested this yet.

    for i in range(tris.shape[0]): # i = **0-indexed** element number
        pt1, pt2, pt3 = [mesh['V'][tris[i][k]] for k in range(3)]
        Area.append(0.5*np.cross(pt2-pt1, pt3-pt1))
        for j in range(3): # j = **0-indexed** local node number
            endpts = np.sort([tris[i,j], tris[i,j+1]])
            faceID = "%d %d"%(endpts[0], endpts[1]) # face ID string serves as hash table key
            if faceID in mesh['B2']:
                # print(f"face ({faceID}) is a boundary face of type {mesh['B2'][faceID]}!")
                # check if periodic boundary
                periodic = False
                for pg in mesh['PG']:
                    if (endpts[0] in pg) and (endpts[1] in pg):
                        # print("this edge is in a periodic group!")
                        periodic = True
                        break
                if periodic:
                    # print("this is a periodic boundary!")
                    other_endpts = np.sort([pg[endpts[0]], pg[endpts[1]]])
                    other_faceID = "%d %d"%(other_endpts[0], other_endpts[1])
                    if other_faceID in visited:
                        # print("second encounter")
                        # face was encountered already
                        elemL, faceL = visited[other_faceID]
                        elemR = i+1 # convert i to 1-indexed element number
                        faceR = 3 if j == 0 else j # convert j to 1-indexed local face number
                        I2E.append([elemL, faceL, elemR, faceR])
                        del visited[other_faceID] # remove the face from the hash table
                        # calculate normal and append to In
                        edge = mesh['V'][endpts[1]] - mesh['V'][endpts[0]]
                        edge /= np.linalg.norm(edge)
                        In.append([edge[1], -edge[0]])
                    else:
                        # print("first encounter")
                        # first time encountering this face
                        elemL = i+1 # convert i to 1-indexed element number
                        faceL = 3 if j == 0 else j # convert j to 1-indexed local face number
                        visited[faceID] = (elemL, faceL) # add face to hash table
                else:
                    # print("this is a true boundary!")
                    elem = i+1 # convert i to 1-indexed element number
                    face = 3 if j == 0 else j # convert j to 1-indexed local face number
                    bgroup = mesh['Bname'].index(mesh['B2'][faceID]) + 1
                    # print(f"{mesh['B2'][faceID]} -> bgroup = {bgroup}")
                    B2E.append([elem, face, bgroup])
                    # calculate normal and append to Bn
                    edge = mesh['V'][tris[i,j+1]] - mesh['V'][tris[i,j]]
                    edge /= np.linalg.norm(edge)
                    Bn.append([edge[1], -edge[0]])
            else:
                # print(f"face ({faceID}) is an interior face!")
                if faceID in visited:
                    # face was encountered already
                    elemL, faceL = visited[faceID]
                    elemR = i+1 # convert i to 1-indexed element number
                    faceR = 3 if j == 0 else j # convert j to 1-indexed local face number
                    I2E.append([elemL, faceL, elemR, faceR])
                    del visited[faceID] # remove the face from the hash table
                    # calculate normal and append to In
                    edge = mesh['V'][endpts[1]] - mesh['V'][endpts[0]]
                    edge /= np.linalg.norm(edge)
                    In.append([edge[1], -edge[0]])
                else:
                    # first time encountering this face
                    elemL = i+1 # convert i to 1-indexed element number
                    faceL = 3 if j == 0 else j # convert j to 1-indexed local face number
                    visited[faceID] = (elemL, faceL) # add face to hash table

    return I2E, B2E, In, Bn, Area   # TODO ask if these should be returned in a single dict/object instead of separately

if __name__ == "__main__":
    args = sys.argv[1:]
    fname = args[0] if len(args) > 0 else "test.gri"
    I2E, B2E, In, Bn, Area = generate_matrices(fname)
    print("=== Output ===")
    print("I2E:")
    print(I2E)
    print("In:")
    print(In)
    print("B2E:")
    print(B2E)
    print("Bn:")
    print(Bn)
    print("Area:")
    print(Area)