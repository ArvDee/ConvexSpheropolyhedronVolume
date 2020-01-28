import math
import numpy

# Calculates the volume of a polyhedron. Needs the positions of the vertices and the indices of the
# faces. Will center vertices around the origin, and faces need not be triangles.
def volumePolyhedron(vertices, faces):
    # center vertices around (0,0,0)
    v_array = numpy.array(vertices)
    center = v_array.mean(axis=0)
    v_array = v_array - center
    # calculate volume by summing volume of tetrahedra of faces with origin
    vol = 0.0
    for f in faces:
        n = len(f)
        v2 = v_array[ f[0] ] # the pivot of the fan
        x2 = v2[0]
        y2 = v2[1]
        z2 = v2[2]
        for i in range(1,n-1): # divide into triangular fan segments
            v0 = v_array[ f[i] ]
            x0 = v0[0]
            y0 = v0[1]
            z0 = v0[2]
            v1 = v_array[ f[i+1] ]
            x1 = v1[0]
            y1 = v1[1]
            z1 = v1[2]
            # Add volume of tetrahedron formed by triangle and origin
            vol += math.fabs(x0 * y1 * z2 + x1 * y2 * z0 \
                         + x2 * y0 * z1 - x0 * y2 * z1 \
                         - x1 * y0 * z2 - x2 * y1 * z0)
    #print(vol/6.0)
    return vol/6.0

# Calculates the volume of a spheropolyhedron. Needs the positions of the polyhedron's vertices
# and the indices of its faces. Will center vertices around the origin, and faces need not be
# triangles.
def volumeSpheropolyhedron(vertices, faces, radius):
    if radius < 0:
        print("Fatal error in volumeSpheropolyhedron: negative rounding radius. Fix that first. Exiting.")
        sys.exit()
    vol_poly = volumePolyhedron(vertices, faces)
    if vol_poly == 0:
        return 4.0 * math.pi / 3.0 * radius**3
    elif vol_poly < 0:
        print("Fatal error in volumeSpheropolyhedron: polyhedron shape has negative volume.")
        print("Undefined behaviour may follow. Exiting instead.")
        sys.exit()
    # print(vol_poly)
    vol_sweep = 0.0

    # ===================== Cylinder segments =====================
    # preparation: find all edges between faces
    edges = []
    for f1_idx, f1 in enumerate(faces):
        for f2 in faces[f1_idx+1:]:
            e = []
            shared_v = 0
            for fv1 in f1:
                for fv2 in f2:
                    if fv1 == fv2:
                        e.append(fv1)
                        shared_v += 1
                if shared_v == 2:
                    edges.append(e)
                    shared_v = 0
                    edge = []
    # print(edges)
    # preparation: find the indices of the faces that share edges
    ef_map = []
    for e_idx, e in enumerate(edges):
        shared_f = 0
        ef_map.append([])
        for f_idx, f in enumerate(faces):
            shared_v = 0
            for fv in f:
                for ev in e:
                    if fv == ev:
                        shared_v += 1
                if shared_v == 2:
                    ef_map[e_idx].append(f_idx)
                    shared_f += 1
                    break
            if shared_f == 2:
                break
    # print(ef_map)
    # Calculate the volume of cylinders around the edges
    for e_idx, e in enumerate(edges):
        e_faces = ef_map[e_idx]
        e_vec = numpy.array(vertices[e[0]]) - numpy.array(vertices[e[1]])
        e_length = numpy.linalg.norm(e_vec)
        # normal of first face
        f1 = faces[ e_faces[0] ]
        f1e1 = numpy.array(vertices[f1[1]]) - numpy.array(vertices[f1[0]])
        f1e2 = numpy.array(vertices[f1[2]]) - numpy.array(vertices[f1[0]])
        n1 = numpy.cross(f1e1,f1e2)
        n1 = n1 / numpy.sqrt(numpy.sum(n1**2)) # normalize
        # normal of second face
        f2 = faces[ e_faces[1] ]
        f2e1 = numpy.array(vertices[f2[1]]) - numpy.array(vertices[f2[0]])
        f2e2 = numpy.array(vertices[f2[2]]) - numpy.array(vertices[f2[0]])
        n2 = numpy.cross(f2e1,f2e2)
        n2 = n2 / numpy.sqrt(numpy.sum(n2**2)) # normalize
		# Angle between the normals divided by 2 pi gives the fraction of
		# the cylinder around the edge. Its volume is then easily calculated.
        angle = math.acos(numpy.dot(n1,n2))
        vol_sweep += math.pi * (radius**2) * e_length * angle / (2.0 * math.pi)
    # print(vol_poly + vol_sweep)

    # ===================== Sphere segments =====================
    # First make a list of the faces of each vertex
    vf_map = []
    for v_idx, v in enumerate(vertices):
        vf_map.append([])
        for f_idx, f in enumerate(faces):
            for fv_idx in f:
                if fv_idx == v_idx:
                    vf_map[v_idx].append(f_idx)
    # print(vf_map)

	# We now generate a list of points for each vertex: for each connected face we make a new point
	# on the new sphereswept surface which is on the border of the spherical segment. This point is
	# the intersection of the cylinder segments we made above, and also one of the vertices of the
	# extended faces.
    for v_idx, v in enumerate(vertices):
        segment_vertices = []
        for f_idx in vf_map[v_idx]:
            f = faces[f_idx]
            fe1 = numpy.array(vertices[f[1]]) - numpy.array(vertices[f[0]])
            fe2 = numpy.array(vertices[f[2]]) - numpy.array(vertices[f[0]])
            n = numpy.cross(fe1,fe2)
            n = n / numpy.sqrt(numpy.sum(n**2)) # normalize
            new_v = numpy.array(vertices[v_idx]) + radius * n
            segment_vertices.append(new_v)

        # We now construct tetrahedra OABC, where O is the old vertex position and ABC are taken from
        # the set of points we just constructed. For each of these (non-overlapping) tetrahedra, we
        # find the solid angle of the ABC triangle and use that to calculate the fraction of a full
        # sphere that this segment occupies. By iterating over all these tetrahedra we get the volume
        # of the entire spherical segment sweeping this vertex.
        v_O = vertices[v_idx]
        v_A = segment_vertices[0] - v_O
        l_A = numpy.linalg.norm(v_A)
        solid_angle = 0.0;
        for idx in range(0,len(segment_vertices)-1):
            v_B = numpy.array(segment_vertices[idx]) - v_O
            l_B = numpy.linalg.norm(v_B)
            v_C = numpy.array(segment_vertices[idx+1]) - v_O
            l_C = numpy.linalg.norm(v_C)
            m = numpy.array([v_A, v_B, v_C])
            triple_product = abs(numpy.linalg.det(m))
            if triple_product == 0.0: continue
            solid_angle = 2.0 * math.atan2(triple_product, (l_A*l_B*l_C +
                                                            l_C * numpy.dot(v_A,v_B) +
                                                            l_A * numpy.dot(v_B,v_C) +
                                                            l_B * numpy.dot(v_C,v_A)))
        vol_sweep += (4.0/3.0) * math.pi * (radius**3) * solid_angle / (4.0 * math.pi)
    # print(vol_poly + vol_sweep)

	# ===================== Face segments =====================
	# Finally, we need to calculate the increase in volume from moving the faces
    # outward. This we can do by first calculating the 2D area of the face
    # (decomposed into triangles), and then simply multiplying by the swept radius.
    for f in faces:
        nv = len(f)
        for i in range(1,nv-1): # divide into triangular fan segments
            e1 = numpy.array(vertices[f[i]]) - numpy.array(vertices[f[0]])
            e2 = numpy.array(vertices[f[i+1]]) - numpy.array(vertices[f[0]])
            f_area = 0.5 * numpy.linalg.norm(numpy.cross(e1,e2))
            vol_sweep += radius * f_area

    # print(vol_poly + vol_sweep)
    return vol_poly + vol_sweep


# Test on some polyhedra:

# A cube
cubeVertices = [
    [-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5], [-0.5, 0.5, -0.5], [-0.5, 0.5, 0.5],
    [0.5, -0.5, -0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5], [0.5, 0.5, 0.5]
]
cubeFaces = [[7, 3, 1, 5], [7, 5, 4, 6], [7, 6, 2, 3], [3, 2, 0, 1], [0, 2, 6, 4], [1, 0, 4, 5]]
cubeVol = volumePolyhedron(vertices=cubeVertices, faces=cubeFaces)
print("Volume of a cube with edge length 1:",cubeVol)
spherocubeVol = volumeSpheropolyhedron(vertices=cubeVertices, faces=cubeFaces, radius=1.0)
print("Volume of a spherocube with edge length 1, rounding radius 1:",spherocubeVol)
# Should be the sum of vols of a unit cube, a radius 1 sphere, 3 radius 1 length 1 cylinders and 6 width 1 length 1 square platelets, i.e.:
print("Should be: ",4.0*math.pi/3.0*(1.0**3.0) + 3 * math.pi*(1.0**2.0)*1.0 + 6*(1.0**2.0)*1.0 + 1.0**3.0)

# A tetrahedron
tetraVertices = [
    [0,  0, math.sqrt(2.0/3.0) - 1.0/(2.0*math.sqrt(6.0))],
    [-1.0/(2.0*math.sqrt(3.0)), -1.0/2.0, -1.0/(2.0*math.sqrt(6.0))],
    [-1.0/(2.0*math.sqrt(3.0)), 1.0/2.0, -1.0/(2.0*math.sqrt(6.0))],
    [1.0/math.sqrt(3.0), 0.0, -1.0/(2.0*math.sqrt(6.0))]
]
tetraFaces = [[1, 2, 3], [2, 1, 0], [3, 0, 1], [0, 3, 2]]
tetraVol = volumePolyhedron(vertices=tetraVertices, faces=tetraFaces)
print("Volume of a tetrahedron with edge length 1:",tetraVol)
spherotetraVol = volumeSpheropolyhedron(vertices=tetraVertices, faces=tetraFaces, radius=1.0)
print("Volume of a spherotetrahedron with edge length 1, rounding radius 1:",spherotetraVol)
