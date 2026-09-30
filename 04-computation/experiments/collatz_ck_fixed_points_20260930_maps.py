#!/usr/bin/env python3
"""collatz_ck_fixed_points_20260930_maps.py -- the map reading of the C_k fixed points (thirteenth note, part B, addendum).
 (a) The paper's edge-removal step at level one: I_2 := C_5(I_1) minus the edge [p,q]; is C_5(I_2) an induced subgraph of
     C_5^2(I_1)?  (Gervacio-Maehara-Ramos, proof of Theorem 3.5.)
 (b) Square tori Cay(Z_n, {+-a, +-b}) = Z^2/L: C_4 = the dual map = the same circulant when the faces are the only induced
     squares; checked for the census list and for the non-example C_9(1,3) (extra triangles) and K_(4,4) = C_8(1,3).
 (c) The honeycomb torus over Z[w]/(13): degree 3, 26 vertices; C_6(honeycomb) = the triangular torus P_13, then fixed:
     the {6,3} map is C_6-preperiodic into the {3,6} Eisenstein fixed point, as the dodecahedron is C_5-preperiodic into I.
 (d) The 4-antiprism C_8(1,2): its eight induced pentagons wind once around the ring; C_5(C_8(1,2)) = C_8(2,3), the image
     under the multiplier 3.
"""
from collatz_pentagon_gadget_20260930 import icosahedron, pentagon_graph, hats, find_icosahedra, induced_subgraph_iso, graph_from_edges
from collatz_ck_fixed_points_20260930 import induced_cycles, cycle_operator, is_isomorphic, circulant, rhombus_graph

# (a)
I = icosahedron()
tad = None
for x in range(12):
    for y in I[x]:
        for z in I[x] & I[y]:
            for t in I[y] - I[x] - I[z] - {x, z}:
                tad = [x, y, z, t]; break
            if tad: break
        if tad: break
    if tad: break
I1 = [set(a) for a in I] + [set()]
for v in tad:
    I1[v].add(12); I1[12].add(v)
C1, pents1 = pentagon_graph(I1)
hat_pents = [i for i, p in enumerate(pents1) if 12 in p]
print("(a) C_5(I_1): %d vertices; pentagons through the hat: %s, adjacent: %s" % (len(C1), hat_pents, hat_pents[1] in C1[hat_pents[0]]))
I2 = [set(a) for a in C1]
p, q = hat_pents
I2[p].discard(q); I2[q].discard(p)
C_I2, _ = pentagon_graph(I2)
C2, _ = pentagon_graph(C1)
print("    C_5(I_2) has %d vertices, C_5^2(I_1) has %d; C_5(I_2) is an induced subgraph of C_5^2(I_1): %s" % (len(C_I2), len(C2), induced_subgraph_iso(C_I2, C2)))
icos = find_icosahedra(C_I2, range(len(C_I2)), max_found=5)
print("    C_5(I_2) contains an induced icosahedron with %s hats (the paper: four)" % ([len(hats(C_I2, Ic)) for Ic in icos][:3]))

# (b)
print("(b) square tori: C_4(Cay(Z_n, {a, b})) equals the same circulant when the faces are the only induced squares")
for n, a, b in ((7, 1, 2), (12, 1, 4), (13, 1, 5), (14, 1, 4), (15, 1, 4), (16, 1, 6), (17, 1, 4), (9, 1, 3), (8, 1, 3), (10, 1, 3), (11, 1, 3), (20, 1, 4), (25, 1, 7)):
    G = circulant(n, [a, b, n - a, n - b])
    sq = induced_cycles(G, 4, limit=4 * n)
    faces = sorted(min((x, (x + a) % n, (x + a + b) % n, (x + b) % n), (x, (x + b) % n, (x + a + b) % n, (x + a) % n)) for x in range(n))
    # canonical form used by induced_cycles: min vertex first, then smaller neighbour: recompute faces canonically
    def canon(c):
        m = min(c); i = c.index(m); r = c[i:] + c[:i]
        return min(r, (r[0],) + r[:0:-1])
    faces = sorted(set(canon(f) for f in faces))
    C4, cyc = cycle_operator(G, 4, limit=4 * n)
    fixed = (C4 is not None and len(cyc) == n and is_isomorphic(C4, G))
    print("   C_%d(%d,%d): induced squares %s, faces are all of them: %s, C_4-fixed: %s" % (n, a, b, len(sq) if len(sq) <= 4 * n else ">%d" % (4 * n), sorted(sq) == faces if len(sq) <= 4 * n else False, fixed))

# (c) honeycomb over Z[w]/(13): vertices (x, s), s in {0,1}; (x,0) ~ (x,1), (x+1,1), (x+1+w,1) with w -> 3 mod 13
n = 13; w = 3
edges = []
for x in range(n):
    for d in (0, 1, (1 + w) % n):   # A_x ~ B_x, B_(x+1), B_(x+1+w): a triangle of the B-lattice (1, 1+w = -w^2, w are units)
        edges.append((x, n + (x + d) % n))
H = graph_from_edges(2 * n, edges)
print("(c) honeycomb torus on %d vertices: degree %s, induced hexagons %d" % (2 * n, sorted(set(len(a) for a in H)), len(induced_cycles(H, 6))))
C6H, _ = cycle_operator(H, 6)
P13 = circulant(13, [1, 3, 4, 9, 10, 12])
print("    C_6(honeycomb) iso P_13: %s; C_6^2 iso P_13: %s -> the {6,3} torus is C_6-preperiodic into the Eisenstein fixed point" % (is_isomorphic(C6H, P13), is_isomorphic(cycle_operator(C6H, 6)[0], P13)))

# (d)
A8 = circulant(8, [1, 2, 6, 7])
C5A, pents = cycle_operator(A8, 5)
print("(d) 4-antiprism C_8(1,2): induced pentagons %s" % pents)
print("    C_5(C_8(1,2)) iso C_8(1,2): %s; equals C_8(2,3) on the pentagon labels P_i -> i: %s" % (
    is_isomorphic(C5A, A8), all(((j - i) % 8 in {2, 3, 5, 6}) == (j in C5A[i]) for i in range(8) for j in range(8) if i != j)))
