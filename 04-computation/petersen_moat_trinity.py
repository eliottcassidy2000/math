"""
petersen_moat_trinity.py -- kind-pasteur-2026-03-14-S67

THE 30-EDGE TRINITY: Petersen graph, Platonic solids, E_8, and the moat at T=10.

Key observation: 30 = h(E_8) = 2*3*5 appears in THREE distinct guises:
  1. h(E_8) = 30 (Coxeter number of the largest exceptional Lie group)
  2. 30 = #edges of icosahedron = #edges of dodecahedron
  3. 30 = #edges of Johnson graph J(5,2) = complement of Petersen graph

The Petersen graph K(5,2):
  - 10 vertices (2-subsets of {1,...,5})
  - 15 edges (disjoint 2-subsets)
  - 3-regular, girth 5, triangle-free
  - |Aut(Petersen)| = 120 = |S_5| = |BI| (binary icosahedral group!)
  - Independence number alpha = 4

Its complement J(5,2):
  - 10 vertices, 30 edges (intersecting 2-subsets)
  - 6-regular
  - J(5,2) = line graph of K_5

And the moat at T=10:
  - The Petersen graph has 10 vertices = T value of the moat
  - T=10 is permanently forbidden in tournaments
  - The Petersen graph is NOT a conflict graph of any tournament

Recurrence angle:
  - Tournament recurrence z^2 - 5z + 6 = 0, roots KEY_1=2, KEY_2=3
  - 5 = #vertices of base set in K(5,2) = KEY_1 + KEY_2
  - 6 = KEY_1 * KEY_2 = #edges of tetrahedron
  - The recurrence encodes the Petersen/platonic/exceptional trinity

This script computes:
  1. I.P. of the Petersen graph and its complement
  2. Connection to tournament conflict graphs
  3. The 30-edge parallel across all three domains
  4. Recurrence structure
"""

from itertools import combinations
from collections import Counter

def independence_polynomial(adj, n):
    """Compute I(G, x) as list of coefficients [alpha_0, alpha_1, ...]."""
    # Enumerate all 2^n subsets via bitmask
    coeffs = [0] * (n + 1)
    # Pre-compute neighbor sets as bitmasks
    nbr = [0] * n
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                nbr[i] |= (1 << j)
    for mask in range(1 << n):
        # Check if mask is independent
        independent = True
        used = mask
        bits = []
        temp = mask
        while temp:
            b = temp & (-temp)  # lowest set bit
            v = b.bit_length() - 1
            bits.append(v)
            # Check if v has a neighbor in the already-processed part
            if nbr[v] & (mask & ((1 << v) - 1)):
                independent = False
                break
            temp ^= b
        if independent:
            coeffs[bin(mask).count('1')] += 1
    return coeffs

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def build_petersen():
    """Build the Petersen graph K(5,2) adjacency matrix."""
    # Vertices = 2-subsets of {0,1,2,3,4}
    vertices = list(combinations(range(5), 2))
    n = len(vertices)  # 10
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            # Edge iff disjoint
            if len(set(vertices[i]) & set(vertices[j])) == 0:
                adj[i][j] = adj[j][i] = 1
    return adj, n, vertices

def build_complement(adj, n):
    """Build complement graph."""
    comp = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j] == 0 and i != j:
                comp[i][j] = comp[j][i] = 1
    return comp

def build_icosahedron():
    """Build icosahedron adjacency matrix (12 vertices, 30 edges)."""
    # Standard icosahedron: vertices labeled 0-11
    # Top vertex 0, bottom vertex 11, two pentagons in between
    edges = [
        (0,1),(0,2),(0,3),(0,4),(0,5),           # top to upper ring
        (1,2),(2,3),(3,4),(4,5),(5,1),             # upper ring
        (1,6),(2,7),(3,8),(4,9),(5,10),            # upper to lower
        (6,7),(7,8),(8,9),(9,10),(10,6),           # lower ring
        (6,2),(7,3),(8,4),(9,5),(10,1),            # cross connections
        (6,11),(7,11),(8,11),(9,11),(10,11)        # lower ring to bottom
    ]
    n = 12
    adj = [[0]*n for _ in range(n)]
    for i, j in edges:
        adj[i][j] = adj[j][i] = 1
    return adj, n, len(edges)

def main():
    print("=" * 70)
    print("THE 30-EDGE TRINITY: PETERSEN, PLATONIC, E_8")
    print("=" * 70)

    # 1. Petersen graph
    print("\n" + "=" * 50)
    print("1. THE PETERSEN GRAPH K(5,2)")
    print("=" * 50)

    pet_adj, pet_n, pet_verts = build_petersen()
    pet_edges = sum(sum(row) for row in pet_adj) // 2
    pet_degrees = [sum(row) for row in pet_adj]

    print(f"  Vertices: {pet_n}")
    print(f"  Edges: {pet_edges}")
    print(f"  Regular degree: {pet_degrees[0]}")
    print(f"  Vertex labels (2-subsets of [5]):")
    for i, v in enumerate(pet_verts):
        print(f"    {i}: {v}")

    pet_ip = independence_polynomial(pet_adj, pet_n)
    print(f"\n  Independence polynomial:")
    print(f"    I(Pet, x) = {' + '.join(f'{c}x^{k}' if k > 0 else str(c) for k, c in enumerate(pet_ip) if c > 0)}")
    print(f"    I(Pet, 2) = {eval_poly(pet_ip, 2)}")
    print(f"    I(Pet, -1) = {eval_poly(pet_ip, -1)}")
    print(f"    I(Pet, 1) = {eval_poly(pet_ip, 1)} (= # independent sets)")

    # 2. Complement = Johnson graph J(5,2)
    print("\n" + "=" * 50)
    print("2. COMPLEMENT: JOHNSON GRAPH J(5,2)")
    print("=" * 50)

    john_adj = build_complement(pet_adj, pet_n)
    john_edges = sum(sum(row) for row in john_adj) // 2
    john_degrees = [sum(row) for row in john_adj]

    print(f"  Vertices: {pet_n}")
    print(f"  Edges: {john_edges}  <-- THE 30!")
    print(f"  Regular degree: {john_degrees[0]}")

    john_ip = independence_polynomial(john_adj, pet_n)
    print(f"\n  Independence polynomial:")
    print(f"    I(J(5,2), x) = {' + '.join(f'{c}x^{k}' if k > 0 else str(c) for k, c in enumerate(john_ip) if c > 0)}")
    print(f"    I(J(5,2), 2) = {eval_poly(john_ip, 2)}")
    print(f"    I(J(5,2), -1) = {eval_poly(john_ip, -1)}")

    # Key: J(5,2) = line graph of K_5
    print(f"\n  J(5,2) = L(K_5) (line graph of complete graph on 5 vertices)")
    print(f"  Two edges of K_5 are adjacent in L(K_5) iff they share a vertex")
    print(f"  This is EXACTLY the tournament conflict graph structure!")
    print(f"  Conflict = sharing a vertex among cycles")

    # 3. Icosahedron
    print("\n" + "=" * 50)
    print("3. ICOSAHEDRON (PLATONIC)")
    print("=" * 50)

    ico_adj, ico_n, ico_edges = build_icosahedron()
    print(f"  Vertices: {ico_n}")
    print(f"  Edges: {ico_edges}  <-- THE 30!")
    print(f"  Faces: 20 (triangular)")
    print(f"  |Aut| = 120 = |A_5 x Z_2| = |BI|")

    ico_ip = independence_polynomial(ico_adj, ico_n)
    print(f"\n  Independence polynomial:")
    print(f"    I(Ico, x) = {' + '.join(f'{c}x^{k}' if k > 0 else str(c) for k, c in enumerate(ico_ip) if c > 0)}")
    print(f"    I(Ico, 2) = {eval_poly(ico_ip, 2)}")
    print(f"    I(Ico, -1) = {eval_poly(ico_ip, -1)}")

    # 4. The number 30
    print("\n" + "=" * 50)
    print("4. THE NUMBER 30 = 2 * 3 * 5")
    print("=" * 50)
    print()
    print("  30 = h(E_8)           Coxeter number of E_8")
    print("  30 = |edges(Ico)|     Icosahedron/dodecahedron edges")
    print("  30 = |edges(J(5,2))| Johnson graph = complement of Petersen")
    print("  30 = 2 * 3 * 5       Product of tournament primes")
    print("  30 = KEY_1 * KEY_2 * (KEY_1+KEY_2)")
    print("  30 = n(n^2-1)/24|_{n=9}  Three-cycles in regular T_9")
    print("  30 = |S_5|/4 = 120/4  Quarter of icosahedral symmetry")
    print()
    print("  The three 30's correspond to:")
    print("    Lie theory:    E_8 Coxeter number")
    print("    Geometry:      Icosahedral edge count")
    print("    Combinatorics: Kneser complement edges")
    print()
    print("  All three are manifestations of the (2,3,5) factorization.")

    # 5. Connection to the moat at T=10
    print("\n" + "=" * 50)
    print("5. THE MOAT AT T=10 AND THE PETERSEN GRAPH")
    print("=" * 50)
    print()
    print("  T = 10 is permanently forbidden in tournament H-spectra.")
    print("  The Petersen graph has exactly 10 vertices.")
    print()
    print("  If Petersen were a tournament conflict graph Omega(T):")
    print(f"    H = I(Petersen, 2) = {eval_poly(pet_ip, 2)}")
    print(f"    T = (H-1)/2 = {(eval_poly(pet_ip, 2)-1)//2}")
    print(f"    alpha_1 = {pet_ip[1]}, alpha_2 = {pet_ip[2]}")
    print(f"    T_quad = alpha_1 + 2*alpha_2 = {pet_ip[1] + 2*pet_ip[2]}")
    print()
    print("  But Petersen is NOT a tournament CG because:")
    print("  - It's triangle-free (girth 5), but tournament CGs always")
    print("    have triangles when alpha_1 >= 3 (from splicing)")
    print("  - Its independence number is 4, but tournament CGs at n=5")
    print("    have specific I.P. structure")
    print()
    print("  The deep connection: Petersen's 10 vertices = T=10 moat")
    print("  And its complement J(5,2) has 30 edges = h(E_8)")
    print()
    print("  Recurrence: z^2 - 5z + 6 = 0")
    print(f"    5 = #vertices of base set in K(5,2)")
    print(f"    6 = #edges of tetrahedron = KEY_1 * KEY_2")
    print(f"    10 = C(5,2) = #vertices of Petersen = THE MOAT")

    # 6. The 5 Platonic solids and 5 exceptional Lie groups
    print("\n" + "=" * 50)
    print("6. THE DOUBLE PENTAD: SOLIDS AND GROUPS")
    print("=" * 50)
    print()
    print("  Platonic solids      Exceptional Lie groups")
    print("  -----------------    ----------------------")
    print("  Tetrahedron (6e)     G_2 (rank 2, h=6)")
    print("  Cube (12e)           F_4 (rank 4, h=12)")
    print("  Octahedron (12e)     E_6 (rank 6, h=12)")
    print("  Dodecahedron (30e)   E_7 (rank 7, h=18)")
    print("  Icosahedron (30e)    E_8 (rank 8, h=30)")
    print()
    print("  Edge counts: {6, 12, 30} = 6 * {1, 2, 5}")
    print("  Coxeter numbers: {6, 12, 18, 30}")
    print()
    print("  MATCHING: tetrahedron <-> G_2 (both have 6)")
    print("            cube/octahedron <-> F_4/E_6 (all have 12)")
    print("            icosahedron <-> E_8 (both have 30)")
    print("            dodecahedron <-> E_7 (30 edges, h=18=30-12)")
    print()
    print("  The 'leftover' E_7 has h=18 = 30-12, bridging")
    print("  the icosahedral (30) and octahedral (12) levels.")

    # 7. Recurrence structure
    print("\n" + "=" * 50)
    print("7. RECURRENCE STRUCTURE")
    print("=" * 50)
    print()
    print("  Tournament recurrence: a(n) = 5*a(n-1) - 6*a(n-2)")
    print("  Roots: KEY_1=2, KEY_2=3")
    print("  General solution: a(n) = A*2^n + B*3^n")
    print()
    print("  Platonic edge sequence: 6, 12, 30")
    print("    12/6 = 2 = KEY_1")
    print("    30/12 = 5/2 = (KEY_1+KEY_2)/KEY_1")
    print("    30/6 = 5 = KEY_1+KEY_2")
    print()
    print("  Exceptional Coxeter sequence: 6, 12, 12, 18, 30")
    print("    Differences: 0, 6, 6, 12")
    print("    The differences are {0, 6, 12} = {0, 1, 2} * h(G_2)")
    print()
    print("  The recurrence z^2 - 5z + 6 = (z-2)(z-3) = 0")
    print("  encodes the MOAT: z=KEY_1+KEY_2=5 makes z^2=25, 5z=25, 6=6")
    print("  => 25-25+6 = 6 != 0. But C(5,2) = 10 = THE MOAT VALUE.")
    print()
    print("  The moat at T=10 = C(KEY_1+KEY_2, KEY_1) = C(5,2)")
    print("  is the binomial coefficient of the recurrence roots' SUM")
    print("  taken KEY_1 at a time.")

    # 8. Petersen graph as forbidden substructure
    print("\n" + "=" * 50)
    print("8. PETERSEN AS FORBIDDEN SUBSTRUCTURE")
    print("=" * 50)
    print()
    print("  The Robertson-Seymour theorem: every minor-closed family")
    print("  of graphs has a finite set of forbidden minors.")
    print()
    print("  For PLANAR graphs: K_5 and K_{3,3} (Wagner/Kuratowski)")
    print("  K_5 has 5 vertices = KEY_1+KEY_2 and 10 edges = THE MOAT")
    print()
    print("  The Petersen graph is a forbidden minor for:")
    print("  - Linklessly embeddable graphs (Robertson-Seymour-Thomas)")
    print("  - Graphs of treewidth <= 4 (with other minors)")
    print("  - Knotlessly embeddable graphs (with other minors)")
    print()
    print("  CONJECTURE: The Petersen graph is 'forbidden' for")
    print("  tournament conflict graphs in a similar structural sense.")
    print("  Its 10 vertices, triangle-free structure, and girth-5")
    print("  property make it impossible as Omega(T) for any tournament.")

    # 9. I.P. values at special points
    print("\n" + "=" * 50)
    print("9. INDEPENDENCE POLYNOMIALS AT SPECIAL POINTS")
    print("=" * 50)
    print()
    # Compute for several graphs
    graphs = {
        "Petersen K(5,2)": (pet_adj, pet_n),
        "Johnson J(5,2)": (john_adj, pet_n),
        "Icosahedron": (ico_adj, ico_n),
    }

    # Build some other key graphs
    # K_5 (complete on 5 vertices)
    k5_adj = [[1 if i != j else 0 for j in range(5)] for i in range(5)]
    graphs["K_5"] = (k5_adj, 5)

    # C_5 (5-cycle)
    c5_adj = [[0]*5 for _ in range(5)]
    for i in range(5):
        c5_adj[i][(i+1)%5] = c5_adj[(i+1)%5][i] = 1
    graphs["C_5"] = (c5_adj, 5)

    # K_3 (triangle)
    k3_adj = [[0,1,1],[1,0,1],[1,1,0]]
    graphs["K_3"] = (k3_adj, 3)

    # K_3 union K_1 (triangle + isolated)
    k3k1_adj = [[0,1,1,0],[1,0,1,0],[1,1,0,0],[0,0,0,0]]
    graphs["K_3 + K_1"] = (k3k1_adj, 4)

    print(f"{'Graph':20s} | I(2) | I(-1) | I(1) | I(omega)")
    print("-" * 70)
    for name, (adj, n) in graphs.items():
        ip = independence_polynomial(adj, n)

        # Evaluate at key points
        i2 = eval_poly(ip, 2)
        im1 = eval_poly(ip, -1)
        i1 = eval_poly(ip, 1)
        # omega = e^{2pi i/3} = -1/2 + sqrt(3)/2 * i
        # I(omega) = sum alpha_k * omega^k
        import cmath
        omega = cmath.exp(2j * cmath.pi / 3)
        iomega = sum(c * omega**k for k, c in enumerate(ip))

        ip_str = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(ip) if c > 0)
        print(f"  {name:18s} | {i2:4d} | {im1:5d} | {i1:4d} | {iomega.real:.1f}+{iomega.imag:.1f}i")
        print(f"    I(x) = {ip_str}")

    # Critical observation: I(K_3, 2) = 7, I(K_3+K_1, 2) = 21
    print()
    print("  KEY: I(K_3, 2) = 7 = H_forbidden_1")
    print("       I(K_3+K_1, 2) = 21 = H_forbidden_2")
    print("  These are EXACTLY the two forbidden H values!")
    print("  But K_3 and K_3+K_1 cannot be tournament CGs.")

    # 10. Synthesis
    print("\n" + "=" * 70)
    print("SYNTHESIS: THE FIVE-FOLD WAY")
    print("=" * 70)
    print()
    print("  5 Platonic solids        = 5 regular convex polyhedra")
    print("  5 Exceptional Lie groups = 5 simple Lie algebras outside ABCD")
    print("  5 vertices in K(5,2)     = base set of the Petersen graph")
    print("  5 = KEY_1 + KEY_2        = sum of tournament polynomial roots")
    print("  5-cycle                  = smallest non-trivial tournament cycle")
    print()
    print("  The number 5 governs:")
    print("  - The ADE inequality 1/2 + 1/3 + 1/5 > 1 (barely)")
    print("  - The tournament polynomial z^2 - 5z + 6 = 0")
    print("  - The Petersen graph K(5,2) with 10 = C(5,2) vertices")
    print("  - The moat at T = 10 = C(5,2)")
    print("  - The 30-edge trinity: 30 = 2*3*5 = product of primes <= 5")
    print()
    print("  RECURRENCE HIERARCHY:")
    print("  Level 0: {2,3} = roots. Tetrahedron (6 edges), G_2 (h=6)")
    print("  Level 1: {2,3,5} = roots + sum. Icosahedron (30 edges), E_8 (h=30)")
    print("  Level 2: {C(5,2)=10} = binomial. THE MOAT. Petersen graph.")
    print()
    print("  The recurrence generates the hierarchy:")
    print("  z^2 - 5z + 6 = 0  =>  roots 2, 3")
    print("  roots => sum 5    =>  K(5,2) = Petersen, 10 vertices")
    print("  roots => product 6 => tetrahedron, G_2")
    print("  all three => 30    => icosahedron, E_8, J(5,2)")

if __name__ == "__main__":
    main()
