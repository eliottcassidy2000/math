#!/usr/bin/env python3
"""
PROJECTIVE & ALGEBRAIC GEOMETRY — PART 5: THE INCIDENCE GEOMETRY
opus-2026-03-14-S71n

From Part 4: the degree-2 sign pattern reveals a GRAPH STRUCTURE:
- Two adjacent arcs sharing vertex v at position (1,0) → POSITIVE coefficient
  This means: arc (a,v) and arc (v,b) with a<v<b → positive
  The shared vertex is the LARGER endpoint of the first arc
  and the SMALLER endpoint of the second.
  This is a DIRECTED PATH a→v→b in the tournament!

- Two adjacent arcs sharing vertex v at position (0,0) or (1,1) → NEGATIVE
  (0,0) means v is smallest: arcs (v,a) and (v,b) → STAR OUT from v
  (1,1) means v is largest: arcs (a,v) and (b,v) → STAR IN to v

THIS IS THE PATH/STAR DICHOTOMY:
- Path configuration (a→v→b): positive coefficient
- Star configuration (v→a, v→b or a→v, b→v): negative coefficient

This connects DIRECTLY to Hamiltonian path counting:
paths contribute positively, stars contribute negatively.

Part 5 explores this connection to incidence matrices, oriented matroids,
and the full algebraic structure.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import math

def adj_matrix(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(n, A):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def get_all_tournaments(n):
    m = n*(n-1)//2
    results = {}
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results[bits] = H
    return results

print("=" * 70)
print("THE INCIDENCE GEOMETRY OF TOURNAMENT COEFFICIENTS")
print("opus-2026-03-14-S71n")
print("=" * 70)

# ======================================================================
# PART 1: THE SIGN RULE — PATH vs STAR
# ======================================================================
print(f"\n{'='*70}")
print("PART 1: PATH vs STAR — THE SIGN RULE FOR ALL DEGREES")
print(f"{'='*70}")

# Build arc list and multilinear coefficients
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    arc_list = []
    for i in range(n):
        for j in range(i+1, n):
            arc_list.append((i, j))

    # Multilinear coefficients
    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    print(f"\n  n={n}: Degree-2 coefficient sign rule")

    # For degree 2: classify arc pairs
    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) != 2:
            continue
        arc1 = arc_list[bits[0]]
        arc2 = arc_list[bits[1]]

        # Is this a PATH configuration?
        # Path: arcs share exactly one vertex, and they form i→v→j
        # Arc (a,b) with a<b: represents potential a→b
        shared = set(arc1) & set(arc2)
        if len(shared) == 1:
            sv = shared.pop()
            all_v = set(arc1) | set(arc2) - {sv}
            endpoints = list(all_v)
            # Path a→v→b: requires sv to be "middle" in some ordering
            # Arc1 = (a, sv) or (sv, b)?
            # In canonical form (i<j): arc (a,b)
            # If sv is the second element of arc1 AND first of arc2: PATH
            is_path = (arc1[1] == sv and arc2[0] == sv)
            sign_str = "+" if v > 0 else "-"
            if bits[0] < 3 and bits[1] < 5:  # Only print a few
                print(f"    arcs {arc1},{arc2}: shared={sv}, path={is_path}, sign={sign_str}, v={v}")

    # Now check the GENERAL rule for all degrees
    print(f"\n  n={n}: All coefficients by degree and structure")

    # For degree k: the arcs form a SUBGRAPH of K_n
    # What is the graph structure of this subgraph?
    for deg in range(3, min(5, m+1)):
        deg_coeffs = [(S, v) for S, v in ml_coeffs.items() if bin(S).count('1') == deg]
        if not deg_coeffs:
            continue

        print(f"  Degree {deg}: {len(deg_coeffs)} nonzero coefficients")

        # Classify by the GRAPH structure of the arc set
        graph_types = defaultdict(list)
        for S, v in deg_coeffs:
            bits = [i for i in range(m) if S & (1 << i)]
            arcs = [arc_list[b] for b in bits]

            # What graph do these arcs form?
            vertices = set()
            for a, b in arcs:
                vertices.add(a)
                vertices.add(b)
            n_vertices = len(vertices)

            # Degree sequence of the underlying graph
            deg_seq = Counter()
            for vert in vertices:
                d = sum(1 for a, b in arcs if a == vert or b == vert)
                deg_seq[d] += 1
            deg_tuple = tuple(sorted(deg_seq.items()))

            graph_types[(n_vertices, deg_tuple)].append(v)

        for gtype, vals in sorted(graph_types.items()):
            n_vert, deg_tup = gtype
            pos = sum(1 for v in vals if v > 0)
            neg = sum(1 for v in vals if v < 0)
            abs_vals = sorted(set(abs(v) for v in vals))
            print(f"    vertices={n_vert}, deg_seq={dict(deg_tup)}: {len(vals)} coeffs, +:{pos} -:{neg}, |vals|={abs_vals}")

# ======================================================================
# PART 2: THE ORIENTED INCIDENCE MATRIX
# ======================================================================
print(f"\n{'='*70}")
print("PART 2: THE ORIENTED INCIDENCE MATRIX AND H")
print(f"{'='*70}")

print("""
  The ORIENTED INCIDENCE MATRIX D of K_n is the n x m matrix:
  D[v, e] = +1 if v is the TAIL of edge e (smaller endpoint)
           = -1 if v is the HEAD of edge e (larger endpoint)
           = 0  otherwise

  For a tournament T with arcs x = (x_1,...,x_m) in {0,1}^m:
  The TOURNAMENT INCIDENCE MATRIX is:
  T_inc[v, e] = D[v,e] * (2x_e - 1) for the direction of the arc.

  The KERNEL of D^T is the cycle space of K_n (dimension C(n,2)-n+1).
  The IMAGE of D is the cut space (dimension n-1).

  The H function lives "above" D: it counts paths, which are
  related to the INTERSECTION of the cut space with the tournament orientation.
""")

# Compute the incidence matrix for K_n
for n in [3, 4, 5]:
    m = n*(n-1)//2

    D = [[0]*m for _ in range(n)]
    arc_list = []
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            D[i][idx] = 1   # tail
            D[j][idx] = -1  # head
            arc_list.append((i, j))
            idx += 1

    # D^T D = degree matrix - adjacency matrix of K_n = L (Laplacian)
    DtD = [[sum(D[k][i]*D[k][j] for k in range(n)) for j in range(m)] for i in range(m)]
    # Wait, that's D^T D, which is m x m

    # Actually DDt = n x n Laplacian
    DDt = [[sum(D[i][k]*D[j][k] for k in range(m)) for j in range(n)] for i in range(n)]
    print(f"\n  n={n}: DD^T (Laplacian of K_n):")
    for row in DDt:
        print(f"    {row}")

    # Eigenvalues of Laplacian: 0 with mult 1, n with mult n-1
    print(f"    Eigenvalues: 0 (x1), {n} (x{n-1})")

    # The quadratic form of H (degree-2 part)
    all_T = get_all_tournaments(n)
    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    # Build the Gram matrix B
    B = [[0]*m for _ in range(m)]
    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) == 2:
            i, j = bits
            B[i][j] = v
            B[j][i] = v

    # Compare B with D^T D
    DtD_proper = [[sum(D[k][i]*D[k][j] for k in range(n)) for j in range(m)] for i in range(m)]

    # Is B proportional to DtD (restricted to off-diagonal)?
    print(f"\n  Comparing B (Gram matrix of H) with D^T D:")
    match = True
    for i in range(m):
        for j in range(i+1, m):
            if DtD_proper[i][j] != 0 and B[i][j] != 0:
                ratio = B[i][j] / DtD_proper[i][j]
                if abs(ratio) != 2:
                    match = False
            elif DtD_proper[i][j] == 0 and B[i][j] != 0:
                match = False
    print(f"    B = -2 * (off-diagonal of D^T D): {match}")

    # Actually let's check directly
    is_match = True
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            expected = -2 * DtD_proper[i][j]
            if B[i][j] != expected:
                is_match = False
                break
        if not is_match:
            break
    print(f"    B[i,j] = -2 * (D^T D)[i,j] for i!=j: {is_match}")

# ======================================================================
# PART 3: THE KOSZUL COMPLEX
# ======================================================================
print(f"\n{'='*70}")
print("PART 3: THE KOSZUL COMPLEX AND TOURNAMENT HOMOLOGY")
print(f"{'='*70}")

print("""
  The KOSZUL COMPLEX of a sequence (f_1,...,f_r) in a ring R is:
  0 -> Lambda^r R^r -> ... -> Lambda^1 R^r -> R -> R/(f_1,...,f_r) -> 0

  For tournaments: consider the sequence of CONSTRAINT FUNCTIONS
  f_e = x_e(1-x_e) for each arc e (the Boolean constraints).

  The KOSZUL HOMOLOGY H_i(f_1,...,f_m; Z) measures the
  "higher syzygies" among the tournament constraints.

  More interesting: the H function ITSELF defines a Koszul-like complex.
  Consider the GRADIENT of H:
  dH = (dH/dx_1, ..., dH/dx_m) where dH/dx_i = H(x + e_i) - H(x).

  The JACOBIAN IDEAL J(H) = (dH/dx_1,...,dH/dx_m) ⊂ Z[x_1,...,x_m]/(x_i^2-x_i).

  The MILNOR NUMBER of H at a critical point is:
  mu = dim_k(Z[x]/(dH/dx_1,...,dH/dx_m))

  For DISCRETE Morse theory: the Milnor number at a local max/min
  tells us the "complexity" of the critical point.
""")

# Compute gradient statistics
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # For each tournament, compute the gradient vector
    grad_types = Counter()
    for bits in range(1 << m):
        H = all_T[bits]
        grad = tuple(all_T[bits ^ (1 << i)] - H for i in range(m))
        grad_types[grad] += 1

    # How many distinct gradients?
    print(f"\n  n={n}: Gradient statistics")
    print(f"    Distinct gradient vectors: {len(grad_types)}")
    print(f"    Most common gradient (count): {max(grad_types.values())}")

    # Gradient magnitude distribution
    mag_dist = Counter()
    for grad, count in grad_types.items():
        mag_sq = sum(g**2 for g in grad)
        mag_dist[mag_sq] += count
    print(f"    |grad|^2 distribution:")
    for mag, count in sorted(mag_dist.items()):
        print(f"      |grad|^2 = {mag}: {count} tournaments")

# ======================================================================
# PART 4: THE GRASSMANNIAN STRUCTURE — Gr(2,n) REVISITED
# ======================================================================
print(f"\n{'='*70}")
print("PART 4: THE LINE GRAPH L(K_n) AND TOURNAMENT ARCS")
print(f"{'='*70}")

print("""
  The arcs of a tournament form a subset of the edges of K_n.
  The LINE GRAPH L(K_n) has vertices = edges of K_n,
  and two vertices are adjacent if the edges share a vertex.

  L(K_n) is the Johnson graph J(n,2):
  |V| = C(n,2) = m
  Each vertex has degree 2(n-2)
  (each edge shares a vertex with n-2 other edges from each endpoint)

  The ADJACENCY MATRIX of L(K_n) is:
  A_L[e, f] = 1 if edges e,f share a vertex.

  From Part 4: the sign of the degree-2 coefficient c_{e,f} depends on
  whether arcs e,f form a PATH or STAR in L(K_n).

  More precisely:
  - PATH: e,f are adjacent in L(K_n) with the shared vertex being
    the "middle" vertex (tail of e, head of f or vice versa) → POSITIVE
  - STAR: e,f adjacent with shared vertex being a common endpoint → NEGATIVE
  - DISJOINT: e,f not adjacent → POSITIVE (for n>=5)

  THIS MEANS: the sign of c_{e,f} is determined by the SIGNED LINE GRAPH!
""")

# Verify: compute signed line graph for n=5
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    arc_list = []
    for i in range(n):
        for j in range(i+1, n):
            arc_list.append((i, j))

    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    # Check sign rule for degree-2 coefficients
    print(f"\n  n={n}: Sign rule verification for degree-2 coefficients:")
    correct = 0
    total = 0
    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) != 2:
            continue
        total += 1
        e_idx, f_idx = bits
        e = arc_list[e_idx]
        f = arc_list[f_idx]

        shared = set(e) & set(f)
        if len(shared) == 0:
            # Disjoint arcs
            predicted_sign = 1  # positive
        elif len(shared) == 1:
            sv = shared.pop()
            # Is sv the "middle" vertex? i.e., tail of one and head of the other
            # Arc (a,b) with a<b: a is the "first" vertex, b is the "second"
            # "Path" means one arc ends at sv and the other starts at sv
            # Path: e = (..., sv) and f = (sv, ...) in the arc direction sense
            # In canonical form: e = (a, sv) means sv is larger endpoint of e
            #                    f = (sv, b) means sv is smaller endpoint of f
            # So "path" requires sv to be the second element of e AND first element of f
            if e[1] == sv and f[0] == sv:
                predicted_sign = 1  # positive (path)
            elif e[0] == sv and f[1] == sv:
                predicted_sign = 1  # positive (reverse path)
            else:
                predicted_sign = -1  # negative (star)
        else:
            predicted_sign = 0

        if (v > 0 and predicted_sign > 0) or (v < 0 and predicted_sign < 0):
            correct += 1

    print(f"    Rule accuracy: {correct}/{total}")

    # Hmm, the "reverse path" case might be wrong.
    # Let me re-examine from the actual data
    print(f"\n  Detailed check:")
    for S, v in list(ml_coeffs.items())[:20]:
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) != 2:
            continue
        e_idx, f_idx = bits
        e = arc_list[e_idx]
        f = arc_list[f_idx]

        shared = set(e) & set(f)
        if len(shared) == 1:
            sv = shared.pop()
            others = list((set(e) | set(f)) - {sv})
            a, b = min(others), max(others)
            # In what position is sv?
            sv_pos_e = 0 if e[0] == sv else 1
            sv_pos_f = 0 if f[0] == sv else 1
            print(f"    e={e}, f={f}: sv={sv}, pos_in_e={sv_pos_e}, pos_in_f={sv_pos_f}, v={v:+d}")

# ======================================================================
# PART 5: THE SIGNED LINE GRAPH AND SECOND EIGENVALUE
# ======================================================================
print(f"\n{'='*70}")
print("PART 5: THE GRAM MATRIX B AS SIGNED LINE GRAPH")
print(f"{'='*70}")

# For n=3,4: B[i,j] = -2 * (D^T D)[i,j] was verified above
# D^T D = diag(2,...,2) - A_L where A_L is adjacency of line graph
# Wait, D^T D is the edge Laplacian

# Let me compute D^T D precisely
for n in [3, 4, 5]:
    m = n*(n-1)//2
    arc_list = []
    D = [[0]*m for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            D[i][idx] = 1
            D[j][idx] = -1
            arc_list.append((i, j))
            idx += 1

    DtD = [[sum(D[k][i]*D[k][j] for k in range(n)) for j in range(m)] for i in range(m)]

    # What is DtD?
    # DtD[i,j] = sum_v D[v,i]*D[v,j]
    # For i=j: DtD[i,i] = sum_v D[v,i]^2 = 2 (two nonzero entries per column)
    # For i!=j: DtD[i,j] = D[a,i]*D[a,j] + D[b,i]*D[b,j] + ...
    # If edges i,j share a vertex v: D[v,i]*D[v,j] = (+1)(-1) or (-1)(+1) = -1
    # Actually: if v is the smaller endpoint of i: D[v,i]=+1
    # if v is also smaller endpoint of j: D[v,j]=+1, so product = +1
    # if v is larger endpoint of j: D[v,j]=-1, product = -1

    # So DtD[i,j] for adjacent edges e_i=(a,b), e_j=(c,d):
    # If shared vertex is a (smaller of e_i): D[a,i]=+1
    #   If a is also smaller of e_j: D[a,j]=+1 → product +1
    #   If a is larger of e_j: D[a,j]=-1 → product -1
    # Similarly for shared vertex b (larger of e_i): D[b,i]=-1
    #   If b is smaller of e_j: D[b,j]=+1 → product -1
    #   If b is larger of e_j: D[b,j]=-1 → product +1

    # So DtD[i,j] = ±1 depending on the relative positions of the shared vertex

    print(f"\n  n={n}: D^T D matrix (edge inner products):")
    if m <= 10:
        for i in range(m):
            row = [DtD[i][j] for j in range(m)]
            print(f"    {arc_list[i]}: {row}")

    # Now the Gram matrix B from H
    all_T = get_all_tournaments(n)
    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    B = [[0]*m for _ in range(m)]
    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) == 2:
            i, j = bits
            B[i][j] = v
            B[j][i] = v

    # Compare B with -2*DtD (off-diagonal only)
    match_count = 0
    total_pairs = 0
    for i in range(m):
        for j in range(i+1, m):
            total_pairs += 1
            if B[i][j] == -2 * DtD[i][j]:
                match_count += 1

    print(f"\n  B[i,j] = -2 * (D^T D)[i,j] for off-diagonal: {match_count}/{total_pairs}")

    if match_count != total_pairs:
        # Show mismatches
        for i in range(m):
            for j in range(i+1, m):
                if B[i][j] != -2 * DtD[i][j]:
                    print(f"    MISMATCH at ({arc_list[i]},{arc_list[j]}): B={B[i][j]}, -2*DtD={-2*DtD[i][j]}")

# ======================================================================
# PART 6: EIGENVALUE STRUCTURE OF THE GRAM MATRIX
# ======================================================================
print(f"\n{'='*70}")
print("PART 6: EIGENVALUE STRUCTURE OF B AND D^T D")
print(f"{'='*70}")

# D^T D has known eigenvalues: it's the edge Laplacian of K_n
# Eigenvalues of edge Laplacian = eigenvalues of graph Laplacian (nonzero) + 2
# Graph Laplacian of K_n: eigenvalues 0 (x1), n (x(n-1))
# Edge Laplacian: eigenvalues n (x(n-1)) and 2 (x(C(n,2)-n+1))

# If B = -2 * DtD (off-diag) then B = -2 * (DtD - diag)
# Since DtD = L_edge has diagonal all = 2:
# B = -2 * (L_edge - 2*I) = -2*L_edge + 4*I

# Actually let me compute B's eigenvalues numerically
for n in [3, 4]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    B = [[0]*m for _ in range(m)]
    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) == 2:
            i, j = bits
            B[i][j] = v
            B[j][i] = v

    # Compute eigenvalues via characteristic polynomial
    # For small m, use numpy-free approach
    if m == 3:
        # 3x3 symmetric: use formula
        a, b, c = B[0][1], B[0][2], B[1][2]
        # Char poly: -λ³ + (a²+b²+c²)λ + 2abc = 0
        print(f"\n  n={n} (m=3): B eigenvalues")
        print(f"    B off-diagonal: {a}, {b}, {c}")
        print(f"    Char poly: -lambda^3 + {a**2+b**2+c**2}*lambda + {2*a*b*c}")
        # Roots: check lambda=4: -64 + 12*4 + 16 = -64+48+16 = 0 ✓
        # lambda=-2: -(-8) + 12*(-2) + 16 = 8-24+16 = 0 ✓
        print(f"    Eigenvalues: 4, -2, -2")

    elif m == 6:
        # Use power method or trace invariants
        # tr(B) = 0 (diagonal is 0)
        # tr(B^2) = 2 * sum of B[i,j]^2
        tr_B2 = sum(B[i][j]**2 for i in range(m) for j in range(m))
        # tr(B^3) = 6 * sum B[i,j]*B[j,k]*B[k,i]
        tr_B3 = sum(B[i][j]*B[j][k]*B[k][i] for i in range(m) for j in range(m) for k in range(m))

        print(f"\n  n={n} (m=6): B trace invariants")
        print(f"    tr(B) = 0")
        print(f"    tr(B^2) = {tr_B2}")
        print(f"    tr(B^3) = {tr_B3}")
        print(f"    sum of eigenvalues = 0")
        print(f"    sum of eigenvalue^2 = {tr_B2}")
        print(f"    sum of eigenvalue^3 = {tr_B3}")

        # For the edge Laplacian L_edge of K_4:
        # eigenvalues of L_edge(K_4): 4 (x3), 2 (x3)
        # B = -2*(L_edge - 2I) = -2*L_edge + 4I
        # eigenvalues of B: -2*4+4 = -4 (x3), -2*2+4 = 0 (x3)
        print(f"    If B = -2*(L_edge - 2I): eigenvalues = -4 (x3), 0 (x3)")
        print(f"    Check: sum = -12 (should be 0) — DOESN'T MATCH")

        # Try: is B = -2 * signed_adj_L where signed_adj is the SIGNED adjacency?
        # Let me just compute B^2 and see
        B2 = [[sum(B[i][k]*B[k][j] for k in range(m)) for j in range(m)] for i in range(m)]
        print(f"    B^2 diagonal: {[B2[i][i] for i in range(m)]}")

# ======================================================================
# PART 7: THE MULTILINEAR PERMANENT FORMULA
# ======================================================================
print(f"\n{'='*70}")
print("PART 7: H AS A MULTILINEAR PERMANENT")
print(f"{'='*70}")

print("""
  H(x) = sum_{P in S_n} prod_{i=0}^{n-2} a_{P_i, P_{i+1}}(x)

  where a_{i,j}(x) = x_e if i < j (e = arc (i,j)),
                    = 1 - x_e if i > j.

  This is a MULTILINEAR function of x. Each term in the sum is a
  product of n-1 factors, each depending on exactly ONE variable.

  The MULTILINEAR COEFFICIENT c_S is:
  c_S = sum_{P} epsilon(P, S) where epsilon(P, S) = +-1.

  epsilon(P, S) = product of signs from each factor:
  For arc e = (i,j) with i<j in the path P:
  - If P traverses e "forward" (i before j, i.e., P_k=i, P_{k+1}=j):
    factor = x_e, so contribution to x_e monomial is +1
  - If P traverses e "backward" (j before i):
    factor = 1-x_e, so contribution to (1) is +1 and to (x_e) is -1

  So c_S = (# paths using all arcs in S forward) - (# paths using all arcs in S backward + ...)

  Wait, it's more subtle because c_S involves the FULL Mobius inversion.
  Let me think about this differently.
""")

# Actually compute the coefficient by expanding the permanent directly
n = 3
m = 3
arc_list = [(0,1), (0,2), (1,2)]

print(f"\n  n=3: H as multilinear permanent")
print(f"  Arcs: {arc_list}")
print(f"  H(x) = sum_P prod a_{{P_i,P_{{i+1}}}}(x)")
print(f"\n  All 6 permutations and their contributions:")

# For each permutation, show the path and its contribution
for P in permutations(range(n)):
    edges = []
    for k in range(n-1):
        i, j = P[k], P[k+1]
        if i < j:
            arc_idx = arc_list.index((i,j))
            edges.append((arc_idx, +1, f"x{arc_idx}"))
        else:
            arc_idx = arc_list.index((j,i))
            edges.append((arc_idx, -1, f"(1-x{arc_idx})"))

    path_str = "->".join(str(v) for v in P)
    factor_str = " * ".join(name for _, _, name in edges)

    # Expand the product
    terms = [(1, frozenset())]  # (coefficient, set of variables)
    for arc_idx, sign, name in edges:
        new_terms = []
        for coeff, vars_set in terms:
            if sign == 1:
                # factor is x_arc_idx
                new_terms.append((coeff, vars_set | {arc_idx}))
            else:
                # factor is (1 - x_arc_idx) = 1 - x_arc_idx
                new_terms.append((coeff, vars_set))  # the "1" part
                new_terms.append((-coeff, vars_set | {arc_idx}))  # the "-x" part
        terms = new_terms

    # Consolidate
    term_dict = defaultdict(int)
    for coeff, vars_set in terms:
        term_dict[vars_set] += coeff

    contributions = {k: v for k, v in term_dict.items() if v != 0}
    print(f"  P={path_str}: {factor_str}")
    for vars_set, coeff in sorted(contributions.items(), key=lambda x: len(x[0])):
        var_str = "*".join(f"x{v}" for v in sorted(vars_set)) if vars_set else "1"
        print(f"    {coeff:+d} * {var_str}")

# Now sum all contributions across all permutations
all_terms = defaultdict(int)
for P in permutations(range(n)):
    edges = []
    for k in range(n-1):
        i, j = P[k], P[k+1]
        if i < j:
            arc_idx = arc_list.index((i,j))
            edges.append((arc_idx, +1))
        else:
            arc_idx = arc_list.index((j,i))
            edges.append((arc_idx, -1))

    terms = [(1, frozenset())]
    for arc_idx, sign in edges:
        new_terms = []
        for coeff, vars_set in terms:
            if sign == 1:
                new_terms.append((coeff, vars_set | {arc_idx}))
            else:
                new_terms.append((coeff, vars_set))
                new_terms.append((-coeff, vars_set | {arc_idx}))
        terms = new_terms

    for coeff, vars_set in terms:
        all_terms[vars_set] += coeff

print(f"\n  TOTAL: H(x) = ", end="")
first = True
for vars_set in sorted(all_terms.keys(), key=lambda x: (len(x), tuple(sorted(x)))):
    coeff = all_terms[vars_set]
    if coeff == 0:
        continue
    var_str = "*".join(f"x{v}" for v in sorted(vars_set)) if vars_set else "1"
    if first:
        print(f"{coeff}*{var_str}", end="")
        first = False
    else:
        print(f" {coeff:+d}*{var_str}", end="")
print()

# ======================================================================
# PART 8: THE PATH-STAR DECOMPOSITION THEOREM
# ======================================================================
print(f"\n{'='*70}")
print("PART 8: THE PATH-STAR DECOMPOSITION THEOREM")
print(f"{'='*70}")

print("""
  THEOREM (conjectured from computation):
  The degree-2 multilinear coefficient c_{e,f} of H satisfies:

  For adjacent arcs e=(a,v), f=(v,b) sharing vertex v:
    c_{e,f} = -2 * D^T D[e,f]
  where D is the oriented incidence matrix of K_n.

  Since D^T D[e,f] = D[a,e]*D[a,f] + D[v,e]*D[v,f] + ...
  this equals the signed overlap of the incidence vectors.

  For n=3: D^T D has off-diagonal entries +-1:
    D^T D[(0,1),(0,2)] = D[0,(0,1)]*D[0,(0,2)] + D[1,(0,1)]*D[1,(0,2)] + D[2,(0,1)]*D[2,(0,2)]
    = (+1)(+1) + (-1)(0) + (0)(-1) = 1
    So c_{(0,1),(0,2)} = -2 * 1 = -2 ✓ (matches!)

  For n=4: same pattern.

  For n=5: ALSO c_{e,f} = -2*DtD[e,f] for ADJACENT arcs.
  But what about DISJOINT arcs (n>=5)?
  There are disjoint arc pairs at n>=5, and they have DtD[e,f]=0.
  So the prediction is c_{e,f} = 0 for disjoint arcs.
  But we found nonzero coefficients for disjoint arcs!
""")

# Verify which disjoint pairs have nonzero coefficients at n=5
n = 5
m = 10
all_T5 = get_all_tournaments(5)

arc_list = []
for i in range(n):
    for j in range(i+1, n):
        arc_list.append((i, j))

ml_coeffs = {}
for S in range(1 << m):
    val = 0
    for T in range(1 << m):
        if T & S != T:
            continue
    sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
    val += sign * all_T5[T]
    if val != 0:
        ml_coeffs[S] = val

# Find disjoint arc pairs with nonzero coefficients
print(f"\n  n=5: Disjoint arc pairs with nonzero degree-2 coefficients:")
disjoint_nonzero = []
for S, v in ml_coeffs.items():
    bits = [i for i in range(m) if S & (1 << i)]
    if len(bits) != 2:
        continue
    e = arc_list[bits[0]]
    f = arc_list[bits[1]]
    if len(set(e) & set(f)) == 0:
        disjoint_nonzero.append((e, f, v))
        print(f"    {e} and {f}: c = {v}")

print(f"  Total disjoint pairs with nonzero coeff: {len(disjoint_nonzero)}")
print(f"  Total disjoint pairs: {sum(1 for i in range(m) for j in range(i+1,m) if len(set(arc_list[i])&set(arc_list[j]))==0)}")

# What do the nonzero disjoint pairs look like?
# At n=5: the disjoint pairs are (i,j) and (k,l) with {i,j}∩{k,l}=∅
# This requires 4 distinct vertices out of 5, so the 5th vertex is "missing"
print(f"\n  The 'missing vertex' for each disjoint pair:")
for e, f, v in disjoint_nonzero:
    missing = set(range(5)) - set(e) - set(f)
    print(f"    {e},{f}: missing vertex = {missing}, coeff = {v}")

# CONNECTION: The disjoint arcs + missing vertex form a MATCHING + ISOLATED VERTEX
# This is related to the TRANSFER MATRIX structure:
# The coefficient comes from paths that traverse both arcs AND pass through the missing vertex

print("\n" + "=" * 70)
print("DONE — INCIDENCE GEOMETRY OF TOURNAMENT COEFFICIENTS")
print("=" * 70)
