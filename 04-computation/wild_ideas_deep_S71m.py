#!/usr/bin/env python3
"""
WILD IDEAS — DEEP DIVES
opus-2026-03-14-S71m

Following up on the most striking findings from the initial wild ideas run:
1. corr(Z, H) = 1.0 — H is an eigenfunction of the hypercube adjacency?
2. Entropy ~0.27 bits/arc — is this a universal constant?
3. Normalized multiplicities sum to 2^(m-3) — combinatorial structure?
4. H mod 7 = 0 forbidden at n=5 — extends to larger n?
5. Sublevel sets always connected — topological proof?
6. NEW: H as spherical function on the symmetric group
7. NEW: Galois group of H's minimal polynomial
8. NEW: H and the Tutte polynomial connection
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
from fractions import Fraction

def adj_matrix(n, bits):
    """Tournament adjacency matrix from bit encoding."""
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
    """Count Hamiltonian paths via DP bitmask."""
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

def compute_F_polynomial(n, A):
    """F-polynomial: F_k = number of HPs with k ascents."""
    F = [0] * n
    for perm in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if not A[perm[k]][perm[k+1]]:
                valid = False
                break
        if valid:
            asc = sum(1 for k in range(n-1) if perm[k] < perm[k+1])
            F[asc] += 1
    return F

def get_all_tournaments(n):
    """Compute all tournaments and their H values."""
    m = n * (n-1) // 2
    results = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results.append((bits, A, H))
    return results

print("=" * 70)
print("DEEP DIVE: WILD IDEAS FOLLOW-UP")
print("opus-2026-03-14-S71m")
print("=" * 70)

# ======================================================================
# DEEP DIVE 1: IS H AN EIGENFUNCTION OF THE HYPERCUBE ADJACENCY?
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 1: H AS EIGENFUNCTION OF HYPERCUBE ADJACENCY")
print("=" * 70)

for n in [3, 4, 5]:
    m = n * (n-1) // 2
    tournaments = get_all_tournaments(n)
    H_vals = {bits: H for bits, A, H in tournaments}

    # For each tournament, compute Z = sum of H values of all arc-flip neighbors
    Z_vals = {}
    for bits, A, H in tournaments:
        Z = 0
        for arc in range(m):
            neighbor = bits ^ (1 << arc)
            Z += H_vals[neighbor]
        Z_vals[bits] = Z

    # Check if Z = alpha * H + beta (affine eigenfunction)
    # Using two points to find alpha, beta
    H_list = [H for _, _, H in tournaments]
    Z_list = [Z_vals[bits] for bits, _, _ in tournaments]

    # Linear regression: Z = alpha * H + beta
    n_pts = len(H_list)
    sum_H = sum(H_list)
    sum_Z = sum(Z_list)
    sum_HH = sum(h*h for h in H_list)
    sum_HZ = sum(h*z for h, z in zip(H_list, Z_list))

    denom = n_pts * sum_HH - sum_H**2
    if denom != 0:
        alpha = (n_pts * sum_HZ - sum_H * sum_Z) / denom
        beta = (sum_Z - alpha * sum_H) / n_pts
    else:
        alpha, beta = 0, sum_Z / n_pts

    # Check residuals
    max_resid = max(abs(Z_vals[bits] - (alpha * H + beta))
                    for bits, _, H in tournaments)

    print(f"\n  n={n} (m={m}):")
    print(f"    Z = sum of H(neighbors)")
    print(f"    Best fit: Z = {alpha:.6f} * H + {beta:.6f}")
    print(f"    Max residual: {max_resid:.10f}")

    if max_resid < 1e-6:
        print(f"    *** EXACT: Z = {alpha} * H + {beta} ***")
        # This means H is an eigenfunction of the adjacency operator!
        # With eigenvalue alpha (shifted by beta)
        # Since the adjacency of Q_m has eigenvalues m-2k for k=0,...,m
        # we need alpha = m - 2k for some k
        print(f"    Eigenvalue interpretation: lambda = {alpha}")
        print(f"    m = {m}, so m - lambda = {m - alpha}, (m-lambda)/2 = {(m-alpha)/2}")

        # Also check: is H a linear combination of Walsh characters?
        # Z = alpha*H + beta means AH = alpha*H + beta*1
        # where 1 is the all-ones vector and A is the adjacency matrix
        # This means H - c*1 is an eigenvector of A with eigenvalue alpha
        # where c = beta/(alpha - m) if alpha != m
        if abs(alpha - m) > 1e-6:
            c = beta / (m - alpha)  # Note: A*1 = m*1
            print(f"    H - {c:.6f}*1 is eigenvector of Q_{m} adjacency with eigenvalue {alpha}")
    else:
        print(f"    NOT exactly affine. Checking quadratic...")

        # Maybe Z = alpha * H^2 + beta * H + gamma?
        # No, let's check if it's piecewise affine by H value
        H_to_Z = defaultdict(list)
        for bits, _, H in tournaments:
            H_to_Z[H].append(Z_vals[bits])

        print(f"    Z values by H class:")
        for H_val in sorted(H_to_Z.keys()):
            Z_set = sorted(set(H_to_Z[H_val]))
            if len(Z_set) <= 5:
                print(f"      H={H_val:4d}: Z in {Z_set}")
            else:
                print(f"      H={H_val:4d}: Z in {Z_set[:3]}...{Z_set[-2:]}, {len(Z_set)} distinct")

# ======================================================================
# DEEP DIVE 2: ENTROPY PER ARC — UNIVERSAL CONSTANT?
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 2: ENTROPY PER ARC")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    if n <= 5:
        tournaments = get_all_tournaments(n)
        H_counter = Counter(H for _, _, H in tournaments)
    else:
        # n=6: sample or compute
        tournaments = get_all_tournaments(n)
        H_counter = Counter(H for _, _, H in tournaments)

    total = sum(H_counter.values())
    entropy = 0
    for count in H_counter.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    spectrum_size = len(H_counter)
    max_entropy = math.log2(spectrum_size)

    print(f"\n  n={n} (m={m}):")
    print(f"    |H-spectrum| = {spectrum_size}")
    print(f"    Shannon entropy = {entropy:.6f} bits")
    print(f"    Entropy per arc = {entropy/m:.6f} bits/arc")
    print(f"    Efficiency = {entropy/max_entropy:.6f}")

    # Also compute: mutual information I(arc_1 ; H)
    # How much does knowing ONE arc tell you about H?
    # I(e; H) = H(H) - H(H|e=0) * P(e=0) - H(H|e=1) * P(e=1)
    # By symmetry P(e=0) = P(e=1) = 1/2 and H(H|e=0) = H(H|e=1) by relabeling
    # Actually not quite — depends on which arc
    if n <= 5:
        # Compute for arc 0 (between vertex 0 and vertex 1)
        H_given_0 = Counter()  # H values when arc 0 = 0
        H_given_1 = Counter()  # H values when arc 0 = 1
        for bits, _, H in tournaments:
            if bits & 1:
                H_given_1[H] += 1
            else:
                H_given_0[H] += 1

        def ent(counter):
            total = sum(counter.values())
            e = 0
            for c in counter.values():
                p = c / total
                if p > 0:
                    e -= p * math.log2(p)
            return e

        H_H = entropy
        H_H_given_arc = 0.5 * ent(H_given_0) + 0.5 * ent(H_given_1)
        MI = H_H - H_H_given_arc

        print(f"    I(arc_0 ; H) = {MI:.6f} bits")
        print(f"    I/H(H) = {MI/entropy:.6f}")

    if n == 6:
        break  # n=6 takes a while, don't go further

# ======================================================================
# DEEP DIVE 3: NORMALIZED MULTIPLICITIES — COMBINATORIAL STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 3: MULTIPLICITY STRUCTURE")
print("=" * 70)

for n in [3, 4, 5]:
    m = n * (n-1) // 2
    tournaments = get_all_tournaments(n)
    H_counter = Counter(H for _, _, H in tournaments)

    counts = [H_counter[h] for h in sorted(H_counter.keys())]
    g = counts[0]
    for c in counts[1:]:
        g = math.gcd(g, c)

    normalized = [c // g for c in counts]

    print(f"\n  n={n} (m={m}):")
    print(f"    H-values: {sorted(H_counter.keys())}")
    print(f"    Multiplicities: {counts}")
    print(f"    GCD = {g}")
    print(f"    Normalized: {normalized}")
    print(f"    Sum of normalized = {sum(normalized)}")
    print(f"    Is power of 2? {sum(normalized)} = 2^{math.log2(sum(normalized)):.2f}")
    print(f"    GCD = {g} = n! * ? = {g}/{math.factorial(n)} = {Fraction(g, math.factorial(n))}")

    # The GCD should be related to the automorphism group
    # For n=5: GCD=8, n!=120, 120/8=15, and 2^10/128=8
    # Hmm, 2^m / sum(normalized) = 2^m / 2^?
    print(f"    2^m / sum(norm) = {2**m} / {sum(normalized)} = {2**m / sum(normalized)}")
    print(f"    2^m / GCD = {2**m // g}")

# ======================================================================
# DEEP DIVE 4: H mod 7 = 0 FORBIDDEN — WHAT ABOUT n=6, 7?
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 4: FORBIDDEN RESIDUES")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    tournaments = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in tournaments))
    H_counter = Counter(H for _, _, H in tournaments)

    print(f"\n  n={n}: H-spectrum = {H_set}")

    for p in [3, 5, 7, 11, 13]:
        residues = sorted(set(h % p for h in H_set))
        missing = sorted(set(range(p)) - set(residues))
        if missing:
            print(f"    H mod {p}: missing residues = {missing}")
        else:
            # Check if any residue is "rare"
            res_counts = Counter(h % p for h in H_set)
            min_res = min(res_counts.values())
            if min_res <= 1:
                rare = [r for r, c in res_counts.items() if c <= 1]
                print(f"    H mod {p}: rare residues = {rare} (count={min_res})")

    # Specifically check: is H ≡ 0 mod 7 ever achievable?
    h_mod7_0 = [h for h in H_set if h % 7 == 0]
    if h_mod7_0:
        print(f"    H ≡ 0 mod 7: achievable! Values = {h_mod7_0}")
    else:
        print(f"    H ≡ 0 mod 7: FORBIDDEN (no H divisible by 7)")

# ======================================================================
# DEEP DIVE 5: SUBLEVEL SET CONNECTIVITY — ALGEBRAIC EXPLANATION
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 5: SUBLEVEL SET CONNECTIVITY")
print("=" * 70)

print("""
  The sublevel sets X_h = {T : H(T) <= h} are always connected
  in the hypercube graph Q_m.

  WHY? Consider two tournaments T, T' with H(T), H(T') <= h.
  Can we find a path T = T_0, T_1, ..., T_k = T' in Q_m where
  each T_i has H(T_i) <= h?

  CLAIM: The path through the TRANSITIVE tournament works.
  Every tournament can be "sorted" one arc flip at a time,
  reaching the transitive tournament while (mostly) increasing H.

  But wait — if H INCREASES along the path, the sublevel set condition
  might be violated (we need H <= h, but if we increase H we go above h).

  ALTERNATIVE: Maybe use the LEXICOGRAPHIC path.
  Start at T, flip arcs one by one to match T'.
  Each flip either increases or decreases H.

  The question is: can we always order the arc flips so that
  H never exceeds max(H(T), H(T'))?

  Let's check: for each pair of tournaments in the SAME H class,
  is there a path staying within that H class or below?
""")

# Check for n=4: within each H class, are all tournaments connected
# via paths that stay at that H level or below?
n = 4
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)
H_map = {bits: H for bits, _, H in tournaments}

for h_target in sorted(set(H for _, _, H in tournaments)):
    # BFS within the sublevel set X_{h_target}
    in_level = {bits for bits, _, H in tournaments if H == h_target}
    at_or_below = {bits for bits, _, H in tournaments if H <= h_target}

    if not in_level:
        continue

    # Check connectivity of at_or_below
    start = next(iter(at_or_below))
    visited = {start}
    queue = [start]
    while queue:
        curr = queue.pop(0)
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr in at_or_below and nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)

    print(f"  n=4, X_{{{h_target}}}: {len(at_or_below)} tournaments, "
          f"{len(visited)} reachable = {'CONNECTED' if len(visited) == len(at_or_below) else 'DISCONNECTED'}")

    # Also: within the EXACT level set {T : H(T) = h_target}
    if len(in_level) > 1:
        start = next(iter(in_level))
        visited = {start}
        queue = [start]
        while queue:
            curr = queue.pop(0)
            for arc in range(m):
                nbr = curr ^ (1 << arc)
                if nbr in in_level and nbr not in visited:
                    visited.add(nbr)
                    queue.append(nbr)
        print(f"          Level set {{H={h_target}}}: {len(in_level)} tournaments, "
              f"{len(visited)} reachable = {'CONNECTED' if len(visited) == len(in_level) else 'DISCONNECTED'}")

# ======================================================================
# DEEP DIVE 6: THE LAPLACIAN SPECTRUM AND H
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 6: H IN THE LAPLACIAN EIGENBASIS")
print("=" * 70)

print("""
  The hypercube Q_m has known eigenvalues: m - 2k for k = 0, ..., m.
  The eigenvectors are Walsh functions chi_S(T) = (-1)^{|T cap S|}.

  H = sum_S hat{H}[S] * chi_S where hat{H}[S] is the Walsh coefficient.

  If Z = alpha * H + beta, then H lives in a single eigenspace!
  This means hat{H}[S] = 0 except for |S| in one particular value.

  But we KNOW H has Walsh components at degrees 0, 2, 4, ..., n-1.
  So Z = alpha * H + beta CAN'T hold exactly if H has multiple degrees.

  The n=4 finding (corr = 1.0) must be because n=4 H is special.
  Let's verify what Walsh degrees n=4 H has.
""")

n = 4
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)
H_vals = {bits: H for bits, _, H in tournaments}

# Walsh transform of H
walsh_H = {}
for S in range(1 << m):
    coeff = 0
    for bits, _, H in tournaments:
        # chi_S(bits) = (-1)^{popcount(S & bits)}
        sign = (-1) ** bin(S & bits).count('1')
        coeff += sign * H
    walsh_H[S] = coeff / (1 << m)

# Group by Walsh degree (popcount of S)
degree_energy = defaultdict(float)
for S, c in walsh_H.items():
    deg = bin(S).count('1')
    degree_energy[deg] += c**2

print(f"  n=4 Walsh spectrum of H by degree:")
total_energy = sum(degree_energy.values())
for deg in sorted(degree_energy.keys()):
    if degree_energy[deg] > 1e-10:
        print(f"    Degree {deg}: energy = {degree_energy[deg]:.6f} "
              f"({100*degree_energy[deg]/total_energy:.2f}%)")

# For n=3 and n=5, check which degrees
for n_check in [3, 5]:
    m_check = n_check * (n_check-1) // 2
    if m_check > 10:
        print(f"\n  n={n_check}: skipping (2^{m_check} too large for Walsh transform)")
        continue

    t_check = get_all_tournaments(n_check)
    H_check = {bits: H for bits, _, H in t_check}

    degree_energy_check = defaultdict(float)
    for S in range(1 << m_check):
        coeff = 0
        for bits, _, H in t_check:
            sign = (-1) ** bin(S & bits).count('1')
            coeff += sign * H
        coeff /= (1 << m_check)
        deg = bin(S).count('1')
        degree_energy_check[deg] += coeff**2

    print(f"\n  n={n_check} Walsh spectrum of H by degree:")
    total = sum(degree_energy_check.values())
    for deg in sorted(degree_energy_check.keys()):
        if degree_energy_check[deg] > 1e-10:
            print(f"    Degree {deg}: energy = {degree_energy_check[deg]:.6f} "
                  f"({100*degree_energy_check[deg]/total:.2f}%)")

# Now understand WHY Z = alpha*H + beta at n=4
# The adjacency operator A on Q_m maps chi_S to (m - 2|S|) * chi_S
# So A*H = sum_S hat{H}[S] * (m - 2|S|) * chi_S
# If H only has degrees 0 and d, then:
# A*H = m * hat{H}[0] * chi_0 + (m-2d) * sum_{|S|=d} hat{H}[S] * chi_S
#     = m * c + (m-2d) * (H - c)
#     = (m-2d) * H + 2d * c
# where c = hat{H}[0] = mean(H)

print("\n  Checking if n=4 H is purely degree 0 + degree d:")
print(f"  n=4 has degrees with nonzero energy:")
active_degrees = [d for d, e in degree_energy.items() if e > 1e-10]
print(f"    Active degrees: {active_degrees}")

if len(active_degrees) == 2 and 0 in active_degrees:
    d = [x for x in active_degrees if x != 0][0]
    mean_H = sum(H for _, _, H in tournaments) / len(tournaments)
    alpha_pred = m - 2*d
    beta_pred = 2*d * mean_H
    print(f"    d = {d}, m = {m}")
    print(f"    Predicted: Z = {alpha_pred} * H + {beta_pred}")
    print(f"    = ({m} - 2*{d}) * H + 2*{d}*{mean_H}")

# ======================================================================
# DEEP DIVE 7: THE CHARACTERISTIC POLYNOMIAL OF H AS OPERATOR
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 7: MINIMAL POLYNOMIAL OF THE H OPERATOR")
print("=" * 70)

print("""
  Treat H as a diagonal matrix on the tournament hypercube.
  Its minimal polynomial is prod_{h in spectrum} (x - h).

  For n=5: spectrum = {1, 3, 5, 9, 11, 13, 15}
  MinPoly = (x-1)(x-3)(x-5)(x-9)(x-11)(x-13)(x-15)

  Interesting: the H values at n=5 are {1,3,5,9,11,13,15}.
  These are all ODD (Redei). But more: they avoid 7!

  The minimal polynomial factors over Q. But what if we look at it mod p?
""")

for n in [3, 4, 5]:
    tournaments = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in tournaments))

    print(f"\n  n={n}: H-spectrum = {H_set}")

    # Factors of the minimal polynomial mod small primes
    for p in [2, 3, 5, 7]:
        residues = sorted(set(h % p for h in H_set))
        print(f"    mod {p}: roots at {residues}, "
              f"{'SPLITS COMPLETELY' if len(residues) == len(H_set) else 'HAS COLLISIONS'}")
        if len(residues) < len(H_set):
            # Show which H values collide
            collisions = defaultdict(list)
            for h in H_set:
                collisions[h % p].append(h)
            for r in sorted(collisions.keys()):
                if len(collisions[r]) > 1:
                    print(f"      Collision at {r}: {collisions[r]}")

# ======================================================================
# DEEP DIVE 8: H AND GRAPH COLORING
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 8: CHROMATIC CONNECTIONS")
print("=" * 70)

print("""
  The conflict graph Omega(T) encodes the 3-cycle structure.
  I(Omega, 2) = H counts independent sets of size <= 2.

  But the CHROMATIC POLYNOMIAL chi(Omega, k) also involves
  independent sets (proper k-colorings = k independent sets covering V).

  For k=2: chi(Omega, 2) = 2 * I_exclusive(Omega, 2)
  where I_exclusive counts independent sets of size EXACTLY |V|/2
  using exactly 2 colors.

  Actually chi(Omega, k) = sum over proper k-colorings.
  This is the CHROMATIC polynomial, not the independence polynomial.

  The INDEPENDENCE polynomial is I(Omega, x) = sum_k i_k * x^k
  where i_k = number of independent sets of size k.
  I(Omega, 1) = total number of independent sets.

  H = I(Omega, 2) = sum_k i_k * 2^k = independence polynomial at x=2.

  What about I(Omega, 1) = total independent sets?
  And I(Omega, -1) = signed count?
""")

for n in [3, 4, 5]:
    tournaments = get_all_tournaments(n)

    # For each tournament, compute I(Omega, x) for x = -1, 0, 1, 2, 3
    I_values = defaultdict(list)

    for bits, A, H in tournaments:
        # Build Omega(T): vertices = edges of complete graph,
        # edge iff they form a 3-cycle
        edges = []
        for i in range(n):
            for j in range(i+1, n):
                edges.append((i, j))

        m_edges = len(edges)
        # Two edges (i,j) and (k,l) form a 3-cycle if they share a vertex
        # and the third pair completes a directed 3-cycle
        omega_adj = [[False]*m_edges for _ in range(m_edges)]

        for e1_idx in range(m_edges):
            for e2_idx in range(e1_idx+1, m_edges):
                i, j = edges[e1_idx]
                k, l = edges[e2_idx]
                # Check if these edges are in a common 3-cycle
                shared = set([i,j]) & set([k,l])
                if len(shared) == 1:
                    v = shared.pop()
                    others = [x for x in [i,j] if x != v] + [x for x in [k,l] if x != v]
                    a, b = others
                    # Check if {v, a, b} form a 3-cycle
                    cyc = (A[v][a] + A[a][b] + A[b][v])
                    anti = (A[a][v] + A[b][a] + A[v][b])
                    if cyc == 3 or anti == 3:
                        omega_adj[e1_idx][e2_idx] = True
                        omega_adj[e2_idx][e1_idx] = True

        # Count independent sets by size
        indep_counts = [0] * (m_edges + 1)
        for mask in range(1 << m_edges):
            vertices = [i for i in range(m_edges) if mask & (1 << i)]
            is_indep = True
            for i in range(len(vertices)):
                for j in range(i+1, len(vertices)):
                    if omega_adj[vertices[i]][vertices[j]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                indep_counts[len(vertices)] += 1

        # Compute I(Omega, x) for various x
        for x in [-1, 1, 2, 3]:
            I_x = sum(indep_counts[k] * x**k for k in range(m_edges + 1))
            I_values[(n, x)].append(I_x)

        # Verify H = I(Omega, 2)
        I_2 = sum(indep_counts[k] * 2**k for k in range(m_edges + 1))
        if I_2 != H:
            print(f"  BUG: I(Omega,2) = {I_2} != H = {H}")

    print(f"\n  n={n}:")
    for x in [-1, 1, 2, 3]:
        vals = sorted(set(I_values[(n, x)]))
        print(f"    I(Omega, {x:2d}): spectrum = {vals}")

# ======================================================================
# DEEP DIVE 9: PRODUCT STRUCTURE OF H
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 9: MULTIPLICATIVE vs ADDITIVE STRUCTURE")
print("=" * 70)

print("""
  H is always odd. Write H = 2k+1 for some k >= 0.
  k = (H-1)/2 = number of "excess" HPs above the minimum.

  For n=5: k values = {0, 1, 2, 4, 5, 6, 7}
  Missing: k=3 (i.e., H=7).

  Is k always avoiding 3? Or is 3 special?
  At n=3: k in {0, 1}, missing everything > 1
  At n=4: k in {0, 1, 2}, missing nothing? Actually H in {1,3,5}, k in {0,1,2}
  At n=5: k in {0,1,2,4,5,6,7}, missing k=3

  Deeper: is the SET of achievable k values a NUMERICAL SEMIGROUP?
  (closed under addition for large enough values?)
""")

for n in [3, 4, 5, 6]:
    tournaments = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in tournaments))
    k_set = sorted(set((h-1)//2 for h in H_set))

    # Check if k_set forms a numerical semigroup
    # A numerical semigroup contains 0, is closed under addition for large values
    # Its complement in N is finite (the "gaps")
    max_k = max(k_set)
    gaps = sorted(set(range(max_k + 1)) - set(k_set))

    print(f"\n  n={n}: k = (H-1)/2 values = {k_set}")
    print(f"    Gaps in [0, {max_k}]: {gaps}")

    # Check closure under addition (for small values)
    if len(k_set) > 1:
        sums = set()
        for a in k_set:
            for b in k_set:
                if a + b <= max_k:
                    sums.add(a + b)
        missing_sums = sums - set(k_set)
        if missing_sums:
            print(f"    NOT closed under addition: {sorted(missing_sums)} achievable as sums but not as k values")
        else:
            print(f"    Closed under pairwise sums up to {max_k}")

# ======================================================================
# DEEP DIVE 10: TOURNAMENT ENTROPY — H AS PARTITION FUNCTION
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 10: PARTITION FUNCTION INTERPRETATION")
print("=" * 70)

print("""
  In statistical mechanics, Z = sum_states exp(-beta * E(state)).
  H = sum_P 1 counts Hamiltonian paths with uniform weight.

  Define Z(T, beta) = sum_P exp(-beta * cost(P))
  where cost(P) = number of "backward" arcs (descents) in P.

  Z(T, 0) = H(T) (uniform weight).
  Z(T, infinity) = F_0(T) (counts ONLY monotone HPs).

  The FREE ENERGY F(beta) = -log Z(T, beta) / beta
  and the ENTROPY S(beta) = beta^2 * dF/dbeta

  Question: at what "temperature" 1/beta does the tournament
  undergo a "phase transition" between ordered and disordered?
""")

import cmath

n = 5
tournaments = get_all_tournaments(n)

# For a few representative tournaments, compute Z(beta)
H_to_bits = defaultdict(list)
for bits, A, H in tournaments:
    H_to_bits[H].append(bits)

for H_target in [1, 5, 9, 15]:
    bits = H_to_bits[H_target][0]
    A = adj_matrix(n, bits)
    F = compute_F_polynomial(n, A)

    print(f"\n  H={H_target}: F-polynomial = {F}")
    print(f"    F(1) = {sum(F)} = H")

    # Z(beta) = sum_k F_k * exp(-beta * k)
    # At beta=0: Z=H
    # At beta -> inf: Z -> F_{n-1} (if nonzero) or F_0
    # Actually cost = descents = (n-1) - ascents
    # Z(beta) = sum_k F_k * exp(-beta * (n-1-k)) = exp(-beta*(n-1)) * F(exp(beta))

    # Phase transition: look for zeros of Z(beta) near real axis
    # Lee-Yang theorem: zeros of partition function determine phase transitions
    # Z(q) = F(q) as a polynomial in q = exp(beta)

    # Zeros of F-polynomial
    if sum(1 for f in F if f > 0) > 1:
        # Find roots numerically
        # F(q) = sum_k F_k * q^k
        coeffs = F  # F_0 + F_1*q + ... + F_{n-1}*q^{n-1}

        # Use numpy-free root finding? Let's just evaluate
        print(f"    F(omega) = F(e^{{2pi i/3}}) = ", end="")
        omega = cmath.exp(2j * cmath.pi / 3)
        F_omega = sum(F[k] * omega**k for k in range(n))
        print(f"{F_omega:.4f}, |F(omega)|^2 = {abs(F_omega)**2:.4f}")

        # F(-1) = alternating sum
        F_neg1 = sum(F[k] * (-1)**k for k in range(n))
        print(f"    F(-1) = {F_neg1}")

        # F(i) = evaluation at i
        F_i = sum(F[k] * (1j)**k for k in range(n))
        print(f"    F(i) = {F_i:.4f}, |F(i)|^2 = {abs(F_i)**2:.4f}")

# ======================================================================
# DEEP DIVE 11: SYMMETRY BREAKING — H UNDER VERTEX RELABELING
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 11: SYMMETRY BREAKING PATTERN")
print("=" * 70)

n = 5
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)

# Group tournaments by isomorphism class (by H value AND score sequence)
from collections import Counter

def score_seq(n, A):
    """Sorted out-degree sequence."""
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

iso_classes = defaultdict(list)
for bits, A, H in tournaments:
    ss = score_seq(n, A)
    iso_classes[(H, ss)].append(bits)

print(f"  n={n}: {len(iso_classes)} (H, score-seq) classes")
print(f"  Score sequences by H:")
for (H_val, ss), members in sorted(iso_classes.items()):
    print(f"    H={H_val:4d}, score={ss}: {len(members)} tournaments")

# ======================================================================
# DEEP DIVE 12: NEW — THE MÖBIUS FUNCTION ON H-SPECTRUM
# ======================================================================
print("\n" + "=" * 70)
print("DEEP DIVE 12: MÖBIUS FUNCTION ON H-SPECTRUM ORDER")
print("=" * 70)

print("""
  Define a PARTIAL ORDER on the H-spectrum:
  h1 <= h2 if there exists a sequence of arc flips
  T_1 -> T_2 -> ... -> T_k where H(T_1) = h1, H(T_k) = h2,
  and H is non-decreasing along the path.

  This defines the "H-reachability order."

  For n=4: Can we reach H=5 from H=1 via non-decreasing H path?
  For n=5: What's the Hasse diagram of this order?
""")

n = 4
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)
H_map = {bits: H for bits, _, H in tournaments}
H_values = sorted(set(H for _, _, H in tournaments))

# For each pair (h1, h2) with h1 < h2, check if there's a non-decreasing path
reachable = defaultdict(set)
for h1 in H_values:
    # BFS from all tournaments with H=h1
    starts = {bits for bits, _, H in tournaments if H == h1}
    visited = set(starts)
    queue = list(starts)
    reached_H = {h1}

    while queue:
        curr = queue.pop(0)
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr not in visited and H_map[nbr] >= H_map[curr]:
                visited.add(nbr)
                queue.append(nbr)
                reached_H.add(H_map[nbr])

    reachable[h1] = reached_H

print(f"\n  n=4: H-reachability:")
for h1 in H_values:
    print(f"    From H={h1}: can reach {sorted(reachable[h1])}")

# Same for n=5
n = 5
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)
H_map = {bits: H for bits, _, H in tournaments}
H_values = sorted(set(H for _, _, H in tournaments))

reachable = defaultdict(set)
for h1 in H_values:
    starts = {bits for bits, _, H in tournaments if H == h1}
    visited = set(starts)
    queue = list(starts)
    reached_H = {h1}

    while queue:
        curr = queue.pop(0)
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr not in visited and H_map[nbr] >= H_map[curr]:
                visited.add(nbr)
                queue.append(nbr)
                reached_H.add(H_map[nbr])

    reachable[h1] = reached_H

print(f"\n  n=5: H-reachability:")
for h1 in H_values:
    print(f"    From H={h1}: can reach {sorted(reachable[h1])}")

# Check: is the reachability order TOTAL? (every pair comparable)
is_total = True
for h1 in H_values:
    for h2 in H_values:
        if h1 < h2 and h2 not in reachable[h1]:
            print(f"    INCOMPARABLE: cannot reach H={h2} from H={h1} via non-decreasing path!")
            is_total = False

if is_total:
    print(f"    Reachability order is TOTAL")

print("\n" + "=" * 70)
print("DONE — DEEP DIVES COMPLETE")
print("=" * 70)
