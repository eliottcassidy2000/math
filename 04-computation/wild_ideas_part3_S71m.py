#!/usr/bin/env python3
"""
WILD IDEAS — PART 3: DEEP INVESTIGATIONS
opus-2026-03-14-S71m

Following up on three major discoveries:
1. H ≡ 0 mod 7 is FORBIDDEN for n=3,4,5,6 — why?
2. I(single arc; H) = 0 — every arc is independently uninformative
3. H-reachability is a total order — the landscape has no dead ends

Plus new ideas:
4. H mod 3 and the tribonacci connection
5. Cluster structure in tournament hypercube
6. The "H gradient flow" on Q_m
7. Lee-Yang zeros of the F-polynomial
8. Double counting: H via matrix permanent
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
from fractions import Fraction
import cmath

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
    m = n * (n-1) // 2
    results = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results.append((bits, A, H))
    return results

print("=" * 70)
print("WILD IDEAS — PART 3: DEEP INVESTIGATIONS")
print("opus-2026-03-14-S71m")
print("=" * 70)

# ======================================================================
# INVESTIGATION 1: WHY H ≡ 0 mod 7 IS FORBIDDEN
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 1: WHY H ≡ 0 mod 7 IS FORBIDDEN")
print("=" * 70)

print("""
  H = I(Omega(T), 2) = independence polynomial of Omega(T) at x=2.
  H = 1 + 2*|V(Omega)| + 4*|ind pairs| + 8*|ind triples| + ...
  H = 1 + 2*C(n,2) + 4*i_2 + 8*i_3 + ...
  where i_k = number of independent sets of size k in Omega(T).

  Wait, Omega(T) has C(n,2) vertices (the edges of K_n).
  An independent set in Omega(T) = set of edges with no two in a 3-cycle.

  Actually: H = I(Omega, 2) where I is the independence polynomial
  evaluated at x=2. But OCF says H = I(Omega(T), 2) where Omega(T)
  is the ODD CYCLE COLLECTION. Let's use the simpler formula.

  H ≡ 1 mod 2 always (Rédei).
  H mod 3: H = 1 + 2*C(n,2) + 4*i_2 + 8*i_3 + ...
           ≡ 1 + 2*C(n,2) + i_2 + 2*i_3 + i_4 + 2*i_5 + ... mod 3

  H mod 7: H ≡ 1 + 2*C(n,2) + 4*i_2 + i_3 + 2*i_4 + 4*i_5 + i_6 + ... mod 7

  For H ≡ 0 mod 7, we need this sum to vanish mod 7.

  But WAIT: Φ_3(2) = 4 + 2 + 1 = 7. And the powers of 2 mod 7 cycle:
  2^0 = 1, 2^1 = 2, 2^2 = 4, 2^3 = 1, 2^4 = 2, 2^5 = 4, ...
  Period 3! Because ord_7(2) = 3.

  So H mod 7 = sum_k i_k * 2^k mod 7
             = sum_{k≡0 mod 3} i_k * 1 + sum_{k≡1 mod 3} i_k * 2 + sum_{k≡2 mod 3} i_k * 4 mod 7
             = A + 2B + 4C mod 7
  where A = sum_{k≡0(3)} i_k, B = sum_{k≡1(3)} i_k, C = sum_{k≡2(3)} i_k.

  For H ≡ 0 mod 7: A + 2B + 4C ≡ 0 mod 7.

  But I(Omega, omega) where omega = e^{2pi i/3} has the property:
  I(Omega, omega) = A*omega^0 + B*omega^1 + C*omega^2 (grouping by residue mod 3)

  Wait, this isn't quite right. Let me think more carefully.

  I(Omega, x) = sum_k i_k * x^k.
  At x = 2: I = H.
  At x = 2*omega where omega^3 = 1:
  I(Omega, 2*omega) = sum_k i_k * (2*omega)^k = sum_k i_k * 2^k * omega^k

  The three values I(Omega, 2), I(Omega, 2*omega), I(Omega, 2*omega^2)
  determine A + 2B + 4C, A + 2B*omega + 4C*omega^2, etc.

  H mod 7 = 0 iff I(Omega, 2) mod 7 = 0 iff 7 | H.

  Since Φ_3(2) = 7, this is related to the factorization of I(Omega, x) mod Φ_3(x).

  KEY INSIGHT: 7 | H iff H ≡ 0 mod Φ_3(2).
  This connects the forbidden value 7 to the cyclotomic polynomial!

  If Φ_3(x) | I(Omega, x) in Z[x], then 7 = Φ_3(2) | I(Omega, 2) = H.
  But Φ_3(x) NEVER divides I(Omega, x) for any tournament!

  QUESTION: Is it true that I(Omega, omega) ≠ 0 for every tournament?
  If so, then Φ_3(x) ∤ I(Omega, x), and the 7-divisibility question
  becomes more nuanced.
""")

# Verify: compute I(Omega, x) at roots of unity
# Actually we should just check H mod 7 and H mod 3 patterns

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    tournaments = get_all_tournaments(n)
    H_vals = [H for _, _, H in tournaments]

    # H mod 7 distribution
    mod7 = Counter(h % 7 for h in H_vals)
    print(f"\n  n={n}: H mod 7 distribution:")
    for r in range(7):
        count = mod7.get(r, 0)
        pct = 100 * count / len(H_vals)
        print(f"    H ≡ {r} mod 7: {count:6d} ({pct:5.1f}%)")

    # H mod 3 distribution
    mod3 = Counter(h % 3 for h in H_vals)
    print(f"  H mod 3 distribution:")
    for r in range(3):
        count = mod3.get(r, 0)
        pct = 100 * count / len(H_vals)
        print(f"    H ≡ {r} mod 3: {count:6d} ({pct:5.1f}%)")

    # H mod 21 = mod (3*7) distribution
    mod21 = Counter(h % 21 for h in H_vals)
    achievable_mod21 = sorted(mod21.keys())
    missing_mod21 = sorted(set(range(21)) - set(achievable_mod21))
    print(f"  H mod 21: {len(achievable_mod21)} achievable residues, "
          f"missing {len(missing_mod21)}: {missing_mod21[:10]}{'...' if len(missing_mod21)>10 else ''}")

    if n == 6:
        break

# ======================================================================
# INVESTIGATION 2: H AND ORD_p(2) — THE CYCLOTOMIC CONNECTION
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 2: H AND THE ORDER OF 2 mod p")
print("=" * 70)

print("""
  ord_7(2) = 3 (since 2^3 = 8 ≡ 1 mod 7).
  Φ_3(2) = 7.

  For which primes p does ord_p(2) = 3?
  Equivalently, which primes divide Φ_3(2) = 7?
  Answer: only p = 7.

  More generally, Φ_d(2) is divisible by primes p with ord_p(2) = d.
  By Zsygmondy's theorem, Φ_d(2) has a "primitive prime divisor" for d ≥ 7.

  Primes p | Φ_d(2):
  d=1: Φ_1(2) = 1 (no prime)
  d=2: Φ_2(2) = 3
  d=3: Φ_3(2) = 7
  d=4: Φ_4(2) = 5
  d=5: Φ_5(2) = 31
  d=6: Φ_6(2) = 3 (note: 3 | Φ_2 and Φ_6)
  d=7: Φ_7(2) = 127
  d=8: Φ_8(2) = 17

  So:
  H ≡ 0 mod 3: forbidden? Let's check.
  H ≡ 0 mod 7: forbidden at n ≤ 6.
  H ≡ 0 mod 5: forbidden?
  H ≡ 0 mod 31: forbidden?
""")

for n in [3, 4, 5, 6]:
    tournaments = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in tournaments))

    print(f"\n  n={n}: H-spectrum mod various Φ_d(2) primes:")
    for p, d in [(3, 2), (7, 3), (5, 4), (31, 5), (127, 7), (17, 8)]:
        residues = set(h % p for h in H_set)
        if 0 not in residues:
            print(f"    H ≡ 0 mod {p:3d} (Φ_{d}(2)): FORBIDDEN")
        else:
            h_div = [h for h in H_set if h % p == 0]
            print(f"    H ≡ 0 mod {p:3d} (Φ_{d}(2)): achievable ({len(h_div)} values: {h_div[:5]})")

    if n == 6:
        break

# ======================================================================
# INVESTIGATION 3: PAIRWISE MUTUAL INFORMATION I(arc_i; arc_j | H=h)
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 3: PAIRWISE ARC INFORMATION GIVEN H")
print("=" * 70)

print("""
  We found I(single arc; H) = 0 — a single arc is independent of H.
  This makes sense: the H value doesn't depend on the orientation of
  any particular arc, by the symmetry of the uniform distribution.

  But what about I(arc_i, arc_j; H) — the JOINT information of two arcs?
  Or equivalently, I(arc_i; H | arc_j) — does knowing one arc reveal
  information about H ONCE YOU KNOW another arc?

  For arcs that form a TRIANGLE: (i,j), (j,k), (i,k):
  Knowing two determines whether there's a 3-cycle.
  Since H depends on 3-cycles, knowing two arcs in a triangle
  should give information about H.
""")

n = 5
m = n * (n-1) // 2
tournaments = get_all_tournaments(n)

# Compute I(arc_pair; H) for all pairs of arcs
# An arc pair can be: sharing a vertex (triangle-adjacent) or disjoint

def entropy(counter):
    total = sum(counter.values())
    e = 0
    for c in counter.values():
        p = c / total
        if p > 0:
            e -= p * math.log2(p)
    return e

H_counter = Counter(H for _, _, H in tournaments)
H_ent = entropy(H_counter)

# Arc pairs that share a vertex vs disjoint
arc_list = []
for i in range(n):
    for j in range(i+1, n):
        arc_list.append((i, j))

# Map arc to bit position
arc_to_bit = {}
idx = 0
for i in range(n):
    for j in range(i+1, n):
        arc_to_bit[(i,j)] = idx
        idx += 1

# Compute MI for each pair of arcs
MI_pairs = {}
for a1 in range(m):
    for a2 in range(a1+1, m):
        # Joint distribution of (arc_a1, arc_a2, H)
        joint = Counter()
        for bits, _, H in tournaments:
            v1 = (bits >> a1) & 1
            v2 = (bits >> a2) & 1
            joint[(v1, v2, H)] += 1

        # Marginal of (arc_pair)
        pair_counter = Counter()
        for (v1, v2, H), c in joint.items():
            pair_counter[(v1, v2)] += c

        # H(H | pair) = sum over pair values p * H(H|pair=p)
        cond_ent = 0
        total = sum(pair_counter.values())
        for (v1, v2), count in pair_counter.items():
            p = count / total
            h_given_pair = Counter()
            for (vv1, vv2, H), c in joint.items():
                if vv1 == v1 and vv2 == v2:
                    h_given_pair[H] += c
            cond_ent += p * entropy(h_given_pair)

        MI = H_ent - cond_ent
        MI_pairs[(a1, a2)] = MI

# Check: do adjacent arcs (sharing vertex) have more MI than disjoint?
adjacent_MIs = []
disjoint_MIs = []

for (a1, a2), mi in MI_pairs.items():
    i1, j1 = arc_list[a1]
    i2, j2 = arc_list[a2]
    shared = set([i1,j1]) & set([i2,j2])
    if shared:
        adjacent_MIs.append(mi)
    else:
        disjoint_MIs.append(mi)

print(f"\n  n={n}: I(arc_pair; H)")
print(f"  Adjacent pairs (share vertex): n={len(adjacent_MIs)}")
print(f"    mean MI = {sum(adjacent_MIs)/len(adjacent_MIs):.6f} bits")
print(f"    max MI = {max(adjacent_MIs):.6f} bits")
print(f"    min MI = {min(adjacent_MIs):.6f} bits")
print(f"  Disjoint pairs: n={len(disjoint_MIs)}")
print(f"    mean MI = {sum(disjoint_MIs)/len(disjoint_MIs):.6f} bits")
print(f"    max MI = {max(disjoint_MIs):.6f} bits")
print(f"    min MI = {min(disjoint_MIs):.6f} bits")

# Are all adjacent MIs the same? (by symmetry they should be)
adj_set = sorted(set(round(mi, 10) for mi in adjacent_MIs))
dis_set = sorted(set(round(mi, 10) for mi in disjoint_MIs))
print(f"  Distinct adjacent MI values: {len(adj_set)}")
print(f"  Distinct disjoint MI values: {len(dis_set)}")
if len(adj_set) == 1:
    print(f"  ALL adjacent pairs have IDENTICAL MI = {adj_set[0]:.10f}")
if len(dis_set) == 1:
    print(f"  ALL disjoint pairs have IDENTICAL MI = {dis_set[0]:.10f}")

# ======================================================================
# INVESTIGATION 4: THE GRADIENT FLOW ON Q_m
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 4: H GRADIENT FLOW ON TOURNAMENT HYPERCUBE")
print("=" * 70)

print("""
  Define the "steepest ascent" flow: from tournament T, flip the arc
  that increases H the most (or stay if all flips decrease H).

  This creates a DIRECTED GRAPH on tournaments: the gradient flow.

  Questions:
  - How many SINKS (local maxima)?
  - How many BASINS OF ATTRACTION?
  - What is the typical path length from a random tournament to a sink?
  - Do all paths lead to the same maximum?
""")

for n in [3, 4, 5]:
    m_val = n * (n-1) // 2
    t = get_all_tournaments(n)
    H_map = {bits: H for bits, _, H in t}

    # Find gradient flow: for each tournament, the arc flip that
    # increases H the most
    gradient = {}  # bits -> next_bits (or None if local max)
    for bits, _, H in t:
        best_nbr = None
        best_delta = 0
        for arc in range(m_val):
            nbr = bits ^ (1 << arc)
            delta = H_map[nbr] - H
            if delta > best_delta:
                best_delta = delta
                best_nbr = nbr
        gradient[bits] = best_nbr

    # Find sinks (local maxima)
    sinks = {bits for bits, nbr in gradient.items() if nbr is None}

    # Find basins of attraction
    def follow_gradient(start):
        """Follow gradient flow to sink."""
        curr = start
        path_len = 0
        while gradient[curr] is not None:
            curr = gradient[curr]
            path_len += 1
            if path_len > 1000:
                return curr, path_len  # safety
        return curr, path_len

    basins = defaultdict(list)
    path_lengths = []
    for bits, _, H in t:
        sink, plen = follow_gradient(bits)
        basins[sink].append(bits)
        path_lengths.append(plen)

    print(f"\n  n={n} (m={m_val}, {len(t)} tournaments):")
    print(f"    Local maxima (sinks): {len(sinks)}")
    print(f"    H values at sinks: {sorted(set(H_map[s] for s in sinks))}")
    print(f"    Basins of attraction: {len(basins)}")
    for sink in sorted(basins.keys(), key=lambda s: -H_map[s]):
        print(f"      Sink H={H_map[sink]}: basin size = {len(basins[sink])}")
    print(f"    Path lengths: min={min(path_lengths)}, max={max(path_lengths)}, "
          f"mean={sum(path_lengths)/len(path_lengths):.2f}")

# ======================================================================
# INVESTIGATION 5: LEE-YANG ZEROS OF THE F-POLYNOMIAL
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 5: ZEROS OF THE F-POLYNOMIAL")
print("=" * 70)

print("""
  The F-polynomial F(T, x) = sum_k F_k x^k counts HPs by ascent count.
  Its roots (Lee-Yang zeros) in the complex plane reveal structure.

  For the transitive tournament: F = x^{n-1}, all roots at x=0.
  For the regular tournament: F has interesting roots.

  Question: Do the roots cluster on a specific curve?
  Lee-Yang theorem: for ferromagnetic Ising, roots are on unit circle.
  What about tournament F-polynomials?
""")

n = 5
tournaments = get_all_tournaments(n)

# Group by isomorphism class (H, score seq)
def score_seq(n, A):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

representatives = {}
for bits, A, H in tournaments:
    ss = score_seq(n, A)
    key = (H, ss)
    if key not in representatives:
        representatives[key] = (bits, A, H)

print(f"\n  n=5 F-polynomial zeros:")
for key in sorted(representatives.keys()):
    bits, A, H = representatives[key]

    # Compute F-polynomial
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

    # Find nonzero degree
    deg = max(k for k in range(n) if F[k] > 0)

    if deg == 0:
        print(f"    H={H:4d}, ss={key[1]}: F = {F[0]} (constant, no roots)")
        continue

    # Find roots using companion matrix eigenvalues (without numpy)
    # For degree <= 4, use explicit formulas or just evaluate

    # Evaluate F on unit circle
    print(f"    H={H:4d}, ss={key[1]}: F = {F[:deg+1]}")

    # Check F at key points
    F_neg1 = sum(F[k] * (-1)**k for k in range(n))
    omega = cmath.exp(2j * cmath.pi / 3)
    F_omega = sum(F[k] * omega**k for k in range(n))
    F_i = sum(F[k] * (1j)**k for k in range(n))

    # F(x) = 0 on unit circle? Check at many points
    zero_count = 0
    near_zero_angles = []
    for theta_deg in range(360):
        theta = theta_deg * cmath.pi / 180
        z = cmath.exp(1j * theta)
        Fz = sum(F[k] * z**k for k in range(n))
        if abs(Fz) < 0.5:
            near_zero_angles.append(theta_deg)

    if near_zero_angles:
        print(f"           Near-zero on unit circle at angles: {near_zero_angles}")

    print(f"           F(-1)={F_neg1}, |F(ω)|²={abs(F_omega)**2:.0f}, |F(i)|²={abs(F_i)**2:.0f}")

# ======================================================================
# INVESTIGATION 6: H AS MATRIX PERMANENT
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 6: H VIA MATRIX PERMANENT")
print("=" * 70)

print("""
  H(T) counts Hamiltonian paths. A Hamiltonian path is a permutation
  sigma such that A[sigma(1)][sigma(2)] * A[sigma(2)][sigma(3)] * ... = 1.

  Define the PATH MATRIX P where P[i][j] = A[i][j] for the tournament T.
  Then H = sum over all n! permutations sigma of prod_{k=1}^{n-1} P[sigma(k)][sigma(k+1)].

  This is NOT a standard permanent (which would be prod P[i][sigma(i)]).
  Instead, it's a "path permanent" or "chain permanent."

  But we CAN relate it to the permanent via:
  H = sum_{start v} sum_{end w} perm of the transfer matrix...

  Actually, the HP count from vertex v to vertex w is given by
  the (v,w) entry of the matrix (A)^{n-1} where we sum over
  paths avoiding repeated vertices (the permanent, not the power).

  More precisely: define the "HP matrix" M[v][w] where
  M[v][w] = number of Hamiltonian paths from v to w.
  Then H = sum_{v,w} M[v][w].

  M can be computed via inclusion-exclusion (permanent of submatrices).

  Question: Is det(M) related to H in a simple way?
  What about the eigenvalues of M?
""")

for n in [3, 4, 5]:
    m_val = n * (n-1) // 2
    t = get_all_tournaments(n)

    # Compute M[v][w] for each tournament
    print(f"\n  n={n}:")

    # Pick representatives
    H_seen = set()
    for bits, A, H in t:
        if H in H_seen:
            continue
        H_seen.add(H)

        # Compute M[v][w]
        M = [[0]*n for _ in range(n)]
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp[mask][v] == 0:
                    continue
                if mask == (1 << n) - 1:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]

        full = (1 << n) - 1
        # M[v][w]: need to track start vertex too
        # Recompute with start tracking
        # dp[mask][v][start]
        M = [[0]*n for _ in range(n)]
        for start in range(n):
            dp2 = [[0]*n for _ in range(1 << n)]
            dp2[1 << start][start] = 1
            for mask in range(1, 1 << n):
                for v in range(n):
                    if not (mask & (1 << v)):
                        continue
                    if dp2[mask][v] == 0:
                        continue
                    if mask == full:
                        continue
                    for u in range(n):
                        if mask & (1 << u):
                            continue
                        if A[v][u]:
                            dp2[mask | (1 << u)][u] += dp2[mask][v]
            for w in range(n):
                M[start][w] = dp2[full][w]

        # Properties of M
        trace_M = sum(M[i][i] for i in range(n))
        sum_M = sum(M[i][j] for i in range(n) for j in range(n))

        # Is M symmetric?
        sym = all(M[i][j] == M[j][i] for i in range(n) for j in range(n))

        # Row sums and column sums
        row_sums = [sum(M[i][j] for j in range(n)) for i in range(n)]
        col_sums = [sum(M[i][j] for i in range(n)) for j in range(n)]

        # Determinant (for small n)
        if n == 3:
            det_M = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                    -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                    +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
        elif n == 4:
            # Cofactor expansion
            det_M = 0
            for j in range(4):
                minor = [[M[i][k] for k in range(4) if k != j] for i in range(1, 4)]
                det_minor = (minor[0][0]*(minor[1][1]*minor[2][2]-minor[1][2]*minor[2][1])
                           -minor[0][1]*(minor[1][0]*minor[2][2]-minor[1][2]*minor[2][0])
                           +minor[0][2]*(minor[1][0]*minor[2][1]-minor[1][1]*minor[2][0]))
                det_M += ((-1)**j) * M[0][j] * det_minor
        else:
            det_M = "skipped"

        print(f"    H={H:4d}: tr(M)={trace_M}, sum(M)={sum_M}=H, "
              f"sym={'Y' if sym else 'N'}, det(M)={det_M}")
        if n <= 4:
            print(f"           row_sums={row_sums}, col_sums={col_sums}")

# ======================================================================
# INVESTIGATION 7: HÖLDER MEANS OF H
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 7: POWER MEANS AND MOMENTS OF H")
print("=" * 70)

print("""
  The power mean M_p = (E[H^p])^{1/p} for various p.
  M_1 = E[H] (arithmetic mean)
  M_2 = sqrt(E[H^2]) (quadratic mean)
  M_{-1} = 1/E[1/H] (harmonic mean)
  M_0 = exp(E[log H]) (geometric mean)

  Question: Does M_p / M_1 approach a simple constant as n grows?
""")

for n in [3, 4, 5, 6]:
    t = get_all_tournaments(n)
    H_vals = [H for _, _, H in t]
    N = len(H_vals)

    M1 = sum(H_vals) / N
    M2 = (sum(h**2 for h in H_vals) / N) ** 0.5
    M3 = (sum(h**3 for h in H_vals) / N) ** (1/3)
    M_neg1 = N / sum(1/h for h in H_vals)
    M0 = math.exp(sum(math.log(h) for h in H_vals) / N)

    print(f"\n  n={n}:")
    print(f"    M_1 (arithmetic) = {M1:.4f}")
    print(f"    M_0 (geometric)  = {M0:.4f}, ratio to M_1 = {M0/M1:.6f}")
    print(f"    M_{'{-1}'} (harmonic)  = {M_neg1:.4f}, ratio to M_1 = {M_neg1/M1:.6f}")
    print(f"    M_2 (quadratic)  = {M2:.4f}, ratio to M_1 = {M2/M1:.6f}")
    print(f"    M_3 (cubic)      = {M3:.4f}, ratio to M_1 = {M3/M1:.6f}")

    # Coefficient of variation
    var_H = sum((h - M1)**2 for h in H_vals) / N
    cv = var_H**0.5 / M1
    print(f"    CV = {cv:.6f}")

    # Skewness and kurtosis
    m3 = sum((h - M1)**3 for h in H_vals) / N
    m4 = sum((h - M1)**4 for h in H_vals) / N
    skew = m3 / var_H**1.5
    kurt = m4 / var_H**2 - 3
    print(f"    Skewness = {skew:.6f}, Excess kurtosis = {kurt:.6f}")

    if n == 6:
        break

# ======================================================================
# INVESTIGATION 8: H mod (2^k - 1) — MERSENNE CONNECTION
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 8: H mod MERSENNE NUMBERS")
print("=" * 70)

print("""
  Mersenne numbers M_k = 2^k - 1.
  M_1 = 1, M_2 = 3, M_3 = 7, M_4 = 15, M_5 = 31, M_6 = 63, M_7 = 127.

  These are related to 2^k via 2^k ≡ 1 mod M_k.

  Since H = I(Omega, 2) = sum i_k * 2^k:
  H mod M_k = sum i_j * 2^j mod (2^k - 1)
            = sum i_j * (2^j mod M_k)
            = sum i_j * 2^{j mod k}  (since 2^k ≡ 1 mod M_k)
            = sum_{r=0}^{k-1} 2^r * sum_{j≡r mod k} i_j

  This is a "DFT-like" projection of the independence number sequence.
""")

for n in [3, 4, 5, 6]:
    t = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in t))

    print(f"\n  n={n}:")
    for k in range(2, 8):
        Mk = 2**k - 1
        residues = set(h % Mk for h in H_set)
        if 0 not in residues:
            print(f"    H mod M_{k}={Mk:3d}: 0 is FORBIDDEN ({len(residues)}/{Mk} residues)")
        else:
            print(f"    H mod M_{k}={Mk:3d}: 0 is achievable ({len(residues)}/{Mk} residues)")

    if n == 6:
        break

# ======================================================================
# INVESTIGATION 9: DISTANCE BETWEEN H-CLASSES IN HYPERCUBE
# ======================================================================
print("\n" + "=" * 70)
print("INVESTIGATION 9: INTER-CLASS DISTANCES IN HYPERCUBE")
print("=" * 70)

n = 5
m_val = n * (n-1) // 2
t = get_all_tournaments(n)
H_map = {bits: H for bits, _, H in t}

# Group by H
H_groups = defaultdict(list)
for bits, _, H in t:
    H_groups[H].append(bits)

H_vals = sorted(H_groups.keys())

# For each pair of H classes, compute minimum Hamming distance
print(f"\n  n={n}: minimum Hamming distance between H classes:")
for i, h1 in enumerate(H_vals):
    for h2 in H_vals[i+1:]:
        # Sample if too many
        g1 = H_groups[h1]
        g2 = H_groups[h2]

        min_dist = m_val + 1
        for b1 in g1[:50]:  # sample
            for b2 in g2[:50]:
                d = bin(b1 ^ b2).count('1')
                if d < min_dist:
                    min_dist = d

        if min_dist <= 2:
            print(f"    d(H={h1}, H={h2}) = {min_dist}")

# Also: what is the minimum number of arc flips to change H by exactly 2?
print(f"\n  Minimum flips to change H by 2:")
for h in H_vals:
    h_plus_2 = h + 2
    if h_plus_2 not in H_groups:
        continue
    min_flips = m_val + 1
    for b1 in H_groups[h][:100]:
        for b2 in H_groups[h_plus_2][:100]:
            d = bin(b1 ^ b2).count('1')
            min_flips = min(min_flips, d)
    print(f"    H={h} -> H={h_plus_2}: min {min_flips} flips")

print("\n" + "=" * 70)
print("DONE — PART 3 DEEP INVESTIGATIONS")
print("=" * 70)
