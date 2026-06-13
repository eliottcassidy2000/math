"""
tropical_fundamental.py — opus-2026-04-04-S7

THINKING TROPICALLY, SEEING WHAT IS MOST FUNDAMENTAL.

The core equation: H(T) = I(Omega(T), 2)
  H = Hamiltonian path count
  Omega = conflict graph (vertices = odd cycles, edges = shared vertices)
  I(G, x) = independence polynomial at x = 2

Everything flows from this. The multilinear polynomial, the sign rule,
the recursive preservation, the degree cap — all are shadows of the
conflict graph's structure as a FUNCTOR from tilings to graphs.

This script explores:
1. Newton polytope of H(t)
2. Tropical H under 2-adic valuation
3. Conflict graph evolution with tile flips
4. The dominant term / tropical permanent
5. Free energy and Legendre transform
6. The conflict graph as the fundamental object
"""
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v

# ==============================================================
# 1. NEWTON POLYTOPE
# ==============================================================
def newton_polytope(n):
    """Compute and analyze the Newton polytope of H(t)."""
    print(f"\n{'='*60}")
    print(f"NEWTON POLYTOPE of H(t), n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # Support = set of exponent vectors with nonzero coefficient
    support = []
    for S, c in coeffs.items():
        if c != 0:
            vec = tuple((S >> i) & 1 for i in range(m))
            support.append(vec)

    support = list(set(support))
    print(f"  Support size: {len(support)} / {2**m} possible monomials")

    # The Newton polytope NP(H) = conv(support) in R^m
    # For small m, compute its properties

    # Dimension: rank of the affine span
    if len(support) > 1:
        arr = np.array(support, dtype=float)
        centered = arr - arr[0]
        rank = np.linalg.matrix_rank(centered)
        print(f"  Affine dimension of NP(H): {rank}")
        print(f"  Ambient dimension: {m}")
        print(f"  Codimension: {m - rank}")
    else:
        print(f"  Trivial (single point)")
        return

    # Does NP(H) = the full hypercube [0,1]^m?
    # This requires: the origin (0,...,0) and all unit vectors are in the support
    # AND all vertices of [0,1]^m... no, just need the support to be "full"
    has_origin = (0,)*m in support
    has_all_ones = (1,)*m in support
    has_all_unit = all(tuple(1 if j==i else 0 for j in range(m)) in support
                       for i in range(m))
    print(f"  Contains origin: {has_origin}")
    print(f"  Contains (1,...,1): {has_all_ones}")
    print(f"  Contains all unit vectors: {has_all_unit}")

    # Check: is the support a DOWNWARD CLOSED set (simplicial complex)?
    # We already know it's NOT (from S5), but let's check the "upward" version
    # Is NP(H) = the Boolean lattice minus some anti-chain?

    # More useful: what is the MINIMUM DEGREE vertex set?
    # I.e., what are the degree-0 and degree-1 exponent vectors in the support?
    by_deg = defaultdict(list)
    for vec in support:
        by_deg[sum(vec)].append(vec)

    print(f"\n  Support by degree:")
    max_deg = max(by_deg.keys())
    for d in range(max_deg + 1):
        count = len(by_deg.get(d, []))
        from math import comb
        total = comb(m, d) if d <= m else 0
        print(f"    d={d}: {count}/{total} = {count/total:.3f}" if total > 0 else f"    d={d}: {count}")

    # The MAXIMAL DEGREE is the degree cap 2*floor((n-1)/2)
    max_support_deg = max(sum(v) for v in support)
    print(f"\n  Max degree in support: {max_support_deg}")
    print(f"  Expected degree cap: {2*((n-1)//2)}")

    # VERTICES of the Newton polytope (extreme points)
    # For Boolean polytopes, the vertices are exactly the support points that
    # are NOT convex combinations of other support points.
    # This is hard to compute in general, but for small m we can check.
    if m <= 10:
        # A support point v is a vertex iff it's NOT in conv(support \ {v})
        # Quick heuristic: check if removing it changes the polytope
        # For now, just count vertices as support points on the "boundary"
        # of the polytope — those that are maximal in some direction.
        pass

    # FACES of the Newton polytope
    # The face lattice encodes how terms interact tropically
    # For each direction w in R^m, the face F_w = argmax_{v in support} w·v
    # is a face of NP(H).

    # Test a few interesting directions:
    directions = {
        "all-ones": np.ones(m),
        "neg-all-ones": -np.ones(m),
        "skip-weighted": np.array([tiles[i][0]-tiles[i][1] for i in range(m)], dtype=float),
        "first-tile": np.array([1 if i==0 else 0 for i in range(m)], dtype=float),
        "last-tile": np.array([1 if i==m-1 else 0 for i in range(m)], dtype=float),
    }

    arr = np.array(support, dtype=float)
    print(f"\n  Face analysis (vertices maximizing w·v):")
    for name, w in directions.items():
        dots = arr @ w
        max_dot = np.max(dots)
        face_indices = np.where(np.abs(dots - max_dot) < 1e-10)[0]
        face_vecs = [support[i] for i in face_indices]
        face_degs = [sum(v) for v in face_vecs]
        print(f"    w={name}: face has {len(face_vecs)} vertices, "
              f"degrees {sorted(set(face_degs))}")

    return support, by_deg

# ==============================================================
# 2. TROPICAL H — piecewise-linear tropicalization
# ==============================================================
def tropical_H(n):
    """
    Tropicalize H(t) using the 2-adic valuation.

    trop(H)(t) = min_S (v_2(c_S) + sum_{i in S} t_i)  [min convention]

    This is a piecewise-linear function on R^m. The "tropical variety"
    is the locus where the minimum is achieved by 2+ terms.

    For our BINARY domain t in {0,1}^m:
    trop(H)(t) = min_S (v_2(c_S) + |S ∩ supp(t)|)
    where supp(t) = {i : t_i = 1}.
    """
    print(f"\n{'='*60}")
    print(f"TROPICAL H(t) via 2-ADIC VALUATION, n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)
    N = 1 << m

    # Compute v_2(c_S) for each nonzero coefficient
    trop_coeffs = {}
    for S, c in coeffs.items():
        if c != 0:
            trop_coeffs[S] = v2(abs(c))

    print(f"  Nonzero terms: {len(trop_coeffs)}")

    # For each tiling t, compute trop(H)(t) = min_S (v_2(c_S) + popcount(S & t))
    # and also v_2(H(t))
    trop_vals = []
    actual_v2 = []
    agreement = 0

    for t in range(N):
        # Tropical value
        min_val = float('inf')
        min_count = 0
        for S, v in trop_coeffs.items():
            val = v + bin(S & t).count('1')  # v_2(c_S * t^S) = v_2(c_S) + #{i in S : t_i=1}
            # Wait, this isn't right. In tropical (min, +):
            # trop(c_S * prod t_i) = val(c_S) + sum val(t_i)
            # For t_i in {0,1}: val(0) = infinity, val(1) = 0
            # So the tropical monomial is val(c_S) + 0*(#ones in S cap t) + inf*(#zeros in S\t)
            # = val(c_S) if S subset supp(t), else infinity
            # So: trop(H)(t) = min_{S subset supp(t)} v_2(c_S)

            if S & ~t:  # S has bits not in t → infinity (t_i=0 for some i in S)
                continue  # this term is infinity
            # S subset supp(t)
            val = v
            if val < min_val:
                min_val = val
                min_count = 1
            elif val == min_val:
                min_count += 1

        trop_vals.append((min_val, min_count))
        actual_v2.append(v2(int(H[t])))

        if min_val == v2(int(H[t])):
            agreement += 1

    print(f"  trop(H)(t) = v_2(H(t))? Agreement: {agreement}/{N} "
          f"= {agreement/N:.3f}")

    # The tropical value should equal v_2(H(t)) when there's NO cancellation.
    # When it doesn't match, there's 2-adic cancellation among terms.

    # Distribution of tropical values
    trop_dist = Counter(tv for tv, _ in trop_vals)
    print(f"\n  Tropical value distribution:")
    for v in sorted(trop_dist.keys()):
        actual_with_this = sum(1 for i in range(N) if actual_v2[i] == v)
        print(f"    v_2 = {v}: trop predicts {trop_dist[v]}, actual {actual_with_this}")

    # Multiplicity: how many tilings have the tropical min achieved by 2+ terms?
    multi = sum(1 for _, c in trop_vals if c >= 2)
    print(f"\n  Tropical variety (min achieved by 2+ terms): {multi}/{N} "
          f"= {multi/N:.3f}")

    # The DOMINANT TERM for each tiling: the S with minimum v_2(c_S) among S subset supp(t)
    # For the transitive (t=0): only S=empty contributes → trop = v_2(c_0) = v_2(1) = 0
    # For all-backward (t=all 1): min over ALL S of v_2(c_S)
    print(f"\n  Dominant term analysis:")
    print(f"    Transitive (t=0): trop={trop_vals[0][0]}, "
          f"v_2(H)={actual_v2[0]}, mult={trop_vals[0][1]}")
    print(f"    All-backward (t=all 1): trop={trop_vals[N-1][0]}, "
          f"v_2(H)={actual_v2[N-1]}, mult={trop_vals[N-1][1]}")

    return trop_vals, actual_v2

# ==============================================================
# 3. CONFLICT GRAPH EVOLUTION
# ==============================================================
def conflict_graph_evolution(n):
    """
    THE MOST FUNDAMENTAL VIEW: track Omega(T) as tiles flip.

    Start from transitive (t=0, all forward). Omega(transitive) = empty graph
    (no odd cycles). As tiles flip to backward, odd cycles appear.

    Track:
    - Number of odd cycles (vertices of Omega)
    - Number of edges in Omega (pairs sharing a vertex)
    - Independence number alpha(Omega) (since H = I(Omega, 2))
    - The structure of Omega itself
    """
    print(f"\n{'='*60}")
    print(f"CONFLICT GRAPH EVOLUTION n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)

    def find_odd_cycles(T_matrix, n):
        """Find all directed odd cycles in tournament T."""
        cycles = []
        # Brute force for small n: check all subsets of size 3, 5, ...
        for k in range(3, n+1, 2):
            for verts in combinations(range(n), k):
                # Check if these vertices form a directed cycle
                # Try all cyclic orderings
                from itertools import permutations
                for perm in permutations(verts):
                    is_cycle = True
                    for i in range(k):
                        u, v = perm[i], perm[(i+1)%k]
                        if T_matrix[u][v] != 1:
                            is_cycle = False
                            break
                    if is_cycle:
                        # Canonical: min rotation
                        rotations = [perm[i:]+perm[:i] for i in range(k)]
                        canon = min(rotations)
                        if perm == canon:
                            cycles.append(set(perm))
                        break  # found a valid ordering for these vertices
        return cycles

    def build_conflict_graph(cycles):
        """Build conflict graph: edge iff cycles share a vertex."""
        nc = len(cycles)
        edges = []
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i] & cycles[j]:
                    edges.append((i, j))
        return edges

    def independence_number(nc, edges):
        """Compute max independent set size (brute force for small graphs)."""
        if nc == 0:
            return 0
        adj = [set() for _ in range(nc)]
        for i, j in edges:
            adj[i].add(j)
            adj[j].add(i)

        max_is = 0
        for size in range(nc, 0, -1):
            for subset in combinations(range(nc), size):
                if all(j not in adj[i] for i in subset for j in subset if j > i):
                    return size
        return 0

    # Start from transitive, flip tiles one by one
    print(f"  Flipping tiles from transitive (skip-ordered: largest first):")

    # Order tiles by skip (largest first — these create the most cycles)
    tile_order = sorted(range(m), key=lambda i: tiles[i][0]-tiles[i][1], reverse=True)

    current_tiling = 0
    for step, tile_idx in enumerate(tile_order):
        current_tiling |= (1 << tile_idx)
        T_matrix = tiling_to_tournament(n, tiles, current_tiling)

        # Count odd cycles
        if n <= 6:  # Only feasible for small n
            cycles = find_odd_cycles(T_matrix, n)
            edges = build_conflict_graph(cycles)
            nc = len(cycles)
            ne = len(edges)

            # Independence polynomial at x=2: I(Omega, 2) = sum_{k=0}^{alpha} alpha_k * 2^k
            # For small Omega, compute directly
            H_val = count_hamiltonian_paths(T_matrix, n)

            print(f"    Step {step+1}: flip tile {tiles[tile_idx]} (skip={tiles[tile_idx][0]-tiles[tile_idx][1]})")
            print(f"      Omega: {nc} cycles, {ne} edges, "
                  f"H={H_val}")

            if step < 3 and nc <= 10:
                # Show the cycles
                for i, cyc in enumerate(cycles):
                    print(f"        cycle {i}: {sorted(cyc)} (len {len(cyc)})")
        else:
            H_val = count_hamiltonian_paths(tiling_to_tournament(n, tiles, current_tiling), n)
            print(f"    Step {step+1}: flip {tiles[tile_idx]} → H={H_val}")

        if step >= 8:
            print(f"    ... (truncated)")
            break

    # Also: evolution along a RANDOM path
    print(f"\n  Random flip sequence (seed=42):")
    np.random.seed(42)
    perm = np.random.permutation(m)
    current_tiling = 0
    h_sequence = [1]  # H(transitive)
    for step, tile_idx in enumerate(perm):
        current_tiling |= (1 << tile_idx)
        T_matrix = tiling_to_tournament(n, tiles, current_tiling)
        H_val = count_hamiltonian_paths(T_matrix, n)
        h_sequence.append(H_val)

    print(f"    H sequence (first 15): {h_sequence[:15]}")
    print(f"    H diffs: {[h_sequence[i+1]-h_sequence[i] for i in range(min(14,len(h_sequence)-1))]}")

    # Monotonicity analysis: is H always increasing?
    increasing = all(h_sequence[i+1] >= h_sequence[i] for i in range(len(h_sequence)-1))
    print(f"    Always increasing? {increasing}")
    if not increasing:
        decreases = [(i, h_sequence[i], h_sequence[i+1])
                     for i in range(len(h_sequence)-1) if h_sequence[i+1] < h_sequence[i]]
        print(f"    Decreases at steps: {decreases[:5]}")

# ==============================================================
# 4. DOMINANT HAMILTONIAN PATH — the tropical permanent
# ==============================================================
def tropical_permanent(n):
    """
    For each tournament, find the "dominant" Hamiltonian path —
    the one that contributes most to H.

    In the OCF picture: H = sum over independent cycle sets.
    The dominant term is the largest independent set weighted by 2^|I|.

    For the tiling polynomial: the "dominant monomial" is the one with
    the largest |c_S| * (product of t_i) term.
    """
    print(f"\n{'='*60}")
    print(f"DOMINANT TERMS AND TROPICAL PERMANENT n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)
    N = 1 << m

    # For each tiling t, find the DOMINANT coefficient:
    # the c_S with S subset supp(t) that has the largest |c_S|
    dominant_S_dist = Counter()
    for t in range(N):
        best_S = None
        best_abs = 0
        for S, c in coeffs.items():
            if S & ~t:  # S not subset of supp(t)
                continue
            if abs(c) > best_abs:
                best_abs = abs(c)
                best_S = S
        if best_S is not None:
            dominant_S_dist[best_S] += 1

    print(f"  Most frequent dominant terms:")
    for S, count in dominant_S_dist.most_common(10):
        tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
        c = coeffs.get(S, 0)
        deg = bin(S).count('1')
        print(f"    S={tile_set} (deg {deg}): dominant for {count}/{N} tilings, c={c}")

    # What fraction of H(t) comes from the constant term (c_0 = 1)?
    const_fraction = [1.0 / int(H[t]) if H[t] > 0 else 0 for t in range(N)]
    print(f"\n  Constant term fraction 1/H(t):")
    print(f"    mean = {np.mean(const_fraction):.4f}")
    print(f"    min = {np.min(const_fraction):.6f} (at H={np.max(H)})")
    print(f"    max = {np.max(const_fraction):.6f} (at H={np.min(H)})")

    # What fraction comes from linear terms?
    lin_fraction = []
    for t in range(N):
        lin_sum = 1  # constant
        for k in range(m):
            if t & (1 << k):
                S = 1 << k
                lin_sum += coeffs.get(S, 0)
        lin_fraction.append(lin_sum / int(H[t]) if H[t] > 0 else 0)

    print(f"\n  Linear approximation (1 + sum c_k*t_k) / H(t):")
    print(f"    mean = {np.mean(lin_fraction):.4f}")
    print(f"    min = {np.min(lin_fraction):.4f}")
    print(f"    max = {np.max(lin_fraction):.4f}")

# ==============================================================
# 5. FREE ENERGY AND LEGENDRE TRANSFORM
# ==============================================================
def free_energy_legendre(n):
    """
    F(lambda) = log H(e^{lambda_1}, ..., e^{lambda_m})
    = log sum_t H(t) * exp(lambda . t)

    This is the cumulant generating function (free energy) of the
    Boltzmann distribution. Its Legendre transform gives the entropy:
    S(mu) = inf_lambda (lambda . mu - F(lambda))

    We compute F along the diagonal: F(lambda, ..., lambda) = log H(e^lambda, ..., e^lambda)
    using the univariate polynomial H(p) = sum_d S_d p^d.
    """
    print(f"\n{'='*60}")
    print(f"FREE ENERGY AND LEGENDRE TRANSFORM n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # Univariate H(p) = sum_d S_d p^d
    by_deg = defaultdict(int)
    for S, c in coeffs.items():
        by_deg[bin(S).count('1')] += c

    max_d = max(by_deg.keys())
    S_d = [by_deg.get(d, 0) for d in range(max_d + 1)]

    print(f"  Degree sums S_d = {S_d}")

    # F(lambda) = log(sum_d S_d * exp(d*lambda))
    # This is the log-partition function along the diagonal

    lambdas = np.linspace(-3, 3, 200)
    F_vals = []
    for lam in lambdas:
        val = sum(S_d[d] * np.exp(d * lam) for d in range(max_d + 1) if S_d[d] != 0)
        if val > 0:
            F_vals.append(np.log(val))
        else:
            F_vals.append(float('nan'))

    F_vals = np.array(F_vals)

    # Key points
    F_at_0 = np.log(sum(S_d))  # = log(H(1,...,1))
    F_at_neg_inf = np.log(abs(S_d[0]))  # → log(1) = 0

    print(f"\n  F(0) = log(H(1)) = log({sum(S_d)}) = {F_at_0:.4f}")
    print(f"  F(-inf) = log(c_0) = log(1) = 0")

    # The derivative F'(lambda) = E[degree] under the tilted distribution
    # F'(0) = sum_d d*S_d / sum_d S_d = E[degree] under uniform-on-tilings
    dF_at_0 = sum(d * S_d[d] for d in range(max_d + 1)) / sum(S_d)
    print(f"  F'(0) = E[degree] = {dF_at_0:.4f}")
    print(f"    (Expected number of backward tiles per tiling)")
    print(f"    Uniform expectation: m/2 = {m/2}")

    # Second derivative = variance of degree
    d2F_at_0 = (sum(d**2 * S_d[d] for d in range(max_d+1)) / sum(S_d)) - dF_at_0**2
    print(f"  F''(0) = Var[degree] = {d2F_at_0:.4f}")

    # The Legendre transform: S(mu) = sup_lambda (lambda*mu - F(lambda))
    # For the diagonal: mu = E[degree] under the tilted distribution
    # mu ranges from 0 (lambda → -inf, only transitive) to max_d (lambda → +inf)

    print(f"\n  Legendre transform S(mu):")
    mus = np.linspace(0.01, max_d - 0.01, 50)
    for mu in [0.5, 1.0, m/4, m/2, 3*m/4, m-1]:
        if mu > max_d:
            continue
        # Find lambda* such that F'(lambda*) = mu
        # Numerical: bisection
        lo, hi = -10, 10
        for _ in range(100):
            mid = (lo + hi) / 2
            # F'(mid) = sum d*S_d*exp(d*mid) / sum S_d*exp(d*mid)
            num = sum(d * S_d[d] * np.exp(d * mid) for d in range(max_d+1) if S_d[d] != 0)
            den = sum(S_d[d] * np.exp(d * mid) for d in range(max_d+1) if S_d[d] != 0)
            if den == 0:
                break
            fprime = num / den
            if fprime < mu:
                lo = mid
            else:
                hi = mid
        lam_star = (lo + hi) / 2
        F_star = np.log(max(1e-300, sum(S_d[d]*np.exp(d*lam_star) for d in range(max_d+1) if S_d[d]!=0)))
        S_val = lam_star * mu - F_star
        print(f"    mu={mu:.1f}: lambda*={lam_star:.3f}, S(mu)={S_val:.4f}")

    # The entropy at mu = m/2 (uniform) should be related to log(H(1)/2^m)
    # H(1)/2^m = average H value / 2^m = fraction of HPs that survive
    print(f"\n  Key entropy values:")
    print(f"    S(0+) = log(c_0) - F(-inf) = 0 (only transitive)")
    print(f"    S(m/2) ≈ {m/2}*0 - log({sum(S_d)}) = {-F_at_0:.4f} bits")
    print(f"    Boltzmann entropy: {np.log2(sum(S_d)):.4f} bits")

# ==============================================================
# 6. THE IRREDUCIBLE CORE
# ==============================================================
def irreducible_core(n):
    """
    What is the SIMPLEST description of H(t)?

    Candidate 1: The conflict graph Omega(T(t)).
    H(t) = I(Omega(t), 2). If we know Omega, we know H.

    Candidate 2: The cycle set.
    Each tile flip creates/destroys odd cycles. The cycles are determined
    by the tiling. H is determined by the cycle packing structure.

    Candidate 3: The 3-cycle count + correction.
    At small n, H depends mostly on the number of 3-cycles (c_3).
    Higher cycles add corrections.

    Let's test: how well does c_3(T) predict H(T)?
    """
    print(f"\n{'='*60}")
    print(f"THE IRREDUCIBLE CORE n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m

    # Compute H and c_3 for each tiling
    h_vals = []
    c3_vals = []

    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        h = count_hamiltonian_paths(T, n)
        # Count 3-cycles
        c3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if T[a][b] + T[b][c] + T[c][a] == 3 or T[a][c] + T[c][b] + T[b][a] == 3:
                        c3 += 1
        h_vals.append(h)
        c3_vals.append(c3)

    h_vals = np.array(h_vals)
    c3_vals = np.array(c3_vals)

    corr = np.corrcoef(h_vals, c3_vals)[0,1]
    print(f"  Corr(H, c_3) = {corr:.4f}")

    # Linear regression H = a*c_3 + b
    A = np.column_stack([c3_vals, np.ones(N)])
    beta = np.linalg.lstsq(A, h_vals, rcond=None)[0]
    pred = A @ beta
    r2 = 1 - np.sum((h_vals - pred)**2) / np.sum((h_vals - np.mean(h_vals))**2)
    print(f"  H ≈ {beta[0]:.3f}*c_3 + {beta[1]:.3f}, R² = {r2:.4f}")

    # Score sequence as predictor
    # The score sequence s = (s_0, ..., s_{n-1}) where s_v = #{u : v beats u}
    score_diffs = []
    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        scores = sorted([sum(T[v]) for v in range(n)])
        score_diffs.append(scores[-1] - scores[0])  # score spread

    score_diffs = np.array(score_diffs)
    corr_score = np.corrcoef(h_vals, score_diffs)[0,1]
    print(f"\n  Corr(H, score_spread) = {corr_score:.4f}")

    # Number of backward tiles (Hamming weight) as predictor
    hw = np.array([bin(t).count('1') for t in range(N)])
    corr_hw = np.corrcoef(h_vals, hw)[0,1]
    print(f"  Corr(H, #backward_tiles) = {corr_hw:.4f}")

    # Weighted backward count (sum of skip values for backward tiles)
    wt_hw = np.array([sum(tiles[k][0]-tiles[k][1] for k in range(m) if t & (1<<k)) for t in range(N)])
    corr_wt = np.corrcoef(h_vals, wt_hw)[0,1]
    print(f"  Corr(H, weighted_backward) = {corr_wt:.4f}")

    # Linear model with BOTH c_3 and weighted backward
    A3 = np.column_stack([c3_vals, wt_hw, np.ones(N)])
    beta3 = np.linalg.lstsq(A3, h_vals, rcond=None)[0]
    pred3 = A3 @ beta3
    r2_3 = 1 - np.sum((h_vals - pred3)**2) / np.sum((h_vals - np.mean(h_vals))**2)
    print(f"\n  H ≈ {beta3[0]:.3f}*c_3 + {beta3[1]:.3f}*wt_hw + {beta3[2]:.3f}")
    print(f"  R² = {r2_3:.4f}")

    # THE FUNDAMENTAL QUESTION: does (c_3, alpha(Omega)) determine H?
    # alpha = max independent set in conflict graph
    # H = sum_{k=0}^{alpha} alpha_k * 2^k

    # Actually, alpha_1 = number of odd cycles = |V(Omega)|
    # This relates to H mod 4: H ≡ 1 + 2*alpha_1 (mod 4)
    # So alpha_1 determines H mod 4.

    # Does alpha_1 predict H well?
    # alpha_1 = number of odd directed cycles = f(tournament)
    # This is related to c_3 + c_5 + c_7 + ...

    # For n=5: all odd cycles are 3-cycles and 5-cycles
    # c_5 at n=5: each 5-tournament has at most 1 directed 5-cycle... actually
    # a tournament on 5 vertices has C(5,5)*4!/5 = 24/5... no. Number of directed
    # 5-cycles = (n-1)!/2 per vertex set (if tournament is a 5-cycle).

    print(f"\n  THE FUNDAMENTAL HIERARCHY:")
    print(f"    Level 0: H(t) mod 2 = 1 always (Redei)")
    print(f"    Level 1: H(t) mod 4 = 1 + 2*alpha_1(t) mod 4")
    print(f"    Level 2: H(t) mod 8 depends on alpha_1, alpha_2")
    print(f"    Level k: H(t) mod 2^(k+1) depends on alpha_0,...,alpha_k")
    print(f"    Full: H(t) = sum_k alpha_k * 2^k")
    print(f"")
    print(f"  The OCF gives H in terms of the INDEPENDENCE SEQUENCE of Omega.")
    print(f"  The multilinear polynomial expresses this in tile coordinates.")
    print(f"  The sign rule, recursive preservation, etc. are all consequences")
    print(f"  of how the conflict graph Omega depends on tile flips.")

# ==============================================================
# MAIN
# ==============================================================
if __name__ == "__main__":
    for n in [5, 6, 7]:
        newton_polytope(n)

    for n in [5, 6]:
        tropical_H(n)

    for n in [5]:
        conflict_graph_evolution(n)

    for n in [5, 6]:
        tropical_permanent(n)

    for n in [5, 6, 7]:
        free_energy_legendre(n)

    for n in [5, 6]:
        irreducible_core(n)
