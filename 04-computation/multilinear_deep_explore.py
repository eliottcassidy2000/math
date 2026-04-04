"""
multilinear_deep_explore.py — opus-2026-04-04-S3

Deep creative exploration of the multilinear polynomial H(t₁,...,tₘ).

INVESTIGATIONS:
1. Walsh/Fourier spectrum (convert to ±1 variables)
2. Quadratic form matrix Q[i,j] = c_{ij} — eigenvalues, rank
3. Coefficient support topology — simplicial complex?
4. 2-adic valuation patterns by degree
5. Gradient structure — ∂H/∂t_k
6. Tile interaction graph from coefficient support
7. Coefficient statistics and distribution
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

# Import the base computation
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v

def walsh_transform(H_values, m):
    """
    Compute Walsh-Hadamard transform.

    H(t) as function on {0,1}^m → convert to f(x) on {±1}^m via x_i = 1-2t_i.
    Walsh coefficient ĥ_S = (1/2^m) Σ_x f(x) Π_{i∈S} x_i

    Equivalently: ĥ_S = (1/2^m) Σ_t H(t) Π_{i∈S} (1-2t_i)
    """
    N = 1 << m
    walsh = {}
    for S in range(N):
        val = 0
        for t in range(N):
            # Π_{i∈S} (1-2t_i) = (-1)^{popcount(S&t)}
            sign = (-1) ** bin(S & t).count('1')
            val += H_values[t] * sign
        walsh[S] = Fraction(val, N)
    return walsh

def multilinear_to_walsh(coeffs, m):
    """
    Convert multilinear coefficients to Walsh coefficients.

    H(t) = Σ_S c_S Π_{i∈S} t_i

    Substituting t_i = (1 - x_i)/2:
    H = Σ_S c_S Π_{i∈S} (1-x_i)/2
      = Σ_S c_S / 2^|S| Π_{i∈S} (1-x_i)

    The Walsh coefficients are:
    ĥ_T = Σ_{S ⊇ T} c_S (-1)^|T| / 2^|S| · C(|S|-|T|, ...) ...

    Actually let me just use the direct formula. More reliable.
    """
    pass  # Use walsh_transform directly

def analyze_walsh_spectrum(walsh, m, tiles):
    """Analyze the Walsh spectrum."""
    print(f"\n{'='*60}")
    print(f"WALSH SPECTRUM ANALYSIS (m={m} tiles)")
    print(f"{'='*60}")

    # Group by order (weight of S)
    by_order = defaultdict(list)
    for S, val in walsh.items():
        if val != 0:
            order = bin(S).count('1')
            by_order[order].append((S, val))

    total_nonzero = sum(len(v) for v in by_order.values())
    print(f"Total nonzero Walsh coefficients: {total_nonzero} / {1 << m}")

    for order in sorted(by_order.keys()):
        entries = by_order[order]
        vals = [v for _, v in entries]
        print(f"\n  Order {order}: {len(entries)} nonzero")
        # Show distinct values
        distinct = sorted(set(vals), key=lambda x: float(x))
        print(f"    Distinct values: {distinct[:20]}")

        # Check if all have same denominator
        denoms = set(v.denominator for v in vals)
        print(f"    Denominators: {denoms}")

        # Check 2-adic structure
        nums = [v.numerator for v in vals]
        min_v2_num = min(v2(abs(n)) for n in nums if n != 0) if nums else 0
        max_v2_denom = max(v2(d) for d in denoms)
        print(f"    Numerator min v₂: {min_v2_num}, Denominator max v₂: {max_v2_denom}")

    # Parseval's identity: Σ ĥ_S² = (1/2^m) Σ H(t)²
    parseval_walsh = sum(float(v)**2 for v in walsh.values())
    print(f"\n  Parseval check: Σ ĥ² = {parseval_walsh:.6f}")

    return by_order

def analyze_quadratic_form(coeffs, m, tiles, n):
    """Build and analyze the quadratic coefficient matrix Q."""
    print(f"\n{'='*60}")
    print(f"QUADRATIC FORM MATRIX (n={n}, m={m})")
    print(f"{'='*60}")

    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        i, j = bits[0], bits[1]
        Q[i][j] = c
        Q[j][i] = c  # Symmetric

    # Also put linear coefficients on diagonal (half of them for the quadratic form)
    # Actually no — the multilinear polynomial separate linear and quadratic terms.
    # Q[i,j] for i≠j is the pairwise interaction. Diagonal is 0 (no t_i² in multilinear).

    print(f"  Q matrix ({m}×{m}):")
    if m <= 15:
        # Print tile labels
        tile_labels = [f"({x},{y})" for x, y in tiles]
        print(f"  {'':>8s}", end="")
        for lbl in tile_labels:
            print(f" {lbl:>6s}", end="")
        print()
        for i in range(m):
            print(f"  {tile_labels[i]:>8s}", end="")
            for j in range(m):
                val = int(Q[i][j])
                if val == 0:
                    print(f"     .", end="")
                else:
                    print(f" {val:5d}", end="")
            print()

    # Eigenvalues
    eigvals = np.linalg.eigvalsh(Q)
    eigvals = sorted(eigvals, reverse=True)
    print(f"\n  Eigenvalues: {[f'{e:.4f}' for e in eigvals]}")

    # Rank
    rank = np.linalg.matrix_rank(Q)
    print(f"  Rank: {rank} / {m}")

    # Trace
    trace = np.trace(Q)
    print(f"  Trace: {trace}")

    # Determinant (if small enough)
    if m <= 15:
        det = np.linalg.det(Q)
        print(f"  Determinant: {det:.6f}")

    # Frobenius norm
    frob = np.linalg.norm(Q, 'fro')
    print(f"  Frobenius norm: {frob:.4f}")

    # Positive/negative eigenvalue count
    pos = sum(1 for e in eigvals if e > 1e-10)
    neg = sum(1 for e in eigvals if e < -1e-10)
    zero = m - pos - neg
    print(f"  Signature: ({pos}+, {neg}-, {zero} zero)")

    # Check: does Q have nice eigenvalues (powers of 2)?
    print(f"\n  Eigenvalue analysis:")
    for e in eigvals:
        if abs(e) > 1e-10:
            ratio_to_2 = abs(e) / 2.0
            print(f"    λ = {e:12.6f}, |λ|/2 = {ratio_to_2:.6f}, log₂|λ| = {np.log2(abs(e)):.4f}")

    return Q, eigvals

def analyze_support_topology(coeffs, m, tiles, n):
    """Analyze the support of nonzero coefficients as a simplicial complex."""
    print(f"\n{'='*60}")
    print(f"COEFFICIENT SUPPORT TOPOLOGY (n={n})")
    print(f"{'='*60}")

    # The support is {S : c_S ≠ 0}
    support = set()
    for S, c in coeffs.items():
        if c != 0 and S != 0:
            support.add(S)

    # Group by degree
    by_degree = defaultdict(set)
    for S in support:
        by_degree[bin(S).count('1')].add(S)

    max_deg = max(by_degree.keys()) if by_degree else 0
    print(f"  Support sizes by degree:")
    total_possible = 0
    for d in range(1, max_deg + 1):
        from math import comb
        possible = comb(m, d)
        actual = len(by_degree.get(d, set()))
        total_possible += possible
        frac = actual / possible if possible > 0 else 0
        print(f"    d={d}: {actual}/{possible} = {frac:.3f}")

    # Is it a simplicial complex? (closed under subsets)
    print(f"\n  Simplicial complex check:")
    violations = 0
    for S in support:
        deg = bin(S).count('1')
        if deg <= 1:
            continue
        # Check all codimension-1 faces
        for i in range(m):
            if S & (1 << i):
                face = S ^ (1 << i)  # Remove bit i
                if face != 0 and face not in support:
                    if violations < 5:
                        face_tiles = [tiles[j] for j in range(m) if face & (1 << j)]
                        S_tiles = [tiles[j] for j in range(m) if S & (1 << j)]
                        print(f"    VIOLATION: {S_tiles} in support but face {face_tiles} missing")
                    violations += 1

    if violations == 0:
        print(f"    ✓ Support IS a simplicial complex!")
    else:
        print(f"    ✗ {violations} violations — NOT a simplicial complex")

    # Check if the POSITIVE and NEGATIVE supports separately form complexes
    pos_support = {S for S, c in coeffs.items() if c > 0 and S != 0}
    neg_support = {S for S, c in coeffs.items() if c < 0 and S != 0}

    for name, supp in [("POSITIVE", pos_support), ("NEGATIVE", neg_support)]:
        viol = 0
        for S in supp:
            deg = bin(S).count('1')
            if deg <= 1:
                continue
            for i in range(m):
                if S & (1 << i):
                    face = S ^ (1 << i)
                    if face != 0 and face not in supp and face not in support:
                        viol += 1
        print(f"    {name} support ({len(supp)} sets): {'complex' if viol == 0 else f'{viol} face violations'}")

    # f-vector (face numbers)
    print(f"\n  f-vector: f = ({', '.join(str(len(by_degree.get(d, set()))) for d in range(1, max_deg+1))})")

    # Euler characteristic of support
    euler = sum((-1)**(d+1) * len(by_degree.get(d, set())) for d in range(1, max_deg+1))
    print(f"  Euler characteristic: χ = {euler}")

    return support, by_degree

def analyze_2adic(coeffs, m, tiles, n):
    """Analyze 2-adic valuation patterns."""
    print(f"\n{'='*60}")
    print(f"2-ADIC VALUATION STRUCTURE (n={n})")
    print(f"{'='*60}")

    by_degree = defaultdict(list)
    for S, c in coeffs.items():
        if c != 0:
            deg = bin(S).count('1')
            by_degree[deg].append((S, c))

    for deg in sorted(by_degree.keys()):
        entries = by_degree[deg]
        vals = [c for _, c in entries]
        v2_vals = [v2(abs(c)) for c in vals]
        min_v2 = min(v2_vals)
        max_v2 = max(v2_vals)

        # Distribution of v2
        v2_dist = defaultdict(int)
        for v in v2_vals:
            v2_dist[v] += 1

        print(f"\n  Degree {deg}: {len(entries)} coefficients")
        print(f"    v₂ range: [{min_v2}, {max_v2}]")
        print(f"    v₂ distribution: {dict(sorted(v2_dist.items()))}")

        # All divided by 2^min_v2
        reduced = sorted(set(c // (2**min_v2) for c in vals))
        print(f"    After dividing by 2^{min_v2}: odd parts = {reduced[:20]}")

        # Check: does v₂(c_S) depend on any simple function of S?
        if deg == 2 and len(entries) <= 100:
            print(f"    Detailed v₂ by tile pair:")
            for S, c in sorted(entries, key=lambda x: v2(abs(x[1]))):
                bits = [i for i in range(m) if S & (1 << i)]
                tile_set = [tiles[i] for i in bits]
                skips = [tiles[i][0] - tiles[i][1] for i in bits]
                print(f"      {tile_set} skips={skips}: c={c}, v₂={v2(abs(c))}")

def analyze_gradient(coeffs, m, tiles, n):
    """Analyze the gradient ∂H/∂t_k."""
    print(f"\n{'='*60}")
    print(f"GRADIENT STRUCTURE ∂H/∂t_k (n={n})")
    print(f"{'='*60}")

    for k in range(m):
        # ∂H/∂t_k = Σ_{S∋k} c_S · Π_{i∈S\{k}} t_i
        # This is a multilinear polynomial in the OTHER m-1 variables
        grad_k_coeffs = {}
        for S, c in coeffs.items():
            if S & (1 << k):
                S_minus_k = S ^ (1 << k)
                grad_k_coeffs[S_minus_k] = c

        if not grad_k_coeffs:
            continue

        # Evaluate at t=0 (all forward = transitive)
        val_at_0 = grad_k_coeffs.get(0, 0)

        # Evaluate at t=1 (all backward)
        val_at_1 = sum(grad_k_coeffs.values())

        # Count positive/negative terms
        pos = sum(1 for c in grad_k_coeffs.values() if c > 0)
        neg = sum(1 for c in grad_k_coeffs.values() if c < 0)

        # Sum of all coefficients
        total = sum(grad_k_coeffs.values())

        tile = tiles[k]
        skip = tile[0] - tile[1]

        # Max degree of gradient
        max_deg = max(bin(S).count('1') for S in grad_k_coeffs.keys()) if grad_k_coeffs else 0

        print(f"  tile {tile} skip={skip}: ∂H/∂t_k has {len(grad_k_coeffs)} terms, "
              f"deg={max_deg}, at t=0: {val_at_0}, at t=1: {val_at_1}, "
              f"sum={total}, (+{pos}/-{neg})")

def analyze_interaction_graph(coeffs, m, tiles, n):
    """Build the tile interaction graph from quadratic coefficients."""
    print(f"\n{'='*60}")
    print(f"TILE INTERACTION GRAPH (n={n})")
    print(f"{'='*60}")

    # Vertices = tiles. Edge (i,j) iff c_{ij} ≠ 0.
    edges = []
    pos_edges = []
    neg_edges = []
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        i, j = bits
        edges.append((i, j, c))
        if c > 0:
            pos_edges.append((i, j))
        else:
            neg_edges.append((i, j))

    print(f"  Vertices: {m} tiles")
    print(f"  Edges: {len(edges)} ({len(pos_edges)} positive, {len(neg_edges)} negative)")

    total_possible = m * (m - 1) // 2
    print(f"  Density: {len(edges)}/{total_possible} = {len(edges)/total_possible:.3f}")

    # Degree sequence
    deg_seq = [0] * m
    for i, j, _ in edges:
        deg_seq[i] += 1
        deg_seq[j] += 1

    print(f"  Degree sequence: {sorted(deg_seq, reverse=True)}")
    print(f"  Min degree: {min(deg_seq)}, Max degree: {max(deg_seq)}")

    # Chromatic structure: is the negative subgraph special?
    print(f"\n  Negative edges (same-vertex pairs):")
    for i, j in neg_edges:
        t1, t2 = tiles[i], tiles[j]
        shared = []
        if t1[0] == t2[0]: shared.append(f"upper={t1[0]}")
        if t1[1] == t2[1]: shared.append(f"lower={t1[1]}")
        if t1[0] == t2[1]: shared.append(f"x1=y2")
        if t1[1] == t2[0]: shared.append(f"y1=x2")
        if len(neg_edges) <= 40:
            print(f"    {t1}-{t2}: {','.join(shared) if shared else 'DISJOINT!'}")

    # Check: are ALL negative edges same-vertex?
    neg_same = 0
    neg_cross = 0
    for i, j in neg_edges:
        t1, t2 = tiles[i], tiles[j]
        verts1 = {t1[0], t1[1]}
        verts2 = {t2[0], t2[1]}
        if verts1 & verts2:
            neg_same += 1
        else:
            neg_cross += 1
    print(f"\n  Negative edges: {neg_same} same-vertex, {neg_cross} cross-vertex")

    # Check: are ALL positive edges cross-end?
    pos_same = 0
    pos_cross_end = 0
    pos_disjoint = 0
    for i, j in pos_edges:
        t1, t2 = tiles[i], tiles[j]
        verts1 = {t1[0], t1[1]}
        verts2 = {t2[0], t2[1]}
        if not (verts1 & verts2):
            pos_disjoint += 1
        elif t1[0] == t2[1] or t1[1] == t2[0]:
            pos_cross_end += 1
        else:
            pos_same += 1
    print(f"  Positive edges: {pos_same} same-end, {pos_cross_end} cross-end, {pos_disjoint} disjoint")

    # Graph Laplacian
    A = np.zeros((m, m))
    for i, j, c in edges:
        A[i][j] = c
        A[j][i] = c
    D = np.diag(A.sum(axis=1))
    L = D - A
    lap_eigvals = sorted(np.linalg.eigvalsh(L))
    print(f"\n  Laplacian eigenvalues: {[f'{e:.4f}' for e in lap_eigvals]}")
    print(f"  Algebraic connectivity (λ₂): {lap_eigvals[1]:.6f}")

def coefficient_statistics(coeffs, m, tiles, n):
    """Comprehensive coefficient statistics."""
    print(f"\n{'='*60}")
    print(f"COEFFICIENT STATISTICS (n={n})")
    print(f"{'='*60}")

    all_vals = [c for c in coeffs.values() if c != 0]
    all_vals_nonconst = [c for S, c in coeffs.items() if c != 0 and S != 0]

    print(f"  Total nonzero: {len(all_vals)}")
    print(f"  Sum of all coefficients: {sum(all_vals)}")
    print(f"  Sum of |coefficients|: {sum(abs(c) for c in all_vals)}")
    print(f"  Sum of positive: {sum(c for c in all_vals if c > 0)}")
    print(f"  Sum of negative: {sum(c for c in all_vals if c < 0)}")

    # The sum of all coefficients = H(1,1,...,1) = H(all-backward tiling)
    # This is the anti-transitive tournament's H
    print(f"  H(all 1s) = sum(c_S) = {sum(all_vals)}")

    # Alternating sum = ?
    alt_sum = sum(c * (-1)**bin(S).count('1') for S, c in coeffs.items() if c != 0)
    print(f"  Alternating sum Σ (-1)^|S| c_S = {alt_sum}")
    # This should be H(all 0s) = H(transitive) = 1... wait, that's c_∅ = 1.
    # Actually H(0,...,0) = Σ_S c_S · 0^... only S=∅ contributes = c_∅ = 1. Good.
    # The alternating sum is H(-1,-1,...,-1) which doesn't make sense for binary.
    # Actually H evaluated at t_i = -1 for all i... hmm, not defined for binary.
    # Let me compute H at t = (1/2, 1/2, ..., 1/2) instead — the "center" of the cube.

    center_val = sum(c / (2**bin(S).count('1')) for S, c in coeffs.items())
    print(f"  H(1/2, ..., 1/2) = {center_val:.6f}")
    print(f"  = average of H over all tilings? No, that's sum H / 2^m = ", end="")
    # Actually H at (1/2,...,1/2) IS the average of H over all tilings (by multilinearity)

    # Moments
    by_degree = defaultdict(list)
    for S, c in coeffs.items():
        if c != 0:
            by_degree[bin(S).count('1')].append(c)

    print(f"\n  By-degree sums:")
    for d in sorted(by_degree.keys()):
        vals = by_degree[d]
        s = sum(vals)
        sa = sum(abs(c) for c in vals)
        print(f"    d={d}: count={len(vals)}, sum={s}, |sum|={sa}, "
              f"mean={np.mean(vals):.3f}")

def main():
    for n in [5, 6, 7]:
        print(f"\n\n{'#'*70}")
        print(f"# n = {n}")
        print(f"{'#'*70}")

        print(f"Computing H for all tilings...")
        H, tiles, m = compute_all_H(n)
        print(f"Done. m={m}, H range=[{H.min()}, {H.max()}]")

        # Multilinear coefficients
        coeffs = mobius_inversion(H, m)

        # 1. Walsh spectrum (only for n≤6 due to cost)
        if n <= 6:
            walsh = walsh_transform(H, m)
            analyze_walsh_spectrum(walsh, m, tiles)

        # 2. Quadratic form matrix
        Q, eigvals = analyze_quadratic_form(coeffs, m, tiles, n)

        # 3. Support topology
        analyze_support_topology(coeffs, m, tiles, n)

        # 4. 2-adic structure
        analyze_2adic(coeffs, m, tiles, n)

        # 5. Gradient
        analyze_gradient(coeffs, m, tiles, n)

        # 6. Interaction graph
        analyze_interaction_graph(coeffs, m, tiles, n)

        # 7. Statistics
        coefficient_statistics(coeffs, m, tiles, n)

    # =============================================
    # SPECIAL: n=7 Walsh spectrum (integer version)
    # =============================================
    print(f"\n\n{'#'*70}")
    print(f"# n=7 WALSH SPECTRUM (integer-scaled)")
    print(f"{'#'*70}")

    n = 7
    H, tiles, m = compute_all_H(n)
    N = 1 << m

    # Instead of dividing by 2^m, just compute the unnormalized Walsh coefficients
    # ĥ_S * 2^m = Σ_t H(t) (-1)^{popcount(S&t)}
    # This is always an integer.
    print(f"Computing Walsh transform for n=7 (m={m}, N={N})...")

    # Walsh transform via fast method (in-place butterfly)
    # Start with H values, apply Hadamard transform
    f = np.array(H, dtype=np.int64)
    for i in range(m):
        step = 1 << i
        for j in range(N):
            if j & step:
                a = f[j ^ step]
                b = f[j]
                f[j ^ step] = a + b
                f[j] = a - b
    # Now f[S] = Σ_t H(t) (-1)^{popcount(S&t)} = ĥ_S * 2^m

    # Analyze
    by_order = defaultdict(list)
    for S in range(N):
        if f[S] != 0:
            order = bin(S).count('1')
            by_order[order].append((S, int(f[S])))

    print(f"\nNonzero unnormalized Walsh coefficients by order:")
    for order in sorted(by_order.keys()):
        entries = by_order[order]
        vals = [v for _, v in entries]
        print(f"  Order {order}: {len(entries)} nonzero")
        distinct = sorted(set(vals))
        print(f"    Distinct values: {distinct[:15]}{'...' if len(distinct) > 15 else ''}")

        # v₂ of the values (to determine actual Walsh coefficient v₂)
        v2_vals = [v2(abs(v)) for v in vals if v != 0]
        if v2_vals:
            min_v2 = min(v2_vals)
            print(f"    Min v₂ of unnorm: {min_v2} (actual Walsh v₂ = {min_v2} - {m})")
            # So actual Walsh coeff has v₂ = min_v2 - m

    # What fraction of Walsh spectrum is zero?
    nonzero_count = sum(len(v) for v in by_order.values())
    print(f"\n  Sparsity: {nonzero_count}/{N} nonzero = {nonzero_count/N:.4f}")

    # Parseval: Σ (ĥ_S * 2^m)² = 2^m * Σ H(t)²
    parseval_lhs = sum(int(f[S])**2 for S in range(N))
    parseval_rhs = N * sum(int(H[t])**2 for t in range(N))
    print(f"  Parseval: LHS={parseval_lhs}, RHS={parseval_rhs}, match={parseval_lhs == parseval_rhs}")

    # Check: is the Walsh spectrum related to the multilinear coefficients?
    # Relationship: c_S = sum_{T supset S} (-1)^{|T|-|S|} h_T ... verify.
    # multilinear coefficients in {0,1} basis vs Walsh in {+/-1} basis.
    # The conversion is: h_S = sum_{T supset S} c_T / (-2)^|T|
    # Or: c_S = sum_{T supset S} (-2)^|T| h_T ... verify computationally.

    coeffs = mobius_inversion(H, m)

    print(f"\n  Checking multilinear ↔ Walsh conversion:")
    # c_S = sum_{T supset S} (-2)^{|T|-|S|} h_T    where h_T = f[T] / 2^m
    # So c_S = sum_{T supset S} (-2)^{|T|-|S|} f[T] / 2^m
    # = (1/2^m) sum_{T supset S} (-1)^{|T|-|S|} * 2^{|T|-|S|} f[T]

    mismatches = 0
    for S in range(min(N, 100)):  # Check first 100
        if S not in coeffs and not any(S == k for k in coeffs):
            c_S_expected = 0
        else:
            c_S_expected = coeffs.get(S, 0)

        # Compute from Walsh
        c_S_from_walsh = 0
        # Iterate over all supersets T of S
        T = S
        complement = ((1 << m) - 1) ^ S
        sub = complement
        while True:
            T_full = S | sub
            diff = bin(T_full ^ S).count('1')  # |T| - |S|
            c_S_from_walsh += ((-2)**diff) * f[T_full]
            if sub == 0:
                break
            sub = (sub - 1) & complement
        c_S_from_walsh //= N  # divide by 2^m

        if c_S_from_walsh != c_S_expected:
            mismatches += 1
            if mismatches <= 5:
                print(f"    MISMATCH at S={S}: expected {c_S_expected}, got {c_S_from_walsh}")

    if mismatches == 0:
        print(f"    ✓ PERFECT MATCH for first {min(N, 100)} coefficients")
    else:
        print(f"    {mismatches} mismatches — conversion formula wrong, trying alternatives...")

        # Try: c_S = Σ_{T⊇S} (-1)^{|T\S|} 2^{|T|} ĥ_T ... different normalization
        mismatches2 = 0
        for S in range(min(N, 20)):
            c_S_expected = coeffs.get(S, 0)
            c_S_v2 = 0
            complement = ((1 << m) - 1) ^ S
            sub = complement
            while True:
                T_full = S | sub
                diff = bin(T_full ^ S).count('1')
                c_S_v2 += ((-1)**diff) * f[T_full]
                if sub == 0:
                    break
                sub = (sub - 1) & complement
            c_S_v2 //= N

            if c_S_v2 != c_S_expected:
                mismatches2 += 1

        if mismatches2 == 0:
            print(f"    Alt formula works: c_S = (1/2^m) sum_{{T supset S}} (-1)^|T\\S| f[T]")

if __name__ == "__main__":
    main()
