"""
spectral_H_formula.py — Search for the spectral formula H = f({lambda_k})
for circulant tournaments.

CORE IDEA:
  H(T) = I(Omega(T), 2) = sum of Hamiltonian path counts.
  For circulant tournaments on Z_p, the adjacency matrix A has eigenvalues
    lambda_k = sum_{s in S} omega^{ks},  omega = e^{2*pi*i/p}
  with Re(lambda_k) = -1/2 and |lambda_k|^2 = 1/4 + y_k^2.

  Known spectral quantities:
    - det(A) = product of lambda_k  (trivially)
    - tr(A^m) = sum lambda_k^m = m * c_m  (directed m-cycle count)
    - permanent(A) = H(T)  ... but perm is NOT a spectral invariant in general!

  HOWEVER: for CIRCULANT matrices, the permanent IS a function of eigenvalues!
  Minc's theorem: perm(circ(c_0,...,c_{n-1})) = product_{k=0}^{n-1} (sum_j c_j omega^{jk})
                 ... NO, that's the determinant!

  Actually: for circulant matrices, the permanent does NOT simplify to a product.
  BUT: H(T) = I(Omega(T), 2), and the independence polynomial can be computed
  from the clique polynomial of the complement, which for structured graphs
  may have spectral formulas.

  EMPIRICAL APPROACH:
  Compute (H, {y_k^2}) for ALL circulant tournaments at primes p=3,5,7,11,13,17
  and search for the functional relationship.

  KEY INSIGHT FROM LADDER ANALYSIS:
  The y_k^2 values live on a simplex: sum y_k^2 = p(p-1)/8 (from sum |lambda_k|^2).
  The question is: on this simplex, where is H maximized?

  For p=7 (mod4=3): center of simplex (flat, y_k^2 = p/4 = 7/4 for all k)
  For p=13 (mod4=1): NOT the center; instead the extremal Dirichlet kernel point.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import time
import cmath
import math
from itertools import combinations
from collections import defaultdict

sys.path.insert(0, '04-computation')

# ======================================================================
# CORE FUNCTIONS
# ======================================================================

def circulant_eigenvalues(n, S):
    """Eigenvalues lambda_k for k=0,...,n-1."""
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def spectral_y_squared(n, S):
    """Return list of y_k^2 for k=1,...,(n-1)/2.
    lambda_k = -1/2 + i*y_k, so y_k = Im(lambda_k).
    y_{n-k} = -y_k, so we only need half."""
    eigs = circulant_eigenvalues(n, S)
    m = (n - 1) // 2
    return [eigs[k].imag ** 2 for k in range(1, m + 1)]


def all_circulant_tournaments(n):
    """Generate all circulant tournament connection sets on Z_n (odd n).
    S must have |S| = (n-1)/2 and S, -S partition {1,...,n-1}."""
    k = (n - 1) // 2
    elts = list(range(1, n))
    # Pair each element with its complement
    pairs = []
    used = set()
    for a in elts:
        if a not in used:
            b = n - a
            if a == b:
                return []  # even n / self-paired
            pairs.append((a, b))
            used.add(a)
            used.add(b)
    # Choose one from each pair
    results = []
    for bits in range(2 ** len(pairs)):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits >> i) & 1 else b)
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    """Hamiltonian path count for circulant tournament on Z_n using Held-Karp DP."""
    # Build adjacency: i->j iff (j-i) mod n in S
    adj = [[False] * n for _ in range(n)]
    S_set = set(S)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i][j] = True

    # Held-Karp DP
    full_mask = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]

    return sum(dp[full_mask][v] for v in range(n))


def qr_set(p):
    """Quadratic residues mod p."""
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})


def cyclic_interval(n):
    """S = {ceil(n/2), ..., n-1}."""
    k = (n - 1) // 2
    return tuple(range(n - k, n))


# ======================================================================
# ANALYSIS: POWER SUM DETERMINATION
# ======================================================================

def power_sums(y2_list, max_pow=6):
    """Compute S_j = sum y_k^{2j} for j=1,...,max_pow."""
    return [sum(y ** j for y in y2_list) for j in range(1, max_pow + 1)]


def analyze_prime(p, verbose=True):
    """Full spectral analysis at prime p."""
    if verbose:
        print(f"\n{'=' * 70}")
        print(f"PRIME p = {p}  (p mod 4 = {p % 4}, p mod 8 = {p % 8})")
        print(f"{'=' * 70}")

    m = (p - 1) // 2
    all_S = all_circulant_tournaments(p)
    n_tours = len(all_S)

    if verbose:
        print(f"  {n_tours} circulant tournaments, {m} spectral parameters")

    # Compute H and spectral data for all
    t0 = time.time()
    data = []
    for S in all_S:
        if p <= 13:
            H = ham_count_dp(p, S)
        else:
            H = ham_count_dp(p, S)
        y2 = spectral_y_squared(p, S)
        ps = power_sums(y2, max_pow=m)
        data.append({
            'S': S,
            'H': H,
            'y2': y2,
            'power_sums': ps,
        })
    elapsed = time.time() - t0
    if verbose:
        print(f"  Computed in {elapsed:.1f}s")

    # Sort by H descending
    data.sort(key=lambda d: d['H'], reverse=True)

    # Distinct H values
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    if verbose:
        print(f"  {len(H_vals)} distinct H values: {H_vals}")
        print(f"  H_max = {H_vals[0]}, H_min = {H_vals[-1]}")
        print(f"  Max H tournament: S = {data[0]['S']}")

    # Check QR
    qr = tuple(sorted(qr_set(p)))
    minus1 = p - 1
    qr_valid = minus1 not in qr
    if verbose:
        print(f"  QR_{p} = {qr}, valid = {qr_valid}")
    qr_H = None
    if qr_valid:
        for d in data:
            if d['S'] == qr:
                qr_H = d['H']
                break
        if verbose:
            print(f"  H(QR) = {qr_H}, is max: {qr_H == H_vals[0]}")

    # Cyclic interval
    ci = cyclic_interval(p)
    ci_H = None
    for d in data:
        if d['S'] == ci:
            ci_H = d['H']
            break
    if verbose:
        print(f"  Cyclic interval S = {ci}, H = {ci_H}, is max: {ci_H == H_vals[0]}")

    # Spectral center distance
    center_y2 = p / 4.0
    for d in data:
        d['dist_center'] = sum((y - center_y2) ** 2 for y in d['y2'])
        d['spread'] = max(d['y2']) - min(d['y2']) if d['y2'] else 0
        d['sum_y4'] = d['power_sums'][1] if len(d['power_sums']) > 1 else d['power_sums'][0] ** 2  # sum y_k^4

    # Print table
    if verbose:
        print(f"\n  {'rank':>4} {'H':>12} {'spread':>8} {'dist_ctr':>10} {'sum_y4':>12} {'S'}")
        print(f"  {'-' * 4} {'-' * 12} {'-' * 8} {'-' * 10} {'-' * 12} {'-' * 20}")
        for i, d in enumerate(data):
            marker = ""
            if d['S'] == qr and qr_valid:
                marker = " <-- QR"
            elif d['S'] == ci:
                marker = " <-- CI"
            if i < 10 or i >= len(data) - 3 or d['S'] == qr or d['S'] == ci:
                print(f"  {i + 1:>4} {d['H']:>12} {d['spread']:>8.4f} "
                      f"{d['dist_center']:>10.4f} {d['sum_y4']:>12.4f} {list(d['S'])}{marker}")
            elif i == 10:
                print(f"  {'...':>4}")

    # KEY TEST: Is H monotone in distance from center?
    # Group by H value and check
    if verbose:
        print(f"\n  SPECTRAL MONOTONICITY TEST:")
        H_to_dist = defaultdict(list)
        for d in data:
            H_to_dist[d['H']].append(d['dist_center'])
        # Check: is H decreasing in distance from center (p=3 mod 4)
        # or increasing (p=1 mod 4)?
        prev_H, prev_dist = None, None
        monotone_inc = True
        monotone_dec = True
        for H in H_vals:
            avg_dist = sum(H_to_dist[H]) / len(H_to_dist[H])
            if prev_H is not None:
                if avg_dist <= prev_dist:
                    monotone_inc = False
                if avg_dist >= prev_dist:
                    monotone_dec = False
            prev_H, prev_dist = H, avg_dist
            print(f"    H = {H:>12}: avg dist from center = {avg_dist:.4f}, "
                  f"count = {len(H_to_dist[H])}")
        if p % 4 == 3:
            print(f"    p = 3 mod 4: H decreasing in dist? {monotone_dec}")
        else:
            print(f"    p = 1 mod 4: H increasing in dist? {monotone_inc}")

    # KEY TEST: Newton's identities — can we reconstruct H from power sums?
    if verbose:
        print(f"\n  POWER SUM DETERMINATION:")
        # Group data by power sum signature
        ps_to_H = defaultdict(set)
        for d in data:
            ps_key = tuple(round(x, 6) for x in d['power_sums'][:m])
            ps_to_H[ps_key].add(d['H'])
        all_determined = all(len(v) == 1 for v in ps_to_H.values())
        print(f"    {len(ps_to_H)} distinct power sum signatures")
        print(f"    {len(H_vals)} distinct H values")
        print(f"    H determined by power sums? {all_determined}")

    # NEW: Multivariate regression H = f(S2, S3, ..., Sm)
    if verbose and len(H_vals) > 1 and m >= 2:
        print(f"\n  POLYNOMIAL FIT (H as function of power sums S_j = sum y_k^{{2j}}):")
        print(f"    Note: S_1 = sum y_k^2 = p(p-1)/8 = {p * (p - 1) / 8} is CONSTANT")

        # Collect non-constant power sums: S_2, S_3, ..., S_m
        n_pts = len(data)
        H_list = [d['H'] for d in data]

        # Try linear in S2 alone
        if m >= 2:
            S2_vals = [d['power_sums'][1] for d in data]
            sx = sum(S2_vals)
            sy = sum(H_list)
            sxx = sum(x ** 2 for x in S2_vals)
            sxy = sum(x * y for x, y in zip(S2_vals, H_list))
            denom = n_pts * sxx - sx ** 2
            if abs(denom) > 1e-10:
                b = (n_pts * sxy - sx * sy) / denom
                a = (sy - b * sx) / n_pts
                residuals = [H - (a + b * s2) for H, s2 in zip(H_list, S2_vals)]
                max_res = max(abs(r) for r in residuals)
                ss_res = sum(r ** 2 for r in residuals)
                ss_tot = sum((H - sy / n_pts) ** 2 for H in H_list)
                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1
                print(f"    Linear in S2: H = {a:.4f} + ({b:.6f}) * S2")
                print(f"    R^2 = {r2:.8f}, max |residual| = {max_res:.4f}")
                if max_res < 0.01:
                    print(f"    *** EXACT LINEAR FORMULA FOUND ***")
                    # Check if coefficients are rational
                    # a should be (H_flat - b * S2_flat)
                    print(f"    Interpretation: H is AFFINE in sum(y_k^4)")

        # Try multiple regression: H = a + b*S2 + c*S3
        if m >= 3 and len(H_vals) >= 3:
            S2_vals = [d['power_sums'][1] for d in data]
            S3_vals = [d['power_sums'][2] for d in data]
            # Solve via Gauss elimination
            # [n  sx2  sx3] [a]   [sy]
            # [sx2 s22 s23] [b] = [sy2]
            # [sx3 s23 s33] [c]   [sy3]
            A_mat = [
                [n_pts, sum(S2_vals), sum(S3_vals)],
                [sum(S2_vals), sum(x ** 2 for x in S2_vals), sum(x * y for x, y in zip(S2_vals, S3_vals))],
                [sum(S3_vals), sum(x * y for x, y in zip(S2_vals, S3_vals)), sum(x ** 2 for x in S3_vals)]
            ]
            rhs = [
                sum(H_list),
                sum(h * s2 for h, s2 in zip(H_list, S2_vals)),
                sum(h * s3 for h, s3 in zip(H_list, S3_vals))
            ]
            # Gauss elimination (3x3)
            for col in range(3):
                # Pivot
                max_row = max(range(col, 3), key=lambda r: abs(A_mat[r][col]))
                A_mat[col], A_mat[max_row] = A_mat[max_row], A_mat[col]
                rhs[col], rhs[max_row] = rhs[max_row], rhs[col]
                if abs(A_mat[col][col]) < 1e-15:
                    continue
                for row in range(col + 1, 3):
                    factor = A_mat[row][col] / A_mat[col][col]
                    for j in range(col, 3):
                        A_mat[row][j] -= factor * A_mat[col][j]
                    rhs[row] -= factor * rhs[col]
            # Back substitution
            coeffs = [0, 0, 0]
            for i in range(2, -1, -1):
                if abs(A_mat[i][i]) > 1e-15:
                    coeffs[i] = (rhs[i] - sum(A_mat[i][j] * coeffs[j] for j in range(i + 1, 3))) / A_mat[i][i]
            a, b, c = coeffs
            residuals = [H - (a + b * s2 + c * s3) for H, s2, s3 in zip(H_list, S2_vals, S3_vals)]
            max_res = max(abs(r) for r in residuals)
            ss_res = sum(r ** 2 for r in residuals)
            ss_tot = sum((H - sum(H_list) / n_pts) ** 2 for H in H_list)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1
            print(f"\n    Bilinear in S2, S3: H = {a:.4f} + ({b:.6f})*S2 + ({c:.6f})*S3")
            print(f"    R^2 = {r2:.8f}, max |residual| = {max_res:.4f}")

        # Exact integer test
        print(f"\n    Exact data points:")
        for H_val in H_vals[:8]:
            matches = [d for d in data if d['H'] == H_val]
            ps_example = matches[0]['power_sums']
            count = len(matches)
            print(f"    H={H_val:>12} (x{count:>3}): S_j = {[round(x, 4) for x in ps_example[:min(5, m)]]}")

    return data, H_vals


# ======================================================================
# SECTION 2: VARIATIONAL PRINCIPLE DEEP TEST
# ======================================================================

def variational_test(p, data):
    """Test the spectral variational principle at prime p."""
    print(f"\n{'=' * 70}")
    print(f"VARIATIONAL PRINCIPLE TEST — p = {p}")
    print(f"{'=' * 70}")

    m = (p - 1) // 2
    center_y2 = p / 4.0
    total_y2 = p * (p - 1) / 8.0  # sum of all y_k^2

    print(f"  Spectral simplex: sum y_k^2 = {total_y2:.1f}")
    print(f"  Center point: y_k^2 = {center_y2:.4f} for all k")
    print(f"  m = {m} spectral parameters")

    # Is center achievable?
    # Center achievable iff Paley/QR is valid (p = 3 mod 4)
    qr = tuple(sorted(qr_set(p)))
    qr_valid = (p - 1) not in qr
    center_achievable = qr_valid and p % 4 == 3

    if center_achievable:
        # Check Paley is actually flat
        y2_paley = spectral_y_squared(p, qr)
        is_flat = all(abs(y - center_y2) < 1e-8 for y in y2_paley)
        print(f"  Center achievable: {center_achievable} (Paley is flat: {is_flat})")
    else:
        print(f"  Center achievable: {center_achievable}")

    # For each tournament, compute distance^2 from center
    for d in data:
        d['dist2'] = sum((y - center_y2) ** 2 for y in d['y2'])
        # Also compute the "4th moment" sum y_k^4 = concentration measure
        d['conc'] = sum(y ** 2 for y in d['y2'])
        # By Cauchy-Schwarz: sum y_k^4 >= (sum y_k^2)^2 / m
        # Equality iff flat. So conc >= total_y2^2 / m.

    min_conc = total_y2 ** 2 / m
    print(f"  Minimum sum y^4 (flat) = {min_conc:.4f}")
    print(f"  Achieved by maximizer: {data[0]['conc']:.4f}")

    # Correlation between H and distance from center
    H_vals = [d['H'] for d in data]
    dist_vals = [d['dist2'] for d in data]
    conc_vals = [d['conc'] for d in data]

    n = len(data)
    mean_H = sum(H_vals) / n
    mean_d = sum(dist_vals) / n
    mean_c = sum(conc_vals) / n

    cov_Hd = sum((H - mean_H) * (d - mean_d) for H, d in zip(H_vals, dist_vals)) / n
    var_H = sum((H - mean_H) ** 2 for H in H_vals) / n
    var_d = sum((d - mean_d) ** 2 for d in dist_vals) / n

    cov_Hc = sum((H - mean_H) * (c - mean_c) for H, c in zip(H_vals, conc_vals)) / n
    var_c = sum((c - mean_c) ** 2 for c in conc_vals) / n

    if var_H > 0 and var_d > 0:
        corr_Hd = cov_Hd / (var_H ** 0.5 * var_d ** 0.5)
    else:
        corr_Hd = 0

    if var_H > 0 and var_c > 0:
        corr_Hc = cov_Hc / (var_H ** 0.5 * var_c ** 0.5)
    else:
        corr_Hc = 0

    print(f"\n  corr(H, dist_from_center) = {corr_Hd:.6f}")
    print(f"  corr(H, sum_y^4) = {corr_Hc:.6f}")

    if p % 4 == 3:
        print(f"  p = 3 mod 4: expect negative correlation (flat=max)")
        print(f"  Maximizer at center? {data[0]['dist2'] < 1e-8}")
    else:
        print(f"  p = 1 mod 4: expect positive correlation (concentrated=max)")
        print(f"  Maximizer farthest from center? ", end="")
        max_dist = max(d['dist2'] for d in data)
        print(f"{abs(data[0]['dist2'] - max_dist) < 1e-6}")

    return corr_Hd, corr_Hc


# ======================================================================
# SECTION 3: GROUP-THEORETIC H FORMULA
# ======================================================================

def group_theoretic_analysis(p, data):
    """Analyze H through the lens of D_{2p} representation theory."""
    print(f"\n{'=' * 70}")
    print(f"D_{{{2 * p}}} REPRESENTATION THEORY — p = {p}")
    print(f"{'=' * 70}")

    m = (p - 1) // 2

    # D_{2p} irreducible representations (p odd prime):
    # - Two 1-dimensional: trivial (rho_0) and sign (rho_s)
    # - (p-1)/2 two-dimensional: rho_k for k=1,...,m
    #   where rho_k(r) = diag(omega^k, omega^{-k}), rho_k(s) = [[0,1],[1,0]]

    print(f"  D_{{{2 * p}}} has {2 + m} irreps:")
    print(f"    2 one-dimensional (trivial, sign)")
    print(f"    {m} two-dimensional (rho_1, ..., rho_{m})")

    # For circulant tournament adjacency matrix A:
    # DFT decomposes C^p = bigoplus V_k
    # V_0 = span(1,...,1) with eigenvalue m (= lambda_0 = |S|)
    # V_k = V_{p-k} for k=1,...,m, each 2-dimensional
    # lambda_k and lambda_{p-k} = conj(lambda_k) are eigenvalues in this block

    # The key: H = perm(A). For circulant matrix,
    # perm factors somehow through the DFT...

    # Actually, let's think about it differently.
    # H(T) = I(Omega(T), 2) where Omega(T) is the conflict graph.
    # For circulant T, Omega(T) is also vertex-transitive.
    # I(G, x) = sum_k alpha_k x^k where alpha_k = # independent sets of size k.

    # For vertex-transitive G on N vertices with I(G,x) = sum a_k x^k:
    # a_k = (N/k) * # ind sets of size k containing vertex 0

    # For circulant T on Z_p: Omega is vertex-transitive under Z_p action
    # (since Omega is defined by distances, i.e. it's also a circulant-type graph)

    # The independence polynomial I(Omega, x) is related to the
    # clique polynomial of the complement bar{Omega}.

    # APPROACH: compute the alpha_k (independence numbers) and see
    # how they relate to eigenvalues.

    print(f"\n  INDEPENDENCE POLYNOMIAL DECOMPOSITION:")
    for d in data[:5]:
        S = d['S']
        H = d['H']
        # Compute I(Omega(T), x) = 1 + a_1*x + a_2*x^2 + ...
        # where a_k = # independent sets of size k in Omega(T)
        # Omega(T) has vertex set = {(i,j,k): i->j->k, k->i}
        # Two vertices conflict if they share a vertex

        # For small p, compute Omega explicitly
        if p <= 13:
            S_set = set(S)
            # Build adjacency
            adj = {}
            for i in range(p):
                for j in range(p):
                    if i != j:
                        adj[(i, j)] = (j - i) % p in S_set

            # Find all directed 3-cycles
            cycles_3 = []
            for i in range(p):
                for j in range(p):
                    if i == j:
                        continue
                    if not adj[(i, j)]:
                        continue
                    for k in range(p):
                        if k == i or k == j:
                            continue
                        if adj[(j, k)] and adj[(k, i)]:
                            cycles_3.append((i, j, k))

            # Conflict graph: two cycles conflict if they share a vertex
            c3 = len(cycles_3) // 3  # each cycle counted 3 times
            cycle_sets = []
            seen = set()
            for c in cycles_3:
                key = tuple(sorted(c))
                if key not in seen:
                    seen.add(key)
                    cycle_sets.append(set(c))

            # Build Omega adjacency
            omega_n = len(cycle_sets)
            omega_adj = [[False] * omega_n for _ in range(omega_n)]
            for a in range(omega_n):
                for b in range(a + 1, omega_n):
                    if cycle_sets[a] & cycle_sets[b]:  # share a vertex
                        omega_adj[a][b] = True
                        omega_adj[b][a] = True

            # Count independent sets by size
            max_ind = 0
            ind_counts = [0] * (omega_n + 1)
            ind_counts[0] = 1

            # For small Omega, enumerate via bitmask
            if omega_n <= 20:
                for mask in range(1, 1 << omega_n):
                    verts = [i for i in range(omega_n) if mask & (1 << i)]
                    is_ind = True
                    for a in range(len(verts)):
                        for b in range(a + 1, len(verts)):
                            if omega_adj[verts[a]][verts[b]]:
                                is_ind = False
                                break
                        if not is_ind:
                            break
                    if is_ind:
                        ind_counts[len(verts)] += 1

                # Verify I(Omega, 2) = H
                I_at_2 = sum(ind_counts[k] * (2 ** k) for k in range(omega_n + 1))
                ip = [ind_counts[k] for k in range(omega_n + 1) if ind_counts[k] > 0]
                print(f"    S={list(S)}: H={H}, I(Omega,2)={I_at_2}, "
                      f"IP={ip}, |Omega|={omega_n}")
                if I_at_2 != H:
                    print(f"    *** MISMATCH ***")

    return


# ======================================================================
# SECTION 4: TRACE FORMULA CONNECTION
# ======================================================================

def trace_formula_analysis(p, data):
    """Explore H in terms of traces (cycle counts) and eigenvalues."""
    print(f"\n{'=' * 70}")
    print(f"TRACE FORMULA — CYCLE COUNTS vs EIGENVALUES — p = {p}")
    print(f"{'=' * 70}")

    # For circulant A with eigenvalues lambda_k:
    # tr(A^m) = sum_k lambda_k^m = p * c_m  (for m odd, directed m-cycle count times p)
    # Actually tr(A^m) = sum of (i1->i2->...->im->i1) walks, each CLOSED walk counted once
    # For m=1: tr(A) = 0 (no self-loops)
    # For m=2: tr(A^2) = # reciprocal edges = 0 (tournament)
    # For m=3: tr(A^3) = 3 * c_3 (3-cycles counted with rotation)
    # ...

    # The cycle counts c_k are Newton power sums in the eigenvalues.
    # By Newton's identities, they determine the characteristic polynomial.
    # But H = permanent, not determinant!

    # KEY OBSERVATION: For CIRCULANT matrices, trace formulas give us
    # ALL cycle counts c_k for k=3,5,7,...,p (odd only; even are 0 for tournaments).
    # Then OCF gives: H = sum_{S subset of cycles, disjoint} 2^|S|
    # i.e., H = I(Omega, 2) where Omega is the conflict graph of cycles.

    # So: eigenvalues -> {c_k} -> Omega -> I(Omega, 2) = H.
    # But eigenvalues DON'T directly determine Omega's structure!
    # They determine the cycle COUNTS, not the conflict structure.

    # HOWEVER: for circulant tournaments, the conflict graph is also structured.
    # Two cycles conflict iff they share a vertex.
    # Since the tournament is circulant, Omega inherits Z_p symmetry.

    # Let's verify: are the cycle counts sufficient to determine H?
    # At p=7: 2 distinct H values, but c_3 is the same for all regular tournaments!
    # So c_k alone do NOT determine H for circulant tournaments.

    print(f"\n  Cycle counts for each circulant tournament:")
    for d in data[:min(10, len(data))]:
        S = d['S']
        H = d['H']
        eigs = circulant_eigenvalues(p, S)

        # Compute cycle counts from traces
        cycles = {}
        for k in range(3, p + 1, 2):  # odd cycle lengths only
            tr_k = sum(e ** k for e in eigs)
            c_k = round(tr_k.real / k)  # tr(A^k) = k * c_k for prime cycles
            cycles[k] = c_k

        # Also compute even traces (should be 0 for tournaments... no, not zero for walks)
        # Actually tr(A^2) = sum |lambda_k|^2 = number of edges from i to j and back = 0
        # No, tr(A^2) = sum_k lambda_k^2 = sum_{i,j} a_{ij} a_{ji}
        # For tournament: a_{ij} a_{ji} = 0 since exactly one of (i,j), (j,i) is an edge
        # So tr(A^2) = 0. Let's verify.
        tr2 = sum(e ** 2 for e in eigs)

        print(f"    S={str(list(S)):>20}: H={H:>10}, "
              f"c3={cycles.get(3, 0)}, c5={cycles.get(5, 0)}, "
              f"c7={cycles.get(7, 0)}, tr2={tr2.real:.6f}")

    # Check if cycle tuple determines H
    print(f"\n  Do cycle counts determine H?")
    cycle_to_H = defaultdict(set)
    for d in data:
        eigs = circulant_eigenvalues(p, d['S'])
        cycle_tuple = tuple(
            round(sum(e ** k for e in eigs).real / k)
            for k in range(3, p + 1, 2)
        )
        cycle_to_H[cycle_tuple].add(d['H'])

    n_cycle_classes = len(cycle_to_H)
    all_determined = all(len(v) == 1 for v in cycle_to_H.values())
    print(f"    {n_cycle_classes} distinct cycle-count tuples")
    print(f"    {len(set(d['H'] for d in data))} distinct H values")
    print(f"    Cycle counts determine H? {all_determined}")

    if not all_determined:
        print(f"    FAILURE CASES:")
        for ct, Hs in cycle_to_H.items():
            if len(Hs) > 1:
                print(f"      cycles {ct} -> H in {Hs}")


# ======================================================================
# SECTION 5: SPECTRAL VARIATIONAL — HIGHER PRIMES
# ======================================================================

def test_higher_primes():
    """Test the variational principle at p=17 (mod4=1) and p=19 (mod4=3)."""
    print(f"\n{'=' * 70}")
    print(f"HIGHER PRIME TESTS — p = 17, 19")
    print(f"{'=' * 70}")

    for p in [17, 19]:
        print(f"\n--- p = {p} (mod 4 = {p % 4}) ---")
        all_S = all_circulant_tournaments(p)
        n_tours = len(all_S)
        m = (p - 1) // 2
        print(f"  {n_tours} circulant tournaments, {m} spectral params")

        if n_tours > 1000:
            print(f"  TOO MANY for exhaustive H computation. Sampling...")
            # Sample a subset
            import random
            random.seed(42)
            sample = random.sample(all_S, min(200, n_tours))
        else:
            sample = all_S

        # Key tournaments to always include
        qr = tuple(sorted(qr_set(p)))
        ci = cyclic_interval(p)
        qr_valid = (p - 1) not in qr

        special = [ci]
        if qr_valid:
            special.append(qr)
        for s in special:
            if s not in sample:
                sample.append(s)

        t0 = time.time()
        results = []
        for i, S in enumerate(sample):
            H = ham_count_dp(p, S)
            y2 = spectral_y_squared(p, S)
            center_y2 = p / 4.0
            dist2 = sum((y - center_y2) ** 2 for y in y2)
            spread = max(y2) - min(y2)
            sum_y4 = sum(y ** 2 for y in y2)
            results.append({
                'S': S, 'H': H, 'y2': y2, 'dist2': dist2,
                'spread': spread, 'sum_y4': sum_y4
            })
            if (i + 1) % 50 == 0:
                print(f"    Computed {i + 1}/{len(sample)} ({time.time() - t0:.1f}s)")

        elapsed = time.time() - t0
        print(f"  Total: {len(results)} tournaments in {elapsed:.1f}s")

        results.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(r['H'] for r in results), reverse=True)
        print(f"  {len(H_vals)} distinct H values in sample")
        print(f"  H range: [{H_vals[-1]}, {H_vals[0]}]")
        print(f"  Max H: S = {results[0]['S']}, H = {results[0]['H']}")
        print(f"  y^2 = {[round(y, 4) for y in results[0]['y2']]}")
        print(f"  spread = {results[0]['spread']:.4f}, dist_center = {results[0]['dist2']:.4f}")

        # QR analysis
        if qr_valid:
            qr_result = next((r for r in results if r['S'] == qr), None)
            if qr_result:
                rank = sorted([r['H'] for r in results], reverse=True).index(qr_result['H']) + 1
                print(f"  QR: H = {qr_result['H']}, rank = {rank}/{len(results)}, "
                      f"spread = {qr_result['spread']:.4f}")
                print(f"    IS MAX: {qr_result['H'] == H_vals[0]}")

        # CI analysis
        ci_result = next((r for r in results if r['S'] == ci), None)
        if ci_result:
            rank = sorted([r['H'] for r in results], reverse=True).index(ci_result['H']) + 1
            print(f"  CI: H = {ci_result['H']}, rank = {rank}/{len(results)}, "
                  f"spread = {ci_result['spread']:.4f}")
            print(f"    IS MAX: {ci_result['H'] == H_vals[0]}")

        # Flattest
        flattest = min(results, key=lambda r: r['spread'])
        flat_rank = sorted([r['H'] for r in results], reverse=True).index(flattest['H']) + 1
        print(f"  Flattest: H = {flattest['H']}, rank = {flat_rank}/{len(results)}, "
              f"spread = {flattest['spread']:.4f}")
        print(f"    IS MAX: {flattest['H'] == H_vals[0]}")

        # Correlation
        H_list = [r['H'] for r in results]
        dist_list = [r['dist2'] for r in results]
        n = len(results)
        mean_H = sum(H_list) / n
        mean_d = sum(dist_list) / n
        cov = sum((h - mean_H) * (d - mean_d) for h, d in zip(H_list, dist_list)) / n
        var_H = sum((h - mean_H) ** 2 for h in H_list) / n
        var_d = sum((d - mean_d) ** 2 for d in dist_list) / n
        if var_H > 0 and var_d > 0:
            corr = cov / (var_H ** 0.5 * var_d ** 0.5)
        else:
            corr = 0
        print(f"  corr(H, dist_from_center) = {corr:.6f}")

        # Winner classification
        if flattest['H'] == H_vals[0]:
            print(f"  ==> FLATTEST WINS (p mod 4 = {p % 4})")
        elif qr_valid and qr_result and qr_result['H'] == H_vals[0]:
            print(f"  ==> QR/PALEY WINS (p mod 4 = {p % 4})")
        elif ci_result and ci_result['H'] == H_vals[0]:
            print(f"  ==> CYCLIC INTERVAL WINS (p mod 4 = {p % 4})")
        else:
            print(f"  ==> OTHER WINS (p mod 4 = {p % 4})")


# ======================================================================
# MAIN
# ======================================================================

def main():
    print("=" * 70)
    print("SPECTRAL H FORMULA — SEARCH FOR H = f(eigenvalues)")
    print("=" * 70)

    all_results = {}
    for p in [3, 5, 7, 11, 13]:
        data, H_vals = analyze_prime(p)
        all_results[p] = (data, H_vals)

        if p >= 7:
            variational_test(p, data)

        if p <= 7:
            group_theoretic_analysis(p, data)

        if p >= 7:
            trace_formula_analysis(p, data)

    # Grand summary
    print(f"\n{'=' * 70}")
    print(f"GRAND SUMMARY: SPECTRAL DETERMINATION OF H")
    print(f"{'=' * 70}")
    print(f"\n  {'p':>3} {'mod4':>4} {'n_H':>4} {'ps_det':>7} {'flat=max':>8} {'winner':<15}")
    print(f"  {'-' * 3} {'-' * 4} {'-' * 4} {'-' * 7} {'-' * 8} {'-' * 15}")

    for p in [3, 5, 7, 11, 13]:
        data, H_vals = all_results[p]
        m = (p - 1) // 2
        # Power sum determination
        ps_to_H = defaultdict(set)
        for d in data:
            ps_key = tuple(round(x, 6) for x in d['power_sums'][:m])
            ps_to_H[ps_key].add(d['H'])
        ps_det = all(len(v) == 1 for v in ps_to_H.values())

        # Flat = max?
        flattest = min(data, key=lambda d: d['spread'])
        flat_max = flattest['H'] == H_vals[0]

        qr = tuple(sorted(qr_set(p)))
        ci = cyclic_interval(p)
        qr_valid = (p - 1) not in qr

        if flat_max:
            if qr_valid and any(d['S'] == qr and d['H'] == H_vals[0] for d in data):
                winner = "PALEY/QR"
            else:
                winner = "FLATTEST"
        else:
            winner = "CONCENTRATED"

        print(f"  {p:>3} {p % 4:>4} {len(H_vals):>4} {str(ps_det):>7} {str(flat_max):>8} {winner:<15}")

    # Higher primes
    print(f"\n{'=' * 70}")
    print(f"HIGHER PRIME TESTS")
    print(f"{'=' * 70}")
    test_higher_primes()


if __name__ == '__main__':
    main()
