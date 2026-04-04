"""
h_ie_deep.py — opus-2026-04-04-S3

DEEP INCLUSION-EXCLUSION SESSION.

Goal 1: Compute log(H(t)) and study its multilinear expansion.
  If H = exp(something simple), the structure of that "something" is fundamental.

Goal 2: Analyze the ZERO PATTERN: which tile subsets S have c_S = 0?
  This should have a combinatorial explanation via cycle packings.

Goal 3: Express c_S via a matrix formula (permanent/determinant).

Goal 4: Understand the full inclusion-exclusion cancellation mechanism.
"""
import numpy as np
from math import log, exp, comb
from itertools import combinations
from collections import defaultdict
from fractions import Fraction

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=int)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def compute_all_H(n):
    tiles = make_tiles(n)
    m = len(tiles)
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)
    return H, tiles, m

def mobius_inversion(H, m):
    coeffs = {}
    for S in range(1 << m):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            coeffs[S] = val
    return coeffs

def goal1_log_structure(n):
    """Compute log(H(t)) and its multilinear expansion."""
    H, tiles, m = compute_all_H(n)

    print(f"\n{'='*60}")
    print(f"GOAL 1: LOG STRUCTURE of H(t) at n={n}")
    print(f"{'='*60}")

    # Compute log(H(t)) for each tiling
    logH = np.log(H.astype(float))

    # Möbius inversion of log(H)
    log_coeffs = {}
    for S in range(1 << m):
        val = 0.0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * logH[T]
            if T == 0: break
            T = (T - 1) & S
        if abs(val) > 1e-10:
            log_coeffs[S] = val

    print(f"Nonzero log-coefficients: {len(log_coeffs)} (vs {len(mobius_inversion(H,m))} for H)")

    # Group by order
    by_order = defaultdict(list)
    for S, c in log_coeffs.items():
        by_order[bin(S).count('1')].append((S, c))

    for order in sorted(by_order.keys()):
        entries = by_order[order]
        values = [c for _, c in entries]
        print(f"\n  Order {order}: {len(entries)} nonzero")
        if order == 0:
            print(f"    log(H(transitive)) = log(1) = {values[0]:.6f}")
        elif order <= 2:
            for S, c in sorted(entries, key=lambda x: abs(x[1]), reverse=True)[:8]:
                ts = [tiles[i] for i in range(m) if S & (1 << i)]
                print(f"    {ts}: {c:.6f}")

    # KEY QUESTION: Is log(H) simpler than H?
    # Compare max degrees
    h_coeffs = mobius_inversion(H, m)
    h_max_deg = max(bin(S).count('1') for S in h_coeffs.keys())
    log_max_deg = max(bin(S).count('1') for S in log_coeffs.keys())
    print(f"\n  H max degree: {h_max_deg}")
    print(f"  log(H) max degree: {log_max_deg}")

    # Are log coefficients related to SINGLE CYCLE indicators?
    # The cluster expansion says: log I(G, λ) = Σ connected induced subgraphs...
    # For the independence polynomial: log I(G, x) = Σ_k (-1)^{k+1} n_k x^k / k
    # where n_k relates to connected induced subgraphs of complement.
    # But this is for the static graph, not the tiling-dependent one.

    # Instead, check: is log(H(t)) ≈ 2 * Σ_C χ_C(t)?
    # i.e., does log I(Ω, 2) ≈ 2 * (number of cycles)?
    # More precisely, log(1 + 2x + 4x^2 + ...) where x tracks cycles

    return log_coeffs

def goal2_zero_pattern(n):
    """Analyze which tile subsets have c_S = 0."""
    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    print(f"\n{'='*60}")
    print(f"GOAL 2: ZERO PATTERN at n={n}")
    print(f"{'='*60}")

    total_subsets = 1 << m
    nonzero = len(coeffs)
    zero = total_subsets - nonzero
    print(f"Total subsets: {total_subsets}")
    print(f"Nonzero: {nonzero} ({100*nonzero/total_subsets:.1f}%)")
    print(f"Zero: {zero} ({100*zero/total_subsets:.1f}%)")

    # By order
    for k in range(m+1):
        total_k = comb(m, k)
        nz_k = sum(1 for S in coeffs if bin(S).count('1') == k)
        z_k = total_k - nz_k
        if total_k > 0:
            print(f"  Order {k}: {nz_k}/{total_k} nonzero ({100*nz_k/total_k:.0f}%)")

    # CHARACTERIZE THE ZEROS: For each zero subset S, what prevents cycle packings?
    # A tile subset S has c_S = 0 iff the signed cycle packings cancel.
    # But if NO cycle packing exists for S, then c_S = 0 trivially.

    # For a cycle packing to exist for S:
    # Need vertex-disjoint odd cycles C₁,...,C_r covering all tiles in S
    # such that B(C_i) ⊆ S and S \ B(C_i) ⊆ F(C_i)

    # The simplest necessary condition: every tile in S must belong to some odd cycle
    # that uses it (forward or backward), and the cycle's other backward tiles must also be in S.

    # Actually, the simplest condition is:
    # S is "reachable" by cycle packings iff there exist vertex-disjoint cycles
    # whose tile sets cover S, with each cycle's backward tiles ⊆ S.

    # Let's check: among zero subsets, how many are "unreachable" (no packing exists)?

    # For this, I need to enumerate valid cycles and their tile requirements.
    # Let me use the OCF machinery.

    from h_ocf_multilinear_v2 import enumerate_directed_cycles, cycle_indicator
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    max_cyc_len = n if n % 2 == 1 else n - 1
    all_cycles = []
    for k in range(3, max_cyc_len + 1, 2):
        for cyc in enumerate_directed_cycles(n, k):
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    # For each cycle, identify which tile subsets it can "serve"
    # A cycle C can serve tile subset S_C iff B(C) ⊆ S_C ⊆ F(C)∪B(C)
    cycle_servable = []
    for cyc, chi in all_cycles:
        # B(C) = tiles used backward, F(C) = tiles used forward
        B_C = set()
        F_C = set()
        L = len(cyc)
        for i in range(L):
            u, v = cyc[i], cyc[(i+1) % L]
            if abs(u - v) == 1: continue  # base path
            high, low = max(u,v), min(u,v)
            tile_label = (high+1, low+1)
            if tile_label not in tile_to_idx: continue
            k = tile_to_idx[tile_label]
            if u > v:  # forward
                F_C.add(k)
            else:  # backward
                B_C.add(k)
        cycle_servable.append((frozenset(B_C), frozenset(F_C), frozenset(cyc)))

    print(f"\nCycle servability analysis:")
    print(f"  Total valid cycles: {len(cycle_servable)}")

    # For zero subsets of small size, check if any cycle packing covers them
    zero_subsets = [S for S in range(total_subsets) if S not in coeffs and S != 0]
    unreachable = 0
    reachable_but_cancel = 0

    for S_mask in zero_subsets[:200]:  # Check first 200 zeros
        S = frozenset(i for i in range(m) if S_mask & (1 << i))
        if not S: continue

        # Check if any cycle can serve some subset of S
        # A packing for S must cover ALL tiles in S
        has_packing = False

        # For small S, try all single-cycle coverings
        if len(S) <= 5:
            for B_C, F_C, verts in cycle_servable:
                if B_C <= S and S <= (F_C | B_C):
                    has_packing = True
                    break

            if not has_packing and len(S) >= 2:
                # Try 2-cycle packings
                for i1, (B1, F1, V1) in enumerate(cycle_servable):
                    for i2, (B2, F2, V2) in enumerate(cycle_servable):
                        if i2 <= i1: continue
                        if V1 & V2: continue  # vertex disjoint
                        S1_options = [s for s in range(1, len(S)) if True]  # TODO
                        # Simplified: check if union of servable sets covers S
                        combined_B = B1 | B2
                        combined_FuB = (F1 | B1) | (F2 | B2)
                        if combined_B <= S and S <= combined_FuB:
                            has_packing = True
                            break
                    if has_packing: break

        if has_packing:
            reachable_but_cancel += 1
        else:
            unreachable += 1

    checked = min(200, len(zero_subsets))
    print(f"\n  Zero subsets analyzed: {checked}")
    print(f"  Unreachable (no cycle packing exists): {unreachable}")
    print(f"  Reachable but cancel: {reachable_but_cancel}")
    print(f"  (Unchecked: {len(zero_subsets) - checked})")

    return coeffs

def goal3_matrix_formula(n):
    """Express c_S via a matrix. The key insight:
    H(t) = Σ_σ Π_{arcs of σ} M[u][v]
    where M is the adjacency matrix with entries depending on t.

    This is NOT a permanent (it sums over paths, not matchings).
    But it IS the trace of a transfer matrix product.

    Let's compute the transfer matrix formulation explicitly.
    """
    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    print(f"\n{'='*60}")
    print(f"GOAL 3: MATRIX FORMULA at n={n}")
    print(f"{'='*60}")

    # The adjacency matrix M(t) of T(t):
    # M[u][v] = 1 if u→v is a base-path arc
    # M[u][v] = 1-t_k if u→v is forward tile k
    # M[u][v] = t_k if u→v is backward tile k
    # M[u][v] = 0 otherwise (u=v or reverse base path)

    # H(t) = Σ_{all n! perms σ of {0,...,n-1}} Π_{i=0}^{n-2} M[σ(i)][σ(i+1)]
    # But only permutations that are valid paths count.

    # More efficiently: H(t) = sum over starting vertices of the (n-1)-step walk count
    # H(t) = 1^T · M^{n-1} · 1 (NOT exactly — need to avoid revisiting vertices)

    # The CORRECT formula uses inclusion-exclusion on visited vertices:
    # H(t) = Σ_{S ⊆ V, |S|=n} Σ_{σ: path on S} Π M[σ_i][σ_{i+1}]
    # But S = V always (Hamiltonian paths use all vertices).

    # Actually, the Held-Karp DP can be expressed as:
    # H(t) = Σ_v Σ_{mask = full} dp[mask][v]
    # dp[{v}][v] = 1
    # dp[mask][v] = Σ_u dp[mask\v][u] · M[u][v]

    # The matrix form: H = trace of product of "transfer matrices"?
    # Not quite standard. But there's a beautiful identity:

    # Let A be the symbolic adjacency matrix. Then:
    # H(t) = coefficient of x₁x₂...x_n in Permanent(diag(x) · A)
    # where diag(x) is the diagonal matrix with x₁,...,x_n.

    # This is the "permanent with vertex markers" — standard in enumeration.
    # But the permanent is #P-hard to compute, so this doesn't simplify things.

    # ALTERNATIVE: The "transfer matrix" approach works for special graph structures.
    # For tournaments on a totally ordered vertex set, the transfer matrix is upper/lower
    # triangular + perturbations from tiles.

    # Let me compute: for the BASE PATH (transitive tournament), the adjacency matrix
    # is strictly upper triangular. All eigenvalues are 0. So det(I - xA) = 1.
    # The characteristic polynomial det(λI - A) = λ^n.

    # When tiles are flipped, we add off-diagonal entries. The "perturbation" from each
    # tile flip changes the matrix.

    # KEY INSIGHT: H(t) = permanent-like quantity.
    # Permanent(A) = Σ_σ Π A[i][σ(i)] sums over bijections.
    # But Hamiltonian paths are ORDERED bijections (the order matters).

    # There IS a connection: H(T) equals the permanent of a specific matrix.
    # For a tournament T, define the n×n matrix P where:
    # P[i][j] = T[i][j] if j is the "next step" from i in a valid path.
    # But this is the adjacency matrix restricted to adjacent positions.

    # Actually, the SIMPLEST matrix formula:
    # H(T) = Σ_{i} (A^{n-1})_{i,*} . 1? No, A^{n-1} counts walks of length n-1,
    # not paths. The inclusion-exclusion to correct for revisits is exactly Held-Karp.

    # Let me try a different matrix approach: the IMMANANT.
    # H(T) relates to the trace of A in the symmetric group representation.
    # Specifically, H(T) is the permanent of A minus diagonal corrections.

    # For n=4, let me compute the symbolic adjacency matrix and its permanent.
    print(f"Symbolic adjacency matrix M(t) for n={n}:")
    for u in range(n):
        row = []
        for v in range(n):
            if u == v:
                row.append('0')
            elif u == v + 1:
                row.append('1')  # base path
            elif v == u + 1:
                row.append('0')  # reverse base path (never)
            else:
                high, low = max(u,v), min(u,v)
                tile_label = (high+1, low+1)
                if tile_label in tile_to_idx:
                    k = tile_to_idx[tile_label]
                    if u > v:
                        row.append(f'(1-t{k})')
                    else:
                        row.append(f't{k}')
                else:
                    row.append('?')
        if n <= 6:
            print(f"  {' '.join(f'{x:>8}' for x in row)}")

    # The permanent of this matrix sums over ALL bijections, including ones
    # that don't form paths. But for TOURNAMENTS, every bijection corresponds
    # to a permutation, and the path structure comes from the ordering.

    # Actually, permanent of the n×n tournament adjacency matrix counts
    # the number of "directed perfect matchings" — which for tournaments
    # is related to H but not identical.

    # Let me instead compute something more useful:
    # The DETERMINANT of (I - x*M) = Σ_k (-x)^k trace_k(M)
    # where trace_k involves cycles of length k.

    # The cycle structure is exactly what OCF captures!
    # det(I - x*M) = exp(-Σ_{k≥1} tr(M^k) x^k / k)

    # And tr(M^k) counts directed walks of length k (with revisits).
    # For k ≤ n, some of these are cycles (Hamiltonian or not).

    # This connection between determinant and cycle counting IS the
    # inclusion-exclusion behind OCF!

    # Let me verify: det(I - 2M) should relate to H somehow.
    # For the independence polynomial: I(G, x) = Σ_{k} α_k x^k
    # There's a matrix formula: I(G, x) = det(I + x*D - x*A) where D is degree matrix?
    # No, that's for the matching polynomial.

    # Actually for independence polynomial on a graph G:
    # I(G, x) can be computed as the permanent of (I + x·J - x·A)?
    # This isn't standard.

    # Let me try a DIRECT matrix computation.
    # For small n, compute det(I - x·M(t)) symbolically and see how it relates to H.

    if n <= 5:
        # Numerical computation at specific tilings
        for test_bits in [0, 1, (1<<m)-1]:
            T = tiling_to_tournament(n, tiles, test_bits)
            A = T.astype(float)
            # Compute det(I - 2A) and compare to... something
            det_val = np.linalg.det(np.eye(n) - 2*A)
            h_val = count_hp(T, n)
            # Compute trace of A^k for k=1,...,n
            traces = []
            Ak = np.eye(n)
            for k in range(1, n+1):
                Ak = Ak @ A
                traces.append(np.trace(Ak))

            print(f"\n  Tiling {test_bits}: H={h_val}, det(I-2A)={det_val:.1f}")
            print(f"    Traces: {[int(t) for t in traces]}")
            # Note: tr(A^k)/k counts directed k-cycles (with multiplicity from starting vertex)
            cycle_counts = [traces[k-1]/k for k in range(1, n+1)]
            print(f"    Cycle counts (tr(A^k)/k): {[f'{c:.1f}' for c in cycle_counts]}")

            # The Ihara-Bass formula: det(I - uA + u^2(D-I)) = ...
            # Or the simpler: det(I - uA) = exp(-Σ c_k u^k / k)
            # where c_k = number of directed closed walks of length k
            # (= tr(A^k))

            # So det(I-2A) = exp(-Σ tr(A^k) 2^k / k) = exp(-2*Σ c_k 2^{k-1}/k)
            # where c_k counts k-walks from vertex i back to i.

            log_det = -sum(traces[k-1] * 2**(k) / (k+1) for k in range(n) if k > 0)
            # Wait, det(I-uA) = exp(-Σ_{k≥1} tr(A^k) u^k / k)
            # So log det(I-2A) = -Σ_{k=1}^{∞} tr(A^k) 2^k / k
            # But for finite matrix, only k up to n matters for eigenvalues.
            log_det_approx = sum(-traces[k-1] * 2**k / k for k in range(1, n+1))
            print(f"    log det(I-2A): exact={np.log(abs(det_val)):.4f}, "
                  f"from traces={log_det_approx:.4f}")

def goal4_cancellation(n):
    """Understand the full cancellation structure.

    For each multilinear coefficient c_S, decompose into:
    - Positive contributions (from cycle packings with even total forward count)
    - Negative contributions (odd total forward count)
    - Net = c_S

    The CANCELLATION RATIO = (positive + |negative|) / |c_S| measures how much
    work the inclusion-exclusion does.
    """
    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    print(f"\n{'='*60}")
    print(f"GOAL 4: CANCELLATION ANALYSIS at n={n}")
    print(f"{'='*60}")

    # Compute the "raw" sum (without cancellation):
    # |c_S| via Möbius: c_S = Σ_{T⊆S} (-1)^{|S\T|} H(1_T)
    # Positive part: Σ_{T: |S\T| even} H(1_T)
    # Negative part: Σ_{T: |S\T| odd} H(1_T)

    by_order = defaultdict(list)
    for S in range(1 << m):
        k = bin(S).count('1')
        if k == 0: continue

        pos = neg = 0
        T = S
        while True:
            diff = bin(S ^ T).count('1')
            if diff % 2 == 0:
                pos += int(H[T])
            else:
                neg += int(H[T])
            if T == 0: break
            T = (T - 1) & S

        net = pos - neg  # = c_S
        raw = pos + neg
        cancel_ratio = raw / abs(net) if net != 0 else float('inf')

        by_order[k].append({
            'S': S, 'net': net, 'pos': pos, 'neg': neg,
            'raw': raw, 'ratio': cancel_ratio
        })

    for k in sorted(by_order.keys()):
        entries = by_order[k]
        nonzero = [e for e in entries if e['net'] != 0]
        zero = [e for e in entries if e['net'] == 0]
        ratios = [e['ratio'] for e in nonzero]

        print(f"\n  Order {k}: {len(nonzero)} nonzero, {len(zero)} zero")
        if ratios:
            print(f"    Cancel ratio: min={min(ratios):.1f}, "
                  f"max={max(ratios):.1f}, mean={sum(ratios)/len(ratios):.1f}")

        # Show extremes
        if nonzero:
            most_cancel = max(nonzero, key=lambda e: e['ratio'])
            least_cancel = min(nonzero, key=lambda e: e['ratio'])
            ts_most = [tiles[i] for i in range(m) if most_cancel['S'] & (1 << i)]
            ts_least = [tiles[i] for i in range(m) if least_cancel['S'] & (1 << i)]
            print(f"    Most cancellation: {ts_most} ratio={most_cancel['ratio']:.1f} "
                  f"({most_cancel['pos']}-{most_cancel['neg']}={most_cancel['net']})")
            print(f"    Least cancellation: {ts_least} ratio={least_cancel['ratio']:.1f} "
                  f"({least_cancel['pos']}-{least_cancel['neg']}={least_cancel['net']})")

        # For zero-coefficient subsets: pos = neg always
        if zero:
            all_equal = all(e['pos'] == e['neg'] for e in zero)
            print(f"    Zero subsets: pos=neg always: {all_equal}")

if __name__ == "__main__":
    for n in [4, 5, 6]:
        goal1_log_structure(n)

    goal2_zero_pattern(5)
    goal2_zero_pattern(6)

    goal3_matrix_formula(4)
    goal3_matrix_formula(5)

    goal4_cancellation(5)
    goal4_cancellation(6)
