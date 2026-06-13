"""
cycle_dominance_mechanism.py — Why does Paley/interval maximize H?

The OCF says H = I(Omega(T), 2) = sum_k 2^k * alpha_k
where alpha_k = # independent k-sets in the conflict graph Omega.

ROBUST MECHANISM HYPOTHESIS:
  The H-maximizer has the MOST total directed odd cycles (alpha_1 = |V(Omega)|).
  The 2*alpha_1 term dominates the independence polynomial at lambda=2.
  Alpha_1 is maximized by spectral flatness for p=3 mod 4 (Paley)
  and by spectral concentration for p=1 mod 4 (cyclic interval).

This script:
  1. Decomposes H into contributions from each cycle length
  2. Checks which cycle lengths the maximizer dominates at
  3. Tests whether alpha_1 alone predicts H ranking (Spearman correlation)
  4. Investigates the spectral mechanism for cycle count maximization

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
from collections import defaultdict
from itertools import combinations

import mpmath
mpmath.mp.dps = 30

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i] |= (1 << j)
    full_mask = (1 << n) - 1
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        if bin(mask).count('1') >= n:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            candidates = adj[v] & ~mask
            w = 0
            while candidates:
                if candidates & 1:
                    dp[(mask | (1 << w), w)] += cnt
                candidates >>= 1
                w += 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))


def directed_cycles_by_length(n, S):
    """Count directed cycles by length in circulant tournament.
    Returns {length: count} for all odd lengths 3, 5, ..., n.
    """
    S_set = set(S)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S_set:
                adj[i][j] = True

    cycle_counts = {}
    for k in range(3, n + 1, 2):  # odd lengths only
        # Count directed k-cycles through vertex subsets
        count = 0
        for subset in combinations(range(n), k):
            # Count Hamiltonian directed cycles in the induced subtournament
            # on this subset
            nodes = list(subset)
            node_idx = {v: i for i, v in enumerate(nodes)}
            # DP for Hamiltonian cycles: paths from node 0, returning to 0
            sub_adj = [[False]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    if a != b:
                        sub_adj[a][b] = adj[nodes[a]][nodes[b]]

            # Count directed Hamiltonian cycles starting and ending at node 0
            # using Held-Karp DP
            full = (1 << k) - 1
            dp_c = [[0]*k for _ in range(1 << k)]
            dp_c[1 << 0][0] = 1  # start at node 0

            for mask in range(1, 1 << k):
                for v in range(k):
                    if dp_c[mask][v] == 0 or not (mask & (1 << v)):
                        continue
                    for w in range(k):
                        if mask & (1 << w):
                            continue
                        if sub_adj[v][w]:
                            dp_c[mask | (1 << w)][w] += dp_c[mask][v]

            # Cycles = paths that return to node 0
            n_cycles = sum(dp_c[full][v] for v in range(1, k) if sub_adj[v][0])
            # Each cycle is counted k times (once per starting vertex within the subset)
            # but we fixed start=0, so each cycle is counted once per starting vertex
            # Actually: we start at node 0 of the SUBSET. A directed cycle visits
            # all k nodes. Starting at different nodes gives different orderings.
            # Since we start at the SMALLEST node (index 0), and cycles have k rotations,
            # but we only count starting at 0. So each DIRECTED cycle is counted
            # once per rotation = k times if we sum over all starting vertices.
            # Since we fix start=0 and the subset has nodes in sorted order,
            # each directed cycle is counted exactly once.
            # WAIT: n_cycles counts paths 0 -> ... -> v -> 0 for each v adjacent to 0.
            # This counts each directed Hamiltonian cycle exactly once (fix start, direction fixed).
            # But directed cycles = both orientations of undirected cycles.
            # In tournaments, exactly one direction works, so n_cycles is correct.

            count += n_cycles

        # Divide by n (rotational symmetry of circulant: each directed cycle
        # appears in n rotations of vertex labeling, but we counted over all subsets)
        # Actually no: we enumerated ALL C(n,k) subsets. A directed cycle uses
        # k specific vertices. Each directed cycle is in exactly 1 subset.
        # So count is correct as-is.
        # BUT: circulant symmetry means the total is always divisible by n/k...
        # Actually, for a general graph, the number of directed k-cycles is just count.
        cycle_counts[k] = count

    return cycle_counts


def trace_cycle_count(n, S):
    """Count k-cycles using eigenvalue traces (faster for large k)."""
    S_set = set(S)
    # Build adjacency matrix
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S_set:
                A[i][j] = 1

    # Compute A^k by matrix multiplication
    traces = {}
    Ak = [[1 if i == j else 0 for j in range(n)] for i in range(n)]  # identity

    for k in range(1, n + 1):
        # Ak = Ak * A
        new_Ak = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new_Ak[i][j] += Ak[i][l] * A[l][j]
        Ak = new_Ak
        traces[k] = sum(Ak[i][i] for i in range(n))

    return traces


def multiplicative_orbit(p, S):
    S_set = frozenset(S)
    orbit = {S_set}
    for a in range(2, p):
        if math.gcd(a, p) != 1:
            continue
        orbit.add(frozenset((a * s) % p for s in S))
    return orbit


def find_orbit_reps(p):
    all_S = all_circulant_tournaments(p)
    seen = set()
    reps = []
    for S in all_S:
        fs = frozenset(S)
        if fs in seen:
            continue
        orbit = multiplicative_orbit(p, S)
        reps.append(S)
        for member in orbit:
            seen.add(member)
    return reps


def main():
    print("=" * 70)
    print("CYCLE DOMINANCE MECHANISM FOR H-MAXIMIZATION")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n{'=' * 60}")
        print(f"p = {p} (mod 4 = {p % 4}), m = {m}")
        print(f"{'=' * 60}")

        reps = find_orbit_reps(p)
        data = []

        for S in reps:
            H = ham_count_dp(p, S)

            # Get cycle counts via trace
            traces = trace_cycle_count(p, S)

            # For odd k, tr(A^k)/k gives directed k-cycles
            # (exact for k <= some threshold due to trace simplicity)
            cycle_counts = {}
            for k in range(3, p + 1, 2):
                cycle_counts[k] = traces[k] // k

            total_cycles = sum(cycle_counts[k] for k in cycle_counts)

            # Spectral data
            pi = mpmath.pi
            y2 = [float(sum(mpmath.sin(2 * pi * k * s / p) for s in S) ** 2)
                  for k in range(1, m + 1)]
            dist_flat = sum((y - p/4)**2 for y in y2)

            # Check if Paley or interval
            qr = set(pow(a, 2, p) % p for a in range(1, p)) - {0}
            nqr = set(range(1, p)) - qr
            is_paley = set(S) == qr or set(S) == nqr
            is_interval = set(S) == set(range(p - m, p))

            orbit_size = len(multiplicative_orbit(p, S))

            data.append({
                'S': S, 'H': H, 'cycles': cycle_counts,
                'total_cycles': total_cycles, 'dist_flat': dist_flat,
                'is_paley': is_paley, 'is_interval': is_interval,
                'orbit_size': orbit_size
            })

        data.sort(key=lambda d: d['H'], reverse=True)

        # Print cycle decomposition table
        cycle_lengths = sorted(set(k for d in data for k in d['cycles']))
        header = f"  {'H':>10} {'total':>8}"
        for k in cycle_lengths:
            header += f" {'c_'+str(k):>7}"
        header += " {'flat':>6} {'paley':>6}"
        print(header)

        for d in data:
            row = f"  {d['H']:>10} {d['total_cycles']:>8}"
            for k in cycle_lengths:
                row += f" {d['cycles'].get(k, 0):>7}"
            row += f" {d['dist_flat']:>8.2f}"
            row += f" {'*' if d['is_paley'] else '':>3}"
            row += f" {'I' if d['is_interval'] else '':>3}"
            print(row)

        # Check: does the maximizer have max total_cycles?
        max_H = data[0]
        max_cycles = max(data, key=lambda d: d['total_cycles'])
        print(f"\n  H-maximizer: H={max_H['H']}, total_cycles={max_H['total_cycles']}")
        print(f"  Cycle-maximizer: H={max_cycles['H']}, total_cycles={max_cycles['total_cycles']}")
        if max_H['total_cycles'] == max_cycles['total_cycles']:
            print(f"  *** H-maximizer = cycle-maximizer! ***")
        else:
            print(f"  H-maximizer != cycle-maximizer")

        # Spearman rank correlation between H and total_cycles
        h_rank = {d['H']: i for i, d in enumerate(sorted(data, key=lambda d: d['H'], reverse=True))}
        c_rank = {d['total_cycles']: i for i, d in enumerate(sorted(data, key=lambda d: d['total_cycles'], reverse=True))}

        n = len(data)
        d_sq = sum((h_rank[d['H']] - c_rank[d['total_cycles']])**2 for d in data)
        rho = 1 - 6 * d_sq / (n * (n**2 - 1))
        print(f"  Spearman rank corr(H, total_cycles) = {rho:.4f}")

        # Per-length analysis: which cycle lengths does the maximizer dominate at?
        print(f"\n  Per-length dominance (maximizer vs others):")
        for k in cycle_lengths:
            max_ck = max(d['cycles'].get(k, 0) for d in data)
            maximizer_ck = max_H['cycles'].get(k, 0)
            rank_ck = sorted(set(d['cycles'].get(k, 0) for d in data), reverse=True).index(maximizer_ck) + 1
            print(f"    c_{k}: maximizer={maximizer_ck}, max={max_ck}, "
                  f"rank={rank_ck}/{len(set(d['cycles'].get(k, 0) for d in data))}")

    # ================================================================
    # SECTION 2: Trace formula and spectral mechanism
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: SPECTRAL MECHANISM FOR CYCLE MAXIMIZATION")
    print(f"{'=' * 60}")

    print("""
    For circulant tournament on Z_p with connection set S:
      tr(A^k) = sum_{j=0}^{p-1} lambda_j^k
              = m^k + sum_{j=1}^{p-1} lambda_j^k

    For Paley (p = 3 mod 4):
      lambda_j = (1/2)(chi(j)*g - 1) where g = Gauss sum, |g|^2 = p
      |lambda_j|^2 = (p+1)/4 - Re(chi(j)*g)/2
      Actually: lambda_j for Paley has |lambda_j| = sqrt(p)/2 for j != 0

    The k-th power sum:
      sum |lambda_j|^k = sum_{j=1}^{p-1} |lambda_j|^k
                       = (p-1) * (sqrt(p)/2)^k for Paley
                       <= (p-1) * max |lambda_j|^{k-2} * sum |lambda_j|^2 for others

    By Parseval: sum |lambda_j|^2 = p*(p-1)/4 (constant for all regular circulants)
    So sum |lambda_j|^k is bounded by (p-1) * max_j |lambda_j|^{k-2} * p/4.

    For Paley: max |lambda_j| = sqrt(p)/2 (minimum possible by averaging).
    For others: max |lambda_j| > sqrt(p)/2 (some eigenvalue is larger).

    Key: for k >= 3, sum |lambda_j|^k is LARGER when the eigenvalue magnitudes
    are more spread (by convexity of x^k for k > 2). This seems to CONTRADICT
    Paley having the most cycles!

    BUT: we're summing lambda_j^k, not |lambda_j|^k. The complex phases matter!
    tr(A^k) = sum lambda_j^k = sum |lambda_j|^k * exp(i*k*arg(lambda_j))
    The REAL PART of this sum is what matters (tr(A^k) is real).

    The phases exp(i*k*arg(lambda_j)) can cause CANCELLATION. Paley's
    "well-distributed" phases may cause LESS cancellation at odd k,
    giving a larger real sum.

    This is the spectral mechanism: Paley maximizes cycle counts not through
    eigenvalue MAGNITUDE uniformity, but through eigenvalue PHASE coherence
    at odd powers.
    """)

    # Verify: compare |sum lambda^k| vs sum |lambda|^k
    for p in [7, 11]:
        m = (p - 1) // 2
        reps = find_orbit_reps(p)
        print(f"\n  p={p}: Phase coherence analysis")
        pi = mpmath.pi

        for S in reps[:4]:  # first few reps
            H = ham_count_dp(p, S)

            # Eigenvalues
            eigs = []
            for j in range(p):
                lam = sum(mpmath.exp(2j * pi * j_idx * s / p)
                          for s in S for j_idx in [j])
                # Actually need to fix: j is loop variable, not index
                pass

            # Use proper eigenvalue computation
            lam_list = []
            for j in range(p):
                lam = sum(mpmath.exp(2j * mpmath.pi * mpmath.mpf(1j if False else 1) * j * s / p)
                          for s in S)
                # This is getting confused with complex j. Let me fix.
                pass

        # Simpler: just use the trace directly
        for S in reps:
            H = ham_count_dp(p, S)
            traces = trace_cycle_count(p, S)

            qr = set(pow(a, 2, p) % p for a in range(1, p)) - {0}
            nqr = set(range(1, p)) - qr
            is_paley = set(S) == qr or set(S) == nqr

            # Compute eigenvalue magnitudes
            omega = [complex(math.cos(2*math.pi*k/p), math.sin(2*math.pi*k/p))
                     for k in range(p)]
            eigs = []
            for j in range(p):
                lam = sum(omega[(j*s) % p] for s in S)
                eigs.append(lam)

            mags = [abs(e) for e in eigs[1:]]  # skip j=0
            max_mag = max(mags)
            min_mag = min(mags)
            mag_spread = max_mag - min_mag

            label = "Paley" if is_paley else f"S={list(S)[:3]}..."
            print(f"    {label:>20}: H={H:>6}, "
                  f"|lam| range=[{min_mag:.3f}, {max_mag:.3f}], "
                  f"spread={mag_spread:.3f}, "
                  f"tr3={traces[3]}, tr5={traces[5]}")

    # ================================================================
    # SECTION 3: Direct alpha_k enumeration at p=7
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: INDEPENDENCE POLYNOMIAL DECOMPOSITION (p=7)")
    print(f"{'=' * 60}")

    p = 7
    m = 3
    reps = find_orbit_reps(p)

    for S in reps:
        H = ham_count_dp(p, S)

        # Enumerate all directed odd cycles
        S_set = set(S)
        adj = [[False]*p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j-i) % p in S_set:
                    adj[i][j] = True

        cycles = []
        for k in range(3, p + 1, 2):
            for subset in combinations(range(p), k):
                nodes = list(subset)
                sub_adj = [[adj[nodes[a]][nodes[b]] for b in range(k)] for a in range(k)]
                # Find Hamiltonian cycles starting at node 0
                full = (1 << k) - 1
                dp_c = [[0]*k for _ in range(1 << k)]
                dp_c[1][0] = 1
                for mask in range(1, 1 << k):
                    for v in range(k):
                        if dp_c[mask][v] == 0 or not (mask & (1 << v)):
                            continue
                        for w in range(k):
                            if mask & (1 << w):
                                continue
                            if sub_adj[v][w]:
                                dp_c[mask | (1 << w)][w] += dp_c[mask][v]
                n_cyc = sum(dp_c[full][v] for v in range(1, k) if sub_adj[v][0])
                if n_cyc > 0:
                    for _ in range(n_cyc):
                        cycles.append(frozenset(subset))

        # Deduplicate cycles (each directed cycle counted once)
        unique_cycles = list(set(cycles))
        alpha_1 = len(unique_cycles)

        # Build conflict graph and count independent sets
        cycle_list = unique_cycles
        n_cyc = len(cycle_list)

        # Conflict: two cycles share a vertex
        conflict = [[False]*n_cyc for _ in range(n_cyc)]
        for i in range(n_cyc):
            for j in range(i + 1, n_cyc):
                if cycle_list[i] & cycle_list[j]:
                    conflict[i][j] = conflict[j][i] = True

        # Count independent sets by size
        alpha = [0] * (n_cyc + 1)
        for size in range(n_cyc + 1):
            if size == 0:
                alpha[0] = 1
                continue
            count = 0
            for subset in combinations(range(n_cyc), size):
                indep = True
                for a in range(len(subset)):
                    for b in range(a + 1, len(subset)):
                        if conflict[subset[a]][subset[b]]:
                            indep = False
                            break
                    if not indep:
                        break
                if indep:
                    count += 1
            alpha[size] = count
            if count == 0:
                break

        # I(Omega, 2) = sum alpha_k * 2^k
        I_omega = sum(alpha[k] * (2**k) for k in range(len(alpha)))

        print(f"\n  S={list(S)}: H={H}")
        print(f"    {n_cyc} directed odd cycles (unique vertex sets)")
        print(f"    alpha = {[alpha[k] for k in range(len(alpha)) if alpha[k] > 0]}")
        print(f"    I(Omega, 2) = {I_omega}")

        # Verify OCF
        if I_omega == H:
            print(f"    OCF VERIFIED: I(Omega, 2) = H")
        else:
            print(f"    *** OCF MISMATCH: I(Omega,2)={I_omega} != H={H} ***")

        # Decompose H into contributions
        print(f"    Decomposition: ", end="")
        for k in range(len(alpha)):
            if alpha[k] > 0:
                print(f"2^{k}*{alpha[k]}={alpha[k]*2**k}", end=" + ")
        print()


if __name__ == '__main__':
    main()
