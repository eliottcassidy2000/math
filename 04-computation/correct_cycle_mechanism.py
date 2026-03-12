"""
correct_cycle_mechanism.py — Fixed cycle enumeration for H-maximization analysis.

BUG IN PREVIOUS: trace formula tr(A^k)/k is NOT c_k for k >= 7 (non-simple walks).
Also: deduplicating cycles by vertex set loses multiplicity.

CORRECT APPROACH: Direct enumeration of directed Hamiltonian cycles on each
k-element subset, keeping each directed cycle as a separate vertex of Omega.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
from collections import defaultdict
from itertools import combinations

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


def enumerate_directed_cycles(n, adj_bool, k):
    """Count directed Hamiltonian cycles on each k-element subset.

    Returns list of (vertex_set, n_directed_cycles) pairs.
    adj_bool[i][j] = True if i->j.
    """
    results = []
    for subset in combinations(range(n), k):
        nodes = list(subset)
        # Count directed Hamiltonian cycles starting at smallest node
        # A directed Hamiltonian cycle visits all k nodes and returns to start
        sub_adj = [[adj_bool[nodes[a]][nodes[b]] for b in range(k)]
                   for a in range(k)]

        full = (1 << k) - 1
        dp_c = [[0] * k for _ in range(1 << k)]
        dp_c[1][0] = 1  # start at node 0

        for mask in range(1, 1 << k):
            for v in range(k):
                if dp_c[mask][v] == 0 or not (mask & (1 << v)):
                    continue
                for w in range(k):
                    if mask & (1 << w):
                        continue
                    if sub_adj[v][w]:
                        dp_c[mask | (1 << w)][w] += dp_c[mask][v]

        # Count paths that return to node 0
        n_cycles = sum(dp_c[full][v] for v in range(1, k) if sub_adj[v][0])

        if n_cycles > 0:
            results.append((frozenset(subset), n_cycles))

    return results


def build_omega_and_compute(n, S):
    """Build conflict graph Omega and compute I(Omega, 2).

    Omega vertices = directed odd cycles (each cycle is a separate vertex).
    Omega edges = cycles that share at least one tournament vertex.
    """
    S_set = set(S)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i][j] = True

    # Enumerate ALL directed odd cycles
    all_cycles = []  # list of (vertex_set, index_within_set)
    cycle_vertex_sets = []
    by_length = {}

    for k in range(3, n + 1, 2):
        cycle_data = enumerate_directed_cycles(n, adj, k)
        count = 0
        for vset, n_cyc in cycle_data:
            for idx in range(n_cyc):
                all_cycles.append(vset)
                count += 1
        by_length[k] = count

    n_cycles = len(all_cycles)

    # Build conflict matrix
    # Two cycles conflict if their vertex sets intersect
    # Note: cycles on the same vertex set ALWAYS conflict
    conflict = [[False] * n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i + 1, n_cycles):
            if all_cycles[i] & all_cycles[j]:  # vertex sets intersect
                conflict[i][j] = conflict[j][i] = True

    # Count independent sets by size using inclusion
    alpha = [0] * (n_cycles + 1)
    alpha[0] = 1

    # For efficiency, only compute up to the independence number
    for size in range(1, n_cycles + 1):
        count = 0
        for subset in combinations(range(n_cycles), size):
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

    # I(Omega, 2)
    I_omega = sum(alpha[k] * (2 ** k) for k in range(len(alpha)))

    return {
        'n_cycles': n_cycles,
        'by_length': by_length,
        'alpha': [alpha[k] for k in range(len(alpha)) if alpha[k] > 0 or k == 0],
        'I_omega': I_omega,
        'alpha_1': n_cycles
    }


def main():
    print("=" * 70)
    print("CORRECT CYCLE DOMINANCE MECHANISM")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n{'=' * 60}")
        print(f"p = {p} (mod 4 = {p % 4}), m = {m}")
        print(f"{'=' * 60}")

        reps = find_orbit_reps(p)
        data = []

        for S in reps:
            H = ham_count_dp(p, S)
            omega_data = build_omega_and_compute(p, S)

            qr = set(pow(a, 2, p) % p for a in range(1, p)) - {0}
            nqr = set(range(1, p)) - qr
            is_paley = set(S) == qr or set(S) == nqr
            is_interval = set(S) == set(range(p - m, p))

            data.append({
                'S': S, 'H': H,
                'total_cycles': omega_data['n_cycles'],
                'by_length': omega_data['by_length'],
                'alpha': omega_data['alpha'],
                'I_omega': omega_data['I_omega'],
                'is_paley': is_paley,
                'is_interval': is_interval
            })

            # Verify OCF
            if omega_data['I_omega'] != H:
                print(f"  *** OCF MISMATCH at S={list(S)}: "
                      f"I(Omega,2)={omega_data['I_omega']} != H={H} ***")

        data.sort(key=lambda d: d['H'], reverse=True)

        # Print table
        cycle_lengths = sorted(set(k for d in data for k in d['by_length']))
        header = f"  {'H':>7} {'alpha_1':>7}"
        for k in cycle_lengths:
            header += f" {'c'+str(k):>5}"
        header += f" {'alpha':>20} {'P':>2} {'I':>2}"
        print(header)

        for d in data:
            row = f"  {d['H']:>7} {d['total_cycles']:>7}"
            for k in cycle_lengths:
                row += f" {d['by_length'].get(k, 0):>5}"
            row += f" {str(d['alpha']):>20}"
            row += f" {'*' if d['is_paley'] else '':>2}"
            row += f" {'I' if d['is_interval'] else '':>2}"
            print(row)

        # H-max vs alpha_1-max
        max_H = data[0]
        max_alpha = max(data, key=lambda d: d['total_cycles'])
        print(f"\n  H-maximizer: H={max_H['H']}, alpha_1={max_H['total_cycles']}")
        print(f"  Alpha_1-maximizer: H={max_alpha['H']}, alpha_1={max_alpha['total_cycles']}")
        if max_H['total_cycles'] == max_alpha['total_cycles']:
            print(f"  *** H-max = alpha_1-max confirmed! ***")

        # Per-length dominance
        print(f"\n  Per-length dominance of H-maximizer:")
        for k in cycle_lengths:
            max_ck = max(d['by_length'].get(k, 0) for d in data)
            maximizer_ck = max_H['by_length'].get(k, 0)
            print(f"    c_{k}: maximizer={maximizer_ck}, max={max_ck}, "
                  f"{'WIN' if maximizer_ck == max_ck else 'LOSE'}")

        # Spearman rank correlation
        h_vals = [(d['H'], d['total_cycles']) for d in data]
        h_vals.sort(key=lambda x: -x[0])
        h_ranks = {h: i for i, (h, _) in enumerate(h_vals)}
        c_vals = sorted(h_vals, key=lambda x: -x[1])
        c_ranks = {c: i for i, (_, c) in enumerate(c_vals)}

        n = len(data)
        if n > 1:
            d_sq = sum((h_ranks[d['H']] - c_ranks[d['total_cycles']]) ** 2 for d in data)
            rho = 1 - 6 * d_sq / (n * (n ** 2 - 1))
            print(f"  Spearman corr(H, alpha_1) = {rho:.4f}")

    # ================================================================
    # SECTION 2: The alpha_1 dominance mechanism
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: WHY ALPHA_1 DETERMINES H")
    print(f"{'=' * 60}")
    print("""
    H = I(Omega, 2) = sum_{k=0}^{alpha(G)} alpha_k * 2^k

    = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

    The question: does the LINEAR TERM 2*alpha_1 dominate?

    If alpha_2 << alpha_1, then H ~ 1 + 2*alpha_1.
    This holds when Omega is DENSE (many conflicts between cycles),
    so few disjoint pairs exist.

    For tournaments: odd cycles tend to share vertices (pigeonhole),
    making Omega dense. The independence number alpha(Omega) is small.
    """)

    for d in data:
        contributions = [(k, d['alpha'][k] * 2**k if k < len(d['alpha']) else 0)
                        for k in range(len(d['alpha']))]
        total = sum(c for _, c in contributions)
        print(f"  S={list(d['S'])}: H={d['H']}")
        for k, c in contributions:
            if c > 0:
                pct = 100.0 * c / total
                print(f"    2^{k} * alpha_{k} = {c:>6} ({pct:>5.1f}%)")

    # ================================================================
    # SECTION 3: Paley cycle count advantage — spectral explanation
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: SPECTRAL EXPLANATION OF PALEY CYCLE DOMINANCE")
    print(f"{'=' * 60}")

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        reps = find_orbit_reps(p)
        for S in reps:
            H = ham_count_dp(p, S)
            qr = set(pow(a, 2, p) % p for a in range(1, p)) - {0}
            nqr = set(range(1, p)) - qr
            is_paley = set(S) == qr or set(S) == nqr

            # Eigenvalues
            import cmath
            omega = cmath.exp(2j * cmath.pi / p)
            eigs = [sum(omega ** (j * s) for s in S) for j in range(p)]

            # Compute sum lambda^k for odd k
            for k in [3, 5]:
                sum_lam_k = sum(e ** k for e in eigs)
                label = "Paley" if is_paley else f"S={list(S)[:3]}..."
                print(f"    {label:>15}: sum lam^{k} = {sum_lam_k.real:>10.2f} "
                      f"(imag: {sum_lam_k.imag:.2e})")

            # Also: for flat spectrum, sum |lam|^k * cos(k*arg(lam))
            # Paley: |lam_j| = sqrt(p)/2 for j != 0
            # So sum lam^k = m^k + (sqrt(p)/2)^k * sum_{j=1}^{p-1} exp(i*k*arg(lam_j))
            if is_paley:
                phases = [cmath.phase(e) for e in eigs[1:]]
                for k in [3, 5]:
                    phase_sum = sum(cmath.exp(1j * k * ph) for ph in phases)
                    mag_factor = (p**0.5 / 2) ** k
                    print(f"    Paley phase analysis at k={k}:")
                    print(f"      (sqrt(p)/2)^k = {mag_factor:.4f}")
                    print(f"      sum exp(ik*theta) = {phase_sum.real:.4f} + {phase_sum.imag:.4f}i")
                    print(f"      Contribution: {mag_factor * phase_sum.real:.4f}")
                    print(f"      Total (with m^k={m**k}): {m**k + mag_factor * phase_sum.real:.4f}")


if __name__ == '__main__':
    main()
