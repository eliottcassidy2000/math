#!/usr/bin/env python3
"""
Complete positivity analysis of U_T for small tournaments.
Instance: opus-2026-03-06-S10

Tests: p-positivity (OCF), h-positivity, e-positivity, Schur-positivity.
Tracks correlation with tiling grid symmetry (SC, SF, SC+SF kernel).

Key finding from h_positivity_test.py: h-positivity FAILS for all non-transitive
tournaments. Only the transitive tournament (unique total order) is h-positive.
This means the (3+1)-free poset → h-positivity pathway does NOT apply to tournaments
(since non-transitive tournaments don't correspond to posets).

Now testing: does e-positivity also fail? The e-basis is related to the h-basis
by the omega involution: h_lambda = omega(e_lambda). So e-positivity of U_T
is equivalent to h-positivity of omega(U_T) = U_{T-bar} = U_{T^op} (for tournaments).
Since U_T = U_{T^op}, e-positivity of U_T ⟺ h-positivity of U_T.
Therefore e-positivity also FAILS for non-transitive tournaments!

This means: for the Schur expansion, we need to check directly.
Irving-Omar (T073) showed s-positivity fails for general digraphs.
But does it hold for TOURNAMENTS specifically?
"""

from itertools import permutations
from math import factorial
from fractions import Fraction
from collections import defaultdict

def tournament_from_bits(bits, n):
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                adj[(i,j)] = 0
    for i in range(n-1):
        adj[(i,i+1)] = 1
        adj[(i+1,i)] = 0
    idx = 0
    for gap in range(2, n):
        for i in range(n - gap):
            j = i + gap
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(j,i)] = 1
                adj[(i,j)] = 0
            idx += 1
    return adj

def get_cycles(perm):
    n = len(perm)
    visited = [False]*n
    cycles = []
    for start in range(n):
        if visited[start]: continue
        cycle = []
        v = start
        while not visited[v]:
            visited[v] = True
            cycle.append(v)
            v = perm[v]
        if len(cycle) > 0:
            cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle_in(cycle, adj):
    k = len(cycle)
    if k == 1: return True
    for i in range(k):
        if adj.get((cycle[i], cycle[(i+1)%k]), 0) != 1:
            return False
    return True

def compute_U_T_power_sum(n, adj):
    adj_op = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                adj_op[(i,j)] = adj.get((j,i), 0)

    coeffs = defaultdict(int)
    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple)
        cycles = get_cycles(perm)
        phi = 0
        valid = True
        for cycle in cycles:
            if len(cycle) == 1: continue
            in_T = is_directed_cycle_in(cycle, adj)
            in_Top = is_directed_cycle_in(cycle, adj_op)
            if not in_T and not in_Top:
                valid = False; break
            if in_Top and not in_T:
                phi += len(cycle) - 1
        if not valid: continue
        cycle_lens = sorted([len(c) for c in cycles], reverse=True)
        partition = tuple(cycle_lens)
        coeffs[partition] += (-1)**phi
    return dict(coeffs)

def generate_partitions(n):
    if n == 0:
        yield ()
        return
    def helper(n, max_part):
        if n == 0:
            yield ()
            return
        for k in range(min(n, max_part), 0, -1):
            for rest in helper(n-k, k):
                yield (k,) + rest
    yield from helper(n, n)

def z_lambda(partition):
    """z_lambda = prod k^{m_k} * m_k! where m_k = multiplicity of k."""
    from collections import Counter
    cnt = Counter(partition)
    result = 1
    for k, mk in cnt.items():
        result *= (k ** mk) * factorial(mk)
    return result

def character_table(n):
    """
    Compute character table of S_n.
    chi[lambda][mu] = character of irrep lambda evaluated at conjugacy class mu.
    Uses the Murnaghan-Nakayama rule for small n.
    For n <= 5, we can hardcode or compute.
    """
    partitions = list(generate_partitions(n))
    # For small n, use a simple recursive computation
    # chi^lambda(mu) via the Murnaghan-Nakayama rule

    # For n=3:
    if n == 3:
        # Partitions: (3), (2,1), (1,1,1)
        # Conjugacy classes (same partitions): (3), (2,1), (1,1,1)
        # Character table:
        #        (1,1,1)  (2,1)  (3)
        # (3)       1       1     1
        # (2,1)     2       0    -1
        # (1,1,1)   1      -1     1
        return {
            ((3,), (1,1,1)): 1, ((3,), (2,1)): 1, ((3,), (3,)): 1,
            ((2,1), (1,1,1)): 2, ((2,1), (2,1)): 0, ((2,1), (3,)): -1,
            ((1,1,1), (1,1,1)): 1, ((1,1,1), (2,1)): -1, ((1,1,1), (3,)): 1,
        }

    if n == 4:
        # Partitions: (4), (3,1), (2,2), (2,1,1), (1,1,1,1)
        # Character table (rows=irreps, cols=conjugacy classes ordered as partitions):
        #              (1^4) (2,1^2) (2^2) (3,1) (4)
        # (4)            1     1      1     1    1
        # (3,1)          3     1     -1     0   -1
        # (2,2)          2     0      2    -1    0
        # (2,1,1)        3    -1     -1     0    1
        # (1,1,1,1)      1    -1      1     1   -1
        parts = [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]
        classes = [(1,1,1,1), (2,1,1), (2,2), (3,1), (4,)]
        table = [
            [1, 1, 1, 1, 1],
            [3, 1, -1, 0, -1],
            [2, 0, 2, -1, 0],
            [3, -1, -1, 0, 1],
            [1, -1, 1, 1, -1],
        ]
        result = {}
        for i, lam in enumerate(parts):
            for j, mu in enumerate(classes):
                result[(lam, mu)] = table[i][j]
        return result

    if n == 5:
        parts = [(5,), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)]
        classes = [(1,1,1,1,1), (2,1,1,1), (2,2,1), (3,1,1), (3,2), (4,1), (5,)]
        table = [
            [1, 1, 1, 1, 1, 1, 1],       # (5)
            [4, 2, 0, 1, -1, 0, -1],      # (4,1)
            [5, 1, 1, -1, 1, -1, 0],      # (3,2)
            [6, 0, -2, 0, 0, 0, 1],       # (3,1,1)
            [5, -1, 1, -1, -1, 1, 0],     # (2,2,1)
            [4, -2, 0, 1, 1, 0, -1],      # (2,1,1,1)
            [1, -1, 1, 1, -1, -1, 1],     # (1,1,1,1,1)
        ]
        result = {}
        for i, lam in enumerate(parts):
            for j, mu in enumerate(classes):
                result[(lam, mu)] = table[i][j]
        return result

    return None

def p_to_schur(n, p_coeffs):
    """Convert p-expansion to Schur expansion using character table.
    s_lambda = sum_mu (chi^lambda(mu) / z_mu) p_mu
    So [s_lambda]f = sum_mu (chi^lambda(mu) / z_mu) [p_mu]f
    """
    chi = character_table(n)
    if chi is None:
        return None

    partitions = list(generate_partitions(n))
    s_coeffs = {}
    for lam in partitions:
        coeff = Fraction(0)
        for mu in partitions:
            if mu not in p_coeffs:
                continue
            key = (lam, mu)
            if key not in chi:
                continue
            coeff += Fraction(chi[key]) * Fraction(p_coeffs[mu]) / Fraction(z_lambda(mu))
        s_coeffs[lam] = coeff
    return s_coeffs

def count_ham_paths(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj.get((v, u), 0):
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                    if i < j and i < k:
                        count += 1
    return count

def is_sc_necessary(n, adj):
    """Necessary condition for self-converse: palindromic score sequence."""
    scores = sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)])
    return scores == [n-1-s for s in reversed(scores)]

def is_self_flip(bits, n):
    """Check if flipping all non-path arcs gives isomorphic tournament."""
    m = n*(n-1)//2 - (n-1)
    flipped = bits ^ ((1 << m) - 1)
    # Check if same iso class (approximate: same scores and H)
    adj1 = tournament_from_bits(bits, n)
    adj2 = tournament_from_bits(flipped, n)
    s1 = sorted([sum(adj1.get((i,j),0) for j in range(n) if j!=i) for i in range(n)])
    s2 = sorted([sum(adj2.get((i,j),0) for j in range(n) if j!=i) for i in range(n)])
    h1 = count_ham_paths(n, adj1)
    h2 = count_ham_paths(n, adj2)
    return s1 == s2 and h1 == h2  # Necessary, not sufficient

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def main():
    print("=== COMPLETE POSITIVITY ANALYSIS OF U_T ===\n")

    for n in range(3, 6):
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        m = n*(n-1)//2 - (n-1)

        # Collect data
        seen = {}
        results = []

        for bits in range(2**m):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
            c3 = count_3cycles(n, adj)

            key = (scores, H, c3)
            if key in seen:
                continue
            seen[key] = bits

            p_coeffs = compute_U_T_power_sum(n, adj)
            s_coeffs = p_to_schur(n, p_coeffs)

            sc = is_sc_necessary(n, adj)
            sf = is_self_flip(bits, n)

            # Check p-positivity (OCF guarantees this)
            p_pos = all(v >= 0 for v in p_coeffs.values())

            # Check s-positivity
            s_pos = s_coeffs is not None and all(v >= 0 for v in s_coeffs.values())

            results.append({
                'bits': bits, 'H': H, 'scores': scores, 'c3': c3,
                'sc': sc, 'sf': sf,
                'p_coeffs': p_coeffs, 's_coeffs': s_coeffs,
                'p_pos': p_pos, 's_pos': s_pos,
            })

        # Print summary table
        print(f"\n{'Bits':>6} {'H':>3} {'c3':>3} {'SC':>3} {'SF':>3} {'p+':>3} {'s+':>3} | Schur expansion")
        print('-'*70)

        for r in sorted(results, key=lambda x: x['H']):
            s_str = ''
            if r['s_coeffs']:
                for p in sorted(r['s_coeffs'].keys()):
                    c = r['s_coeffs'][p]
                    if c != 0:
                        s_str += f"{c}*s{partition_str(p)} "

            print(f"{r['bits']:>6} {r['H']:>3} {r['c3']:>3} {'Y' if r['sc'] else 'N':>3} {'Y' if r['sf'] else 'N':>3} "
                  f"{'Y' if r['p_pos'] else 'N':>3} {'Y' if r['s_pos'] else 'N':>3} | {s_str[:60]}")

        # Summary statistics
        total = len(results)
        p_pos_count = sum(1 for r in results if r['p_pos'])
        s_pos_count = sum(1 for r in results if r['s_pos'])

        print(f"\n  p-positive: {p_pos_count}/{total} ({100*p_pos_count/total:.0f}%)")
        print(f"  s-positive: {s_pos_count}/{total} ({100*s_pos_count/total:.0f}%)")

        # Check correlation: is s-positivity correlated with SC or H?
        if total > 1:
            sc_and_s_pos = sum(1 for r in results if r['sc'] and r['s_pos'])
            nsc_and_s_pos = sum(1 for r in results if not r['sc'] and r['s_pos'])
            sc_total = sum(1 for r in results if r['sc'])
            nsc_total = sum(1 for r in results if not r['sc'])
            if sc_total > 0:
                print(f"  SC tournaments: {sc_and_s_pos}/{sc_total} s-positive")
            if nsc_total > 0:
                print(f"  NSC tournaments: {nsc_and_s_pos}/{nsc_total} s-positive")

if __name__ == '__main__':
    main()
