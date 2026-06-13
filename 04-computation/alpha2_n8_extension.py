import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
alpha2_n8_extension.py
kind-pasteur-2026-03-07-S39b

At n=8, alpha_2 can include:
- (3-cycle, 3-cycle) disjoint pairs: 3+3=6 <= 8 vertices
- (5-cycle, 3-cycle) disjoint pairs: 5+3=8 vertices exactly
- (5-cycle, 5-cycle) disjoint pairs: 5+5=10 > 8, IMPOSSIBLE

So alpha_2(Omega) = alpha_2(3-cycles) + alpha_{2,cross}
where alpha_{2,cross} = # (5-cycle, 3-cycle) vertex-disjoint pairs.

Question: How significant is alpha_{2,cross}?
If it's always 0 or always small, the trace formula is still useful.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from tournament_fast import c3_from_score, c5_fast, alpha2_from_trace
from itertools import combinations
from math import comb
import random


def find_all_odd_cycles_with_vsets(T):
    """Find all directed odd cycles and return their vertex sets."""
    n = len(T)
    cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            v = list(verts)
            kk = len(v)
            dp = [[0]*kk for _ in range(1 << kk)]
            dp[1][0] = 1
            for mask in range(1, 1 << kk):
                for last in range(kk):
                    if dp[mask][last] == 0 or not (mask & (1 << last)):
                        continue
                    for nxt in range(1, kk):
                        if mask & (1 << nxt):
                            continue
                        if T[v[last]][v[nxt]]:
                            dp[mask | (1 << nxt)][nxt] += dp[mask][last]
            full = (1 << kk) - 1
            for last in range(1, kk):
                if T[v[last]][v[0]]:
                    cnt = dp[full][last]
                    for _ in range(cnt):
                        cycles.append((k, frozenset(verts)))
    return cycles


print("=" * 70)
print("alpha_2 at n=8: contribution from (5-cycle, 3-cycle) pairs")
print("=" * 70)

n = 8
m = n * (n - 1) // 2
random.seed(42)
sample = [random.randint(0, (1 << m) - 1) for _ in range(200)]

stats = {
    'alpha2_3only': [],
    'alpha2_cross': [],
    'alpha2_full': [],
    'fraction_cross': [],
}

for idx, bits in enumerate(sample):
    T = tournament_from_bits(n, bits)

    # Find all odd cycles
    all_cycles = find_all_odd_cycles_with_vsets(T)

    # Separate 3-cycles and 5-cycles (7-cycles exist too but 7+3=10 > 8)
    cycles_3 = [c for c in all_cycles if c[0] == 3]
    cycles_5 = [c for c in all_cycles if c[0] == 5]
    cycles_7 = [c for c in all_cycles if c[0] == 7]

    # alpha_2 from 3-cycles only
    a2_3 = 0
    for i in range(len(cycles_3)):
        for j in range(i+1, len(cycles_3)):
            if not (cycles_3[i][1] & cycles_3[j][1]):
                a2_3 += 1

    # alpha_2 cross: (3-cycle, 5-cycle) disjoint pairs
    a2_cross = 0
    for c3 in cycles_3:
        for c5 in cycles_5:
            if not (c3[1] & c5[1]):
                a2_cross += 1

    # Full alpha_2 (all disjoint pairs of odd cycles)
    m_cyc = len(all_cycles)
    a2_full = 0
    for i in range(m_cyc):
        for j in range(i+1, m_cyc):
            if not (all_cycles[i][1] & all_cycles[j][1]):
                a2_full += 1

    # Check: alpha_2_full = a2_3 + a2_cross (since 5+5>8 and 7+3>8 and 7+5>8)
    assert a2_full == a2_3 + a2_cross, \
        f"bits={bits}: full={a2_full}, 3only={a2_3}, cross={a2_cross}"

    # alpha_2 from trace formula (3-cycles only)
    a2_trace = alpha2_from_trace(T)
    assert a2_trace == a2_3, f"bits={bits}: trace={a2_trace}, direct={a2_3}"

    stats['alpha2_3only'].append(a2_3)
    stats['alpha2_cross'].append(a2_cross)
    stats['alpha2_full'].append(a2_full)
    if a2_full > 0:
        stats['fraction_cross'].append(a2_cross / a2_full)

    if idx < 5 or (a2_cross > 0 and idx < 50):
        print(f"  bits={bits}: a2_3={a2_3}, a2_cross={a2_cross}, a2_full={a2_full}, "
              f"c3={len(cycles_3)}, c5={len(cycles_5)}, c7={len(cycles_7)}")


print(f"\nSummary (200 tournaments at n=8):")
print(f"  alpha_2(3-cycles only): min={min(stats['alpha2_3only'])}, "
      f"max={max(stats['alpha2_3only'])}, avg={sum(stats['alpha2_3only'])/len(stats['alpha2_3only']):.1f}")
print(f"  alpha_2(cross 3x5): min={min(stats['alpha2_cross'])}, "
      f"max={max(stats['alpha2_cross'])}, avg={sum(stats['alpha2_cross'])/len(stats['alpha2_cross']):.1f}")
print(f"  alpha_2(full): min={min(stats['alpha2_full'])}, "
      f"max={max(stats['alpha2_full'])}, avg={sum(stats['alpha2_full'])/len(stats['alpha2_full']):.1f}")

nonzero_cross = sum(1 for x in stats['alpha2_cross'] if x > 0)
print(f"  Tournaments with nonzero cross contribution: {nonzero_cross}/200 ({nonzero_cross/2:.1f}%)")

if stats['fraction_cross']:
    avg_frac = sum(stats['fraction_cross']) / len(stats['fraction_cross'])
    max_frac = max(stats['fraction_cross'])
    print(f"  When nonzero: avg fraction = {avg_frac:.4f}, max = {max_frac:.4f}")


# ============================================================
# Impact on H(T): how much does ignoring cross terms affect H?
# ============================================================
print("\n" + "=" * 70)
print("H(T) ERROR from ignoring cross (5,3) pairs at n=8")
print("=" * 70)

from tournament_lib import hamiltonian_path_count

h_errors = []
for idx, bits in enumerate(sample[:100]):
    T = tournament_from_bits(n, bits)
    h_true = hamiltonian_path_count(T)

    all_cycles = find_all_odd_cycles_with_vsets(T)
    alpha_1 = len(all_cycles)
    alpha_2_full = 0
    m_cyc = len(all_cycles)
    for i in range(m_cyc):
        for j in range(i+1, m_cyc):
            if not (all_cycles[i][1] & all_cycles[j][1]):
                alpha_2_full += 1

    # Also need alpha_3
    alpha_3 = 0
    # At n=8: 3 disjoint 3-cycles need 9 > 8 vertices. So alpha_3 from
    # (3,3,3) = 0. What about (3,3,other)? 3+3+3=9>8, 3+3+5=11>8. All impossible.
    # So alpha_3 = 0 at n=8.

    h_ocf = 1 + 2*alpha_1 + 4*alpha_2_full + 8*alpha_3

    # Using 3-cycle-only alpha_2
    a2_3only = alpha2_from_trace(T)
    h_approx = 1 + 2*alpha_1 + 4*a2_3only

    error = h_true - h_approx
    h_errors.append(error)

    if abs(error) > 0 and idx < 10:
        print(f"  bits={bits}: H_true={h_true}, H_ocf={h_ocf}, H_approx={h_approx}, "
              f"error={error}")

    # Also verify OCF
    if h_true != h_ocf:
        print(f"  OCF MISMATCH at bits={bits}: true={h_true}, ocf={h_ocf}")

nonzero_errors = sum(1 for e in h_errors if e != 0)
print(f"\n  Tournaments with H error: {nonzero_errors}/100 ({nonzero_errors}%)")
if h_errors:
    max_err = max(abs(e) for e in h_errors)
    print(f"  Max |error|: {max_err}")
    avg_err = sum(abs(e) for e in h_errors) / len(h_errors)
    print(f"  Avg |error|: {avg_err:.2f}")


# ============================================================
# Can we compute alpha_{2,cross} from trace data?
# ============================================================
print("\n" + "=" * 70)
print("FORMULA for alpha_{2,cross}: (5-cycle, 3-cycle) disjoint pairs")
print("=" * 70)

# Observation: a (5-cycle, 3-cycle) pair at n=8 uses ALL 8 vertices.
# The 5-cycle occupies 5 vertices, the 3-cycle the remaining 3.
# So alpha_{2,cross} = sum over all C(8,5)=56 vertex 5-subsets:
#   (# directed 5-cycles on that subset) * (# directed 3-cycles on complement)

# And we can compute both using trace formulas!
# For the 5-subset S: c_5(T[S]) = tr(A_S^5)/5 where A_S is submatrix
# For the 3-complement S^c: c_3(T[S^c]) = 0 or 1 (only 1 triple)

# So: alpha_{2,cross} = sum over C(8,5) 5-subsets S:
#   c_5(T[S]) * c_3(T[S^c])
# where c_3 on 3 vertices is 0 (transitive) or 1 (cyclic).

# c_5(T[S]) via trace: O(5^3) = O(125) per subset
# c_3(T[S^c]): check one triple = O(1)
# Total: O(C(8,5) * 125) = O(7000) which is O(n^3) overall!

print("Testing: alpha_{2,cross} = sum over 5-subsets of c5(S)*c3(S^c)")

match = 0
for idx, bits in enumerate(sample[:100]):
    T = tournament_from_bits(n, bits)
    nn = len(T)

    # Direct computation
    all_cycles = find_all_odd_cycles_with_vsets(T)
    cycles_3 = [c for c in all_cycles if c[0] == 3]
    cycles_5 = [c for c in all_cycles if c[0] == 5]
    a2_cross_direct = 0
    for c3 in cycles_3:
        for c5 in cycles_5:
            if not (c3[1] & c5[1]):
                a2_cross_direct += 1

    # Formula: sum over 5-subsets
    a2_cross_formula = 0
    for subset in combinations(range(nn), 5):
        S = list(subset)
        Sc = [v for v in range(nn) if v not in subset]
        assert len(Sc) == 3

        # c_5 on S via trace
        T_sub = [[T[S[i]][S[j]] for j in range(5)] for i in range(5)]
        A2 = [[sum(T_sub[i][k]*T_sub[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A4 = [[sum(A2[i][k]*A2[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A5 = [[sum(A4[i][k]*T_sub[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        tr5 = sum(A5[i][i] for i in range(5))
        c5_S = tr5 // 5

        # c_3 on S^c (complement): just check if it's cyclic
        a, b, c = Sc
        c3_Sc = 0
        if T[a][b] and T[b][c] and T[c][a]:
            c3_Sc = 1
        elif T[a][c] and T[c][b] and T[b][a]:
            c3_Sc = 1

        a2_cross_formula += c5_S * c3_Sc

    if a2_cross_direct == a2_cross_formula:
        match += 1
    elif idx < 5:
        print(f"  MISMATCH bits={bits}: direct={a2_cross_direct}, formula={a2_cross_formula}")

print(f"  Match: {match}/100 ({'PASS' if match == 100 else 'FAIL'})")


print("\nDone.")
