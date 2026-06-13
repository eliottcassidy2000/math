import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
alpha_structure_n9.py
kind-pasteur-2026-03-07-S39b

Analyze the independence polynomial structure at n=9 where NEW phenomena appear:
1. alpha_3 can be nonzero (3+3+3 = 9 vertices exactly)
2. alpha_2 has cross terms: (5-cycle, 3-cycle) with 5+3=8 < 9 spare vertex
3. c_7 contributes to alpha_1 (7-cycles exist)
4. c_9 contributes to alpha_1 (full Hamiltonian cycles)

Goal: Compute exact I(Omega(T), x) = 1 + alpha_1*x + alpha_2*x^2 + alpha_3*x^3
for sampled n=9 tournaments and verify H(T) = I(Omega(T), 2).
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from tournament_fast import c3_from_score, c5_fast, alpha2_from_trace
from itertools import combinations
from math import comb
import random
import time


def count_directed_cycles_subset(T, verts):
    """Count directed Hamiltonian cycles on the subtournament T[verts].
    Uses subset DP. Returns count of directed cycles (each direction counted separately).
    """
    k = len(verts)
    if k < 3:
        return 0
    V = list(verts)

    # Subset DP: count Hamiltonian paths starting at V[0]
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for last in range(k):
            if dp[mask][last] == 0 or not (mask & (1 << last)):
                continue
            for nxt in range(1, k):
                if mask & (1 << nxt):
                    continue
                if T[V[last]][V[nxt]]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]

    full = (1 << k) - 1
    cycles = 0
    for last in range(1, k):
        if T[V[last]][V[0]]:
            cycles += dp[full][last]
    return cycles


def c7_count(T):
    """Count total directed 7-cycles in T."""
    n = len(T)
    if n < 7:
        return 0
    total = 0
    for subset in combinations(range(n), 7):
        total += count_directed_cycles_subset(T, subset)
    return total


def c9_count(T):
    """Count directed 9-cycles (Hamiltonian cycles) in T."""
    n = len(T)
    if n != 9:
        return count_directed_cycles_subset(T, range(n)) if n >= 3 else 0
    return count_directed_cycles_subset(T, range(9))


def alpha2_full(T):
    """Compute full alpha_2 at n=9: disjoint pairs of odd cycles.

    Possible pair types at n=9:
    - (3, 3): 3+3=6 <= 9, spare 3 vertices
    - (3, 5): 3+5=8 <= 9, spare 1 vertex
    - (3, 7): 3+7=10 > 9, IMPOSSIBLE
    - (5, 5): 5+5=10 > 9, IMPOSSIBLE

    So alpha_2 = alpha_2(3,3) + alpha_2(3,5)
    """
    n = len(T)

    # alpha_2(3,3): vertex-disjoint 3-cycle pairs
    a2_33 = alpha2_from_trace(T)  # Uses THM-097

    # alpha_2(3,5): vertex-disjoint (3-cycle, 5-cycle) pairs
    # For each 5-subset S: c_5(T[S]) * c_3(T[V\S]) summed over all
    # But at n=9, V\S has 4 vertices, so c_3(V\S) = c_3(T[4-vertex subtournament])
    a2_35 = 0
    for subset in combinations(range(n), 5):
        S = list(subset)
        Sc = [v for v in range(n) if v not in subset]

        # c_5 on S via trace
        Ts = [[T[S[i]][S[j]] for j in range(5)] for i in range(5)]
        A2 = [[sum(Ts[i][k]*Ts[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A4 = [[sum(A2[i][k]*A2[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A5 = [[sum(A4[i][k]*Ts[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        c5_S = sum(A5[i][i] for i in range(5)) // 5

        if c5_S == 0:
            continue

        # c_3 on S^c (4 vertices) via Moon's formula
        T_sc = [[T[Sc[i]][Sc[j]] for j in range(len(Sc))] for i in range(len(Sc))]
        scores_sc = [sum(row) for row in T_sc]
        c3_Sc = comb(len(Sc), 3) - sum(comb(s, 2) for s in scores_sc)

        a2_35 += c5_S * c3_Sc

    return a2_33, a2_35


def alpha3_full(T):
    """Compute alpha_3 at n=9: vertex-disjoint triples of odd cycles.

    Only possible type at n=9: (3, 3, 3) using all 9 vertices.
    3+3+3 = 9 = n, so we need three vertex-disjoint 3-cycles covering all vertices.

    alpha_3 = number of ordered triples / 6 (unordered).
    But actually alpha_3 counts independent sets of size 3 in Omega(T).
    """
    n = len(T)
    if n < 9:
        return 0

    # Find all 3-cycles
    cycles_3 = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a]):
            cycles_3.append(frozenset(triple))

    # Count triples of mutually vertex-disjoint 3-cycles
    count = 0
    nc = len(cycles_3)
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles_3[i] & cycles_3[j]:
                continue
            for k in range(j+1, nc):
                if (not (cycles_3[i] & cycles_3[k])) and (not (cycles_3[j] & cycles_3[k])):
                    count += 1

    return count


# ============================================================
# Main verification
# ============================================================
print("=" * 70)
print("ALPHA STRUCTURE AT n=9: I(Omega(T), x) decomposition")
print("=" * 70)

n = 9
m = n * (n - 1) // 2  # 36
random.seed(42)
sample = [random.randint(0, (1 << m) - 1) for _ in range(50)]

results = []
t0 = time.time()

for idx, bits in enumerate(sample):
    T = tournament_from_bits(n, bits)

    # Compute all cycle counts
    c3 = c3_from_score(T)
    c5 = c5_fast(T)
    c7 = c7_count(T)
    c9 = c9_count(T)

    # alpha_1 = total directed odd cycles
    alpha_1 = c3 + c5 + c7 + c9

    # alpha_2 decomposition
    a2_33, a2_35 = alpha2_full(T)
    alpha_2 = a2_33 + a2_35

    # alpha_3
    alpha_3 = alpha3_full(T)

    # OCF prediction
    h_ocf = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

    # True H
    h_true = hamiltonian_path_count(T)

    match = (h_ocf == h_true)

    results.append({
        'bits': bits, 'c3': c3, 'c5': c5, 'c7': c7, 'c9': c9,
        'alpha_1': alpha_1, 'a2_33': a2_33, 'a2_35': a2_35,
        'alpha_2': alpha_2, 'alpha_3': alpha_3,
        'h_ocf': h_ocf, 'h_true': h_true, 'match': match
    })

    if idx < 10 or not match:
        status = "OK" if match else "MISMATCH"
        print(f"  [{status}] bits={bits}: c3={c3} c5={c5} c7={c7} c9={c9}")
        print(f"    alpha_1={alpha_1}, alpha_2={a2_33}+{a2_35}={alpha_2}, alpha_3={alpha_3}")
        print(f"    H_ocf={h_ocf}, H_true={h_true}")

t1 = time.time()

matches = sum(r['match'] for r in results)
print(f"\nn=9 ({len(results)} sampled): OCF verification: {matches}/{len(results)} "
      f"({'PASS' if matches == len(results) else 'FAIL'})")
print(f"Time: {t1-t0:.1f}s ({(t1-t0)/len(results):.2f}s per tournament)")

# Statistics
print(f"\nStatistics:")
print(f"  c3: {min(r['c3'] for r in results)}-{max(r['c3'] for r in results)}, "
      f"avg={sum(r['c3'] for r in results)/len(results):.1f}")
print(f"  c5: {min(r['c5'] for r in results)}-{max(r['c5'] for r in results)}, "
      f"avg={sum(r['c5'] for r in results)/len(results):.1f}")
print(f"  c7: {min(r['c7'] for r in results)}-{max(r['c7'] for r in results)}, "
      f"avg={sum(r['c7'] for r in results)/len(results):.1f}")
print(f"  c9: {min(r['c9'] for r in results)}-{max(r['c9'] for r in results)}, "
      f"avg={sum(r['c9'] for r in results)/len(results):.1f}")

print(f"\n  alpha_1: {min(r['alpha_1'] for r in results)}-{max(r['alpha_1'] for r in results)}, "
      f"avg={sum(r['alpha_1'] for r in results)/len(results):.1f}")
print(f"  alpha_2(3,3): avg={sum(r['a2_33'] for r in results)/len(results):.1f}")
print(f"  alpha_2(3,5): avg={sum(r['a2_35'] for r in results)/len(results):.1f}")
print(f"  alpha_2: avg={sum(r['alpha_2'] for r in results)/len(results):.1f}")

nonzero_a3 = sum(1 for r in results if r['alpha_3'] > 0)
print(f"  alpha_3: nonzero in {nonzero_a3}/{len(results)} "
      f"({100*nonzero_a3/len(results):.0f}%)")
if nonzero_a3 > 0:
    a3_vals = [r['alpha_3'] for r in results if r['alpha_3'] > 0]
    print(f"    When nonzero: min={min(a3_vals)}, max={max(a3_vals)}, "
          f"avg={sum(a3_vals)/len(a3_vals):.1f}")

# H contribution analysis
print(f"\nH contribution analysis:")
h_vals = [r['h_true'] for r in results]
a1_contrib = [2*r['alpha_1'] for r in results]
a2_contrib = [4*r['alpha_2'] for r in results]
a3_contrib = [8*r['alpha_3'] for r in results]

print(f"  H: min={min(h_vals)}, max={max(h_vals)}, avg={sum(h_vals)/len(h_vals):.0f}")
print(f"  2*alpha_1 contribution: avg={sum(a1_contrib)/len(a1_contrib):.0f} "
      f"({100*sum(a1_contrib)/sum(h_vals):.1f}% of H)")
print(f"  4*alpha_2 contribution: avg={sum(a2_contrib)/len(a2_contrib):.0f} "
      f"({100*sum(a2_contrib)/sum(h_vals):.1f}% of H)")
print(f"  8*alpha_3 contribution: avg={sum(a3_contrib)/len(a3_contrib):.0f} "
      f"({100*sum(a3_contrib)/sum(h_vals):.1f}% of H)")

# Newton inequality check: does I(Omega, x) have all real roots at n=9?
print(f"\n" + "=" * 70)
print(f"NEWTON INEQUALITY CHECK at n=9")
print(f"=" * 70)

newton_fail = 0
disc_fail = 0
for r in results:
    a1, a2, a3 = r['alpha_1'], r['alpha_2'], r['alpha_3']

    # For degree-3 polynomial 1 + a1*x + a2*x^2 + a3*x^3:
    # Newton inequalities: a_k^2 >= a_{k-1} * a_{k+1} * (k+1)/(k) ... actually
    # For the sequence (1, a1, a2, a3):
    # Newton at k=1: a1^2 >= 2*1*a2 = 2*a2
    # Newton at k=2: a2^2 >= (3/2)*a1*a3

    if a2 > 0:
        newton1 = a1**2 - 2*a2
        if newton1 < 0:
            newton_fail += 1
            if newton_fail <= 3:
                print(f"  Newton-1 FAIL: bits={r['bits']}, a1={a1}, a2={a2}, a3={a3}")

    if a3 > 0 and a1 > 0:
        # For the sequence to have all real roots, Newton inequalities must hold
        # a_k^2 >= a_{k-1}*a_{k+1} * C(d,k)^2 / (C(d,k-1)*C(d,k+1))
        # where d = degree = 3
        # At k=1: a1^2 >= (C(3,1)^2/(C(3,0)*C(3,2))) * 1 * a2 = (9/3)*a2 = 3*a2
        # At k=2: a2^2 >= (C(3,2)^2/(C(3,1)*C(3,3))) * a1*a3 = (9/3)*a1*a3 = 3*a1*a3
        newton2 = a2**2 - 3*a1*a3
        if newton2 < 0 and newton_fail <= 5:
            print(f"  Newton-2 FAIL: bits={r['bits']}, a1={a1}, a2={a2}, a3={a3}, "
                  f"margin={newton2}")

    # Discriminant of cubic I(Omega, x) = a3*x^3 + a2*x^2 + a1*x + 1
    # For all real roots, discriminant >= 0
    if a3 > 0:
        D = 18*a3*a2*a1*1 - 4*a2**3*1 + a2**2*a1**2 - 4*a3*a1**3 - 27*a3**2*1**2
        if D < 0:
            disc_fail += 1
            if disc_fail <= 3:
                print(f"  Disc < 0: bits={r['bits']}, a1={a1}, a2={a2}, a3={a3}, D={D}")

total_cubic = sum(1 for r in results if r['alpha_3'] > 0)
print(f"\n  Tournaments with alpha_3 > 0 (cubic I.P.): {total_cubic}/{len(results)}")
print(f"  Newton-1 failures: {newton_fail}")
print(f"  Discriminant < 0 (complex roots): {disc_fail}")

print("\nDone.")
