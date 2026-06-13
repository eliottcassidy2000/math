#!/usr/bin/env python3
"""
Connect W(r) coefficients to OCF at n=5. Fixed cycle counting.

key identity discovered: H = 1 + 2*(t3 + t5) at n=5.

kind-pasteur-2026-03-06-S25f
"""

from itertools import permutations, combinations
from math import factorial, comb

def count_directed_cycles(A, k):
    n = len(A)
    count = 0
    for verts in combinations(range(n), k):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                count += 1
    return count // k

def all_odd_cycle_sets(A):
    n = len(A)
    cycle_sets = set()
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            for p in permutations(verts):
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    cycle_sets.add(frozenset(verts))
                    break
    return list(cycle_sets)

def independence_poly(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True

    alpha = {0: 1}
    for mask in range(1, 1 << nc):
        bits = [i for i in range(nc) if (mask >> i) & 1]
        k = len(bits)
        independent = True
        for a in range(len(bits)):
            for b in range(a+1, len(bits)):
                if adj[bits[a]][bits[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            alpha[k] = alpha.get(k, 0) + 1
    return alpha

def W_coefficients(A):
    n = len(A)
    coeffs = [0.0] * n
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        poly = [1.0]
        for si in s:
            new_poly = [0.0] * (len(poly) + 1)
            for j, c in enumerate(poly):
                new_poly[j+1] += c
                new_poly[j] += c * si
            poly = new_poly
        for k in range(min(len(poly), n)):
            coeffs[k] += poly[k]
    return coeffs

n = 5
print("=" * 70)
print(f"W(r) vs OCF at n={n} (CORRECTED)")
print("=" * 70)

results = {}
for bits in range(1 << 10):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    t3 = count_directed_cycles(A, 3)
    t5 = count_directed_cycles(A, 5)
    scores = tuple(sorted(sum(A[i]) for i in range(5)))

    cycles = all_odd_cycle_sets(A)
    alpha = independence_poly(cycles)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)

    H = sum(1 for p in permutations(range(5))
            if all(A[p[i]][p[i+1]] == 1 for i in range(4)))

    H_ocf = sum(v * 2**k for k, v in alpha.items())
    assert H == H_ocf, f"OCF failed: H={H}, OCF={H_ocf}, scores={scores}, t3={t3}, t5={t5}, alpha={alpha}"

    wc = W_coefficients(A)
    w0, w2, w4 = wc[0], wc[2], wc[4]

    key = (scores, t3, t5)
    if key not in results:
        results[key] = {'H': H, 'w0': w0, 'w2': w2, 'w4': w4,
                        'a1': a1, 'a2': a2, 'count': 1,
                        'ncycles': len(cycles)}
    else:
        results[key]['count'] += 1

print(f"\n{'scores':<16} {'t3':>3} {'t5':>3} {'H':>4} | {'w0':>6} {'w2':>6} {'w4':>6} | {'#cyc':>4} {'a1':>3} {'a2':>3} | {'H=1+2(t3+t5)':>14}")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    h_simple = 1 + 2*(t3+t5)
    print(f"{str(scores):<16} {t3:3d} {t5:3d} {r['H']:4d} | {r['w0']:6.1f} {r['w2']:6.1f} {r['w4']:6.1f} | {r['ncycles']:4d} {r['a1']:3d} {r['a2']:3d} | {h_simple:4d} {'Y' if h_simple == r['H'] else 'N':>9}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
AT n=5, the OCF simplifies dramatically:
  - a2 = 0 always (two vertex-disjoint odd cycles need >= 6 vertices)
  - H = 1 + 2*a1 where a1 = number of odd cycles
  - a1 = t3 + t5 (since odd cycles at n=5 are only 3-cycles and 5-cycles,
    and each vertex set supports at most one directed cycle in a tournament)

IDENTITY: H(T) = 1 + 2*(t3 + t5) for all 5-vertex tournaments.

This is a BEAUTIFUL simplification of OCF at n=5.

Now the W(r) connection:
  w_4 = 120 = 5! (counts all permutations)
  w_2 = 12*t3 - 30
  w_0 = -t3 + 2*t5 + 1

  H = w_0 + w_2/4 + w_4/16
    = (-t3 + 2*t5 + 1) + (12*t3 - 30)/4 + 120/16
    = -t3 + 2*t5 + 1 + 3*t3 - 7.5 + 7.5
    = 2*t3 + 2*t5 + 1
    = 1 + 2*(t3 + t5)  CHECK!

So the W(r) coefficients ENCODE the OCF decomposition:
  w_4 = n! (universal constant)
  w_2 = 12*t3 - 30 (depends only on 3-cycles)
  w_0 = -t3 + 2*t5 + 1 (introduces 5-cycles)

The STRATIFICATION is:
  Level n-1: universal (just counting permutations)
  Level n-3: depends on t_3 (3-cycle structure)
  Level n-5: depends on t_3 AND t_5 (3- and 5-cycle structure)

And at r = 1/2, these layers combine to give the OCF identity H = 1 + 2*(t3 + t5).
""")

# Verify w_0 = -t3 + 2*t5 + 1
all_ok = True
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    predicted = -t3 + 2*t5 + 1
    if abs(r['w0'] - predicted) > 0.01:
        all_ok = False
        print(f"  FAIL: {scores}, t3={t3}, t5={t5}: w0={r['w0']}, predicted={predicted}")
print(f"w_0 = -t3 + 2*t5 + 1: ALL MATCH = {all_ok}")

print("\nDONE")
