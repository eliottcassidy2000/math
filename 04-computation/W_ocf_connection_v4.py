#!/usr/bin/env python3
"""
W(r) vs OCF at n=5. Fixed: each DIRECTED odd cycle is a separate vertex in Omega.

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

def all_directed_odd_cycles(A):
    """Return all directed odd cycles as (vertex_set, canonical_ordering) pairs."""
    n = len(A)
    cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            for p in permutations(verts):
                # Canonical: start at min vertex, go in direction where second vertex < last
                if p[0] != min(p):
                    continue
                if k > 2 and p[1] > p[-1]:
                    continue
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    cycles.append((frozenset(verts), tuple(p)))
    return cycles

def independence_poly(cycles):
    nc = len(cycles)
    # Two cycles conflict if they share a vertex
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][0] & cycles[j][0]:  # vertex sets intersect
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
print(f"W(r) vs OCF at n={n} (v4 - per directed cycle)")
print("=" * 70)

results = {}
ocf_ok = 0
total = 0
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

    cycles = all_directed_odd_cycles(A)
    alpha = independence_poly(cycles)

    H = sum(1 for p in permutations(range(5))
            if all(A[p[i]][p[i+1]] == 1 for i in range(4)))

    H_ocf = sum(v * 2**k for k, v in alpha.items())
    total += 1
    if H == H_ocf:
        ocf_ok += 1

    wc = W_coefficients(A)
    w0, w2, w4 = wc[0], wc[2], wc[4]

    key = (scores, t3, t5)
    if key not in results:
        results[key] = {'H': H, 'H_ocf': H_ocf, 'w0': w0, 'w2': w2, 'w4': w4,
                        'alpha': dict(alpha), 'count': 1, 'ncycles': len(cycles)}
    else:
        results[key]['count'] += 1

print(f"\nOCF check: {ocf_ok}/{total} pass")
print(f"\n{'scores':<16} {'t3':>3} {'t5':>3} {'H':>4} {'Hocf':>4} | {'w0':>6} {'w2':>6} {'w4':>6} | {'#cyc':>4} {'alpha':>20}")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    alpha_str = str({k: v for k, v in sorted(r['alpha'].items()) if k > 0})
    match = 'Y' if r['H'] == r['H_ocf'] else 'N'
    print(f"{str(scores):<16} {t3:3d} {t5:3d} {r['H']:4d} {r['H_ocf']:4d}{match} | {r['w0']:6.1f} {r['w2']:6.1f} {r['w4']:6.1f} | {r['ncycles']:4d} {alpha_str}")

# Summary relationships
print(f"\nKEY: At n=5, directed cycles =")
print(f"  Each 3-vertex set supports exactly 1 directed 3-cycle (in a tournament)")
print(f"  Each 5-vertex set can support 0, 1, 2, or 3 directed 5-cycles")
print(f"  So #odd_cycles = t3 + (5-cycle count per vertex set summed)")
print(f"  And a2 can be nonzero when two vertex-disjoint 3-cycles exist... ")
print(f"  Wait: at n=5, two disjoint 3-cycles need 6 vertices. Impossible!")
print(f"  But two disjoint cycles where one is a 3-cycle on {a,b,c} and")
print(f"  a 5-cycle on {a,b,c,d,e} overlap! So a2=0 always at n=5.")
print(f"  Actually: two disjoint 3-cycles need 6 > 5. Two disjoint cycles of")
print(f"  ANY odd size need 3+3=6 > 5. So a2 = 0 always at n=5.")
print(f"  => H = 1 + 2*a1 where a1 = total number of directed odd cycles")

# Check H = 1 + 2*a1
print(f"\nH = 1 + 2*a1 check:")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    a1 = r['alpha'].get(1, 0)
    predicted = 1 + 2*a1
    ok = r['H'] == predicted
    print(f"  {str(scores):<16} H={r['H']:4d}, a1={a1:3d}, 1+2*a1={predicted:4d}, match={ok}")

# Check a1 = t3 + t5
print(f"\na1 = t3 + t5 check:")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    a1 = r['alpha'].get(1, 0)
    total_cycles = t3 + t5
    print(f"  {str(scores):<16} a1={a1:3d}, t3+t5={total_cycles:3d}, #cycles={r['ncycles']:3d}, match={a1 == r['ncycles']}")

# The number of directed odd cycles should be ncycles
# a1 = ncycles (each cycle is an independent set of size 1)
# So H = 1 + 2*ncycles

print(f"\nH = 1 + 2*(#directed_odd_cycles) at n=5:")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    predicted = 1 + 2*r['ncycles']
    ok = r['H'] == predicted
    print(f"  {str(scores):<16} H={r['H']:4d}, ncycles={r['ncycles']:3d}, 1+2*nc={predicted:4d}, match={ok}")

print(f"\nFINAL IDENTITY: At n=5, H(T) = 1 + 2*(total directed odd cycles)")
print("DONE")
