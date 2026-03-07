#!/usr/bin/env python3
"""
Connect W(r) coefficients to OCF at n=5.

Fix: properly count odd cycles and build conflict graph.

kind-pasteur-2026-03-06-S25f
"""

from itertools import permutations, combinations
from math import factorial, comb

def count_directed_cycles(A, k):
    """Count directed k-cycles."""
    n = len(A)
    count = 0
    for verts in combinations(range(n), k):
        for p in permutations(verts):
            if p[0] != min(p):
                continue
            if k > 2 and p[1] > p[-1]:
                continue
            if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                count += 1
    return count

def all_odd_cycles(A):
    """Return all directed odd cycles as frozensets of their vertex sets."""
    n = len(A)
    cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            # Check if there exists ANY directed k-cycle on these vertices
            found = False
            for p in permutations(verts):
                if p[0] != min(p):
                    continue
                if k > 2 and p[1] > p[-1]:
                    continue
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    cycles.append(frozenset(verts))
                    found = True
            # Note: multiple cycles can exist on same vertex set
    return cycles

def independence_poly_from_cycles(cycles):
    """Compute independence polynomial of conflict graph."""
    nc = len(cycles)
    # Build adjacency (conflict = shared vertex)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:  # shared vertex
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size
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
print(f"W(r) vs OCF at n={n}")
print("=" * 70)

# Generate all labeled tournaments
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

    cycles = all_odd_cycles(A)
    alpha = independence_poly_from_cycles(cycles)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)

    H = sum(1 for p in permutations(range(5))
            if all(A[p[i]][p[i+1]] == 1 for i in range(4)))

    H_ocf = 1 + a1*2 + a2*4
    assert H == H_ocf, f"OCF failed: H={H}, OCF={H_ocf}, t3={t3}, t5={t5}"

    wc = W_coefficients(A)
    w0, w2, w4 = wc[0], wc[2], wc[4]

    key = (scores, t3, t5)
    if key not in results:
        results[key] = {'H': H, 'w0': w0, 'w2': w2, 'w4': w4,
                        'a1': a1, 'a2': a2, 'count': 1,
                        'cycles': len(cycles)}
    else:
        results[key]['count'] += 1
        # Verify all same
        assert results[key]['H'] == H

print(f"\n{'scores':<16} {'t3':>3} {'t5':>3} {'H':>4} | {'w0':>6} {'w2':>6} {'w4':>6} | {'#cyc':>4} {'a1':>3} {'a2':>3} | {'cnt':>4}")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    print(f"{str(scores):<16} {t3:3d} {t5:3d} {r['H']:4d} | {r['w0']:6.1f} {r['w2']:6.1f} {r['w4']:6.1f} | {r['cycles']:4d} {r['a1']:3d} {r['a2']:3d} | {r['count']:4d}")

# Check identity: H = 1 + 2*a1 + 4*a2
print(f"\nAll OCF checks passed!" if True else "FAILED")

# Check: w0 = -t3 + 2*t5 + 1
print(f"\nVerify w0 = -t3 + 2*t5 + 1:")
all_ok = True
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    predicted = -t3 + 2*t5 + 1
    ok = abs(r['w0'] - predicted) < 0.01
    if not ok:
        all_ok = False
        print(f"  FAIL: {scores}, t3={t3}, t5={t5}: w0={r['w0']}, predicted={predicted}")
print(f"  All match: {all_ok}")

# Derive relationship between W and OCF coefficients
print(f"\n{'='*70}")
print("RELATIONSHIP BETWEEN W(r) AND OCF")
print(f"{'='*70}")
print(f"""
  H = W(1/2) = w0 + w2/4 + w4/16
  H = I(Omega,2) = 1 + 2*a1 + 4*a2

  w4 = n! = 120 (always)
  w2 = 12*t3 - 30 (proven, opus S27)
  w0 = -t3 + 2*t5 + 1 (verified above)

  So: (-t3 + 2*t5 + 1) + (12*t3 - 30)/4 + 120/16 = 1 + 2*a1 + 4*a2
      -t3 + 2*t5 + 1 + 3*t3 - 7.5 + 7.5 = 1 + 2*a1 + 4*a2
      2*t3 + 2*t5 + 1 = 1 + 2*a1 + 4*a2
      2*(t3 + t5) = 2*a1 + 4*a2
      t3 + t5 = a1 + 2*a2
""")

print("Verify t3 + t5 = a1 + 2*a2:")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    lhs = t3 + t5
    rhs = r['a1'] + 2*r['a2']
    ok = lhs == rhs
    print(f"  {str(scores):<16} t3+t5={lhs}, a1+2*a2={rhs}, match={ok}")

print(f"""
Interpretation:
  a1 = total number of odd cycles (= t3 + t5 at n=5)
  a2 = number of vertex-DISJOINT pairs of odd cycles

  So: t3 + t5 = (t3 + t5) + 2*a2
  => a2 = 0 always at n=5!

  This makes sense: at n=5 with 5 vertices, two vertex-disjoint
  odd cycles would need 3+3 = 6 vertices. Impossible!
  (A 3-cycle + 5-cycle shares at least 1 vertex too, since 3+5=8 > 5)

  So a2 = 0 at n=5. The OCF simplifies to H = 1 + 2*(t3 + t5).
""")

# Verify H = 1 + 2*(t3+t5)
print("Verify H = 1 + 2*(t3+t5):")
for key in sorted(results.keys()):
    scores, t3, t5 = key
    r = results[key]
    predicted = 1 + 2*(t3 + t5)
    ok = r['H'] == predicted
    print(f"  {str(scores):<16} H={r['H']}, 1+2*(t3+t5)={predicted}, match={ok}")

print("\nDONE")
