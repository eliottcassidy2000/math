#!/usr/bin/env python3
"""
Connect W(r) coefficients to OCF: H = I(Omega(T), 2) = sum alpha_k * 2^k.

At n=5:
  W(r) = w_0 + w_2*r^2 + w_4*r^4  (odd coefficients vanish)
  H = W(1/2) = w_0 + w_2/4 + w_4/16

  w_4 = H (number of permutations = n! at r^{n-1}... wait, w_{n-1} = H)
  Actually w_4 = H * (1/2)^0 ... no.

Let me be more careful about the indexing.

W(r) = sum_{k=0}^{n-1} w_k * r^k

At r = 1/2: W(1/2) = sum_P prod(1/2 + s_e)
  prod(1/2 + s_e) = prod_e T(v_i, v_{i+1}) since (1/2 + s_e) = T[u,v] when s_e = T-1/2
  Wait: s_e = A[u,v] - 1/2, so r + s_e = 1/2 + A[u,v] - 1/2 = A[u,v].
  So W(1/2) = sum_P prod A[v_i,v_{i+1}] = number of DIRECTED Hamiltonian paths = H(T).

At r = -1/2: W(-1/2) = sum_P prod(-1/2 + s_e) = sum_P prod(A[v_i,v_{i+1}] - 1)
  = sum_P prod(-(1-A[v_i,v_{i+1}])) = sum_P (-1)^{n-1} prod A^{op}[v_i,v_{i+1}]
  = (-1)^{n-1} H(T^op) = (-1)^{n-1} H(T)  (since H(T) = H(T^op))

So W(-1/2) = -H (at even n-1) or +H (at odd n-1).
For n=5: n-1=4, (-1)^4 = 1, so W(-1/2) = H. But W has only even powers, so W(-1/2) = W(1/2) = H. Consistent!

Now: H = W(1/2) = w_0*(1/2)^0 + w_2*(1/2)^2 + w_4*(1/2)^4
Wait, W(r) = w_0 + w_1*r + w_2*r^2 + w_3*r^3 + w_4*r^4
With w_1 = w_3 = 0:
  H = w_0 + w_2/4 + w_4/16

Also: H = I(Omega,2) = 1 + alpha_1*2 + alpha_2*4

These are different decompositions of H!

kind-pasteur-2026-03-06-S25f
"""

from itertools import permutations, combinations
from math import factorial, comb

def tournament_all_n5():
    """Generate all non-isomorphic tournaments on 5 vertices (12 classes), returning one per class."""
    # Just generate all 2^10 = 1024 labeled tournaments
    all_T = []
    for bits in range(1024):
        A = [[0]*5 for _ in range(5)]
        idx = 0
        for i in range(5):
            for j in range(i+1, 5):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        all_T.append(A)
    return all_T

def count_kcycles(A, k):
    n = len(A)
    count = 0
    for verts in combinations(range(n), k):
        for p in permutations(verts):
            is_cycle = all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k))
            if is_cycle:
                count += 1
    return count // k

def conflict_graph(A, n):
    """Build conflict graph Omega(T): vertices = odd cycles, edges = shared vertex."""
    cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            for p in permutations(verts):
                # Normalize: start from smallest vertex, go in canonical direction
                if p[0] != min(p):
                    continue
                if k > 1 and p[1] > p[-1]:
                    continue
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    cycles.append(frozenset(verts))
    # Remove duplicates
    cycles = list(set(cycles))
    return cycles

def independence_poly(cycles):
    """Compute independence polynomial of conflict graph."""
    n_cycles = len(cycles)
    # Build adjacency (conflict = shared vertex)
    adj = [[False]*n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True

    # Enumerate independent sets
    alpha = {}  # alpha[k] = number of independent sets of size k
    for mask in range(1 << n_cycles):
        bits = [i for i in range(n_cycles) if (mask >> i) & 1]
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

print("=" * 70)
print("W(r) vs OCF DECOMPOSITION at n=5")
print("=" * 70)

all_T = tournament_all_n5()

# Group by isomorphism class (use score sequence + cycle counts as proxy)
seen = set()
representatives = []
for A in all_T:
    t3 = count_kcycles(A, 3)
    t5 = count_kcycles(A, 5)
    scores = tuple(sorted(sum(A[i]) for i in range(5)))
    key = (scores, t3, t5)
    if key not in seen:
        seen.add(key)
        representatives.append((A, key))

print(f"\nFound {len(representatives)} distinct (score,t3,t5) classes")

print(f"\n{'scores':<16} {'t3':>3} {'t5':>3} {'H':>4} | {'w0':>6} {'w2':>6} {'w4':>6} | {'a0':>3} {'a1':>3} {'a2':>3} | {'OCF':>4} {'W-check':>7}")

for A, key in sorted(representatives, key=lambda x: x[1]):
    scores, t3, t5 = key
    n = 5

    # W coefficients
    wc = W_coefficients(A)
    w0, w2, w4 = wc[0], wc[2], wc[4]

    # H from permutation count
    H = int(round(w4))  # w_{n-1} = H? No... let me check
    # Actually w_4 is coefficient of r^4 = r^{n-1}
    # W(r) = sum_P prod(r + s_e) where each factor is degree 1
    # So W has degree n-1 = 4
    # Coefficient of r^4 = sum_P 1 = n! (total permutations)
    # Wait, that's wrong. W sums over ALL permutations, not just Ham paths.
    # Each permutation contributes prod(r + s_e) regardless of arc directions.
    # At r=1/2: prod(1/2 + s_e) = prod(A[u,v]) = 1 if directed path, 0 otherwise.
    # So w_4 = sum_P 1 = n! = 120

    H_actual = sum(1 for p in permutations(range(5))
                   if all(A[p[i]][p[i+1]] == 1 for i in range(4)))

    # OCF
    cycles = conflict_graph(A, 5)
    alpha = independence_poly(cycles)
    a0 = alpha.get(0, 1)  # always 1
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)

    H_ocf = a0 + a1*2 + a2*4

    # W check: H = w0 + w2/4 + w4/16
    H_W = w0 + w2/4 + w4/16

    print(f"{str(scores):<16} {t3:3d} {t5:3d} {H_actual:4d} | {w0:6.1f} {w2:6.1f} {w4:6.1f} | {a0:3d} {a1:3d} {a2:3d} | {H_ocf:4d} {H_W:7.1f}")

print(f"\n{'':>39} Note: w4 = n! = 120 always (sum over all perms)")

# Verify w0 = -t3 + 2*t5 + 1
print(f"\nVerify w0 = -t3 + 2*t5 + 1:")
for A, key in sorted(representatives, key=lambda x: x[1]):
    scores, t3, t5 = key
    wc = W_coefficients(A)
    w0 = wc[0]
    predicted = -t3 + 2*t5 + 1
    print(f"  {str(scores):<16} t3={t3}, t5={t5}: w0={w0:.1f}, predicted={predicted}, match={'Y' if abs(w0-predicted)<0.01 else 'N'}")

# Now the key connection:
# H = w0 + w2/4 + w4/16
# H = 1 + 2*a1 + 4*a2
#
# So: w0 + w2/4 + w4/16 = 1 + 2*a1 + 4*a2
# With w4 = 120, w2 = 12*t3 - 30, w0 = -t3 + 2*t5 + 1:
# (-t3 + 2*t5 + 1) + (12*t3 - 30)/4 + 120/16 = 1 + 2*a1 + 4*a2
# -t3 + 2*t5 + 1 + 3*t3 - 7.5 + 7.5 = 1 + 2*a1 + 4*a2
# 2*t3 + 2*t5 + 1 = 1 + 2*a1 + 4*a2
# 2*t3 + 2*t5 = 2*a1 + 4*a2
# t3 + t5 = a1 + 2*a2

print(f"\nDERIVED IDENTITY: t3 + t5 = a1 + 2*a2")
print(f"Where a1 = # odd cycles, a2 = # independent pairs of odd cycles")
print(f"\nVerification:")
for A, key in sorted(representatives, key=lambda x: x[1]):
    scores, t3, t5 = key
    cycles = conflict_graph(A, 5)
    alpha = independence_poly(cycles)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    lhs = t3 + t5
    rhs = a1 + 2*a2
    print(f"  {str(scores):<16} t3+t5={lhs}, a1+2*a2={rhs}, match={'Y' if lhs==rhs else 'N'}")

print(f"\nBut wait: a1 = total number of odd cycles = t3 + t5 (at n=5, only 3-cycles and 5-cycles)")
print(f"So: t3 + t5 = (t3 + t5) + 2*a2 => a2 = 0?!")
print(f"Check:")
for A, key in sorted(representatives, key=lambda x: x[1]):
    scores, t3, t5 = key
    cycles = conflict_graph(A, 5)
    alpha = independence_poly(cycles)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    print(f"  {str(scores):<16} a1={a1}, t3+t5={t3+t5}, a2={a2}")

print("\n" + "=" * 70)
print("DONE")
