#!/usr/bin/env python3
"""
What happens to W(r) under GS tiling flip?

THM-060: t3 parity flips for GS tilings at odd n.
THM-059: W(r) = C_0(r) + C_t3(r)*t3 + C_t5(r)*t5 + ...
THM-058: w_{n-3} = (n-2)! * [2*t3 - C(n,3)/2]

Question: does the FULL W(r) change in a predictable way under flip?
If we know W(r) for tiling T and W(r) for flip(T), what's the relationship?

Also: the Eulerian polynomial F_f(r) has the property F_f(1/2) = 1.
Does it have a symmetry at r = -1/2? This would be relevant since
flip reverses edge weights: s_e -> -s_e for non-backbone edges,
which maps r -> -r (heuristically).

Actually: backbone edges contribute (r + 1/2) always. Non-backbone edges
contribute (r + s_e) where s_e = +1/2 or -1/2. Flip changes s_e -> -s_e
for non-backbone edges, so those factors become (r - s_e). But backbone
factors stay as (r + 1/2).

So W_flip(r) is NOT simply W(-r). It's more subtle.

Let's compute W(r) and W_flip(r) for GS tilings at n=5 and look for patterns.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations
from fractions import Fraction
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def compute_W_poly(A, n):
    """Compute W(r) = sum over Ham paths P of prod(r + s_e).
    Returns dict: power -> coefficient (exact Fraction)."""
    W = defaultdict(lambda: Fraction(0))

    def ham_paths(mask, last):
        if mask == (1 << n) - 1:
            yield [last]
            return
        for v in range(n):
            if not (mask & (1 << v)) and A[last][v]:
                for path in ham_paths(mask | (1 << v), v):
                    yield [last] + path
            elif not (mask & (1 << v)) and A[v][last]:
                # v beats last, but we go last -> v anyway (this is wrong direction)
                pass

    # Actually, W(r) sums over ALL orderings, not just paths following edges
    # W(r) = sum over all permutations P of prod_{i=0}^{n-2} (r + A[P[i]][P[i+1]] - 1/2)

    for perm in permutations(range(n)):
        product = Fraction(1)
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            # (r + s) contributes to the polynomial
            # We need to track powers of r
            pass

    # Better: use polynomial multiplication
    # Start with constant 1, multiply by (r + s_i) for each step

    result = defaultdict(lambda: Fraction(0))

    for perm in permutations(range(n)):
        # Product of (r + s_i) for i=0..n-2
        # s_i = A[perm[i]][perm[i+1]] - 1/2
        poly = {0: Fraction(1)}  # polynomial: power -> coeff

        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            # Multiply by (r + s)
            new_poly = {}
            for power, coeff in poly.items():
                # coeff * r -> power+1
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                # coeff * s -> power
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly

        for power, coeff in poly.items():
            result[power] += coeff

    return dict(result)

def poly_str(poly, var='r'):
    """Pretty print polynomial."""
    terms = []
    for power in sorted(poly.keys(), reverse=True):
        c = poly[power]
        if c == 0:
            continue
        if power == 0:
            terms.append(str(c))
        elif power == 1:
            terms.append(f"{c}*{var}")
        else:
            terms.append(f"{c}*{var}^{power}")
    return " + ".join(terms) if terms else "0"

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs, fixed, seen = [], [], set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen: continue
        ti, tj = n-1-j, n-1-i
        if ti > tj: ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx); seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx)); seen.add(idx); seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

from itertools import combinations
def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

print(f"W(r) UNDER GS FLIP at n={n}")
print(f"{'='*70}")
print(f"  {len(gs_tilings)} GS tilings, m={m}")

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)

    if idx > flip_idx:
        continue  # Already printed the pair

    A = tournament_from_tiling(n, bits)
    Af = tournament_from_tiling(n, flipped)

    W = compute_W_poly(A, n)
    Wf = compute_W_poly(Af, n)

    t3_before = count_3cycles(A, n)
    t3_after = count_3cycles(Af, n)

    # Compute W + Wf and W - Wf
    W_sum = {}
    W_diff = {}
    all_powers = set(W.keys()) | set(Wf.keys())
    for p in all_powers:
        W_sum[p] = W.get(p, Fraction(0)) + Wf.get(p, Fraction(0))
        W_diff[p] = W.get(p, Fraction(0)) - Wf.get(p, Fraction(0))

    print(f"\n  GS pair: free={idx:04b} <-> free={flip_idx:04b}")
    print(f"    t3: {t3_before} <-> {t3_after}")
    print(f"    W(r):      {poly_str(W)}")
    print(f"    W_flip(r): {poly_str(Wf)}")
    print(f"    W+Wf:      {poly_str(W_sum)}")
    print(f"    W-Wf:      {poly_str(W_diff)}")

    # Check: evaluate at r=1/2
    H = sum(c * Fraction(1,2)**p for p, c in W.items())
    Hf = sum(c * Fraction(1,2)**p for p, c in Wf.items())
    print(f"    H(T)={H}, H(flip)={Hf}")
