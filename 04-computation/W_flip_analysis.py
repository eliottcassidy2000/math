#!/usr/bin/env python3
"""
Deeper analysis of W(r) under GS tiling flip.

From W_poly_gs_flip.py output at n=5:
- W_4 = 120 = 5! always (universal)
- W_2 coefficient tracks t3 via THM-059: w_2 = -1 + 12*t3  (wait, let me check)
  Actually THM-059 says w_{n-3} = (n-2)! * [2*t3 - C(n,3)/2]... let me verify

Key question: Does W(r) - W_flip(r) have a clean form in terms of delta_t3?

Also: what happens to the Eulerian polynomial F_f(r) evaluation at r = -1/2?
F_f(1/2) = 1 for all f. What about F_f(-1/2)?

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict
from math import comb

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
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(n)):
        poly = {0: Fraction(1)}
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            new_poly = {}
            for power, coeff in poly.items():
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly
        for power, coeff in poly.items():
            result[power] += coeff
    return dict(result)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

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

# Compute F_f(r) master polynomials via Eulerian numbers
def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with exactly k descents."""
    if n == 0:
        return 1 if k == 0 else 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def master_poly_F(f):
    """Compute F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k.
    Returns dict: power -> Fraction coefficient."""
    result = defaultdict(lambda: Fraction(0))
    p = Fraction(1, 2)  # r + 1/2 constant part
    q = Fraction(-1, 2)  # r - 1/2 constant part

    for k in range(f+1):
        Ank = eulerian_number(f+1, k)
        # (r+1/2)^{f-k} * (r-1/2)^k
        # Expand each binomially
        poly_term = {0: Fraction(Ank)}
        for _ in range(f-k):
            # Multiply by (r + 1/2)
            new = {}
            for pw, c in poly_term.items():
                new[pw+1] = new.get(pw+1, Fraction(0)) + c
                new[pw] = new.get(pw, Fraction(0)) + c * p
            poly_term = new
        for _ in range(k):
            # Multiply by (r - 1/2)
            new = {}
            for pw, c in poly_term.items():
                new[pw+1] = new.get(pw+1, Fraction(0)) + c
                new[pw] = new.get(pw, Fraction(0)) + c * q
            poly_term = new
        for pw, c in poly_term.items():
            result[pw] += c
    return dict(result)

def poly_str(poly, var='r'):
    terms = []
    for power in sorted(poly.keys(), reverse=True):
        c = poly[power]
        if c == 0: continue
        if power == 0:
            terms.append(str(c))
        elif power == 1:
            terms.append(f"{c}*{var}")
        else:
            terms.append(f"{c}*{var}^{power}")
    return " + ".join(terms) if terms else "0"

def poly_eval(poly, r):
    return sum(c * r**p for p, c in poly.items())

print("="*70)
print("PART 1: F_f(-1/2) values — symmetry under flip?")
print("="*70)
for f in range(9):
    Ff = master_poly_F(f)
    val_half = poly_eval(Ff, Fraction(1,2))
    val_neg_half = poly_eval(Ff, Fraction(-1,2))
    val_zero = poly_eval(Ff, Fraction(0))
    print(f"  F_{f}(r) = {poly_str(Ff)}")
    print(f"    F_{f}(1/2) = {val_half},  F_{f}(-1/2) = {val_neg_half},  F_{f}(0) = {val_zero}")

print()
print("="*70)
print("PART 2: W(r) decomposition for GS flip pairs at n=5")
print("="*70)

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# THM-059 says: W(r) = C_0(r) + C_t3(r)*t3 + C_t5(r)*t5
# C_0(r) = F_4(r) (since n-1=4, free=4)
# C_t3(r) = 2*F_2(r) (parts=1, |pi|=1, free=4-2=2)
# C_t5(r) = 2*F_0(r) (parts=1, |pi|=2, free=4-4=0)

F4 = master_poly_F(4)
F2 = master_poly_F(2)
F0 = master_poly_F(0)
print(f"\n  F_4(r) = {poly_str(F4)}")
print(f"  F_2(r) = {poly_str(F2)}")
print(f"  F_0(r) = {poly_str(F0)}")
print(f"  C_0(r) = F_4(r) = {poly_str(F4)}")
C_t3 = {p: 2*c for p, c in F2.items()}
C_t5 = {p: 2*c for p, c in F0.items()}
print(f"  C_t3(r) = 2*F_2(r) = {poly_str(C_t3)}")
print(f"  C_t5(r) = 2*F_0(r) = {poly_str(C_t5)}")

print()
print("VERIFICATION: W(r) = C_0 + C_t3*t3 + C_t5*t5")
print("-"*60)

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)
    if idx > flip_idx:
        continue

    A = tournament_from_tiling(n, bits)
    Af = tournament_from_tiling(n, flipped)
    W = compute_W_poly(A, n)
    Wf = compute_W_poly(Af, n)
    t3 = count_3cycles(A, n)
    t3f = count_3cycles(Af, n)

    # Count 5-cycles
    def count_5cycles(A, n):
        count = 0
        for combo in combinations(range(n), 5):
            # n=5: only one 5-vertex set = all vertices
            sub = list(combo)
            for perm in permutations(sub):
                is_cycle = all(A[perm[i]][perm[(i+1)%5]] for i in range(5))
                if is_cycle:
                    count += 1
            # Divide by 5 for rotational equivalence
        return count // 5

    if n == 5:
        t5 = count_5cycles(A, n)
        t5f = count_5cycles(Af, n)
    else:
        t5 = t5f = 0

    # Predicted W using THM-059
    W_pred = {}
    for p in set(list(F4.keys()) + list(F2.keys()) + list(F0.keys())):
        W_pred[p] = F4.get(p, Fraction(0)) + Fraction(2)*F2.get(p, Fraction(0))*t3 + Fraction(2)*F0.get(p, Fraction(0))*t5

    # W - Wf using THM-059:
    # W - Wf = C_t3*(t3-t3f) + C_t5*(t5-t5f)
    delta_t3 = t3 - t3f
    delta_t5 = t5 - t5f if n == 5 else 0
    W_diff_pred = {}
    for p in set(list(F2.keys()) + list(F0.keys())):
        W_diff_pred[p] = Fraction(2)*F2.get(p, Fraction(0))*delta_t3 + Fraction(2)*F0.get(p, Fraction(0))*delta_t5

    # Actual W - Wf
    all_powers = set(W.keys()) | set(Wf.keys())
    W_diff = {p: W.get(p, Fraction(0)) - Wf.get(p, Fraction(0)) for p in all_powers}

    match = all(W_diff.get(p, 0) == W_diff_pred.get(p, 0) for p in all_powers | set(W_diff_pred.keys()))

    print(f"  pair {idx:04b}<->{flip_idx:04b}: t3={t3}->{t3f} (dt3={delta_t3}), t5={t5}->{t5f} (dt5={delta_t5})")
    print(f"    W-Wf actual:    {poly_str(W_diff)}")
    print(f"    W-Wf predicted: {poly_str(W_diff_pred)}")
    print(f"    MATCH: {match}")

print()
print("="*70)
print("PART 3: What is W(-r) vs W_flip(r)?")
print("="*70)
print("  If flip maps s_e -> -s_e for non-backbone, then W_flip(r) != W(-r)")
print("  But via THM-059: W(r) = F_{n-1}(r) + sum 2*F_f(r)*I(T)")
print("  F_f has parity f, so F_f(-r) = (-1)^f * F_f(r)")
print()

for idx in range(min(4, len(gs_tilings))):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    A = tournament_from_tiling(n, bits)
    W = compute_W_poly(A, n)

    # Compute W(-r): negate odd powers
    W_neg = {}
    for p, c in W.items():
        W_neg[p] = c * ((-1)**p)

    Af = tournament_from_tiling(n, flipped)
    Wf = compute_W_poly(Af, n)

    print(f"  tiling {idx:04b}:")
    print(f"    W(r)  = {poly_str(W)}")
    print(f"    W(-r) = {poly_str(W_neg)}")
    print(f"    Wflip = {poly_str(Wf)}")
    print(f"    W(-r) == Wflip? {W_neg == Wf}")

print()
print("="*70)
print("PART 4: The 'anti' evaluation W(-1/2)")
print("="*70)
print("  W(1/2) = H(T) counts Ham paths. What is W(-1/2)?")

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    A = tournament_from_tiling(n, bits)
    W = compute_W_poly(A, n)
    t3 = count_3cycles(A, n)

    H = poly_eval(W, Fraction(1,2))
    anti = poly_eval(W, Fraction(-1,2))

    print(f"  tiling {idx:04b}: t3={t3}, H=W(1/2)={H}, W(-1/2)={anti}, H+anti={H+anti}, H-anti={H-anti}")

print()
print("="*70)
print("PART 5: t3 parity and W(0)")
print("="*70)
print("  W(0) = constant term = F_{n-1}(0) + 2*F_{n-3}(0)*t3 + 2*F_{n-5}(0)*t5")
print(f"  F_4(0) = {poly_eval(F4, Fraction(0))}")
print(f"  F_2(0) = {poly_eval(F2, Fraction(0))}")
print(f"  F_0(0) = {poly_eval(F0, Fraction(0))}")
print(f"  So W(0) = {poly_eval(F4, Fraction(0))} + {2*poly_eval(F2, Fraction(0))}*t3 + {2*poly_eval(F0, Fraction(0))}*t5")
print(f"         = 1 - t3 + 2*t5")
print()

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    A = tournament_from_tiling(n, bits)
    W = compute_W_poly(A, n)
    t3 = count_3cycles(A, n)

    def count_5cycles_simple(A, n):
        count = 0
        for perm in permutations(range(n)):
            is_cycle = all(A[perm[i]][perm[(i+1)%n]] for i in range(n))
            if is_cycle:
                count += 1
        return count // n

    t5 = count_5cycles_simple(A, n) if n == 5 else 0
    W0 = poly_eval(W, Fraction(0))
    pred = 1 - t3 + 2*t5
    print(f"  tiling {idx:04b}: t3={t3}, t5={t5}, W(0)={W0}, predicted={pred}, match={W0==pred}")
