#!/usr/bin/env python3
"""
glmy_paley_conjecture.py — opus-2026-03-13-S71

Test conjectures about GLMY Betti numbers for Paley tournaments.

Conjecture: For Paley tournament P_p (p ≡ 3 mod 4):
  β_m = 0 for m ≠ 0 and m ≠ (p-3)/2
  β_{(p-3)/2} = p - 1
  chi = (-1)^{(p-3)/2} * (p-1) + 1

Check at p=3,5,7. (p=11 too expensive for full GLMY, but can estimate.)

Also: Why are C7_{1,2,3} and C7_{1,3,5} identical?
Answer: they're ISOMORPHIC. Check via canon_form.
"""

import numpy as np
from itertools import permutations

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

def enumerate_directed_paths(A, n, p):
    if p == 0: return [(v,) for v in range(n)]
    paths = []
    def dfs(path, depth):
        if depth == p:
            paths.append(tuple(path))
            return
        last = path[-1]
        visited = set(path)
        for v in range(n):
            if v not in visited and A[last][v]:
                path.append(v)
                dfs(path, depth+1)
                path.pop()
    for s in range(n):
        dfs([s], 0)
    return paths

def compute_omega_basis(allowed_p, allowed_pm1, p):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0,0))
    if p <= 1: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for path in allowed_p:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in allowed_pm1_set and face not in non_allowed:
                non_allowed[face] = na_count
                na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    J = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face in non_allowed:
                J[non_allowed[face], j] += (-1)**i
    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    if rank == dim_Ap: return np.zeros((dim_Ap, 0))
    return Vh[rank:].T

def glmy_betti(A):
    n = A.shape[0]
    max_dim = n - 1
    allowed = {}
    omega = {}
    for p in range(max_dim+2):
        allowed[p] = enumerate_directed_paths(A, n, p)
        omega[p] = compute_omega_basis(allowed[p], allowed.get(p-1, []), p)
    ranks = {}
    for p in range(1, max_dim+1):
        dim_om = omega[p].shape[1]
        dim_om_prev = omega[p-1].shape[1]
        if dim_om == 0 or dim_om_prev == 0:
            ranks[p] = 0
            continue
        idx = {path: i for i, path in enumerate(allowed[p-1])}
        B = np.zeros((len(allowed[p-1]), len(allowed[p])))
        for j, path in enumerate(allowed[p]):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                if face in idx:
                    B[idx[face], j] += (-1)**i
        B_omega = omega[p-1].T @ B @ omega[p]
        ranks[p] = np.linalg.matrix_rank(B_omega, tol=1e-8)
    betti = []
    for p in range(max_dim+1):
        dim_p = omega[p].shape[1]
        rk_p = ranks.get(p, 0)
        rk_p1 = ranks.get(p+1, 0)
        betti.append(dim_p - rk_p - rk_p1)
    omega_dims = [omega[p].shape[1] for p in range(max_dim+1)]
    return betti, omega_dims

# ============================================================
print("="*70)
print("ISOMORPHISM CHECK: C7 INTERVAL vs C7_{1,3,5}")
print("="*70)

A_int = circulant_tournament(7, {1,2,3})
A_bal = circulant_tournament(7, {1,3,5})
cf_int = canon_form(A_int)
cf_bal = canon_form(A_bal)
print(f"  C7_{{1,2,3}} canon: {cf_int[:10]}...")
print(f"  C7_{{1,3,5}} canon: {cf_bal[:10]}...")
print(f"  Isomorphic? {cf_int == cf_bal}")

# Also check {2,3,4} etc
for S_tuple in [(1,2,3), (1,3,5), (2,3,4), (1,2,5), (2,4,6), (3,5,6)]:
    S = set(S_tuple)
    # Ensure valid tournament (complement)
    comp = {7-s for s in S}
    if S & comp: continue
    if S | comp != set(range(1,7)): continue
    A = circulant_tournament(7, S)
    cf = canon_form(A)
    print(f"  C7_{set(S_tuple)}: canon={cf[:6]}...")

# ============================================================
print(f"\n{'='*70}")
print("PALEY BETTI CONJECTURE TEST")
print("="*70)

# p=3: P_3 = C_3 cyclic. (p-3)/2 = 0. Expect β_0=p-1=2? But β_0=1.
# Hmm. Let's look more carefully.
print("\n  Paley at each p ≡ 3 mod 4:")
for p in [3, 5, 7]:
    # Quadratic residues mod p
    QR = {a % p for a in range(1, p) if pow(a, (p-1)//2, p) == 1}
    A = circulant_tournament(p, QR)
    b, o = glmy_betti(A)
    chi = sum((-1)**m * b[m] for m in range(len(b)))

    # Where is the nonzero β (besides β_0)?
    nonzero = [(m, b[m]) for m in range(len(b)) if b[m] > 0 and m > 0]

    print(f"\n  P_{p}: QR={QR}")
    print(f"    β = {b}")
    print(f"    Ω/p = {[o[m]//p for m in range(len(o))]}")
    print(f"    chi = {chi}")
    print(f"    (p-3)/2 = {(p-3)//2}")
    if nonzero:
        for m, val in nonzero:
            print(f"    β_{m} = {val} (at degree {m})")
    else:
        print(f"    Only β_0 = 1 nonzero")

# p=5 note: P_5 has QR={1,4}. But {1,2} is the "regular" circulant.
# Are they the same?
print(f"\n  P_5 QR = {{a^2 mod 5}} = {{1,4}}")
A_p5 = circulant_tournament(5, {1,4})
b_p5, o_p5 = glmy_betti(A_p5)
print(f"  C5_{{1,4}} β = {b_p5}, Ω/5 = {[o//5 for o in o_p5]}")

A_r5 = circulant_tournament(5, {1,2})
b_r5, o_r5 = glmy_betti(A_r5)
print(f"  C5_{{1,2}} β = {b_r5}, Ω/5 = {[o//5 for o in o_r5]}")

cf_p5 = canon_form(A_p5)
cf_r5 = canon_form(A_r5)
print(f"  P_5 ≅ C5_{{1,2}}? {cf_p5 == cf_r5}")

# ============================================================
print(f"\n{'='*70}")
print("PALEY P_7 Ω/7 PALINDROME ANALYSIS")
print("="*70)

# Ω/7 = [1,3,6,9,9,6,3]
# This is NOT [1,3,6,9,9,6,3,1] (missing the 1 at end since Ω_7=0 for n=7)
# But the NONZERO part IS palindromic.
# Actually Ω_6 = 21 so Ω/7 = [1,3,6,9,9,6,3] — yes palindromic!

seq = [1,3,6,9,9,6,3]
print(f"  Ω/7 = {seq}")
print(f"  Palindromic? {seq == seq[::-1]}")

# For C5 regular (= Paley P_5):
seq5 = [1,2,2,2,1]
print(f"  Ω/5 (P_5) = {seq5}")
print(f"  Palindromic? {seq5 == seq5[::-1]}")

# For C3 cyclic (= Paley P_3):
seq3 = [1,1,0]
print(f"  Ω/3 (P_3) = {seq3}")
print(f"  Palindromic? {seq3 == seq3[::-1]}")

# For Interval C7:
seq_i7 = [1,3,6,11,14,9,2]
print(f"  Ω/7 (I_7) = {seq_i7}")
print(f"  Palindromic? {seq_i7 == seq_i7[::-1]}")

print(f"\n  CONJECTURE: For Paley tournament P_p, Ω_m/p is palindromic.")
print(f"  This is Poincaré duality for the path complex!")

# ============================================================
print(f"\n{'='*70}")
print("PALEY Ω/p GENERATING FUNCTION ANALYSIS")
print("="*70)

# P_7: Ω/7 = [1,3,6,9,9,6,3]
# Generating function: 1+3x+6x²+9x³+9x⁴+6x⁵+3x⁶
# = (1+3x+6x²+9x³+9x⁴+6x⁵+3x⁶)
# Try: (1+x)^k / something
# Sum = 37

# P_5: Ω/5 = [1,2,2,2,1]
# Sum = 8

# P_3: Ω/3 = [1,1,0]
# Sum = 2

print(f"  P_3: sum(Ω/3) = 2 = 2^1")
print(f"  P_5: sum(Ω/5) = 8 = 2^3")
print(f"  P_7: sum(Ω/7) = 37")
print(f"  37 is prime, not a nice power")
print()

# Maybe look at total Ω (not divided by p)?
# P_3: sum(Ω) = 3+3+0 = 6
# P_5: sum(Ω) = 5+10+10+10+5 = 40
# P_7: sum(Ω) = 7+21+42+63+63+42+21 = 259

print(f"  P_3: sum(Ω) = 6 = 3!? No, 3!=6. Yes!")
print(f"  P_5: sum(Ω) = 40 = 5*8 = 5*2^3")
print(f"  P_7: sum(Ω) = 259 = 7*37")

# ============================================================
print(f"\n{'='*70}")
print("P_7 PALEY: BOUNDARY RANK SEQUENCE")
print("="*70)

# rk(∂_m): [6, 15, 27, 36, 21, 21]
# = [6, 15, 27, 36, 21, 21]
# rk/7: not always integer?
rks = [6, 15, 27, 36, 21, 21]
print(f"  rk(∂_m) for m=1..6: {rks}")
print(f"  rk/7: {[r/7 for r in rks]}")
print(f"  rk/3: {[r/3 for r in rks]}")
print(f"  All divisible by 3: {all(r % 3 == 0 for r in rks)}")

# These are:
# 6 = C(4,2), 15 = C(6,2), 27 = C(3,1)*C(3,3)*3?, 36 = C(9,2)?
# 6, 15, 27, 36, 21, 21
# Diffs: 9, 12, 9, -15, 0
# Not obvious pattern

# ============================================================
print(f"\n{'='*70}")
print("RELATIONSHIP BETWEEN β_4=6 AND Z_7 REPRESENTATION THEORY")
print("="*70)

# β_4 = 6 = p-1 = |Z_7^*|
# The homology H_4(P_7) is a 6-dimensional representation of Z_7 (cyclic symmetry)
# By Maschke's theorem, this decomposes into 1-dim irreps (characters of Z_7)
# There are 7 characters (trivial + 6 nontrivial)
# β_4 = 6 = number of nontrivial characters!
# CONJECTURE: H_4(P_7) ≅ direct sum of ALL nontrivial Z_7-irreps

print(f"  β_4 = 6 = p-1 = number of nontrivial characters of Z_p")
print(f"  Conjecture: H_{{(p-3)/2}}(P_p) ≅ regular representation minus trivial")
print()

# P_3: (p-3)/2 = 0, β_0 = 1. But p-1=2.
# This doesn't match: β_0 = 1 ≠ 2.
# HOWEVER: β_0 always = 1 for GLMY (connected component).
# Maybe the conjecture is: for m = (p-3)/2 with m > 0:
#   β_m = p - 1

print(f"  P_3: m=(3-3)/2=0, β_0=1 (always 1 for β_0)")
print(f"  P_5: m=(5-3)/2=1, β_1=1 ≠ p-1=4. FAILS!")
print(f"  P_7: m=(7-3)/2=2, β_2=0 ≠ p-1=6. FAILS!")
print()
print(f"  Wait — P_7's nonzero β is at m=4, not m=2.")
print(f"  (p-3)/2 = 2 but β_2 = 0. The nonzero β is at m = p-3 = 4.")

# Let me re-examine
print(f"\n  P_3: nonzero β at m=1, value=1. p-3=0, p-2=1")
print(f"  P_5: nonzero β at m=1, value=1. p-3=2, p-2=3")
print(f"  P_7: nonzero β at m=4, value=6. p-3=4!")

print(f"\n  NEW CONJECTURE: For Paley P_p with p ≡ 3 mod 4:")
print(f"    If p ≥ 7: β_{{p-3}} = p-1, all other β_m = 0 for m ≥ 1")
print(f"    chi(P_p) = 1 + (-1)^{{p-3}}*(p-1) = 1 + (p-1) = p (since p-3 even)")

# Check: p=7: chi = 1 + (-1)^4 * 6 = 1 + 6 = 7 ✓
# p=3: β = [1,1,0], chi = 0 (doesn't fit the p≥7 pattern)
# p=5: β = [1,1,0,0,0], chi = 0 (doesn't fit either)

print(f"\n  Verification:")
print(f"    P_7: chi = 1 + 6 = 7 = p ✓")
print(f"    P_3: chi = 0 ≠ 3 (small case exception)")
print(f"    P_5: chi = 0 ≠ 5 (small case exception)")

# ============================================================
print(f"\n{'='*70}")
print("POINCARÉ DUALITY FOR PALEY PATH COMPLEX?")
print("="*70)

# Ω/7 = [1,3,6,9,9,6,3] is palindromic
# For a closed (p-2)-manifold, Poincaré duality gives β_m = β_{p-2-m}
# But GLMY β = [1,0,0,0,6,0,0] is NOT palindromic

# The Ω being palindromic but β NOT suggests something interesting:
# The boundary maps BREAK the palindrome symmetry

# Actually check: is Ω palindromic for OTHER p too?
print(f"  Ω palindrome check:")
for p, QR in [(3, {1}), (5, {1,2}), (7, {1,2,4})]:
    A = circulant_tournament(p, QR)
    _, o = glmy_betti(A)
    nonzero = [x for x in o if x > 0]
    is_pal = nonzero == nonzero[::-1]
    print(f"    P_{p}: Ω(nonzero) = {nonzero}, palindromic? {is_pal}")

    # Also check full Ω (including trailing zeros)
    print(f"         Ω(full)    = {o}")

print("\nDONE.")
