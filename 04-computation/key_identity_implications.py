#!/usr/bin/env python3
"""
Explore implications of the Key Identity (THM-030) for the broader project.

KEY IDENTITY: B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)

Questions to investigate:
1. What does T = r*Sigma mean for T evaluated at r=1 (the actual tournament case)?
   At r = 1/2 (the 0/1 case)? At c = 2r?
2. Does the proof give new information about Hamiltonian path counts H(T)?
3. Can we extract the independence polynomial from the transfer matrix?
4. What is the relationship between Sigma and H(T)?
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Rational, factor, collect
import sys

MAX_N = 6
r = symbols('r')
s = {}
for i in range(MAX_N):
    for j in range(i+1, MAX_N):
        s[(i,j)] = symbols(f's{i}{j}')
        s[(j,i)] = -s[(i,j)]

def t(i, j):
    return r + s[(i,j)]

def path_weight(path):
    w = 1
    for i in range(len(path)-1):
        w *= t(path[i], path[i+1])
    return expand(w)

def B_v(v, W):
    total = 0
    for p in permutations(W):
        if p[0] == v:
            total += path_weight(p)
    return expand(total)

def E_v(v, W):
    total = 0
    for p in permutations(W):
        if p[-1] == v:
            total += path_weight(p)
    return expand(total)

def T_total(W):
    total = 0
    for p in permutations(W):
        total += path_weight(p)
    return expand(total)

def M_entry(a, b, W):
    W = tuple(sorted(W))
    U = [v for v in W if v != a and v != b]
    total = 0
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            ea = E_v(a, tuple(sorted(list(S) + [a])))
            bb = B_v(b, tuple(sorted(R + [b])))
            total += ((-1)**k) * ea * bb
    return expand(total)

def Sigma(W):
    total = 0
    for a in W:
        for b in W:
            if a != b:
                total += M_entry(a, b, W)
    return expand(total)

def diag_sum(W):
    """Sum of diagonal entries M[a,a]."""
    total = 0
    for a in W:
        total += M_entry(a, a, W)
    return expand(total)

print("=" * 70)
print("IMPLICATION 1: T = r*Sigma at r = c/2 for c-tournaments")
print("=" * 70)
print()
print("For a {0,1}-tournament, c=1 so r=1/2.")
print("T(W)|_{r=1/2} = H(T_W) = number of Hamiltonian paths through W.")
print("Sigma(W)|_{r=1/2} = sum_{a!=b} M[a,b]|_{r=1/2}.")
print()
print("For even m: H = (1/2) * Sigma|_{r=1/2}")
print("For odd m: Sigma|_{r=1/2} = 0")
print()

# Let's check at specific tournaments
# Transitive tournament on {0,1,...,m-1}: t(i,j) = 1 if i<j, 0 if i>j
# So s_{ij} = 1/2 if i<j, -1/2 if i>j

print("--- Transitive tournament (t_{ij} = 1 if i<j, 0 else) ---")
for m in [3, 4, 5]:
    W = tuple(range(m))
    # Substitute r=1/2, s_{ij} = 1/2 for i<j
    subs = {r: Rational(1, 2)}
    for i in range(m):
        for j in range(i+1, m):
            subs[s[(i,j)]] = Rational(1, 2)

    tt = T_total(W).subs(subs)
    sig = Sigma(W).subs(subs)
    print(f"  m={m}: H = T|_{{r=1/2}} = {tt}")
    print(f"        Sigma|_{{r=1/2}} = {sig}")
    if m % 2 == 0:
        print(f"        H = r*Sigma? {tt} = {Rational(1,2)*sig}? {tt == Rational(1,2)*sig}")
    else:
        print(f"        Sigma = 0? {sig == 0}")

print()
print("=" * 70)
print("IMPLICATION 2: Diagonal entries of M")
print("=" * 70)
print()
print("M[a,a] = sum_S (-1)^|S| E_a(S+a) B_a(R+a)")
print("This is the 'self-transfer' entry. What is its relationship to H?")
print()

for m in [3, 4, 5]:
    W = tuple(range(m))
    print(f"  m={m}:")
    ds = diag_sum(W)
    sig = Sigma(W)
    tt = T_total(W)

    # Total sum (including diagonal) = sum of ALL M[a,b]
    total = expand(sig + ds)
    print(f"    diag_sum = {ds}")
    print(f"    off_diag_sum (Sigma) = {sig}")
    print(f"    total_sum = {total}")
    print(f"    T(W) = {tt}")

    # Check: is total_sum = T?
    print(f"    total_sum = T? {expand(total - tt) == 0}")
    # Check: is diag_sum = trace relation from THM-027?
    # THM-027 says trace = H for odd n, 0 for even n
    # But here n = m+? No, W IS the full set.
    # For odd m: trace should be... let me check
    subs = {r: Rational(1, 2)}
    for i in range(m):
        for j in range(i+1, m):
            subs[s[(i,j)]] = Rational(1, 2)  # transitive
    ds_num = ds.subs(subs)
    print(f"    diag_sum (transitive) = {ds_num}")

print()
print("=" * 70)
print("IMPLICATION 3: Full matrix sum = T (always?)")
print("=" * 70)
print()
print("Checking: sum_{a,b} M[a,b] (including diagonal) = T(W)?")
for m in [2, 3, 4, 5]:
    W = tuple(range(m))
    total = 0
    for a in W:
        for b in W:
            total += M_entry(a, b, W)
    total = expand(total)
    tt = T_total(W)
    print(f"  m={m}: sum_{{a,b}} M[a,b] = T? {expand(total - tt) == 0}")
    if expand(total - tt) != 0:
        print(f"    diff = {expand(total - tt)}")

print()
print("=" * 70)
print("IMPLICATION 4: Relationship between trace(M) and H(T)")
print("=" * 70)
print()
print("THM-027: tr(M) = H for odd |W|, = 0 for even |W|.")
print("Combined with Sigma identity: T = tr(M) + Sigma.")
print()
print("For odd m: T = tr(M) + 0 = tr(M). So T = tr(M) = H (at r=1/2).")
print("Wait, T = H always (T is total path weight = H at r=1/2).")
print("So for odd m: tr(M) = T. This IS THM-027.")
print()
print("For even m: T = tr(M) + Sigma. And T = r*Sigma.")
print("So: r*Sigma = tr(M) + Sigma => tr(M) = (r-1)*Sigma.")
print()

for m in [4, 5, 6]:
    W = tuple(range(m))
    tr_M = diag_sum(W)
    sig = Sigma(W)
    tt = T_total(W)

    # Check T = tr(M) + Sigma
    check1 = expand(tt - tr_M - sig)
    print(f"  m={m}: T = tr(M) + Sigma? {check1 == 0}")

    if m % 2 == 0:
        # For even m: tr(M) = (r-1)*Sigma
        check2 = expand(tr_M - (r-1)*sig)
        print(f"    tr(M) = (r-1)*Sigma? {check2 == 0}")

print()
print("=" * 70)
print("IMPLICATION 5: New identity — tr(M) = (r-1)*Sigma for even m")
print("=" * 70)
print()
print("This means at r = 1/2 (standard tournament):")
print("  tr(M) = -1/2 * Sigma = -H (since H = r*Sigma = Sigma/2)")
print("  So tr(M)|_{r=1/2} = -H for even m.")
print()
print("And at r = 1 (c=2 tournament):")
print("  tr(M) = 0 for even m.")
print("  This makes sense — THM-027 says tr(M) = H for odd, 0 for even.")
print("  But wait, THM-027 was for {0,1} tournaments. Need to check general c.")
print()

# At r=0 (pure skew):
print("At r = 0 (pure skew tournament):")
for m in [3, 4, 5]:
    W = tuple(range(m))
    tr_M = diag_sum(W)
    sig = Sigma(W)
    tt = T_total(W)

    tr_0 = tr_M.subs(r, 0)
    sig_0 = sig.subs(r, 0)
    tt_0 = tt.subs(r, 0)

    print(f"  m={m}: T(0) = {tt_0}, tr(M)(0) = {tr_0}, Sigma(0) = {sig_0}")
    if m % 2 == 0:
        check = expand(tr_0 + sig_0)
        print(f"    tr(M) = -Sigma? {check == 0}")
    else:
        check = expand(tr_0 - tt_0)
        print(f"    tr(M) = T? {check == 0}")

print()
print("=" * 70)
print("IMPLICATION 6: Full structure of the total matrix sum")
print("=" * 70)
print()
print("Define Tau = sum_{a,b in W} M[a,b] (full sum including diagonal).")
print("Then Tau = tr(M) + Sigma.")
print()
print("From T = Tau (if true), we get T = tr(M) + Sigma.")
print("For even m: T = r*Sigma, so tr(M) = (r-1)*Sigma, and Tau = T.")
print("For odd m: Sigma = 0, so Tau = tr(M) = T.")
print("So in BOTH cases: Tau = T.")
print()
print("This means: sum_{a,b} M[a,b] = T(W) for ALL m.")
print("This is a CLEAN IDENTITY connecting the transfer matrix to total path weight.")
print()

# Verify once more
for m in [2, 3, 4, 5]:
    W = tuple(range(m))
    tau = 0
    for a in W:
        for b in W:
            tau += M_entry(a, b, W)
    tau = expand(tau)
    tt = T_total(W)
    print(f"  m={m}: Tau = T? {expand(tau - tt) == 0}")

print()
print("=" * 70)
print("SUMMARY OF NEW IDENTITIES FROM THM-030")
print("=" * 70)
print("""
1. KEY IDENTITY: B_b + (-1)^m E_b = 2r * col_sum(b)    [THM-030]
2. SIGMA IDENTITY: T[1+(-1)^m] = 2r * Sigma            [Corollary]
3. FULL SUM: sum_{a,b} M[a,b] = T(W)                    [Verified]
4. TRACE IDENTITY:
   - Odd m: tr(M) = T                                    [THM-027]
   - Even m: tr(M) = (r-1) * Sigma = (1 - 1/r) * T      [NEW]
5. AT r = 1/2 (standard tournament, c=1):
   - Even m: tr(M) = -H, Sigma = 2H
   - Odd m: tr(M) = H, Sigma = 0
6. M[a,b] = M[b,a] (symmetry)                            [Corollary]
7. M[a,b] has only even powers of r                      [Corollary]
""")
