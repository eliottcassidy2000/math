# opus-2026-07-17-S343 -- HYP-7270: THE TWO-VARIABLE FOLDED IDENTITY.
# Claim (g = gcd(a,b); fold_M(r) = r*(M - r), r taken mod M):
#   mu(D_a cap D_b) = [4ab + fold_{14g}((a+b) mod 14g)
#                          - fold_{14g}((b-a) mod 14g)] / (196 ab)   ... (I)
# derived by min(2a, X) = X - (X - 2a)+ telescoping the sawtooth sum into
# F(S) - F(b-a), 14g*F(S) = S^2 + fold_{14g}(S mod 14g).
# Consequences: boxeph's consecutive form is the b = a+1 slice; the ENTIRE
# floor table becomes analytic: mu >= 1/49 - fold-corrections/(196ab).
from fractions import Fraction
from math import gcd
import random

F = Fraction
LAM = F(1, 14)

def mu_direct(a, b):
    def teeth(x):
        w = LAM / x
        out = []
        for j in range(x):
            lo, hi = F(j, x) - w, F(j, x) + w
            out.append((max(lo, F(0)), min(hi, F(1))))
            if lo < 0: out.append((lo + 1, F(1)))
            if hi > 1: out.append((F(0), hi - 1))
        return sorted(out)
    u, v = teeth(a), teeth(b)
    tot, i, j = F(0), 0, 0
    while i < len(u) and j < len(v):
        lo, hi = max(u[i][0], v[j][0]), min(u[i][1], v[j][1])
        if lo < hi: tot += hi - lo
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return tot

def fold(r, M): return r * (M - r)

def mu_formula(a, b):
    g = gcd(a, b)
    M = 14 * g
    S, D = a + b, b - a
    val = 4 * a * b + fold(S % M, M) - fold(D % M, M)
    return F(val, 196 * a * b)

print("(I) THE TWO-VARIABLE FOLDED IDENTITY, exact battery:")
random.seed(343)
ok = bad = 0
worst = None
for _ in range(400):
    a = random.randint(2, 150)
    b = random.randint(a, 13 * a)
    if mu_formula(a, b) == mu_direct(a, b): ok += 1
    else:
        bad += 1
        if bad <= 3: print(f"   MISMATCH at (a,b)=({a},{b}) g={gcd(a,b)}: "
                           f"formula {mu_formula(a,b)} vs {mu_direct(a,b)}")
print(f"   400 random pairs (ratios <= 13, all gcds): {ok} exact, {bad} mismatches")

print()
print("(II) consequences if (I) holds:")
print("   mu = 1/49 + [fold(S%14g) - fold(D%14g)]/(196ab); folds in [0, 49g^2]")
print("   FLOOR: mu >= 1/49 - fold(D%14g)_max-effect/(196ab) >= 1/49 - g^2/(4ab)")
print("        = 1/49 - 1/(4 a' b' * 196/...)  [a'=a/g, b'=b/g: 1/49 - 49g^2/(196ab)]")
# re-derive the S341 table minima from the formula
print("   re-deriving the S341 floor table from (I):")
for rnum in (13, 8, 5, 3, 2):
    best = None
    for a in range(2, 121):
        for b in range(a, rnum * a + 1):
            m = mu_formula(a, b)
            if best is None or m < best[0]: best = (m, a, b)
    print(f"   ratio <= {rnum}: min = {best[0]} at {best[1:]} "
          f"(S341 table check)")
print()
print("(III) the fold is the triangle: fold_14(r) = r(14-r) = the staircase")
print("   leg product at the lam = 1/14 modulus; boxeph's consecutive form is")
print("   the b = a+1 slice (fold_14 with D%14 = 1: fold=13 constant).")
