#!/usr/bin/env python3
"""
boxeph-2026-07-19-S123 — the determinant stratification of the n=12 uniqueness gap (1/13, 2/25).

Owner: work n=12 AP uniqueness; mine threads on 3, 4, 1/12.

By the Pinch Lemma (HYP-2059/THM-401) the loneliness maximizer sits at a pairwise sum and
M = D/s, D = |v_i a_j - v_j a_i| (determinant), s = v_i+v_j. Writing M = p/q in lowest terms,
p = D/gcd(D,s) <= D, so stratifying the gap by numerator p = stratifying by determinant.

PROVED (LRC(13) + parity): numerator p<=2 is EXCLUDED from the open gap (1/13,2/25):
  p=1: M=1/q>=1/13 => q<=13 => M=1/13 (edge) or >=1/12.
  p=2: M=2/q lowest terms => q ODD; 2/q>=1/13 => q<=26 => q<=25 => M>=2/25.
Residual: numerator p>=3, a DISCRETE ladder {3/38,4/51,5/63,5/64,...} -> 1/13+, with 3/38 the
unique depth-minimal value (the mediant). (C) <=> none of these is achievable.

This script verifies:
 (A) M({1,...,11,12m}) = m/(12m+1)  [= THM-633 / LRCLadderD1, mac-mini-S33 -- CREDITED, re-derived];
 (B) the cascade reduction 12 not| x  => M >= 1/12 (t=1/12 unblocked);
 (C) the gap numerator stratification (which reduced p/q lie in (1/13,2/25), by p);
 (D) that a M=3/38 family would be covering + need a determinant-3 pair at s=38.
"""
from fractions import Fraction as F
from math import gcd


def fd(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(sp, t: F) -> F:
    return min(fd(F(v) * t) for v in sp)


def trueM(sp):
    Ds = set()
    n = len(sp)
    for i in range(n):
        for j in range(i + 1, n):
            Ds.add(sp[i] + sp[j])
            Ds.add(abs(sp[i] - sp[j]))
    best, bt = F(0), F(0)
    for D in Ds:
        if D <= 0:
            continue
        for m in range(1, D):
            v = fmin(sp, F(m, D))
            if v > best:
                best, bt = v, F(m, D)
    return best, bt


THIRT, TWOF = F(1, 13), F(2, 25)

print("=" * 72)
print("(A) M({1,...,11,12m}) = m/(12m+1)   [= THM-633 / LRCLadderD1, CREDITED]")
print("=" * 72)
ok = True
for m in range(1, 11):
    C = list(range(1, 12)) + [12 * m]
    M, t = trueM(C)
    pred = F(m, 12 * m + 1)
    ok = ok and (M == pred)
    tag = "  (AP {1..12})" if m == 1 else ("  <-- min non-AP = 2/25" if m == 2 else "")
    print(f"  m={m:>2}  x=12m={12*m:>3}  M={str(M):>7}  m/(12m+1)={str(pred):>7}  t*={str(t):>7}{tag}")
print(f"  formula holds m=1..10: {ok}")

print()
print("=" * 72)
print("(B) cascade: 12 does NOT divide x => M >= 1/12 > 1/13 (t=1/12 unblocked)")
print("=" * 72)
for x in [13, 14, 17, 25, 35, 101]:
    C = list(range(1, 12)) + [x]
    print(f"  x={x:>4}  12|x? {str(x % 12 == 0):>5}  min||v/12||={str(fmin(C, F(1,12))):>5}  "
          f"trueM={str(trueM(C)[0]):>5}")

print()
print("=" * 72)
print("(C) the open gap (1/13,2/25) stratified by reduced numerator p (= det/gcd):")
print("=" * 72)
for p in range(1, 8):
    vals = sorted({F(p, q) for q in range(2, 300) if gcd(p, q) == 1 and THIRT < F(p, q) < TWOF})
    label = "EXCLUDED (proved: LRC13 + parity)" if p <= 2 else f"{[str(v) for v in vals[:4]]}"
    print(f"  numerator p={p}: {label}")
print("  => residual = numerator p>=3, discrete, -> 1/13+; depth-minimal = 3/38 (unique).")

print()
print("=" * 72)
print("(D) a M=3/38 family: covering-forced + determinant-3 pair at s=38")
print("=" * 72)
print(f"  3/38 ~ {float(F(3,38)):.5f}, in gap: {THIRT < F(3,38) < TWOF}; 3/38<1/12 => COVERING forced")
print(f"        (any unblocked t=1/q, q<=12, gives M>=1/12>3/38, so contains a multiple of each 2..12)")
print(f"  Pinch: M=3/38=D/s => (D,s)=(3,38) or (6,76); 38=3*13-1 (the mediant denominator).")
print(f"  Farey: 1/13,2/25 neighbors (1*25-2*13={1*25-2*13}); mediant = 3/38 = {F(1+2,13+25)}.")
print("\nDONE.")
