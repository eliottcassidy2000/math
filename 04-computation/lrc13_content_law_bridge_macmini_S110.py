#!/usr/bin/env python3
"""
The content law — the n=12 completeness bridge as one invariant (THM-1006)
=========================================================================
mac-mini-2026-07-18-S110.  Owner directive: attempt the completeness bridge
invariant (codex-S64 sec.6).

Engine: klein THM-1002 writes M(A) = val/q on the pair-sum ruler, with q | v_i+v_j
for an active pair, hence q <= 2*max(A).

Facts verified here:
  (A) DILATION LAW:  M(cA)=M(A), q(cA)=c q(A), val(cA)=c val(A).
      => val(A) = gcd(A) * val(A/gcd(A))  => val >= gcd ALWAYS (free direction).
  (B) TIGHT <=> q = 13*val  (debt d := 13 val - q; tight <=> d=0).
      With q <= 2 max A:  tight => max(A) >= 13 val / 2.
  (C) val IS THM-769's sheet number s (Q = 13 s and q = 13 val).
      val=1 <=> shallow (q=13); val>=2 <=> deep.  Dilates c*{1..12} realize val=c.
  (D) NEW BOUND: tight => val <= 4/(169 * delta(A\\{max})), delta = widest arc with
      phi > 1/13 (THM-1001).  Sharp within ~3.4 on the dilates.
  (E) THE BRIDGE (conjecture, EQUIVALENT to the open problem HYP-6800/6820):
      tight => val = gcd(A), i.e. A = gcd(A)*{1..12}.  Since val>=gcd is free,
      the whole content is  val <= gcd  on the tight locus.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

ONE13 = F(1, 13)

def M_full(S):
    """exact M with a realizing (val, q) on the pair-sum ruler (+ self-cusps, MISTAKE-144)"""
    S = sorted(set(S)); best = F(0); bq = None; bval = None
    dens = {2 * a for a in S}
    for a, b in combinations(S, 2):
        dens.add(a + b); dens.add(b - a)
    dens.discard(0)
    for q in dens:
        for m in range(1, q):
            num = q
            for v in S:
                r = (v * m) % q; d = min(r, q - r)
                if d < num: num = d
                if num == 0: break
            if num > 0:
                c = F(num, q)
                if c > best: best = c; bq = q; bval = num
    return best, bval, bq

def delta_largest(C, L=ONE13):
    C = sorted(set(C)); bnds = set()
    for v in C:
        for a in range(v):
            bnds.add((F(a) + L) / v % 1); bnds.add((F(a) - L) / v % 1)
    B = sorted(bnds); best = F(0)
    for i in range(len(B)):
        lo = B[i]; hi = B[(i + 1) % len(B)]; wid = (hi - lo) % 1
        mid = (lo + wid / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) > L for v in C):
            if wid > best: best = wid
    return best

print("=" * 74)
print("(A) DILATION LAW and val >= gcd (free direction)")
print("=" * 74)
for base in ([1,2,3,4,5,6,7,8,9,10,11,13], list(range(1,13))):
    M0, v0, q0 = M_full(base)
    for c in (1, 2, 3):
        A = [c * x for x in base]; M, val, q = M_full(A)
        ok = (M == M0 and val == c * v0 and q == c * q0)
        print(f"  base max={max(base):2d} c={c}: M={str(M):8s} val={val:3d}(=c*{v0}) q={q:4d}(=c*{q0})  law: {ok}")
print("  => val(A) = gcd(A)*val(A/gcd) hence val >= gcd for every set.")

print()
print("=" * 74)
print("(B)+(C) tight <=> q=13val ; val = sheet number ; dilates realize every val")
print("=" * 74)
for c in range(1, 8):
    A = [c * i for i in range(1, 13)]
    M, val, q = M_full(A); g = reduce(gcd, A); d = 13 * val - q
    print(f"  c={c}: A=c*{{1..12}}  M={M} val={val:2d} q={q:3d} debt={d} gcd={g:2d} "
          f"val==gcd:{val==g} tight:{M==ONE13} max>=13val/2:{max(A) >= 13*val/2}")

print()
print("=" * 74)
print("(D) NEW BOUND  val <= 4/(169*delta(A\\max))  [THM-1001 + q<=2max]")
print("=" * 74)
for c in (1, 2, 3):
    A = [c * i for i in range(1, 13)]
    M, val, q = M_full(A); C = A[:-1]; dl = delta_largest(C)
    bd = F(4, 169) / dl
    print(f"  c={c}: val={val} delta={dl} bound={float(bd):.2f} holds:{val <= bd} "
          f"slack factor={float(bd)/val:.2f}")

print()
print("=" * 74)
print("(E) CONTENT LAW — the bridge (equivalent restatement, NOT proved)")
print("=" * 74)
print("  tight => val = gcd(A)  <=>  A = gcd(A)*{1..12}")
print("  val >= gcd is FREE (A); so the entire content is  val <= gcd  on the tight locus,")
print("  i.e. the primitive quotient has val=1 = the deep branch is empty for primitive sets.")
print("  sporadic emptiness <=> [content law] + [shallow rigidity].")
print()
print("  POSITIVE CONTROL (must fail at n=13): Goddyn-Wong {1..11,13,24}")
GW = [1,2,3,4,5,6,7,8,9,10,11,13,24]
M, val, q = M_full(GW); g = reduce(gcd, GW)
print(f"    GW: M={M} val={val} q={q} gcd={g} val==gcd:{val==g} tight(1/14):{M==F(1,14)}")
print("    -> GW satisfies val=gcd yet is NOT {1..13}: at n=13 it is the SHALLOW RIGIDITY half")
print("       that fails, not val<=gcd. Any proof of (E) must be n=12-specific.")
