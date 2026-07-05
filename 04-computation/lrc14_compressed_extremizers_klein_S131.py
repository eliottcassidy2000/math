#!/usr/bin/env python3
"""
klein-2026-07-04-S131 (HYP-4093) - CHARACTERIZING THE COMPRESSED 1/13-EXTREMIZERS (dilated deep-wells).

CORRECTION FIRST: my initial compressed-floor pass (minimal-tightener + blocks) reported min 7/89 --
WRONG, the same minimal-tightener-scope error mac-mini-S46 corrected. The TRUE compressed floor is
**1/13** (mac-mini-S46), attained by DILATED DEEP-WELLS like 3*{1..12} u {182} that a minimal-tightener
enumeration misses. This script confirms the 1/13 floor INCLUDING dilated families and characterizes
the extremizers (which the compressed-peel proof must handle).

Offset-forcer/free-rider (mac-mini-S46): at the base optimum t*=1/(13c), the killer's phase is
integer when DOMINANT (deep well 182/13=14 -> kills optimum -> 14/183) but NON-integer when
COMPRESSED (dilated 182/39=14/3, ||.||=1/3 -> free rider -> M=M(base)=1/13). The 13x line = the
integer/non-integer line. So dilated deep-wells c*{1..12} u {killer} (c>=2, compressed) attain 1/13.

Test: (a) confirm min over compressed covering families INCLUDING dilated deep-wells is 1/13 (none
below); (b) characterize ALL compressed families attaining exactly 1/13 -- are they exactly the
dilated deep-wells c*{1..12} u {182d} (compressed)?  Exact.
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations

N = 14
ONE13 = F(1, 13); ONE14 = F(1, 14); DW = F(14, 183)

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

def cov(S): return all(any(v % q == 0 for v in S) for q in range(2, N + 1))
def comp(S): s = sorted(S); return s[-1] <= 13 * s[-2]
def prim(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1
def qc(S): return min(2 * max(S) + 2, 500)

print(f"1/14={float(ONE14):.6f}  14/183={float(DW):.6f}  1/13={float(ONE13):.6f}")
print("=" * 86)
print("(a) DILATED DEEP-WELLS c*{1..12} u {182d}: compressed? covering? M?")
print(f"{'c':>3} {'killer':>7} {'compressed':>11} {'covering':>9} {'primitive':>10} {'M':>10}")
attain = []
for c in range(1, 8):
    base = [c * i for i in range(1, 13)]
    for d in range(1, 6):
        killer = 182 * d
        S = sorted(base + [killer])
        if len(set(S)) != 13: continue
        cp = comp(S); cv = cov(S); pr = prim(S)
        if cv:
            M = Mval(S, qc(S))
            tag = ""
            if cp and M == ONE13: tag = " <-- 1/13 extremizer"; attain.append(S)
            if c <= 5 and d <= 3:
                print(f"{c:>3} {killer:>7} {str(cp):>11} {str(cv):>9} {str(pr):>10} {str(M):>10}{tag}")

print("=" * 86)
print("(b) broad compressed covering search INCLUDING dilations c*U for small U -- any M < 1/13?")
# dilations of the S128 minimal-tightener families + dilated deep wells
from itertools import combinations as comb
def miss(A): return [q for q in range(2, N + 1) if not any(a % q == 0 for a in A)]
def parts(items, k):
    items = list(items)
    if k <= 0: return
    if k == 1: yield [items]; return
    if len(items) < k: return
    f, rest = items[0], items[1:]
    for p in parts(rest, k - 1): yield [[f]] + p
    for p in parts(rest, k):
        for i in range(len(p)): yield p[:i] + [[f] + p[i]] + p[i + 1:]
def st(block, used):
    L = 1
    for q in block: L = lcm(L, q)
    x = ((14 + L - 1) // L) * L
    while x in used: x += L
    return x
AP = list(range(1, 14)); base_fams = set()
for d in range(0, 4):
    for drop in comb(AP, d):
        A = [v for v in AP if v not in drop]; Qm = miss(A)
        if len(Qm) < d: continue
        for part in parts(Qm, d):
            used = set(A); T = []
            for bl in part: t = st(bl, used); T.append(t); used.add(t)
            S = tuple(sorted(A + T))
            if len(S) == 13 and not miss(list(S)): base_fams.add(S)

below13 = []; at13 = []; nchk = 0
for S0 in base_fams:
    for c in range(1, 5):
        S = [c * v for v in S0]
        if len(set(S)) != 13 or not cov(S) or not comp(S): continue
        nchk += 1
        M = Mval(S, qc(S))
        if M < ONE13: below13.append((M, S))
        elif M == ONE13: at13.append(S)
print(f"  checked {nchk} (dilated) compressed covering families; below 1/13: {len(below13)}; at exactly 1/13: {len(at13)}")
for M, S in sorted(below13)[:8]: print(f"    BELOW 1/13: M={M} (~{float(M):.6f}) {S} <-- !!!")

print("=" * 86)
print(f"1/13-extremizers found (dilated deep-wells c*{{1..12}} u {{182d}}): {len(attain)}")
for S in attain[:12]:
    c = S[0]  # smallest element = c*1 = c
    print(f"    c={c}: {S}")
print()
print("READING: if NO compressed covering family is below 1/13 and the 1/13-extremizers are exactly")
print("the dilated deep-wells c*{1..12} u {182d} (compressed, c>=2), then mac-mini-S46's tight target")
print("compressed => M >= 1/13 is confirmed AND the extremizers are characterized -- the compressed-")
print("peel proof must (only) handle the dilated deep-wells (killer = free rider at t*=1/(13c)).")
print("DONE")
