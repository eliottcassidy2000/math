#!/usr/bin/env python3
"""
klein-2026-07-06-S150 - IS THE ONE OPEN NODE OF (C) A FINITE CHECK?

kps-S49 clean (C) skeleton: the sole open node =
   COMPRESSED (max<=13*min), NON-TRANSLATE (not {m..m+11}), NON-AP mod-25 BLOCKER (M_25<2/25)
   => M >= 2/25 (cleared at some q<=37).
Verified (klein-S144: q<=31 to height 650k; mac-mini-S35: q up to 37, Q0 saturates; kps-S49: 615
families clear q<=29). DECISIVE question for formalization strategy: are these families BOUNDED in
height (=> finite check, native_decide) or UNBOUNDED (=> need uniform covering)?

This script: enumerate compressed non-translate non-AP mod-25-blockers by min m, over increasing m,
and see whether (a) the COUNT saturates, (b) min is BOUNDED, (c) they clear at bounded q. Also test
whether they reduce to finitely many RESIDUE CLASSES mod lcm of the covering (the height-uniform
finiteness that would make the node a finite residue check).
Exact.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

TWO25 = F(2, 25)
def cd(a, q): r = a % q; return min(r, q - r)
def Mq(W, q):
    b = 0
    for c in range(1, q // 2 + 1):
        if gcd(c, q) != 1: continue
        m = min(cd(v * c, q) for v in W)
        if m > b: b = m
    return F(b, q)
def is_blocker25(W): return Mq(W, 25) < TWO25          # obstructed at 25
def clearq(W, qmax=60):
    for q in range(2, qmax + 1):
        if Mq(W, q) >= TWO25: return q
    return None
def tg(W):
    g = 0
    for v in W: g = gcd(g, v)
    return g
def is_translate(W): W=sorted(W); return W == list(range(W[0], W[0]+12))
AP = list(range(1, 13))

print("Open node: COMPRESSED (max<=13min) NON-TRANSLATE NON-AP mod-25-BLOCKER (M_25<2/25) => cleared.")
print("Q: bounded height? finite? clearing-q bound?")
print("="*84)

# Enumerate by min m: all 12-subsets with min=m, max<=13m, that are non-translate non-AP blockers.
# For each m, this is subsets of [m, 13m] containing m, size 12 -> too many for large m; SAMPLE for
# large m, EXHAUST for small m. Track: count of blockers, their clearing-q, whether min stays bounded.
from random import seed, sample
seed(150)
print(f"{'min m':>6} {'window':>10} {'#blockers(exh/samp)':>20} {'maxclearq':>10} {'example blocker (clear-q)':>30}")
overall_maxclear = 0
blocker_mins = []
for m in [1,2,3,4,5,6,8,10,15,20,30,50,100,300,1000]:
    hi = 13*m
    pool = list(range(m, hi+1))
    n_block = 0; maxclear = 0; ex = None; sampled = False
    # exhaustive if C(|pool|-1, 11) small, else sample
    from math import comb
    total = comb(len(pool)-1, 11)  # subsets containing m
    if total <= 300000 and m <= 6:
        it = ([m]+list(c) for c in combinations([x for x in pool if x!=m], 11))
    else:
        sampled = True
        it = ([m]+sorted(sample([x for x in pool if x!=m], 11)) for _ in range(200000))
    seen=set()
    for W in it:
        W = tuple(sorted(W))
        if W in seen: continue
        seen.add(W)
        if tg(list(W))!=1: continue
        if list(W)==AP: continue
        if is_translate(list(W)): continue
        if not is_blocker25(list(W)): continue
        n_block += 1
        cq = clearq(list(W), 60)
        if cq is None:
            print(f"   *** UNCLEARED (q>60) blocker at m={m}: {list(W)}"); continue
        if cq > maxclear: maxclear = cq; ex = (list(W), cq)
    if n_block: blocker_mins.append(m)
    overall_maxclear = max(overall_maxclear, maxclear)
    tag = "samp" if sampled else "exh"
    print(f"{m:>6} {f'[{m},{hi}]':>10} {f'{n_block} ({tag})':>20} {maxclear:>10}  {str(ex):>30}")

print("="*84)
print(f"largest min m with a compressed non-translate non-AP blocker found: {max(blocker_mins) if blocker_mins else 'NONE'}")
print(f"overall max clearing-q: {overall_maxclear}")
print()
print("READING: if blockers PERSIST at large min (unbounded height) -> the node is NOT a finite family")
print("check; it needs the height-uniform covering (finite RESIDUE classes). If blockers VANISH above")
print("some m0 -> the node is a FINITE check (native_decide over m<=m0). Either way max-clear-q bounds")
print("the covering set. This decides the formalization strategy.")
print("DONE")
