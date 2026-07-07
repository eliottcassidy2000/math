#!/usr/bin/env python3
r"""
lrc_intersection_ledger_loo27_kps_S60.py   (kind-pasteur-2026-07-07-S60, HYP-4847 part 2)

THE 2/7 LEAVE-ONE-OUT INTERSECTION LEDGER: extends the S60 intersection ledger past
raw bounded diameter to "k-1 tame + one arbitrarily far element" cluster shapes.

POINTWISE CHAIN (all one-line, rigorous):
  (i)  bisection (klein-S154 C2 mechanism): removing one point can at most bisect a gap,
       so maxgap_E(x) >= maxgap_{E\e_j}(x) / 2 for every x
       ==> {maxgap_{E\e_j} > 2/7} ⊆ {maxgap_E > 1/7}.
  (ii) subset lemma (kps-S59) on E\e_j at theta = 2/7:
       {roof_{D_j+1} > 2/7} ⊆ {maxgap_{E\e_j} > 2/7},  D_j = diam(E \ e_j) (raw).
  Intersecting with G_P:
       rho*_{1/7}(P,E) >= meas(G_P ∩ {roof_{D_j+1} > 2/7}) =: ILedger27(P, D_j+1).

So a shape (P,E) clears the hlarge bar whenever SOME leave-one-out of the cluster has
diameter within the 2/7-bite: coverage of far-element clusters at EVERY leg, exactly.
(P = ∅ case = klein-S154 part 3's n2* = 37 crossing; here the G_P-intersected legs.)

OUTPUT: per |P| = 0..5, the largest n with min over all P of ILedger27(P,n) >= m_P,
plus the composite two-level coverage table (raw-diam bite ⊕ LOO-2/7 bite).
"""
from fractions import Fraction as F
from itertools import combinations

TH27 = F(2, 7)
MP = F(14249, 252252)

def merge(iv):
    iv = sorted((a, b) for a, b in iv if b > a)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return [(a, b) for a, b in out]

def complement01(iv):
    iv = merge(iv); out = []; cur = F(0)
    for a, b in iv:
        if a > cur: out.append((cur, a))
        cur = max(cur, b)
    if cur < 1: out.append((cur, F(1)))
    return out

def inter_measure(A, B):
    tot = F(0); i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: tot += hi - lo
        if b1 < b2: i += 1
        else: j += 1
    return tot

def gp_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            c = F(j, p)
            bad.append((max(F(0), c - w), min(F(1), c + w)))
    return complement01(bad)

_farey_cache = {}
def farey(n):
    if n not in _farey_cache:
        fr = set()
        for q in range(1, n + 1):
            for p in range(0, q + 1):
                fr.add(F(p, q))
        _farey_cache[n] = sorted(fr)
    return _farey_cache[n]

_roof_cache = {}
def roof_superlevel(n, theta):
    key = (n, theta)
    if key not in _roof_cache:
        Fs = farey(n); out = []
        for a, b in zip(Fs[:-1], Fs[1:]):
            q, qp = a.denominator, b.denominator
            vl, vr = F(1, q), F(1, qp)
            if vl <= theta and vr <= theta: continue
            if vl > theta and vr > theta:
                out.append((a, b)); continue
            xc = a + (theta - vl) * (b - a) / (vr - vl)
            out.append((a, xc) if vl > theta else (xc, b))
        _roof_cache[key] = merge(out)
    return _roof_cache[key]

print("=" * 96)
print("2/7 LOO INTERSECTION LEDGER: min_P meas(G_P cap {roof_n > 2/7}) vs m_P")
print("=" * 96)
print(f"  m_P = {float(MP):.6f};  sanity P=empty crossing should be 37 (klein-S154 part 3)")
composite = {}
for s in range(0, 6):
    k = 13 - s
    Ps = list(combinations(range(1, 14), s)) if s else [()]
    n = max(k - 1, 3)          # LOO subset has k-1 points => diam >= k-2 => n >= k-1
    last_ok = None; firstfail = None
    while n <= 45:
        worst = None; worstP = None
        for P in Ps:
            v = inter_measure(gp_intervals(P), roof_superlevel(n, TH27))
            if worst is None or v < worst:
                worst, worstP = v, P
        if worst >= MP:
            last_ok = n; n += 1
        else:
            firstfail = (worstP, worst); break
    composite[s] = last_ok
    ff = f"first-fail P={firstfail[0]} val={float(firstfail[1]):.5f}" if firstfail else "no fail through n=45"
    if last_ok:
        print(f"  |P|={s} (k={k:2d}): ILedger27 >= m_P through n = {last_ok}  =>  LOO-bite: some 12-of-{k} subset diam <= {last_ok-1}   [{ff}]")
    else:
        print(f"  |P|={s} (k={k:2d}): NO 2/7 bite (fails already at n={max(k-1,3)}: {ff})")

print()
print("=" * 96)
print("COMPOSITE TWO-LEVEL COVERAGE (raw-diam bite [S60 part 1]  OR  LOO-2/7 bite [this part])")
print("=" * 96)
raw_bites = {0: 75, 1: 34, 2: 21, 3: 17, 4: 11, 5: 11}   # k=13..8 (k=13 primitive)
for s in range(0, 6):
    k = 13 - s
    lo = composite[s]
    lo_txt = f"OR some (k-1)-subset diam <= {lo-1}" if lo else "(no LOO bite)"
    prim = " (primitive)" if s == 0 else " (raw)"
    print(f"  k={k:2d}: cluster diam <= {raw_bites[s]}{prim}  {lo_txt}   ==> G2 >= m_P")
print()
print("  Residual per leg: diam > raw-bite AND every leave-one-out diam > LOO-bite")
print("  (= genuinely two-level-spread clusters; THM-527-D says spread raises mu).")
print("DONE.")
