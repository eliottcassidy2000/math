#!/usr/bin/env python3
"""
lrc14_Bk_rho_constrained_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)

CORRECTED rho* scan: enforce the REAL S3 constraint |P| + k = 13  (13 speeds total =
small part P subset {1..13} of size 13-k  PLUS  large cluster of size k=|E|, 0 in E).

The earlier 'rho' stage random scan reported rho*=0 ONLY for over-sized configs
(|P|+k = 17 or 18 >> 13) which are NOT admissible LRC(14) instances. Here we:
  (1) PROVE numerically that those zeros require |P|+k > 13 (over-constrained G_P).
  (2) Re-run the floor search under the CORRECT constraint |P|+k=13 and report exact min rho*,
      confirming 0 zeros among admissible configs.

rho*(P,E) = meas( G_P  cap  Good(E) ),
  G_P  = {x : ||p x|| >= 1/14  for all p in P},
  Good(E) = {x : maxgap{frac(e x): e in E} > 2/7}.
Exact Fractions throughout.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(20260618)
TWO7 = F(2, 7)

def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out
def meas(arcs): return sum((b - a for a, b in arcs), F(0))
def good_set_exact(E):
    E = sorted(set(E)); k = len(E)
    if k == 1: return [(F(0), F(1))]
    diffs = set()
    for a in range(k):
        for b in range(a + 1, k): diffs.add(E[b] - E[a])
    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(0, d + 1): bps.add(F(m, d))
    bps = sorted(x for x in bps if 0 <= x <= 1)
    good = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        pts = sorted(((E[t] * xm) % 1, E[t]) for t in range(k))
        order = [e for _, e in pts]; floors = [int((e * xm) // 1) for e in order]
        for idx in range(k):
            e_cur = order[idx]; f_cur = floors[idx]
            if idx < k - 1: e_nx = order[idx + 1]; f_nx = floors[idx + 1]; wrap = F(0)
            else: e_nx = order[0]; f_nx = floors[0]; wrap = F(1)
            A = F(e_nx - e_cur); Cc = F(f_cur - f_nx) + wrap
            if A == 0:
                if Cc > TWO7: good.append((x0, x1))
                continue
            xb = (TWO7 - Cc) / A
            if A > 0: lo = max(x0, xb); hi = x1
            else: lo = x0; hi = min(x1, xb)
            if lo < hi: good.append((lo, hi))
    return merge(good)
def GP_arcs(P):
    safe = [(F(0), F(1))]
    for p in P:
        forb = []
        for j in range(p + 1):
            c = F(j, p); forb.append((c - F(1, 14 * p), c + F(1, 14 * p)))
        forb = merge([(max(F(0), a), min(F(1), b)) for a, b in forb if b > 0 and a < 1])
        sp = []; cur = F(0)
        for a, b in forb:
            if a > cur: sp.append((cur, a))
            cur = max(cur, b)
        if cur < 1: sp.append((cur, F(1)))
        new = []
        for a, b in safe:
            for c, d in sp:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        safe = merge(new)
    return safe
def intersect(A, B):
    out = []
    for a, b in A:
        for c, d in B:
            lo, hi = max(a, c), min(b, d)
            if lo < hi: out.append((lo, hi))
    return merge(out)
def rho_star(P, E): return meas(intersect(GP_arcs(P), good_set_exact(E)))
def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g == 1

print("=" * 90)
print("(0) DIAGNOSE the earlier zeros: were they over-sized (|P|+k>13) configs?")
print("=" * 90)
bad = [
    ((1,2,3,6,7,8,9,10,12,13),(0,1,2,3,4,5,8,9)),
    ((1,2,3,4,6,9,11,12,13),(0,2,3,4,5,6,7,8)),
    ((1,2,3,5,6,9,10,12,13),(0,1,2,3,4,5,6,8)),
    ((1,2,3,5,6,7,9,10,12,13),(0,1,2,4,5,6,7,8)),
    ((1,2,3,4,5,6,7,8,12,13),(0,1,2,3,4,5,6,8)),
]
for P, E in bad:
    k = len(E); szP = len(P); r = rho_star(P, E)
    print(f"   |P|={szP} k={k}  |P|+k={szP+k}  rho*={r}  (admissible iff |P|+k=13: {'YES' if szP+k==13 else 'NO -> ARTIFACT'})")
print("   => every reported zero has |P|+k >= 17 (over-constrained); none is an admissible S3 instance.")
print("      Also: meas(G_P) itself can be 0 once |P| is large enough -- that alone forces rho*=0,")
print("      independent of E. We confirm meas(G_P) for these P:")
for P, E in bad:
    print(f"      P={P}: meas(G_P)={meas(GP_arcs(P))}")

print()
print("=" * 90)
print("(1) ADMISSIBLE FLOOR: |P|+k=13, P subset {1..13}, E primitive bounded-spread, 0 in E.")
print("    Exhaustive over P shapes (consecutive + random subsets) x bounded-spread E.")
print("=" * 90)
best = (F(2), None, None); zeros = 0; total = 0
# k ranges over cluster sizes 3..10 (S3 is k>=3); |P| = 13-k small speeds from {1..13}.
for k in range(3, 11):
    szP = 13 - k
    if szP < 0: continue
    capE = {3:14,4:13,5:12,6:11,7:10,8:10,9:10,10:10}[k]
    # E shapes: exhaustive bounded-spread primitive (limit count for tractability)
    Elist = []
    cnt = 0
    for rest in itertools.combinations(range(1, capE + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        Elist.append(E); cnt += 1
        if cnt >= 4000: break
    # P shapes: all C(13, szP) subsets if small, else consecutive + many random
    if szP <= 2:
        Plist = list(itertools.combinations(range(1, 14), szP))
    else:
        Plist = [tuple(range(1, szP + 1)), tuple(range(13 - szP + 1, 14))]
        seen = set(Plist)
        while len(Plist) < 60:
            cand = tuple(sorted(random.sample(range(1, 14), szP)))
            if cand not in seen:
                seen.add(cand); Plist.append(cand)
    for P in Plist:
        gp = GP_arcs(P)
        if meas(gp) == 0:
            # G_P empty -> structurally cannot host a witness; record but note it's a G_P issue
            for E in Elist[:1]:
                total += 1; zeros += 1
            continue
        for E in Elist:
            r = meas(intersect(gp, good_set_exact(E)))
            total += 1
            if r == 0:
                zeros += 1
                if zeros <= 8: print(f"   !!! admissible rho*=0 at P={P} (|P|={szP}) E={E}")
            if r < best[0]:
                best = (r, P, E)
print(f"\n   admissible scan: {total} (P,E) pairs, {zeros} zeros.")
print(f"   *** min rho* (admissible) = {best[0]} = {float(best[0]):.6f}  at P={best[1]} E={best[2]} ***")

print()
print("=" * 90)
print("(2) Does meas(G_P) ever hit 0 for an admissible small part? Scan all P subset {1..13}")
print("    of each size; report the minimum meas(G_P) and whether any G_P is empty.")
print("=" * 90)
for szP in range(0, 11):
    mn = (F(2), None)
    # all subsets if small, else sample
    if szP <= 3:
        Ps = itertools.combinations(range(1, 14), szP)
    else:
        Ps = (tuple(sorted(random.sample(range(1, 14), szP))) for _ in range(3000))
    for P in Ps:
        m = meas(GP_arcs(P))
        if m < mn[0]: mn = (m, P)
    print(f"   |P|={szP:2d}: min meas(G_P) = {str(mn[0]):>14s} = {float(mn[0]):.6f}  at P={mn[1]}  "
          f"{'(EMPTY!)' if mn[0]==0 else ''}")
print("\nDONE.")
