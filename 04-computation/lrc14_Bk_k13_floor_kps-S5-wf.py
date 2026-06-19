#!/usr/bin/env python3
"""
LRC(14) ANGLE subtorus-relation-lattice -- k up to 13 FLOOR (kps-S5).

Two deliverables at the ACTUAL LRC(14) sizes k=3..13:
 (A) the consecutive-AP value  a(k) = mu({0,1,...,k-1})  EXACT  (a natural floor
     candidate and one of the subtorus limits: the all-singletons shared-scale torus).
 (B) the bounded-shape minimum at moderate caps for k=8..13 via a TARGETED search over
     perforated near-APs (full exhaustion is infeasible for k>=9, so we search the
     family the smaller-k extremizers live in: {0,...,k} minus a few interior points,
     plus small one-sided spread), which contains the true min for k<=7.
 (C) re-confirm at k=13 that the dangerous shared-scale AP-of-runs subtori
     (a few representative run-compositions) do not undercut a(13).

This pins the candidate uniform floor  c0 = min_k a(k) or the perforated-min, an
EXPLICIT positive rational, and verifies no tested subtorus type beats it.

EXACT Fractions for mu.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from functools import reduce

G0 = F(2, 7)

def mu_exact(E):
    E = sorted(set(int(e) for e in E)); k = len(E)
    if k == 1: return F(1)
    diffs = {E[i]-E[j] for i in range(k) for j in range(k) if E[i]-E[j] > 0}
    bps = {F(0), F(1)}
    for d in diffs:
        for t in range(0, d+1): bps.add(F(t, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        fr = [F(E[i])*mid - (F(E[i])*mid).__floor__() for i in range(k)]
        order = sorted(range(k), key=lambda i: fr[i])
        n = [(F(E[i])*mid).__floor__() for i in range(k)]
        cross = {a, b}
        for r in range(k):
            i1, i2 = order[r], order[(r+1)%k]; wrap = 1 if r == k-1 else 0
            slope = E[i2]-E[i1]; const = -n[i2]+n[i1]+wrap
            if slope != 0:
                xc = (G0-const)/slope
                if a < xc < b: cross.add(xc)
        cross = sorted(cross)
        for u, v in zip(cross, cross[1:]):
            if u == v: continue
            mm = (u+v)/2
            P = sorted(F(E[i])*mm - n[i] for i in range(k))
            gaps = [P[r+1]-P[r] for r in range(k-1)] + [P[0]+1-P[-1]]
            if max(gaps) > G0: total += (v-u)
    return total

def ap_of_runs(run_sizes, q):
    E = []
    for g, s in enumerate(run_sizes):
        for j in range(s):
            E.append(g*q + j)
    return sorted(set(E))

if __name__ == "__main__":
    print("="*70)
    print("(A) consecutive-AP value a(k) = mu({0,1,...,k-1}), EXACT")
    print("="*70)
    aP = {}
    for k in range(3, 14):
        v = mu_exact(list(range(k)))
        aP[k] = v
        print(f"  k={k:2d}: a(k) = {str(v):>22s} = {float(v):.6f}", flush=True)

    print("\n" + "="*70)
    print("(B) TARGETED perforated near-AP min for k=8..13")
    print("    search {0,1,...,W} choose k, W = k-1 .. k+3 (bounded spread <= k+3)")
    print("="*70)
    pmin = {}
    for k in range(8, 14):
        best, bestE = None, None
        for W in range(k-1, k+4):           # window width; spread bounded by k+3
            for combo in combinations(range(1, W+1), k-1):
                E = (0,)+combo
                if reduce(gcd, E) != 1: continue
                m = mu_exact(list(E))
                if best is None or m < best:
                    best, bestE = m, E
        pmin[k] = best
        below = best < aP[k]
        print(f"  k={k:2d}: perforated-min mu = {str(best):>14s} = {float(best):.6f} "
              f"at {bestE}  {'(< AP!)' if below else '(>= AP)'}", flush=True)

    print("\n" + "="*70)
    print("CANDIDATE UNIFORM FLOOR c0 = min over k=3..13 of the bounded minimum")
    print("="*70)
    floor = None; fk = None
    for k in range(3, 14):
        cand = aP[k]
        if k in pmin and pmin[k] < cand: cand = pmin[k]
        if floor is None or cand < floor: floor, fk = cand, k
    print(f"  c0 (pure mu floor candidate) = {floor} = {float(floor):.6f}  (at k={fk})")
    print(f"  NOTE: pure-mu floor; ROUGH lower bound for the rho*=mu intersect G_P route.")

    print("\n" + "="*70)
    print("(C) k=13 dangerous shared-scale subtori do not undercut a(13)=%s" % str(aP[13]))
    print("="*70)
    a13 = aP[13]
    for comp in [(12,1),(11,2),(7,6),(6,6,1),(1,)*13,(2,)*6+(1,),(11,1,1)]:
        if sum(comp) != 13:
            print(f"  (skip {comp}: sums to {sum(comp)})"); continue
        m97 = mu_exact(ap_of_runs(comp, 97))
        below = m97 < a13
        print(f"  k=13 runs {comp}: mu(q=97) = {float(m97):.5f}  "
              f"{'<<< BELOW a(13)!' if below else 'OK (>= a(13))'}", flush=True)
    print("\nDONE.")
