#!/usr/bin/env python3
"""
LRC(14) ANGLE subtorus-relation-lattice -- CRUX DIAGNOSTIC (kps-S5).

The whole reduction rests on:  is there an EFFECTIVE spread bound D(k) s.t.
   mu(E) >= mu_min_bounded(k)   for ALL primitive E with spread > D(k) ?
Equivalently: does mu(E) ever DIP BELOW the bounded-shape minimum at INTERMEDIATE
spread (between the small extremizer and the spread->inf torus limit)?

If mu is "monotone up" past the bounded extremizer, B(k) is a finite check.
If mu DIPS below mu_min_bounded at some intermediate spread, the bounded-shape
minimum is NOT the global infimum and the angle FAILS as stated.

TEST.  For each k, take the bounded extremizer E* (the min-mu bounded shape) and
perturb ONE coordinate outward to ever-larger values, tracking exact mu.  Also do a
RANDOMIZED wide search over spreads up to large caps.  Report the running min and
whether anything undercuts mu_min_bounded.

EXACT Fractions.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from functools import reduce
import random

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

def bounded_min_mu(k, cap):
    best, bestE = None, None
    for combo in combinations(range(1, cap+1), k-1):
        E = (0,)+combo
        if reduce(gcd, E) != 1: continue
        m = mu_exact(list(E))
        if best is None or m < best: best, bestE = m, E
    return best, bestE

if __name__ == "__main__":
    random.seed(11)
    # bounded minima (small caps -- enough to find the known extremizers)
    bm = {}
    for k, cap in [(4,7),(5,9),(6,11),(7,13)]:
        b, E = bounded_min_mu(k, cap)
        bm[k] = (b, E)
        print(f"k={k}: bounded-min mu = {b} = {float(b):.6f} at {E}")
    print()

    # ---- TEST 1: outward perturbation of one coordinate of the extremizer ----
    print("TEST 1: push ONE coord of the bounded extremizer outward; track exact mu.")
    for k in (5, 6, 7):
        bmin, E0 = bm[k]
        E0 = list(E0)
        # replace the largest coordinate with M growing
        base = E0[:-1]
        print(f"  k={k}, base={base}, pushing last coord M:")
        running = bmin
        for M in [E0[-1], E0[-1]+1, E0[-1]+2, 10, 20, 50, 100, 300, 1000, 5000]:
            if M <= base[-1]: continue
            E = base + [M]
            if reduce(gcd, E) != 1:  # keep primitive; if not, skip
                # make primitive by adding 0 already there; gcd should be 1 since 0 in set->gcd over nonzero
                pass
            m = mu_exact(E)
            if m < running: running = m
            flag = "  *** BELOW bounded min ***" if m < bmin else ""
            print(f"     M={M:6d}: mu={float(m):.5f} = {m}{flag}")
        print(f"     running min over this family: {float(running):.5f}")

    # ---- TEST 2: randomized wide search at growing spread caps ----
    print("\nTEST 2: randomized primitive shapes at growing spread caps; running min mu.")
    for k in (5, 6, 7):
        bmin, _ = bm[k]
        gmin = None; gE = None
        for cap in [20, 60, 200, 1000, 8000]:
            local = None; localE = None
            trials = 4000 if k <= 6 else 2500
            for _ in range(trials):
                combo = tuple(sorted(random.sample(range(1, cap+1), k-1)))
                E = (0,)+combo
                if reduce(gcd, E) != 1: continue
                m = mu_exact(list(E))
                if local is None or m < local: local, localE = m, E
            if gmin is None or local < gmin: gmin, gE = local, localE
            below = local < bmin
            print(f"  k={k} cap={cap:5d}: min mu found = {float(local):.5f}  "
                  f"{'<<< UNDERCUTS bounded min' if below else 'ok (>= bounded min)'}")
        print(f"  k={k}: GLOBAL random min = {float(gmin):.6f} = {gmin} at {gE}  "
              f"(bounded min {float(bmin):.6f})  "
              f"{'FAIL' if gmin < bmin else 'PASS'}")
    print("\nDONE.")
