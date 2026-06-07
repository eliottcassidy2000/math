#!/usr/bin/env python3
"""
unit_distance_3n_crossover_focus_s4.py
monad-explorer-2026-06-07-S4  (focused follow-up to the family sweep)

The family sweep showed sqrt7 (eisen t=7) and sqrt13 (eisen t=13) are the two
single-norm lattice families closest to beating 3N at small N, sitting at
deficit -1 (u=83 vs 84) at N=28 -- WITHIN search noise.  AMP only PROVES
u(n)<=3n up to n=24; for n in {25,26,27} it is OPEN.  So a sqrt7/sqrt13 patch
with u > 3N at N in {25,26,27} would LOWER the OPEN-Q-057 ceiling below 28.

Here we hammer the densest-k-subset for these two families at N=24..34 with a
stronger search (peel-from-large + multi-restart anneal) and EXACTLY certify
the best patch's unit-distance count (re-counted from raw coordinates).
"""
import random
random.seed(424242)

def eisen_offsets(t, rad=16):
    return [(a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)
            if a*a + a*b + b*b == t]

def eisen_Q(p):
    return p[0]*p[0] + p[0]*p[1] + p[1]*p[1]

def deg_in(p, Sset, V):
    return sum((p[0]+v[0], p[1]+v[1]) in Sset for v in V)

def edges(Sset, V):
    return sum(deg_in(p, Sset, V) for p in Sset) // 2

def big_patch(V, size, rad=16):
    cells = sorted(((a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)),
                   key=lambda p: (eisen_Q(p), p))
    return set(cells[:size])

def peel_to(Sset, V, N, jitter=0.0):
    """remove min-degree vertices until size N (random tie-break via jitter)."""
    S = set(Sset)
    while len(S) > N:
        # vertex with fewest internal neighbours leaves
        worst = min(S, key=lambda p: (deg_in(p, S, V) + jitter*random.random(),
                                       eisen_Q(p)))
        S.discard(worst)
    return S

def anneal(Sset, V, Vset, iters):
    S = set(Sset)
    E = edges(S, V)
    best = E; bestS = set(S)
    N = len(S)
    for it in range(iters):
        T = max(0.05, 3.0*(1 - it/iters))
        u = random.choice(tuple(S))
        p = random.choice(tuple(S)); v = random.choice(V)
        w = (p[0]+v[0], p[1]+v[1])
        if w in S or w == u:
            continue
        du = deg_in(u, S, V)
        dw = deg_in(w, S, V) - (1 if (w[0]-u[0], w[1]-u[1]) in Vset else 0)
        delta = dw - du
        if delta >= 0 or random.random() < pow(2.718281828, delta/T):
            S.remove(u); S.add(w); E += delta
            if E > best:
                best = E; bestS = set(S)
    return best, bestS

def densest(V, N, big, iters=120000, restarts=14):
    Vset = set(V)
    best = -1; bestS = None
    # start 1: disk peel (no jitter)
    cands = []
    cands.append(peel_to(big, V, N, jitter=0.0))
    for r in range(restarts):
        cands.append(peel_to(big, V, N, jitter=1.5))
    for c in cands:
        e, S = anneal(c, V, Vset, iters)
        if e > best:
            best, bestS = e, S
    return best, bestS

def exact_recount(S, V):
    """independent exact unit-distance recount from raw coords."""
    Sset = set(S)
    Vset = set(V)
    cnt = 0
    L = list(S)
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            d = (L[i][0]-L[j][0], L[i][1]-L[j][1])
            if d in Vset:
                cnt += 1
    return cnt

def main():
    print("="*70)
    print("FOCUSED densest-k-subset: sqrt7 (t=7) and sqrt13 (t=13) Eisenstein")
    print("target window N=24..34; exact recount; OPEN range for ceiling = {25,26,27}")
    print("="*70)
    for t in (7, 13):
        V = eisen_offsets(t)
        big = big_patch(V, 160)
        print(f"\n### eisen t={t} (degree {len(V)}) ###")
        print(f"   {'N':>3} {'best':>5} {'3N':>5} {'u-3N':>6}  recount  status")
        for N in range(24, 35):
            best, S = densest(V, N, big)
            rc = exact_recount(S, V)
            assert rc == best, f"recount mismatch {rc} vs {best}"
            d = best - 3*N
            status = ""
            if d > 0:
                status = "BEATS 3N"
                if N <= 27:
                    status += "  *** LOWERS CEILING (N*<=%d) ***" % N
            elif d == 0:
                status = "ties 3N"
            print(f"   {N:>3} {best:>5} {3*N:>5} {d:>6}  {rc:>6}   {status}")
    print("\n" + "="*70)
    print("AMP-proven: u(n)<=3n for n<=24. So a >3N here at N=25,26,27 is a NEW")
    print("ceiling; at N=28 it only re-certifies (exact coords) the known N*<=28.")

if __name__ == "__main__":
    main()
