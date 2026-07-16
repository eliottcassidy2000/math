#!/usr/bin/env python3
"""death-star-2026-07-16-S21 (HYP-7027): THE T1540 CUT+CYCLE PROBE + the Moore-bound census.

PROBE A (the movie): for a cluster E, the wall-crossing event stream around the circle:
  states sigma(x) = (floor(7{fx}))_f in Z_7^6 between events; colors = walls (f, w).
  - RETURN SPECTRUM: repeated states; per return-pair: the gap Dx and the palette
    (per-runner odd-crossed wall counts); trivial (palette 0 = silent) vs expressive.
  - CYCLE RANK: movie graph cycle dimension = #events - #states + #components.
  - RELATION SPECTRUM: short vectors of {k in Z^6 : sum k_i f_i = 0} (support <= 4,
    |k_i| <= 8, direct search) — THM-890's coherence side.
  CLAIM UNDER TEST (transference reading): return-richness <-> short relations
  (coherent species recur; incoherent species do not) — cut side (states/returns) vs
  cycle side (palettes/relations) are dual.

PROBE B (2607.14068-inspired census): does COVERING-saturation force shorter additive
  relations among 13-sets? Sample covering vs non-covering primitive 13-sets from the same
  magnitude window; compare shortest-relation statistics (min |k|_1 over nonzero relations
  sum k_i v_i = 0 with support <= 3, coefficients |k| <= 6). The hypergraph Moore bound is
  the density => short-even-cover forcing principle; this asks for its LRC face.
"""
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations, product as iproduct
import random, sys, time

def movie(E):
    """Event stream: sorted section boundaries of all moving runners. Returns
    (states list per cell, events list of (runner index, wall) colors)."""
    mov = [f for f in E if f > 0]
    evs = []
    for i, f in enumerate(mov):
        for w in range(7 * f):
            evs.append((Fr(w, 7 * f), i, w))
    evs.sort()
    # dedupe positions (collisions: multiple runners at once = one event with several colors)
    cells = []   # state per cell
    events = []  # list of sets of (i, w)
    pos = []
    i = 0
    while i < len(evs):
        p = evs[i][0]
        cols = set()
        while i < len(evs) and evs[i][0] == p:
            cols.add((evs[i][1], evs[i][2]))
            i += 1
        pos.append(p); events.append(cols)
    n = len(pos)
    for j in range(n):
        a = pos[j]
        b = pos[(j + 1) % n] if j + 1 < n else Fr(1)
        mid = (a + b) / 2
        st = tuple(int((f * mid % 1) * 7) for f in mov)
        cells.append(st)
    return pos, events, cells, mov

def probe_A(E, label):
    pos, events, cells, mov = movie(E)
    n = len(cells)
    # returns: state -> list of cell indices
    from collections import defaultdict
    where = defaultdict(list)
    for j, st in enumerate(cells):
        where[st].append(j)
    nstates = len(where)
    ncomp = 1  # single closed walk => connected on visited states
    cycrank = n - nstates + ncomp
    # adjacent-return palettes: for states visited >= 2, take consecutive visit pairs
    triv = expr = 0
    best_gap = None
    for st, js in where.items():
        if len(js) < 2: continue
        for a, b in zip(js, js[1:]):
            # palette: per-runner count of crossings in (a, b] mod 2 per wall — summarize
            # as the multiset of walls crossed an odd number of times
            cnt = defaultdict(int)
            for j in range(a + 1, b + 1):
                for (i, w) in events[j % n]:
                    cnt[(i, w)] ^= 1
            odd = [k for k, v in cnt.items() if v]
            if odd: expr += 1
            else: triv += 1
            gap = pos[b] - pos[a]
            if best_gap is None or gap < best_gap: best_gap = gap
    # relation spectrum: support <= 4, |k| <= 8
    rels = []
    idx = range(len(mov))
    for sup in range(2, 5):
        for comb in combinations(idx, sup):
            for ks in iproduct(range(-8, 9), repeat=sup):
                if any(k == 0 for k in ks): continue
                if sum(k * mov[i] for k, i in zip(ks, comb)) == 0:
                    l1 = sum(abs(k) for k in ks)
                    rels.append((l1, sup, tuple(zip(comb, ks))))
        if rels: break
    rels.sort()
    kmin = rels[0][0] if rels else None
    nret = sum(max(0, len(js) - 1) for js in where.values())
    print(f"  {label}: events={n} states={nstates} cycrank={cycrank} "
          f"returns={nret} (silent={triv}, expressive={expr}) "
          f"min-relation-L1={kmin} (found {len(rels)} at min support)")
    return dict(events=n, states=nstates, cycrank=cycrank, returns=nret,
                triv=triv, expr=expr, kmin=kmin)

def shortest_relation_L1(S, supmax=3, kmax=6):
    best = None
    n = len(S)
    for sup in range(2, supmax + 1):
        for comb in combinations(range(n), sup):
            vals = [S[i] for i in comb]
            for ks in iproduct(range(-kmax, kmax + 1), repeat=sup):
                if any(k == 0 for k in ks): continue
                if sum(k * v for k, v in zip(ks, vals)) == 0:
                    l1 = sum(abs(k) for k in ks)
                    if best is None or l1 < best: best = l1
    return best

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def probe_B(trials=400, lo=1, hi=200, seed=20260716):
    rnd = random.Random(seed)
    stats = {True: [], False: []}
    t0 = time.time()
    while (len(stats[True]) < trials or len(stats[False]) < trials) and time.time() - t0 < 240:
        S = sorted(rnd.sample(range(lo, hi + 1), 13))
        g = 0
        for v in S: g = gcd(g, v)
        if g != 1: continue
        cov = is_covering(S)
        if len(stats[cov]) >= trials: continue
        r = shortest_relation_L1(S)
        stats[cov].append(r if r is not None else 99)
    import statistics as st
    for cov in [True, False]:
        xs = stats[cov]
        if xs:
            print(f"  covering={cov}: n={len(xs)} min-relation-L1: mean={st.mean(xs):.3f} "
                  f"median={st.median(xs)} min={min(xs)} frac(L1<=4)={sum(1 for x in xs if x<=4)/len(xs):.3f}")

if __name__ == "__main__":
    t0 = time.time()
    print("PROBE A: the movie's return spectrum vs the relation spectrum (cut vs cycle)")
    probe_A([0, 50, 51, 52, 53, 54, 55], "consecutive c=50 (coherent)")
    probe_A([0, 50, 69, 96, 142, 207, 294], "generic c=50 (incoherent)")
    probe_A([0, 48, 69, 96, 144, 207, 294], "planted 48+96=144 (one relation)")
    probe_A([0, 1, 2, 3, 4, 5, 50], "far bank t=50")
    sys.stdout.flush()
    print("\nPROBE B: does covering-saturation shorten additive relations? (2607.14068 face)")
    probe_B()
    print(f"[total {time.time()-t0:.1f}s]")
