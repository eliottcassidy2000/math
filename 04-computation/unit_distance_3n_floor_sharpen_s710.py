#!/usr/bin/env python3
"""
unit_distance_3n_floor_sharpen_s710.py
monad-explorer-2026-06-07-S710

Two jobs, all EXACT integer arithmetic:

(1) SHARPEN THM-421's [17,32] gap for N* = smallest N with u(N) > 3N.
    Uses the proven values of the Erdos unit-distance maximum u(n) from
    Alexeev-Mixon-Parshall (arXiv:2412.11914, "The Erdos unit distance problem
    for small point sets", Thm 1):
       exact   u(n), n<=21
       bounds  u(22)<=61, u(23)<=66, u(24)<=72  (proven upper bounds)
               u(28)>=85 (Engel/Schade realizable lower bound)
    => u(n) <= 3n PROVEN for all n<=24 ; u(28)>3*28 PROVEN.  N* in [25,28].

(2) INDEPENDENT exact-integer search in the sqrt(7) Eisenstein unit-distance
    graph (the construction lane of THM-421): for each k, find the densest
    k-subset we can (greedy disk + multistart swap hill-climb) and report the
    smallest k for which that family beats 3k. This is a LOWER-bound witness
    machine; it certifies u(k) >= (count) for the k it finds.

Eisenstein lattice Z[w], w = e^{i pi/3}. Point (a,b) = a + b*w.
Squared length form  Q(a,b) = a^2 + a*b + b^2  (exact integer).
"unit^2 = 7" graph: edges are pairs with Q(difference) = 7 (12 vectors).
"unit^2 = 1" graph: the triangular lattice / penny graph (6 vectors).
"""
import itertools, sys

def Q(a, b):
    return a*a + a*b + b*b

def shell(R, rad=8):
    """all integer vectors with Q = R, |coords| <= rad"""
    out = []
    for a in range(-rad, rad+1):
        for b in range(-rad, rad+1):
            if Q(a, b) == R:
                out.append((a, b))
    return out

def box(rad):
    return [(a, b) for a in range(-rad, rad+1) for b in range(-rad, rad+1)]

def count_edges(pts, R):
    """exact count of pairs in pts at squared-distance R"""
    S = set(pts)
    sh = shell(R)
    c = 0
    for (a, b) in pts:
        for (da, db) in sh:
            n = (a+da, b+db)
            if n in S and (n > (a, b)):   # unordered, dedup
                c += 1
    return c

def degree_in(pt, S, sh):
    a, b = pt
    return sum(1 for (da, db) in sh if (a+da, b+db) in S)

def best_ksubset(R, k, centers, rad=8, restarts=40, seed=0):
    """greedy disk around each center + swap hill-climb; returns (best_count,best_set)"""
    sh = shell(R, rad)
    allpts = box(rad)
    # frozen pseudo-random via linear congruential (Date/random banned in workflow,
    # but this is a plain script: use a fixed deterministic shuffle per restart)
    best_c, best_set = -1, None
    # candidate seeds: disk by squared dist to each center
    seeds = []
    for cx, cy in centers:
        ordered = sorted(allpts, key=lambda p: ( (2*p[0]-2*cx)**2 + (2*p[0]-2*cx)*(2*p[1]-2*cy)
                                                 + (2*p[1]-2*cy)**2, p))
        seeds.append(ordered[:k])
    # deterministic perturbed seeds
    for r in range(restarts):
        base = seeds[r % len(seeds)]
        cur = set(base)
        # hill climb: try swapping each in-point for a neighbor-rich out-point
        improved = True
        guard = 0
        while improved and guard < 2000:
            guard += 1
            improved = False
            cur_c = sum(degree_in(p, cur, sh) for p in cur)//2
            # outside candidates = points adjacent to current set
            outside = set()
            for p in cur:
                a, b = p
                for (da, db) in sh:
                    q = (a+da, b+db)
                    if q not in cur:
                        outside.add(q)
            best_swap = None
            for pin in list(cur):
                for pout in outside:
                    new = set(cur); new.discard(pin); new.add(pout)
                    nc = sum(degree_in(p, new, sh) for p in new)//2
                    if nc > cur_c + (0 if best_swap is None else best_swap[0]-cur_c):
                        best_swap = (nc, pin, pout)
            if best_swap and best_swap[0] > cur_c:
                _, pin, pout = best_swap
                cur.discard(pin); cur.add(pout)
                improved = True
        c = count_edges(list(cur), R)
        if c > best_c:
            best_c, best_set = c, set(cur)
    return best_c, best_set


def main():
    print("="*70)
    print("PART 1: N* in [25,28]  (AMP 2412.11914 + Engel/Schade)")
    print("="*70)
    exact = {16:41,17:43,18:46,19:50,20:54,21:57}
    bnd = {22:(60,61),23:(64,66),24:(68,72),25:(72,78),26:(76,84),
           27:(81,90),28:(85,96),29:(89,103),30:(93,110)}
    lo = dict(exact); up = dict(exact)
    for n,(l,u) in bnd.items(): lo[n]=l; up[n]=u
    print(f"{'n':>3} {'lo':>4} {'up':>4} {'3n':>4}  verdict (vs strictly beating 3N)")
    for n in range(16,31):
        t=3*n
        v = "cannot beat (up<=3n)" if up[n]<=t else ("BEATS (lo>3n)" if lo[n]>t else "undetermined")
        print(f"{n:>3} {lo[n]:>4} {up[n]:>4} {t:>4}  {v}")
    floor = max(n for n in range(16,31) if up[n] <= 3*n) + 1
    ceil  = min(n for n in range(16,31) if lo[n] > 3*n)
    print(f"\n  PROVEN: u(n)<=3n for all n<=24  => N* >= {floor}")
    print(f"  PROVEN: u(28)>3*28              => N* <= {ceil}")
    print(f"  ==> N* in [{floor},{ceil}]   (was THM-421: [17,32])")

    print()
    print("="*70)
    print("PART 2: sqrt(7) Eisenstein construction lane (exact integer counts)")
    print("="*70)
    centers = [(0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(1,2),(2,1),(2,2),(3,1)]
    # off-lattice centers (use doubled coords) handled crudely via lattice centers here
    print(f"{'k':>3} {'bestU(sqrt7)':>12} {'3k':>4} {'U-3k':>5}  beats?")
    first_beat = None
    for k in range(13, 40):
        c, st = best_ksubset(7, k, centers, rad=7, restarts=len(centers))
        beat = c > 3*k
        if beat and first_beat is None:
            first_beat = (k, c, st)
        print(f"{k:>3} {c:>12} {3*k:>4} {c-3*k:>5}  {'YES' if beat else ''}")
    if first_beat:
        k,c,st = first_beat
        print(f"\n  sqrt(7) lane first beats 3N at k={k}: U={c} > {3*k}")
        print(f"  witness vertex set (Eisenstein (a,b)): {sorted(st)}")
        # exact re-verify
        print(f"  EXACT recount Q=7 edges: {count_edges(list(st),7)}")
    else:
        print("\n  sqrt(7) lane did not beat 3N in tested range")

if __name__ == "__main__":
    main()
