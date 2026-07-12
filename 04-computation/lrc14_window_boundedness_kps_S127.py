# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont47: is the clearing window BOUNDED over spread DC families? (the crux of whether
# Route B is a genuinely FINITE check, per klein-S258, or unbounded/diameter-growing).
#
# opus-S241 + my cont.46: clearing at q needs (a) q nmid every v_i AND (b) coprime sub-family misses a
# fold-class. opus-S241: "no fixed bounded window is a shortcut" -- an adversary blocks a modulus by carrying a
# multiple of it. BUT that was general; the DC+primitive+spread constraints may BOUND how much can be blocked.
# THIS resolves it: can a spread DC primitive 13-family block ALL non-14 q in [15,W] (=> min-clear-q > W)?
# Construction: to block prime p in [15,W] we need p|v_i; combine with small factors to stay DC. Push W up.
import random
from math import gcd, lcm
from functools import reduce

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def lrun(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def valid(v): return len(set(v)) == 13 and prim(v) and is_DC(v) and lrun(v) <= 7
def clears_at(v, q):
    lo = -(-q // 14)
    return any(all(lo <= (vi * p) % q <= q - lo for vi in v) for p in range(1, q))
def min_clear_q(v, hi=200):
    for q in range(15, hi):
        if q % 14 and clears_at(v, q): return q
    return None
def blocked(v, W):  # non-14 moduli in [15,W] that some v_i is divisible by
    return [q for q in range(15, W+1) if q % 14 and any(x % q == 0 for x in v)]

def main():
    rng = random.Random(20260714)
    Wtargets = [31, 37, 43, 50, 60]
    print("Q1: can a spread DC primitive family block ALL non-14 q in [15,W]? (=> min-clear-q > W)")
    for W in Wtargets:
        allq = [q for q in range(15, W+1) if q % 14]
        best_unblocked = len(allq); best_v = None
        # annealing on (# unblocked moduli in [15,W]), allowing large speeds
        for restart in range(40):
            vmax = rng.choice([200, 500, 1500, 5000, 20000])
            v = None
            for _ in range(800):
                w = sorted(rng.sample(range(1, vmax+1), 13))
                if valid(w): v = w; break
            if not v: continue
            cur = len([q for q in allq if q not in blocked(v, W)])
            for _ in range(2500):
                w = v[:]; i = rng.randrange(13)
                # bias mutation toward multiples of an unblocked modulus
                unb = [q for q in allq if not any(x % q == 0 for x in w)]
                if unb and rng.random() < 0.6:
                    q0 = rng.choice(unb); mult = rng.randint(1, max(2, vmax//q0)); w[i] = q0*mult
                else:
                    w[i] = rng.randint(1, vmax)
                w = sorted(set(w))
                if len(w) != 13 or not valid(w): continue
                u = len([q for q in allq if not any(x % q == 0 for x in w)])
                if u <= cur:
                    v, cur = w, u
                    if u < best_unblocked: best_unblocked, best_v = u, w
        mcq = min_clear_q(best_v) if best_v else None
        print(f"  W={W}: min #unblocked-in-[15,{W}] found = {best_unblocked} of {len(allq)}; that family min-clear-q = {mcq}")
        if best_unblocked == 0:
            print(f"    *** FULLY BLOCKED [15,{W}] by a spread DC family => min-clear-q > {W} => window UNBOUNDED for DC ***")
        else:
            print(f"    could NOT fully block [15,{W}] (>=1 modulus always free) -- consistent with a bounded window")

    # Q2: relation between blocking-count and Vmax needed (does full blocking require huge diameter?)
    print("\nQ2: max min-clear-q found across all restarts (the empirical window ceiling for DC):")
    ceil = 0; cv = None
    for _ in range(4000):
        vmax = rng.choice([100, 400, 2000, 10000])
        w = sorted(rng.sample(range(1, vmax+1), 13))
        if valid(w):
            q = min_clear_q(w) or 0
            if q > ceil: ceil, cv = q, w
    print(f"  random high-Vmax: max min-clear-q = {ceil}")
    print(f"  => if the ceiling keeps rising with Vmax, the window is UNBOUNDED (Route B is not a fixed finite check);")
    print(f"     if it saturates, there is a bounded window. (Structural: blocking q needs a multiple => diameter >= q.)")

if __name__ == "__main__":
    main()
