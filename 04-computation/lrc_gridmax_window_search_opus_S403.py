"""
opus-2026-07-19-S403 (H1 / HYP-7975): FIRST DATA ON THE G*0(11) WINDOW.

By THM-1289/HYP-7930-UPDATE, accumulation of the 13-speed spectrum inside
(1/14, 1/13) can only sit AT grid-loneliness values of finite proper subgroups
of the 11-torus.  Cyclic case: gridmax(a; q) = max_j min_i ||j a_i / q||, an
11-tuple a mod q.  A value g/q lies in the open window iff 13g < q < 14g,
and gridmax(a; q) = g/q iff the radius-g dilated intervals of the a_i COVER
Z/q while the radius-(g-1) ones do NOT -- a discrete dilated-interval covering
problem (the danger-comb geometry, discretized).

This script: for g = 2..7 and every q in (13g, 14g), search for realizing
11-tuples via (i) structured seeds (dilates of {1..11}, {1..10,x}, APs),
(ii) randomized tuples, (iii) greedy set-cover with bitmask coverage.
Collect the exact gridmax VALUES found in the window; report per-(g,q)
FOUND/NOT-FOUND-with-effort and the minimal window value found overall
(= the lowest currently-known candidate accumulation site above 1/14).
HONEST: NOT-FOUND is search evidence only (no exhaustion; detection floor
applies); FOUND values are exact certificates (tuple printed).
"""
import random
from math import gcd
from fractions import Fraction

random.seed(20260719)

def masks_for(q, g):
    """mask[a] = bitset of j in [0,q) with ||j*a||_q <= g."""
    ms = [0]*(q)
    for a in range(1, q):
        m = 0
        for j in range(q):
            r = (j*a) % q
            if min(r, q-r) <= g:
                m |= 1 << j
        ms[a] = m
    return ms

def gridmax_exact(tup, q):
    best = 0
    for j in range(q):
        m = min(min((j*a) % q, q-(j*a) % q) for a in tup)
        if m > best: best = m
    return best

def proper(tup, q):
    for j in range(1, q):
        if all((j*a) % q for a in tup): return True
    return False

def try_tuple(tup, q, g, found):
    gm = gridmax_exact(tup, q)
    if gm == g and proper(tup, q):
        found.append((g, q, tuple(sorted(tup))))
        return True
    return False

def search(q, g, budget_random=20000, greedy_restarts=120):
    full = (1 << q) - 1
    ms = masks_for(q, g)
    ms1 = masks_for(q, g-1)
    found = []
    # structured seeds
    for c in range(1, q):
        tup = [(c*i) % q for i in range(1, 12)]
        if 0 in tup: continue
        if try_tuple(tup, q, g, found): return found, 0
    for c in range(1, q):
        base = [(c*i) % q for i in range(1, 11)]
        if 0 in base: continue
        for x in range(1, q):
            tup = base + [x]
            cov = 0
            for a in tup: cov |= ms[a]
            if cov != full: continue
            cov1 = 0
            for a in tup: cov1 |= ms1[a]
            if cov1 == full: continue
            if try_tuple(tup, q, g, found): return found, 0
    # greedy + random
    best_deficit = q
    for _ in range(greedy_restarts):
        tup = []
        cov = 0
        cands = list(range(1, q))
        random.shuffle(cands)
        for _ in range(11):
            bestc, bestgain = None, -1
            for a in cands[:60]:
                gain = bin(ms[a] | cov).count('1')
                if gain > bestgain:
                    bestgain, bestc = gain, a
            tup.append(bestc)
            cov |= ms[bestc]
        deficit = q - bin(cov).count('1')
        best_deficit = min(best_deficit, deficit)
        if deficit == 0:
            cov1 = 0
            for a in tup: cov1 |= ms1[a]
            if cov1 != full and try_tuple(tup, q, g, found):
                return found, 0
    for _ in range(budget_random):
        tup = random.sample(range(1, q), 11)
        cov = 0
        for a in tup:
            cov |= ms[a]
            if cov == full: break
        deficit = q - bin(cov).count('1')
        best_deficit = min(best_deficit, deficit)
        if cov == full:
            cov1 = 0
            for a in tup: cov1 |= ms1[a]
            if cov1 != full and try_tuple(tup, q, g, found):
                return found, 0
    return found, best_deficit

if __name__ == "__main__":
    print("G*0(11)-window search (cyclic subgroups) -- opus-S403")
    window_vals = []
    for g in range(2, 8):
        for q in range(13*g + 1, 14*g):
            found, deficit = search(q, g)
            if found:
                gg, qq, tup = found[0]
                print(f"g={g} q={q}: FOUND gridmax = {gg}/{qq} = "
                      f"{gg/qq:.6f} at {tup}")
                window_vals.append(Fraction(gg, qq))
            else:
                print(f"g={g} q={q}: not found (min radius-{g} cover deficit "
                      f"over search: {deficit} positions)")
    if window_vals:
        print(f"\nwindow values found: {sorted(set(window_vals))}")
        print(f"MINIMAL window value found (lowest known candidate "
              f"accumulation site above 1/14 = {1/14:.6f}): "
              f"{min(window_vals)} = {float(min(window_vals)):.6f}")
    else:
        print("\nNO cyclic G*0(11) values found in (1/14, 1/13) at g <= 7 "
              "(search evidence, not exhaustion)")
