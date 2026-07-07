"""
lrc_palindromy_perD_exact_opus_S137.py

EXACT per-diameter palindromy sweeps (kps-S62's named next step; owner worklist).

kps-S62/HYP-4877 + mac-mini-S43/THM-639: gap laws are invariant under step-reversal
(E -> max+min-E), so mu-minimizers come in mirror pairs unless palindromic; every known
exact record family is a palindromic step word (AP = 1^12, prim-sat = 2^5 1^2 2^5,
parity record = 2^4 1^4 2^4).  kps's PALINDROMIC-EXTREMIZER conjecture: per stratum, the
mu_{1/7}-minimizer is palindromic.  Their descents were noise-limited; the S136 exact
order-cell engine removes the noise entirely.

THIS SWEEP (k=13, exact rationals): for each primitive diameter D = 12..19, enumerate ALL
compositions of D into 12 positive steps with gcd 1, up to reversal (= all normalized
13-sets with that diameter, up to mirror); compute exact mu_{1/7}; report the exact
minimizer(s), whether palindromic, the runner-up, and ties.  For D = 20..24, palindromic-
restricted exhaustive + matched-size random non-palindromic sample (exact) as evidence.

Also reports per-D minima of E[maxgap] is NOT attempted here (mean sidecar; kps/klein lane).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random, time, sys

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import mu_exact

def compositions(total, parts):
    """Yield all compositions of `total` into `parts` positive integers."""
    if parts == 1:
        yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest

def steps_to_E(steps):
    E = [1]
    for s in steps:
        E.append(E[-1] + s)
    return E

def is_palin(steps):
    return tuple(steps) == tuple(reversed(steps))

def gcd_all(steps):
    g = 0
    for s in steps: g = gcd(g, s)
    return g

def main():
    K = 13; NSTEPS = 12
    print("=" * 100)
    print("EXACT per-D palindromy sweep, k=13 (all compositions up to reversal, gcd 1)")
    print("=" * 100)
    for D in range(12, 20):
        t0 = time.time()
        seen = set()
        best = None; runner = None; nclass = 0; ties = []
        for steps in compositions(D, NSTEPS):
            if gcd_all(steps) != 1: continue
            canon = min(tuple(steps), tuple(reversed(steps)))
            if canon in seen: continue
            seen.add(canon)
            nclass += 1
            m = mu_exact(steps_to_E(steps))
            if best is None or m < best[0]:
                runner = best
                best = (m, canon)
                ties = [canon]
            elif m == best[0]:
                ties.append(canon)
            elif runner is None or m < runner[0]:
                runner = (m, canon)
        pal = all(is_palin(t) for t in ties)
        print(f"\n D={D:2d}: classes={nclass}")
        print(f"   min mu = {best[0]} = {float(best[0]):.6f} at {ties if len(ties)<=3 else ties[:3]}"
              f"{' (+' + str(len(ties)-3) + ' ties)' if len(ties) > 3 else ''}")
        print(f"   minimizer(s) palindromic: {pal}"
              f"{'   *** NON-PALINDROMIC MINIMIZER ***' if not pal else ''}")
        if runner:
            print(f"   runner-up = {runner[0]} = {float(runner[0]):.6f} at {runner[1]}"
                  f" (palin: {is_palin(runner[1])})")
        print(f"   [{time.time()-t0:.0f}s]")

    print("\n" + "=" * 100)
    print("D = 20..24: palindromic-exhaustive vs random non-palindromic sample (exact)")
    print("=" * 100)
    rng = random.Random(4913)
    for D in range(20, 25):
        t0 = time.time()
        # palindromic words: first 6 steps free, mirrored (12 = 2*6): sum = 2*half
        palbest = None
        if D % 2 == 0:
            half = D // 2
            for h in compositions(half, 6):
                steps = h + tuple(reversed(h))
                if gcd_all(steps) != 1: continue
                m = mu_exact(steps_to_E(steps))
                if palbest is None or m < palbest[0]: palbest = (m, steps)
            note = ""
        else:
            note = " (odd D: no even-length palindrome; skip)"
            palbest = None
        # random non-palindromic sample, exact
        sbest = None
        for _ in range(300):
            cuts = sorted(rng.sample(range(1, D), 11))
            steps = tuple(b - a for a, b in zip([0] + cuts, cuts + [D]))
            if gcd_all(steps) != 1 or is_palin(steps): continue
            m = mu_exact(steps_to_E(steps))
            if sbest is None or m < sbest[0]: sbest = (m, steps)
        if palbest:
            print(f" D={D}: palindromic-exhaustive min = {float(palbest[0]):.6f} at {palbest[1]}")
        else:
            print(f" D={D}:{note}")
        if sbest:
            cmp = ""
            if palbest:
                cmp = ("  [sample >= pal-min OK]" if sbest[0] >= palbest[0]
                       else "  *** RANDOM SAMPLE BEAT PALINDROMIC MIN ***")
            print(f"        random non-pal sample (300) min = {float(sbest[0]):.6f} at {sbest[1]}{cmp}")
        print(f"   [{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
