#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 7a (ORIGINAL): Schur-convexity of measS7 in the residue-LEG magnitude vector.

Huffer-Shepp (1987): circle coverage probability is SCHUR-CONVEX in the arc-length vector.
Repo: the 6 nonzero residues r=1..6 each have a "leg" = the minimal-magnitude speed of that
residue (the binding clock). consec = all legs minimal = leg-magnitude vector (1,2,3,4,5,6).
Foster conservation: sum_r g(leg_r mod 7) = 112 (conserved over the residue transversal).

ORIGINAL CLAIM: measS7 (or its complement) is a SCHUR-MONOTONE function of the leg-magnitude
vector L = (|leg_1|,...,|leg_6|). consec has the componentwise-MINIMAL L (1,2,3,4,5,6) -- it
MAJORIZES-from-below every other full-residue config's leg vector. If measS7 is Schur-CONCAVE
decreasing (smaller, more-balanced legs => more coverage), consec is the unique maximizer.

This is a CLEANER object than the c-vector (which failed monotonicity): the LEG vector ignores
the parallel-clock corrections. Test:
  (1) is measS7 a function of L alone? (collisions)  [HYP-2760 says residue-legs determine measS7]
  (2) does consec minorize every L componentwise (L(consec)=(1,..,6) <= L(E))?  -> YES by def of min rep
  (3) MAJORIZATION/Schur test: order L descending; does measS7 increase when L is made MORE
      balanced / smaller in the majorization order? i.e. is measS7 Schur-monotone in L?
  (4) the DECISIVE test: among full-residue configs, is consec (min L) the UNIQUE measS7 max,
      and is the ordering of measS7 consistent with majorization of L?
"""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        if len(s & set(range(1, 7))) == 6: tot += hi - lo
    return tot

def legvec(E):
    leg = {}
    for e in E:
        e = int(e)
        if e == 0: continue
        r = e % 7
        if r == 0: continue
        a = abs(e)
        if r not in leg or a < leg[r]: leg[r] = a
    return tuple(leg.get(r, None) for r in range(1, 7))

def full_residue(E):
    return all(legvec(E)[i] is not None for i in range(6))

def majorizes(a, b):
    """does sorted-desc a majorize b (a more spread)? both same length, same sum assumed-ish."""
    sa = sorted(a, reverse=True); sb = sorted(b, reverse=True)
    pa = pb = 0
    for i in range(len(sa)):
        pa += sa[i]; pb += sb[i]
        if pa < pb - 1e-9: return False
    return True

def rand_full(k, hi, seed, maxtry=500):
    random.seed(seed)
    for _ in range(maxtry):
        E = sorted(set([0] + random.sample(range(1, hi), k - 1)))
        if len(E) == k and reduce(gcd, E) == 1 and full_residue(E):
            return E
    return None

def main():
    print("=== ANGLE 7a: Schur-monotonicity of measS7 in the residue-LEG vector ===\n")

    print("(1) does L determine measS7? (collisions = same L, different measS7)")
    for k in (8, 9):
        seen = {}; coll = 0; tot = 0; exs = []
        for s in range(500):
            E = rand_full(k, 16, 20000 + s)
            if E is None: continue
            L = legvec(E); m = measS7(E); tot += 1
            if L in seen:
                if seen[L][0] != m: coll += 1; exs.append((L, seen[L][1], E))
            else:
                seen[L] = (m, E)
        print(f"  k={k}: {coll} L-collisions w/ different measS7 out of {tot}")
        if exs[:1]: print(f"     ex: L={exs[0][0]} from {exs[0][1]} and {exs[0][2]}")
    print()

    print("(2) consec L = (1,2,3,4,5,6) is componentwise-minimal => every other L >= it.")
    for k in (8, 9):
        consec = list(range(k)); Lc = legvec(consec)
        cw_dom = 0; n = 0
        for s in range(200):
            E = rand_full(k, 16, 22000 + s)
            if E is None: continue
            L = legvec(E); n += 1
            if all(L[i] >= Lc[i] for i in range(6)): cw_dom += 1
        print(f"  k={k} consec L={Lc}: {cw_dom}/{n} configs have L >= consec componentwise")
    print()

    print("(3) Schur test: order configs by measS7; check measS7 ANTI-monotone with L-majorization")
    # collect (measS7, L) pairs; for pairs where L1 majorizes L2 (L1 more spread/larger), is measS7_1 <= measS7_2?
    for k in (8, 9):
        data = []
        for s in range(160):
            E = rand_full(k, 14, 24000 + s)
            if E is None: continue
            data.append((measS7(E), legvec(E), E))
        viol = 0; comp = 0
        for i in range(len(data)):
            for j in range(len(data)):
                if i == j: continue
                m1, L1, _ = data[i]; m2, L2, _ = data[j]
                if majorizes(L1, L2) and not majorizes(L2, L1):
                    comp += 1
                    # L1 more spread => expect LESS coverage: m1 <= m2
                    if m1 > m2 + F(1, 10**6): viol += 1
        print(f"  k={k}: {viol}/{comp} majorization-comparable pairs VIOLATE 'more-spread L => less coverage'")
    print()
    print("(4) DECISIVE: is consec the UNIQUE measS7 max, and L the right order parameter?")
    for k in (8, 9, 10):
        consec = list(range(k)); cm = measS7(consec)
        beat = 0; tie = 0; n = 0
        for s in range(150):
            E = rand_full(k, 14, 26000 + s)
            if E is None: continue
            n += 1; m = measS7(E)
            if m > cm: beat += 1
            elif m == cm and legvec(E) != legvec(consec): tie += 1
        print(f"  k={k} consec measS7={float(cm):.4f}: {beat} beat, {tie} tie(diff L), out of {n}")
    print("\nVERDICT: low majorization-violation + consec unique max => Schur/HS port is the wall proof.")

if __name__ == "__main__":
    main()
