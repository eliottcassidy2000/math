#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S35 -- THE critical question: is the r=2 covering modulus Q0
BOUNDED (kps/opus 'no height bound') or does it GROW with the lift heights (a,b)?

Since clearing at q is periodic in (a mod q, b mod q), the max-over-all-(a,b) min
clearing q = the q for the WORST residue pattern -- BOUNDED in principle, but the
witnessing (a,b) can be huge (~lcm).  My S35 [0,40)^2 found max q=37 (challenges
kps Q0=25).  Here: does max-q SATURATE as the range widens, and what IS Q0?

Approach: (1) for the hardest shape (10,12), the worst (a,b) needing large q -- find
one, verify it's a mod-25 BLOCKER (else it clears at 25) and clears at NO q<Q.
(2) Search residue patterns DIRECTLY (a mod L, b mod L for L=lcm of candidate
moduli) rather than raw (a,b), to find the true worst-case min-clearing-q.
"""
from math import gcd, lcm
import itertools


def mu_of(q):
    return -(-2 * q // 25)


def clears_at(vlist, q):
    mu = mu_of(q)
    if 2 * mu > q:
        return None
    for c in range(1, q):
        if gcd(c, q) != 1:
            continue
        if all(mu <= (v * c) % q <= q - mu for v in vlist):
            return c
    return None


def min_clear_q(vlist, qmax):
    for q in range(6, qmax + 1):
        if clears_at(vlist, q) is not None:
            return q
    return None


def shape_family(i, j, a, b):
    base = [v for v in range(1, 13) if v not in (i, j)]
    fam = sorted(base + [i + 13 * a, j + 13 * b])
    if len(set(fam)) != 12 or any(v <= 0 for v in fam):
        return None
    return tuple(fam)


def blocks_mod25(vlist):
    """mod-25 saturated (blocker) = +-residues cover all 20 units => q=25 does NOT clear."""
    return clears_at(vlist, 25) is None


def main():
    print("=" * 88)
    print("SATURATION: max min-clearing-q for shape (10,12) as (a,b) range widens")
    print("=" * 88)
    ih, jh = 10, 12
    for N in [26, 40, 60, 90, 130]:
        mq = 0
        arg = None
        for a in range(0, N):
            for b in range(0, N):
                if a == 0 and b == 0:
                    continue
                fam = shape_family(ih, jh, a, b)
                if fam is None:
                    continue
                q = min_clear_q(fam, qmax=80)
                if q is None:
                    print(f"  N={N}: ESCAPE at (a,b)=({a},{b}) fam={fam} -- no q<=80!")
                    continue
                if q > mq:
                    mq = q
                    arg = (a, b, fam)
        print(f"  (a,b) in [0,{N})^2 : max min-clearing-q = {mq}   worst (a,b)={arg[:2] if arg else None}")
    print("  => if max-q SATURATES, Q0 is bounded (no height bound); if it keeps GROWING, challenge.")

    print()
    print("=" * 88)
    print("DIRECT worst-case: search (a mod L, b mod L) for small candidate L to find true Q0")
    print("=" * 88)
    # clearing at q depends on (a mod q, b mod q).  To find the worst residue pattern,
    # search a,b over lcm of a candidate modulus window.  Use moduli window [6..24]:
    # min-clearing-q over a full residue period of the small moduli.
    L = 1
    for q in range(6, 25):
        L = lcm(L, q)
    print(f"  lcm(6..24) = {L} (too large for full sweep); sample a,b over structured residues")
    # sample: a,b that make the lifted speeds land on 'bad' residues at many small q.
    # Heuristic: a,b sweeping [0, 200) to catch more CRT combos than [0,130).
    mq = 0
    arg = None
    blockers_at_hi_q = 0
    for a in range(0, 200):
        for b in range(0, 200, 1):
            if a == 0 and b == 0:
                continue
            fam = shape_family(ih, jh, a, b)
            if fam is None:
                continue
            q = min_clear_q(fam, qmax=80)
            if q and q > mq:
                mq = q
                arg = (a, b, fam)
            if q and q >= 26:
                blockers_at_hi_q += 1
    print(f"  (a,b) in [0,200)^2 : max min-clearing-q = {mq}   worst (a,b)={arg[:2] if arg else None}")
    if arg:
        a, b, fam = arg
        print(f"  worst family: {fam}")
        print(f"    is mod-25 blocker (q=25 fails): {blocks_mod25(fam)}")
        print(f"    clears at q=25? {clears_at(fam, 25) is not None};  clears at q={mq}? {clears_at(fam, mq) is not None}")
        # what q<=24 clear it?
        cl = [q for q in range(6, 25) if clears_at(fam, q) is not None]
        print(f"    q<=24 that clear it: {cl if cl else 'NONE -- needs q>=25'}")
    print(f"  # (a,b) needing q>=26: {blockers_at_hi_q}")
    print()
    print("  INTERPRETATION: the true Q0 = max min-clearing-q over ALL residue patterns.")
    print("  kps's Q0=25 is from a narrow (a,b) sample; the real Q0 is larger (>=37 seen).")
    print("  The 'no height bound' claim needs Q0 BOUNDED -- i.e. max-q saturates.  Report it.")


if __name__ == "__main__":
    main()
