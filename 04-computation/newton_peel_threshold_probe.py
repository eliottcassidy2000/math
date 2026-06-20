#!/usr/bin/env python3
"""Probe the TRUE vanishing threshold of Delta_S in the seven-sector p0 model.

The algebraic Newton identity p0(BuF)=sum_S Delta_S is a tautology (Mobius
inversion) -- true for ANY p0. The substantive claim is Delta_S=0 for |S|>6.

This probe asks: what actually drives Delta_S=0? Is 6 the real threshold, or is
the vanishing far stronger (kicking in at |S|>=2 for independent runners) /
or weaker (some |S|=7 nonzero exists that the earlier probe missed)?

Mechanism hypothesis: Delta_S measures the FAILURE of p0 to be additive in the
"survival deficits" of the runners in S. If the survival events of distinct
runners were independent, p0 would be MULTIPLICATIVE, and the second mixed
difference Delta_{a,b} would generally be NONZERO (multiplicative != additive).
So vanishing of high-order Delta is a real structural fact, not automatic.
"""
from fractions import Fraction
from itertools import combinations

HALF = Fraction(1, 14)


def _danger(c):
    if c == 0:
        return [(Fraction(0), Fraction(1))]
    hw = Fraction(1, 14 * c)
    arcs = []
    for a in range(c):
        center = Fraction(a, c)
        lo, hi = center - hw, center + hw
        if lo < 0:
            arcs.append((Fraction(0), hi)); arcs.append((lo + 1, Fraction(1)))
        elif hi > 1:
            arcs.append((lo, Fraction(1))); arcs.append((Fraction(0), hi - 1))
        else:
            arcs.append((lo, hi))
    return arcs


def _merge(iv):
    out = []
    for a, b in sorted(iv):
        if a >= b: continue
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1][1] = b
        else: out.append([a, b])
    return out


def p0(runners):
    danger = []
    for c in runners: danger.extend(_danger(c))
    bad = sum((b - a for a, b in _merge(danger)), Fraction(0))
    return Fraction(1) - bad


def delta_S(B, S):
    S = list(S); tot = Fraction(0)
    for k in range(len(S) + 1):
        sgn = (-1) ** (len(S) - k)
        for T in combinations(S, k):
            tot += sgn * p0(tuple(B) + T)
    return tot


def first_nonzero_size(B, F):
    """Smallest |S|>=1 over subsets S of F with Delta_S != 0; and the largest."""
    F = list(F)
    sizes_nonzero = []
    for k in range(1, len(F) + 1):
        any_nz = any(delta_S(B, S) != 0 for S in combinations(F, k))
        sizes_nonzero.append((k, any_nz))
    return sizes_nonzero


def main():
    print("Threshold probe: at which |S| does Delta_S stop being able to be nonzero?\n")

    cases = [
        ((1, 2, 3, 4, 5, 6), (17, 19, 23, 29, 31, 37, 41, 43)),  # r=8 far
        ((1, 2, 3, 4, 5, 7), (17, 19, 23, 29, 31, 37, 41, 43)),
        ((2,), (17, 19, 23, 29, 31, 37, 41, 43)),                # tiny base
        ((1, 2), (3, 4, 5, 6, 7, 8, 9, 10)),                     # contiguous far
    ]
    for B, F in cases:
        print(f"B={B}, F={F}")
        for k, nz in first_nonzero_size(B, F):
            print(f"  |S|={k}: some Delta_S nonzero? {nz}")
        print()

    # Key structural test: are danger-sets multiplicatively independent so that
    # p0 is MULTIPLICATIVE across runners?  If so, Delta_{a,b} = p0(B)*[...] and
    # high order differences DO vanish only under additivity, not multiplicativity.
    print("Multiplicativity check: does p0({a,b}) == p0({a})*p0({b})?")
    for a, b in [(17, 19), (17, 23), (19, 23), (2, 3), (3, 5)]:
        lhs = p0((a, b)); rhs = p0((a,)) * p0((b,))
        print(f"  a={a},b={b}: p0(a,b)={lhs} ; p0(a)p0(b)={rhs} ; equal? {lhs==rhs}"
              f" ; ratio-1={float(lhs/rhs - 1):.3e}" if rhs else "")
    print()

    # If p0 were multiplicative, what is Delta_{a,b}? For multiplicative f over a
    # base where adding a runner multiplies survival by a factor q_a INDEPENDENT
    # of B, Delta_{a,b}(B) = p0(B)*(q_a - 1)*(q_b - 1) -- generally NONZERO.
    # So vanishing of Delta_{a,b} would mean ADDITIVE survival deficits.
    print("Second mixed difference Delta_{a,b}(B) for far a,b (expect NONZERO):")
    B = (1, 2, 3, 4, 5, 6)
    for a, b in [(17, 19), (17, 23), (19, 23)]:
        d = delta_S(B, (a, b))
        print(f"  Delta_{{{a},{b}}}(B) = {d}  nonzero? {d != 0}")
    print()

    # And third: Delta_{a,b,c} -- the claim allows nonzero up to |S|=6.
    print("Third mixed difference Delta_{a,b,c}(B) (claim: can be nonzero, <=6):")
    for trip in [(17, 19, 23), (17, 19, 29), (19, 23, 29)]:
        d = delta_S(B, trip)
        print(f"  Delta_{trip}(B) = {d}  nonzero? {d != 0}")
    print()

    # Now the decisive question: find the LARGEST |S| for which ANY Delta_S != 0
    # over a rich far-set, to see if the true ceiling is 6 or lower.
    print("LARGEST |S| with a nonzero Delta_S, scanning a rich far set:")
    for B in [(1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 7), (1, 2, 3, 4, 5, 6, 7)]:
        F = (17, 19, 23, 29, 31, 37, 41, 43)
        largest = 0
        for k in range(1, len(F) + 1):
            if any(delta_S(B, S) != 0 for S in combinations(F, k)):
                largest = k
        print(f"  B={B} (|B|={len(B)}): largest nonzero |S| = {largest}")
    print()
    print("Interpretation: largest nonzero |S| reveals the TRUE order; the claim")
    print("says it must be <= 6.  If it is consistently < 6 here, '6' is a safe")
    print("over-estimate (an upper bound), still correct as stated ('Delta_S=0 for |S|>6').")


if __name__ == "__main__":
    main()
