#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_consec_vs_oddAP_IE_macmini.py  (mac-mini 2026-06-21, THREAD D)

PINPOINT what separates consec [1..8] (p0=0.345) from the equal-additive-energy
primitive odd-AP [1,3,5,...,15] (p0=0.309).

p0 = P(all 7 sectors hit) = sum_{j=0}^{7} (-1)^j S_j,  S_j = measure{ v : >= j specified
sectors empty } summed appropriately. The cleanest exact decomposition:
  Let A_t = { v : sector t is EMPTY } = { v : no e in E has sector(e v)=t }.
  p0 = 1 - meas(union_t A_t) = sum_{T subset Z/7} (-1)^|T| meas( intersect_{t in T} A_t ).
We compute S_j = sum_{|T|=j} meas(cap_{t in T} A_t) for j=0..7, EXACT, for both sets.
The signed alternating sum reconstructs p0. Comparing S_j term-by-term shows WHERE
the two sets diverge — i.e. which order of "sectors-simultaneously-empty" correlation
is responsible for consec's advantage.

This is the inclusion-exclusion fingerprint THM-534 / HYP-2726 (Delsarte) acts on:
S_j corresponds to the depth-law / distance distribution. If consec dominates at a
SPECIFIC even j (Krawtchouk band), that pins the Delsarte mechanism concretely.
"""
import itertools
from fractions import Fraction as Fr
from collections import defaultdict

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def empty_pattern_law(E):
    """Return dict: frozenset(empty sectors) -> Lebesgue measure, exact."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E)
    law = defaultdict(lambda: Fr(0))
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        hit = set(sector((e * mid) % 1) for e in E)
        empty = frozenset(set(range(P)) - hit)
        law[empty] += (b - a)
    return law

def S_terms(E):
    """S_j = sum over |T|=j of meas(all sectors in T empty) = sum_{patterns sup T} meas."""
    law = empty_pattern_law(E)
    S = [Fr(0)] * (P + 1)
    for empty, m in law.items():
        # this pattern has 'empty' as the empty set; it contributes to all T subset empty
        e = len(empty)
        for j in range(e + 1):
            # number of T of size j contained in 'empty'
            from math import comb
            S[j] += comb(e, j) * m
    return S

def measS7_from_S(S):
    # p0 = sum_{j} (-1)^j S_j  (inclusion-exclusion on the 7 "sector empty" events)
    return sum((-1)**j * S[j] for j in range(P + 1))

def main():
    print("#"*80)
    print("# IE FINGERPRINT: consec vs equal-AE odd-AP")
    print("#"*80)

    sets = {
        "consec[1..8]":      list(range(1, 9)),
        "oddAP[1,3,..,15]":  [1,3,5,7,9,11,13,15],
    }
    Sdata = {}
    for name, E in sets.items():
        S = S_terms(E)
        Sdata[name] = S
        p0 = measS7_from_S(S)
        print(f"\n{name}:  p0 (reconstructed) = {p0} = {float(p0):.6f}")
        print("  S_j (measure with >= the j named sectors empty):")
        for j in range(P + 1):
            print(f"    S_{j} = {str(S[j]):>14s} = {float(S[j]):.6f}")

    print("\n=== TERM-BY-TERM DIFFERENCE (consec - oddAP), signed contribution (-1)^j*dS_j ===")
    Sc = Sdata["consec[1..8]"]; So = Sdata["oddAP[1,3,..,15]"]
    cum = Fr(0)
    for j in range(P + 1):
        dS = Sc[j] - So[j]
        signed = (-1)**j * dS
        cum += signed
        print(f"  j={j}: dS_j={float(dS):+.6f}  (-1)^j dS_j={float(signed):+.6f}  "
              f"running p0-gap={float(cum):+.6f}")
    print(f"\n  total p0 gap (consec - oddAP) = {float(Sc[0]-So[0] + measS7_from_S(Sc)-measS7_from_S(So) - (Sc[0]-So[0])):.6f}")
    print(f"  (sanity) p0_consec - p0_oddAP = {float(measS7_from_S(Sc)-measS7_from_S(So)):+.6f}")

    # WHICH single j dominates the gap?
    print("\n=== which j carries the consec advantage? ===")
    contribs = [( (-1)**j*(Sc[j]-So[j]), j) for j in range(P+1)]
    for val, j in sorted(contribs, key=lambda t: -float(t[0])):
        print(f"  j={j}: signed gap contribution = {float(val):+.6f}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
