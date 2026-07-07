#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S29 -- toward the FULL crux: does the parity/competing-
denominator mechanism (THM-632) extend from the CANONICAL family to the whole
BORDERED-AP candidate class at N=12?

THM-632: F(12)={1..10,12,33} clears at 3/35 (competitor denominator Q=3N-1=35,
odd) > 2/25 => not a gap member.  The RESIDUAL for the full crux (opus-S118 fleet
sweep empty): show EVERY candidate gap family at N=12 has such a competing
clearance above the gap.

Gap members are BORDERED APs (my S25 / opus-S117): an AP core + doubled endpoints
+ possibly a far resonator.  This script enumerates the bordered-AP candidate
class at N=12 (12 speeds) targeting the gap (1/13, 2/25), and for each records:
  - M and its denominator q_M,
  - whether M is in the gap, above (loose), or below (tight),
  - the PARITY of q_M and whether q_M relates to 3N-1=35 or the odd competitors.
GOAL: is the gap EMPTY over this class, and is the escape ALWAYS via an odd /
factorable competing denominator (the THM-632 mechanism generalized)?
"""
import itertools
from fractions import Fraction as F
from math import gcd

LO, HI = F(1, 13), F(2, 25)  # the second gap


def exact_M_and_q(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0)
    seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    W2 = tuple(sorted(set(v // g for v in W))) if g > 1 else tuple(sorted(set(W)))
    return W2


def bordered_ap_families(N, dmax=6, xmax=None):
    """Bordered-AP candidate families with N speeds:
       AP core {a, a+d, ..., a+(L-1)d} + doubled endpoints (a-d? / +d borders)
       + far resonator(s), all normalized to N distinct primitive speeds."""
    if xmax is None:
        xmax = 3 * N + 4
    fams = set()
    for d in range(1, dmax + 1):
        for a in range(1, d + 2):
            for L in range(max(2, N - 4), N + 1):
                core = [a + i * d for i in range(L)]
                if len(set(core)) != L:
                    continue
                # borders: doubled endpoints (endpoint +-1), classic bordered-AP
                border_opts = [
                    [], [core[0] - 1], [core[-1] + 1], [core[0] - 1, core[-1] + 1],
                    [core[0] + 1], [core[-1] - 1],
                ]
                for borders in border_opts:
                    base = [v for v in core + borders if v >= 1]
                    base = sorted(set(base))
                    need = N - len(base)
                    if need < 0:
                        continue
                    if need == 0:
                        W = primitive(tuple(base))
                        if len(W) == N:
                            fams.add(W)
                    else:
                        # add `need` far resonators near 3(N-1) (the mediant far element)
                        far_center = 3 * (N - 1)
                        cand = [far_center + o for o in range(-4, 5)] + \
                               list(range(max(base) + 1, xmax))
                        cand = [c for c in cand if c not in base and c >= 1]
                        for combo in itertools.combinations(sorted(set(cand)), need):
                            W = primitive(tuple(sorted(set(base) | set(combo))))
                            if len(W) == N:
                                fams.add(W)
    return fams


def main():
    N = 12
    print("=" * 88)
    print(f"BORDERED-AP CANDIDATE CLASS at N={N}: is the gap (1/13,2/25) empty, and")
    print(f"is every escape via an odd/competing denominator (THM-632 generalized)?")
    print("=" * 88)
    fams = bordered_ap_families(N)
    print(f"  enumerated {len(fams)} bordered-AP candidate families (12 speeds each)")
    in_gap = []
    above = 0
    below = 0
    q35 = 0            # escapes at denominator 35 = 3N-1 (the THM-632 competitor)
    odd_q = 0          # escapes at an ODD denominator
    even_q = 0
    above_examples = []
    for W in fams:
        M = exact_M_and_q(W)
        if LO < M < HI:
            in_gap.append((M, W))
        elif M >= HI:
            above += 1
            q = M.denominator
            if q == 35:
                q35 += 1
            if q % 2 == 1:
                odd_q += 1
            else:
                even_q += 1
            if len(above_examples) < 8 and M < F(1, 11):  # just above the gap
                above_examples.append((M, q, W))
        else:
            below += 1
    print(f"  IN GAP (1/13,2/25):  {len(in_gap)}")
    print(f"  ABOVE gap (loose, M >= 2/25): {above}")
    print(f"  BELOW gap (tight, M <= 1/13): {below}")
    print()
    print(f"  Of the {above} above-gap escapes:")
    print(f"    denominator == 35 (=3N-1, THM-632 competitor): {q35}")
    print(f"    ODD denominator:  {odd_q}")
    print(f"    EVEN denominator: {even_q}")
    print()
    if in_gap:
        print(f"  *** {len(in_gap)} IN-GAP families found -- INVESTIGATE (would be gap members!):")
        for M, W in sorted(in_gap)[:10]:
            print(f"      M={M} ({float(M):.5f})  W={list(W)}")
    else:
        print("  => the bordered-AP candidate class at N=12 is EMPTY in the gap.")
    print()
    print("  Escape denominators just above the gap (M in (2/25, 1/11)):")
    for M, q, W in sorted(above_examples):
        print(f"    M={M} ({float(M):.5f}) q={q} ({'odd' if q%2 else 'even'})  W={list(W)[:8]}..")
    print()
    print("  INTERPRETATION: if the class is gap-empty AND escapes concentrate at")
    print("  odd/small competing denominators (esp. 35=3N-1), THM-632's parity")
    print("  mechanism generalizes -- the gap is protected by odd-denominator")
    print("  competitors, a structural REASON for opus-S118's empty sweep at N=12.")


if __name__ == "__main__":
    main()
