#!/usr/bin/env python3
"""Exact scale-three transport of the inherited antipodal safe-body certificate.

Universe: three named ten-body Haar hostiles; all 504 distinct positive
ternary-unit triples below24 containing an odd speed; named minimality,
strict endpoint, all-even, and sufficiency-not-necessity controls.
No arbitrary-body existence claim and no silent normalization of tails.
"""

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


SPEC = spec_from_file_location("body", Path(__file__).with_name(
    "lrc14_ten_body_haar_hostile_overnight_hexagon_sep05.py"))
BODY = module_from_spec(SPEC)
SPEC.loader.exec_module(BODY)
CHECKS = 0


def need(test, payload):
    global CHECKS
    CHECKS += 1
    if not test:
        raise AssertionError(payload)


def shift(parts, t):
    out = []
    for lo, hi in parts:
        lo += t
        hi += t
        if lo >= 1:
            lo -= 1
            hi -= 1
        if hi <= 1:
            out.append((lo, hi))
        else:
            out.extend(((lo, Q(1)), (Q(0), hi - 1)))
    return sorted(out)


def paired(parts, t):
    # Returned y obeys y in G and y-t in G; this sign is retained explicitly.
    return BODY.meet(parts, shift(parts, t))


def band_components(C):
    """Independent wall-cell construction of the parity-aware middle band."""
    def good(y):
        for c in C:
            d = min((c * y) % 1, (-c * y) % 1)
            if d < Q(1, 14) or (c % 2 and d > Q(3, 7)):
                return False
        return True

    ends = {Q(0), Q(1)}
    for c in C:
        for k in range(c):
            for sign in (-1, 1):
                y = Q(14 * k + sign, 14 * c) % 1
                ends.add(y)
                if c % 2:
                    ends.add((y + Q(1, 2)) % 1)
    ends = sorted(ends)
    pieces = [(x, x) for x in ends if good(x)]
    pieces += [(lo, hi) for lo, hi in zip(ends, ends[1:]) if good((lo + hi) / 2)]
    out = []
    for lo, hi in sorted(pieces):
        if out and lo <= out[-1][1]:
            out[-1] = out[-1][0], max(hi, out[-1][1])
        else:
            out.append((lo, hi))
    return out


def physical_witness(C, T, y):
    labels = {owner for n, owner in BODY.owners(T, y) if owner is not None}
    if labels == {0, 1, 2}:
        return None
    j = min({0, 1, 2} - labels)
    x = (y + j) / 3
    need(BODY.safe(tuple(3 * c for c in C) + T, x), ("literal physical witness", C, T, x))
    return x


def main():
    positives = (
        ((1, 2, 3, 4, 5, 7, 8, 9, 11, 13), Q(15, 182), Q(8, 1001), Q(21514, 315315)),
        ((1, 2, 3, 5, 7, 8, 9, 11, 12, 13), Q(43, 182), Q(1, 52), Q(14249, 252252)),
        ((1, 3, 4, 10, 11, 13, 14, 16, 17, 18), Q(15, 182), Q(1101, 136136), Q(534689, 7796880)),
    )
    tails = [T for T in combinations((w for w in range(1, 24) if w % 3), 3)
             if any(w % 2 for w in T)]
    need(len(tails) == 504, "complete typed 504-tail bank")
    digest = sha256()
    for C, y, expected_pair_mass, expected_body_mass in positives:
        G = BODY.body_intersection(C)
        need(G == BODY.body_sweep(C), ("independent body construction", C))
        pair = paired(G, Q(1, 2))
        need(pair == band_components(C), ("independent parity-band construction", C))
        need(BODY.mass(pair) == expected_pair_mass, ("pair mass", C))
        need(BODY.mass(G) == expected_body_mass < Q(6, 77), ("actual Haar hostile", C))
        need(BODY.safe(C, y) and BODY.safe(C, y + Q(1, 2)), ("fixed positive packet", C))
        for T in tails:
            odd = next(w for w in T if w % 2)
            distances = [min((odd * z) % 1, (-odd * z) % 1) for z in (y, y + Q(1, 2))]
            need(sum(distances) == Q(1, 2), ("actual odd-coordinate complement", C, T))
            z = (y, y + Q(1, 2))[distances.index(max(distances))]
            need(max(distances) >= Q(1, 4) > Q(3, 14), ("one coordinate inactive", C, T))
            x = physical_witness(C, T, z)
            need(x is not None, ("free scale-three sheet", C, T))
            digest.update((repr((C, T, z, x)) + "\n").encode())
        print("HAAR-HOSTILE UNIFORM-ODD-COORDINATE SUPPLIER", C,
              "PAIR", (y, y + Q(1, 2)), "PAIR MASS", expected_pair_mass,
              "BODY MASS", expected_body_mass)

    minimal = (1, 3, 4, 5, 7, 8, 9, 11, 12)
    need(not paired(BODY.body_intersection(minimal), Q(1, 2)), "nine-core antipodal hostile")
    need(not band_components(minimal), "independent nine-core band hostile")
    deletions = []
    for c in minimal:
        smaller = tuple(v for v in minimal if v != c)
        P = paired(BODY.body_intersection(smaller), Q(1, 2))
        need(P == band_components(smaller) and bool(P), ("inclusion minimality", c))
        deletions.append((c, P[0][0]))
    print("INCLUSION-MINIMAL NINE-CORE NO-HALF-PAIR", minimal)
    print("SINGLE-DELETION RESTORATION PHASES", deletions)

    C = tuple(sorted(minimal + (2,)))
    T = (1, 11, 13)
    G = BODY.body_intersection(C)
    need(BODY.mass(G) == Q(16573, 194040) > Q(6, 77), "Haar passes while fixed pair fails")
    need(all(not paired(G, Q(1, 2 * w)) for w in T), "all canonical per-tail shifts fail")
    need(physical_witness(C, T, Q(1, 14)) == Q(5, 14), "pair test is only sufficient")
    positive_shift = Q(3, 26)
    need(paired(G, positive_shift) == [(Q(3, 14), Q(3, 14)), (Q(82, 91), Q(82, 91))],
         "noncanonical antiphase shift has only singleton intersections")
    y = Q(3, 14)
    need(BODY.safe(C, y) and BODY.safe(C, y - positive_shift), "closed singleton pair")
    need((13 * positive_shift) % 1 == Q(1, 2), "actual antiphase target")
    print("TEN-BODY CANONICAL-SHIFT HOSTILE", C, "T", T, "BODY MASS", BODY.mass(G),
          "PHYSICAL", Q(5, 14), "REPAIRED SHIFT", positive_shift,
          "SINGLETON PAIR", (y, y - positive_shift))

    T = (2, 10, 22)
    y = Q(1, 11)
    need(BODY.spoiled(T, y) and BODY.spoiled(T, y + Q(1, 2)), "all-even tail cannot be silently normalized")
    print("ALL-EVEN HALF-SHIFT HOSTILE", T, "SPOILED BOTH", (y, y + Q(1, 2)))

    C = tuple(range(1, 14))
    T = (1, 5, 11)
    G = BODY.body_intersection(C)
    need(G == [(Q(k, 14), Q(k, 14)) for k in (1, 3, 5, 9, 11, 13)], "exact larger-body singleton control")
    need(all(not paired(G, Q(2 * k + 1, 2 * w)) for w in T for k in range(w)),
         "all antiphase shifts can fail outside ten-body scope")
    need(physical_witness(C, T, Q(1, 14)) == Q(5, 14), "full pair test still not necessary")
    print("OUTSIDE-TEN-BODY FULL-ANTIPHASE NONNECESSITY", C, "T", T,
          "ALL17 SHIFTS FAIL; PHYSICAL", Q(5, 14))
    print("COMPLETE TAILS", len(tails), "POSITIVE PHYSICAL RECOVERIES", 3 * len(tails))
    print("SEMANTIC SHA256", digest.hexdigest())
    print("CHECKS", CHECKS)
    print("PROVED inherited half-shift transport and adaptive antiphase sufficiency; FINITE-EXACT named boundaries")
    print("OPEN arbitrary-row entry; adaptive suppliers remain conditional; no tail-only normalization")


if __name__ == "__main__":
    main()
