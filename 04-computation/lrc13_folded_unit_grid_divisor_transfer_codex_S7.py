#!/usr/bin/env python3
"""Exact finite shell certificate for THM-772.

This script does not search speed sets.  It enumerates the residue/colour
options forced at every unit fraction a/m when a quotient core U has no
m-multiple in THM-769's s=2 or s=3 equality packet.
"""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd, lcm


def units(m):
    return tuple(a for a in range(1, m) if gcd(a, m) == 1)


def inverse_mod(a, m):
    return pow(a % m, -1, m)


def colour_word(s, m, w):
    """Return colours at all unit fractions, or None if ever ineligible."""
    assert gcd(w, s) == 1
    ans = []
    for a in units(m):
        z = a * w
        r = z % m
        if 2 * r < m:
            signed_error = r
        elif 2 * r > m:
            signed_error = r - m
        else:
            # The eligibility radius s/13 is <1/2 for s=2,3.
            return None
        distance = abs(signed_error)
        if 13 * distance > s * m:
            return None
        nearest = (z - signed_error) // m
        colour = (-inverse_mod(w, s) * nearest) % s
        ans.append(colour)
    return tuple(ans)


def local_shell(s, m):
    """Eligible residues modulo sm and all colour-surjective s-packets."""
    eligible = []
    for w in range(s * m):
        if gcd(w, s) != 1:
            continue
        word = colour_word(s, m, w)
        if word is not None:
            eligible.append((w, word))

    packets = []
    target = set(range(s))
    for indices in combinations_with_replacement(range(len(eligible)), s):
        if all(
            {eligible[i][1][column] for i in indices} == target
            for column in range(len(units(m)))
        ):
            packets.append(tuple(eligible[i][0] for i in indices))
    return tuple(eligible), tuple(packets)


def joint_shell(s, moduli):
    """Packets that are eligible/colour-surjective at every listed modulus."""
    period = lcm(*(s * m for m in moduli))
    eligible = []
    for w in range(period):
        if gcd(w, s) != 1:
            continue
        words = tuple(colour_word(s, m, w) for m in moduli)
        if all(word is not None for word in words):
            eligible.append((w, words))

    packets = []
    target = set(range(s))
    for indices in combinations_with_replacement(range(len(eligible)), s):
        good = True
        for modulus_index, m in enumerate(moduli):
            for column in range(len(units(m))):
                if {
                    eligible[i][1][modulus_index][column] for i in indices
                } != target:
                    good = False
                    break
            if not good:
                break
        if good:
            packets.append(tuple(eligible[i][0] for i in indices))
    return period, tuple(eligible), tuple(packets)


def canonical_signed_class(w, modulus):
    return min(w % modulus, (-w) % modulus)


def tournament_fingerprint(vertices, key):
    """Transitive comparison tournament, with vertex order as the tie path."""
    ordered = sorted(range(len(vertices)), key=lambda i: (key(vertices[i]), i))
    rank = {v: i for i, v in enumerate(ordered)}
    edges = {(i, j) if rank[i] < rank[j] else (j, i)
             for i in range(len(vertices)) for j in range(i + 1, len(vertices))}
    scores = [0] * len(vertices)
    for u, _ in edges:
        scores[u] += 1
    # A total-order tournament is transitive: no directed triangles, singleton
    # SCCs, and a unique Hamiltonian path.
    return {
        "edges": edges,
        "order": tuple(vertices[i][0] for i in ordered),
        "score_hist": tuple(sorted(scores)),
        "directed_3cycles": 0,
        "scc_sizes": (1,) * len(vertices),
        "hamiltonian_paths": 1,
    }


def main():
    shells = {}
    print("THM-772 folded unit-grid divisor-transfer certificate")
    print("strict core-safe / closed exception-danger convention")
    print()

    for s in (2, 3):
        print(f"s={s} local shells")
        for m in range(2, 13):
            eligible, packets = local_shell(s, m)
            shells[s, m] = (eligible, packets)
            print(
                f"m={m:2d} units={units(m)} "
                f"eligible={len(eligible):2d} packets={len(packets):2d} "
                f"residues={tuple(w for w, _ in eligible)}"
            )
        print()

    # Two-sheet conclusion: every missing divisor gives an impossible local
    # two-colour packet.
    assert all(len(shells[2, m][1]) == 0 for m in range(2, 13))
    print("s=2 conclusion: U must contain a multiple of every m=2,...,12")

    # Three-sheet local exceptions are exactly m=6 and m=12.
    feasible_s3 = tuple(m for m in range(2, 13) if shells[3, m][1])
    assert feasible_s3 == (6, 12)
    period, joint_eligible, joint_packets = joint_shell(3, (6, 12))
    assert period == 36
    assert joint_eligible == ()
    assert joint_packets == ()
    print("s=3 locally feasible missing divisors:", feasible_s3)
    print("s=3 joint m=6,12 shell: period=36 eligible=0 packets=0")
    print("therefore U must contain a multiple of 6 (no-6 also means no-12)")

    packets12 = shells[3, 12][1]
    signed_packets12 = {
        tuple(sorted(canonical_signed_class(w, 36) for w in packet))
        for packet in packets12
    }
    assert signed_packets12 == {(2, 10, 14)}
    assert len(packets12) == 8
    required_s3 = tuple(m for m in range(2, 12))
    assert all(m == 6 or not shells[3, m][1] for m in required_s3)
    print("s=3 conclusion: U must contain a multiple of every m=2,...,11")
    print("if U has no 12-multiple, exception shell modulo 36 is")
    print("  one independent sign choice from each of +/-2, +/-10, +/-14")
    print("  exact labelled packets:", packets12)
    print()

    # Safe-interval/tooth-width ratios.  For s exceptions over a (12-s)-speed
    # core, LRC gives mu>=1/(13-s).  Comparing the guaranteed safe interval
    # with one eligibility tooth yields w <= (13-s) max(U).
    for s in (2, 3):
        mu_floor = Fraction(1, 13 - s)
        delta = Fraction(1, 13)
        safe_length_per_inverse_b = 2 * (mu_floor - delta)
        tooth_length_per_inverse_w = Fraction(2 * s, 13)
        ratio = tooth_length_per_inverse_w / safe_length_per_inverse_b
        assert ratio == 13 - s
        print(
            f"s={s} safe-interval ratio: "
            f"L*b>={safe_length_per_inverse_b}, tooth*w={tooth_length_per_inverse_w}, "
            f"hence every exception w<={(13-s)}*max(U)"
        )
    print()

    # Tournament Analysis on modulus obligations.  Pairwise observable:
    # constraint strength from (packet count, eligible count).  Gauge A uses
    # raw counts; gauge B normalizes eligible count by the residue period.
    # The fixed (s,m) order is the tie Hamiltonian path.
    vertices = []
    for s in (2, 3):
        for m in range(2, 13):
            eligible, packets = shells[s, m]
            vertices.append((f"s{s}m{m}", s, m, len(eligible), len(packets)))

    raw = tournament_fingerprint(
        vertices, lambda v: (v[4], v[3])
    )
    normalized = tournament_fingerprint(
        vertices, lambda v: (Fraction(v[4], 1), Fraction(v[3], v[1] * v[2]))
    )
    flips = len(raw["edges"] ^ normalized["edges"]) // 2
    print("Tournament Analysis (vertices = modulus obligations)")
    print("pairwise observable = local obstruction strength")
    print("gauge A = raw (packet_count, eligible_count)")
    print("gauge B = (packet_count, eligible_density)")
    print("raw order:", raw["order"])
    print("normalized order:", normalized["order"])
    print("score histogram:", raw["score_hist"])
    print("directed 3-cycles:", raw["directed_3cycles"])
    print("SCC sizes:", raw["scc_sizes"])
    print("Hamiltonian-path counts:", raw["hamiltonian_paths"], normalized["hamiltonian_paths"])
    print("edge flips between gauges:", flips)
    print("carrier warning: this tournament loses unit-fraction columns and the m=6/12 joint splice")
    print("predicate-preserving carrier: (modulus, unit fraction, sheet colour) obligation hypergraph")
    print()
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
