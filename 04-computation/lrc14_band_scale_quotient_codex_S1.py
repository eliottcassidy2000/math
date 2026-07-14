#!/usr/bin/env python3
"""Exact audit of the THM-755 residual under dilation.

Pairwise/Tournament Analysis uses quotient carriers, not runners.  The LRC
predicate is global under dilation, so a runner tournament would discard the
scale action that this audit is testing.  Each carrier declares which labels it
retains; the pairwise observable is the weighted retention advantage.  Two
gauges (predicate-first and compression-first) expose which orientations are
stable, and a fixed carrier order resolves ties into a Hamiltonian path.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
import math

from lrc14_certificates import LAM, good_intervals, is_covering


PI_HI = F(22, 7)
BASE = tuple(range(1, 13))


def good_stats(speeds):
    intervals = good_intervals(speeds)
    return sum(b - a for a, b in intervals), len(intervals)


def clearance(speeds, t):
    return min(min((v * t) % 1, 1 - (v * t) % 1) for v in speeds)


def near_dilate(c):
    return tuple(c * i for i in range(1, 13)) + (13 * c + 1,)


def near_dilate_covering_condition(c):
    # q <= 12 is carried by qc.  Modulus 13 needs 13|c.  Modulus 14 is
    # carried by the block iff gcd(c,14)>1, or by 13c+1 iff c=1 (mod 14).
    return c % 13 == 0 and (math.gcd(c, 14) > 1 or c % 14 == 1)


def tournament(vertices, features, weights, tie_order):
    rank = {v: i for i, v in enumerate(tie_order)}
    value = {
        v: sum(weights[k] * int(features[v][k]) for k in weights)
        for v in vertices
    }
    edge = {}
    for a, b in combinations(vertices, 2):
        delta = value[a] - value[b]
        if delta == 0:
            winner = a if rank[a] < rank[b] else b
        else:
            winner = a if delta > 0 else b
        edge[a, b] = winner

    def beats(a, b):
        key = (a, b) if (a, b) in edge else (b, a)
        return edge[key] == a

    outdeg = {v: sum(beats(v, w) for w in vertices if w != v) for v in vertices}
    cycles = sum(
        beats(a, b) and beats(b, c) and beats(c, a)
        or beats(b, a) and beats(c, b) and beats(a, c)
        for a, b, c in combinations(vertices, 3)
    )
    reach = {a: {b for b in vertices if a != b and beats(a, b)} for a in vertices}
    for k in vertices:
        for a in vertices:
            if k in reach[a]:
                reach[a] |= reach[k]
    unseen = set(vertices)
    sccs = []
    while unseen:
        a = min(unseen)
        comp = {b for b in unseen if b == a or (b in reach[a] and a in reach[b])}
        sccs.append(sorted(comp))
        unseen -= comp
    hamiltonian = [
        p for p in permutations(vertices)
        if all(beats(p[i], p[i + 1]) for i in range(len(p) - 1))
    ]
    return {
        "value": value,
        "edge": edge,
        "score_hist": dict(sorted(Counter(outdeg.values()).items())),
        "cycles3": cycles,
        "sccs": sccs,
        "hamiltonian_count": len(hamiltonian),
        "tie_path": tuple(tie_order),
    }


def main():
    base_measure, base_components = good_stats(BASE)
    print("HYP-6780 exact scale audit")
    print(f"base P={{1..12}}: |G_P|={base_measure}, r_P={base_components}")
    assert base_measure == F(6617, 194040)
    assert base_components == 12

    print("\nDirect dilation checks:")
    for c in (1, 2, 7, 13, 26, 91):
        measure, components = good_stats(tuple(c * v for v in BASE))
        assert measure == base_measure
        assert components == c * base_components
        vstar_ratio = components / (math.pi * float(measure) * c)
        print(
            f"  c={c:3d}: |G_cP|={measure}, r_cP={components:4d}, "
            f"v*(cP)/c={vstar_ratio:.12f}"
        )

    print("\nNear-dilate covering condition (exhaustive c<=500 check):")
    for c in range(1, 501):
        assert is_covering(near_dilate(c)) == near_dilate_covering_condition(c)
    covering_scales = [c for c in range(1, 200) if near_dilate_covering_condition(c)]
    print(f"  first covering scales: {covering_scales[:8]}")
    assert covering_scales[0] == 26

    print("\nUnbounded uncapped covering ray:")
    for c in (26, 91, 169, 182, 364, 1820):
        body = near_dilate(c)
        v = body[-1]
        t = F(c + 1, 13 * c)
        cl = clearance(body, t)
        assert is_covering(body)
        assert cl == F(1, 13)
        # PI_HI > pi.  PI_HI*v*G < r therefore proves pi*v*G < r,
        # equivalently v < v*(cP), without a floating-point decision.
        r = c * base_components
        assert PI_HI * v * base_measure < r
        far_count = sum(x > 14 for x in body)
        print(
            f"  c={c:4d}: max={v:6d}, f={far_count:2d}, "
            f"clearance@((c+1)/(13c))={cl}, top-peel=UNCAPPED"
        )
    assert F(1, 13) < F(97, 1000)
    print("  Hence f>=4 does not imply the sampled M>=0.097 margin.")
    print("  The ray is nevertheless closed exactly: core upper bound and witness give M=1/13.")

    vertices = (
        "raw_set",
        "gcd_normalized",
        "core_shape",
        "shape_plus_residue",
        "certificate_only",
    )
    features = {
        "raw_set": dict(lrc=1, covering=1, primitive=1, band=1, witness=1, orbit=0, finite=0),
        "gcd_normalized": dict(lrc=1, covering=1, primitive=1, band=1, witness=1, orbit=0, finite=0),
        "core_shape": dict(lrc=0, covering=0, primitive=0, band=1, witness=0, orbit=1, finite=1),
        "shape_plus_residue": dict(lrc=1, covering=1, primitive=1, band=1, witness=1, orbit=1, finite=1),
        "certificate_only": dict(lrc=1, covering=0, primitive=0, band=0, witness=1, orbit=1, finite=1),
    }
    gauges = {
        "predicate_first": dict(lrc=5, covering=4, primitive=3, band=3, witness=3, orbit=1, finite=1),
        "compression_first": dict(lrc=4, covering=2, primitive=2, band=3, witness=2, orbit=5, finite=5),
    }
    tie_path = (
        "shape_plus_residue",
        "raw_set",
        "gcd_normalized",
        "core_shape",
        "certificate_only",
    )
    print("\nTournament Analysis (vertices are quotient carriers, not runners):")
    tours = {}
    for name, weights in gauges.items():
        tours[name] = tournament(vertices, features, weights, tie_path)
        q = tours[name]
        print(
            f"  {name}: values={q['value']}; score_hist={q['score_hist']}; "
            f"3cycles={q['cycles3']}; SCCs={q['sccs']}; "
            f"Hamiltonian_paths={q['hamiltonian_count']}"
        )
    flips = sum(
        tours["predicate_first"]["edge"][e] != tours["compression_first"]["edge"][e]
        for e in tours["predicate_first"]["edge"]
    )
    print(f"  edge_flips_between_gauges={flips}")
    print(f"  declared_tie_Hamiltonian_path={tie_path}")
    print(
        "  Verdict: shape+residue is stable first. Core shape alone preserves the capped ratio "
        "but destroys covering/primitivity; raw runners preserve truth but do not collapse the infinite ray."
    )


if __name__ == "__main__":
    main()
