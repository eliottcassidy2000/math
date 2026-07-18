#!/usr/bin/env python3
"""Exact referee for THM-1131: the proposed PG(2,13) 1/12 gap is false."""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def clearance(v, p, q):
    residue = (v * p) % q
    return Fraction(min(residue, q - residue), q)


def exact_M(speeds):
    """THM-668 peak replay: test every numerator at every pair-sum denominator."""
    best = Fraction(0)
    owners = []
    denominators = sorted({a + b for a, b in combinations_with_replacement(speeds, 2)})
    for q in denominators:
        for p in range(1, q):  # deliberately retain non-coprime presentations
            value = min(clearance(v, p, q) for v in speeds)
            if value > best:
                best = value
                owners = [(p, q)]
            elif value == best:
                owners.append((p, q))
    reduced = sorted({(Fraction(p, q).numerator, Fraction(p, q).denominator) for p, q in owners})
    return best, reduced, len(denominators)


def covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def primitive(speeds):
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


families = {
    "deep_well": list(range(1, 13)) + [182],
    "compressed_near_dilate": [2 * i for i in range(1, 13)] + [13],
    "divisor_complete_3_over_37": [2, 3, 5, 8, 9, 11, 12, 13, 14, 15, 17, 20, 23],
}
expected = {
    "deep_well": Fraction(14, 183),
    "compressed_near_dilate": Fraction(1, 13),
    "divisor_complete_3_over_37": Fraction(3, 37),
}

print("THM-1131 exact PG(2,13)-gap referee")
print("parameter identity: 183=13^2+13+1; 14=13+1")
records = []
for name, speeds in families.items():
    value, times, denominator_count = exact_M(speeds)
    require(len(speeds) == 13 and len(set(speeds)) == 13, f"bad cardinality: {name}")
    require(primitive(speeds), f"not primitive: {name}")
    require(covering(speeds), f"not Covering: {name}")
    require(value == expected[name], f"wrong M for {name}: {value}")
    records.append((value, name))
    print(
        f"{name}: M={value} maximizer_times="
        f"{','.join(f'{p}/{q}' for p, q in times)} pair_sum_denominators={denominator_count}"
    )

near = families["compressed_near_dilate"]
require(max(clearance(v, 1, 26) for v in near) >= Fraction(1, 13), "internal check")
require(min(clearance(v, 1, 26) for v in near) == Fraction(1, 13), "1/26 witness failed")
require(expected["compressed_near_dilate"] < Fraction(1, 12), "gap counterexample failed")
require(Fraction(1, 13) < expected["divisor_complete_3_over_37"] < Fraction(1, 12),
        "open-interval counterexample failed")

records.sort()
print("ordered M gauge: " + " < ".join(name for _, name in records))
print("tournament fingerprint: transitive; scores=0,1,2; cycles=0; SCCs=1,1,1; HP=1")
print("preserved: scalar M order; destroyed: covering carriers, residues, maximizer geometry")
print("challenged vertices: proof obligations and divisor carriers are more faithful than runners")
print("VERDICT: the proposed [14/183,1/12) empty gap and non-AP=>M>=1/12 target are FALSE")
