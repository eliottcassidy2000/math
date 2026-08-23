#!/usr/bin/env python3
"""Exact support-grammar probe for the THM-3795 -> THM-3798 near-connection.

The first two-term nonconstant r-correction is r(g0+g1*e), with Euler
weights {3,0}.  The mandatory minimal nodal core e^2+f0*z contributes
weights {-6,1}.  Against a four-term step-three mate, the actual cell is
therefore 4x4, not the pure 2x4 cell proved in THM-3798.

Incoming, independently audited THM-3803 now closes that affine-linear
carrier.  The probe therefore also checks the next adjacent two-term
correction r(g1*e+g2*e^2), whose weights are {0,-3}.  It remains quadratic
and open at the theorem level.

This probe enumerates every possible scalar bucket and applies only the
proved singleton-bucket sign law.  It is scratch mathematics, not canon.
"""

from collections import defaultdict


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ANCHORS = {"e^2": -6, "f0*z": 1}


def correction(start):
    terms = {}
    for j in (start, start + 1):
        suffix = "" if j == 0 else ("*e" if j == 1 else f"*e^{j}")
        terms[f"r*g{j}{suffix}"] = 3 - 3 * j
    return terms


def carrier(start):
    return {**correction(start), **ANCHORS}


def buckets(components, b):
    out = defaultdict(list)
    for a_name, a_weight in components.items():
        for k in range(4):
            b_weight = b + 3 * k
            out[a_weight + b_weight].append((a_name, a_weight, k, b_weight))
    return dict(sorted(out.items()))


def relative_histogram(components):
    # b=0 records the translation-independent collision grammar.
    return [(weight, len(entries)) for weight, entries in buckets(components, 0).items()]


def singleton_opposite_sign_witness(components, b):
    for total, entries in buckets(components, b).items():
        if len(entries) != 1:
            continue
        a_name, a_weight, k, b_weight = entries[0]
        if a_weight != 0 and b_weight != 0 and a_weight * b_weight < 0:
            return total, a_name, a_weight, k, b_weight
    return None


def probe(start, expected_histogram, expected_offsets):
    components = carrier(start)
    histogram = relative_histogram(components)
    require(histogram == expected_histogram, f"j={start} histogram changed")
    collision_offsets = [q for q, multiplicity in histogram if multiplicity >= 2]
    require(collision_offsets == expected_offsets, f"j={start} unexpected collision offsets {collision_offsets}")

    # A scalar bracket has total Euler weight zero after the bracket shift +2.
    scalar_rows = [(q, -q - 2) for q in collision_offsets]

    witnesses = []
    for q, b in scalar_rows:
        witness = singleton_opposite_sign_witness(components, b)
        require(witness is not None, f"j={start} scalar row q={q}, b={b} survived sign gate")
        total, a_name, a_weight, k, b_weight = witness
        require(total != -2, "hostile singleton accidentally occupies scalar bucket")
        witnesses.append((q, b, total, a_name, a_weight, k, b_weight))

    pure = defaultdict(int)
    for a_weight in correction(start).values():
        for k in range(4):
            pure[a_weight + 3 * k] += 1
    pure_histogram = sorted(pure.items())
    require(histogram != pure_histogram, "anchors failed to change collision grammar")
    return components, histogram, pure_histogram, scalar_rows, witnesses


def print_probe(label, start, result):
    components, histogram, pure_histogram, scalar_rows, witnesses = result
    print(label)
    print(f"r*g correction weights: {sorted(correction(start).values())}")
    print(f"mandatory nodal anchor weights: {sorted(ANCHORS.values())}")
    print(f"actual carrier weights: {sorted(components.values())}")
    print(f"pure projected 2x4 bucket histogram: {pure_histogram}")
    print(f"actual mixed 4x4 bucket histogram: {histogram}")
    print(f"possible scalar rows (relative collision q, mate base b): {scalar_rows}")
    for q, b, total, a_name, a_weight, k, b_weight in witnesses:
        print(
            f"row q={q}, b={b}: singleton hostile total={total}, "
            f"{a_name}[{a_weight}] with mate k={k}[{b_weight}]"
        )


def main():
    affine = probe(
        0,
        [(-6, 1), (-3, 1), (0, 2), (1, 1), (3, 3), (4, 1),
         (6, 2), (7, 1), (9, 2), (10, 1), (12, 1)],
        [0, 3, 6, 9],
    )
    quadratic = probe(
        1,
        [(-6, 1), (-3, 2), (0, 3), (1, 1), (3, 3), (4, 1),
         (6, 2), (7, 1), (9, 1), (10, 1)],
        [-3, 0, 3, 6],
    )

    print("MIXED ANCHORED STEP-THREE SUPPORT AUDIT")
    print_probe("AFFINE-LINEAR CONTROL j=0,1 (now inside THM-3803)", 0, affine)
    print_probe("NEXT LIVE QUADRATIC CONTROL j=1,2", 1, quadratic)
    print("VERDICT: all eight minimal mixed rows fail the singleton opposite-sign gate.")
    print("MAP/LOSS: forgetting e^2 and f0*z recovers the pure step-3 2x4 grammar but does not preserve bracket buckets.")
    print("STATUS FIREWALL: THM-3798 and THM-3803 are audited but do not consume forgotten anchors; THM-3805 remains reserved.")


if __name__ == "__main__":
    main()
