#!/usr/bin/env python3
"""Exact signed-dilation audit for the p=2q-1 three-tooth envelope.

For q>=2 put p=2q-1, delta=1/(4q), and

    B_-=[0,delta),
    B_0=[1/2-delta,1/2+delta),
    B_+=[1-delta,1).

THM-2684 proves that q=7 gives the complete supported THM-2584 rail
envelope.  This companion checks both natural signed scale transporters

    D_+(x)={ p x},             D_-(x)={-p x}.

At positive-length support level, D_+ has identity tooth adjacency and D_-
swaps the endpoint teeth while fixing the central tooth.  For either sign,
the only possible three-tooth words have raw return cylinders obtained by
replacing delta with delta/p^2.  On each cylinder the two intrinsic event
clocks c_q({p x}) and c_q({p^2 x}) agree.  Hence the current stored edge is
constant and no clock-legal three-event rail product exists.

The proof is the endpoint calculation in the module docstring; the exact
q=2..64 sweep is an implementation hostile.  Reflection changes finitely
many half-open endpoint assignments, so all negative-slope claims explicitly
concern positive-length/a.e. support.  This does not declare an arbitrary
signed circle map to be a lawful LRC transporter.  It closes the already
typed reflected candidate rho o D on the THM-2684 envelope.
"""

from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def intersect(left, right):
    """Intersect sorted half-open interval unions, ignoring null contacts."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return merge(out)


def merge(intervals):
    out = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def reflect_support(intervals):
    """Canonical half-open representative of rho(x)={-x}, modulo nulls."""
    out = []
    for left, right in intervals:
        if left == 0:
            out.append((1 - right, Fraction(1)))
        elif right == 1:
            out.append((Fraction(0), 1 - left))
        else:
            out.append((1 - right, 1 - left))
    return merge(out)


def pullback_positive(intervals, slope):
    """Positive-slope pullback x -> {slope*x}."""
    return merge(
        ((branch + left) / slope, (branch + right) / slope)
        for left, right in intervals
        for branch in range(slope)
    )


def pullback_positive_restricted(intervals, slope, windows):
    """Positive pullback branches whose interiors can meet given windows."""
    out = []
    for target_left, target_right in intervals:
        for window_left, window_right in windows:
            first = max(
                0,
                floor_fraction(slope * window_left - target_right) + 1,
            )
            last = min(
                slope - 1,
                ceil_fraction(slope * window_right - target_left) - 1,
            )
            out.extend(
                ((branch + target_left) / slope,
                 (branch + target_right) / slope)
                for branch in range(first, last + 1)
            )
    return merge(out)


def pullback(intervals, slope, sign):
    target = intervals if sign == 1 else reflect_support(intervals)
    return pullback_positive(target, slope)


def measure(intervals):
    return sum((right - left for left, right in intervals), Fraction(0))


def clock_label(x, q):
    x %= 1
    return ((q * x + Fraction(1, 2)).numerator
            // (q * x + Fraction(1, 2)).denominator) % q


def clock_boundaries(q):
    return tuple(Fraction(2 * ell + 1, 2 * q) for ell in range(q))


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def boundary_preimages_inside(boundary, slope, intervals):
    """Only preimages strictly inside the tiny intervals under audit."""
    out = []
    for left, right in intervals:
        first = max(0, floor_fraction(slope * left - boundary) + 1)
        last = min(slope - 1, ceil_fraction(slope * right - boundary) - 1)
        out.extend((branch + boundary) / slope
                   for branch in range(first, last + 1))
    return tuple(out)


def clock_pairs_on(intervals, q, p):
    cuts = {point for interval in intervals for point in interval}
    for boundary in clock_boundaries(q):
        cuts.update(boundary_preimages_inside(boundary, p, intervals))
        cuts.update(boundary_preimages_inside(boundary, p * p, intervals))
    cuts = sorted(cuts)
    pairs = set()
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left + right) / 2
        if not any(a <= midpoint < b for a, b in intervals):
            continue
        pairs.add((
            clock_label(p * midpoint, q),
            clock_label(p * p * midpoint, q),
        ))
    return tuple(sorted(pairs))


def teeth(q):
    delta = Fraction(1, 4 * q)
    return (
        [(Fraction(0), delta)],
        [(Fraction(1, 2) - delta, Fraction(1, 2) + delta)],
        [(Fraction(1) - delta, Fraction(1))],
    )


def signed_atlas(q, sign):
    p = 2 * q - 1
    banks = teeth(q)
    pair_support = {}
    for first in range(3):
        for second in range(3):
            support = intersect(
                banks[first], pullback(banks[second], p, sign)
            )
            pair_support[first, second] = support

    permutation = (0, 1, 2) if sign == 1 else (2, 1, 0)
    require(
        tuple(key for key, value in pair_support.items() if measure(value))
        == tuple((index, permutation[index]) for index in range(3)),
        "signed tooth adjacency changed",
    )

    triples = {}
    pairs = {}
    for first in range(3):
        second = permutation[first]
        third = first
        support = intersect(
            intersect(
                banks[first], pullback(banks[second], p, sign)
            ),
            # Both signed maps square to x -> {p^2*x}.
            pullback_positive_restricted(
                banks[third], p * p, banks[first]
            ),
        )
        triples[first, second, third] = support
        pairs[first] = clock_pairs_on(support, q, p)
    delta2 = Fraction(1, 4 * q * p * p)
    expected = {
        (0, permutation[0], 0): [(Fraction(0), delta2)],
        (1, 1, 1): [
            (Fraction(1, 2) - delta2, Fraction(1, 2) + delta2)
        ],
        (2, permutation[2], 2): [(Fraction(1) - delta2, Fraction(1))],
    }
    require(
        all(measure(triples[key]) == measure(value)
            and triples[key] == value
            for key, value in expected.items()),
        "signed raw triple cylinders changed",
    )
    require(all(all(left == right for left, right in value)
                for value in pairs.values()),
            "a signed triple acquired a nonconstant first clock edge")
    return permutation, triples, pairs


def main():
    controls = 0
    for q in range(2, 65):
        for sign in (1, -1):
            signed_atlas(q, sign)
            controls += 1
    require(controls == 126, "general signed-control universe changed")

    plus = signed_atlas(7, 1)
    minus = signed_atlas(7, -1)
    print("Three-tooth signed p=2q-1 dilation atlas")
    print("scope=positive-length/a.e. support; negative-slope endpoint assignments are null data")
    print("general_exact_controls=q=2..64,sign=+/-; cases=126")
    print("general_plus_tooth_adjacency=identity")
    print("general_minus_tooth_adjacency=endpoint_swap_central_fixed")
    print(f"q7_plus_permutation={plus[0]} raw_triples={plus[1]}")
    print(f"q7_minus_permutation={minus[0]} raw_triples={minus[1]}")
    print(f"q7_plus_clock_pairs={plus[2]}")
    print(f"q7_minus_clock_pairs={minus[2]}")
    print("verdict=PASS: D and typed rho_o_D have no clock-legal three-event product on the three-tooth envelope")
    print("scope_boundary=no arbitrary signed circle map is declared a lawful LRC handoff; other parent carriers and transporters remain open")


if __name__ == "__main__":
    main()
