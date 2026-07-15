#!/usr/bin/env python3
"""Exact global-erosion budget and component-tournament liar for THM-789.

The proved general budget behind this certificate is the following.  Put

    delta = (1/11 - 1/13) / max(U) = 2/(143 max(U)).

The centered interval I_delta is contained in the simultaneous return set
R_U.  Under hypothetical two-sheet tightness,

    E_U + I_delta subset E_U + R_U subset H_(x,y).

If the cyclic gaps of E_U are g_i, one-dimensional thickening is exact:

    mu(E_U + I_delta) = mu(E_U) + sum_i min(g_i, 2 delta).

Kneser-Macbeath also gives

    mu(E_U) + mu(R_U) <= mu(H_(x,y)),

while a 72^10 phase-cell difference packet gives

    mu(R_U) >= max(2 delta, 2 * 72^(-10)).

The computation below exactly replays these quantities on THM-789's core and
then exhibits a limitation of scalar/component Tournament Analysis.  The
pairwise observable is the raw component escape margin; the switch is exact
comparison of those margins, and the tie Hamiltonian path is left-endpoint
order.  Two admissible odd pairs have the same diamond measure, the same full
eroded-diamond measure, the same signed raw component tournament, but
different component-by-component erosion incidence.  Thus the tournament is
telemetry unless its vertices retain signed affine tooth/slope addresses (or
the full eroded margin itself).

The challenged vertex assumption is that runners, raw deep components, or
their scalar ranks suffice.  The faithful carrier is instead a deep component
together with folded-tooth addresses, signed affine slopes, and return-set
incidence.  The quotient preserves raw escape ordering but destroys the
erosion predicate.
"""

from fractions import Fraction as F
from hashlib import sha256


U = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12)
ALPHA = F(1, 13)
BETA = F(1, 11)
GAMMA = BETA - ALPHA
THRESHOLD = F(11, 13)
DELTA = GAMMA / max(U)
EXPECTED_DIGEST = "9c319fb9ea80a9ef193345d85bed9b752a2a7a77cdee570cde8ef4bc00361cd6"


def norm(z: F) -> F:
    z %= 1
    return min(z, 1 - z)


def intersect(left, right, *, allow_points: bool):
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        lo, hi = max(a, c), min(b, d)
        if hi > lo or (allow_points and hi == lo):
            answer.append((lo, hi))
        if b < d:
            i += 1
        elif d < b:
            j += 1
        else:
            i += 1
            j += 1
    return answer


def deep_components(speeds):
    current = [(F(0), F(1))]
    for speed in speeds:
        safe = [
            (
                (F(k) + BETA) / speed,
                (F(k + 1) - BETA) / speed,
            )
            for k in range(speed)
        ]
        current = intersect(current, safe, allow_points=True)
    return tuple(current)


def return_components(speeds):
    """Closure endpoints of R_U in the representative (-1/2,1/2)."""
    current = [(F(-1, 2), F(1, 2))]
    for speed in speeds:
        allowed = sorted(
            (
                (F(k) - GAMMA) / speed,
                (F(k) + GAMMA) / speed,
            )
            for k in range(-speed - 1, speed + 2)
            if (F(k) + GAMMA) / speed > F(-1, 2)
            and (F(k) - GAMMA) / speed < F(1, 2)
        )
        current = intersect(current, allowed, allow_points=False)
    return tuple(current)


def q_value(a: int, b: int, t: F) -> F:
    return norm(a * t) + norm(b * t)


def folded_components(x: int, y: int):
    """Closed components of ||a t||+||b t|| >= 11/13."""
    a, b = (x + y) // 2, (x - y) // 2
    breakpoints = {F(0), F(1)}
    for speed in (a, b):
        breakpoints.update(F(k, 2 * speed) for k in range(2 * speed + 1))
    points = sorted(breakpoints)
    pieces = []
    for left, right in zip(points, points[1:]):
        q_left = q_value(a, b, left)
        q_right = q_value(a, b, right)
        if q_left >= THRESHOLD and q_right >= THRESHOLD:
            pieces.append((left, right))
        elif q_left < THRESHOLD and q_right < THRESHOLD:
            continue
        else:
            root = left + (right - left) * (THRESHOLD - q_left) / (
                q_right - q_left
            )
            pieces.append(
                (root, right) if q_right >= THRESHOLD else (left, root)
            )
    merged = []
    for left, right in pieces:
        if merged and merged[-1][1] == left:
            merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def interval_min_q(a: int, b: int, left: F, right: F) -> F:
    candidates = {left, right}
    for speed in (a, b):
        start = int(2 * speed * left) - 2
        stop = int(2 * speed * right) + 3
        candidates.update(
            point
            for k in range(start, stop + 1)
            if left <= (point := F(k, 2 * speed)) <= right
        )
    return min(q_value(a, b, point) for point in candidates)


def component_profile(components, x: int, y: int, radius: F):
    a, b = (x + y) // 2, (x - y) // 2
    return tuple(
        THRESHOLD - interval_min_q(a, b, left - radius, right + radius)
        for left, right in components
    )


def cyclic_gaps(components):
    return tuple(
        [
            components[i + 1][0] - components[i][1]
            for i in range(len(components) - 1)
        ]
        + [F(1) - components[-1][1] + components[0][0]]
    )


def measure(components) -> F:
    return sum((right - left for left, right in components), F(0))


def eroded_measure(components, radius: F) -> F:
    return sum(
        (max(F(0), right - left - 2 * radius) for left, right in components),
        F(0),
    )


def tournament_order(components, margins):
    """Ascending margins; left endpoints are the fixed tie path."""
    return tuple(
        sorted(
            range(len(components)),
            key=lambda i: (margins[i], components[i][0]),
        )
    )


def format_fraction(value: F) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def format_interval(interval) -> str:
    return f"[{format_fraction(interval[0])},{format_fraction(interval[1])}]"


def format_tuple(values) -> str:
    return "(" + ",".join(format_fraction(value) for value in values) + ")"


def main() -> None:
    deep = deep_components(U)
    returns = return_components(U)
    expected_deep = (
        (F(12, 77), F(7, 44)),
        (F(23, 99), F(21, 88)),
        (F(23, 88), F(32, 121)),
        (F(89, 121), F(65, 88)),
        (F(67, 88), F(76, 99)),
        (F(37, 44), F(65, 77)),
    )
    assert deep == expected_deep
    assert GAMMA == F(2, 143)
    assert DELTA == F(1, 858)
    assert returns == ((-DELTA, DELTA),)

    gaps = cyclic_gaps(deep)
    mu_deep = measure(deep)
    thickened_measure = mu_deep + sum(
        (min(gap, 2 * DELTA) for gap in gaps), F(0)
    )
    assert mu_deep == F(193, 7623)
    assert min(gaps) == F(1, 44) > 2 * DELTA
    assert thickened_measure == F(3895, 99099)
    return_measure_lower = max(2 * DELTA, F(2, 72**10))
    assert return_measure_lower == F(1, 429)

    records = {}
    for pair in ((13, 9), (17, 13)):
        x, y = pair
        folded = folded_components(x, y)
        raw = component_profile(deep, x, y, F(0))
        eroded = component_profile(deep, x, y, DELTA)
        order = tournament_order(deep, raw)
        signs = tuple(value > 0 for value in raw)
        erosion_flags = tuple(value > 0 for value in eroded)
        records[pair] = {
            "folded": folded,
            "raw": raw,
            "eroded": eroded,
            "order": order,
            "signs": signs,
            "erosion_flags": erosion_flags,
            "measure": measure(folded),
            "eroded_measure": eroded_measure(folded, DELTA),
        }

    first = records[(13, 9)]
    second = records[(17, 13)]
    assert first["folded"] == (
        (F(37, 169), F(28, 117)),
        (F(37, 117), F(54, 169)),
        (F(115, 169), F(80, 117)),
        (F(89, 117), F(132, 169)),
    )
    assert second["folded"] == (
        (F(50, 221), F(41, 169)),
        (F(50, 169), F(67, 221)),
        (F(154, 221), F(119, 169)),
        (F(128, 169), F(171, 221)),
    )
    assert first["raw"] == (
        F(159, 572),
        F(-7, 1144),
        F(447, 1573),
        F(447, 1573),
        F(-7, 1144),
        F(159, 572),
    )
    assert first["eroded"] == (
        F(15, 52),
        F(5, 1144),
        F(2825, 9438),
        F(2825, 9438),
        F(5, 1144),
        F(15, 52),
    )
    assert second["raw"] == (
        F(197, 1001),
        F(-59, 1144),
        F(538, 1573),
        F(538, 1573),
        F(-59, 1144),
        F(197, 1001),
    )
    assert second["eroded"] == (
        F(1301, 6006),
        F(-125, 3432),
        F(3415, 9438),
        F(3415, 9438),
        F(-125, 3432),
        F(1301, 6006),
    )
    expected_order = (1, 4, 0, 5, 2, 3)
    expected_signs = (True, False, True, True, False, True)
    assert first["order"] == second["order"] == expected_order
    assert first["signs"] == second["signs"] == expected_signs
    assert first["erosion_flags"] == (True,) * 6
    assert second["erosion_flags"] == expected_signs
    assert first["measure"] == second["measure"] == F(8, 169)
    assert (
        first["eroded_measure"]
        == second["eroded_measure"]
        == F(212, 5577)
    )
    assert thickened_measure < F(8, 169)

    expanded_middle = (deep[1][0] - DELTA, deep[1][1] + DELTA)
    assert expanded_middle == (F(595, 2574), F(823, 3432))
    assert expanded_middle[1] - first["folded"][0][1] == F(5, 10296)
    assert expanded_middle[0] - second["folded"][0][0] == F(215, 43758)
    assert second["folded"][0][1] - expanded_middle[1] == F(125, 44616)

    # A second, pointwise liar: the same deep anchor has identical unsigned
    # odd errors, folded margin, parity opposition, and sharp determinant, but
    # a different signed affine slope makes the natural return interval escape.
    t0 = F(4, 17)
    q_13_9 = q_value(11, 2, t0)
    q_43_13 = q_value(28, 15, t0)
    assert q_13_9 == q_43_13 == F(15, 17)
    assert tuple(sorted((norm(13 * t0), norm(9 * t0)))) == (
        F(1, 17),
        F(2, 17),
    )
    assert tuple(sorted((norm(43 * t0), norm(13 * t0)))) == (
        F(1, 17),
        F(2, 17),
    )
    assert 13 * 2 - 9 * 3 == 43 * 3 - 13 * 10 == -1
    assert F(1, 1000) < DELTA
    assert q_value(28, 15, t0 + F(1, 1000)) - THRESHOLD == F(-1503, 221000)

    canonical = []
    canonical.append("U=" + ",".join(map(str, U)))
    canonical.append("E=" + ";".join(format_interval(row) for row in deep))
    canonical.append("R=" + ";".join(format_interval(row) for row in returns))
    canonical.append("gaps=" + format_tuple(gaps))
    canonical.append(
        "budget="
        + ",".join(
            map(
                format_fraction,
                (mu_deep, thickened_measure, return_measure_lower),
            )
        )
    )
    for pair in ((13, 9), (17, 13)):
        row = records[pair]
        canonical.append(
            f"pair={pair[0]},{pair[1]}"
            + "|H="
            + ";".join(format_interval(component) for component in row["folded"])
            + "|measures="
            + format_tuple((row["measure"], row["eroded_measure"]))
            + "|raw="
            + format_tuple(row["raw"])
            + "|eroded="
            + format_tuple(row["eroded"])
            + "|order="
            + ",".join(map(str, row["order"]))
            + "|signs="
            + "".join("1" if value else "0" for value in row["signs"])
            + "|erosion_flags="
            + "".join("1" if value else "0" for value in row["erosion_flags"])
        )
    canonical.append(
        "anchor_liar="
        + format_tuple(
            (
                t0,
                q_13_9,
                q_43_13,
                q_value(28, 15, t0 + F(1, 1000)) - THRESHOLD,
            )
        )
    )
    digest = sha256(("\n".join(canonical) + "\n").encode()).hexdigest()
    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert digest == EXPECTED_DIGEST

    print("THM-789 global erosion budget + component-tournament liar")
    print(f"U={U}")
    print("E_components=" + " ".join(format_interval(row) for row in deep))
    print(f"R_U=(-{DELTA},{DELTA}) gaps_min={min(gaps)}")
    print(
        f"mu(E)={mu_deep} mu(E+I_delta)={thickened_measure} "
        f"mu(R)_uniform_lower={return_measure_lower}"
    )
    print("general_budget: mu(E)+sum_i min(g_i,2delta) <= mu(H)")
    print("kneser_budget: mu(E)+mu(R) <= mu(H), mu(R)>=max(2delta,2*72^-10)")
    print()
    for pair in ((13, 9), (17, 13)):
        row = records[pair]
        print(f"pair={pair}")
        print("  H_components=" + " ".join(format_interval(c) for c in row["folded"]))
        print(
            f"  mu(H)={row['measure']} mu(H-minus-R)={row['eroded_measure']}"
        )
        print(f"  raw_component_margins={format_tuple(row['raw'])}")
        print(f"  eroded_component_margins={format_tuple(row['eroded'])}")
        print(f"  tournament_order_zero_based={row['order']}")
        print(f"  raw_escape_signs={row['signs']}")
        print(f"  eroded_escape_signs={row['erosion_flags']}")
    print()
    print(f"expanded_C2={format_interval(expanded_middle)}")
    print("pair_13_9_right_overshoot=5/10296")
    print("pair_17_13_containment_slacks=(215/43758,125/44616)")
    print(
        "tournament_fingerprint: score_histogram=(0,1,2,3,4,5) "
        "directed_cycles=0 scc_sizes=(1,1,1,1,1,1) hp_count=1 edge_flips=0"
    )
    print("preserved=raw_component_escape_order_and_sign")
    print("destroyed=signed_affine_tooth_slope_and_return-incidence")
    print()
    print(
        "anchor_slope_liar: t0=4/17 pairs=(13,9),(43,13) "
        "Q=15/17 odd_error_multiset=(1/17,2/17) determinant=-1"
    )
    print(
        "anchor_slope_liar_escape: d=1/1000 in R, "
        "Q_(43,13)(t0+d)-11/13=-1503/221000"
    )
    print(f"canonical_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
