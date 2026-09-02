#!/usr/bin/env python3
"""Exact probe for the minority-anchor owner-renewal trichotomy.

For S={2h} union W, each component of the anchor-safe set is a closed
interval.  The open 1/14-danger teeth of the odd tails either miss a point,
one tooth spans the component, or a shortest covering chain has a proper
crossing.  This script constructs shortest greedy chains with Fraction
arithmetic and checks the endpoint/transition slack numerators.

This is a probe of a deterministic interval lemma, not an LRC(14) census.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, floor, gcd


@dataclass(frozen=True)
class Tooth:
    w: int
    n: int
    left: Fraction
    right: Fraction


def anchor_component(h: int, k: int) -> tuple[Fraction, Fraction]:
    assert 0 <= k < 2 * h
    return Fraction(14 * k + 1, 28 * h), Fraction(14 * k + 13, 28 * h)


def meeting_teeth(w: int, left: Fraction, right: Fraction) -> list[Tooth]:
    # (14n-1)/(14w) < right and left < (14n+1)/(14w).
    lo = floor(w * left - Fraction(1, 14))
    hi = ceil(w * right + Fraction(1, 14))
    ans: list[Tooth] = []
    for n in range(lo - 1, hi + 2):
        a = Fraction(14 * n - 1, 14 * w)
        b = Fraction(14 * n + 1, 14 * w)
        if a < right and left < b:
            ans.append(Tooth(w, n, a, b))
    return ans


def shortest_cover_chain(
    h: int, k: int, speeds: tuple[int, ...]
) -> tuple[list[Tooth] | None, Fraction | None]:
    left, right = anchor_component(h, k)
    teeth = [t for w in speeds for t in meeting_teeth(w, left, right)]

    chain: list[Tooth] = []
    cursor = left
    first = True
    while True:
        eligible = [
            t
            for t in teeth
            if t.left < cursor < t.right
            and (first or not chain or t != chain[-1])
        ]
        if not eligible:
            return None, cursor
        nxt = max(eligible, key=lambda t: (t.right, -t.left, -t.w, -t.n))
        if chain and nxt.right <= chain[-1].right:
            # No interval crossing cursor extends the covered prefix.
            return None, cursor
        chain.append(nxt)
        cursor = nxt.right
        first = False
        if cursor > right:
            return chain, None


def independently_covered(h: int, k: int, speeds: tuple[int, ...]) -> bool:
    """Cell-by-cell open-cover decision, independent of greedy reachability."""
    left, right = anchor_component(h, k)
    teeth = [t for w in speeds for t in meeting_teeth(w, left, right)]
    walls = {left, right}
    for tooth in teeth:
        if left < tooth.left < right:
            walls.add(tooth.left)
        if left < tooth.right < right:
            walls.add(tooth.right)
    ordered = sorted(walls)
    tests = list(ordered)
    tests.extend((a + b) / 2 for a, b in zip(ordered, ordered[1:]))
    return all(any(t.left < x < t.right for t in teeth) for x in tests)


def mass_profile(
    h: int, speeds: tuple[int, ...]
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return safe, load, collision-excess, and anchor-tail overlap masses."""
    union_mass = Fraction(0)
    load_mass = Fraction(0)
    collision_excess = Fraction(0)
    for k in range(2 * h):
        left, right = anchor_component(h, k)
        teeth = [t for w in speeds for t in meeting_teeth(w, left, right)]
        walls = {left, right}
        for tooth in teeth:
            walls.add(max(left, tooth.left))
            walls.add(min(right, tooth.right))
        ordered = sorted(x for x in walls if left <= x <= right)
        for a, b in zip(ordered, ordered[1:]):
            if a == b:
                continue
            midpoint = (a + b) / 2
            multiplicity = sum(t.left < midpoint < t.right for t in teeth)
            width = b - a
            load_mass += multiplicity * width
            if multiplicity:
                union_mass += width
                collision_excess += (multiplicity - 1) * width
    anchor_safe_mass = Fraction(6, 7)
    safe_mass = anchor_safe_mass - union_mass
    anchor_tail_overlap = Fraction(len(speeds), 7) - load_mass
    assert safe_mass == anchor_safe_mass - load_mass + collision_excess
    assert safe_mass == (
        Fraction(6 - len(speeds), 7)
        + anchor_tail_overlap
        + collision_excess
    )
    return safe_mass, load_mass, collision_excess, anchor_tail_overlap


def check_chain(h: int, k: int, chain: list[Tooth]) -> tuple[int, ...]:
    left, right = anchor_component(h, k)
    assert chain
    assert chain[0].left < left < chain[0].right
    assert chain[-1].left < right < chain[-1].right

    left_num = (left - chain[0].left) * (28 * h * chain[0].w)
    right_num = (chain[-1].right - right) * (28 * h * chain[-1].w)
    assert left_num.denominator == right_num.denominator == 1
    assert left_num == (
        chain[0].w * (14 * k + 1) - 2 * h * (14 * chain[0].n - 1)
    )
    assert right_num == (
        2 * h * (14 * chain[-1].n + 1)
        - chain[-1].w * (14 * k + 13)
    )
    assert left_num >= 1 and right_num >= 1
    assert int(left_num) % 2 == int(right_num) % 2 == 1

    qs: list[int] = []
    for one, two in zip(chain, chain[1:]):
        assert one.left < two.left < one.right < two.right
        determinant = two.n * one.w - one.n * two.w
        assert determinant > 0
        q = one.w + two.w - 14 * determinant
        overlap = one.right - two.left
        assert overlap == Fraction(q, 14 * one.w * two.w)
        assert q >= 2 and q % 2 == 0
        pair_gcd = gcd(one.w, two.w)
        assert q % (2 * pair_gcd) == 0
        assert overlap >= Fraction(1, 7 * (one.w * two.w // pair_gcd))
        qs.append(q)

    reconstructed_width = sum(
        (t.right - t.left for t in chain), start=Fraction(0)
    )
    reconstructed_width -= left - chain[0].left
    reconstructed_width -= chain[-1].right - right
    reconstructed_width -= sum(
        (a.right - b.left for a, b in zip(chain, chain[1:])),
        start=Fraction(0),
    )
    assert reconstructed_width == right - left == Fraction(3, 7 * h)

    # Every active tail on the kth physical anchor component has the same
    # MC2 owner bit: 0 on the first h components, 1 on the second h.
    epsilon = int(k >= h)
    for tooth in chain:
        owner_integer = 2 * tooth.n - epsilon * tooth.w
        assert owner_integer % 2 == epsilon
    return tuple(qs)


def profile(h: int, speeds: tuple[int, ...], label: str) -> dict[str, int]:
    assert all(w > 0 and w % 2 == 1 for w in speeds)
    counts = {
        "components": 2 * h,
        "missing": 0,
        "spanning": 0,
        "transition": 0,
        "repeated_owner_chain": 0,
        "max_chain": 0,
        "min_q": 0,
    }
    q_values: list[int] = []
    span_hist: dict[int, int] = {}
    first_missing: tuple[int, Fraction] | None = None
    first_span: tuple[int, int, int] | None = None
    first_transition: tuple[int, tuple[tuple[int, int], ...], tuple[int, ...]] | None = None
    for k in range(2 * h):
        chain, gap = shortest_cover_chain(h, k, speeds)
        if chain is None:
            counts["missing"] += 1
            if first_missing is None:
                assert gap is not None
                first_missing = (k, gap)
            continue
        qs = check_chain(h, k, chain)
        counts["max_chain"] = max(counts["max_chain"], len(chain))
        if len(chain) == 1:
            counts["spanning"] += 1
            # A full anchor component has length 3/(7h), while a tooth has
            # length 1/(7w); strict endpoint coverage forces 3w<h.
            assert 3 * chain[0].w < h
            span_hist[chain[0].w] = span_hist.get(chain[0].w, 0) + 1
            if first_span is None:
                first_span = (k, chain[0].w, chain[0].n)
        else:
            counts["transition"] += 1
            q_values.extend(qs)
            if len({t.w for t in chain}) < len(chain):
                counts["repeated_owner_chain"] += 1
            if first_transition is None:
                first_transition = (
                    k,
                    tuple((t.w, t.n) for t in chain),
                    qs,
                )
    if q_values:
        counts["min_q"] = min(q_values)
    print(label)
    print(" h", h)
    print(" speeds", speeds)
    print(" counts", counts)
    print(" span_hist", dict(sorted(span_hist.items())))
    safe_mass, load_mass, collision_excess, anchor_tail_overlap = mass_profile(
        h, speeds
    )
    print(" safe_mass", safe_mass)
    print(" tail_load_on_anchor_safe", load_mass)
    print(" collision_excess", collision_excess)
    print(" anchor_tail_overlap_sum", anchor_tail_overlap)
    if first_missing is not None:
        missing_k, witness = first_missing
        full_speeds = (2 * h,) + speeds
        distances = tuple(
            min((v * witness) % 1, 1 - (v * witness) % 1)
            for v in full_speeds
        )
        clearance = min(distances)
        binding = tuple(v for v, d in zip(full_speeds, distances) if d == clearance)
        assert clearance >= Fraction(1, 14)
        print(" first_missing", (missing_k, witness, clearance, binding))
    else:
        print(" first_missing", None)
    print(" first_span", first_span)
    print(" first_transition", first_transition)
    return counts


def find_transition_control() -> tuple[int, tuple[int, ...], int, list[Tooth]]:
    # Find a small literal component which is covered, but not by one tooth.
    for h in range(2, 31):
        odds = tuple(range(1, 6 * h, 2))
        for i, u in enumerate(odds):
            for v in odds[i + 1 :]:
                speeds = (u, v)
                for k in range(2 * h):
                    chain, _ = shortest_cover_chain(h, k, speeds)
                    if chain is not None and len(chain) >= 2:
                        check_chain(h, k, chain)
                        return h, speeds, k, chain
    raise AssertionError("transition control not found")


def small_independent_audit() -> int:
    checks = 0
    for h in range(2, 15):
        odds = tuple(range(1, 4 * h, 2))
        banks = [(w,) for w in odds]
        banks.extend((u, v) for i, u in enumerate(odds) for v in odds[i + 1 :])
        for speeds in banks:
            for k in range(2 * h):
                chain, _ = shortest_cover_chain(h, k, speeds)
                assert (chain is not None) == independently_covered(h, k, speeds)
                if chain is not None:
                    check_chain(h, k, chain)
                checks += 1
    return checks


def main() -> None:
    print("small_independent_audit_components", small_independent_audit())
    h0, speeds0, k0, chain0 = find_transition_control()
    print(
        "transition_control",
        h0,
        speeds0,
        k0,
        tuple((t.w, t.n) for t in chain0),
        "q",
        check_chain(h0, k0, chain0),
    )
    profile(10, (1,), "spanning_control")

    common = tuple(11 + 1680 * k for k in range(7)) + (525, 945, 1365, 1575)
    c1287 = profile(420, common + (1287,), "joint_hostile_P1287")
    c9009 = profile(420, common + (9009,), "joint_hostile_P9009")

    # These are safe rows, so at least one component must expose a missing
    # owner bit.  The probe also records how often the other components need
    # a spanning state versus a collision-backed renewal.
    assert c1287["missing"] > 0 and c9009["missing"] > 0
    print("PASS")


if __name__ == "__main__":
    main()
