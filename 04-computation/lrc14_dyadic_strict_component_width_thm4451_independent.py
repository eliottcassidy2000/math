"""Independent exact audit of the proposed THM-4451 component-width result.

This file is deliberately self-contained.  It does not import either of the
author's discovery/audit programs.  All geometry is recomputed from the
defining strict inequalities at an exact rational wall decomposition.
"""

from fractions import Fraction as Q
from hashlib import sha256


ZERO = Q(0)
ONE = Q(1)
HALF = Q(1, 2)
RADIUS = Q(1, 14)


def floor_q(x):
    return x.numerator // x.denominator


def mod_one(x):
    return x - floor_q(x)


def is_allowed(n, ternary):
    return n > 0 and n % 2 == 1 and (not ternary or n % 3 != 0)


def in_danger(n, x):
    """Strict membership in D_n at a rational circle point x."""
    y = mod_one(n * x)
    return y < RADIUS or y > ONE - RADIUS


def in_failure(tails, x):
    first = any(in_danger(n, x) for n in tails)
    second = any(in_danger(n, mod_one(x - HALF)) for n in tails)
    return first and second


def strict_walls(tails, quotient=False):
    """All defining walls, plus zero, in physical or doubled coordinate."""
    walls = {ZERO}
    for n in tails:
        for k in range(n):
            for sign in (-1, 1):
                z = mod_one(Q(14 * k + sign, 14 * n))
                walls.add(mod_one(2 * z) if quotient else z)
                z = mod_one(z + HALF)
                walls.add(mod_one(2 * z) if quotient else z)
    return sorted(walls)


def component_widths(tails, quotient=False, essential=False):
    """Exact circular component widths.

    strict mode retains a deleted defining wall.  essential mode fills all
    isolated wall holes and is used only for the P* carrier geometry.
    """
    walls = strict_walls(tails, quotient=quotient)
    count = len(walls)
    lengths = []
    live_cell = []
    live_wall = []

    def pred(y):
        x = mod_one(y / 2) if quotient else y
        return in_failure(tails, x)

    for i, left in enumerate(walls):
        right = walls[i + 1] if i + 1 < count else ONE
        lengths.append(right - left)
        live_cell.append(pred((left + right) / 2))
        live_wall.append(pred(left))

    live_indices = [i for i, yes in enumerate(live_cell) if yes]
    if not live_indices:
        return []
    if all(live_cell):
        return [ONE]

    # Break the circle at an inactive cell.  A union-find implementation makes
    # the strict/essential distinction explicit at every intervening wall.
    parent = list(range(count))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    for wall_i in range(count):
        left_i = (wall_i - 1) % count
        right_i = wall_i
        if not (live_cell[left_i] and live_cell[right_i]):
            continue
        if essential or live_wall[wall_i]:
            union(left_i, right_i)

    totals = {}
    for i in live_indices:
        root = find(i)
        totals[root] = totals.get(root, ZERO) + lengths[i]
    return sorted(totals.values(), reverse=True)


def carrier_geometry(a, b):
    """Largest essential P* component and least positive circular gap."""
    tails = (a, b)
    walls = strict_walls(tails)
    count = len(walls)
    lengths = []
    active = []
    for i, left in enumerate(walls):
        right = walls[i + 1] if i + 1 < count else ONE
        lengths.append(right - left)
        active.append(in_failure(tails, (left + right) / 2))

    if not any(active):
        return ZERO, ONE
    if all(active):
        return ONE, ZERO

    def cyclic_run_sums(value):
        # Rotate after a cell of the opposite value, so no desired run wraps.
        pivot = next(i for i, bit in enumerate(active) if bit != value)
        ordered = [(pivot + 1 + j) % count for j in range(count)]
        sums = []
        run = ZERO
        for i in ordered:
            if active[i] == value:
                run += lengths[i]
            elif run:
                sums.append(run)
                run = ZERO
        if run:
            sums.append(run)
        return sums

    widths = cyclic_run_sums(True)
    gaps = cyclic_run_sums(False)
    return max(widths), min(gaps)


def cap(n, length):
    x = 2 * n * length
    m = floor_q(x)
    r = x - m
    return (Q(2 * m, 7) + min(r, Q(2, 7))) / (2 * n)


def surplus(n, length):
    return cap(n, length) - Q(2, 7) * length


def enumerate_reduction(length, ternary, safe_thresholds):
    """Rebuild the inclusive capacity boxes plus finite star prefixes."""
    delta = Q(8, 7) * length
    a_bound = Q(15, 49) / delta
    finite = set()
    unbounded = []
    bounded_pairs = []

    for a in range(1, floor_q(a_bound) + 1):
        if not is_allowed(a, ternary):
            continue
        rem_after_a = delta - surplus(a, length)
        assert rem_after_a > 0
        b_bound = Q(10, 49) / rem_after_a
        for b in range(a + 1, floor_q(b_bound) + 1):
            if not is_allowed(b, ternary):
                continue
            rem_after_b = rem_after_a - surplus(b, length)
            if rem_after_b <= 0:
                unbounded.append((a, b))
                continue
            c_bound = Q(5, 49) / rem_after_b
            bounded_pairs.append((a, b, c_bound))
            for c in range(b + 1, floor_q(c_bound) + 1):
                if is_allowed(c, ternary):
                    finite.add((a, b, c))

    assert set(unbounded) == set(safe_thresholds)
    star_prefix = set()
    for (a, b), threshold in safe_thresholds.items():
        w, g = carrier_geometry(a, b)
        for c in range(b + 1, threshold):
            if is_allowed(c, ternary):
                row = (a, b, c)
                finite.add(row)
                star_prefix.add(row)

        # The threshold itself is the first claimed tail in the infinite
        # region.  The inequality then improves monotonically with c.
        assert is_allowed(threshold, ternary)
        ell = Q(1, 7 * threshold)
        assert ell < g
        # Empty carriers have no star centre: the only components are teeth.
        star_bound = ell if w == 0 else max(ell, w + 2 * ell)
        assert star_bound < length

    return {
        "a_bound": a_bound,
        "unbounded": sorted(unbounded),
        "bounded_pairs": bounded_pairs,
        "rows": sorted(finite),
        "prefix_count": len(star_prefix),
    }


def audit_domain(name, length, ternary, thresholds, expected_rows,
                 expected_winner, expected_runner_up):
    reduction = enumerate_reduction(length, ternary, thresholds)
    rows = reduction["rows"]
    assert len(rows) == expected_rows

    values = []
    for tails in rows:
        widths = component_widths(tails)
        assert widths
        physical = widths[0]
        quotient_widths = component_widths(tails, quotient=True)
        quotient = quotient_widths[0]
        assert quotient == 2 * physical
        values.append((physical, tails, widths, quotient_widths))

    values.sort(key=lambda row: (row[0], row[1]), reverse=True)
    maximum = values[0][0]
    leaders = sorted(t for value, t, _, _ in values if value == maximum)
    distinct_values = sorted({row[0] for row in values}, reverse=True)
    runner_up = distinct_values[1]
    assert (maximum, leaders) == expected_winner
    assert runner_up == expected_runner_up
    assert maximum == length

    winning_widths = next(w for v, t, w, _ in values if t == leaders[0])
    quotient_winning_widths = next(q for v, t, _, q in values if t == leaders[0])

    print(name)
    print("  target", length, "a_bound", reduction["a_bound"])
    print("  unbounded_pairs", reduction["unbounded"])
    print("  finite_rows", len(rows), "star_prefix_rows", reduction["prefix_count"])
    print("  maximum", maximum, "leaders", leaders)
    print("  runner_up", runner_up)
    print("  winning_physical_widths", winning_widths)
    print("  winning_quotient_widths", quotient_winning_widths)
    return reduction, values


def main():
    all_thresholds = {
        (1, 3): 7,
        (1, 5): 7,
        (1, 7): 21,
        (3, 5): 15,
        (3, 7): 21,
        (5, 7): 21,
    }
    ternary_thresholds = {
        (1, 5): 11,
        (1, 7): 35,
        (5, 7): 35,
    }

    all_reduction, all_values = audit_domain(
        "all odd",
        Q(17, 693),
        False,
        all_thresholds,
        123,
        (Q(17, 693), [(1, 9, 11)]),
        Q(1, 42),
    )
    ternary_reduction, ternary_values = audit_domain(
        "odd 3-unit",
        Q(19, 1001),
        True,
        ternary_thresholds,
        209,
        (Q(19, 1001), [(1, 11, 13)]),
        Q(29, 1547),
    )

    # Exact carrier table from the proof, recomputed from definitions.
    expected_carriers = {
        (1, 3): (ZERO, ONE),
        (1, 5): (ZERO, ONE),
        (1, 7): (Q(1, 98), Q(6, 49)),
        (3, 5): (Q(1, 210), Q(5, 42)),
        (3, 7): (Q(1, 98), Q(19, 98)),
        (5, 7): (Q(1, 98), Q(1, 14)),
    }
    for pair, expected in expected_carriers.items():
        actual = carrier_geometry(*pair)
        assert actual == expected, (pair, actual, expected)
    print("carrier_geometry", expected_carriers)

    # The endpoint-hole control must split the apparent 1/49 interval.
    control = component_widths((1, 7, 13))
    assert control.count(Q(1, 91)) == 4
    assert control.count(Q(1, 98)) == 8
    assert len(control) == 12
    assert not in_failure((1, 7, 13), Q(1, 14))
    print("endpoint_control_(1,7,13)", control)

    # Winning multiplicities and the report's complete width lists.
    all_w = component_widths((1, 9, 11))
    assert all_w.count(Q(17, 693)) == 4
    assert all_w.count(Q(13, 1386)) == 4
    assert len(all_w) == 8
    ter_w = component_widths((1, 11, 13))
    assert ter_w.count(Q(19, 1001)) == 4
    assert ter_w.count(Q(17, 2002)) == 4
    assert ter_w.count(Q(3, 2002)) == 4
    assert len(ter_w) == 12

    # Capacity algebra and the claimed integer candidate sets.
    for n in range(1, 400, 2):
        for length in (Q(17, 693), Q(19, 1001), Q(1, 100), Q(2, 17)):
            x = 2 * n * length
            r = x - floor_q(x)
            expected = (min(r, Q(2, 7)) - Q(2, 7) * r) / (2 * n)
            assert surplus(n, length) == expected
            assert ZERO <= expected <= Q(5, 49 * n)

    # Report-text correction: (7,11) has a cutoff above b, but no admissible
    # odd c lies below it.  It is therefore harmless yet should not be called
    # a pair whose cutoff is at most b.
    extra = [row for row in all_reduction["bounded_pairs"] if row[:2] == (7, 11)]
    assert extra == [(7, 11, Q(165, 14))]
    assert not any(t[:2] == (7, 11) for t in all_reduction["rows"])
    print("wording_correction_(7,11)_cutoff", Q(165, 14))

    # Digest the exact evaluated row/value table, making this run auditable
    # without inflating stdout with 332 rows.
    payload = "\n".join(
        f"{tag}:{a},{b},{c}:{value}"
        for tag, values in (("A", all_values), ("T", ternary_values))
        for value, (a, b, c), _, _ in sorted(values, key=lambda z: z[1])
    ).encode()
    print("row_value_sha256", sha256(payload).hexdigest())
    print("PASS")


if __name__ == "__main__":
    main()
