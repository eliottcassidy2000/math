"""Clean-room exact audit for the q=4 one-v2 component proposal.

No author code is imported.  Components are reconstructed from the original
strict danger predicate on an exact rational wall decomposition.
"""

from fractions import Fraction as Q
from hashlib import sha256


ZERO = Q(0)
ONE = Q(1)
QUARTER = Q(1, 4)
HALF = Q(1, 2)
RADIUS = Q(1, 14)


def need(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def floor_q(x):
    return x.numerator // x.denominator


def mod_one(x):
    return x - floor_q(x)


def dangerous(n, x):
    y = mod_one(n * x)
    return y < RADIUS or y > ONE - RADIUS


def q4_failure(r, a, b, x):
    tails = (2 * r, a, b)
    return all(
        any(dangerous(n, mod_one(x + j * QUARTER)) for n in tails)
        for j in range(4)
    )


def pair_failure(a, b, x):
    return (
        (dangerous(a, x) or dangerous(b, x))
        and (
            dangerous(a, mod_one(x + HALF))
            or dangerous(b, mod_one(x + HALF))
        )
    )


def decomposed_failure(r, a, b, x):
    first = dangerous(2 * r, x) and pair_failure(
        a, b, mod_one(x + QUARTER)
    )
    second = dangerous(2 * r, mod_one(x + QUARTER)) and pair_failure(a, b, x)
    return first or second


def walls_for_shifted_dangers(speeds, shifts, multiplier=1):
    walls = {ZERO}
    for n in speeds:
        for shift in shifts:
            for k in range(n):
                for sign in (-1, 1):
                    z = mod_one(Q(14 * k + sign, 14 * n) - shift)
                    walls.add(mod_one(multiplier * z))
    return sorted(walls)


def component_widths(r, a, b, quotient=False):
    multiplier = 4 if quotient else 1
    walls = walls_for_shifted_dangers(
        (2 * r, a, b),
        (ZERO, QUARTER, HALF, 3 * QUARTER),
        multiplier,
    )
    size = len(walls)

    def predicate(y):
        x = mod_one(y / 4) if quotient else y
        return q4_failure(r, a, b, x)

    cell_lengths = []
    live_cells = []
    live_walls = []
    for i, left in enumerate(walls):
        right = walls[i + 1] if i + 1 < size else ONE
        cell_lengths.append(right - left)
        live_cells.append(predicate((left + right) / 2))
        live_walls.append(predicate(left))

    live = [i for i, value in enumerate(live_cells) if value]
    if not live:
        return []
    if len(live) == size:
        return [ONE]

    parent = list(range(size))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    for wall_i in range(size):
        left_i = (wall_i - 1) % size
        right_i = wall_i
        if live_cells[left_i] and live_cells[right_i] and live_walls[wall_i]:
            union(left_i, right_i)

    totals = {}
    for i in live:
        root = find(i)
        totals[root] = totals.get(root, ZERO) + cell_lengths[i]
    return sorted(totals.values(), reverse=True)


def identity_points(r, a, b):
    walls = walls_for_shifted_dangers(
        (2 * r, a, b), (ZERO, QUARTER, HALF, 3 * QUARTER)
    )
    points = list(walls)
    for i, left in enumerate(walls):
        right = walls[i + 1] if i + 1 < len(walls) else ONE
        points.append((left + right) / 2)
    return points


def allowed(n, ternary):
    return n % 2 == 1 and (not ternary or n % 3 != 0)


def audit_box(name, target, r_values, tail_values, ternary,
              predicted_equalities):
    rows = []
    for r in r_values:
        need(allowed(r, ternary), ("disallowed r", r, ternary))
        for i, a in enumerate(tail_values):
            if not allowed(a, ternary):
                continue
            for b in tail_values[i + 1 :]:
                if not allowed(b, ternary):
                    continue

                # The decomposition is checked pointwise on every defining
                # wall and every complementary open cell.
                for x in identity_points(r, a, b):
                    need(
                        q4_failure(r, a, b, x)
                        == decomposed_failure(r, a, b, x),
                        ("decomposition", r, a, b, x),
                    )

                physical_widths = component_widths(r, a, b)
                quotient_widths = component_widths(r, a, b, quotient=True)
                physical = physical_widths[0] if physical_widths else ZERO
                quotient = quotient_widths[0] if quotient_widths else ZERO
                need(quotient == 4 * physical, ("quotient", r, a, b))

                analytic_bound = min(Q(1, 14 * r), Q(1, 7 * max(a, b)))
                need(physical <= analytic_bound, ("bound", r, a, b, physical))
                rows.append(
                    (physical, r, a, b, physical_widths, quotient_widths)
                )

    maximum = max(row[0] for row in rows)
    equality = sorted(
        (r, (a, b)) for value, r, a, b, _, _ in rows if value == target
    )
    above = [(value, r, a, b) for value, r, a, b, _, _ in rows if value > target]
    need(not above, ("above target", name, above))
    need(maximum == target, ("maximum", name, maximum, target))
    need(
        equality == sorted(predicted_equalities),
        ("equality", name, equality, predicted_equalities),
    )

    distinct = sorted({row[0] for row in rows}, reverse=True)
    runner_up = distinct[1]
    print(name)
    print("  rows", len(rows), "target", target, "maximum", maximum)
    print("  equality", equality)
    print("  runner_up", runner_up)
    for equality_row in equality:
        er, (ea, eb) = equality_row
        row = next(row for row in rows if row[1:4] == (er, ea, eb))
        print("  equality_widths", equality_row, row[4], row[5])
    return rows


def main():
    all_rows = audit_box(
        "all odd",
        Q(1, 98),
        (1, 3, 5, 7),
        (1, 3, 5, 7, 9, 11, 13),
        False,
        ((1, (7, 9)), (3, (7, 13)), (5, (3, 7))),
    )
    ternary_rows = audit_box(
        "odd 3-unit",
        Q(1, 110),
        (1, 5, 7),
        (1, 5, 7, 11, 13),
        True,
        ((5, (1, 11)),),
    )

    # The finite lists are exactly the inclusive residuals forced by
    # min(1/(14r),1/(7M)) >= target.
    need(len(all_rows) == 4 * (7 * 6 // 2) == 84, "all-row count")
    need(len(ternary_rows) == 3 * (5 * 4 // 2) == 30, "3-unit-row count")

    payload = "\n".join(
        f"{tag}:{r}:{a},{b}:{value}"
        for tag, rows in (("A", all_rows), ("T", ternary_rows))
        for value, r, a, b, _, _ in sorted(rows, key=lambda z: z[1:4])
    ).encode()
    print("row_value_sha256", sha256(payload).hexdigest())
    print("PASS")


if __name__ == "__main__":
    main()
