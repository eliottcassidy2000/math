#!/usr/bin/env python3
"""Exact R=8 integer-entry certificate for AMM 12592.

The script rebuilds the halved sparse junk-flow polytope with exact integer
arithmetic, verifies an integer point, and compares its active constraints
with the denominator-103 vertex from the cross-frontier hostile.  It also
checks the THM-3332 one-row entry conditions on the state entering the first
feed-free row and verifies that the next clamp captures immediately.

No numerical solver and no non-standard package is used.

Reproduction:
  python 04-computation/amm12592_r8_integer_entry_section_kps_s178.py
  python -O 04-computation/amm12592_r8_integer_entry_section_kps_s178.py
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path


R = 8


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fib_pair(n: int) -> tuple[int, int]:
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def compare_five_power_to_phi2(d: int, m: int) -> int:
    """Return sign(5**d - phi**(2*m)) using integer arithmetic."""
    f, f_next = fib_pair(2 * m)
    lucas = 2 * f_next - f
    a = 2 * 5**d - lucas
    if a <= 0:
        return -1
    delta = a * a - 5 * f * f
    return (delta > 0) - (delta < 0)


def floor_gamma_star(m: int) -> int:
    d = 3 * m // 5
    while compare_five_power_to_phi2(d + 1, m) <= 0:
        d += 1
    while compare_five_power_to_phi2(d, m) > 0:
        d -= 1
    return d


def two_g_coeffs(r: int) -> list[int]:
    values = [2]
    b = 1
    for j in range(1, r):
        b = b * (r - j) // j
        values.append((-b if j & 1 else b) - 1)
    return values


def t4_row_load(r: int, d: int) -> list[int]:
    return [
        (-1) ** (d - t) * comb(r - 2 - t, d - t)
        - comb(d + 1, t + 1)
        + 2 * comb(d, t)
        for t in range(d + 1)
    ]


def dot(row: tuple[int, ...], point: tuple[Fraction, ...]) -> Fraction:
    return sum((Fraction(a) * x for a, x in zip(row, point)), Fraction(0))


def rank(rows: list[tuple[int, ...]]) -> int:
    """Exact row rank by rational elimination."""
    if not rows:
        return 0
    a = [[Fraction(x) for x in row] for row in rows]
    pivot_row = 0
    for column in range(len(a[0])):
        pivot = next(
            (i for i in range(pivot_row, len(a)) if a[i][column]), None
        )
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        q = a[pivot_row][column]
        a[pivot_row] = [x / q for x in a[pivot_row]]
        for i in range(len(a)):
            if i != pivot_row and a[i][column]:
                q = a[i][column]
                a[i] = [x - q * y for x, y in zip(a[i], a[pivot_row])]
        pivot_row += 1
        if pivot_row == len(a):
            break
    return pivot_row


def build_polytope():
    """Build y=j/2 sparse constraints with explicit row names."""
    profile = [floor_gamma_star(R + i) for i in range(R)]
    degrees = profile[:-1]
    offsets = [0]
    for d in degrees:
        offsets.append(offsets[-1] + d)

    def variable(i: int, t: int) -> int:
        return offsets[i] + t

    meta: list[tuple[int, int] | None] = [None] * offsets[-1]
    for i, d in enumerate(degrees):
        for t in range(d):
            meta[variable(i, t)] = (i, t)

    inequalities: list[tuple[tuple[int, ...], int, tuple]] = []
    bounds: list[tuple[int | None, int | None]] = [
        (None, None) for _ in meta
    ]

    d0 = degrees[0]
    initial = t4_row_load(R, d0)
    require(all(x % 2 == 0 for x in initial), "initial parity")
    for t in range(d0):
        load = initial[t] // 2
        lower_u = -comb(d0 - 1, t)
        upper_u = comb(d0 - 1, t - 1) if t else 0
        bounds[variable(0, t)] = (load - upper_u, load - lower_u)
    require(0 <= initial[d0] // 2 <= 1, "initial top-cell gate")

    feed = two_g_coeffs(R)
    require(all(x % 2 == 0 for x in feed), "feed parity")
    for i in range(1, len(degrees)):
        d, previous_degree = degrees[i], degrees[i - 1]
        kernel = (1, 1) if d == previous_degree else (1, 2, 1)
        require(d - previous_degree in (0, 1), "degree jump")
        feed0 = feed1 = 0
        if d + i <= R - 1:
            feed0 += feed[d + i] // 2
        if d > previous_degree and d - 1 + i <= R - 1:
            z = feed[d - 1 + i] // 2
            feed0 += z
            feed1 += z
        for t in range(d + 1):
            row = [0] * len(meta)
            for a, coefficient in enumerate(kernel):
                if 0 <= t - a < previous_degree:
                    row[variable(i - 1, t - a)] += coefficient
            if t < d:
                row[variable(i, t)] -= 1
            f = feed0 if t == 0 else (feed1 if t == 1 else 0)
            lower_u = -comb(d - 1, t) if t < d else 0
            upper_u = comb(d - 1, t - 1) if t else 0
            inequalities.append(
                (tuple(row), upper_u - f, ("row-upper", i, t))
            )
            inequalities.append(
                (tuple(-x for x in row), f - lower_u, ("row-lower", i, t))
            )

    previous_degree = degrees[-1]
    last_degree = profile[-1]
    kernel = (
        (1, 1) if last_degree == previous_degree else (1, 2, 1)
    )
    last_row = len(degrees) - 1
    for t in range(last_degree + 1):
        row = [0] * len(meta)
        for a, coefficient in enumerate(kernel):
            if 0 <= t - a < previous_degree:
                row[variable(last_row, t - a)] += coefficient
        inequalities.append(
            (tuple(row), comb(last_degree, t), ("terminal-upper", t))
        )
        inequalities.append(
            (tuple(-x for x in row), 0, ("terminal-lower", t))
        )

    require(all(x is not None for x in meta), "variable metadata")
    return profile, tuple(meta), inequalities, bounds


def active_rows(
    point: tuple[Fraction, ...],
    meta: tuple[tuple[int, int], ...],
    inequalities: list[tuple[tuple[int, ...], int, tuple]],
    bounds: list[tuple[int | None, int | None]],
) -> dict[tuple, tuple[int, ...]]:
    active: dict[tuple, tuple[int, ...]] = {}
    for row, rhs, name in inequalities:
        value = dot(row, point)
        require(value <= rhs, f"inequality violation: {name}")
        if value == rhs:
            active[name] = row
    for j, ((lower, upper), value) in enumerate(zip(bounds, point)):
        unit = tuple(1 if k == j else 0 for k in range(len(point)))
        if lower is not None:
            require(value >= lower, f"lower-bound violation: {meta[j]}")
            if value == lower:
                active[("bound-lower",) + meta[j]] = unit
        if upper is not None:
            require(value <= upper, f"upper-bound violation: {meta[j]}")
            if value == upper:
                active[("bound-upper",) + meta[j]] = unit
    return active


def first_feed_free_row(r: int, profile: list[int]) -> int:
    feed = two_g_coeffs(r)
    for i in range(1, r - 1):
        d, previous_degree = profile[i], profile[i - 1]
        feed0 = d + i <= r - 1 and feed[d + i] != 0
        feed1 = (
            d > previous_degree
            and d - 1 + i <= r - 1
            and feed[d - 1 + i] != 0
        )
        if not feed0 and not feed1:
            return i
    raise RuntimeError("no feed-free row")


def entry_certificate(profile: list[int], witness_rows: tuple[tuple[int, ...], ...]):
    i_pf = first_feed_free_row(R, profile)
    d = profile[i_pf]
    # Variables are y=j/2.  The state entering row i_pf is the preceding
    # junk row, restored here to the original even junk coordinates.
    entry = tuple(2 * x for x in witness_rows[i_pf - 1])
    support = {t: x for t, x in enumerate(entry) if x}
    require(support, "empty entry state")
    m = max(support)
    require(all(x < 0 for x in support.values()), "F1 negativity")
    require(m + 2 < d, "F1 support front")
    a = {t: -x for t, x in support.items()}
    require(a.get(0, 0) <= d - 1, "F2 edge debt")
    margins = {
        t: 2 * comb(d - 1, t)
        - (2 * a.get(t - 1, 0) + a.get(t - 2, 0))
        for t in range(2, m + 3)
    }
    require(all(x >= 0 for x in margins.values()), "F3 marginal surface")

    previous_degree = profile[i_pf - 1]
    kernel = (1, 1) if d == previous_degree else (1, 2, 1)
    load = tuple(
        sum(
            kernel[q] * entry[t - q]
            for q in range(len(kernel))
            if 0 <= t - q < len(entry)
        )
        for t in range(d + 1)
    )
    lower = tuple(-2 * comb(d - 1, t) if t < d else 0 for t in range(d + 1))
    upper = tuple(
        2 * comb(d - 1, t - 1) if t else 0 for t in range(d + 1)
    )
    require(all(lo <= x <= hi for x, lo, hi in zip(load, lower, upper)),
            "entry load is not fully clampable")
    require(all(x == 0 for x in witness_rows[i_pf]), "next junk is not zero")
    capture_clock = i_pf + (a.get(0, 0) + 1) // 2
    require(capture_clock <= R - 2, "F4 capture budget")
    return i_pf, d, entry, load, margins, capture_clock


def main() -> None:
    profile, meta, inequalities, bounds = build_polytope()
    require(profile == [4, 5, 5, 6, 7, 7, 8, 8], "profile")
    require(len(meta) == 42 and len(inequalities) == 106, "polytope size")

    witness_rows = (
        (6, -3, 2, 1),
        (-1, 0, 0, 0, 0),
        (-1, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0),
    )
    witness = tuple(Fraction(x) for row in witness_rows for x in row)
    require(len(witness) == len(meta), "witness length")
    witness_active = active_rows(witness, meta, inequalities, bounds)

    fractional_raw = (
        "7 -7 1 1 0 -490/103 0 0 0 "
        "-1 -169/103 128/103 -2 1 "
        "-1 140/103 -828/103 911/103 -3 1 "
        "-1 -169/103 894/103 -2150/103 751/103 7 0 "
        "0 -375/103 2270/103 -2801/103 146/103 -73/103 1 "
        "1 346/103 799/103 4969/103 419/103 -419/103 8 0"
    )
    fractional = tuple(Fraction(x) for x in fractional_raw.split())
    require(len(fractional) == len(meta), "fractional vertex length")
    fractional_active = active_rows(fractional, meta, inequalities, bounds)

    shared = set(witness_active) & set(fractional_active)
    require(len(witness_active) == 25, "integer active count")
    require(len(fractional_active) == 45, "fractional active count")
    require(len(shared) == 8, "shared active count")
    require(rank(list(fractional_active.values())) == 42, "vertex rank")
    require(rank(list(witness_active.values())) == 22, "integer active rank")
    require(rank([fractional_active[name] for name in shared]) == 7,
            "shared active rank")
    require(
        all(
            dot(fractional_active[name], witness)
            < next(rhs for row, rhs, n in inequalities if n == name)
            for name in set(fractional_active) - shared
            if name[0] not in {"bound-lower", "bound-upper"}
        ),
        "a hostile active inequality did not move inward",
    )

    i_pf, d, entry, clamp_load, margins, capture_clock = entry_certificate(
        profile, witness_rows
    )
    require(i_pf == 3 and d == 6, "entry row")
    require(entry == (-2, 0, 0, 0, 0), "entry state")
    require(clamp_load == (-2, -4, -2, 0, 0, 0, 0), "capture clamp")
    require(margins == {2: 18}, "F3 margin")
    require(capture_clock == 4, "F4 clock")

    profile512 = [floor_gamma_star(512 + i) for i in range(512)]
    i_pf512 = first_feed_free_row(512, profile512)
    require(i_pf512 == 130, "R512 first feed-free row")
    require(profile512[129] == profile512[130] == 383, "R512 entry degrees")
    prefix_variables512 = sum(profile512[:i_pf512])
    full_variables512 = sum(profile512[:-1])
    prefix_causal512 = 2 * sum(
        profile512[i] + 1 for i in range(1, i_pf512)
    )
    require(prefix_variables512 == 44750, "R512 prefix variable count")
    require(full_variables512 == 234117, "R512 full variable count")
    require(prefix_causal512 == 89146, "R512 prefix causal count")

    witness_digest = sha256(
        ";".join(
            f"{i},{t}:{value.numerator}/{value.denominator}"
            for (i, t), value in zip(meta, witness)
        ).encode("ascii")
    ).hexdigest()
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")

    print("AMM12592 R8 INTEGER ENTRY-SECTION EXACT AUDIT")
    print("status=PROVED+FINITE-EXACT")
    print(
        f"polytope=R={R};D0=0;profile={profile};variables={len(meta)};"
        f"inequalities={len(inequalities)}"
    )
    print(
        "integer_witness="
        + ";".join(f"y{i}={list(row)}" for i, row in enumerate(witness_rows))
    )
    print(
        f"integer_active_set=active={len(witness_active)};"
        f"active_rank={rank(list(witness_active.values()))};"
        f"witness_sha256={witness_digest}"
    )
    print(
        f"hostile_separation=fractional_active={len(fractional_active)};"
        f"fractional_active_rank={rank(list(fractional_active.values()))};"
        f"shared_active={len(shared)};shared_rank=7;"
        f"hostile_rows_made_strict={len(fractional_active)-len(shared)};"
        f"new_integer_active_rows={len(witness_active)-len(shared)}"
    )
    print(
        f"entry=i_pf={i_pf};degree={d};junk={list(entry)};support_max=0;"
        f"F1=true;F2=2<={d-1};F3_margin_t2={margins[2]};"
        f"F4_capture_clock={capture_clock}<=R_minus_2={R-2}"
    )
    print(
        f"immediate_capture=next_load={list(clamp_load)};"
        "load_is_inside_clamp_box=true;next_junk_is_zero=true"
    )
    print(
        f"R512_entry_reduction=last_feed_row=129;i_pf={i_pf512};"
        f"entry_degree=383;prefix_variables={prefix_variables512};"
        f"full_variables={full_variables512};"
        f"prefix_causal_inequalities={prefix_causal512}"
    )
    print(
        "R512_half_junk_cone=y_t<=0_for_0<=t<=382;"
        "y_381=y_382=0;-2*y_0<=382;"
        "-2*y_(t-1)-y_(t-2)<=C(382,t)_for_2<=t<=382"
    )
    print(
        "conditional_compiler=an_integer_feasible_R512_prefix_landing_in_"
        "this_cone_extends_by_THM3332_to_an_exact_floor_witness"
    )
    print(
        "conclusion=the_integer_certificate_is_selected_by_a_feed_free_"
        "entry_section_and_not_by_the_denominator_103_vertex_active_face"
    )
    print(
        "nonconsequence=no_integrality_or_total_unimodularity_theorem;"
        "no_R512_integer_point;the_entry_cone_condition_is_sufficient_not_"
        "necessary_for_arbitrary_witnesses;no_new_AMM_deadline_bound"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
