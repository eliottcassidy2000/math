#!/usr/bin/env python3
"""All-level closure of THM-2941's canonical reflected-stalk ``k=1`` family.

For a six-body root ``E``, put ``L=14*lcm(E)``.  This extends the canonical
level-one packet to every reflected level

    Z_(E,q) = {qL-e : e in E},                         q >= 1.

Only one aligned tail remains.  On a body-safe cell ``t=(j+u)/L``, its six
reflected predecessors have union

    U_j^q = union_(e in E) {u : ||q u-e(j+u)/L|| < 1/14}.

By THM-2941(25i), ``mu(U_j^q)<6/7`` implies that the lossless projected
residual has mass greater than ``1/7`` and therefore cannot fit in the one
remaining normalized danger comb.

There are two scale laws.

TOOTHPICK SUBDIVISION.  Write ``u=(r+x)/q`` and ``J=qj+r``.  Then

    (qL-e)(j+u)/L = integer + x-e(J+x)/(qL),

so a level-q coarse cell is the average of q level-one fine cells at ruler
``qL``:

    mu(U_j^q) = (1/q) sum_(0<=r<q) mu(U_(qj+r)^1; qL).  (A)

DIAGONAL DILATION.  Put

    B_j(u)=union_e ({e(j+u)/L}+(-1/14,1/14)).

On a safe cell these are ordinary non-wrapping intervals.  On subdivision
branch r, freeze one such interval as ``I=(a,b)`` at ``u=r/q`` and put
``beta=e/(qL)``.  Its moving diagonal section is exactly

    {x : x in I+beta*x} = (a/(1-beta), b/(1-beta)) intersect [0,1].

Its symmetric difference from I has mass at most ``2beta/(1-beta)``.
The exact finite base-envelope audit below gives

    max_(E,u) mu(B_j(u)) = 265/336,
    max_E sum(E)/L       = 1/6.

Since ``e/L<=1/14``, this proves for every selected cell and every q

    mu(U_j^q)
      <= 265/336 + sum_e 2e/(qL-e)
      <= 265/336 + 14/[3(14q-1)].                       (B)

At q=5 the right side is below 6/7 by 19/23184, and it decreases with q.
Thus (B) closes the infinite tail q>=5.  Exact interval clauses close the
finite head q=1,2,3,4.  As a deliberately redundant hostile control, the
same exact engine also audits every root at all levels 1<=q<=30.

The selected cell is the first safe component's leftmost cell for 2,997
roots and the second component's leftmost cell for the same six exceptions
found at level one.  All arithmetic is integer or Fraction exact.  Every
assertion is a RuntimeError check and remains active under ``python -O``.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as Q
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.out"
)

H = Q(1, 14)
THRESHOLD = Q(6, 7)
SCOUT_Q_MAX = 30
HEAD_Q_MAX = 4
SELECTOR_EXCEPTIONS = (
    (1, 3, 5, 7, 9, 11),
    (1, 3, 5, 7, 9, 12),
    (1, 3, 5, 7, 10, 12),
    (1, 3, 5, 8, 10, 12),
    (1, 3, 6, 8, 10, 12),
    (1, 4, 6, 8, 10, 12),
)
EXPECTED_BASE_MAXIMUM = (
    Q(265, 336),
    Q(1),
    (1, 3, 4, 6, 8, 12),
    336,
    24,
    0,
)
EXPECTED_RATIO_MAXIMUM = (Q(1, 6), (1, 2, 3, 4, 6, 12), 168)
EXPECTED_HEAD_MAXIMA = (
    (1, Q(10028643748, 12527514945)),
    (2, Q(83066098364, 104645906565)),
    (3, Q(23225895503112, 29341517507375)),
    (4, Q(219047327076728, 277109491304595)),
)
EXPECTED_MAXIMUM_BODY = (1, 3, 4, 6, 8, 12)
EXPECTED_SEMANTIC_SHA256 = (
    "c5898923940813859f7ce7401e227cdfb5c1d223b322974ceee76142b10f25fd"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def safe_cell_ranges(E: tuple[int, ...]) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Exact positive-length body-safe ranges on the minimal integer ruler."""
    L = 14 * lcm(*E)
    danger = []
    for e in E:
        radius = L // (14 * e)
        period = L // e
        require(14 * e * radius == L, ("body ruler mismatch", E, e, L))
        for k in range(e + 1):
            center = k * period
            danger.append(
                (max(0, center - radius), min(L, center + radius))
            )
    danger.sort()
    merged = []
    for left, right in danger:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    safe = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < L:
        safe.append((cursor, L))
    require(safe, ("empty six-body carrier", E, L))
    return L, tuple(safe)


def clip(left: Q, right: Q) -> tuple[Q, Q] | None:
    left = max(Q(0), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def direct_multiplier_arcs(L: int, a: int, j: int) -> tuple[tuple[Q, Q], ...]:
    """Pull back ``||a(j+u)/L||<1/14`` directly; ``a`` may exceed L."""
    require(a > 0, ("nonpositive multiplier", L, a, j))
    lo = (a * j) // L - 1
    hi = (a * (j + 1)) // L + 1
    arcs = []
    for n in range(lo, hi + 1):
        left = Q(L * (14 * n - 1), 14 * a) - j
        right = Q(L * (14 * n + 1), 14 * a) - j
        piece = clip(left, right)
        if piece is not None:
            arcs.append(piece)
    return tuple(sorted(set(arcs)))


def reflected_level_arcs(
    L: int, e: int, q: int, j: int
) -> tuple[tuple[Q, Q], ...]:
    """Independent cell law ``(qL-e)u-ej`` modulo the L-ruler."""
    z = q * L - e
    r = (e * j) % L
    lo = (-r) // L - 1
    hi = (z - r) // L + 1
    radius = Q(L, 14 * z)
    arcs = []
    for k in range(lo, hi + 1):
        center = Q(r + k * L, z)
        piece = clip(center - radius, center + radius)
        if piece is not None:
            arcs.append(piece)
    return tuple(sorted(set(arcs)))


def merge_intervals(arcs: tuple[tuple[Q, Q], ...]) -> tuple[tuple[Q, Q], ...]:
    merged = []
    for left, right in sorted(arcs):
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def interval_mass(intervals: tuple[tuple[Q, Q], ...]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def slab_sweep_mass(arcs: tuple[tuple[Q, Q], ...]) -> Q:
    endpoints = sorted({Q(0), Q(1), *(x for arc in arcs for x in arc)})
    total = Q(0)
    for left, right in zip(endpoints, endpoints[1:]):
        middle = (left + right) / 2
        if any(a < middle < b for a, b in arcs):
            total += right - left
    return total


def body_cell_is_safe(L: int, E: tuple[int, ...], j: int) -> bool:
    return all(not direct_multiplier_arcs(L, e, j) for e in E)


def level_union(
    L: int,
    E: tuple[int, ...],
    q: int,
    j: int,
    independent: bool,
) -> tuple[Q, tuple[tuple[Q, Q], ...], int]:
    arcs = []
    comparisons = 0
    for e in E:
        z = q * L - e
        direct = direct_multiplier_arcs(L, z, j)
        if independent:
            reflected = reflected_level_arcs(L, e, q, j)
            require(direct == reflected, ("level clause mismatch", E, L, q, e, j))
            comparisons += 1
        clause_mass = interval_mass(merge_intervals(direct))
        require(
            clause_mass == Q(q * L, 7 * z),
            ("safe-cell individual mass law failed", E, L, q, e, j, clause_mass),
        )
        arcs.extend(direct)
    arc_tuple = tuple(arcs)
    union = merge_intervals(arc_tuple)
    mass = interval_mass(union)
    if independent:
        require(
            mass == slab_sweep_mass(arc_tuple),
            ("head union route mismatch", E, L, q, j, mass),
        )
    return mass, union, comparisons


def fine_subdivision_mass(L: int, E: tuple[int, ...], q: int, j: int) -> Q:
    """Right side of the exact toothpick identity (A)."""
    fine_L = q * L
    total = Q(0)
    for r in range(q):
        J = q * j + r
        arcs = tuple(
            arc
            for e in E
            for arc in direct_multiplier_arcs(fine_L, fine_L - e, J)
        )
        total += interval_mass(merge_intervals(arcs))
    return total / q


def fractional_part(value: Q) -> Q:
    return value - value.numerator // value.denominator


def base_centers(E: tuple[int, ...], L: int, j: int, u: Q) -> tuple[Q, ...]:
    centers = tuple(sorted(fractional_part(Q(e, L) * (j + u)) for e in E))
    require(
        all(H <= center <= 1 - H for center in centers),
        ("base interval wraps on selected safe cell", E, L, j, u, centers),
    )
    return centers


def base_union_gap_mass(centers: tuple[Q, ...]) -> Q:
    """Independent equal-interval union formula on the non-wrapping line."""
    return 2 * H + sum(
        (min(right - left, 2 * H) for left, right in zip(centers, centers[1:])),
        Q(0),
    )


def base_union_merge_mass(centers: tuple[Q, ...]) -> Q:
    arcs = tuple((center - H, center + H) for center in centers)
    return interval_mass(merge_intervals(arcs))


def base_envelope(
    E: tuple[int, ...], L: int, j: int
) -> tuple[Q, Q, int, str]:
    """Exact piecewise-linear maximum of ``mu(B_j(u))`` on 0<=u<=1."""
    midpoint = Q(1, 2)
    floors = {
        e: (Q(e, L) * (j + midpoint)).numerator
        // (Q(e, L) * (j + midpoint)).denominator
        for e in E
    }
    events = {Q(0), Q(1)}
    for e, f in combinations(E, 2):
        for difference in (-2 * H, Q(0), 2 * H):
            u = Q(L, f - e) * (
                difference + floors[f] - floors[e]
            ) - j
            if 0 <= u <= 1:
                events.add(u)
    ordered = sorted(events)
    tests = set(ordered)
    tests.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
    digest = hashlib.sha256()
    values = []
    for u in sorted(tests):
        centers = base_centers(E, L, j, u)
        route_a = base_union_gap_mass(centers)
        route_b = base_union_merge_mass(centers)
        require(route_a == route_b, ("base union route mismatch", E, L, j, u))
        digest.update(f"{u}|{centers}|{route_a}\n".encode())
        values.append((route_a, u))
    maximum = max(values)
    return maximum[0], maximum[1], len(events), digest.hexdigest()


def selected_component(E: tuple[int, ...]) -> int:
    return 1 if E in SELECTOR_EXCEPTIONS else 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3003, ("body universe changed", len(bodies)))
    require(len(SELECTOR_EXCEPTIONS) == 6, "selector exception count changed")

    body_data = []
    base_digest = hashlib.sha256()
    base_rows = []
    ratio_rows = []
    safe_checks = 0
    base_event_count = 0
    for E in bodies:
        L, safe = safe_cell_ranges(E)
        component = selected_component(E)
        require(component < len(safe), ("missing selected safe component", E, safe))
        j = safe[component][0]
        require(body_cell_is_safe(L, E, j), ("selected body cell unsafe", E, L, j))
        safe_checks += len(E)
        base_mass, base_u, event_count, event_digest = base_envelope(E, L, j)
        base_event_count += event_count
        base_row = (base_mass, base_u, E, L, j, component)
        base_rows.append(base_row)
        ratio_rows.append((Q(sum(E), L), E, L))
        base_digest.update(f"{base_row}|{event_count}|{event_digest}\n".encode())
        body_data.append((E, L, j, component))

    base_maximum = max(base_rows)
    ratio_maximum = max(ratio_rows)
    require(
        base_maximum == EXPECTED_BASE_MAXIMUM,
        ("base envelope maximum changed", base_maximum),
    )
    require(
        ratio_maximum == EXPECTED_RATIO_MAXIMUM,
        ("body ratio maximum changed", ratio_maximum),
    )

    tail_correction_q5 = Q(14, 3 * (14 * 5 - 1))
    tail_gap = THRESHOLD - base_maximum[0] - tail_correction_q5
    require(tail_gap == Q(19, 23184), ("analytic tail gap changed", tail_gap))
    require(tail_gap > 0, ("analytic tail no longer closes", tail_gap))

    exact_correction_max = max(
        (
            sum((Q(2 * e, 5 * L - e) for e in E), Q(0)),
            E,
            L,
        )
        for E, L, _j, _component in body_data
    )
    require(
        exact_correction_max[0] <= tail_correction_q5,
        ("coarse tail correction failed", exact_correction_max),
    )

    scout_digest = hashlib.sha256()
    q_records = []
    head_route_comparisons = 0
    head_self_similarity_checks = 0
    exact_clause_mass_checks = 0
    for q in range(1, SCOUT_Q_MAX + 1):
        maximum = None
        for E, L, j, component in body_data:
            independent = q <= HEAD_Q_MAX
            mass, union, comparisons = level_union(L, E, q, j, independent)
            exact_clause_mass_checks += len(E)
            head_route_comparisons += comparisons
            if independent:
                fine_mass = fine_subdivision_mass(L, E, q, j)
                require(
                    mass == fine_mass,
                    ("toothpick subdivision mismatch", E, L, q, j, mass, fine_mass),
                )
                head_self_similarity_checks += q
            require(
                mass < THRESHOLD,
                ("finite reflected-level scout survivor", E, L, q, j, mass),
            )
            if q == 1:
                first_j = safe_cell_ranges(E)[1][0][0]
                first_mass = level_union(L, E, 1, first_j, False)[0]
                require(
                    (first_mass >= THRESHOLD) == (E in SELECTOR_EXCEPTIONS),
                    ("level-one selector classification changed", E, first_mass),
                )
            row = (mass, E, L, j, component)
            if maximum is None or row > maximum:
                maximum = row
            scout_digest.update(f"{q}|{row}|{union}\n".encode())
        require(maximum is not None, ("empty q scout", q))
        require(
            maximum[1] == EXPECTED_MAXIMUM_BODY,
            ("q-level maximum body changed", q, maximum),
        )
        q_records.append((q, *maximum))

    require(
        tuple((q, mass) for q, mass, _E, _L, _j, _component in q_records[:4])
        == EXPECTED_HEAD_MAXIMA,
        ("finite head maxima changed", q_records[:4]),
    )
    finite_head_gap = min(
        THRESHOLD - mass
        for q, mass, _E, _L, _j, _component in q_records
        if q <= HEAD_Q_MAX
    )
    require(
        finite_head_gap == Q(4964583434, 87692604615),
        ("finite head gap changed", finite_head_gap),
    )

    semantic_payload = (
        SELECTOR_EXCEPTIONS,
        tuple(base_rows),
        ratio_maximum,
        exact_correction_max,
        tuple(q_records),
        safe_checks,
        base_event_count,
        head_route_comparisons,
        head_self_similarity_checks,
        exact_clause_mass_checks,
        base_digest.hexdigest(),
        scout_digest.hexdigest(),
        tail_gap,
        finite_head_gap,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest changed", semantic_sha256),
        )

    lines = [
        "LRC14 THM-2941 canonical reflected-stalk all-level k=1 mass closure",
        (
            "scope=all six-body roots E subset {1,...,14};"
            "Z_(E,q)={qL-e:e in E} for every integer q>=1;"
            "one final aligned comb;not the whole j=7 sector"
        ),
        (
            "exact_cell_law=(qL-e)(j+u)/L=q*u-e(j+u)/L mod1;"
            "U_j^q=union_e {u:||q*u-e(j+u)/L||<1/14}"
        ),
        (
            "toothpick_subdivision=u=(r+x)/q,J=qj+r,L'=qL;"
            "mu(U_j^q)=(1/q)sum_(r=0..q-1)mu(U_J^1 at ruler L')"
        ),
        (
            "selector=first safe-component left cell for 2997 roots;"
            "second safe-component left cell for 6 frozen exceptions"
        ),
        f"selector_exceptions={SELECTOR_EXCEPTIONS}",
        f"body_roots={len(bodies)};direct_body_safe_checks={safe_checks}",
        (
            f"base_envelope_max={qtext(base_maximum[0])};u={qtext(base_maximum[1])};"
            f"body={base_maximum[2]};L={base_maximum[3]};j={base_maximum[4]};"
            f"safe_component_index={base_maximum[5]}"
        ),
        (
            f"base_event_rows={base_event_count};base_union_routes="
            "equal_interval_gap_formula+exact_interval_merge;all_equal"
        ),
        (
            f"max_sumE_over_L={qtext(ratio_maximum[0])};"
            f"body={ratio_maximum[1]};L={ratio_maximum[2]}"
        ),
        (
            "diagonal_dilation_law={x:x in (a,b)+(e/(qL))*x}="
            "(a/(1-e/(qL)),b/(1-e/(qL))) intersect [0,1]"
        ),
        (
            "all_q_bound=mu(U_j^q)<=265/336+sum_e 2e/(qL-e)"
            "<=265/336+14/[3(14q-1)]"
        ),
        (
            f"q_ge_5_uniform_gap_6/7_minus_bound={qtext(tail_gap)};"
            f"exact_q5_correction_max={qtext(exact_correction_max[0])};"
            f"correction_max_body={exact_correction_max[1]}"
        ),
        (
            f"finite_head_uniform_gap={qtext(finite_head_gap)};"
            f"head_route_comparisons={head_route_comparisons};"
            f"toothpick_fine_cell_checks={head_self_similarity_checks}"
        ),
        (
            f"redundant_exact_scout=q=1..{SCOUT_Q_MAX};"
            f"body_level_rows={len(bodies)*SCOUT_Q_MAX};"
            f"individual_clause_mass_checks={exact_clause_mass_checks};survivors=0"
        ),
    ]
    for q, mass, E, L, j, component in q_records:
        lines.append(
            f"q={q};maximum_selected_union={qtext(mass)};"
            f"gap={qtext(THRESHOLD-mass)};body={E};L={L};j={j};"
            f"safe_component_index={component}"
        )
    lines.extend(
        (
            (
                "projection_consequence=mu(P_(E,Z_(E,q)))>1/7 for every q>=1;"
                "therefore P cannot lie in one normalized danger comb"
            ),
            (
                "endpoint_status=open danger intervals retained;"
                "strict positive finite-head and analytic-tail gaps make seams irrelevant"
            ),
            f"base_digest_sha256={base_digest.hexdigest()}",
            f"scout_digest_sha256={scout_digest.hexdigest()}",
            f"semantic_sha256={semantic_sha256}",
            (
                "conclusion=canonical reflected-level k=1 family is empty for all "
                "3003 body roots and every integer level q>=1"
            ),
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
