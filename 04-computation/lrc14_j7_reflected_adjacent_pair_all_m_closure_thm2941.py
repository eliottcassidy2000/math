#!/usr/bin/env python3
r"""Universal adjacent-pair closure and the peeled reflected D=1 tail.

Let ``E`` be any six-element subset of ``{1,...,14}``, put

    L=14*lcm(E),

and give its reflected labels arbitrary positive integer levels ``q_e``.
Write ``m=min q_e``.  This referee proves:

> If two consecutive labels ``a,a+1`` have a common level ``Q``, then the
> packet closes on one explicit body-safe cell.  The other four levels are
> arbitrary and ``Q`` need not equal ``m``.

The cells form a two-chart atlas.  If ``r=max(E)<=12``, use ``j=L/14``.  If
``r>=13``, let ``s`` be the bottom of the terminal consecutive block ending
at ``r`` and use

    j=15L/(14s).                                             (1)

Because ``s in E``, (1) is integral.  Six labels force ``s>=8``.  The exact
3003-body census below checks that these whole open cells are body-safe; the
terminal-block definition is the coordinate that makes the near-danger
labels enlarge ``L`` exactly when extra endpoint room is needed.

Put

    h=floor(aj/L),
    C=1/7-j/L-1/(2L),
    B=h+(3-a)/7.                                             (2)

Lifting the wrapped endpoint pair so its centres differ by ``j`` gives the
repeated word

    aL_k < (a+1)L_k < aR_k < (a+1)R_k,

with successive teeth disjoint.  Summing the ``Q`` overlaps gives

    P(E,a,Q)
      = LQ [ Q(L/7-j-1/2) + h+(3-a)/7 ]
        /[(LQ-a)(LQ-a-1)]
      = LQ(QLC+B)/[(LQ-a)(LQ-a-1)].                         (3)

This formula includes phase wrap through ``h``.  If ``X=LQ`` and
``D=(X-a)(X-a-1)<X^2``, then ``XC+B>0`` and

    P > C + min(0,B)/(LQ)
      >= C + min(0,B)/(Lm).                                 (4)

The singleton excess above ``6/7`` is

    epsilon=sum_e e/[7(q_eL-e)]
            <= sum(E)/[7m(L-r)].                            (5)

For every one of the 6435 labelled body/consecutive-pair rows, exact rational
arithmetic verifies the compatible margin

    C + min(0,B)/L - sum(E)/[7(L-r)] > 0.                   (6)

If the coefficient of ``1/m`` in (4)--(5) is negative, (6) is its worst
case ``m=1``; if nonnegative, positivity is immediate from ``C>0``.  Thus
one-edge Bonferroni proves the theorem for every integer level word in the
stated scope.  Direct multiplier pullback, the independent reflected-cell
formula, and (3) agree for every body/pair at ``Q=1,2,3,7``.

For the global ``D=1`` branch, a nontrivial mask gives its vertices level
``m`` and the others level ``m+1``.  The theorem removes exactly the masks
with a monochromatic consecutive pair.  Rebuilding all ``3003*62=186186``
exact joint spanning-tree rows removes 133200 and leaves 52986.  The largest
remaining tree threshold is below 241, hence the combined result closes
every reflected ``D=1`` packet for ``m>=241``; certificate failure forces
``m<=240``.

As a redundant hostile control, the former unique worst ray

    E=(1,2,3,4,6,12), q=(m,m,m,m,m,m+1)

is swept exactly on ``j=12`` for ``1<=m<=296``.  The theorem itself is
all-level and does not depend on that finite sweep.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve()
BASE_RELATIVE = Path("04-computation") / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ROOT = next(parent for parent in HERE.parents if (parent / BASE_RELATIVE).is_file())
BASE_SOURCE = ROOT / BASE_RELATIVE
TREE_SOURCE = ROOT / "04-computation" / "lrc14_j7_reflected_hostile_carrier_tree_allq_bridge_thm2941.py"
OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_adjacent_pair_all_m_closure_thm2941.out"

EXPECTED_BASE_SOURCE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_TREE_SOURCE_SHA256 = "39ba7ba13197669d654e66cb2043dda619e9644b2f4b2667a11e420a65aa0876"
EXPECTED_SEMANTIC_SHA256 = "7478ca70d6178fda71f1cf463bbe7ac88303ff0494a62b3dfbdc181f14953a04"

PAIR_AUDIT_LEVELS = (1, 2, 3, 7)
EXPECTED_CHART_COUNTS = ((0, 924), (8, 1), (9, 8), (10, 36), (11, 120), (12, 330), (13, 792), (14, 792))
EXPECTED_PAIR_MARGIN_MINIMUM = (
    F(453354931, 50853085920),
    (8, 9, 10, 11, 12, 13),
    720720,
    96525,
    8,
    12,
)
EXPECTED_D1_COVERAGE = 133200
EXPECTED_ORIGINAL_D1_MAXIMUM = (
    F(551653043364, 1857935729), (1, 2, 3, 4, 6, 12), 31
)
EXPECTED_PEELED_D1_MAXIMUM = (
    F(122608946603, 509473633), (1, 2, 3, 4, 6, 12), 5
)

HOSTILE_BODY = (1, 2, 3, 4, 6, 12)
HOSTILE_J = 12
HOSTILE_FINITE_MAX_M = 296
EXPECTED_HOSTILE_SAFE_RANGES = (
    (12, 13), (15, 26), (30, 39), (45, 52), (60, 69), (71, 78),
    (90, 97), (99, 108), (116, 123), (129, 138), (142, 153), (155, 156),
)
EXPECTED_HOSTILE_WORST = (F(7636804, 12540055), 2, 6)
EXPECTED_HOSTILE_PAIR_MARGIN_MINIMUM = (F(1845934, 35374185), 1)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "reflected interval engine changed")
require(sha256(TREE_SOURCE) == EXPECTED_TREE_SOURCE_SHA256, "joint-tree source changed")
R = load("adjacent_pair_interval_engine", BASE_SOURCE)
TREE = load("adjacent_pair_joint_tree_engine", TREE_SOURCE)


def body_chart(body: tuple[int, ...]):
    ruler = 14 * lcm(*body)
    maximum = body[-1]
    if maximum <= 12:
        terminal_start = 0
        j = ruler // 14
    else:
        terminal_start = maximum
        while terminal_start - 1 in body:
            terminal_start -= 1
        require(terminal_start >= 8, ("terminal block too long", body, terminal_start))
        require((ruler // 14) % terminal_start == 0,
                ("terminal chart not integral", body, ruler, terminal_start))
        j = 15 * (ruler // 14) // terminal_start
    require(0 < j < ruler // 7, ("chart outside overlap range", body, ruler, j))
    return ruler, maximum, terminal_start, j


def pair_parameters(body: tuple[int, ...], first_label: int):
    ruler, maximum, terminal_start, j = body_chart(body)
    require(first_label in body and first_label + 1 in body,
            ("pair absent from body", body, first_label))
    h = (first_label * j) // ruler
    c = F(1, 7) - F(j, ruler) - F(1, 2 * ruler)
    b = F(h) + F(3 - first_label, 7)
    require(c > 0, ("nonpositive chart floor", body, ruler, j, c))
    return ruler, maximum, terminal_start, j, h, c, b


def pair_overlap_formula(body: tuple[int, ...], first_label: int, level: int) -> F:
    require(level >= 1, ("nonpositive pair level", level))
    ruler, _maximum, _terminal_start, j, h, _c, _b = pair_parameters(body, first_label)
    numerator = F(ruler * level) * (
        F(level) * (F(ruler, 7) - j - F(1, 2))
        + h
        + F(3 - first_label, 7)
    )
    denominator = (ruler * level - first_label) * (ruler * level - first_label - 1)
    return numerator / denominator


def compatible_margin_at_one(body: tuple[int, ...], first_label: int):
    ruler, maximum, terminal_start, j, h, c, b = pair_parameters(body, first_label)
    floor = c + min(F(0), b) / ruler
    epsilon_bound = F(sum(body), 7 * (ruler - maximum))
    return floor - epsilon_bound, floor, epsilon_bound, ruler, j, terminal_start, h, b, c


def interval_intersection_mass(first, second) -> F:
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            total += right - left
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def endpoint_event_mass(arcs) -> F:
    events: dict[F, int] = {}
    for left, right in arcs:
        events[left] = events.get(left, 0) + 1
        events[right] = events.get(right, 0) - 1
    active = 0
    total = F(0)
    previous = None
    for endpoint in sorted(events):
        if previous is not None and active > 0:
            total += endpoint - previous
        active += events[endpoint]
        require(active >= 0, ("negative endpoint multiplicity", endpoint, active))
        previous = endpoint
    require(active == 0, ("endpoint sweep did not close", active))
    return total


def monochromatic_consecutive_witness(body: tuple[int, ...], mask: int):
    positions = {label: index for index, label in enumerate(body)}
    for first_label in range(1, 14):
        if first_label in positions and first_label + 1 in positions:
            first = positions[first_label]
            second = positions[first_label + 1]
            low = bool(mask & (1 << first))
            if low == bool(mask & (1 << second)):
                return first_label, low
    return None


def hostile_epsilon(m: int) -> F:
    levels = (m, m, m, m, m, m + 1)
    return sum((F(e, 7 * (q * 168 - e)) for e, q in zip(HOSTILE_BODY, levels)), F(0))


def hostile_pair_formula(m: int) -> F:
    return F(12 * m * (161 * m + 4), (168 * m - 1) * (168 * m - 2))


def hostile_finite_profile(m: int):
    levels = (m, m, m, m, m, m + 1)
    clauses = []
    singleton_sum = F(0)
    checks = 0
    for e, q in zip(HOSTILE_BODY, levels):
        direct = R.direct_multiplier_arcs(168, q * 168 - e, HOSTILE_J)
        reflected = R.reflected_level_arcs(168, e, q, HOSTILE_J)
        require(direct == reflected, ("hostile direct/reflected mismatch", m, e))
        clause_mass = R.interval_mass(R.merge_intervals(direct))
        require(clause_mass == F(q * 168, 7 * (q * 168 - e)),
                ("hostile singleton law failed", m, e))
        singleton_sum += clause_mass
        clauses.append(direct)
        checks += 1
    arcs = tuple(arc for clause in clauses for arc in clause)
    merged = R.merge_intervals(arcs)
    mass = R.interval_mass(merged)
    require(mass == endpoint_event_mass(arcs), ("hostile merge/event mismatch", m))
    debt = hostile_epsilon(m)
    pair = interval_intersection_mass(clauses[0], clauses[1])
    require(singleton_sum == F(6, 7) + debt, ("hostile singleton identity failed", m))
    require(pair == hostile_pair_formula(m), ("hostile pair formula failed", m, pair))
    require(mass <= singleton_sum - pair and pair > debt and mass < F(6, 7),
            ("hostile finite closure failed", m, mass, pair, debt))
    return m, mass, len(merged), pair, debt, pair - debt, checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3003, ("body universe changed", len(bodies)))

    chart_counter: Counter[int] = Counter()
    safe_checks = 0
    pair_rows = []
    pair_census_digest = hashlib.sha256()
    direct_formula_checks = 0
    for body in bodies:
        ruler, maximum, terminal_start, j = body_chart(body)
        require(R.body_cell_is_safe(ruler, body, j), ("chart cell unsafe", body, ruler, j))
        chart_counter[terminal_start] += 1
        safe_checks += 1
        for first_label in range(1, 14):
            if first_label not in body or first_label + 1 not in body:
                continue
            margin_row = compatible_margin_at_one(body, first_label)
            margin, floor, epsilon_bound, *_rest = margin_row
            require(margin > 0, ("compatible pair margin failed", body, first_label, margin_row))
            pair_record = (margin, body, ruler, j, terminal_start, first_label, floor, epsilon_bound, margin_row[-2], margin_row[-1])
            pair_rows.append(pair_record)
            for level in PAIR_AUDIT_LEVELS:
                direct_first = R.direct_multiplier_arcs(ruler, level * ruler - first_label, j)
                direct_second = R.direct_multiplier_arcs(ruler, level * ruler - first_label - 1, j)
                reflected_first = R.reflected_level_arcs(ruler, first_label, level, j)
                reflected_second = R.reflected_level_arcs(ruler, first_label + 1, level, j)
                require(direct_first == reflected_first and direct_second == reflected_second,
                        ("direct/reflected pair mismatch", body, first_label, level))
                overlap = interval_intersection_mass(direct_first, direct_second)
                formula = pair_overlap_formula(body, first_label, level)
                require(overlap == formula, ("pair formula mismatch", body, first_label, level, overlap, formula))
                row = (body, ruler, j, terminal_start, first_label, level, overlap)
                pair_census_digest.update(f"{row}\n".encode())
                direct_formula_checks += 1
    chart_counts = tuple(sorted(chart_counter.items()))
    require(chart_counts == EXPECTED_CHART_COUNTS, ("chart census changed", chart_counts))
    require(safe_checks == len(bodies), ("safe chart count changed", safe_checks))
    require(len(pair_rows) == 13 * 495 == 6435, ("body-pair count changed", len(pair_rows)))
    require(direct_formula_checks == len(pair_rows) * len(PAIR_AUDIT_LEVELS) == 25740,
            ("direct/formula check count changed", direct_formula_checks))
    pair_margin_minimum = min(pair_rows)
    require(pair_margin_minimum[:6] == EXPECTED_PAIR_MARGIN_MINIMUM,
            ("minimum compatible pair margin changed", pair_margin_minimum))
    require(sum(row[0] == pair_margin_minimum[0] for row in pair_rows) == 1,
            ("minimum compatible pair margin is not unique", pair_margin_minimum))

    # Full exact D=1 spanning-tree ledger, peeled by the universal theorem.
    tree_digest = hashlib.sha256()
    original_maximum = None
    original_equality = 0
    peeled_maximum = None
    peeled_equality = 0
    covered = 0
    tree_rows = 0
    for body in bodies:
        for mask in range(1, 63):
            floor_sum, selected_tax, selected_edges = TREE.maximum_floor_tree(body, mask)
            slope = F(45, 2) / floor_sum
            intercept = (F(9, 2) * selected_tax + F(1, 39)) / floor_sum
            threshold = slope + intercept
            witness = monochromatic_consecutive_witness(body, mask)
            row = (threshold, body, mask, floor_sum, selected_tax, slope, intercept, selected_edges, witness)
            tree_digest.update(f"{row}\n".encode())
            leader = (threshold, body, mask)
            if original_maximum is None or threshold > original_maximum[0]:
                original_maximum = leader
                original_equality = 1
            elif threshold == original_maximum[0]:
                original_equality += 1
            if witness is not None:
                covered += 1
            elif peeled_maximum is None or threshold > peeled_maximum[0]:
                peeled_maximum = leader
                peeled_equality = 1
            elif threshold == peeled_maximum[0]:
                peeled_equality += 1
            tree_rows += 1
    require(tree_rows == 3003 * 62 == 186186, ("tree row count changed", tree_rows))
    require(covered == EXPECTED_D1_COVERAGE and tree_rows - covered == 52986,
            ("D=1 adjacent coverage changed", covered, tree_rows - covered))
    require(original_maximum == EXPECTED_ORIGINAL_D1_MAXIMUM and original_equality == 1,
            ("original D=1 maximum changed", original_maximum, original_equality))
    require(peeled_maximum == EXPECTED_PEELED_D1_MAXIMUM and peeled_equality == 1,
            ("peeled D=1 maximum changed", peeled_maximum, peeled_equality))
    require(F(240) < peeled_maximum[0] < F(241), ("peeled D=1 cutoff changed", peeled_maximum))

    hostile_ruler, hostile_safe = R.safe_cell_ranges(HOSTILE_BODY)
    require(hostile_ruler == 168 and hostile_safe == EXPECTED_HOSTILE_SAFE_RANGES,
            ("hostile safe carrier changed", hostile_ruler, hostile_safe))
    hostile_rows = tuple(hostile_finite_profile(m) for m in range(1, HOSTILE_FINITE_MAX_M + 1))
    hostile_worst = max(
        (mass, m, count)
        for m, mass, count, _pair, _debt, _margin, _checks in hostile_rows
    )
    hostile_pair_margin_minimum = min(
        (margin, m)
        for m, _mass, _count, _pair, _debt, margin, _checks in hostile_rows
    )
    require(hostile_worst == EXPECTED_HOSTILE_WORST,
            ("hostile finite worst changed", hostile_worst))
    require(hostile_pair_margin_minimum == EXPECTED_HOSTILE_PAIR_MARGIN_MINIMUM,
            ("hostile finite pair margin changed", hostile_pair_margin_minimum))
    require(sum(row[-1] for row in hostile_rows) == 6 * HOSTILE_FINITE_MAX_M,
            "hostile finite route count changed")

    semantic_payload = (
        chart_counts,
        safe_checks,
        len(pair_rows),
        direct_formula_checks,
        pair_margin_minimum,
        pair_census_digest.hexdigest(),
        tree_rows,
        covered,
        original_maximum,
        original_equality,
        peeled_maximum,
        peeled_equality,
        tree_digest.hexdigest(),
        hostile_rows,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))

    lines = [
        "LRC14 universal reflected adjacent-pair all-level closure and D=1 peel exact referee",
        f"base_source_sha256={sha256(BASE_SOURCE)}",
        f"tree_source_sha256={sha256(TREE_SOURCE)}",
        "THEOREM=arbitrary positive levels; any same-level consecutive labels close; other four levels arbitrary; common Q need not be minimum",
        "CHART_A=max(E)<=12:j=L/14;CHART_B=max(E)>=13:s=terminal-block-start,j=15L/(14s)",
        f"chart_counts={chart_counts};whole_safe_cells={safe_checks}",
        "pair_formula=LQ[Q(L/7-j-1/2)+floor(aj/L)+(3-a)/7]/[(LQ-a)(LQ-a-1)]",
        f"body_pair_rows={len(pair_rows)};direct_reflected_formula_checks={direct_formula_checks};audit_levels={PAIR_AUDIT_LEVELS}",
        f"pair_census_digest_sha256={pair_census_digest.hexdigest()}",
        f"minimum_compatible_margin={qtext(pair_margin_minimum[0])};E={pair_margin_minimum[1]};L={pair_margin_minimum[2]};j={pair_margin_minimum[3]};s={pair_margin_minimum[4]};pair={pair_margin_minimum[5]},{pair_margin_minimum[5]+1};unique=1",
        f"D1_LEDGER;rows={tree_rows};adjacent_pair_covered={covered};remaining={tree_rows-covered};digest_sha256={tree_digest.hexdigest()}",
        f"D1_ORIGINAL_MAX;threshold={qtext(original_maximum[0])};E={original_maximum[1]};mask={original_maximum[2]};unique={original_equality}",
        f"D1_PEELED_MAX;threshold={qtext(peeled_maximum[0])};E={peeled_maximum[1]};mask={peeled_maximum[2]};unique={peeled_equality}",
        "D1_GLOBAL_CONSEQUENCE=all reflected D=1 packets close for m>=241;combined certificate failure forces m<=240",
        f"HOSTILE_CONTROL;E={HOSTILE_BODY};q=(m,m,m,m,m,m+1);j={HOSTILE_J};m=1..{HOSTILE_FINITE_MAX_M};survivors=0",
        f"HOSTILE_WORST;mass={qtext(hostile_worst[0])};m={hostile_worst[1]};gap={qtext(F(6,7)-hostile_worst[0])};components={hostile_worst[2]}",
        f"HOSTILE_PAIR_MARGIN_MINIMUM={qtext(hostile_pair_margin_minimum[0])};m={hostile_pair_margin_minimum[1]}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
