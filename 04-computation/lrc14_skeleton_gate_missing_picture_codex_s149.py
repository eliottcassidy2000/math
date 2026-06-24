#!/usr/bin/env python3
"""S149: compose the LRC14 boundary skeleton, gate, C27, and unital lanes.

This script is a small finite "missing picture" scout.  It does not try to
prove LRC14.  It makes executable the synthesis suggested by recent work:

    Haar/Baire boundary skeleton says AP and GW are endpoint-only.
    Jacobsthal gate says only the 12-site can be a hidden single acceleration.
    C27/unital labels distinguish the legal tight 12->24 branch from the
    open K33 12->36 branch.
    Covering/divisibility catches residue impostors that endpoint residues miss.

Tournament Analysis declaration
-------------------------------
Vertices are proof lenses, not runners:

    covering_divisibility, Haar_Baire_boundary_skeleton, Jacobsthal_gate,
    C27_unital_transfer, derived_relative_profile, K33_state_lift,
    raw_residue_multiset

Pairwise observable:

    detects residue impostors, separates open fronts from boundary-only rows,
    forces the movable 12-site, distinguishes D3 from D9/unit transfers,
    feeds the state-lift endpoint, remains finite-checkable, resists scalar
    collapse.

Switch/gauge:

    A lens beats B if it wins more coordinates in the observable vector.
    Ties follow the displayed Hamiltonian path.

Assumption challenge:

    The quotient is not runner vertices.  Alternative vertices considered are
    runners, residues, gaps, boundary points, interval fronts, cover arcs,
    C27 shell transfers, q=3 unital blocks, derived-profile classes, K33
    obligations, and proof lenses.  This scout chooses proof lenses plus
    labelled row packets.  It preserves the LRC predicate "open witness versus
    endpoint-only residual" and the transfer branch that must feed a proof.  It
    destroys full speed identity and therefore keeps covering and C27 labels as
    external checks.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


N = 14
DELTA = Fraction(1, N)
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])
U = (1, 3, 5, 9, 11, 13)
RESIDUE_TIE_PATH = (1, 13, 3, 11, 5, 9, 7, 2, 12, 4, 10, 6, 8, 0)
RESIDUE_RANK = {r: i for i, r in enumerate(RESIDUE_TIE_PATH)}


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circular_distance_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def split_circle_interval(a: Fraction, b: Fraction) -> list[tuple[Fraction, Fraction]]:
    while a < 0:
        a += 1
        b += 1
    while a >= 1:
        a -= 1
        b -= 1
    if b <= 1:
        return [(a, b)]
    return [(a, Fraction(1)), (Fraction(0), b - 1)]


UNSAFE_CACHE: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {}
ENDPOINT_CACHE: dict[int, tuple[Fraction, ...]] = {}


def unsafe_intervals_for_speed(v: int) -> tuple[tuple[Fraction, Fraction], ...]:
    if v in UNSAFE_CACHE:
        return UNSAFE_CACHE[v]
    radius = DELTA / v
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        center = Fraction(m, v)
        out.extend(split_circle_interval(center - radius, center + radius))
    ans = tuple(sorted(out))
    UNSAFE_CACHE[v] = ans
    return ans


def endpoint_candidates_for_speed(v: int) -> tuple[Fraction, ...]:
    if v in ENDPOINT_CACHE:
        return ENDPOINT_CACHE[v]
    out: list[Fraction] = []
    for m in range(v):
        out.append(frac_part((Fraction(m) - DELTA) / v))
        out.append(frac_part((Fraction(m) + DELTA) / v))
    ans = tuple(sorted(set(out)))
    ENDPOINT_CACHE[v] = ans
    return ans


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[list[Fraction]] = [[intervals[0][0], intervals[0][1]]]
    for a, b in intervals[1:]:
        if a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def complement_intervals(covered: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not covered:
        return [(Fraction(0), Fraction(1))]
    out: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for a, b in covered:
        if cursor < a:
            out.append((cursor, a))
        if b > cursor:
            cursor = b
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def safe_open_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(unsafe_intervals_for_speed(v))
    return complement_intervals(merge_intervals(intervals))


def interval_measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), start=Fraction(0))


def threshold_safe_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    candidates = {Fraction(0)}
    for v in speeds:
        candidates.update(endpoint_candidates_for_speed(v))
    good = []
    for t in sorted(candidates):
        if all(circular_distance_to_integer(Fraction(v) * t) >= DELTA for v in speeds):
            good.append(t)
    return tuple(good)


def active_owners(speeds: tuple[int, ...], t: Fraction) -> tuple[str, ...]:
    owners: list[str] = []
    for v in speeds:
        y = frac_part(Fraction(v) * t)
        if y == DELTA:
            owners.append(f"{v}L")
        elif y == 1 - DELTA:
            owners.append(f"{v}R")
    return tuple(owners)


def midpoint_slack(speeds: tuple[int, ...], a: Fraction, b: Fraction) -> Fraction:
    mid = (a + b) / 2
    return min(circular_distance_to_integer(Fraction(v) * mid) for v in speeds) - DELTA


def clock_distance_num(speed: int, unit: int) -> int:
    r = (speed * unit) % N
    return min(r, N - r)


def endpoint_units_ok(row: tuple[int, ...]) -> bool:
    return all(min(clock_distance_num(s, u) for s in row) == 1 for u in U)


def full_unit_shell(row: tuple[int, ...]) -> bool:
    return set(U).issubset({s % N for s in row})


def has_7_apex(row: tuple[int, ...]) -> bool:
    return any(s % N == 7 for s in row)


def max_empty_gap(row: tuple[int, ...]) -> int:
    occupied = sorted({s % N for s in row})
    if not occupied:
        return N
    return max(b - a for a, b in zip(occupied, occupied[1:] + [occupied[0] + N]))


def covering_ok(row: tuple[int, ...]) -> bool:
    return all(any(s % d == 0 for s in row) for d in range(1, 14)) and all(s % N != 0 for s in row)


def distinct_speeds(row: tuple[int, ...]) -> bool:
    return len(set(row)) == len(row)


def passes_basic(row: tuple[int, ...]) -> bool:
    return (
        distinct_speeds(row)
        and covering_ok(row)
        and endpoint_units_ok(row)
        and full_unit_shell(row)
        and has_7_apex(row)
        and max_empty_gap(row) <= 3
    )


def jacobsthal_window(v: int) -> tuple[int, int]:
    return (N - v, 2 * N - 1 - 2 * v)


def gate_admissible(v: int) -> bool:
    lo, hi = jacobsthal_window(v)
    if lo > hi:
        return False
    return all(gcd(k, v) > 1 for k in range(lo, hi + 1))


def gate_profile(v: int) -> tuple[int, int, int, int]:
    lo, hi = jacobsthal_window(v)
    window = [] if lo > hi else list(range(lo, hi + 1))
    primes = prime_divisors(v)
    hit = {p for p in primes if any(k % p == 0 for k in window)}
    return (int(gate_admissible(v)), len(hit), -len(window), v)


def prime_divisors(n: int) -> set[int]:
    out: set[int] = set()
    d = 2
    x = n
    while d * d <= x:
        while x % d == 0:
            out.add(d)
            x //= d
        d += 1
    if x > 1:
        out.add(x)
    return out


def gate_order() -> list[int]:
    return sorted(range(1, 14), key=gate_profile, reverse=True)


def shell_fold(x: int, modulus: int = 27) -> int:
    r = x % modulus
    return min(r, modulus - r)


def removed_added(row: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    s = set(row)
    return tuple(sorted(set(AP) - s)), tuple(sorted(s - set(AP)))


def transfer_label(row: tuple[int, ...]) -> str:
    holes, adds = removed_added(row)
    if not holes and not adds:
        return "AP"
    htxt = ",".join(f"H{h}:g{gcd(h, 27)}" for h in holes) or "-"
    atxt = ",".join(f"D{shell_fold(a)}:g{gcd(shell_fold(a), 27)}@{a}" for a in adds) or "-"
    return f"{htxt} -> {atxt}"


KNOWN_COMPONENTS = {
    (12, 24): "GW:D3-tight",
    (12, 36): "K33:D9-open",
    (10, 20): "P10:unit-petal",
    (13, 26): "P13:unit-petal",
    (8, 16): "P8:gate-fail",
}


def component_tags(row: tuple[int, ...]) -> tuple[str, ...]:
    holes, adds = removed_added(row)
    tags: list[str] = []
    used_holes: set[int] = set()
    used_adds: set[int] = set()
    for h, a in KNOWN_COMPONENTS:
        if h in holes and a in adds:
            tags.append(KNOWN_COMPONENTS[(h, a)])
            used_holes.add(h)
            used_adds.add(a)
    for h in holes:
        if h in used_holes:
            continue
        candidates = [a for a in adds if a not in used_adds]
        if candidates:
            a = candidates[0]
            used_adds.add(a)
            if a == 2 * h:
                tags.append(f"P{h}:doubled-gate-{'ok' if gate_admissible(h) else 'fail'}")
            else:
                tags.append(f"H{h}->D{shell_fold(a)}:unknown")
        else:
            tags.append(f"H{h}:unmatched")
    for a in adds:
        if a not in used_adds:
            tags.append(f"D{shell_fold(a)}@{a}:unmatched")
    return tuple(tags)


def doubled_gate_sites(row: tuple[int, ...]) -> tuple[tuple[int, tuple[int, int], bool], ...]:
    holes, adds = removed_added(row)
    out: list[tuple[int, tuple[int, int], bool]] = []
    for h in holes:
        if 2 * h in adds:
            out.append((h, jacobsthal_window(h), gate_admissible(h)))
    return tuple(out)


def boundary_skeleton(row: tuple[int, ...]) -> tuple[tuple[Fraction, tuple[str, ...]], ...]:
    if interval_measure(safe_open_components(row)) != 0:
        return ()
    return tuple((t, active_owners(row, t)) for t in threshold_safe_points(row))


AP_BOUNDARY_SKELETON = boundary_skeleton(AP)


def boundary_skeleton_key(row: tuple[int, ...]) -> str:
    skel = boundary_skeleton(row)
    if not skel:
        return "open-or-empty"
    return "AP/GW-six-pair" if skel == AP_BOUNDARY_SKELETON else short_hash(skel)


def first_fronts(row: tuple[int, ...], limit: int = 2) -> tuple[str, ...]:
    comps = safe_open_components(row)
    out: list[str] = []
    for a, b in comps[:limit]:
        out.append(
            f"({a},{b}) {active_owners(row, a)}->{active_owners(row, b)} "
            f"slack={midpoint_slack(row, a, b)}"
        )
    return tuple(out)


def short_hash(obj: object) -> str:
    return hashlib.sha1(repr(obj).encode("utf-8")).hexdigest()[:12]


def relative_profile(outer: int, inner: int) -> tuple[tuple[tuple[int, int, int], int], ...]:
    counts = Counter(
        (
            clock_distance_num(outer, u),
            clock_distance_num(inner, u),
            clock_distance_num(inner - outer, u),
        )
        for u in U
    )
    return tuple(sorted(counts.items()))


def inner_profile_class(row: tuple[int, ...], outer: int) -> tuple[tuple[tuple[tuple[int, int, int], int], int], ...]:
    return tuple(sorted(Counter(relative_profile(outer, inner) for inner in row).items()))


def square_profile_hash_histogram(row: tuple[int, ...]) -> tuple[tuple[str, int], ...]:
    hist = Counter(short_hash(inner_profile_class(row, outer)) for outer in row)
    return tuple(sorted(hist.items()))


@dataclass(frozen=True)
class RowCase:
    name: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class RowAudit:
    case: RowCase
    strict_mass: Fraction
    components: int
    closed_points: int
    covering_ok: bool
    endpoint_units_ok: bool
    full_unit_shell: bool
    has_7_apex: bool
    max_empty_gap: int
    basic_ok: bool
    boundary_key: str
    transfer: str
    tags: tuple[str, ...]
    gates: tuple[tuple[int, tuple[int, int], bool], ...]
    route: str
    first_failure: str
    square_profile_hashes: tuple[tuple[str, int], ...]


def classify_route(row: tuple[int, ...], mass: Fraction, basic: bool) -> str:
    tags = component_tags(row)
    gates = doubled_gate_sites(row)
    if row == AP:
        return "tight_AP_boundary"
    if row == GW:
        return "tight_GW_hidden_D3_boundary"
    if not distinct_speeds(row):
        return "duplicate_speed_not_candidate"
    if not covering_ok(row):
        return "residue_liar_covering_failure"
    if mass == 0:
        if boundary_skeleton(row) == AP_BOUNDARY_SKELETON:
            return "AP_GW_skeleton_unresolved"
        return "other_boundary_skeleton_obligation"
    if any("K33" in tag for tag in tags):
        return "open_K33_state_lift_lane"
    if len(tags) > 1:
        return "open_splice_discharge"
    if any("unit-petal" in tag for tag in tags):
        return "open_unit_petal_discharge"
    if gates and not all(ok for _, _, ok in gates):
        return "open_gate_fail_impostor"
    if basic:
        return "open_basic_survivor_needs_next_label"
    return "open_general_discharge"


def first_failure(row: tuple[int, ...], mass: Fraction) -> str:
    checks = [
        ("distinct_speeds", distinct_speeds(row)),
        ("covering_ok", covering_ok(row)),
        ("endpoint_units_ok", endpoint_units_ok(row)),
        ("full_unit_shell", full_unit_shell(row)),
        ("has_7_apex", has_7_apex(row)),
        ("max_empty_gap<=3", max_empty_gap(row) <= 3),
    ]
    gates = doubled_gate_sites(row)
    if gates:
        checks.append(("all_doubled_sites_gate_admissible", all(ok for _, _, ok in gates)))
    if mass > 0:
        checks.append(("strict_Haar_zero", False))
    for name, ok in checks:
        if not ok:
            return name
    return "none_or_boundary_endpoint"


def audit_case(case: RowCase) -> RowAudit:
    comps = safe_open_components(case.speeds)
    mass = interval_measure(comps)
    pts = threshold_safe_points(case.speeds)
    basic = passes_basic(case.speeds)
    return RowAudit(
        case=case,
        strict_mass=mass,
        components=len(comps),
        closed_points=len(pts),
        covering_ok=covering_ok(case.speeds),
        endpoint_units_ok=endpoint_units_ok(case.speeds),
        full_unit_shell=full_unit_shell(case.speeds),
        has_7_apex=has_7_apex(case.speeds),
        max_empty_gap=max_empty_gap(case.speeds),
        basic_ok=basic,
        boundary_key=boundary_skeleton_key(case.speeds),
        transfer=transfer_label(case.speeds),
        tags=component_tags(case.speeds),
        gates=doubled_gate_sites(case.speeds),
        route=classify_route(case.speeds, mass, basic),
        first_failure=first_failure(case.speeds, mass),
        square_profile_hashes=square_profile_hash_histogram(case.speeds),
    )


def accelerate(*sites: int) -> tuple[int, ...]:
    values = list(AP)
    for v in sites:
        values.remove(v)
    values.extend(2 * v for v in sites)
    return tuple(sorted(values))


NAMED_CASES = [
    RowCase("AP", AP, "baseline tight boundary"),
    RowCase("GW_12_to_24", GW, "Goddyn-Wong hidden D3 transfer"),
    RowCase("false_AP_residue_12_to_26", tuple(list(range(1, 12)) + [13, 26]), "same AP residues, lost multiple of 12"),
    RowCase("near_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "first D9/K33 Farey child"),
    RowCase("loose_gate_8_to_16", accelerate(8), "basic filters pass but Jacobsthal gate fails"),
    RowCase("petal_10_to_20", accelerate(10), "unit C27 petal"),
    RowCase("petal_13_to_26", accelerate(13), "unit C27 petal"),
    RowCase("splice_10_12_to_20_24", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]), "P10 plus tight GW"),
    RowCase("splice_10_12_to_20_36", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]), "P10 plus K33"),
]


def print_gate_table() -> None:
    print("[1] Jacobsthal acceleration gate")
    print(f"  gate_order={tuple(gate_order())}")
    print(f"  admissible_sites={[v for v in range(1, 14) if gate_admissible(v)]}")
    for v in range(1, 14):
        print(f"  v={v:2d} W={jacobsthal_window(v)} profile={gate_profile(v)} gate={gate_admissible(v)}")


def print_named_cases() -> None:
    print()
    print("[2] Named skeleton-gate row atlas")
    for audit in map(audit_case, NAMED_CASES):
        print(f"  {audit.case.name}:")
        print(f"    row={audit.case.speeds}")
        print(f"    note={audit.case.note}")
        print(
            f"    strict_mass={audit.strict_mass} comps={audit.components} "
            f"closed_pts={audit.closed_points} boundary_key={audit.boundary_key}"
        )
        print(
            f"    covering={audit.covering_ok} endpoint_units={audit.endpoint_units_ok} "
            f"unit_shell={audit.full_unit_shell} apex7={audit.has_7_apex} "
            f"gap={audit.max_empty_gap} basic={audit.basic_ok}"
        )
        print(f"    transfer={audit.transfer}")
        print(f"    tags={audit.tags or ('-',)} gates={audit.gates or ('-',)}")
        print(f"    first_failure={audit.first_failure} route={audit.route}")
        print(f"    square_profile_hashes={audit.square_profile_hashes}")
        for front in first_fronts(audit.case.speeds):
            print(f"      front {front}")


def print_family_audit() -> None:
    print()
    print("[3] Finite single/double acceleration family audit")
    single_audits = {v: audit_case(RowCase(f"single_{v}_to_{2*v}", accelerate(v), "")) for v in range(1, 14)}
    double_audits = {
        sites: audit_case(RowCase(f"double_{sites[0]}_{sites[1]}", accelerate(*sites), ""))
        for sites in combinations(range(1, 14), 2)
    }
    single_basic = [v for v, a in single_audits.items() if a.basic_ok]
    single_gate = [v for v, a in single_audits.items() if a.basic_ok and all(ok for _, _, ok in a.gates)]
    double_basic = [sites for sites, a in double_audits.items() if a.basic_ok]
    single_routes = Counter(a.route for a in single_audits.values())
    double_routes = Counter(a.route for a in double_audits.values())
    print(f"  single_accelerations={len(single_audits)} route_hist={dict(sorted(single_routes.items()))}")
    print(f"  single_basic_pass_sites={single_basic}")
    print(f"  single_basic_and_gate_pass_sites={single_gate}")
    for v in single_basic:
        a = single_audits[v]
        print(
            f"    single {v:2d}->{2*v:2d}: mass={str(a.strict_mass):>8} "
            f"gate={a.gates} route={a.route}"
        )
    print(f"  double_accelerations={len(double_audits)} route_hist={dict(sorted(double_routes.items()))}")
    print(f"  double_basic_pass_sites={double_basic}")
    for sites in double_basic:
        a = double_audits[sites]
        print(
            f"    double {sites}: mass={str(a.strict_mass):>8} gates={a.gates or ('-',)} "
            f"tags={a.tags or ('-',)} route={a.route}"
        )
    boundary_double = [sites for sites, a in double_audits.items() if a.strict_mass == 0 and a.closed_points]
    print(f"  double_boundary_only_sites={boundary_double}")


@dataclass(frozen=True)
class Lens:
    name: str
    vector: tuple[int, int, int, int, int, int, int]
    note: str


LENSES = [
    Lens("covering_divisibility", (5, 1, 2, 1, 1, 5, 5), "catches false AP residues such as 12->26"),
    Lens("Haar_Baire_boundary_skeleton", (3, 5, 2, 2, 2, 5, 5), "splits endpoint-only rows from open fronts"),
    Lens("Jacobsthal_gate", (2, 2, 5, 2, 2, 5, 4), "forces the only movable single site v=12"),
    Lens("C27_unital_transfer", (2, 3, 4, 5, 4, 4, 5), "separates D3, D9, and unit-petal transfer branches"),
    Lens("derived_relative_profile", (2, 2, 3, 4, 3, 4, 4), "detects first-fold profile classes after coarse shells fail"),
    Lens("K33_state_lift", (1, 3, 2, 4, 5, 3, 5), "routes D9/Farey child packets to HYP-2908/THM-572"),
    Lens("raw_residue_multiset", (0, 1, 1, 1, 0, 5, 1), "useful ledger, but magnitude-blind and easy to spoof"),
]


def print_lens_tournament() -> None:
    print()
    print("[4] Proof-lens tournament")
    names = [lens.name for lens in LENSES]
    tie = {name: i for i, name in enumerate(names)}
    adj = [[False] * len(LENSES) for _ in LENSES]
    for i, a in enumerate(LENSES):
        for j, b in enumerate(LENSES):
            if i == j:
                continue
            wins = sum(x > y for x, y in zip(a.vector, b.vector))
            losses = sum(x < y for x, y in zip(a.vector, b.vector))
            adj[i][j] = tie[a.name] < tie[b.name] if wins == losses else wins > losses
    scores = [sum(row) for row in adj]
    c3 = 0
    cycle_examples: list[tuple[str, str, str]] = []
    for a, b, c in combinations(range(len(LENSES)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
            if len(cycle_examples) < 3:
                if adj[a][b] and adj[b][c] and adj[c][a]:
                    cycle_examples.append((names[a], names[b], names[c]))
                else:
                    cycle_examples.append((names[a], names[c], names[b]))
    order = [name for _score, name in sorted(zip(scores, names), reverse=True)]
    hp_count = sum(
        all(adj[path[i]][path[i + 1]] for i in range(len(path) - 1))
        for path in permutations(range(len(LENSES)))
    )
    reach = [[adj[i][j] or i == j for j in range(len(LENSES))] for i in range(len(LENSES))]
    for k in range(len(LENSES)):
        for i in range(len(LENSES)):
            for j in range(len(LENSES)):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    seen: set[int] = set()
    sccs: list[tuple[str, ...]] = []
    for i in range(len(LENSES)):
        if i in seen:
            continue
        comp = tuple(names[j] for j in range(len(LENSES)) if reach[i][j] and reach[j][i])
        seen.update(j for j in range(len(LENSES)) if reach[i][j] and reach[j][i])
        sccs.append(comp)
    print("  criteria: residue_liar, open_boundary, site12, D3_D9_unit, state_lift, finite_check, anti_scalar")
    for lens, score in zip(LENSES, scores):
        print(f"  {lens.name:30s} vector={lens.vector} score={score} note={lens.note}")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))} directed_3cycles={c3} hp_count={hp_count}")
    print(f"  sccs={sccs}")
    if cycle_examples:
        print(f"  cycle_examples={cycle_examples}")
    print(f"  Hamiltonian_path={' > '.join(order)}")


def print_readout() -> None:
    print()
    print("[5] Missing-picture synthesis")
    print("  The finite data supports a two-switch proof interface:")
    print("    switch A: strict Haar mass > 0 gives an open witness interval.")
    print("    switch B: strict Haar mass = 0 routes to boundary-owner skeletons.")
    print("  AP and GW share the same boundary-owner skeleton, so Haar/Baire alone")
    print("  cannot see the H12->D3 transfer.  The transfer must be carried by")
    print("  C27/unital labels.")
    print("  The Jacobsthal gate then says the only single hidden acceleration site")
    print("  compatible with the AP/GW summit is v=12.  Among the 12-branches,")
    print("  C27/unital labels split D3 (GW tight boundary) from D9 (open K33 lane).")
    print("  Residue impostors such as 12->26 show why covering/divisibility must")
    print("  precede raw residue or apex-tournament conclusions.")
    print()
    print("  Candidate theorem shape:")
    print("    after standard reductions, any primitive LRC14 endpoint atom with no")
    print("    strict Haar witness has the AP/GW six-pair boundary skeleton; the")
    print("    derived gate forces any hidden single acceleration to be 12; C27/unital")
    print("    pair completion permits only the D3 tight branch, while D9 and unit")
    print("    petals open Haar fronts or feed the K33/state-lift route.")
    print()
    print("  What remains global:")
    print("    prove every bad atom enters this skeleton-gate language, or prove the")
    print("    complementary open-front/state-lift discharge for atoms outside it.")


def main() -> None:
    print("LRC14 skeleton-gate missing-picture scout")
    print("=" * 72)
    print_gate_table()
    print_named_cases()
    print_family_audit()
    print_lens_tournament()
    print_readout()


if __name__ == "__main__":
    main()
