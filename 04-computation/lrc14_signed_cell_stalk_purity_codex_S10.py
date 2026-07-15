#!/usr/bin/env python3
"""Exact replay for THM-821, the signed-cell stalk-purity boundary.

The fixed numerical predicate is the THM-803 erosion predicate for
``(x,y)=(13,5)``.  Its folded frequencies are ``(a,b)=(9,4)``.  For every
closed deep component ``E`` and every signed return cell ``R`` we form the
raw circular sum arc ``E+R`` and record

    success(E,R) <=> min_{t in E+R} (||9t||+||4t||) >= 11/13.

The audit reconstructs the 213 deterministic THM-817 random rows from seed
803807 and adds every one of the 438 x 3 stalks of its disconnected U_0 row.
It tests which forgetful signatures have mixed success fibres.  Exact sum
arcs decide this fixed predicate; endpoint owners are retained separately as
ancestry needed to reconstruct cells under later descent operations.

This is a finite exact audit.  It does not assert an all-size purity theorem,
a global row-level separator, or any Cech gluing statement.

Implementation guardrail
------------------------
The max-speed cell and deep-component constructors are imported from the
independently replayed THM-817 and THM-803 artifacts.  This file independently
reconstructs endpoint owners, circular sum arcs, folded selector events,
margins, signature fibres, counterexamples, and the U_n local family stalk.
Its bank and full-stalk digests therefore catch drift in either dependency.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, fields, is_dataclass
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import floor, gcd
from pathlib import Path
from random import Random
import sys
from typing import Callable, Hashable, Iterable, Sequence


HERE = Path(__file__).resolve().parent
X, Y = 13, 5
A, B_FOLD = (X + Y) // 2, abs(X - Y) // 2
BETA = F(1, 11)
GAMMA = F(2, 143)
THRESHOLD = F(11, 13)
EXPECTED_BANK_DIGEST = "303009db5bf61e2b5584f0664e740039aefda134e8ba80cf34de30cd897fcc71"
EXPECTED_CERTIFICATE_DIGEST = "ebfe916e128a632c0e61bc5891df396f054b925747c56ca99ac7d3f252db4f8e"


def load_dependency(module_name: str, filename: str):
    path = HERE / filename
    spec = spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load dependency {path}")
    module = module_from_spec(spec)
    # dataclasses inspect sys.modules while decorating ReturnCell.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


THM817 = load_dependency(
    "thm817_return_cells_for_thm821",
    "lrc13_return_satellites_cell_classifier_codex_S10.py",
)
THM803 = load_dependency(
    "thm803_selector_for_thm821",
    "lrc13_antigrid_all_component_selector_codex_S10.py",
)

assert THM817.GAMMA == GAMMA
assert THM803.BETA == BETA
assert THM803.THRESHOLD == THRESHOLD
assert THM817.SATELLITE_ROW == (1, 2, 3, 4, 7, 500, 503, 504, 505, 506)


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def folded_value(t: F) -> F:
    value = norm(A * t) + norm(B_FOLD * t)
    assert value == THM803.q_value(X, Y, t)
    return value


def canonical(value: object) -> str:
    if isinstance(value, F):
        return fmt(value)
    if isinstance(value, bool):
        return "1" if value else "0"
    if value is None:
        return "none"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if is_dataclass(value):
        return canonical(tuple(getattr(value, field.name) for field in fields(value)))
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{canonical(key)}:{canonical(value[key])}"
            for key in sorted(value, key=lambda item: canonical(item))
        ) + "}"
    if isinstance(value, (tuple, list)):
        return "(" + ",".join(canonical(item) for item in value) + ")"
    raise TypeError(f"no canonical encoding for {type(value)!r}")


def digest(value: object) -> str:
    return sha256(canonical(value).encode()).hexdigest()


@dataclass(frozen=True)
class SelectorMinimum:
    value: F
    argmins: tuple[F, ...]
    events: tuple[tuple[F, tuple[str, ...]], ...]


@dataclass(frozen=True)
class Stalk:
    source: str
    row_index: int
    speeds: tuple[int, ...]
    deep_index: int
    cell_label: int
    return_interval: tuple[F, F]
    return_owners: tuple[tuple[int, ...], tuple[int, ...]]
    deep_interval: tuple[F, F]
    deep_owners: tuple[tuple[int, ...], tuple[int, ...]]
    sum_arc: tuple[F, F]  # (start modulo one, circular width)
    selector: SelectorMinimum
    margin: F            # min Q - 11/13
    success: bool


@dataclass(frozen=True)
class FibreAudit:
    fibres: int
    mixed: int
    largest_mixed: int
    largest_bucket: int


def selector_minimum(left: F, right: F) -> SelectorMinimum:
    """Exact minimum of Q_(9,4) on one lifted interval of width < 1."""

    assert left <= right and right - left < 1
    candidates = {left, right}
    for frequency in (A, B_FOLD):
        for k in range(floor(2 * frequency * left) - 2,
                       floor(2 * frequency * right) + 4):
            point = F(k, 2 * frequency)
            if left <= point <= right:
                candidates.add(point)

    value = min(folded_value(point) for point in candidates)
    minimizing_lifts = tuple(
        sorted(point for point in candidates if folded_value(point) == value)
    )
    events = []
    for point in minimizing_lifts:
        tags = []
        if point == left:
            tags.append("L")
        if point == right:
            tags.append("R")
        if (2 * A * (point % 1)).denominator == 1:
            tags.append("A")
        if (2 * B_FOLD * (point % 1)).denominator == 1:
            tags.append("B")
        events.append((point % 1, tuple(tags)))
    return SelectorMinimum(
        value,
        tuple(sorted({point % 1 for point in minimizing_lifts})),
        tuple(events),
    )


def endpoint_owners(
    speeds: tuple[int, ...], interval: tuple[F, F]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left, right = interval
    owners = (
        tuple(speed for speed in speeds if norm(speed * left) == BETA),
        tuple(speed for speed in speeds if norm(speed * right) == BETA),
    )
    assert owners[0] and owners[1]
    return owners


def make_stalks(
    source: str, row_index: int, speeds: tuple[int, ...]
) -> tuple[tuple[Stalk, ...], int, int]:
    cells = THM817.return_cells(speeds)
    direct = THM817.direct_return_components(speeds)
    assert THM817.cell_intervals(cells) == direct
    deep = THM803.deep_components(speeds)

    answer = []
    for deep_index, deep_interval in enumerate(deep):
        deep_owner_pair = endpoint_owners(speeds, deep_interval)
        for cell in cells:
            left = deep_interval[0] + cell.left
            right = deep_interval[1] + cell.right
            width = right - left
            assert 0 <= width < 1
            selector = selector_minimum(left, right)
            margin = selector.value - THRESHOLD
            answer.append(
                Stalk(
                    source,
                    row_index,
                    speeds,
                    deep_index,
                    cell.label,
                    (cell.left, cell.right),
                    (cell.left_owners, cell.right_owners),
                    deep_interval,
                    deep_owner_pair,
                    (left % 1, width),
                    selector,
                    margin,
                    margin >= 0,
                )
            )
    return tuple(answer), len(deep), len(cells)


def deterministic_random_rows() -> tuple[tuple[int, ...], ...]:
    rng = Random(803807)
    rows = []
    for maximum in range(10, 81):
        for _ in range(3):
            row = tuple(sorted((*rng.sample(range(1, maximum), 9), maximum)))
            if gcd(*row) == 1:
                rows.append(row)
    assert len(rows) == 213
    return tuple(rows)


def audit_signature(
    stalks: Sequence[Stalk], key: Callable[[Stalk], Hashable]
) -> FibreAudit:
    fibres: dict[Hashable, list[bool]] = defaultdict(list)
    for stalk in stalks:
        fibres[key(stalk)].append(stalk.success)
    mixed = [values for values in fibres.values() if len(set(values)) > 1]
    return FibreAudit(
        len(fibres),
        len(mixed),
        max((len(values) for values in mixed), default=0),
        max((len(values) for values in fibres.values()), default=0),
    )


SIGNATURES: tuple[tuple[str, int, Callable[[Stalk], Hashable]], ...] = (
    ("cell_label", 1, lambda s: s.cell_label),
    ("return_interval", 2, lambda s: s.return_interval),
    ("return_owners", 2, lambda s: s.return_owners),
    ("signed_cell_full", 5,
     lambda s: (s.cell_label, s.return_interval, s.return_owners)),
    ("deep_width", 1,
     lambda s: s.deep_interval[1] - s.deep_interval[0]),
    ("deep_interval", 2, lambda s: s.deep_interval),
    ("deep_owners", 2, lambda s: s.deep_owners),
    ("deep_interval_owners", 4,
     lambda s: (s.deep_interval, s.deep_owners)),
    ("cell_label_plus_deep", 3,
     lambda s: (s.cell_label, s.deep_interval)),
    ("return_plus_deep", 4,
     lambda s: (s.return_interval, s.deep_interval)),
    ("full_input_stalk", 9,
     lambda s: (s.cell_label, s.return_interval, s.return_owners,
                s.deep_interval, s.deep_owners)),
    ("sum_width", 1, lambda s: s.sum_arc[1]),
    ("sum_arc", 2, lambda s: s.sum_arc),
    ("selector_event_shape", 1,
     lambda s: tuple(tags for _point, tags in s.selector.events)),
    ("selector_argmin", 1, lambda s: s.selector.argmins),
    ("exact_margin", 1, lambda s: s.margin),
)


EXPECTED_COMBINED_AUDIT = {
    "cell_label": FibreAudit(3, 3, 9098, 9098),
    "return_interval": FibreAudit(74, 61, 438, 438),
    "return_owners": FibreAudit(74, 61, 438, 438),
    "signed_cell_full": FibreAudit(74, 61, 438, 438),
    "deep_width": FibreAudit(2718, 65, 668, 668),
    "deep_interval": FibreAudit(7290, 24, 3, 38),
    "deep_owners": FibreAudit(3519, 234, 144, 144),
    "deep_interval_owners": FibreAudit(7808, 24, 3, 7),
    "cell_label_plus_deep": FibreAudit(8166, 0, 0, 38),
    "return_plus_deep": FibreAudit(9796, 0, 0, 3),
    "full_input_stalk": FibreAudit(9850, 0, 0, 3),
    "sum_width": FibreAudit(4088, 48, 28, 38),
    "sum_arc": FibreAudit(9796, 0, 0, 3),
    "selector_event_shape": FibreAudit(3, 2, 4804, 4804),
    "selector_argmin": FibreAudit(9086, 0, 0, 60),
    "exact_margin": FibreAudit(4397, 0, 0, 126),
}


def find_stalk(
    stalks: Iterable[Stalk],
    *,
    speeds: tuple[int, ...],
    deep_interval: tuple[F, F],
    cell_label: int,
) -> Stalk:
    matches = [
        stalk for stalk in stalks
        if stalk.speeds == speeds
        and stalk.deep_interval == deep_interval
        and stalk.cell_label == cell_label
    ]
    assert len(matches) == 1
    return matches[0]


def local_deep_component(
    speeds: tuple[int, ...], point: F
) -> tuple[tuple[F, F], tuple[tuple[int, ...], tuple[int, ...]]]:
    left, right = F(0), F(1)
    constraints = []
    for speed in speeds:
        scaled = speed * point
        integer = floor(scaled)
        residue = scaled - integer
        assert BETA <= residue <= 1 - BETA
        lo = (F(integer) + BETA) / speed
        hi = (F(integer + 1) - BETA) / speed
        constraints.append((speed, lo, hi))
        left = max(left, lo)
        right = min(right, hi)
    interval = (left, right)
    owners = (
        tuple(speed for speed, lo, _hi in constraints if lo == left),
        tuple(speed for speed, _lo, hi in constraints if hi == right),
    )
    assert left <= point <= right and owners[0] and owners[1]
    return interval, owners


def family_local_audit() -> tuple[tuple[object, ...], ...]:
    """Audit one exact local E-stalk and representative satellites for U_n."""

    point = F(479, 616)
    table = []
    for index in range(6):
        maximum = THM817.family_B(index)
        speeds = THM817.family_row(index)
        pairs = THM817.family_satellite_pairs(maximum)
        assert pairs == 1 + 720 * index
        component_count = 1 + 2 * pairs
        assert component_count == 3 + 1440 * index

        deep, deep_owners = local_deep_component(speeds, point)
        expected_deep = (
            point - F(25, 616 * (maximum - 3)),
            point,
        )
        assert deep == expected_deep
        assert deep_owners == ((maximum - 3,), (maximum - 2,))

        central = (-F(2, 143 * maximum), F(2, 143 * maximum))
        central_selector = selector_minimum(
            deep[0] + central[0], deep[1] + central[1]
        )
        expected_value = F(69, 616) - F(10, 143 * maximum)
        expected_argmin = (point + F(2, 143 * maximum),)
        assert central_selector.value == expected_value
        assert central_selector.argmins == expected_argmin
        central_margin = expected_value - THRESHOLD
        assert central_margin < 0

        satellite_margins = []
        for label in sorted({1, pairs}):
            satellite = THM817.family_positive_interval(maximum, label)
            selector = selector_minimum(
                deep[0] + satellite[0], deep[1] + satellite[1]
            )
            assert selector.value < THRESHOLD
            satellite_margins.append((label, selector.value - THRESHOLD))

        table.append(
            (
                index,
                maximum,
                component_count,
                deep,
                deep_owners,
                central,
                ((maximum,), (maximum,)),
                central_margin,
                tuple(satellite_margins),
            )
        )
    return tuple(table)


def tournament_telemetry(
    audits: dict[str, FibreAudit]
) -> tuple[tuple[str, ...], tuple[str, ...], int]:
    """Two transitive priorities on signature carriers, not on runners."""

    declaration_rank = {name: index for index, (name, _cost, _key) in enumerate(SIGNATURES)}
    cost = {name: payload for name, payload, _key in SIGNATURES}
    names = tuple(declaration_rank)
    purity_order = tuple(sorted(
        names,
        key=lambda name: (
            audits[name].mixed,
            audits[name].largest_mixed,
            cost[name],
            declaration_rank[name],
        ),
    ))
    compression_order = tuple(sorted(
        names,
        key=lambda name: (
            cost[name],
            audits[name].mixed,
            audits[name].largest_mixed,
            declaration_rank[name],
        ),
    ))

    rank_a = {name: index for index, name in enumerate(purity_order)}
    rank_b = {name: index for index, name in enumerate(compression_order)}
    flips = sum(
        (rank_a[left] < rank_a[right]) != (rank_b[left] < rank_b[right])
        for index, left in enumerate(names)
        for right in names[index + 1:]
    )
    return purity_order, compression_order, flips


def main() -> None:
    rows = deterministic_random_rows()
    bank_digest = digest(rows)
    if EXPECTED_BANK_DIGEST != "TO_BE_FILLED":
        assert bank_digest == EXPECTED_BANK_DIGEST

    random_stalks = []
    deep_histogram: Counter[int] = Counter()
    return_histogram: Counter[int] = Counter()
    row_success = []
    for row_index, speeds in enumerate(rows):
        stalks, deep_count, return_count = make_stalks(
            "random", row_index, speeds
        )
        random_stalks.extend(stalks)
        deep_histogram[deep_count] += 1
        return_histogram[return_count] += 1
        row_success.append(all(stalk.success for stalk in stalks))

    assert len(random_stalks) == 8660
    assert Counter(stalk.success for stalk in random_stalks) == {
        False: 8318,
        True: 342,
    }
    assert return_histogram == {1: 213}
    assert not any(row_success)

    satellite_stalks, satellite_deep_count, satellite_return_count = make_stalks(
        "U0", len(rows), THM817.SATELLITE_ROW
    )
    assert satellite_deep_count == 438
    assert satellite_return_count == 3
    assert len(satellite_stalks) == 1314
    satellite_by_label = {
        label: Counter(
            stalk.success
            for stalk in satellite_stalks
            if stalk.cell_label == label
        )
        for label in (-1, 0, 1)
    }
    assert satellite_by_label == {
        -1: {False: 388, True: 50},
        0: {False: 388, True: 50},
        1: {False: 388, True: 50},
    }

    combined = tuple((*random_stalks, *satellite_stalks))
    assert len(combined) == 9974
    assert Counter(stalk.success for stalk in combined) == {
        False: 9482,
        True: 492,
    }

    audits = {
        name: audit_signature(combined, key)
        for name, _payload, key in SIGNATURES
    }
    assert audits == EXPECTED_COMBINED_AUDIT

    # Same full central signed cell, opposite E-components and verdicts.
    central_row = (30, 33, 35, 42, 50, 53, 55, 72, 73, 75)
    return_fail = find_stalk(
        combined,
        speeds=central_row,
        deep_interval=(F(1, 330), F(2, 165)),
        cell_label=0,
    )
    return_pass = find_stalk(
        combined,
        speeds=central_row,
        deep_interval=(F(309, 803), F(106, 275)),
        cell_label=0,
    )
    assert return_fail.return_interval == return_pass.return_interval == (
        -F(2, 10725), F(2, 10725)
    )
    assert return_fail.return_owners == return_pass.return_owners == (
        (75,), (75,)
    )
    assert return_fail.margin == -F(17357, 21450) and not return_fail.success
    assert return_pass.margin == F(12049, 156585) and return_pass.success

    # Same exact E-component, opposite signed satellites and verdicts.
    shared_deep = (F(49, 132), F(2067, 5566))
    deep_fail = find_stalk(
        combined,
        speeds=THM817.SATELLITE_ROW,
        deep_interval=shared_deep,
        cell_label=-1,
    )
    deep_pass = find_stalk(
        combined,
        speeds=THM817.SATELLITE_ROW,
        deep_interval=shared_deep,
        cell_label=1,
    )
    assert deep_fail.margin == -F(557, 12012) and not deep_fail.success
    assert deep_pass.margin == F(281, 53625) and deep_pass.success

    family = family_local_audit()
    tournament = tournament_telemetry(audits)

    certificate = (
        rows,
        tuple(random_stalks),
        tuple(satellite_stalks),
        audits,
        return_fail,
        return_pass,
        deep_fail,
        deep_pass,
        family,
        tournament,
    )
    certificate_digest = digest(certificate)
    if EXPECTED_CERTIFICATE_DIGEST != "TO_BE_FILLED":
        assert certificate_digest == EXPECTED_CERTIFICATE_DIGEST

    print("THM-821 SIGNED RETURN-CELL / DEEP-COMPONENT STALK PURITY")
    print("arithmetic=integer+fractions.Fraction floating_point=none")
    print("fixed_pair=(13,5) folded_frequencies=(9,4) threshold=11/13")
    print("atomic_predicate=min_(t in E_component+return_cell) Q_(9,4)(t)>=11/13")
    print("scope=finite_exact_atomic_audit_not_all_size_not_global_row_separator_not_Cech_gluing")
    print()
    print("RANDOM_THM817_BANK")
    print(f"seed=803807 rows={len(rows)} bank_sha256={bank_digest}")
    print(f"return_component_histogram={dict(sorted(return_histogram.items()))}")
    print(f"deep_component_histogram={dict(sorted(deep_histogram.items()))}")
    print("atoms=8660 success=342 failure=8318 global_row_success=0 global_row_failure=213")
    print("note=global_row_predicate_is_constant_false_so_only_atomic_fibres_are_nontrivial")
    print()
    print("DISCONNECTED_U0")
    print(f"U0={THM817.SATELLITE_ROW} deep_components=438 return_cells=3 atoms=1314")
    for label in (-1, 0, 1):
        counts = satellite_by_label[label]
        print(f"cell_label={label} success={counts[True]} failure={counts[False]}")
    print()
    print("COMBINED_SIGNATURE_FIBRES")
    print("atoms=9974 success=492 failure=9482")
    print("columns=signature,fibres,mixed,largest_mixed,largest_bucket")
    for name, _payload, _key in SIGNATURES:
        row = audits[name]
        print(
            f"{name},{row.fibres},{row.mixed},"
            f"{row.largest_mixed},{row.largest_bucket}"
        )
    print()
    print("MIXED_FIBRE_RETURN_CELL")
    print(f"U={central_row}")
    print(
        "shared_cell="
        f"[{fmt(return_fail.return_interval[0])},{fmt(return_fail.return_interval[1])}] "
        f"owners={return_fail.return_owners}"
    )
    for role, stalk in (("failure", return_fail), ("success", return_pass)):
        print(
            f"{role}_E=[{fmt(stalk.deep_interval[0])},{fmt(stalk.deep_interval[1])}] "
            f"sum_arc=({fmt(stalk.sum_arc[0])},{fmt(stalk.sum_arc[1])}) "
            f"argmin={canonical(stalk.selector.argmins)} margin={fmt(stalk.margin)}"
        )
    print()
    print("MIXED_FIBRE_DEEP_COMPONENT")
    print(f"U0={THM817.SATELLITE_ROW}")
    print(f"shared_E=[{fmt(shared_deep[0])},{fmt(shared_deep[1])}]")
    for role, stalk in (("failure", deep_fail), ("success", deep_pass)):
        print(
            f"{role}_k={stalk.cell_label} "
            f"R=[{fmt(stalk.return_interval[0])},{fmt(stalk.return_interval[1])}] "
            f"owners={stalk.return_owners} "
            f"sum_arc=({fmt(stalk.sum_arc[0])},{fmt(stalk.sum_arc[1])}) "
            f"argmin={canonical(stalk.selector.argmins)} margin={fmt(stalk.margin)}"
        )
    print()
    print("UN_FAMILY_LOCAL_STALK")
    print("s=479/616 E=[s-25/(616*(B-3)),s] E_owners=(B-3,B-2)")
    print("central_R=[-2/(143B),2/(143B)] R_owners=(B,B)")
    print("central_min_Q=69/616-10/(143B)<11/13 checked_n=0..5")
    print("first_and_last_positive_satellite_stalks=failure checked_n=0..5")
    print("family_component_law=N_R=3+1440n")
    print()
    print("TOURNAMENT_TELEMETRY")
    print("vertices=forgetful_signatures pair_observable=(mixed,largest_mixed,payload)")
    print("switch=payload_first tie_path=declaration_order")
    print(f"purity_path={tournament[0]}")
    print(f"compression_path={tournament[1]}")
    print(f"edge_flips={tournament[2]}")
    print("both_tournaments=transitive score_histogram=0..15 cycles=0 SCCs=16_singletons HP=1")
    print("theorem_carrier=exact_return_interval+exact_deep_interval_or_exact_sum_arc")
    print("endpoint_owners=numerically_forgettable_here_but_retained_for_ancestry_and_descent")
    print("Cech_status=not_tested_no_overlap_or_gluing_data_emitted")
    print()
    print(f"certificate_sha256={certificate_digest}")
    print("PASS: cell/R/E/width/event projections mix; exact input/sum/argmin/margin signatures are pure")


if __name__ == "__main__":
    main()
