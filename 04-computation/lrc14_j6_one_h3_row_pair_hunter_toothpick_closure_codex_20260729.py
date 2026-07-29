#!/usr/bin/env python3
r"""Exact closure of the 42 children left by the one-H3-row top-four test.

The input is the locked THM-2912 census.  Its 42 open ordered children
already retain the parent prefix, the first hostile centre, and every
earlier hostile centre as forbidden labels.  For each literal child R this
verifier:

1. seals the exact global allowed pair cap beta_2(R);
2. applies the four-slot Hunter-star envelope

       G_4=max_a [a+sum_(j=2)^4 min(a,q_j,beta_2-a)];

3. when G_4 >= |R|, solves the least hostile level lambda exactly,
   enumerates the finite second-centre core c_R(w)>=lambda, and allocates
   covers to the earliest maximum-coverage centre; and
4. after subtracting that centre, seals the exact top three allowed
   singleton coverages and checks that their sum is below the remaining
   mass.

All tails use the strict THM-735 discrepancy bound.  The second-centre
step is a literal smaller copy of the first ordered hostile-centre step:
the proof object is a finite rooted "toothpick", not an unmarked residual.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = (
    ROOT
    / "04-computation/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py"
)
BASE_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out"
)
THM2911_VERIFY = ROOT / "04-computation/lrc14_j6_ranked_h1_thm2911/verify.py"
THM2911_OUT = ROOT / "05-knowledge/results/lrc14_j6_ranked_h1_thm2911/locked.out"

BASE_SHA256 = "641ed154aabe7eaa156c0ea27c3be5a5af21a56c5818ec18024f5ab069deb9d5"
BASE_OUT_SHA256 = "ba6bbf7d0c778fd398830663ed543ca18ea1f11369fa3a450a1bb61b15886578"
THM2911_VERIFY_SHA256 = (
    "48f279b90955592f34656eb45e6d6dfaf4f88618a157032e58fcbd7beb104aee"
)
THM2911_OUT_SHA256 = (
    "e5c58cc2eb325928c00839c2593450ea7cce8945b3835898ec83c6c5f42fac9b"
)

FIRST_EXTERNAL = 15
PAIR_HORIZON = 2_500
TOP3_HORIZON = 2_500

EXPECTED_COUNTS: tuple[object, ...] | None = (
    11_842,
    11_511,
    3_344,
    210,
    807,
    765,
    42,
    22,
    20,
    10,
    10,
    29,
    29,
    38,
    1,
    37,
    181,
    314,
    351,
    3_081,
    104_320,
    342,
    2_931,
    71_974,
    2,
    5,
    0,
    0,
    0,
)
EXPECTED_FAILED_ROOT_DIGEST: str | None = (
    "62f177cca38cfa9dc8a1aa33c9842efa61ab7c9a4270cc024f248cb3ee5fd37b"
)
EXPECTED_ADDITIVE_ROOT_DIGEST: str | None = (
    "35dd2a9391ec0ba4e4cc95821382dd1ca8f08dc2571ca1c9ae6d25855beb8f54"
)
EXPECTED_LEDGER_SHA256: str | None = (
    "236df90a5497eb85b6536b9cbfebc1e8a7b30990fe9bd24d3a890d748c3bc428"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    # Git may materialize the same text blob with LF or CRLF line endings.
    # Path.read_text() applies universal-newline normalization, so dependency
    # pins describe mathematical text rather than checkout policy.
    return hashlib.sha256(path.read_text().encode()).hexdigest()


def load_base():
    require(file_sha256(BASE_PATH) == BASE_SHA256, "THM2912 source changed")
    spec = importlib.util.spec_from_file_location("thm2912_toothpick", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM2912")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = load_base()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def interval_mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def literal_output_value(path: Path, prefix: str) -> object:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return ast.literal_eval(matches[0])


def dependency_controls() -> tuple[set[tuple[int, ...]], set[tuple[int, ...]]]:
    require(file_sha256(BASE_OUT) == BASE_OUT_SHA256, "THM2912 output changed")
    require(
        file_sha256(THM2911_VERIFY) == THM2911_VERIFY_SHA256,
        "THM2911 verifier changed",
    )
    require(file_sha256(THM2911_OUT) == THM2911_OUT_SHA256, "THM2911 output changed")
    C.check_dependency_outputs()
    require(
        next(
            line.removeprefix("mode=")
            for line in BASE_OUT.read_text().splitlines()
            if line.startswith("mode=")
        )
        == "LOCKED"
        and next(
            line.removeprefix("all_exact_controls=")
            for line in BASE_OUT.read_text().splitlines()
            if line.startswith("all_exact_controls=")
        )
        == "PASS",
        "THM2912 replay controls changed",
    )
    thm2911_text = THM2911_OUT.read_text()
    require(
        "mode=LOCKED_FINITE_EXACT" in thm2911_text
        and "all_exact_controls=PASS" in thm2911_text
        and "proved_union=181,residual=3251" in thm2911_text,
        "THM2911 lock controls changed",
    )
    top4_closed = set(literal_output_value(BASE_OUT, "closed_roots="))
    top4_failed = set(literal_output_value(BASE_OUT, "failed_roots="))
    require(
        len(top4_closed) == 172
        and len(top4_failed) == 38
        and not (top4_closed & top4_failed),
        "THM2912 root partition changed",
    )
    return top4_closed, top4_failed


def thm2911_route_roots() -> set[tuple[int, ...]]:
    roots = {
        tuple(map(int, line.split("=", 1)[1].split(",")))
        for line in THM2911_OUT.read_text().splitlines()
        if line.startswith("THREE_ROUTE_ROOT=")
    }
    require(len(roots) == 138, "THM2911 route-root list changed")
    return roots


def compute_parent(row: dict[str, object]) -> dict[str, object]:
    return C.compute_parent(row)


def compute_child(
    item: tuple[dict[str, object], int],
) -> dict[str, object]:
    return C.child_task(item)


def subtract_multi(
    carrier: list[tuple[F, F]],
    labels: tuple[int, ...],
) -> list[tuple[F, F]]:
    residual = carrier
    for label in labels:
        residual = C.H.R.subtract_local(residual, label)
    return residual


def direct_carrier(family: tuple[int, ...]) -> tuple[list[tuple[F, F]], int, F]:
    require(len(set(family)) == len(family), f"repeated label in family {family}")
    return C.H.R.CORE.good_norm(tuple(sorted(family)))


def exact_ranked_head(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
    horizon: int,
    seal_rank: int,
) -> tuple[list[tuple[F, int]], F, int]:
    labels = [
        label
        for label in range(FIRST_EXTERNAL, horizon + 1)
        if label not in forbidden
    ]
    require(len(labels) >= seal_rank, "finite head has too few allowed labels")
    rows = C.H.exact_coverages(carrier, labels)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    mass = interval_mass(carrier)
    tail = mass / 7 + C.H.S2 * len(carrier) / (7 * (horizon + 1))
    require(
        ranked[seal_rank - 1][0] >= tail,
        f"rank {seal_rank} was not sealed at horizon {horizon}",
    )
    return ranked, tail, len(labels)


def exact_pair_cap(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
) -> dict[str, object]:
    ranked, singleton_tail, scanned = exact_ranked_head(
        carrier,
        forbidden,
        PAIR_HORIZON,
        4,
    )
    mass = interval_mass(carrier)
    q1 = ranked[0][0]
    best = F(-1)
    witness: tuple[int, int] | None = None
    paid = 0
    paid_digest = hashlib.sha256(b"LRC14/j6/THM2913/paid-pairs/v1\n")
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= best:
            break
        after_first = C.H.R.subtract_local(carrier, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= best:
                break
            residual = C.H.R.subtract_local(after_first, second)
            union = mass - interval_mass(residual)
            paid += 1
            pair = tuple(sorted((first, second)))
            paid_digest.update(
                f"P={pair[0]},{pair[1]};U={ftext(union)}\n".encode()
            )
            if (union, tuple(-label for label in pair)) > (
                best,
                tuple(-label for label in witness) if witness is not None else (),
            ):
                best = union
                witness = pair
    require(paid > 0 and witness is not None, "empty exact pair search")
    pair_tail = q1 + singleton_tail
    require(
        best > pair_tail,
        "finite pair head did not strictly beat the infinite-tail cap",
    )
    direct = mass - interval_mass(subtract_multi(carrier, witness))
    require(direct == best, "pair winner failed independent residual check")
    return {
        "ranked": ranked,
        "top4": tuple(ranked[:4]),
        "singleton_tail": singleton_tail,
        "scanned": scanned,
        "head": best,
        "tail": pair_tail,
        "tail_gap": best - pair_tail,
        "witness": witness,
        "paid": paid,
        "paid_digest": paid_digest.hexdigest(),
    }


def clipped(value: F, upper: F) -> F:
    return min(max(value, F(0)), upper)


def hunter_data(
    qs: tuple[F, ...],
    pair_cap: F,
    mass: F,
) -> tuple[F, F | None]:
    require(2 <= len(qs) <= 4, "Hunter arity changed")
    require(
        all(qs[index] >= qs[index + 1] for index in range(len(qs) - 1)),
        "singleton ranks are not decreasing",
    )
    require(pair_cap >= qs[0], "pair cap fell below singleton cap")
    upper = min(qs[0], pair_cap)
    breakpoints = {F(0), upper, clipped(pair_cap / 2, upper)}
    for singleton in qs[1:]:
        breakpoints.add(clipped(singleton, upper))
        breakpoints.add(clipped(pair_cap - singleton, upper))
    ordered = sorted(breakpoints)

    def objective(center: F) -> F:
        return center + sum(
            (
                min(center, singleton, pair_cap - center)
                for singleton in qs[1:]
            ),
            F(0),
        )

    envelope = max(objective(point) for point in ordered)
    require(envelope <= 2 * pair_cap, "Hunter envelope exceeded pair partition")
    for left, right in zip(ordered, ordered[1:]):
        require(
            2 * objective((left + right) / 2)
            == objective(left) + objective(right),
            "Hunter breakpoint list is incomplete",
        )
    if envelope < mass:
        return envelope, None
    for left, right in zip(ordered, ordered[1:]):
        left_value = objective(left)
        right_value = objective(right)
        if left_value >= mass:
            return envelope, left
        if left_value < mass <= right_value:
            threshold = (
                left
                + (mass - left_value)
                * (right - left)
                / (right_value - left_value)
            )
            require(
                objective(threshold) == mass,
                "hostile second-centre threshold solve failed",
            )
            return envelope, threshold
    require(objective(ordered[-1]) >= mass, "hostile level disappeared")
    return envelope, ordered[-1]


def second_pivot(
    item: tuple[
        dict[str, object],
        list[tuple[F, F]],
        frozenset[int],
        tuple[tuple[F, int], ...],
        int,
    ],
) -> dict[str, object]:
    child_row, child, forbidden, core, index = item
    coverage, center = core[index]
    earlier = tuple(label for _, label in core[:index])
    second_forbidden = frozenset((*forbidden, center, *earlier))
    grandchild = C.H.R.subtract_local(child, center)
    grandchild_mass = interval_mass(grandchild)
    require(grandchild_mass > 0, "second centre covered the child")
    direct, components, direct_mass = direct_carrier(
        tuple(
            (
                *child_row["body"],
                child_row["apex"],
                child_row["center"],
                center,
            )
        )
    )
    require(
        grandchild == direct
        and len(grandchild) == components
        and grandchild_mass == direct_mass,
        "literal/direct grandchild disagreement",
    )
    ranked, tail, scanned = exact_ranked_head(
        grandchild,
        second_forbidden,
        TOP3_HORIZON,
        3,
    )
    top3 = tuple(ranked[:3])
    margin = grandchild_mass - sum((value for value, _ in top3), F(0))
    return {
        "center": center,
        "coverage": coverage,
        "earlier": earlier,
        "mass": grandchild_mass,
        "components": len(grandchild),
        "scanned": scanned,
        "tail": tail,
        "tail_gap": top3[-1][0] - tail,
        "top3": top3,
        "margin": margin,
        "closed": margin > 0,
    }


def analyze_failed_child(child_row: dict[str, object]) -> dict[str, object]:
    parent, _, _ = direct_carrier(
        tuple((*child_row["body"], child_row["apex"]))
    )
    child = C.H.R.subtract_local(parent, child_row["center"])
    mass = interval_mass(child)
    require(
        mass == child_row["child_mass"]
        and len(child) == child_row["child_components"],
        "THM2912 child reconstruction changed",
    )
    direct, components, direct_mass = direct_carrier(
        tuple((*child_row["body"], child_row["apex"], child_row["center"]))
    )
    require(
        child == direct and len(child) == components and mass == direct_mass,
        "literal/direct child disagreement",
    )
    forbidden = frozenset(
        (*child_row["prefix"], child_row["center"], *child_row["earlier"])
    )
    pair = exact_pair_cap(child, forbidden)
    require(
        pair["top4"] == child_row["top_four"],
        "2500-head and THM2912 top-four ranks disagree",
    )
    pair_margin = mass - 2 * pair["head"]
    qs = tuple(value for value, _ in pair["top4"])
    envelope, threshold = hunter_data(qs, pair["head"], mass)
    hunter_margin = mass - envelope
    second_core: tuple[tuple[F, int], ...] = ()
    second_rows: tuple[dict[str, object], ...] = ()
    cutoff = None
    delta = None
    core_scan = 0
    if threshold is not None:
        delta = threshold - mass / 7
        require(delta > 0, "second-centre threshold is not discrepancy-finite")
        gamma = C.H.S2 * len(child) / 7
        cutoff = C.H.ceiling(gamma / delta) - 1
        require(cutoff >= FIRST_EXTERNAL, "second-centre cutoff is empty")
        labels = [
            label
            for label in range(FIRST_EXTERNAL, cutoff + 1)
            if label not in forbidden
        ]
        coverages = C.H.exact_coverages(child, labels)
        core_scan = len(labels)
        second_core = tuple(
            sorted(
                (
                    (value, label)
                    for value, label in coverages
                    if value >= threshold
                ),
                key=lambda item: (-item[0], item[1]),
            )
        )
        require(second_core, "hostile second-centre core is empty")
        require(
            mass / 7 + gamma / (cutoff + 1) <= threshold,
            "second-centre discrepancy tail did not seal",
        )
        second_rows = tuple(
            second_pivot((child_row, child, forbidden, second_core, index))
            for index in range(len(second_core))
        )
        require(
            all(row["closed"] for row in second_rows),
            "open second-centre pivot survived",
        )
    if pair_margin > 0:
        route = "pair"
    elif hunter_margin > 0:
        route = "hunter"
    else:
        require(threshold is not None and second_rows, "recursive route vanished")
        route = "toothpick"
    return {
        **child_row,
        "pair": pair,
        "pair_margin": pair_margin,
        "hunter": envelope,
        "hunter_margin": hunter_margin,
        "threshold": threshold,
        "delta": delta,
        "cutoff": cutoff,
        "core_scan": core_scan,
        "second_core": second_core,
        "second_rows": second_rows,
        "route": route,
        "closed": True,
    }


def child_line(row: dict[str, object]) -> str:
    pair = row["pair"]
    return (
        f"CHILD;E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"x={row['center']};earlier={','.join(map(str, row['earlier']))};"
        f"h={ftext(row['child_mass'])};r={row['child_components']};"
        f"scan={pair['scanned']};Stail={ftext(pair['singleton_tail'])};"
        + "top4="
        + ",".join(f"{label}:{ftext(value)}" for value, label in pair["top4"])
        + f";B2={ftext(pair['head'])};B2tail={ftext(pair['tail'])};"
        f"B2gap={ftext(pair['tail_gap'])};B2w={pair['witness']};"
        f"paid={pair['paid']};paiddigest={pair['paid_digest']};"
        f"pairmargin={ftext(row['pair_margin'])};"
        f"G4={ftext(row['hunter'])};G4margin={ftext(row['hunter_margin'])};"
        f"lambda={None if row['threshold'] is None else ftext(row['threshold'])};"
        f"delta={None if row['delta'] is None else ftext(row['delta'])};"
        f"N={row['cutoff']};corescan={row['core_scan']};"
        + "H1="
        + ",".join(
            f"{label}:{ftext(value)}" for value, label in row["second_core"]
        )
        + f";route={row['route']}\n"
    )


def second_line(parent: dict[str, object], row: dict[str, object]) -> str:
    return (
        f"SECOND;E={','.join(map(str, parent['body']))};rank={parent['rank']};"
        f"a={parent['apex']};x={parent['center']};y={row['center']};"
        f"earlier={','.join(map(str, row['earlier']))};"
        f"h={ftext(row['mass'])};r={row['components']};scan={row['scanned']};"
        f"tail={ftext(row['tail'])};tailgap={ftext(row['tail_gap'])};"
        + "top3="
        + ",".join(f"{label}:{ftext(value)}" for value, label in row["top3"])
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}\n"
    )


def boundary_text(row: dict[str, object], margin_key: str) -> str:
    return (
        f"{row['body']};rank={row['rank']};a={row['apex']};"
        f"x={row['center']};margin={ftext(row[margin_key])}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    top4_closed_artifact, top4_failed_artifact = dependency_controls()

    inputs = C.H.survivor_inputs()
    context = mp.get_context("spawn")
    if args.workers == 1:
        parent_rows = list(map(compute_parent, inputs))
    else:
        with context.Pool(args.workers) as pool:
            parent_rows = pool.map(compute_parent, inputs)
    parent_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    one_rows, _ = C.residual_one_rows(parent_rows)
    tasks = [
        (row, index)
        for row in one_rows
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]
    require(len(tasks) == 807, "ordered child universe changed")
    if args.workers == 1:
        children = list(map(compute_child, tasks))
    else:
        with context.Pool(args.workers) as pool:
            children = pool.map(compute_child, tasks)
    children.sort(
        key=lambda row: (row["body"], row["rank"], row["apex"], row["center"])
    )
    base_closed = [row for row in children if row["closed"]]
    base_failed = [row for row in children if not row["closed"]]
    require(
        len(base_closed) == 765 and len(base_failed) == 42,
        "THM2912 child partition changed",
    )

    if args.workers == 1:
        repaired = list(map(analyze_failed_child, base_failed))
    else:
        with context.Pool(args.workers) as pool:
            repaired = pool.map(analyze_failed_child, base_failed)
    repaired.sort(
        key=lambda row: (row["body"], row["rank"], row["apex"], row["center"])
    )
    failed_roots = tuple(sorted({row["body"] for row in repaired}))
    children_by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in children:
        children_by_body[row["body"]].append(row)
    computed_top4_closed = {
        body
        for body, body_rows in children_by_body.items()
        if all(row["closed"] for row in body_rows)
    }
    require(set(failed_roots) == top4_failed_artifact, "failed-root list changed")
    require(
        computed_top4_closed == top4_closed_artifact,
        "closed-root list changed",
    )

    prior88 = C.current_proved_union()
    route2911 = thm2911_route_roots()
    baseline181 = prior88 | route2911
    top4_union = baseline181 | top4_closed_artifact
    additive_roots = tuple(sorted(set(failed_roots) - top4_union))
    full_one_row = top4_closed_artifact | set(failed_roots)
    final_union = baseline181 | full_one_row
    require(
        len(prior88) == 88
        and len(baseline181) == 181
        and len(top4_closed_artifact & baseline181) == 39
        and len(top4_union) == 314
        and len(set(failed_roots) & baseline181) == 1
        and set(failed_roots) & baseline181 == {(1, 2, 3, 5, 6, 8, 13)}
        and not (set(failed_roots) & top4_closed_artifact)
        and len(additive_roots) == 37
        and len(full_one_row) == 210
        and len(full_one_row & baseline181) == 40
        and len(final_union) == 351
        and 3_432 - len(final_union) == 3_081,
        "live root recomposition changed",
    )

    second_rows = [
        (parent, row)
        for parent in repaired
        for row in parent["second_rows"]
    ]
    counts = (
        len(parent_rows),
        11_511,
        3_344,
        len(one_rows),
        len(children),
        len(base_closed),
        len(repaired),
        sum(row["route"] == "pair" for row in repaired),
        sum(row["pair_margin"] <= 0 for row in repaired),
        sum(row["route"] == "hunter" for row in repaired),
        sum(row["route"] == "toothpick" for row in repaired),
        len(second_rows),
        sum(row["closed"] for _, row in second_rows),
        len(failed_roots),
        len(set(failed_roots) & baseline181),
        len(additive_roots),
        len(baseline181),
        len(top4_union),
        len(final_union),
        3_432 - len(final_union),
        sum(row["pair"]["scanned"] for row in repaired),
        sum(row["pair"]["paid"] for row in repaired),
        sum(row["core_scan"] for row in repaired),
        sum(row["scanned"] for _, row in second_rows),
        min(len(row["second_core"]) for row in repaired if row["second_core"]),
        max(len(row["second_core"]) for row in repaired if row["second_core"]),
        sum(row["pair_margin"] == 0 for row in repaired),
        sum(row["hunter_margin"] == 0 for row in repaired),
        sum(row["margin"] == 0 for _, row in second_rows),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "THM2913 counts changed")
    failed_root_digest = hashlib.sha256(repr(failed_roots).encode()).hexdigest()
    additive_root_digest = hashlib.sha256(repr(additive_roots).encode()).hexdigest()
    if EXPECTED_FAILED_ROOT_DIGEST is not None:
        require(
            failed_root_digest == EXPECTED_FAILED_ROOT_DIGEST,
            "failed roots changed",
        )
    if EXPECTED_ADDITIVE_ROOT_DIGEST is not None:
        require(
            additive_root_digest == EXPECTED_ADDITIVE_ROOT_DIGEST,
            "additive roots changed",
        )

    ledger_lines: list[str] = []
    for row in repaired:
        ledger_lines.append(child_line(row))
        ledger_lines.extend(second_line(row, child) for child in row["second_rows"])
    ledger_hash = hashlib.sha256()
    ledger_hash.update(b"LRC14/j6/one-H3-row/pair-Hunter-toothpick/v1\n")
    for line in ledger_lines:
        ledger_hash.update(line.encode())
    ledger_sha256 = ledger_hash.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "THM2913 ledger changed")
    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 one-H3-row pair-Hunter toothpick ledger\n"
            + "".join(ledger_lines)
            + f"ledger_sha256={ledger_sha256}\n"
            + "scope=42 THM2912-open ordered children;38 roots;not LRC14\n"
        )

    pair_positive = [row for row in repaired if row["pair_margin"] > 0]
    pair_failed = [row for row in repaired if row["pair_margin"] <= 0]
    hunter_positive = [
        row
        for row in pair_failed
        if row["hunter_margin"] > 0
    ]
    hunter_failed = [
        row
        for row in pair_failed
        if row["hunter_margin"] <= 0
    ]
    closest_pair_positive = min(
        pair_positive,
        key=lambda row: (
            row["pair_margin"],
            row["body"],
            row["center"],
        ),
    )
    closest_pair_failure = max(
        pair_failed,
        key=lambda row: (
            row["pair_margin"],
            tuple(-value for value in row["body"]),
            -row["center"],
        ),
    )
    closest_hunter_positive = min(
        hunter_positive,
        key=lambda row: (
            row["hunter_margin"],
            row["body"],
            row["center"],
        ),
    )
    closest_hunter_failure = max(
        hunter_failed,
        key=lambda row: (
            row["hunter_margin"],
            tuple(-value for value in row["body"]),
            -row["center"],
        ),
    )
    closest_second = min(
        second_rows,
        key=lambda item: (
            item[1]["margin"],
            item[0]["body"],
            item[0]["center"],
            item[1]["center"],
        ),
    )

    lines = [
        "LRC14 one-H3-row pair-Hunter toothpick closure",
        f"pair_horizon={PAIR_HORIZON};top3_horizon={TOP3_HORIZON}",
        f"counts={counts}",
        f"closed_roots={failed_roots}",
        f"additive_roots={additive_roots}",
        f"closed_root_digest={failed_root_digest}",
        f"additive_root_digest={additive_root_digest}",
        "overlap_THM2911=((1, 2, 3, 5, 6, 8, 13),)",
        f"closest_pair_positive={boundary_text(closest_pair_positive, 'pair_margin')}",
        f"closest_pair_failure={boundary_text(closest_pair_failure, 'pair_margin')}",
        f"closest_hunter_positive={boundary_text(closest_hunter_positive, 'hunter_margin')}",
        f"closest_hunter_failure={boundary_text(closest_hunter_failure, 'hunter_margin')}",
        (
            "closest_second_positive="
            f"{closest_second[0]['body']};rank={closest_second[0]['rank']};"
            f"a={closest_second[0]['apex']};x={closest_second[0]['center']};"
            f"y={closest_second[1]['center']};"
            f"margin={ftext(closest_second[1]['margin'])}"
        ),
        (
            "minimum_pair_tail_gap="
            f"{ftext(min(row['pair']['tail_gap'] for row in repaired))}"
        ),
        (
            "minimum_second_tail_gap="
            f"{ftext(min(row['tail_gap'] for _, row in second_rows))}"
        ),
        (
            "minimum_hostile_delta="
            f"{ftext(min(row['delta'] for row in hunter_failed))}"
        ),
        f"second_core_sizes={tuple(len(row['second_core']) for row in hunter_failed)}",
        f"ledger_sha256={ledger_sha256}",
        (
            "mode=DISCOVERY"
            if any(
                value is None
                for value in (
                    EXPECTED_COUNTS,
                    EXPECTED_FAILED_ROOT_DIGEST,
                    EXPECTED_ADDITIVE_ROOT_DIGEST,
                    EXPECTED_LEDGER_SHA256,
                )
            )
            else "mode=LOCKED"
        ),
        (
            "scope=42 ordered children and 38 roots left open by THM2912;"
            "37 additive roots over the live THM2911+THM2912 union;not LRC14"
        ),
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
