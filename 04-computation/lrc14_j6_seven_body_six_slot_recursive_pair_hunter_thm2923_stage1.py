#!/usr/bin/env python3
"""Exact first/second-level repair of the 4,866 THM-2915 failures.

This proof-bearing stage keeps the full THM-2904 branch identity

    (body, rank, apex, prefix, first centre, earlier first centres)

and reconstructs every child both sequentially and directly.  All singleton,
pair, and recursive top-three heads are sealed against the global THM-735
tail; no finite search cap is treated as a theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_all_open_centre_child_top4_closure_thm2915.ledger.out"
)
ENGINE_PATH = (
    ROOT
    / "04-computation/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.py"
)
FIRST_EXTERNAL = 15
INITIAL_PAIR_HORIZON = 2_000
INITIAL_TOP3_HORIZON = 2_000
EXPECTED_ENGINE_SHA256 = (
    "14e56e124197cd1bdae841efa195a1e7c282e7ea368a610e5f4d56509431858b"
)
EXPECTED_INPUT_SHA256 = (
    "798cd660ab60e2021b28074a1390af3f6b1367c99f2d0ab63a581513f7871071"
)
EXPECTED_PAIR_DIGEST = (
    "2b0580c8e98e13aeb34adc2f4ec3f8da5be30363adb882cd6794da192ce9476c"
)
EXPECTED_RECURSIVE_DIGEST = (
    "c9e4a568f17778355b4ce91ee387add8f67dd6dcca257108ec945a2dfe1bd327"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_engine():
    require(file_sha256(ENGINE_PATH) == EXPECTED_ENGINE_SHA256, "engine changed")
    spec = importlib.util.spec_from_file_location("generalized_thm2913_engine", ENGINE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-2913 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_engine()
H = T.C.H


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def parse_ranked(text: str) -> tuple[tuple[F, int], ...]:
    return tuple(
        (F(value), int(label))
        for label, value in (item.split(":", 1) for item in text.split(","))
    )


def parse_failures(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    all_keys: set[tuple[object, ...]] = set()
    open_parent_keys: set[tuple[object, ...]] = set()
    closed_count = 0
    for line in path.read_text().splitlines():
        if not line.startswith("CHILD;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (
            ints(fields["E"]),
            int(fields["rank"]),
            int(fields["a"]),
            ints(fields["P"]),
            int(fields["x"]),
            ints(fields["earlier"]),
        )
        require(key not in all_keys, f"duplicate THM-2915 child identity {key}")
        all_keys.add(key)
        open_parent_keys.add(key[:4])
        margin = F(fields["margin"])
        top4 = parse_ranked(fields["top4"])
        derived_closed = margin > 0
        require(
            (fields["closed"] == "1") == derived_closed,
            f"{key}: THM-2915 closed bit disagrees with exact margin",
        )
        require(
            F(fields["hR"]) - sum((value for value, _ in top4), F(0)) == margin,
            f"{key}: THM-2915 top-four margin arithmetic changed",
        )
        closed_count += derived_closed
        if derived_closed:
            continue
        row: dict[str, object] = {
            "body": ints(fields["E"]),
            "rank": int(fields["rank"]),
            "apex": int(fields["a"]),
            "prefix": ints(fields["P"]),
            "center": int(fields["x"]),
            "earlier": ints(fields["earlier"]),
            "child_mass": F(fields["hR"]),
            "child_components": int(fields["rR"]),
            "ledger_horizon": int(fields["M"]),
            "ledger_top4": top4,
            "ledger_margin": margin,
        }
        rows.append(row)
    rows.sort(key=identity)
    require(
        len(all_keys) == 51_222
        and len(open_parent_keys) == 11_563
        and closed_count == 46_356,
        "THM-2915 all-child/open-parent identity universe changed",
    )
    require(len(rows) == 4_866, f"expected 4,866 failures, found {len(rows)}")
    require(len({identity(row) for row in rows}) == len(rows), "child identity collision")
    return rows


def parent_key(row: dict[str, object]) -> tuple[object, ...]:
    return row["body"], row["rank"], row["apex"], row["prefix"]


def identity(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["rank"],
        row["apex"],
        row["prefix"],
        row["center"],
        row["earlier"],
    )


def interval_mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def subtract_multi(
    carrier: list[tuple[F, F]],
    labels: Iterable[int],
) -> list[tuple[F, F]]:
    residual = carrier
    for label in labels:
        residual = H.R.subtract_local(residual, label)
    return residual


def direct_carrier(family: tuple[int, ...]) -> tuple[list[tuple[F, F]], int, F]:
    require(len(set(family)) == len(family), f"repeated label in {family}")
    return H.R.CORE.good_norm(tuple(sorted(family)))


def reconstruct_child(
    row: dict[str, object],
) -> tuple[list[tuple[F, F]], frozenset[int]]:
    parent, _, _ = direct_carrier(tuple((*row["body"], row["apex"])))
    child = H.R.subtract_local(parent, row["center"])
    child_family = tuple((*row["body"], row["apex"], row["center"]))
    require(len(child_family) == 9, "child arity is not 9 fixed + 4 remaining")
    direct, components, mass = direct_carrier(child_family)
    require(
        child == direct
        and len(child) == components == row["child_components"]
        and interval_mass(child) == mass == row["child_mass"],
        f"{identity(row)}: direct/sequential child mismatch",
    )
    forbidden = frozenset((*row["prefix"], row["center"], *row["earlier"]))
    require(
        len(forbidden) == len(row["prefix"]) + 1 + len(row["earlier"])
        and row["center"] not in row["prefix"]
        and row["center"] not in row["earlier"]
        and not (set(row["prefix"]) & set(row["earlier"])),
        f"{identity(row)}: inherited forbidden-prefix collision",
    )
    return child, forbidden


def ranked_head(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
    horizon: int,
) -> tuple[list[tuple[F, int]], int]:
    labels = [
        label
        for label in range(FIRST_EXTERNAL, horizon + 1)
        if label not in forbidden
    ]
    ranked = sorted(
        H.exact_coverages(carrier, labels),
        key=lambda item: (-item[0], item[1]),
    )
    return ranked, len(labels)


def pair_head(
    carrier: list[tuple[F, F]],
    ranked: list[tuple[F, int]],
) -> tuple[F, tuple[int, int], int]:
    mass = interval_mass(carrier)
    first, second = ranked[0][1], ranked[1][1]
    witness = tuple(sorted((first, second)))
    best = mass - interval_mass(subtract_multi(carrier, witness))
    paid = 1
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] < best:
            break
        after_first = H.R.subtract_local(carrier, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value < best:
                break
            candidate_pair = tuple(sorted((first, second)))
            union = mass - interval_mass(
                H.R.subtract_local(after_first, second)
            )
            paid += 1
            if (union, tuple(-label for label in candidate_pair)) > (
                best,
                tuple(-label for label in witness),
            ):
                best = union
                witness = candidate_pair
    require(
        mass - interval_mass(subtract_multi(carrier, witness)) == best,
        "pair witness failed direct residual reconstruction",
    )
    return best, witness, paid


def exact_pair_cap(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
    initial_horizon: int,
) -> dict[str, object]:
    mass = interval_mass(carrier)
    gamma = H.S2 * len(carrier) / 7
    horizon = max(INITIAL_PAIR_HORIZON, initial_horizon)
    extensions = 0
    while True:
        ranked, scanned = ranked_head(carrier, forbidden, horizon)
        require(len(ranked) >= 4, "pair head has fewer than four labels")
        best, witness, paid = pair_head(carrier, ranked)
        singleton_tail = mass / 7 + gamma / (horizon + 1)
        pair_tail = ranked[0][0] + singleton_tail
        if best >= pair_tail:
            break
        delta = best - ranked[0][0] - mass / 7
        require(
            delta > 0,
            "pair cap cannot be sealed by one-head plus discrepancy tail",
        )
        required = ceiling(gamma / delta) - 1
        require(required > horizon, "pair horizon extension failed to advance")
        horizon = required
        extensions += 1
        require(extensions <= 4, "pair horizon failed to stabilize")
    top4 = tuple(ranked[:4])
    return {
        "top4": top4,
        "horizon": horizon,
        "scanned": scanned,
        "head": best,
        "witness": witness,
        "paid": paid,
        "singleton_tail": singleton_tail,
        "tail": pair_tail,
        "tail_gap": best - pair_tail,
        "extensions": extensions,
    }


def analyze_pair(row: dict[str, object]) -> dict[str, object]:
    child, forbidden = reconstruct_child(row)
    pair = exact_pair_cap(child, forbidden, int(row["ledger_horizon"]))
    pair_residual = subtract_multi(child, pair["witness"])
    pair_family = tuple(
        (
            *row["body"],
            row["apex"],
            row["center"],
            *pair["witness"],
        )
    )
    require(len(pair_family) == 11, "child B2 did not consume exactly two slots")
    pair_direct, pair_components, pair_mass = direct_carrier(pair_family)
    require(
        pair_residual == pair_direct
        and len(pair_residual) == pair_components
        and interval_mass(pair_residual) == pair_mass
        and row["child_mass"] - pair_mass == pair["head"],
        f"{identity(row)}: child pair winner failed full-family reconstruction",
    )
    ledger_top4 = row["ledger_top4"]
    require(
        tuple(value for value, _ in pair["top4"])
        == tuple(value for value, _ in ledger_top4),
        f"{identity(row)}: sealed top-four values disagree with THM-2915 ledger",
    )
    require(
        row["child_mass"] - sum((value for value, _ in pair["top4"]), F(0))
        == row["ledger_margin"],
        f"{identity(row)}: child top-four margin changed",
    )
    qs = tuple(value for value, _ in pair["top4"])
    envelope, threshold = T.hunter_data(qs, pair["head"], row["child_mass"])
    pair_margin = row["child_mass"] - 2 * pair["head"]
    hunter_margin = row["child_mass"] - envelope
    route = (
        "pair"
        if pair_margin > 0
        else "hunter"
        if hunter_margin > 0
        else "recursive"
    )
    return {
        **row,
        "pair": pair,
        "hunter": envelope,
        "threshold": threshold,
        "pair_margin": pair_margin,
        "hunter_margin": hunter_margin,
        "route": route,
    }


def exact_top_k(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
    rank: int,
    initial_horizon: int,
) -> dict[str, object]:
    mass = interval_mass(carrier)
    gamma = H.S2 * len(carrier) / 7
    horizon = initial_horizon
    extensions = 0
    while True:
        ranked, scanned = ranked_head(carrier, forbidden, horizon)
        require(len(ranked) >= rank, "top-k head too short")
        kth = ranked[rank - 1][0]
        tail = mass / 7 + gamma / (horizon + 1)
        if kth >= tail:
            break
        delta = kth - mass / 7
        require(delta > 0, "top-k rank has no discrepancy-finite excess")
        required = ceiling(gamma / delta) - 1
        require(required > horizon, "top-k extension failed to advance")
        horizon = required
        extensions += 1
        require(extensions <= 4, "top-k horizon failed to stabilize")
    return {
        "top": tuple(ranked[:rank]),
        "horizon": horizon,
        "scanned": scanned,
        "tail": tail,
        "tail_gap": kth - tail,
        "extensions": extensions,
    }


def second_pivot(
    row: dict[str, object],
    child: list[tuple[F, F]],
    forbidden: frozenset[int],
    core: tuple[tuple[F, int], ...],
    index: int,
) -> dict[str, object]:
    coverage, center = core[index]
    earlier = tuple(label for _, label in core[:index])
    second_forbidden = frozenset((*forbidden, center, *earlier))
    require(
        center not in forbidden
        and not (set(earlier) & forbidden)
        and center not in earlier
        and len(second_forbidden) == len(forbidden) + 1 + len(earlier),
        f"{identity(row)};y={center}: second-centre forbidden sidecar collided",
    )
    grandchild = H.R.subtract_local(child, center)
    grandchild_family = tuple(
        (*row["body"], row["apex"], row["center"], center)
    )
    require(
        len(grandchild_family) == 10,
        "second pivot is not 10 fixed + 3 remaining",
    )
    direct, components, mass = direct_carrier(grandchild_family)
    require(
        grandchild == direct
        and len(grandchild) == components
        and interval_mass(grandchild) == mass,
        f"{identity(row)};y={center}: direct/sequential grandchild mismatch",
    )
    top3 = exact_top_k(
        grandchild,
        second_forbidden,
        rank=3,
        initial_horizon=INITIAL_TOP3_HORIZON,
    )
    margin = mass - sum((value for value, _ in top3["top"]), F(0))
    return {
        "center": center,
        "coverage": coverage,
        "earlier": earlier,
        "mass": mass,
        "components": components,
        "top3": top3,
        "margin": margin,
        "closed": margin > 0,
    }


def analyze_recursive(row: dict[str, object]) -> dict[str, object]:
    require(row["route"] == "recursive", "nonrecursive row entered toothpick pass")
    child, forbidden = reconstruct_child(row)
    threshold = row["threshold"]
    require(threshold is not None, "Hunter failure has no hostile threshold")
    mass = row["child_mass"]
    delta = threshold - mass / 7
    require(delta > 0, "recursive threshold is not discrepancy-finite")
    gamma = H.S2 * len(child) / 7
    cutoff = ceiling(gamma / delta) - 1
    require(cutoff >= FIRST_EXTERNAL, "empty second-centre cutoff")
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in forbidden
    ]
    coverages = H.exact_coverages(child, labels)
    core = tuple(
        sorted(
            (
                (value, label)
                for value, label in coverages
                if value >= threshold
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    require(core, "recursive hostile core is empty")
    require(
        mass / 7 + gamma / (cutoff + 1) <= threshold,
        "recursive hostile core tail did not seal",
    )
    pivots = tuple(
        second_pivot(row, child, forbidden, core, index)
        for index in range(len(core))
    )
    return {
        **row,
        "delta": delta,
        "cutoff": cutoff,
        "core": core,
        "pivots": pivots,
        "recursive_closed": all(pivot["closed"] for pivot in pivots),
    }


def key_text(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"x={row['center']};earlier={','.join(map(str, row['earlier']))}"
    )


def pair_line(row: dict[str, object]) -> str:
    pair = row["pair"]
    return (
        "PAIR;"
        + key_text(row)
        + f";h={ftext(row['child_mass'])};r={row['child_components']};"
        + f"M={pair['horizon']};scan={pair['scanned']};"
        + "top4="
        + ",".join(f"{label}:{ftext(value)}" for value, label in pair["top4"])
        + f";B2={ftext(pair['head'])};B2w={pair['witness'][0]},{pair['witness'][1]};"
        + f"B2tail={ftext(pair['tail'])};B2gap={ftext(pair['tail_gap'])};"
        + f"paid={pair['paid']};extensions={pair['extensions']};"
        + f"pairmargin={ftext(row['pair_margin'])};"
        + f"G4={ftext(row['hunter'])};G4margin={ftext(row['hunter_margin'])};"
        + f"lambda={None if row['threshold'] is None else ftext(row['threshold'])};"
        + f"route={row['route']}"
    )


def recursive_lines(row: dict[str, object]) -> list[str]:
    lines = [
        (
            "RECURSIVE;"
            + key_text(row)
            + f";lambda={ftext(row['threshold'])};delta={ftext(row['delta'])};"
            + f"N={row['cutoff']};"
            + "H1="
            + ",".join(f"{label}:{ftext(value)}" for value, label in row["core"])
            + f";closed={int(row['recursive_closed'])}"
        )
    ]
    for pivot in row["pivots"]:
        top3 = pivot["top3"]
        lines.append(
            "SECOND;"
            + key_text(row)
            + f";y={pivot['center']};"
            + f"yearlier={','.join(map(str, pivot['earlier']))};"
            + f"h={ftext(pivot['mass'])};r={pivot['components']};"
            + f"M={top3['horizon']};scan={top3['scanned']};"
            + "top3="
            + ",".join(f"{label}:{ftext(value)}" for value, label in top3["top"])
            + f";tail={ftext(top3['tail'])};tailgap={ftext(top3['tail_gap'])};"
            + f"margin={ftext(pivot['margin'])};closed={int(pivot['closed'])}"
        )
    return lines


def digest_lines(domain: bytes, lines: list[str]) -> str:
    digest = hashlib.sha256(domain + b"\n")
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def root_and_parent_composition(
    all_failures: list[dict[str, object]],
    newly_closed: set[tuple[object, ...]],
) -> dict[str, object]:
    failed_by_parent: dict[tuple[object, ...], set[tuple[object, ...]]] = defaultdict(set)
    failed_by_body: dict[tuple[int, ...], set[tuple[object, ...]]] = defaultdict(set)
    for row in all_failures:
        failed_by_parent[parent_key(row)].add(identity(row))
        failed_by_body[row["body"]].add(identity(row))
    completed_parents = {
        key for key, failures in failed_by_parent.items() if failures <= newly_closed
    }
    completed_roots = {
        body for body, failures in failed_by_body.items() if failures <= newly_closed
    }
    root_digest = hashlib.sha256(b"LRC14/THM2915-failures/one-level-roots/v1\n")
    for body in sorted(completed_roots):
        root_digest.update((",".join(map(str, body)) + "\n").encode())
    return {
        "failed_parents": len(failed_by_parent),
        "failed_roots": len(failed_by_body),
        "completed_parents": len(completed_parents),
        "remaining_parents": len(failed_by_parent) - len(completed_parents),
        "completed_roots": len(completed_roots),
        "remaining_roots": len(failed_by_body) - len(completed_roots),
        "completed_root_sha256": root_digest.hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--pair-ledger", type=Path, required=True)
    parser.add_argument("--recursive-ledger", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--skip-recursive", action="store_true")
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")

    failures = parse_failures(args.input)
    require(
        file_sha256(args.input) == EXPECTED_INPUT_SHA256,
        "THM-2915 child ledger changed",
    )
    selected = failures if args.limit is None else failures[: args.limit]
    if args.workers == 1:
        pair_rows = list(map(analyze_pair, selected))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            pair_rows = pool.map(analyze_pair, selected, chunksize=1)
    pair_rows.sort(key=identity)
    pair_lines = [pair_line(row) for row in pair_rows]
    pair_digest = digest_lines(b"LRC14/THM2915-failures/pair-Hunter/v1", pair_lines)
    args.pair_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.pair_ledger.write_text(
        "LRC14 THM2915-failure exact pair/Hunter finite-exact ledger\n"
        + "\n".join(pair_lines)
        + "\n"
        + f"semantic_sha256={pair_digest}\n"
        + f"scope={len(selected)} of 4866 THM2915 child-top4 failures;"
        + "seven-body/six-external rung;not unrestricted LRC14\n",
        encoding="utf-8",
        newline="\n",
    )

    recursive_inputs = [row for row in pair_rows if row["route"] == "recursive"]
    recursive_rows: list[dict[str, object]] = []
    if recursive_inputs and not args.skip_recursive:
        if args.workers == 1:
            recursive_rows = list(map(analyze_recursive, recursive_inputs))
        else:
            with mp.get_context("spawn").Pool(args.workers) as pool:
                recursive_rows = pool.map(analyze_recursive, recursive_inputs, chunksize=1)
        recursive_rows.sort(key=identity)
        lines = [
            line
            for row in recursive_rows
            for line in recursive_lines(row)
        ]
        recursive_digest = digest_lines(
            b"LRC14/THM2915-failures/one-toothpick/v1", lines
        )
        require(args.recursive_ledger is not None, "--recursive-ledger is required")
        args.recursive_ledger.write_text(
            "LRC14 THM2915-failure one-toothpick finite-exact ledger\n"
            + "\n".join(lines)
            + "\n"
            + f"semantic_sha256={recursive_digest}\n"
            + f"scope={len(recursive_rows)} Hunter-star failures;"
            + "seven-body/six-external rung;not unrestricted LRC14\n",
            encoding="utf-8",
            newline="\n",
        )
    else:
        recursive_digest = None

    recursive_by_id = {identity(row): row for row in recursive_rows}
    closed_ids = {
        identity(row)
        for row in pair_rows
        if row["route"] in {"pair", "hunter"}
        or (
            identity(row) in recursive_by_id
            and recursive_by_id[identity(row)]["recursive_closed"]
        )
    }
    pair_positive = [row for row in pair_rows if row["pair_margin"] > 0]
    pair_failed = [row for row in pair_rows if row["pair_margin"] <= 0]
    hunter_positive = [
        row for row in pair_failed if row["hunter_margin"] > 0
    ]
    hunter_failed = [
        row for row in pair_failed if row["hunter_margin"] <= 0
    ]
    second_rows = [
        (row, pivot)
        for row in recursive_rows
        for pivot in row["pivots"]
    ]
    if len(selected) == 4_866 and not args.skip_recursive:
        require(pair_digest == EXPECTED_PAIR_DIGEST, "pair semantic digest changed")
        require(
            recursive_digest == EXPECTED_RECURSIVE_DIGEST,
            "recursive semantic digest changed",
        )
        require(
            (
                len(pair_positive),
                len(pair_failed),
                len(hunter_positive),
                len(hunter_failed),
                len(recursive_rows),
                sum(row["recursive_closed"] for row in recursive_rows),
                sum(not row["recursive_closed"] for row in recursive_rows),
                len(second_rows),
                sum(pivot["closed"] for _, pivot in second_rows),
                sum(not pivot["closed"] for _, pivot in second_rows),
                sum(row["pair"]["paid"] for row in pair_rows),
            )
            == (
                1_612,
                3_254,
                1_175,
                2_079,
                2_079,
                1_884,
                195,
                6_172,
                5_944,
                228,
                40_753,
            ),
            "full pair/Hunter/second-centre census changed",
        )
    summary: dict[str, object] = {
        "selected": len(selected),
        "pair_positive": len(pair_positive),
        "pair_nonpositive": len(pair_failed),
        "pair_equal": sum(row["pair_margin"] == 0 for row in pair_rows),
        "hunter_positive_after_pair": len(hunter_positive),
        "hunter_nonpositive_after_pair": len(hunter_failed),
        "hunter_equal": sum(row["hunter_margin"] == 0 for row in pair_failed),
        "pair_horizon_extensions": sum(row["pair"]["extensions"] for row in pair_rows),
        "maximum_pair_horizon": max(row["pair"]["horizon"] for row in pair_rows),
        "minimum_pair_tail_gap": min(row["pair"]["tail_gap"] for row in pair_rows),
        "paid_pairs": sum(row["pair"]["paid"] for row in pair_rows),
        "recursive_processed": len(recursive_rows),
        "recursive_closed": sum(row["recursive_closed"] for row in recursive_rows),
        "recursive_open": sum(not row["recursive_closed"] for row in recursive_rows),
        "second_pivots": len(second_rows),
        "second_closed": sum(pivot["closed"] for _, pivot in second_rows),
        "second_open": sum(not pivot["closed"] for _, pivot in second_rows),
        "maximum_second_core": max(
            (len(row["core"]) for row in recursive_rows), default=0
        ),
        "maximum_second_cutoff": max(
            (row["cutoff"] for row in recursive_rows), default=0
        ),
        "maximum_top3_horizon": max(
            (pivot["top3"]["horizon"] for _, pivot in second_rows), default=0
        ),
        "minimum_top3_tail_gap": min(
            (pivot["top3"]["tail_gap"] for _, pivot in second_rows), default=None
        ),
        "closed_children": len(closed_ids),
        "open_children": len(selected) - len(closed_ids),
        "route_histogram": tuple(
            sorted(Counter(row["route"] for row in pair_rows).items())
        ),
        "slot_chain": ((9, 4), (10, 3)),
        "composition": root_and_parent_composition(selected, closed_ids),
        "pair_semantic_sha256": pair_digest,
        "recursive_semantic_sha256": recursive_digest,
    }
    if pair_positive:
        summary["closest_pair_positive"] = pair_line(
            min(pair_positive, key=lambda row: (row["pair_margin"], identity(row)))
        )
    if pair_failed:
        summary["closest_pair_failure"] = pair_line(
            max(pair_failed, key=lambda row: (row["pair_margin"], identity(row)))
        )
    if hunter_positive:
        summary["closest_hunter_positive"] = pair_line(
            min(hunter_positive, key=lambda row: (row["hunter_margin"], identity(row)))
        )
    if hunter_failed:
        summary["closest_hunter_failure"] = pair_line(
            max(hunter_failed, key=lambda row: (row["hunter_margin"], identity(row)))
        )
    if second_rows:
        closest_second_positive = [
            pair for pair in second_rows if pair[1]["margin"] > 0
        ]
        closest_second_failure = [
            pair for pair in second_rows if pair[1]["margin"] <= 0
        ]
        if closest_second_positive:
            row, pivot = min(
                closest_second_positive,
                key=lambda pair: (pair[1]["margin"], identity(pair[0]), pair[1]["center"]),
            )
            summary["closest_second_positive"] = (
                key_text(row)
                + f";y={pivot['center']};margin={ftext(pivot['margin'])}"
            )
        if closest_second_failure:
            row, pivot = max(
                closest_second_failure,
                key=lambda pair: (pair[1]["margin"], identity(pair[0]), pair[1]["center"]),
            )
            summary["closest_second_failure"] = (
                key_text(row)
                + f";y={pivot['center']};margin={ftext(pivot['margin'])}"
            )

    rendered = "\n".join(f"{key}={value}" for key, value in summary.items()) + "\n"
    if args.summary is not None:
        args.summary.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
