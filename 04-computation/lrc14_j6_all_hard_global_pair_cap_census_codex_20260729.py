#!/usr/bin/env python3
r"""Globally sealed pair-cap census on all scalar-hard j=6 suffixes.

The input ledger contains the 14,806 marked suffixes left open by the exact
all-root scalar top-five census.  On each literal apex carrier C of mass h,
write c(w)=|C intersect D_w| and q_j for its allowed singleton order
statistics.  Since every input row satisfies q_1<3h/7, THM-2897 gives the
strict one-slot cutoff

    W_2 = floor(gamma / (3h/7-q_1)) + 1,
    gamma = (99/70)r/7.

Every pair containing w>=W_2 has union strictly below

    h/7 + gamma/W_2 + q_1 < 4h/7.

    The program first computes the exact maximum H_2 over pairs entirely
    below W_2.  If

        eta_2 = H_2-q_1-h/7 > 0,
        R_2 = ceil(gamma/eta_2),

    then every pair with an endpoint w>=R_2 has union strictly below H_2.
    Recomputing the finite head through max(W_2,R_2)-1 therefore gives the
    actual global pair supremum beta_2, not merely an upper cap.  If eta_2
    is nonpositive, the safe fallback

    B_2 = max(H_2, h/7 + gamma/W_2 + q_1)

    remains a rigorous global pair cap but is recorded as analytically
    tail-dominated rather than exact.  Every finite H_2 is maximized exactly:
    singleton sums
are only branch-and-bound upper bounds, while every paid pair union is
evaluated from cached literal comb pieces and checked against direct
residual subtraction.

Three consequences are tallied:

* B_2<4h/7 gives THM-2893's finite heavy H_3 core;
* q_5+2B_2<h directly closes the five-slot parent;
* B_2+q_3<5h/7 is a rank-selective global triple cap and gives a finite
  heavy H_2 core containing at least four vertices of a hypothetical cover.
* 2B_2<6h/7 is a global four-set cap and gives a finite heavy H_1 core
  containing all five vertices of a hypothetical cover.

This is a finite exact/scoped census of THM-2896 suffixes, not LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, floor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = Path(
    os.environ.get(
        "LRC14_HARD_LEDGER_PATH",
        ROOT
        / "05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out",
    )
)
LEDGER_SHA256 = (
    "6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37"
)
RESIDUAL_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
RESIDUAL_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)

FIRST_EXTERNAL = 15
S2 = F(99, 70)

# Filled after discovery and then locked for ordinary/optimized replay.
EXPECTED_COUNTS: tuple[int, ...] | None = None
EXPECTED_GROUPS: tuple[tuple[object, ...], ...] | None = None
EXPECTED_EXTREMA: tuple[tuple[object, ...], ...] | None = None
EXPECTED_QUANTILES: tuple[tuple[str, tuple[tuple[int, int], ...]], ...] | None = None
EXPECTED_FAILURE_DIGESTS: tuple[tuple[str, str], ...] | None = None
EXPECTED_LEDGER_DIGEST: str | None = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_residual():
    require(file_sha256(RESIDUAL_PATH) == RESIDUAL_SHA256, "residual source changed")
    spec = importlib.util.spec_from_file_location("j6_global_pair_residual", RESIDUAL_PATH)
    require(spec is not None and spec.loader is not None, "cannot load residual module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_residual()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def parse_csv_ints(text: str) -> tuple[int, ...]:
    return tuple(map(int, text.split(","))) if text else ()


def parse_hard_ledger() -> list[dict[str, object]]:
    require(LEDGER_PATH.exists(), f"missing hard ledger: {LEDGER_PATH}")
    if LEDGER_SHA256 is not None:
        require(file_sha256(LEDGER_PATH) == LEDGER_SHA256, "hard ledger changed")
    rows: list[dict[str, object]] = []
    for raw_line in LEDGER_PATH.read_text().splitlines():
        if not raw_line.startswith("HARD;"):
            continue
        fields: dict[str, str] = {}
        for field in raw_line.removeprefix("HARD;").split(";"):
            key, value = field.split("=", 1)
            fields[key] = value
        top5 = tuple(
            (int(item.split(":", 1)[0]), F(item.split(":", 1)[1]))
            for item in fields["top5"].split(",")
        )
        require(
            len(top5) == 5
            and all(top5[index][1] >= top5[index + 1][1] for index in range(4)),
            "hard-row singleton ranks changed",
        )
        row = {
            "stratum": fields["S"],
            "body": parse_csv_ints(fields["E"]),
            "K": int(fields["K"]),
            "rank": int(fields["rank"]),
            "apex": int(fields["apex"]),
            "prefix": parse_csv_ints(fields["prefix"]),
            "mass": F(fields["m"]),
            "components": int(fields["r"]),
            "scalar_margin": F(fields["margin"]),
            "top5": top5,
        }
        require(
            fields["closed"] == "0"
            and row["scalar_margin"] <= 0
            and len(row["body"]) == 7,
            "input row is not scalar-hard",
        )
        rows.append(row)
    require(len(rows) == 14_806, f"hard-row count changed: {len(rows)}")
    require(
        len(
            {
                (row["body"], row["rank"], row["apex"], row["prefix"])
                for row in rows
            }
        )
        == len(rows),
        "duplicate marked suffix in hard ledger",
    )
    return rows


def allowed_count(last_inclusive: int, forbidden: frozenset[int]) -> int:
    if last_inclusive < FIRST_EXTERNAL:
        return 0
    return (
        last_inclusive - FIRST_EXTERNAL + 1
        - sum(
            FIRST_EXTERNAL <= label <= last_inclusive
            for label in forbidden
        )
    )


def cutoff_data(row: dict[str, object]) -> dict[str, object]:
    mass = row["mass"]
    components = row["components"]
    q1 = row["top5"][0][1]
    entry_margin = 3 * mass / 7 - q1
    require(entry_margin > 0, "one-slot pair ladder is not eligible")
    gamma = S2 * components / 7
    cutoff = floor(gamma / entry_margin) + 1
    forbidden = frozenset((*row["prefix"], row["apex"]))
    head_count = allowed_count(cutoff - 1, forbidden)
    require(head_count >= 2, "pair head has fewer than two labels")
    return {
        **row,
        "q1": q1,
        "q3": row["top5"][2][1],
        "q5": row["top5"][4][1],
        "gamma": gamma,
        "entry_margin": entry_margin,
        "cutoff": cutoff,
        "head_count": head_count,
        "raw_pairs": comb(head_count, 2),
        "forbidden": forbidden,
    }


class PairLedger:
    def __init__(
        self,
        carrier: list[tuple[F, F]],
        singleton: dict[int, F],
    ) -> None:
        self.carrier = carrier
        self.singleton = singleton
        self.pieces: dict[int, tuple[tuple[F, F], ...]] = {}
        self.edges: dict[tuple[int, int], F] = {}

    def bad_pieces(self, label: int) -> tuple[tuple[F, F], ...]:
        cached = self.pieces.get(label)
        if cached is not None:
            return cached
        pieces = tuple(
            piece
            for left, right in self.carrier
            for piece in R.local_bad_pieces(left, right, label)
        )
        require(
            sum((right - left for left, right in pieces), F(0))
            == self.singleton[label],
            f"piece/singleton mismatch at {label}",
        )
        self.pieces[label] = pieces
        return pieces

    @staticmethod
    def intersection_mass(
        first: tuple[tuple[F, F], ...],
        second: tuple[tuple[F, F], ...],
    ) -> F:
        i = 0
        j = 0
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

    def edge(self, first: int, second: int) -> F:
        key = tuple(sorted((first, second)))
        cached = self.edges.get(key)
        if cached is not None:
            return cached
        value = (
            self.singleton[first]
            + self.singleton[second]
            - self.intersection_mass(
                self.bad_pieces(first),
                self.bad_pieces(second),
            )
        )
        self.edges[key] = value
        return value


def exact_pair_head(
    carrier: list[tuple[F, F]],
    labels: tuple[int, ...],
) -> dict[str, object]:
    rows = R.coverages_many(carrier, list(labels))
    require(len(rows) == len(labels), "pair-head singleton vector changed")
    singleton = {label: value for value, label in rows}
    require(len(singleton) == len(labels), "duplicate pair-head label")
    ranked = tuple(sorted(labels, key=lambda label: (-singleton[label], label)))
    ledger = PairLedger(carrier, singleton)
    warm = tuple(sorted(ranked[:2]))
    best = ledger.edge(*warm)
    witness = warm
    paid = 1
    for first_index in range(len(ranked) - 1):
        first = ranked[first_index]
        if singleton[first] + singleton[ranked[first_index + 1]] <= best:
            break
        for second_index in range(first_index + 1, len(ranked)):
            second = ranked[second_index]
            if singleton[first] + singleton[second] <= best:
                break
            pair = tuple(sorted((first, second)))
            if pair == warm:
                continue
            value = ledger.edge(*pair)
            paid += 1
            if (value, tuple(-label for label in pair)) > (
                best,
                tuple(-label for label in witness),
            ):
                best = value
                witness = pair

    # Independent literal subtraction validates the winning pair union.
    direct_residual = R.subtract_local_multi(carrier, witness)
    require(
        best == interval_mass(carrier) - interval_mass(direct_residual),
        "winning pair reconstruction failed",
    )
    return {
        "head": best,
        "witness": witness,
        "paid": paid,
        "piece_labels": len(ledger.pieces),
    }


def exact_branch(row: dict[str, object]) -> dict[str, object]:
    data = cutoff_data(row)
    body = data["body"]
    apex = data["apex"]
    carrier, components, mass = R.CORE.good_norm(tuple(sorted((*body, apex))))
    require(
        components == data["components"]
        and mass == data["mass"]
        and mass > 0,
        "literal carrier reconstruction changed",
    )
    require(
        R.coverage(carrier, data["top5"][0][0]) == data["q1"],
        "ledger q1/direct singleton mismatch",
    )
    labels = tuple(
        label
        for label in range(FIRST_EXTERNAL, data["cutoff"])
        if label not in data["forbidden"]
    )
    require(len(labels) == data["head_count"], "pair-head count mismatch")
    initial = exact_pair_head(carrier, labels)
    minimal_tail = mass / 7 + data["gamma"] / data["cutoff"] + data["q1"]
    rank_gap = initial["head"] - data["q1"] - mass / 7
    rank_cutoff = None
    exact_cutoff = data["cutoff"]
    global_exact = False
    exact = initial
    if rank_gap > 0:
        rank_cutoff = ceiling(data["gamma"] / rank_gap)
        exact_cutoff = max(data["cutoff"], rank_cutoff)
        if exact_cutoff != data["cutoff"]:
            exact_labels = tuple(
                label
                for label in range(FIRST_EXTERNAL, exact_cutoff)
                if label not in data["forbidden"]
            )
            exact = exact_pair_head(carrier, exact_labels)
        final_rank_gap = exact["head"] - data["q1"] - mass / 7
        require(final_rank_gap >= rank_gap > 0, "pair rank gap decreased")
        final_rank_cutoff = ceiling(data["gamma"] / final_rank_gap)
        require(
            final_rank_cutoff <= exact_cutoff,
            "expanded pair head did not seal its own tail",
        )
        exact_tail = mass / 7 + data["gamma"] / exact_cutoff + data["q1"]
        require(
            exact_tail <= initial["head"] <= exact["head"],
            "exact pair tail did not fall below the finite witness",
        )
        pair_cap = exact["head"]
        global_exact = True
    else:
        exact_tail = minimal_tail
        pair_cap = max(initial["head"], minimal_tail)

    exact_head_count = allowed_count(exact_cutoff - 1, data["forbidden"])
    exact_raw_pairs = comb(exact_head_count, 2)
    paid_total = initial["paid"]
    piece_labels_total = initial["piece_labels"]
    if exact_cutoff != data["cutoff"]:
        paid_total += exact["paid"]
        piece_labels_total += exact["piece_labels"]
    pair_margin = 4 * mass / 7 - pair_cap
    direct_margin = mass - data["q5"] - 2 * pair_cap
    triple_cap = pair_cap + data["q3"]
    triple_margin = 5 * mass / 7 - triple_cap
    quadruple_cap = 2 * pair_cap
    quadruple_margin = 6 * mass / 7 - quadruple_cap
    pair_eligible = pair_margin > 0
    direct_closed = direct_margin > 0
    triple_eligible = triple_margin > 0
    quadruple_eligible = quadruple_margin > 0

    h3_cutoff = None
    h3_count = None
    h3_raw = None
    if pair_eligible:
        beta3 = pair_margin / 3
        h3_cutoff = ceiling(data["gamma"] / beta3) - 1
        h3_count = allowed_count(h3_cutoff, data["forbidden"])
        h3_raw = comb(h3_count, 3) if h3_count >= 3 else 0

    h2_cutoff = None
    h2_count = None
    h2_raw = None
    if triple_eligible:
        beta2 = triple_margin / 2
        h2_cutoff = ceiling(data["gamma"] / beta2) - 1
        h2_count = allowed_count(h2_cutoff, data["forbidden"])
        h2_raw = comb(h2_count, 4) if h2_count >= 4 else 0

    h1_cutoff = None
    h1_count = None
    h1_raw = None
    if quadruple_eligible:
        beta1 = quadruple_margin
        h1_cutoff = ceiling(data["gamma"] / beta1) - 1
        h1_count = allowed_count(h1_cutoff, data["forbidden"])
        h1_raw = comb(h1_count, 5) if h1_count >= 5 else 0

    if direct_closed:
        adaptive = "direct"
        adaptive_raw = 0
    else:
        routes = []
        if pair_eligible:
            routes.append((h3_raw, "H3"))
        if triple_eligible:
            routes.append((h2_raw, "H2"))
        if quadruple_eligible:
            routes.append((h1_raw, "H1"))
        if routes:
            adaptive_raw, adaptive = min(routes)
        else:
            adaptive = "open"
            adaptive_raw = None

    return {
        **data,
        **exact,
        "initial_head": initial["head"],
        "initial_witness": initial["witness"],
        "initial_paid": initial["paid"],
        "minimal_tail": minimal_tail,
        "rank_gap": rank_gap,
        "rank_cutoff": rank_cutoff,
        "exact_cutoff": exact_cutoff,
        "exact_head_count": exact_head_count,
        "exact_raw_pairs": exact_raw_pairs,
        "exact_tail": exact_tail,
        "global_exact": global_exact,
        "paid_total": paid_total,
        "piece_labels_total": piece_labels_total,
        "pair_cap": pair_cap,
        "pair_margin": pair_margin,
        "direct_margin": direct_margin,
        "triple_cap": triple_cap,
        "triple_margin": triple_margin,
        "quadruple_cap": quadruple_cap,
        "quadruple_margin": quadruple_margin,
        "pair_eligible": pair_eligible,
        "direct_closed": direct_closed,
        "triple_eligible": triple_eligible,
        "quadruple_eligible": quadruple_eligible,
        "h3_cutoff": h3_cutoff,
        "h3_count": h3_count,
        "h3_raw": h3_raw,
        "h2_cutoff": h2_cutoff,
        "h2_count": h2_count,
        "h2_raw": h2_raw,
        "h1_cutoff": h1_cutoff,
        "h1_count": h1_count,
        "h1_raw": h1_raw,
        "adaptive": adaptive,
        "adaptive_raw": adaptive_raw,
    }


def result_line(row: dict[str, object]) -> str:
    return (
        f"S={row['stratum']};E={','.join(map(str, row['body']))};"
        f"K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"h={ftext(row['mass'])};r={row['components']};"
        f"q1={ftext(row['q1'])};q3={ftext(row['q3'])};q5={ftext(row['q5'])};"
        f"W2={row['cutoff']};n2={row['head_count']};raw2={row['raw_pairs']};"
        f"H2initial={ftext(row['initial_head'])};"
        f"pair_initial={row['initial_witness']};paid_initial={row['initial_paid']};"
        f"tail_min={ftext(row['minimal_tail'])};"
        f"eta2={ftext(row['rank_gap'])};R2={row['rank_cutoff']};"
        f"X2={row['exact_cutoff']};nx2={row['exact_head_count']};"
        f"rawx2={row['exact_raw_pairs']};"
        f"H2={ftext(row['head'])};pair={row['witness']};"
        f"paid={row['paid_total']};tail_exact={ftext(row['exact_tail'])};"
        f"B2={ftext(row['pair_cap'])};global_exact={int(row['global_exact'])};"
        f"mB2={ftext(row['pair_margin'])};"
        f"mdirect={ftext(row['direct_margin'])};"
        f"mB3={ftext(row['triple_margin'])};"
        f"mB4={ftext(row['quadruple_margin'])};"
        f"H3={row['h3_cutoff']}:{row['h3_count']}:{row['h3_raw']};"
        f"H2core={row['h2_cutoff']}:{row['h2_count']}:{row['h2_raw']};"
        f"H1={row['h1_cutoff']}:{row['h1_count']}:{row['h1_raw']};"
        f"adaptive={row['adaptive']}:{row['adaptive_raw']}\n"
    )


def nearest_rank_quantiles(
    values: list[int],
    percentages: tuple[int, ...] = (0, 25, 50, 75, 90, 95, 99, 100),
) -> tuple[tuple[int, int], ...]:
    require(values, "empty quantile list")
    ordered = sorted(values)
    return tuple(
        (
            percentage,
            ordered[
                0
                if percentage == 0
                else (percentage * len(ordered) + 99) // 100 - 1
            ],
        )
        for percentage in percentages
    )


def row_id(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["rank"],
        row["apex"],
        row["prefix"],
    )


def extremum(
    name: str,
    rows: list[dict[str, object]],
    field: str,
    minimum: bool = True,
) -> tuple[object, ...]:
    require(rows, f"empty extremum: {name}")
    if minimum:
        row = min(rows, key=lambda item: (item[field], row_id(item)))
    else:
        row = max(
            rows,
            key=lambda item: (item[field], tuple(-x for x in item["body"]), -item["rank"]),
        )
    value = row[field]
    value_text = ftext(value) if isinstance(value, F) else value
    return (name, value_text, *row_id(row))


def failure_digest(
    name: str,
    rows: list[dict[str, object]],
) -> tuple[str, str]:
    digest = hashlib.sha256()
    digest.update(f"LRC14/j6/global-pair/{name}/v1\n".encode())
    for row in rows:
        digest.update(result_line(row).encode())
    return name, digest.hexdigest()


def group_summary(
    key: object,
    rows: list[dict[str, object]],
) -> tuple[object, ...]:
    return (
        key,
        len(rows),
        sum(row["pair_eligible"] for row in rows),
        sum(row["direct_closed"] for row in rows),
        sum(row["triple_eligible"] for row in rows),
        sum(row["quadruple_eligible"] for row in rows),
        sum(row["adaptive"] != "open" for row in rows),
        sum(row["adaptive"] == "direct" for row in rows),
        sum(row["adaptive"] == "H1" for row in rows),
        sum(row["adaptive"] == "H2" for row in rows),
        sum(row["adaptive"] == "H3" for row in rows),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
    )
    parser.add_argument(
        "--workload-only",
        action="store_true",
        help="parse the hard ledger and report the cheap W2 workload only",
    )
    parser.add_argument(
        "--emit-rows",
        action="store_true",
        help="emit the complete exact branch ledger",
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    inputs = [cutoff_data(row) for row in parse_hard_ledger()]

    if args.workload_only:
        print("LRC14 j6 all-hard global pair-cap workload scout")
        print(
            "counts="
            f"{(len(inputs), min(row['cutoff'] for row in inputs), max(row['cutoff'] for row in inputs), sum(row['raw_pairs'] for row in inputs))}"
        )
        print(
            "W2_quantiles_nearest_rank="
            f"{nearest_rank_quantiles([row['cutoff'] for row in inputs])}"
        )
        print(
            "head_n_quantiles_nearest_rank="
            f"{nearest_rank_quantiles([row['head_count'] for row in inputs])}"
        )
        print(
            "raw_pair_quantiles_nearest_rank="
            f"{nearest_rank_quantiles([row['raw_pairs'] for row in inputs])}"
        )
        print("scope=cheap workload only;no pair unions evaluated;not LRC14")
        return

    context = mp.get_context("spawn")
    if args.workers == 1:
        rows = list(map(exact_branch, inputs))
    else:
        with context.Pool(args.workers) as pool:
            rows = pool.map(exact_branch, inputs, chunksize=8)
    require(len(rows) == 14_806, "exact result count changed")
    if args.emit_rows:
        for row in rows:
            print("PAIR;" + result_line(row).rstrip())

    groups: list[tuple[object, ...]] = []
    for stratum in ("low", "one", "both"):
        groups.append(
            group_summary(
                f"S:{stratum}",
                [row for row in rows if row["stratum"] == stratum],
            )
        )
    for rank in sorted({row["rank"] for row in rows}):
        groups.append(
            group_summary(
                f"R:{rank}",
                [row for row in rows if row["rank"] == rank],
            )
        )
    groups_tuple = tuple(groups)

    pair_rows = [row for row in rows if row["pair_eligible"]]
    direct_rows = [row for row in rows if row["direct_closed"]]
    triple_rows = [row for row in rows if row["triple_eligible"]]
    quadruple_rows = [row for row in rows if row["quadruple_eligible"]]
    open_rows = [row for row in rows if row["adaptive"] == "open"]
    counts = (
        len(rows),
        len(pair_rows),
        len(rows) - len(pair_rows),
        len(direct_rows),
        len(triple_rows),
        len(quadruple_rows),
        sum(row["direct_closed"] or row["triple_eligible"] for row in rows),
        sum(
            row["direct_closed"]
            or row["triple_eligible"]
            or row["quadruple_eligible"]
            for row in rows
        ),
        len(rows) - len(open_rows),
        len(open_rows),
        sum(row["global_exact"] for row in rows),
        sum(row["paid_total"] for row in rows),
        sum(row["exact_raw_pairs"] for row in rows),
        sum(row["adaptive"] == "direct" for row in rows),
        sum(row["adaptive"] == "H1" for row in rows),
        sum(row["adaptive"] == "H2" for row in rows),
        sum(row["adaptive"] == "H3" for row in rows),
    )

    extrema = (
        extremum("minimum_pair_margin", rows, "pair_margin"),
        extremum("minimum_direct_margin", rows, "direct_margin"),
        extremum("minimum_triple_margin", rows, "triple_margin"),
        extremum("minimum_quadruple_margin", rows, "quadruple_margin"),
        extremum("maximum_W2", rows, "cutoff", minimum=False),
        extremum("maximum_initial_head_n", rows, "head_count", minimum=False),
        extremum("maximum_exact_cutoff", rows, "exact_cutoff", minimum=False),
        extremum("maximum_exact_head_n", rows, "exact_head_count", minimum=False),
        extremum("maximum_exact_raw_pairs", rows, "exact_raw_pairs", minimum=False),
        extremum("maximum_paid_pairs", rows, "paid_total", minimum=False),
        extremum("maximum_H3_cutoff", pair_rows, "h3_cutoff", minimum=False),
        extremum("maximum_H3_raw", pair_rows, "h3_raw", minimum=False),
        extremum("maximum_H2_cutoff", triple_rows, "h2_cutoff", minimum=False),
        extremum("maximum_H2_raw", triple_rows, "h2_raw", minimum=False),
        extremum(
            "maximum_H1_cutoff",
            quadruple_rows,
            "h1_cutoff",
            minimum=False,
        ),
        extremum(
            "maximum_H1_raw",
            quadruple_rows,
            "h1_raw",
            minimum=False,
        ),
    )
    quantiles = (
        ("W2", nearest_rank_quantiles([row["cutoff"] for row in rows])),
        ("head_n", nearest_rank_quantiles([row["head_count"] for row in rows])),
        ("raw_pairs", nearest_rank_quantiles([row["raw_pairs"] for row in rows])),
        ("R2", nearest_rank_quantiles([row["rank_cutoff"] for row in rows if row["rank_cutoff"] is not None])),
        ("exact_cutoff", nearest_rank_quantiles([row["exact_cutoff"] for row in rows])),
        ("exact_head_n", nearest_rank_quantiles([row["exact_head_count"] for row in rows])),
        ("exact_raw_pairs", nearest_rank_quantiles([row["exact_raw_pairs"] for row in rows])),
        ("paid_pairs", nearest_rank_quantiles([row["paid_total"] for row in rows])),
        ("H3_cutoff", nearest_rank_quantiles([row["h3_cutoff"] for row in pair_rows])),
        ("H2_cutoff", nearest_rank_quantiles([row["h2_cutoff"] for row in triple_rows])),
        (
            "H1_cutoff",
            nearest_rank_quantiles(
                [
                    row["h1_cutoff"]
                    for row in quadruple_rows
                ]
            ),
        ),
    )
    failure_digests = (
        failure_digest(
            "pair_failures",
            [row for row in rows if not row["pair_eligible"]],
        ),
        failure_digest(
            "direct_failures",
            [row for row in rows if not row["direct_closed"]],
        ),
        failure_digest(
            "triple_failures",
            [row for row in rows if not row["triple_eligible"]],
        ),
        failure_digest(
            "quadruple_failures",
            [row for row in rows if not row["quadruple_eligible"]],
        ),
        failure_digest("adaptive_open", open_rows),
    )
    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-hard-global-pair-cap/v1\n")
    for row in rows:
        digest.update(result_line(row).encode())
    ledger_digest = digest.hexdigest()

    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "aggregate counts changed")
    if EXPECTED_GROUPS is not None:
        require(groups_tuple == EXPECTED_GROUPS, "group census changed")
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_QUANTILES is not None:
        require(quantiles == EXPECTED_QUANTILES, "quantiles changed")
    if EXPECTED_FAILURE_DIGESTS is not None:
        require(
            failure_digests == EXPECTED_FAILURE_DIGESTS,
            "failure ledgers changed",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(ledger_digest == EXPECTED_LEDGER_DIGEST, "result ledger changed")

    print("LRC14 j6 all-hard globally sealed pair-cap census")
    print(f"counts={counts}")
    print(f"groups={groups_tuple}")
    print(f"extrema={extrema}")
    print(f"quantiles_nearest_rank={quantiles}")
    print(f"failure_digests={failure_digests}")
    print(f"ledger_sha256={ledger_digest}")
    print("mode=DISCOVERY" if EXPECTED_COUNTS is None else "mode=LOCKED")
    print(
        "scope=all 14806 scalar-hard marked THM2896 suffixes;"
        "exact finite pair heads plus strict discrepancy tails;"
        "global pair caps and heavy-core workload only;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
