#!/usr/bin/env python3
"""Batched exact replacement for graph_free_star_peel.py.

The proof algorithm, query order, digests, and output schema are inherited
from the hash-pinned scalar reference.  Only the evaluation of a pivot's
child coverages is changed: prospective scalar queries are evaluated in
guarded int64 batches by the exact ``coverages_many`` primitive, then
consumed in the original order.  Values speculatively computed past an
early stopping point are discarded before the next pivot and never enter
the proof ledger or cache-hit counts.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import importlib.util
import math
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_PATH = (
    ROOT
    / ".scratch/lrc_rank3_h2_clique_descent_20260729/graph_free_star_peel.py"
)
REFERENCE_SHA256 = (
    "ae59a1b28ea1f216b23fa0de58bcd06410f2419287130ac75aa3157c37daef4c"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_reference():
    require(
        hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest()
        == REFERENCE_SHA256,
        "scalar graph-free reference changed",
    )
    spec = importlib.util.spec_from_file_location(
        "rank3_h2_graph_free_scalar_reference",
        REFERENCE_PATH,
    )
    require(spec is not None and spec.loader is not None, "cannot load reference")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = load_reference()


class ExactCoverageKernel:
    """One residual's exact common-denominator tooth primitive.

    Safe rows use the same guarded int64 vector primitive as the canonical
    engine.  Rows outside that guard use the identical integer numerator
    formula with Python arbitrary-precision integers.
    """

    def __init__(self, carrier: list[tuple[F, F]]):
        require(carrier, "empty vector-coverage carrier")
        self.carrier = carrier
        self.denominator = math.lcm(
            *(
                endpoint.denominator
                for interval in carrier
                for endpoint in interval
            )
        )
        self.endpoint_ints = tuple(
            endpoint.numerator
            * (self.denominator // endpoint.denominator)
            for interval in carrier
            for endpoint in interval
        )

    def int64_safe(self, speeds: list[int]) -> bool:
        maximum = max(speeds)
        return (
            4
            * len(self.carrier)
            * self.denominator
            * maximum
            < 2**62
            and 14 * self.denominator * maximum < 2**63
        )

    def many(self, speeds: list[int]) -> tuple[list[tuple[F, int]], str]:
        require(speeds, "empty exact-coverage batch")
        if self.int64_safe(speeds):
            # Reuse the canonical implementation, including both overflow
            # guards.  Its output consists of exact Fraction values.
            return G.P.S.G.T.coverages_many(self.carrier, speeds), "int64"

        denominator = self.denominator

        def primitive(scaled: int) -> int:
            quotient, remainder = divmod(scaled, denominator)
            return (
                2 * denominator * quotient
                + min(14 * remainder, denominator)
                + max(0, 14 * remainder - 13 * denominator)
            )

        out = []
        for speed in speeds:
            numerator = 0
            for index in range(0, len(self.endpoint_ints), 2):
                numerator += primitive(self.endpoint_ints[index + 1] * speed)
                numerator -= primitive(self.endpoint_ints[index] * speed)
            out.append((F(numerator, 14 * denominator * speed), speed))
        return out, "bigint"


class BatchedBranchPeeler(G.BranchPeeler):
    """Reference peeler with speculative vector batches behind pair queries."""

    def __init__(
        self,
        row: dict[str, object],
        coverages: dict[int, F],
        pair_cap: F,
        batch_size: int,
    ):
        super().__init__(row, coverages, pair_cap)
        require(batch_size >= 1, "batch size must be positive")
        self.batch_size = batch_size
        self.raw_batch: dict[tuple[int, int], F] = {}
        self.vector_batches = 0
        self.vector_values = 0
        self.vector_controls = 0
        self.speculative_discards = 0
        self.bigint_batch_size = min(8, batch_size)
        self.coverage_kernels: dict[int, ExactCoverageKernel] = {}
        self.int64_batches = 0
        self.int64_values = 0
        self.bigint_batches = 0
        self.bigint_values = 0

    def clear_raw_batch(self) -> None:
        self.speculative_discards += len(self.raw_batch)
        self.raw_batch.clear()

    def prepare_pair_batch(self, pivot: int, speeds: tuple[int, ...]) -> None:
        missing = []
        seen = set()
        for speed in speeds:
            key = tuple(sorted((pivot, speed)))
            if key in seen or key in self.pairs or key in self.raw_batch:
                continue
            seen.add(key)
            missing.append(speed)
        if not missing:
            return

        carrier = self.reconstruct_carrier()
        if pivot not in self.after:
            self.after[pivot] = G.P.S.R.subtract_local(carrier, pivot)
        if pivot not in self.coverage_kernels:
            self.coverage_kernels[pivot] = ExactCoverageKernel(self.after[pivot])
        kernel = self.coverage_kernels[pivot]
        if not kernel.int64_safe(missing):
            missing = missing[: self.bigint_batch_size]
        child_rows, mode = kernel.many(missing)
        require(
            len(child_rows) == len(missing)
            and tuple(speed for _, speed in child_rows) == tuple(missing),
            "vector child-coverage order changed",
        )
        for child, speed in child_rows:
            key = tuple(sorted((pivot, speed)))
            union = self.coverages[pivot] + child
            require(
                max(self.coverages[pivot], self.coverages[speed]) <= union
                <= min(
                    self.h,
                    self.coverages[pivot] + self.coverages[speed],
                    self.pair_cap,
                ),
                "vector pair-union bounds failed",
            )
            self.raw_batch[key] = union

        # Independent scalar controls on the first, middle, and last entries
        # of the first batch, then one endpoint every 128 batches.  These do
        # not enter the theorem-bearing query ledger.
        control_indices: tuple[int, ...]
        if self.vector_batches == 0:
            control_indices = tuple(
                dict.fromkeys((0, len(missing) // 2, len(missing) - 1))
            )
        elif self.vector_batches % 128 == 0:
            control_indices = (len(missing) - 1,)
        else:
            control_indices = ()
        for index in control_indices:
            speed = missing[index]
            scalar = (
                self.coverages[pivot]
                + G.P.S.G.T.coverage(self.after[pivot], speed)
            )
            require(
                self.raw_batch[tuple(sorted((pivot, speed)))] == scalar,
                f"vector/scalar pair mismatch: {(pivot, speed)}",
            )
            self.vector_controls += 1
        self.vector_batches += 1
        self.vector_values += len(missing)
        if mode == "int64":
            self.int64_batches += 1
            self.int64_values += len(missing)
        else:
            self.bigint_batches += 1
            self.bigint_values += len(missing)

    def pair_union(self, x: int, y: int) -> F:
        key = tuple(sorted((x, y)))
        cached = self.pairs.get(key)
        if cached is not None:
            self.pair_cache_hits += 1
            return cached
        if key not in self.raw_batch:
            self.prepare_pair_batch(x, (y,))
        union = self.raw_batch.pop(key)
        require(
            max(self.coverages[x], self.coverages[y]) <= union
            <= min(self.h, self.coverages[x] + self.coverages[y], self.pair_cap),
            "pair-union bounds failed",
        )
        self.pairs[key] = union
        self.pair_queries += 1
        self.first_pair = key if self.first_pair is None else self.first_pair
        self.last_pair = key
        heavy = union >= self.theta
        self.threshold_ties += union == self.theta
        self.pair_digest.update(
            (
                f"{key[0]},{key[1]}:{G.ftext(union)}:"
                f"{int(heavy)}\n"
            ).encode()
        )
        candidate = (abs(union - self.theta), key)
        if heavy:
            if self.closest_heavy is None or candidate < self.closest_heavy:
                self.closest_heavy = candidate
        elif self.closest_light is None or candidate < self.closest_light:
            self.closest_light = candidate
        return union

    def _test_pivot(self, pivot: int) -> dict[str, object]:
        self.clear_raw_batch()
        cx = self.coverages[pivot]
        heavy_floor = self.theta - cx
        upper_cap = self.pair_cap - cx
        require(upper_cap >= heavy_floor, "pair cap would globally delete heavy edges")

        upper = self.upper_three(pivot)
        if len(upper) < 3:
            return {
                "certified": True,
                "kind": "FEWER_THAN_THREE_VERTICES",
                "margin": None,
                "triple": (),
                "scanned": 0,
            }
        upper_margin = self.h - cx - self.q1 - sum(upper, F(0))
        if upper_margin > 0:
            return {
                "certified": True,
                "kind": "PAIR_CAP_UPPER",
                "margin": upper_margin,
                "triple": (),
                "scanned": 0,
            }

        ordered = tuple(speed for speed in self.order if speed != pivot)
        top: list[tuple[F, int]] = []
        scanned = 0
        cursor = 0
        while cursor < len(ordered):
            speed = ordered[cursor]
            upper_bound = min(self.coverages[speed], upper_cap)
            if len(top) == 3 and upper_bound <= top[-1][0]:
                margin = self.h - cx - self.q1 - sum(
                    (value for value, _ in top), F(0)
                )
                require(margin > 0, "hard exact top three escaped early")
                return {
                    "certified": True,
                    "kind": "EXACT_TOP_THREE",
                    "margin": margin,
                    "triple": tuple(item for _, item in top),
                    "scanned": scanned,
                }
            if upper_bound < heavy_floor:
                break

            if tuple(sorted((pivot, speed))) not in self.pairs:
                prospective = []
                for candidate in ordered[cursor : cursor + self.batch_size]:
                    if min(self.coverages[candidate], upper_cap) < heavy_floor:
                        break
                    prospective.append(candidate)
                self.prepare_pair_batch(pivot, tuple(prospective))

            union = self.pair_union(pivot, speed)
            scanned += 1
            cursor += 1
            if union < self.theta:
                continue
            child = union - cx
            top.append((child, speed))
            top.sort(key=lambda item: (-item[0], item[1]))
            del top[3:]
            if len(top) == 3:
                margin = self.h - cx - self.q1 - sum(
                    (value for value, _ in top), F(0)
                )
                if margin <= 0:
                    return {
                        "certified": False,
                        "kind": "HOSTILE_TRIPLE",
                        "margin": margin,
                        "triple": tuple(item for _, item in top),
                        "scanned": scanned,
                    }

        if len(top) < 3:
            return {
                "certified": True,
                "kind": "FEWER_THAN_THREE_HEAVY",
                "margin": None,
                "triple": tuple(item for _, item in top),
                "scanned": scanned,
            }
        margin = self.h - cx - self.q1 - sum(
            (value for value, _ in top), F(0)
        )
        require(margin > 0, "terminal exact top three was not strict")
        return {
            "certified": True,
            "kind": "EXACT_TOP_THREE",
            "margin": margin,
            "triple": tuple(item for _, item in top),
            "scanned": scanned,
        }

    def test_pivot(self, pivot: int) -> dict[str, object]:
        try:
            return self._test_pivot(pivot)
        finally:
            # Pivot-local carriers and scaled endpoints can be reconstructed
            # if a final independent control needs them.  They never affect a
            # later proof query, so retaining thousands of copies on a large
            # branch only wastes memory.
            self.clear_raw_batch()
            self.coverage_kernels.pop(pivot, None)
            self.after.pop(pivot, None)

    def run(self) -> dict[str, object]:
        row = super().run()
        self.clear_raw_batch()
        row.update(
            {
                "vector_batches": self.vector_batches,
                "vector_values": self.vector_values,
                "vector_controls": self.vector_controls,
                "speculative_discards": self.speculative_discards,
                "int64_batches": self.int64_batches,
                "int64_values": self.int64_values,
                "bigint_batches": self.bigint_batches,
                "bigint_values": self.bigint_values,
            }
        )
        return row


def profile(
    task: tuple[dict[str, object], dict[int, F], F, int],
) -> dict[str, object]:
    row, coverages, pair_cap, batch_size = task
    return BatchedBranchPeeler(
        row,
        coverages,
        pair_cap,
        batch_size,
    ).run()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-core", type=int, default=41)
    parser.add_argument("--max-core", type=int)
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--emit-ledger", action="store_true")
    parser.add_argument("--emit-performance", action="store_true")
    args = parser.parse_args()
    require(args.min_core >= 0, "min-core must be nonnegative")
    require(
        args.max_core is None or args.max_core >= args.min_core,
        "bad core range",
    )
    require(args.workers >= 1, "workers must be positive")
    require(args.batch_size >= 1, "batch size must be positive")

    rows = [
        row
        for row in G.P.TARGET_ROWS
        if len(G.CORES[(row["body"], row["rank"])]) >= args.min_core
        and (
            args.max_core is None
            or len(G.CORES[(row["body"], row["rank"])]) <= args.max_core
        )
    ]
    tasks = tuple(
        (
            row,
            G.CORES[(row["body"], row["rank"])],
            G.PAIR_CAPS[(row["body"], row["rank"])],
            args.batch_size,
        )
        for row in rows
    )
    if args.workers == 1:
        profiled = [profile(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiled = list(pool.imap(profile, tasks, chunksize=1))
    require(
        tuple((row["body"], row["rank"]) for row in profiled)
        == tuple((row["body"], row["rank"]) for row in rows),
        "worker order changed",
    )

    digest = hashlib.sha256(
        (
            f"LRC14/j6/r3/H2/graph-free-star-peel/"
            f"{args.min_core}/{args.max_core}/v1\n"
        ).encode()
    )
    for row in profiled:
        line = G.ledger_line(row)
        digest.update(line.encode())
        if args.emit_ledger:
            print("ROW;" + line.rstrip())

    histogram = tuple(
        sorted(Counter(len(row["residual"]) for row in profiled).items())
    )
    print("LRC14 j6 strongest-r3 graph-free ordered Hunter-star peel")
    print(
        f"core_range={args.min_core}.."
        f"{'-' if args.max_core is None else args.max_core};"
        f"branches={len(profiled)};closed={sum(row['closed'] for row in profiled)};"
        f"survivors={sum(not row['closed'] for row in profiled)}"
    )
    print(
        f"vertices={sum(row['H'] for row in profiled)};"
        f"peeled={sum(len(row['peeled']) for row in profiled)};"
        f"residual={sum(len(row['residual']) for row in profiled)};"
        f"residual_histogram={histogram}"
    )
    print(
        f"pair_queries={sum(row['pair_queries'] for row in profiled)};"
        f"cache_hits={sum(row['pair_cache_hits'] for row in profiled)};"
        f"candidate_scans={sum(row['total_scanned'] for row in profiled)};"
        f"retests={sum(row['retests'] for row in profiled)};"
        f"threshold_ties={sum(row['threshold_ties'] for row in profiled)};"
        f"controls={sum(row['pair_controls'] for row in profiled)}"
    )
    print(
        "residual_max="
        + repr(
            max(
                (
                    len(row["residual"]),
                    row["H"],
                    row["body"],
                    row["rank"],
                )
                for row in profiled
            )
            if profiled
            else None
        )
    )
    print(
        "closed_keys="
        + repr(
            tuple(
                (row["body"], row["rank"], row["H"])
                for row in profiled
                if row["closed"]
            )
        )
    )
    print(f"canonical_ledger_sha256={digest.hexdigest()}")
    print(
        "orientation=acyclic proof-elimination gauge only;"
        "queried pair unions symmetric;not tournament"
    )
    print(
        "scope=exact graph-free Hunter-star peel;"
        "residual size >=4 is unresolved, not a cover witness"
    )
    print("all_exact_controls=PASS")
    if args.emit_performance:
        print(
            "vector_performance="
            f"(batch_size={args.batch_size},"
            f"batches={sum(row['vector_batches'] for row in profiled)},"
            f"values={sum(row['vector_values'] for row in profiled)},"
            f"controls={sum(row['vector_controls'] for row in profiled)},"
            f"speculative_discards={sum(row['speculative_discards'] for row in profiled)},"
            f"int64_batches={sum(row['int64_batches'] for row in profiled)},"
            f"int64_values={sum(row['int64_values'] for row in profiled)},"
            f"bigint_batches={sum(row['bigint_batches'] for row in profiled)},"
            f"bigint_values={sum(row['bigint_values'] for row in profiled)})"
        )


if __name__ == "__main__":
    main()
