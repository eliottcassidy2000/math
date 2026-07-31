#!/usr/bin/env python3
"""Exact graph-free ordered Hunter-star peeling on the ranked-r3 H2 cores.

The script never constructs the theta-heavy graph.  For a pivot x it seeks
only three heavy later neighbors whose child coverages already make the
Hunter-star invoice hard, or else proves the exact top-three child invoice
strictly below h.  A failed pivot watches one explicit hostile triple and
is reconsidered only when a member of that triple is peeled.

This is scratch discovery code until its census is locked and promoted.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import heapq
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/pilot.py"
PILOT_SHA256 = (
    "9ec1b8d8c697f54fdbf2836638fa6c7dc5284f4ea2644849d98d71398c1f3520"
)
CORE_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/core.discovery.out"
CORE_SHA256 = (
    "5569d8d34b59e5eed2bfa82148648dcb5e63515146ff55b95234a002d0c87ba2"
)
PREFILTER_PATH = (
    ROOT
    / ".scratch/lrc_rank3_h2_clique_descent_20260729/prefilter.discovery.out"
)
PREFILTER_SHA256 = (
    "3727e5a285a289ab74fca6f7ce567e60c66f85439e3d13a62da3fa4fa5a743d7"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_pilot():
    require(file_sha256(PILOT_PATH) == PILOT_SHA256, "pilot changed")
    require(file_sha256(CORE_PATH) == CORE_SHA256, "core ledger changed")
    require(file_sha256(PREFILTER_PATH) == PREFILTER_SHA256, "prefilter changed")
    spec = importlib.util.spec_from_file_location("rank3_h2_graph_free_pilot", PILOT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load pilot")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_pilot()


def ftext(value: F | None) -> str:
    return "-" if value is None else f"{value.numerator}/{value.denominator}"


def parse_inputs() -> tuple[dict[tuple[tuple[int, ...], int], dict[int, F]], dict[tuple[tuple[int, ...], int], F]]:
    cores = {}
    for line in CORE_PATH.read_text().splitlines():
        if not line.startswith("H2;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (ast.literal_eval(fields["E"]), int(fields["rank"]))
        cores[key] = {
            int(speed): F(value)
            for speed, value in ast.literal_eval(fields["Hrows"])
        }

    pair_caps = {}
    for line in PREFILTER_PATH.read_text().splitlines():
        if not line.startswith("ROW;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (ast.literal_eval(fields["E"]), int(fields["rank"]))
        pair_caps[key] = F(fields["B2"])

    expected = {(row["body"], row["rank"]) for row in P.TARGET_ROWS}
    require(set(cores) == expected, "core universe changed")
    require(set(pair_caps) == expected, "pair-cap universe changed")
    return cores, pair_caps


CORES, PAIR_CAPS = parse_inputs()


class CoverageOrder:
    """Coverage-descending active linked list with O(1) deletions."""

    def __init__(self, coverages: dict[int, F]):
        self.vertices = tuple(
            sorted(coverages, key=lambda speed: (-coverages[speed], speed))
        )
        self.position = {speed: index for index, speed in enumerate(self.vertices)}
        size = len(self.vertices)
        self.previous = [index - 1 for index in range(size)]
        self.following = [index + 1 for index in range(size)]
        if size:
            self.following[-1] = -1
        self.active = [True] * size
        self.head = 0 if size else -1

    def __iter__(self):
        index = self.head
        while index >= 0:
            require(self.active[index], "inactive vertex remained in linked order")
            yield self.vertices[index]
            index = self.following[index]

    def remove(self, speed: int) -> None:
        index = self.position[speed]
        require(self.active[index], "linked-order vertex removed twice")
        before = self.previous[index]
        after = self.following[index]
        if before >= 0:
            self.following[before] = after
        else:
            self.head = after
        if after >= 0:
            self.previous[after] = before
        self.active[index] = False


class BranchPeeler:
    def __init__(self, row: dict[str, object], coverages: dict[int, F], pair_cap: F):
        self.row = row
        self.coverages = coverages
        self.pair_cap = pair_cap
        self.h = row["h"]
        self.q1 = row["top5"][0][0]
        self.r3 = sum((value for value, _ in row["top5"][:3]), F(0))
        self.theta = self.h - self.r3
        self.order = CoverageOrder(coverages)
        self.remaining = set(coverages)
        self.carrier = None
        self.after = {}
        self.pairs = {}
        self.pair_queries = 0
        self.pair_cache_hits = 0
        self.threshold_ties = 0
        self.first_pair = None
        self.last_pair = None
        self.closest_heavy = None
        self.closest_light = None
        self.pair_digest = hashlib.sha256(
            (
                f"LRC14/j6/r3/H2/graph-free-pairs/{row['body']}/"
                f"{row['rank']}/v1\n"
            ).encode()
        )
        self.event_digest = hashlib.sha256(
            (
                f"LRC14/j6/r3/H2/graph-free-star-events/{row['body']}/"
                f"{row['rank']}/v1\n"
            ).encode()
        )

    def reconstruct_carrier(self):
        if self.carrier is None:
            root_good, _, _ = P.S.G.T.CORE.good_norm(self.row["body"])
            self.carrier = P.S.R.subtract_local(root_good, self.row["apex"])
            require(P.mass(self.carrier) == self.h, "literal carrier changed")
        return self.carrier

    def pair_union(self, x: int, y: int) -> F:
        key = tuple(sorted((x, y)))
        cached = self.pairs.get(key)
        if cached is not None:
            self.pair_cache_hits += 1
            return cached
        carrier = self.reconstruct_carrier()
        if x not in self.after:
            self.after[x] = P.S.R.subtract_local(carrier, x)
        union = self.coverages[x] + P.S.G.T.coverage(self.after[x], y)
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
                f"{key[0]},{key[1]}:{ftext(union)}:"
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

    def upper_three(self, pivot: int) -> tuple[F, ...]:
        cap = self.pair_cap - self.coverages[pivot]
        out = []
        for speed in self.order:
            if speed == pivot:
                continue
            out.append(min(self.coverages[speed], cap))
            if len(out) == 3:
                break
        return tuple(out)

    def test_pivot(self, pivot: int) -> dict[str, object]:
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

        top: list[tuple[F, int]] = []
        scanned = 0
        for speed in self.order:
            if speed == pivot:
                continue
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
                    "triple": tuple(speed for _, speed in top),
                    "scanned": scanned,
                }
            if upper_bound < heavy_floor:
                break

            union = self.pair_union(pivot, speed)
            scanned += 1
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
                        "triple": tuple(speed for _, speed in top),
                        "scanned": scanned,
                    }

        if len(top) < 3:
            return {
                "certified": True,
                "kind": "FEWER_THAN_THREE_HEAVY",
                "margin": None,
                "triple": tuple(speed for _, speed in top),
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
            "triple": tuple(speed for _, speed in top),
            "scanned": scanned,
        }

    def control_pair(self, key: tuple[int, int]) -> None:
        x, y = key
        carrier = self.reconstruct_carrier()
        if y not in self.after:
            self.after[y] = P.S.R.subtract_local(carrier, y)
        reverse = self.coverages[y] + P.S.G.T.coverage(self.after[y], x)
        literal = self.h - P.mass(P.S.R.subtract_local_multi(carrier, key))
        require(
            self.pairs[key] == reverse == literal,
            f"pair scalar/reverse/literal mismatch: {key}",
        )

    def run(self) -> dict[str, object]:
        witnesses: dict[int, tuple[int, int, int]] = {}
        watchers: dict[int, set[int]] = defaultdict(set)
        pending = set(self.remaining)
        queue = [
            (self.coverages[speed], speed)
            for speed in self.remaining
        ]
        heapq.heapify(queue)
        peel = []
        retests = 0
        total_scanned = 0

        while queue:
            _, pivot = heapq.heappop(queue)
            if pivot not in pending:
                continue
            pending.remove(pivot)
            if pivot not in self.remaining:
                continue
            retests += pivot in witnesses
            result = self.test_pivot(pivot)
            total_scanned += result["scanned"]
            self.event_digest.update(
                (
                    f"test={pivot};kind={result['kind']};"
                    f"margin={ftext(result['margin'])};"
                    f"triple={result['triple']};scanned={result['scanned']}\n"
                ).encode()
            )
            if not result["certified"]:
                triple = result["triple"]
                require(len(triple) == 3, "failed pivot lacks hostile triple")
                witnesses[pivot] = triple
                for speed in triple:
                    watchers[speed].add(pivot)
                continue

            peel.append(
                (
                    pivot,
                    result["kind"],
                    ftext(result["margin"]),
                    result["triple"],
                    result["scanned"],
                )
            )
            witnesses.pop(pivot, None)
            self.remaining.remove(pivot)
            self.order.remove(pivot)
            for dependent in tuple(watchers.pop(pivot, ())):
                if (
                    dependent in self.remaining
                    and pivot in witnesses.get(dependent, ())
                ):
                    witnesses.pop(dependent)
                    if dependent not in pending:
                        pending.add(dependent)
                        heapq.heappush(
                            queue,
                            (self.coverages[dependent], dependent),
                        )

        require(not pending, "pending pivot queue did not drain")
        for pivot in self.remaining:
            triple = witnesses.get(pivot)
            require(triple is not None and len(triple) == 3, "residual lacks witness")
            require(all(speed in self.remaining for speed in triple), "stale witness")
            children = []
            for speed in triple:
                union = self.pairs[tuple(sorted((pivot, speed)))]
                require(union >= self.theta, "witness edge is not heavy")
                children.append(union - self.coverages[pivot])
            require(
                self.coverages[pivot] + self.q1 + sum(children, F(0)) >= self.h,
                "residual hostile triple became strict",
            )

        controls = {
            key
            for key in (
                self.first_pair,
                self.last_pair,
                None if self.closest_heavy is None else self.closest_heavy[1],
                None if self.closest_light is None else self.closest_light[1],
            )
            if key is not None
        }
        for key in sorted(controls):
            self.control_pair(key)

        residual = tuple(sorted(self.remaining))
        residual_witnesses = tuple(
            (pivot, witnesses[pivot])
            for pivot in residual
        )
        self.event_digest.update(
            (
                f"residual={residual};witnesses={residual_witnesses};"
                f"pairs={self.pair_queries};ties={self.threshold_ties}\n"
            ).encode()
        )
        return {
            **self.row,
            "H": len(self.coverages),
            "pair_cap": self.pair_cap,
            "peeled": tuple(peel),
            "residual": residual,
            "residual_witnesses": residual_witnesses,
            "closed": len(residual) < 4,
            "pair_queries": self.pair_queries,
            "pair_cache_hits": self.pair_cache_hits,
            "threshold_ties": self.threshold_ties,
            "pair_controls": len(controls),
            "total_scanned": total_scanned,
            "retests": retests,
            "pair_digest": self.pair_digest.hexdigest(),
            "event_digest": self.event_digest.hexdigest(),
        }


def profile(task: tuple[dict[str, object], dict[int, F], F]) -> dict[str, object]:
    row, coverages, pair_cap = task
    return BranchPeeler(row, coverages, pair_cap).run()


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};H={row['H']};"
        f"B2={ftext(row['pair_cap'])};peeled={len(row['peeled'])};"
        f"residual={len(row['residual'])};closed={int(row['closed'])};"
        f"pairs={row['pair_queries']};hits={row['pair_cache_hits']};"
        f"scanned={row['total_scanned']};retests={row['retests']};"
        f"ties={row['threshold_ties']};controls={row['pair_controls']};"
        f"residual_vertices={row['residual']};"
        f"residual_witnesses={row['residual_witnesses']};"
        f"peel={row['peeled']};pair_digest={row['pair_digest']};"
        f"event_digest={row['event_digest']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-core", type=int, default=41)
    parser.add_argument("--max-core", type=int)
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--emit-ledger", action="store_true")
    args = parser.parse_args()
    require(args.min_core >= 0, "min-core must be nonnegative")
    require(args.max_core is None or args.max_core >= args.min_core, "bad core range")
    require(args.workers >= 1, "workers must be positive")

    rows = [
        row
        for row in P.TARGET_ROWS
        if len(CORES[(row["body"], row["rank"])]) >= args.min_core
        and (
            args.max_core is None
            or len(CORES[(row["body"], row["rank"])]) <= args.max_core
        )
    ]
    tasks = tuple(
        (
            row,
            CORES[(row["body"], row["rank"])],
            PAIR_CAPS[(row["body"], row["rank"])],
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
        line = ledger_line(row)
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


if __name__ == "__main__":
    main()
