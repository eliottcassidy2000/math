#!/usr/bin/env python3
"""Scratch census for the THM-2888/2893 heavy-triangle frontier."""

from __future__ import annotations

import importlib.util
import multiprocessing as mp
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PATH = (
    ROOT
    / "04-computation/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.py"
)


def load():
    spec = importlib.util.spec_from_file_location("thm2888_stats", PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import THM-2888")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load()


def profile(body):
    return T.profile_body(body)


def main() -> None:
    with mp.get_context("fork").Pool(6) as pool:
        nested = pool.map(profile, T.BODIES, chunksize=1)
    rows = [row for block in nested for row in block]
    by_key = {(row["body"], row["apex"]): row for row in rows}

    hostile_keys = {
        (body, apex)
        for body, _, apex, _, _, _ in T.EXPECTED_HOSTILES
    }
    normalization_sources = [
        row
        for row in rows
        if row["scalar_class"] == "failure"
        and (
            row["margin"] <= 0
            or (
                row["cutoff"] is not None
                and row["cutoff"] > T.PAIR_HORIZON
            )
        )
    ]
    normalized = {
        (row["body"], row["apex"]): T.normalize_heavy_edge(row)
        for row in normalization_sources
    }

    residual_bodies = []
    residual_keys = []
    for body_index, body in enumerate(T.BODIES):
        body_rows = rows[10 * body_index : 10 * body_index + 10]
        root = T.THM2885.profile_body(body)
        terminals = {
            row["apex"]
            for row in body_rows
            if row["scalar_class"] == "direct"
            or row["pairpair_direct"]
            or (body, row["apex"]) in hostile_keys
        }
        allowed = [
            (value, speed)
            for rank, (value, speed) in enumerate(root["top15"], start=1)
            if rank >= 11 or speed not in terminals
        ]
        margin = root["m"] - sum(value for value, _ in allowed[:5])
        if max(body) >= 13 and margin <= 0:
            residual_bodies.append(body)
            residual_keys.extend(
                (body, row["apex"])
                for row in body_rows
                if row["apex"] not in terminals
            )

    if len(residual_bodies) != 1444 or len(residual_keys) != 10202:
        raise RuntimeError("THM-2888 residual universe changed")

    core_rows = []
    for index, key in enumerate(residual_keys, start=1):
        row = by_key[key]
        norm = normalized.get(key)
        cap = row["global_cap"] if norm is None else norm["deleted_cap"]
        cutoff = row["cutoff"] if norm is None else norm["deleted_cutoff"]
        forbidden = None if norm is None else norm["edge"]
        if cutoff is None:
            raise RuntimeError(f"missing cutoff {key}")
        root_good, _, _ = T.CORE.good_norm(row["body"])
        carrier = T.THM2883.subtract_local(root_good, row["apex"])
        threshold = (row["m"] - cap) / 2
        speeds = [
            speed
            for speed in range(T.FIRST_EXTERNAL, cutoff + 1)
            if speed != row["apex"]
        ]
        coverage_rows = T.THM2885.coverages_many(carrier, speeds)
        core = tuple(
            speed for value, speed in coverage_rows if value >= threshold
        )
        core_rows.append(
            (
                key,
                row["scalar_class"],
                cutoff,
                len(core),
                forbidden,
                core,
            )
        )
        if index % 500 == 0:
            print(f"progress={index}/{len(residual_keys)}", flush=True)

    size_counter = Counter(size for _, _, _, size, _, _ in core_rows)
    class_counter = Counter(kind for _, kind, _, _, _, _ in core_rows)
    cutoff_counter = Counter()
    for _, _, cutoff, _, _, _ in core_rows:
        for bound in (50, 100, 250, 500, 1000, 2500, 5000, 10000):
            if cutoff <= bound:
                cutoff_counter[bound] += 1
                break
        else:
            cutoff_counter[None] += 1

    print(f"residual_bodies={len(residual_bodies)}")
    print(f"residual_branches={len(residual_keys)}")
    print(f"class_counter={dict(sorted(class_counter.items()))}")
    print(f"normalized_branches={sum(key in normalized for key in residual_keys)}")
    print(f"maximum_cutoff={max(row[2] for row in core_rows)}")
    print(f"cutoff_buckets={dict(cutoff_counter)}")
    print(f"core_size_min={min(size_counter)}")
    print(f"core_size_max={max(size_counter)}")
    print(f"core_size_hist={dict(sorted(size_counter.items()))}")
    largest = sorted(
        core_rows,
        key=lambda item: (-item[3], item[0]),
    )[:20]
    for key, kind, cutoff, size, forbidden, core in largest:
        print(
            f"largest={key};class={kind};cutoff={cutoff};"
            f"core_size={size};forbidden={forbidden};core={core}"
        )


if __name__ == "__main__":
    main()
