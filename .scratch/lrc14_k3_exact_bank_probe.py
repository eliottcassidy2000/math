#!/usr/bin/env python3
"""Generic parallel exact ray/status probe for one projected k=3 atlas height."""

from __future__ import annotations

import argparse
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_PATH = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)


def load_ray(first):
    spec = importlib.util.spec_from_file_location(
        f"ray_bank_{first}",
        RAY_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(RAY_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIRST = first
    return module


def atlas_bodies(first):
    pattern = re.compile(r"^row=E=([0-9,]+);")
    bodies = []
    for line in ATLAS_PATH.read_text().splitlines():
        if f";z1={first};" not in line:
            continue
        match = pattern.match(line)
        if match is None:
            raise RuntimeError(line)
        bodies.append(tuple(map(int, match.group(1).split(","))))
    if not bodies:
        raise RuntimeError(("empty atlas height", first))
    return tuple(bodies)


def evaluate(task):
    first, body = task
    ray = load_ray(first)
    stream = ray.Stream(body)
    trials, states, checks, _signs = ray.ray_quotient_states(stream)
    crude, status, residual = ray.common_status_screen(stream, states)
    return (
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        trials,
        checks,
        len(states),
        len(crude),
        len(status),
        len(residual),
        tuple((ds, states[ds]) for ds in sorted(residual)),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, required=True)
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    bodies = atlas_bodies(args.first)
    tasks = tuple((args.first, body) for body in bodies)
    totals = [0, 0, 0, 0]
    residual_summaries = []
    with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
        for record in pool.imap_unordered(evaluate, tasks, chunksize=1):
            body, L, high, first_d, trials, checks, states, crude, status, residual, rows = record
            totals = [
                totals[index] + value
                for index, value in enumerate((states, crude, status, residual))
            ]
            if residual:
                residual_summaries.append(
                    (body, states, crude, status, residual)
                )
            print(
                "BODY",
                body,
                "L",
                L,
                "high",
                high,
                "d1",
                first_d,
                "counts",
                states,
                crude,
                status,
                residual,
                "trials",
                trials,
                "checks",
                checks,
                flush=True,
            )
            if not args.compact:
                for ds, state in rows:
                    print("RES", body, ds, state, flush=True)
    print("RESIDUAL_SUMMARIES", tuple(sorted(residual_summaries)))
    print("TOTAL", totals)


if __name__ == "__main__":
    main()
