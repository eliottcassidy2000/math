#!/usr/bin/env python3
"""Parallel exact ray/status probe of the projected k=3 height-275 atlas bank."""

from __future__ import annotations

import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_PATH = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
FIRST = 275


def load_ray():
    spec = importlib.util.spec_from_file_location("ray275", RAY_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(RAY_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIRST = FIRST
    return module


def atlas_bodies():
    bodies = []
    pattern = re.compile(r"^row=E=([0-9,]+);")
    for line in ATLAS_PATH.read_text().splitlines():
        if f";z1={FIRST};" not in line:
            continue
        match = pattern.match(line)
        if match is None:
            raise RuntimeError(line)
        bodies.append(tuple(map(int, match.group(1).split(","))))
    if len(bodies) != 10:
        raise RuntimeError(bodies)
    return tuple(bodies)


def evaluate(body):
    ray = load_ray()
    stream = ray.Stream(body)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
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
    bodies = atlas_bodies()
    with mp.get_context("spawn").Pool(len(bodies)) as pool:
        records = tuple(pool.map(evaluate, bodies))
    totals = [0, 0, 0, 0]
    for record in records:
        body, L, high, first_d, trials, checks, states, crude, status, residual, rows = record
        totals = [
            totals[index] + value
            for index, value in enumerate((states, crude, status, residual))
        ]
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
        for ds, state in rows:
            print("RES", body, ds, state, flush=True)
    print("TOTAL", totals)


if __name__ == "__main__":
    main()
