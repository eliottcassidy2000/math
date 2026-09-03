#!/usr/bin/env python3
"""Deterministic sampled search for hostile half-turn cubic-depth rows.

This is a scout only.  It searches twelve odd tails at h=420 on a fixed
midpoint mesh of the anchor-safe core.  Any reported candidate must be
recomputed by the separate exact rational wall sweep before it is used.
"""

from __future__ import annotations

import argparse
import random
from math import gcd

import numpy as np


H = 420
G = np.array([1.0, 0.0, 0.0, 0.0, -1.0, -4.0, -10.0])
RICH = np.full((13, 13), np.inf)
for _a in range(13):
    for _b in range(13 - _a):
        RICH[_a, _b] = (
            4.0 / 5.0
            - (_a + _b) / 5.0
            + 7.0 / 15.0 * (_a * (_a - 1) / 2 + _b * (_b - 1) / 2)
            - 2.0 / 5.0 * _a * _b
            - 3.0
            / 5.0
            * (
                _a * (_a - 1) * (_a - 2) / 6
                + _b * (_b - 1) * (_b - 2) / 6
            )
            + 2.0
            / 15.0
            * (
                _a * (_a - 1) / 2 * _b
                + _a * _b * (_b - 1) / 2
            )
        )


def candidate_pool() -> tuple[int, ...]:
    values = set(range(1, 401, 2))
    offsets = (-41, -27, -15, -13, -11, -9, -7, -5, -3, -1,
               1, 3, 5, 7, 9, 11, 13, 15, 27, 41)
    for j in range(1, 18):
        for offset in offsets:
            value = 840 * j + offset
            if value > 0 and value % 2:
                values.add(value)
    # Exact denominator-completion and inherited hostile labels.
    values.update((9, 11, 13, 99, 117, 143, 1287, 9009, 5041, 5881))
    return tuple(sorted(values))


def mesh(points_per_component: int) -> np.ndarray:
    k = np.repeat(np.arange(H, dtype=np.float64), points_per_component)
    j = np.tile(np.arange(points_per_component, dtype=np.float64), H)
    x = 1.0 + 12.0 * (j + 0.5) / points_per_component
    return (14.0 * k + x) / (28.0 * H)


def states(pool: tuple[int, ...], times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    phase = np.remainder(np.asarray(pool, dtype=np.float64)[:, None] * times[None, :], 1.0)
    low = (np.minimum(phase, 1.0 - phase) < 1.0 / 14.0).astype(np.int8)
    high = (np.abs(phase - 0.5) < 1.0 / 14.0).astype(np.int8)
    assert not np.any(low & high)
    return low, high


def score(a: np.ndarray, b: np.ndarray, objective: str) -> float:
    values = G[np.minimum(a, b)] if objective == "cubic" else RICH[a, b]
    return float(values.mean() * 3.0 / 7.0)


def complete_denominators(row: tuple[int, ...]) -> bool:
    speeds = (2 * H,) + row
    return all(any(w % modulus == 0 for w in speeds) for modulus in range(2, 15))


def admissible(row: tuple[int, ...], require_complete: bool) -> bool:
    return gcd(2 * H, *row) == 1 and (
        not require_complete or complete_denominators(row)
    )


def local_search(
    pool: tuple[int, ...],
    low: np.ndarray,
    high: np.ndarray,
    rng: random.Random,
    restarts: int,
    require_complete: bool,
    objective: str,
) -> tuple[float, tuple[int, ...]]:
    best_global = (float("inf"), ())
    n_pool = len(pool)
    all_indices = tuple(range(n_pool))

    for restart in range(restarts):
        while True:
            chosen = set(rng.sample(all_indices, 12))
            row0 = tuple(pool[i] for i in chosen)
            if admissible(row0, require_complete):
                break

        a = low[list(chosen)].sum(axis=0)
        b = high[list(chosen)].sum(axis=0)
        current = score(a, b, objective)

        for _round in range(40):
            improved = False
            order = list(chosen)
            rng.shuffle(order)
            for removed in order:
                a0 = a - low[removed]
                b0 = b - high[removed]
                retained = chosen - {removed}
                candidates = [
                    i
                    for i in all_indices
                    if (i not in chosen or i == removed)
                    and admissible(
                        tuple(pool[j] for j in retained | {i}), require_complete
                    )
                ]
                candidate_array = np.asarray(candidates, dtype=np.int32)
                aa = a0[None, :] + low[candidate_array]
                bb = b0[None, :] + high[candidate_array]
                if objective == "cubic":
                    values = G[np.minimum(aa, bb)]
                else:
                    values = RICH[aa, bb]
                scores = values.mean(axis=1) * 3.0 / 7.0
                position = int(np.argmin(scores))
                proposed = int(candidate_array[position])
                proposed_score = float(scores[position])
                if proposed_score < current - 1e-12:
                    chosen.remove(removed)
                    chosen.add(proposed)
                    a = a0 + low[proposed]
                    b = b0 + high[proposed]
                    current = proposed_score
                    improved = True
                    break
            if not improved:
                break

        row = tuple(sorted(pool[i] for i in chosen))
        assert admissible(row, require_complete)
        candidate = (current, row)
        if candidate < best_global:
            best_global = candidate
            print(
                f"best restart={restart} sampled_{objective}={current:.12f} "
                + "speeds=" + ",".join(map(str, row)),
                flush=True,
            )
    return best_global


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=48)
    parser.add_argument("--restarts", type=int, default=24)
    parser.add_argument("--seed", type=int, default=611420)
    parser.add_argument("--unconstrained", action="store_true")
    parser.add_argument("--objective", choices=("cubic", "rich"), default="cubic")
    args = parser.parse_args()

    pool = candidate_pool()
    times = mesh(args.points)
    low, high = states(pool, times)
    print(f"h={H} pool={len(pool)} mesh={len(times)} seed={args.seed}")
    value, row = local_search(
        pool,
        low,
        high,
        random.Random(args.seed),
        args.restarts,
        not args.unconstrained,
        args.objective,
    )
    print(f"FINAL sampled_{args.objective}={value:.12f}")
    print("FINAL speeds=" + ",".join(map(str, row)))
    print("FINAL denominator_complete=" + str(complete_denominators(row)))
    print("FINAL primitive=" + str(gcd(2 * H, *row) == 1))


if __name__ == "__main__":
    main()
