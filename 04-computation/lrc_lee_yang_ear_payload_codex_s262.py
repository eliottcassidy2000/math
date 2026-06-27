#!/usr/bin/env python3
"""Lee-Yang / ear-payload scout for the LRC14 miss-count PGF.

The S66 signal studies the roots of

    G_E(z) = E[z^N],  N = number of empty 7-sectors after the runners in E.

This scout adds the extension payload for a one-runner "ear".  If F=E+{a} and
A_t = P(N_E=t and runner a lands in a sector that is empty for E), then

    q_F[t] = q_E[t] - A_t + A_{t+1}
    G_F(z) = G_E(z) + (z^-1 - 1) A_E,a(z).

Thus the hidden payload A, not the scalar p0 alone, is the root-motion carrier.
"""

from __future__ import annotations

import math
from collections import defaultdict
from fractions import Fraction as F

import numpy as np


def sector_of(x: F) -> int:
    return int((x % 1) * 7)


def breakpoints(E: tuple[int, ...]) -> list[F]:
    b = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * abs(e) + 1):
            b.add(F(m, 7 * abs(e)))
    return sorted(b)


def hit_set(E: tuple[int, ...], x: F) -> set[int]:
    return {sector_of(F(e) * x) for e in E}


def missdist(E: tuple[int, ...]) -> list[F]:
    q = [F(0)] * 7
    b = breakpoints(E)
    for x0, x1 in zip(b, b[1:]):
        if x1 <= x0:
            continue
        x = (x0 + x1) / 2
        t = 7 - len(hit_set(E, x))
        if 0 <= t <= 6:
            q[t] += x1 - x0
    return q


def ear_payload(base: tuple[int, ...], new: int) -> tuple[list[F], list[F]]:
    """Return (q_full, A), where A_t is the hidden empty-hit payload."""
    full = tuple(sorted(set(base + (new,))))
    q_full = missdist(full)
    A = [F(0)] * 7
    b = breakpoints(tuple(sorted(set(base + (new,)))))
    for x0, x1 in zip(b, b[1:]):
        if x1 <= x0:
            continue
        x = (x0 + x1) / 2
        before = hit_set(base, x)
        t = 7 - len(before)
        s = sector_of(F(new) * x)
        if 0 <= t <= 6 and s not in before:
            A[t] += x1 - x0
    return q_full, A


def roots(q: list[F]) -> list[complex]:
    coeffs = [float(q[t]) for t in range(6, -1, -1)]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
        coeffs = coeffs[1:]
    return sorted(np.roots(coeffs), key=lambda z: (abs(z.imag) < 1e-8, abs(z)))


def root_metrics(q: list[F]) -> dict[str, float | int]:
    r = roots(q)
    if not r:
        return {
            "degree": 0,
            "real_roots": 0,
            "nearest_modulus": float("inf"),
            "negative_interval_distance": float("inf"),
            "axis_gap_deg": 180.0,
            "min_modulus": float("inf"),
            "max_modulus": 0.0,
        }
    real_roots = sum(1 for z in r if abs(z.imag) < 1e-8)
    nearest_modulus = min(abs(z) for z in r)
    min_modulus = min(abs(z) for z in r)
    max_modulus = max(abs(z) for z in r)
    distances = []
    axis_gaps = []
    for z in r:
        # Distance to the closed Lee-Yang danger interval [-1,0].
        x = min(0.0, max(-1.0, z.real))
        distances.append(abs(z - x))
        theta = abs(math.degrees(math.atan2(z.imag, z.real)))
        axis_gaps.append(min(theta, 180.0 - theta))
    return {
        "degree": len(r),
        "real_roots": real_roots,
        "nearest_modulus": nearest_modulus,
        "negative_interval_distance": min(distances),
        "axis_gap_deg": min(axis_gaps),
        "min_modulus": min_modulus,
        "max_modulus": max_modulus,
    }


def fmt_frac(x: F) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_float(x: F | float) -> str:
    return f"{float(x):.6f}"


def payload_stats(A: list[F]) -> dict[str, F | int | str]:
    total = sum(A)
    if total == 0:
        return {
            "total": F(0),
            "mean_level": F(0),
            "support": "-",
            "even_mass": F(0),
            "odd_mass": F(0),
            "edge_mass": F(0),
        }
    mean_level = sum(F(t) * A[t] for t in range(7)) / total
    support = ",".join(str(t) for t, a in enumerate(A) if a)
    even_mass = sum(A[t] for t in range(0, 7, 2))
    odd_mass = sum(A[t] for t in range(1, 7, 2))
    edge_mass = A[0] + A[6]
    return {
        "total": total,
        "mean_level": mean_level,
        "support": support,
        "even_mass": even_mass,
        "odd_mass": odd_mass,
        "edge_mass": edge_mass,
    }


def eval_pgf(q: list[F], z: float) -> float:
    return sum(float(q[t]) * (z**t) for t in range(7))


def print_config(name: str, E: tuple[int, ...]) -> None:
    q = missdist(E)
    m = root_metrics(q)
    print(f"{name:18s} E={E}")
    print(
        "  q="
        + str([fmt_frac(x) for x in q])
        + f"  p0={fmt_float(q[0])}  q6={fmt_float(q[6])}  extreme={fmt_float(q[0]+q[6])}"
    )
    print(
        f"  roots: real={m['real_roots']}/{m['degree']}  "
        f"nearest={m['nearest_modulus']:.4f}  "
        f"dist([-1,0])={m['negative_interval_distance']:.4f}  "
        f"axis_gap={m['axis_gap_deg']:.2f}deg  "
        f"mod_span={m['min_modulus']:.3f}..{m['max_modulus']:.3f}"
    )


def rank_curve(configs: dict[str, tuple[int, ...]]) -> None:
    qs = {name: missdist(E) for name, E in configs.items()}
    wins: dict[str, int] = defaultdict(int)
    sample_z = [i / 20 for i in range(0, 61)]  # 0..3
    flips = []
    last = None
    for z in sample_z:
        order = sorted(qs, key=lambda name: -eval_pgf(qs[name], z))
        wins[order[0]] += 1
        if order[0] != last:
            flips.append((z, order[0]))
            last = order[0]
    print("FUGACITY WINNER PROFILE on z in [0,3] sampled by 1/20:")
    for name, count in sorted(wins.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {name:18s} wins {count:2d}/61 grid points")
    print("  winner changes: " + ", ".join(f"z={z:.2f}->{name}" for z, name in flips))


def print_ear_chain(name: str, chain: tuple[int, ...]) -> None:
    print(f"\nEAR CHAIN: {name}")
    base: tuple[int, ...] = tuple()
    for new in chain:
        q_base = missdist(base) if base else [F(0)] * 6 + [F(1)]
        q_full, A = ear_payload(base, new)
        stats = payload_stats(A)
        m = root_metrics(q_full)
        reconstructed = [F(0)] * 7
        for t in range(7):
            reconstructed[t] = q_base[t] - A[t] + (A[t + 1] if t < 6 else F(0))
        ok = reconstructed == q_full
        print(
            f"  +{new:2d}: |E|={len(set(base+(new,))):2d} "
            f"p0={fmt_float(q_full[0])} real={m['real_roots']}/{m['degree']} "
            f"near={m['nearest_modulus']:.3f} d[-1,0]={m['negative_interval_distance']:.3f} "
            f"A_total={fmt_float(stats['total'])} A_mean={fmt_float(stats['mean_level'])} "
            f"A_support={stats['support']} parity={fmt_float(stats['even_mass'])}/{fmt_float(stats['odd_mass'])} "
            f"recon={'ok' if ok else 'FAIL'}"
        )
        base = tuple(sorted(set(base + (new,))))


def main() -> None:
    configs = {
        "consec_8": tuple(range(8)),
        "even_AP": tuple(2 * i for i in range(8)),
        "top_cluster": (0, 7, 8, 9, 10, 11, 12, 13),
        "break_mid": (0, 1, 5, 7, 8, 9, 11, 13),
        "single_far_21": (0, 1, 2, 3, 4, 5, 6, 21),
        "random_spread": (0, 1, 3, 7, 12, 20, 33, 54),
    }
    print("=" * 96)
    print("LEE-YANG ROOT METRICS FOR SELECTED LRC MISS-COUNT PGFs")
    print("=" * 96)
    for name, E in configs.items():
        print_config(name, E)

    print("\n" + "=" * 96)
    rank_curve(configs)

    print("\n" + "=" * 96)
    print("ONE-RUNNER EAR PAYLOADS")
    print("=" * 96)
    print_ear_chain("nested AP/consec growth", tuple(range(8)))
    print_ear_chain("top-cluster growth", (0, 7, 8, 9, 10, 11, 12, 13))
    print_ear_chain("middle-break growth", (0, 1, 5, 7, 8, 9, 11, 13))
    print_ear_chain("single-far resonance", (0, 1, 2, 3, 4, 5, 6, 21))

    print("\n" + "=" * 96)
    print("NEW SIGNALS SUGGESTED BY THE SCOUT")
    print("=" * 96)
    print("  lee_yang_negative_interval_distance: min distance from a root of G_N to [-1,0].")
    print("  root_axis_gap_deg: minimum angular clearance from the real axis.")
    print("  ear_empty_hit_payload_A: A_t = P(N_base=t and new runner hits an empty sector).")
    print("  ear_payload_mean_level: average miss-count level on which the new ear acts.")
    print("  ear_payload_parity_bias: even/odd mass split of A_t, the odd-ear analogue.")
    print("  root_motion_reconstruction_check: q_full[t] = q_base[t] - A_t + A_{t+1}.")
    print("  fugacity_winner_profile: which configuration maximizes G_N(z) along z in [0,3].")


if __name__ == "__main__":
    main()
