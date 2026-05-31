#!/usr/bin/env python3
"""
lrc_bruhat_tits_descent_s410.py

codex-2026-05-31 S410

Bruhat-Tits descent probe for the Lonely Runner endpoint program.

This is not a full use of buildings.  It extracts the local p-adic tree
shadow that is already present in THM-357, THM-360, and THM-367:

* an endpoint t is a rational boundary point;
* at each prime p | n, the p-adic depth of t is v_p(denominator(t));
* a speed divisible by n is a gate that translates the unit endpoint layer
  back toward the root, but the gate's own endpoints reappear deeper in the
  p-adic tree.

The useful exact identity is:

    if v = n*q and p | n, then every endpoint of a v-interval has
    p-depth 2*v_p(n) + v_p(q),

because endpoints are (n*m +/- 1)/(n*v) and n*m +/- 1 is a p-adic unit.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S391 = SourceFileLoader(
    "lrc_n16_dyadic_endpoint_formula_s391",
    str(ROOT / "04-computation" / "lrc_n16_dyadic_endpoint_formula_s391.py"),
).load_module()


ONE = Fraction(1, 1)


@dataclass(frozen=True)
class PeelDepthLayer:
    step: int
    dead_endpoints: int
    dead_intervals: int
    remaining_endpoints: int
    remaining_intervals: int
    depth_hist: tuple[tuple[tuple[int, ...], int], ...]


def prime_factors(value: int) -> tuple[int, ...]:
    n = value
    out: list[int] = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.append(n)
    return tuple(out)


def vp(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("v_p(0) is not finite in this finite-depth ledger")
    value = abs(value)
    out = 0
    while value % prime == 0:
        out += 1
        value //= prime
    return out


def depth_vector(point: Fraction, primes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(vp(point.denominator, prime) for prime in primes)


def speed_vector(speed: int, primes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(vp(speed, prime) for prime in primes)


def fmt_frac(value: Fraction | None) -> str:
    return S360.fmt_frac(value)


def fmt_depth(depth: tuple[int, ...], primes: tuple[int, ...]) -> str:
    if not primes:
        return "{}"
    return "{" + ",".join(f"{prime}:{exp}" for prime, exp in zip(primes, depth)) + "}"


def fmt_hist(
    hist: Counter[tuple[int, ...]] | dict[tuple[int, ...], int],
    primes: tuple[int, ...],
    limit: int = 6,
) -> str:
    if not hist:
        return "empty"
    items = sorted(hist.items(), key=lambda item: (sum(item[0]), item[0], item[1]))
    parts = [f"{fmt_depth(depth, primes)}:{count}" for depth, count in items[:limit]]
    if len(items) > limit:
        parts.append("...")
    return " ".join(parts)


def endpoint_values(speeds: tuple[int, ...]) -> set[Fraction]:
    return {endpoint.value for endpoint in S360.endpoints(speeds)}


def unprotected_values(speeds: tuple[int, ...]) -> set[Fraction]:
    values = endpoint_values(speeds)
    return {
        value
        for value in values
        if not any(S360.direct_protects(speeds, protector, value) for protector in speeds)
    }


def depth_hist(values: set[Fraction], primes: tuple[int, ...]) -> Counter[tuple[int, ...]]:
    return Counter(depth_vector(value, primes) for value in values)


def unit_points(n: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(a, n) for a in range(1, n) if gcd(a, n) == 1)


def unit_protection_count(n: int, residue: int) -> int:
    threshold = Fraction(1, n)
    count = 0
    for point in unit_points(n):
        distance = S360.circular_distance_to_integer(residue * point)
        if distance < threshold:
            count += 1
    return count


def print_unit_gate_table(n: int) -> None:
    primes = prime_factors(n)
    gate = speed_vector(n, primes)
    units = unit_points(n)
    print(f"Unit layer gate table for n={n}")
    print(f"  primes={primes} gate_depth={fmt_depth(gate, primes)} unit_count={len(units)}")
    print("  residue classes that protect at least one unit endpoint:")
    rows = []
    for residue in range(n):
        count = unit_protection_count(n, residue)
        if count:
            rows.append((residue, count, speed_vector(residue or n, primes)))
    for residue, count, vector in rows:
        print(
            f"    r={residue:>2} protects={count:>2}/{len(units)} "
            f"local_vector={fmt_depth(vector, primes)}"
        )
    if rows == [(0, len(units), gate)]:
        print("  exact THM-360 shadow: the unit layer is killed only by an n-gate.")
    print()


def gate_export_depth(n: int, q: int, primes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(2 * vp(n, prime) + vp(q, prime) for prime in primes)


def print_gate_export_table(n: int, qmax: int = 10) -> None:
    primes = prime_factors(n)
    gate = speed_vector(n, primes)
    print(f"Gate endpoint export depths for n={n}")
    print(
        "  If v=n*q, every endpoint (n*m+/-1)/(n*v) has "
        "depth 2*gate_depth + depth(q)."
    )
    print(f"  gate_depth={fmt_depth(gate, primes)}")
    for q in range(1, qmax + 1):
        print(
            f"    q={q:<2} speed={n*q:<4} exported_depth="
            f"{fmt_depth(gate_export_depth(n, q, primes), primes)}"
        )
    print()


def peel_depth_layers(
    speeds: tuple[int, ...], primes: tuple[int, ...], max_layers: int = 4
) -> tuple[PeelDepthLayer, ...]:
    endpoints, intervals, _owners, protectors, boundary = S362.build_endpoint_system(
        speeds
    )
    remaining_endpoints = set(endpoints)
    remaining_intervals = set(intervals)
    layers: list[PeelDepthLayer] = []
    step = 0
    while step < max_layers:
        dead_endpoints = {
            endpoint
            for endpoint in remaining_endpoints
            if not (protectors[endpoint] & remaining_intervals)
        }
        if not dead_endpoints:
            break
        dead_intervals = {
            interval
            for interval in remaining_intervals
            if boundary[interval] & dead_endpoints
        }
        next_endpoints = remaining_endpoints - dead_endpoints
        next_intervals = remaining_intervals - dead_intervals
        hist = tuple(sorted(depth_hist(dead_endpoints, primes).items()))
        layers.append(
            PeelDepthLayer(
                step=step,
                dead_endpoints=len(dead_endpoints),
                dead_intervals=len(dead_intervals),
                remaining_endpoints=len(next_endpoints),
                remaining_intervals=len(next_intervals),
                depth_hist=hist,
            )
        )
        remaining_endpoints = next_endpoints
        remaining_intervals = next_intervals
        step += 1
    return tuple(layers)


def print_candidate(label: str, raw_speeds: tuple[int, ...]) -> None:
    summary = S360.summarize(list(raw_speeds))
    speeds = summary.speeds
    n = len(speeds) + 1
    primes = prime_factors(n)
    gate = speed_vector(n, primes)
    all_endpoints = endpoint_values(speeds)
    unprotected = unprotected_values(speeds)
    speed_hist = Counter(speed_vector(speed, primes) for speed in speeds)
    n_gates = tuple(speed for speed in speeds if speed % n == 0)
    gap_ratio = summary.max_gap / summary.threshold
    print(f"[{label}]")
    print(f"  n={n} primes={primes} gate_depth={fmt_depth(gate, primes)}")
    print(f"  speeds={speeds}")
    print(
        f"  class={summary.classification} forbidden={fmt_frac(summary.forbidden_length)} "
        f"gap/th={float(gap_ratio):.6f} Q={summary.boundary_modulus}"
    )
    print(
        f"  n_gates={n_gates if n_gates else 'none'} "
        f"speed_depth_hist={fmt_hist(speed_hist, primes)}"
    )
    print(
        f"  endpoint_depths={fmt_hist(depth_hist(all_endpoints, primes), primes)} "
        f"unique={len(all_endpoints)}"
    )
    print(
        f"  exposed_depths={fmt_hist(depth_hist(unprotected, primes), primes)} "
        f"unprotected={len(unprotected)} first={fmt_frac(summary.first_unprotected)}"
    )
    layers = peel_depth_layers(speeds, primes)
    if layers:
        print("  peel-depth first layers:")
        for layer in layers:
            print(
                f"    step={layer.step} deadE={layer.dead_endpoints:<4} "
                f"deadI={layer.dead_intervals:<4} remE={layer.remaining_endpoints:<4} "
                f"depths={fmt_hist(dict(layer.depth_hist), primes)}"
            )
    else:
        print("  peel-depth first layers: terminal core already nonempty")
    print()


def print_pure_dyadic_bt_table() -> None:
    print("Pure dyadic n=16 radial table")
    print(
        "  For owner u=2^k, every endpoint has 2-adic depth 4+k. "
        "A protector p=2^j*q translates the ray upward by j edges."
    )
    for owner in (2, 4, 8, 16, 32, 64):
        k = owner.bit_length() - 1
        by_j: dict[int, Counter[int]] = defaultdict(Counter)
        for protector in range(1, owner):
            by_j[vp(protector, 2)][S391.pure_dyadic_formula(k, protector)] += 1
        lower_covered = len(S391.lower_union(owner))
        endpoint_count = len(S391.endpoint_labels(owner))
        pieces = []
        for j in sorted(by_j):
            hist = ",".join(f"{count}:{freq}" for count, freq in sorted(by_j[j].items()))
            pieces.append(f"j={j}[{hist}]")
        print(
            f"  owner={owner:<3} depth={{2:{4+k}}} lower_union="
            f"{lower_covered:>3}/{endpoint_count:<3} " + " ".join(pieces)
        )
    print(
        "  Read this as a BT-tree warning: local radial descent is real, but "
        "u>=16 has a stable nine-fan that can close one branch locally."
    )
    print()


def print_synthesis() -> None:
    print("Synthesis")
    print(
        "  Bruhat-Tits descent is a sharper language for the existing "
        "Bohr-boundary descent.  Unit endpoints are first-layer boundary rays "
        "at the primes dividing n.  THM-360 says only an n-gate moves all of "
        "those rays back to the root."
    )
    print(
        "  But an n-gate creates its own endpoints at depth 2*v_p(n)+v_p(q). "
        "So gate repair is not deletion; it is debt export into child vertices "
        "of the local p-adic tree."
    )
    print(
        "  The n=16 dyadic theorem is exactly the rank-one BT lab: diagonal "
        "translation by p=2^j controls endpoint counts, while odd residue "
        "classes decide which branches are hit.  The next proof target is a "
        "global product-tree potential whose boundary strictly increases "
        "under every attempted all-protected core."
    )


def main() -> None:
    print("LRC Bruhat-Tits descent probe (codex-2026-05-31 S410)")
    print()
    print_unit_gate_table(14)
    print_unit_gate_table(16)
    print_gate_export_table(14, qmax=8)
    print_gate_export_table(16, qmax=8)

    print("Candidate endpoint-depth ledgers")
    print_candidate("n=14 initial tight", tuple(range(1, 14)))
    print_candidate("n=14 replace 13 by 14", tuple(range(1, 13)) + (14,))
    print_candidate(
        "n=14 S380 14-multiple ladder",
        (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
    )
    print_candidate("n=16 initial tight", tuple(range(1, 16)))
    print_candidate("n=16 drop 15 add 16", tuple(range(1, 15)) + (16,))
    print_candidate(
        "n=16 odd units plus seven 16-gates",
        tuple(sorted(tuple(range(1, 16, 2)) + tuple(16 * q for q in range(1, 8)))),
    )

    print_pure_dyadic_bt_table()
    print_synthesis()


if __name__ == "__main__":
    main()
