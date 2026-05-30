#!/usr/bin/env python3
"""
Fejer/root-sign phase probe.

A circulant tournament on Z/pZ is a sign choice on cyclic root pairs
{d,-d}, d=1,...,(p-1)/2.  If sigma_d=+1 chooses d and sigma_d=-1 chooses -d,
then for every nonzero Fourier character k,

    lambda_k = sum_d cos(2*pi*k*d/p)
               + i sum_d sigma_d sin(2*pi*k*d/p)
             = -1/2 + i S_k,

so

    |lambda_k|^2 = 1/4 + S_k^2.

The interval tournament is the all-one root-sign chamber, and its spectrum is
the Fejer kernel.  This script makes that dictionary measurable and compares
root-sign phase features with H for p <= 13.
"""

from __future__ import annotations

import cmath
import itertools
import math
from dataclasses import dataclass
from typing import Iterable


SMALL_H_PRIMES = (7, 11, 13)
PROFILE_PRIMES = (7, 11, 13, 17, 19, 23, 29, 31)


@dataclass(frozen=True)
class Profile:
    p: int
    signs: tuple[int, ...]
    connection: tuple[int, ...]
    label: str
    top_fraction: float
    low_pair_fraction: float
    ipr: float
    fejer_alignment: float
    sign_changes: int
    chamber_bias: float
    additive_energy: int
    is_interval_orbit: bool
    is_paley: bool
    hamilton_paths: int | None = None


def connection_from_signs(p: int, signs: tuple[int, ...]) -> tuple[int, ...]:
    values: list[int] = []
    for d, sign in enumerate(signs, start=1):
        values.append(d if sign > 0 else p - d)
    return tuple(sorted(values))


def interval_connection(p: int) -> tuple[int, ...]:
    return tuple(range(1, (p - 1) // 2 + 1))


def paley_connection(p: int) -> tuple[int, ...] | None:
    if p % 4 != 3:
        return None
    return tuple(sorted({pow(a, 2, p) for a in range(1, p)}))


def all_signs(p: int) -> Iterable[tuple[int, ...]]:
    m = (p - 1) // 2
    return itertools.product((1, -1), repeat=m)


def eigen_magnitudes_direct(p: int, connection: tuple[int, ...]) -> list[float]:
    omega = cmath.exp(2j * math.pi / p)
    values = []
    for k in range(1, p):
        lam = sum(omega ** (k * s) for s in connection)
        values.append(abs(lam) ** 2)
    return values


def eigen_magnitudes_from_root_signs(p: int, signs: tuple[int, ...]) -> list[float]:
    values = []
    for k in range(1, p):
        sine_sum = 0.0
        for d, sign in enumerate(signs, start=1):
            sine_sum += sign * math.sin(2 * math.pi * k * d / p)
        values.append(0.25 + sine_sum * sine_sum)
    return values


def fejer_magnitudes(p: int) -> list[float]:
    m = (p - 1) // 2
    values = []
    for k in range(1, p):
        numerator = math.sin(math.pi * m * k / p) ** 2
        denominator = math.sin(math.pi * k / p) ** 2
        values.append(numerator / denominator)
    return values


def cosine_alignment(a: list[float], b: list[float]) -> float:
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    return dot / (na * nb)


def additive_energy(connection: tuple[int, ...], p: int) -> int:
    sset = set(connection)
    total = 0
    for a in connection:
        for b in connection:
            for c in connection:
                if (a + b - c) % p in sset:
                    total += 1
    return total


def is_unit_interval_orbit(connection: tuple[int, ...], p: int) -> bool:
    base = set(interval_connection(p))
    target = set(connection)
    for u in range(1, p):
        if {u * x % p for x in base} == target:
            return True
    return False


def classify(p: int, signs: tuple[int, ...], connection: tuple[int, ...]) -> str:
    paley = paley_connection(p)
    if set(connection) == set(interval_connection(p)):
        return "interval"
    if paley is not None and set(connection) == set(paley):
        return "paley"
    if is_unit_interval_orbit(connection, p):
        return "interval_orbit"
    return "other"


def out_masks(p: int, connection: tuple[int, ...]) -> list[int]:
    masks = [0] * p
    for v in range(p):
        mask = 0
        for s in connection:
            mask |= 1 << ((v + s) % p)
        masks[v] = mask
    return masks


def count_hamiltonian_paths(p: int, connection: tuple[int, ...]) -> int:
    outs = out_masks(p, connection)
    full = (1 << p) - 1
    dp = [[0] * p for _ in range(1 << p)]
    for v in range(p):
        dp[1 << v][v] = 1
    for mask in range(1 << p):
        row = dp[mask]
        for v, count in enumerate(row):
            if count == 0:
                continue
            moves = outs[v] & ~mask
            while moves:
                bit = moves & -moves
                w = bit.bit_length() - 1
                dp[mask | bit][w] += count
                moves ^= bit
    return sum(dp[full])


def sign_change_count(signs: tuple[int, ...]) -> int:
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def make_profile(p: int, signs: tuple[int, ...], include_h: bool) -> Profile:
    connection = connection_from_signs(p, signs)
    mags = eigen_magnitudes_from_root_signs(p, signs)
    fejer = fejer_magnitudes(p)
    total = sum(mags)
    low_pair = mags[0] + mags[-1]
    paley = paley_connection(p)
    h = count_hamiltonian_paths(p, connection) if include_h else None
    return Profile(
        p=p,
        signs=signs,
        connection=connection,
        label=classify(p, signs, connection),
        top_fraction=max(mags) / total,
        low_pair_fraction=low_pair / total,
        ipr=sum(x * x for x in mags) / (total * total),
        fejer_alignment=cosine_alignment(mags, fejer),
        sign_changes=sign_change_count(signs),
        chamber_bias=abs(sum(signs)) / len(signs),
        additive_energy=additive_energy(connection, p),
        is_interval_orbit=is_unit_interval_orbit(connection, p),
        is_paley=(paley is not None and set(connection) == set(paley)),
        hamilton_paths=h,
    )


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def fmt_connection(connection: tuple[int, ...], width: int = 24) -> str:
    text = "{" + ",".join(str(x) for x in connection) + "}"
    return text if len(text) <= width else text[: width - 3] + "..."


def verify_identity() -> None:
    print("=" * 72)
    print("IDENTITY: CIRCULANT ROOT SIGNS -> FOURIER/FEJER CHANNELS")
    print("=" * 72)
    for p in PROFILE_PRIMES:
        max_root_err = 0.0
        sample_count = 0
        for signs in itertools.islice(all_signs(p), 0, min(64, 1 << ((p - 1) // 2))):
            connection = connection_from_signs(p, signs)
            direct = eigen_magnitudes_direct(p, connection)
            root = eigen_magnitudes_from_root_signs(p, signs)
            max_root_err = max(max_root_err, max(abs(a - b) for a, b in zip(direct, root)))
            sample_count += 1

        interval_signs = (1,) * ((p - 1) // 2)
        interval_root = eigen_magnitudes_from_root_signs(p, interval_signs)
        fejer = fejer_magnitudes(p)
        fejer_err = max(abs(a - b) for a, b in zip(interval_root, fejer))
        print(
            f"p={p:2d} sampled_signs={sample_count:3d} "
            f"max_root_formula_err={max_root_err:.2e} "
            f"interval_fejer_err={fejer_err:.2e}"
        )


def summarize_profiles() -> dict[int, list[Profile]]:
    print("\n" + "=" * 72)
    print("PHASE PROFILE: INTERVAL CHAMBER VS PALEY FLATNESS")
    print("=" * 72)
    profiles_by_p: dict[int, list[Profile]] = {}
    for p in PROFILE_PRIMES:
        profiles = [make_profile(p, tuple(s), include_h=(p in SMALL_H_PRIMES)) for s in all_signs(p)]
        profiles_by_p[p] = profiles

        interval = next(x for x in profiles if x.label == "interval")
        best_top = max(profiles, key=lambda x: x.top_fraction)
        best_fejer = max(profiles, key=lambda x: x.fejer_alignment)
        paley = next((x for x in profiles if x.is_paley), None)

        print(f"\np={p}, root-pair signs={len(interval.signs)}, circulants={len(profiles)}")
        print(
            "  interval: "
            f"top={interval.top_fraction:.6f} low_pair={interval.low_pair_fraction:.6f} "
            f"ipr={interval.ipr:.6f} fejer={interval.fejer_alignment:.6f} "
            f"changes={interval.sign_changes} E={interval.additive_energy}"
        )
        if paley is not None:
            print(
                "  paley:    "
                f"top={paley.top_fraction:.6f} low_pair={paley.low_pair_fraction:.6f} "
                f"ipr={paley.ipr:.6f} fejer={paley.fejer_alignment:.6f} "
                f"changes={paley.sign_changes} E={paley.additive_energy}"
            )
        print(
            "  best top: "
            f"{best_top.label:14s} top={best_top.top_fraction:.6f} "
            f"S={fmt_connection(best_top.connection)}"
        )
        print(
            "  best Fejer alignment: "
            f"{best_fejer.label:14s} align={best_fejer.fejer_alignment:.6f} "
            f"S={fmt_connection(best_fejer.connection)}"
        )
    return profiles_by_p


def compare_with_h(profiles_by_p: dict[int, list[Profile]]) -> None:
    print("\n" + "=" * 72)
    print("H COMPARISON FOR SMALL PRIME CIRCULANTS")
    print("=" * 72)
    features = [
        ("top_fraction", lambda x: x.top_fraction),
        ("low_pair_fraction", lambda x: x.low_pair_fraction),
        ("ipr", lambda x: x.ipr),
        ("fejer_alignment", lambda x: x.fejer_alignment),
        ("sign_changes", lambda x: float(x.sign_changes)),
        ("chamber_bias", lambda x: x.chamber_bias),
        ("additive_energy", lambda x: float(x.additive_energy)),
    ]
    for p in SMALL_H_PRIMES:
        profiles = profiles_by_p[p]
        hs = [float(x.hamilton_paths) for x in profiles if x.hamilton_paths is not None]
        print(f"\np={p}:")
        for name, getter in features:
            xs = [getter(x) for x in profiles]
            print(f"  corr(H,{name}) = {pearson(xs, hs):+.6f}")

        ranked = sorted(profiles, key=lambda x: (-int(x.hamilton_paths or 0), -x.fejer_alignment))
        print("  top H rows:")
        print("    rank       H  label           top      fejer   changes  E      S")
        for rank, prof in enumerate(ranked[:8], start=1):
            print(
                f"    {rank:>4} {prof.hamilton_paths:>8}  {prof.label:14s} "
                f"{prof.top_fraction:>7.4f}  {prof.fejer_alignment:>7.4f} "
                f"{prof.sign_changes:>7} {prof.additive_energy:>4}  "
                f"{fmt_connection(prof.connection, width=22)}"
            )


def tangent_scan(profiles_by_p: dict[int, list[Profile]]) -> None:
    print("\n" + "=" * 72)
    print("NEW TANGENT CANDIDATES")
    print("=" * 72)
    print("1. Root-sign quotient:")
    print("   The type-A root sign cube collapses under cyclic translation to")
    print("   m=(p-1)/2 binary choices.  Fourier characters read this quotient")
    print("   through sine projections; Fejer is the all-one chamber shadow.")
    print()
    print("2. Fejer alignment is not simply H:")
    print("   At p=7 and p=11 Paley wins H while being anti-Fejer/flat; by p=13")
    print("   the interval orbit wins, but fixed Fejer alignment remains coordinate fragile.")
    print("   This marks a phase transition, not a monotone spectral theorem.")
    print()
    print("3. Interval orbits are the harmonic-analysis extremals:")
    for p in (7, 11, 13, 17):
        profiles = profiles_by_p[p]
        max_top = max(x.top_fraction for x in profiles)
        count = sum(1 for x in profiles if abs(x.top_fraction - max_top) < 1e-12)
        orbit_count = sum(
            1
            for x in profiles
            if abs(x.top_fraction - max_top) < 1e-12 and x.is_interval_orbit
        )
        print(f"   p={p}: top_fraction maximizers={count}, in interval orbit={orbit_count}")
    print()
    print("4. Root-sign roughness is a separate channel:")
    print("   Sign changes correlate weakly with H compared with phase features.")
    print("   The order on root pairs is chamber data; the character projection is")
    print("   the phase data.  They should be kept as separate feature blocks.")


def main() -> None:
    verify_identity()
    profiles_by_p = summarize_profiles()
    compare_with_h(profiles_by_p)
    tangent_scan(profiles_by_p)


if __name__ == "__main__":
    main()
