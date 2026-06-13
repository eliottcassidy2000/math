#!/usr/bin/env python3
"""
lrc_time_burnside_fourier_s571.py

codex-2026-06-03 S571

Finite time-slot Burnside/Fourier atlas for LRC time modeling.

The user proposed discretising time to Z/N and applying Burnside to lonely
time slots.  This script makes that precise and records the correction:

  * For a fixed speed set, the lonely-slot subset X is usually NOT invariant
    under the full time-translation group Z/N.  Therefore |X/G| under the full
    group is not a valid Burnside quotient except in the empty/all cases.
  * The correct finite Burnside group is the stabilizer K of the whole lonely
    time word.  Then X is K-invariant and Burnside gives |X/K|=|X|/|K| because
    nontrivial translations act freely on time slots.
  * The representation-theoretic dual is the discrete Fourier transform of the
    lonely indicator.  K-invariance forces nonzero characters to lie in the
    annihilator of K.  Frequency concentration is the dual shadow of reset
    folding/resonance.

Grid:
  N = lcm(14, 12, 15, 16, 23) = 38640.

This includes the n=14 clock and the main S564/S570 witnessed denominators for
the audited n=14 rows.  It is a finite approximation, not a proof by itself.

Tournament Analysis:
  Vertices: speed-set time words.
  Pairwise observable: (lonely density, stabilizer size, Fourier entropy,
    top nonzero Fourier mass, minimal period).
  Switch/gauge: harder row wins if it has lower density, then larger reset
    folding (stabilizer), then lower Fourier entropy / larger top mode.
  Tie Hamiltonian path: displayed sample order.

Assumption challenge:
  Time vertices could be individual slots, runner states, pair-sum events,
  wall crossings, Fourier modes, subgroup periods, or speed-set time words.
  This script uses whole time words because that preserves the predicate
  "which times are lonely?" and exposes the legal Burnside stabilizer.  It
  destroys interval ownership and exact continuous endpoints; those remain the
  S564/S570 objects.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isclose, lcm

import numpy as np


THRESHOLD_N = 14
GRID_N = lcm(14, 12, 15, 16, 23)


@dataclass(frozen=True)
class Sample:
    label: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class TimeWordReport:
    sample: Sample
    lonely_count: int
    density: F
    minimal_period: int
    stabilizer_size: int
    word_orbit_size: int
    lonely_orbits_under_stabilizer: int
    naive_full_burnside_average: F
    full_action_valid: bool
    exact_period: int
    allowed_frequency_modulus: int
    fourier_support_count: int
    fourier_entropy: float
    top_nonzero_mass: float
    top_frequencies: tuple[tuple[int, float], ...]
    annihilator_violations: int


def normalize(speeds: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for s in speeds:
        g = gcd(g, abs(s))
    return tuple(sorted({abs(s) // g for s in speeds if s}))


def divisors(n: int) -> list[int]:
    out = []
    for d in range(1, int(n**0.5) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def mobius(n: int) -> int:
    x = n
    p = 2
    primes = 0
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            if e > 1:
                return 0
            primes += 1
        p += 1 if p == 2 else 2
    if x > 1:
        primes += 1
    return -1 if primes % 2 else 1


def lonely_word(speeds: tuple[int, ...], grid_n: int = GRID_N) -> np.ndarray:
    slots = np.arange(grid_n, dtype=np.int64)
    ok = np.ones(grid_n, dtype=bool)
    for v in speeds:
        residues = (int(v % grid_n) * slots) % grid_n
        dist = np.minimum(residues, grid_n - residues)
        ok &= dist * THRESHOLD_N >= grid_n
    return ok


def invariant_under_shift(word: np.ndarray, shift: int) -> bool:
    return bool(np.array_equal(word, np.roll(word, -shift)))


def minimal_period(word: np.ndarray) -> int:
    n = len(word)
    for p in divisors(n):
        if invariant_under_shift(word, p):
            return p
    return n


def exact_period_by_mobius(word: np.ndarray) -> int:
    n = len(word)
    fixed = {p: int(invariant_under_shift(word, p)) for p in divisors(n)}
    exact = []
    for p in divisors(n):
        val = sum(mobius(p // d) * fixed[d] for d in divisors(p) if p % d == 0)
        if val:
            exact.append(p)
    if len(exact) != 1:
        return -1
    return exact[0]


def fourier_summary(word: np.ndarray, stabilizer_size: int) -> tuple[int, float, float, tuple[tuple[int, float], ...], int]:
    values = word.astype(float)
    coeffs = np.fft.fft(values)
    power = np.abs(coeffs) ** 2
    total = float(power.sum())
    if total == 0.0:
        return 0, 0.0, 0.0, tuple(), 0

    probs = power / total
    support = np.where(probs > 1e-12)[0]
    entropy = float(-np.sum(probs[support] * np.log2(probs[support]))) if len(support) else 0.0
    nonzero = [(int(f), float(probs[f])) for f in support if f != 0]
    nonzero.sort(key=lambda x: x[1], reverse=True)
    top = tuple(nonzero[:5])
    top_mass = top[0][1] if top else 0.0

    # If the word is invariant under a subgroup K of size k, nonzero Fourier
    # coefficients must be supported on frequencies divisible by k.
    violations = sum(1 for f in support if f % stabilizer_size != 0)
    return len(support), entropy, top_mass, top, violations


def packet(n: int, scale: int, skip: int) -> tuple[int, ...]:
    return (1,) + tuple(scale * q for q in range(1, n) if q != skip)


def samples() -> list[Sample]:
    return [
        Sample("AP_wall", tuple(range(1, 14)), "maximally folded wall; closed n-clock witnesses"),
        Sample("V_star_wall", tuple(list(range(1, 12)) + [13, 24]), "sporadic wall row"),
        Sample("near_AP_apex", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14), "near AP with apex replacement"),
        Sample("S562_packet_n14", packet(14, 7, 6), "HYP-2073 n=14 residual packet"),
        Sample("S562_packet_n14_lift", packet(14, 14, 6), "dyadic lift of n=14 packet"),
        Sample("no_small_pinch_proxy", (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201), "THM-396 proxy"),
        Sample("random_low_resonance", (2, 5, 11, 17, 23, 31, 37, 41, 43, 47, 53, 59, 61), "low-resonance primitive sample"),
    ]


def analyze(sample: Sample) -> TimeWordReport:
    speeds = normalize(sample.speeds)
    word = lonely_word(speeds)
    count = int(word.sum())
    density = F(count, GRID_N)
    period = minimal_period(word)
    stab = GRID_N // period
    full_valid = stab == GRID_N or count in (0, GRID_N)
    support_count, entropy, top_mass, top, violations = fourier_summary(word, stab)
    return TimeWordReport(
        sample=Sample(sample.label, speeds, sample.note),
        lonely_count=count,
        density=density,
        minimal_period=period,
        stabilizer_size=stab,
        word_orbit_size=period,
        lonely_orbits_under_stabilizer=count // stab if stab else count,
        naive_full_burnside_average=density,
        full_action_valid=full_valid,
        exact_period=exact_period_by_mobius(word),
        allowed_frequency_modulus=stab,
        fourier_support_count=support_count,
        fourier_entropy=entropy,
        top_nonzero_mass=top_mass,
        top_frequencies=top,
        annihilator_violations=violations,
    )


def fmt_frac(x: F) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def tournament_fingerprint(reports: list[TimeWordReport]) -> dict[str, object]:
    def key(r: TimeWordReport) -> tuple[F, int, float, float, int]:
        # Harder = lower density, more stabilizer folding, lower entropy,
        # stronger top mode, shorter period.  Negate density by comparing
        # smaller as larger via -float for the secondary tuple only.
        return (-r.density, r.stabilizer_size, -r.fourier_entropy, r.top_nonzero_mass, -r.minimal_period)

    n = len(reports)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(reports):
        for j, right in enumerate(reports):
            if i == j:
                continue
            adj[i][j] = key(left) > key(right) or (key(left) == key(right) and i < j)

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1

    def reach(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v, edge in enumerate(adj[u]):
                if edge and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reach(u)
        comp = {v for v in remaining if v in ru and u in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "vertices": [r.sample.label for r in reports],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": 1 if c3 == 0 and len(set(scores)) == n else "not_counted",
    }


def main() -> None:
    print("S571 finite time Burnside/Fourier atlas")
    print("=" * 78)
    print(f"grid N={GRID_N}; threshold=1/{THRESHOLD_N}")
    print("Full Z/N Burnside on lonely slots is valid only when the lonely word is")
    print("translation-invariant.  Otherwise |X|/N is a density, not an orbit count.")
    print()

    reports = [analyze(sample) for sample in samples()]
    for r in reports:
        print(r.sample.label)
        print(f"  note: {r.sample.note}")
        print(f"  speeds={r.sample.speeds}")
        print(
            f"  lonely_slots={r.lonely_count}/{GRID_N} density={fmt_frac(r.density)} "
            f"({float(r.density):.6f})"
        )
        print(
            f"  time-word period={r.minimal_period}; stabilizer |K|={r.stabilizer_size}; "
            f"word orbit under full shifts={r.word_orbit_size}"
        )
        print(
            f"  Burnside: naive full average={fmt_frac(r.naive_full_burnside_average)} "
            f"valid_full_action={r.full_action_valid}; "
            f"lonely orbits under K={r.lonely_orbits_under_stabilizer}"
        )
        print(
            f"  Mobius exact period={r.exact_period}; Fourier support={r.fourier_support_count}; "
            f"allowed freq multiple={r.allowed_frequency_modulus}; "
            f"annihilator violations={r.annihilator_violations}"
        )
        top = ", ".join(f"{f}:{mass:.4f}" for f, mass in r.top_frequencies[:4]) or "-"
        print(
            f"  Fourier entropy={r.fourier_entropy:.4f}; "
            f"top nonzero mass={r.top_nonzero_mass:.6f}; top freqs={top}"
        )

    print("\nTournament Analysis")
    print("  vertices: speed-set time words")
    print("  observable: density, stabilizer size, Fourier entropy/top mode, period")
    print("  switch: harder = lower density, more folding, lower entropy, stronger mode")
    print(f"  fingerprints: {tournament_fingerprint(reports)}")

    print("\nSynthesis")
    print("  Raw primitive reset time is one lap, so it does not classify integer rows.")
    print("  The useful reset invariant is the finite time word: period, stabilizer,")
    print("  and Fourier character support.  Burnside becomes legal after replacing")
    print("  the full shift group by the stabilizer K of that word.")
    print("  Pair-sum and n-clock witnesses are visible as time-word frequencies;")
    print("  the next proof object should combine S570 witness-or-core with this")
    print("  Fourier/stabilizer ledger on owner-labelled time events.")


if __name__ == "__main__":
    main()
