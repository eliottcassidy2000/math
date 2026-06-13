#!/usr/bin/env python3
"""
lrc_dual_burnside_constraints_s573.py

codex-2026-06-03 S573

Sharper dual Burnside constraints for LRC time words.

HYP-2085 corrected the naive cyclic Burnside count: the lonely subset is not
usually invariant under all translations, so one must use the cyclic stabilizer
of the whole lonely word.  HYP-2086 reframed LRC's open/boundary split as
Burnside orbit/fix duality.

This script tightens both:

  * Every LRC lonely word is invariant under time reversal t -> -t, because
    ||v t|| = ||-v t||.  Thus the legal finite symmetry is not just cyclic; it
    contains a dihedral reflection.
  * If the cyclic stabilizer has size k, the full legal dihedral stabilizer has
    size 2k.  Its reflection terms contribute actual fixed-time slots to
    Burnside, so the quotient is sharper than |X|/k when fixed slots are lonely.
  * On the dual side, reflection invariance kills the sine/odd sector.  In
    Fourier language, all nonzero energy lives in real cosine characters, and
    cyclic k-folding still forces frequencies to be multiples of k.

Tournament Analysis:
  Vertices: speed-set time words.
  Pairwise observable: (density, cyclic fold, dihedral orbit count, cosine
    energy defect, boundary fixed lonely slots, top character mass).
  Switch/gauge: harder row wins by lower density, more cyclic folding, fewer
    dihedral lonely orbits, smaller odd-sector defect, and stronger top mode.
  Tie Hamiltonian path: displayed sample order.

Assumption challenge:
  Vertices could be time slots, translations, reflections, Fourier modes,
  wall-crossing events, or labelled endpoint-core states.  This audit keeps
  whole binary time words.  It preserves lonely/non-lonely slots plus legal
  dihedral symmetry and destroys runner ownership, wall endpoint labels, and
  pair-sum/core provenance.  The challenged assumption is that cyclic shifts
  are the whole legal Burnside group; time reversal is always present.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm

import numpy as np


THRESHOLD_N = 14
GRID_N = lcm(14, 12, 15, 16, 23)


@dataclass(frozen=True)
class Sample:
    label: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class DihedralReport:
    sample: Sample
    lonely_count: int
    density: F
    cyclic_period: int
    cyclic_stabilizer: int
    dihedral_stabilizer: int
    reflection_axes: tuple[int, ...]
    reflection_fixed_lonely_total: int
    dihedral_lonely_orbits: F
    cyclic_annihilator_violations: int
    odd_sector_l2: float
    imag_l2: float
    top_cosine_frequencies: tuple[tuple[int, float], ...]


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


def reflection_fixed_lonely(word: np.ndarray, axis: int) -> int:
    n = len(word)
    # Fixed slots solve 2t = axis mod n.  Brute force is cheap enough here and
    # avoids parity mistakes when n is even.
    fixed = [t for t in range(n) if (2 * t - axis) % n == 0]
    return sum(1 for t in fixed if word[t])


def fourier_constraints(word: np.ndarray, cyclic_stabilizer: int) -> tuple[int, float, float, tuple[tuple[int, float], ...]]:
    values = word.astype(float)
    coeffs = np.fft.fft(values)
    power = np.abs(coeffs) ** 2
    total = float(power.sum())
    if total == 0.0:
        return 0, 0.0, 0.0, tuple()

    support = np.where(power / total > 1e-12)[0]
    annihilator_violations = sum(1 for f in support if f % cyclic_stabilizer != 0)

    imag_l2 = float(np.sum(np.imag(coeffs) ** 2) / total)

    # Reflection t -> -t makes the word even.  The odd part is the anti-fixed
    # projection under reflection; its L2 mass should vanish.
    reflected = values[(-np.arange(len(values))) % len(values)]
    odd = (values - reflected) / 2.0
    odd_l2 = float(np.dot(odd, odd) / max(1.0, np.dot(values, values)))

    probs = power / total
    nonzero = [(int(f), float(probs[f])) for f in support if f != 0]
    # Identify +/-f pairs by their smaller representative; cosine mass is the
    # combined mass of the pair.
    paired: dict[int, float] = {}
    n = len(word)
    for f, mass in nonzero:
        rep = min(f, (-f) % n)
        paired[rep] = paired.get(rep, 0.0) + mass
    top = tuple(sorted(paired.items(), key=lambda kv: kv[1], reverse=True)[:5])
    return annihilator_violations, odd_l2, imag_l2, top


def packet(n: int, scale: int, skip: int) -> tuple[int, ...]:
    return (1,) + tuple(scale * q for q in range(1, n) if q != skip)


def samples() -> list[Sample]:
    return [
        Sample("AP_wall", tuple(range(1, 14)), "maximally folded wall"),
        Sample("V_star_wall", tuple(list(range(1, 12)) + [13, 24]), "sporadic wall row"),
        Sample("near_AP_apex", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14), "near AP with apex"),
        Sample("S562_packet_n14", packet(14, 7, 6), "n=14 residual packet"),
        Sample("S562_packet_n14_lift", packet(14, 14, 6), "dyadic packet lift"),
        Sample("no_small_pinch_proxy", (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201), "THM-396 proxy"),
        Sample("random_low_resonance", (2, 5, 11, 17, 23, 31, 37, 41, 43, 47, 53, 59, 61), "low-resonance sample"),
        Sample("all_odd_probe", (1, 3, 5, 9, 11, 15, 17, 19, 23, 25, 29, 31, 33), "tests reflection fixed half-turn"),
    ]


def analyze(sample: Sample) -> DihedralReport:
    speeds = normalize(sample.speeds)
    word = lonely_word(speeds)
    count = int(word.sum())
    period = minimal_period(word)
    cyclic_stab = GRID_N // period
    # The base reflection r_0: t -> -t always preserves an LRC lonely word.
    # Composing r_0 with each cyclic stabilizer shift gives all reflection
    # symmetries; if another reflection existed, its product with r_0 would be
    # a translation stabilizer.
    axes = tuple((j * period) % GRID_N for j in range(cyclic_stab))
    fixed_total = sum(reflection_fixed_lonely(word, axis) for axis in axes)
    dihedral_stab = cyclic_stab + len(axes)
    dihedral_orbits = F(count + fixed_total, dihedral_stab)
    violations, odd_l2, imag_l2, top = fourier_constraints(word, cyclic_stab)
    return DihedralReport(
        sample=Sample(sample.label, speeds, sample.note),
        lonely_count=count,
        density=F(count, GRID_N),
        cyclic_period=period,
        cyclic_stabilizer=cyclic_stab,
        dihedral_stabilizer=dihedral_stab,
        reflection_axes=axes,
        reflection_fixed_lonely_total=fixed_total,
        dihedral_lonely_orbits=dihedral_orbits,
        cyclic_annihilator_violations=violations,
        odd_sector_l2=odd_l2,
        imag_l2=imag_l2,
        top_cosine_frequencies=top,
    )


def fmt_frac(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def tournament_fingerprint(reports: list[DihedralReport]) -> dict[str, object]:
    def key(r: DihedralReport) -> tuple[F, int, F, float, float]:
        return (
            -r.density,
            r.cyclic_stabilizer,
            -r.dihedral_lonely_orbits,
            -r.odd_sector_l2,
            r.top_cosine_frequencies[0][1] if r.top_cosine_frequencies else 0.0,
        )

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
    print("S573 sharpened dual Burnside constraints for LRC")
    print("=" * 78)
    print(f"grid N={GRID_N}; threshold=1/{THRESHOLD_N}")
    print("Universal constraint: 1_X(t)=1_X(-t).  The legal stabilizer is")
    print("dihedral: cyclic folds plus their reflection cosets.")
    print("Dual constraint: the lonely word has zero odd/sine sector.")
    print()

    reports = [analyze(sample) for sample in samples()]
    for r in reports:
        axes_preview = ",".join(str(a) for a in r.reflection_axes[:6])
        if len(r.reflection_axes) > 6:
            axes_preview += ",..."
        top = ", ".join(f"{f}:{mass:.4f}" for f, mass in r.top_cosine_frequencies[:4]) or "-"
        print(r.sample.label)
        print(f"  note: {r.sample.note}")
        print(f"  speeds={r.sample.speeds}")
        print(
            f"  lonely={r.lonely_count}/{GRID_N} density={fmt_frac(r.density)} "
            f"({float(r.density):.6f})"
        )
        print(
            f"  cyclic: period={r.cyclic_period}; |K_cyc|={r.cyclic_stabilizer}; "
            f"annihilator violations={r.cyclic_annihilator_violations}"
        )
        print(
            f"  dihedral: |K_D|={r.dihedral_stabilizer}; reflection axes={len(r.reflection_axes)} "
            f"[{axes_preview}]; reflection fixed lonely total={r.reflection_fixed_lonely_total}; "
            f"lonely orbits={fmt_frac(r.dihedral_lonely_orbits)}"
        )
        print(
            f"  dual: odd-sector L2={r.odd_sector_l2:.3e}; imag Fourier L2={r.imag_l2:.3e}; "
            f"top cosine pairs={top}"
        )

    print("\nTournament Analysis")
    print("  vertices: speed-set time words")
    print("  observable: density, cyclic fold, dihedral orbits, odd-sector defect, top cosine mode")
    print("  switch: harder = lower density, more fold, fewer quotient orbits, zero odd defect, stronger mode")
    print(f"  fingerprints: {tournament_fingerprint(reports)}")

    print("\nSharper constraints")
    print("  1. Cyclic Burnside is still too small: time reversal is always legal.")
    print("  2. A counterexample/time-boundary obstruction must live in the dihedral")
    print("     fixed sector, hence in the real cosine Fourier subspace.")
    print("  3. Nonzero sine/odd character mass is not merely small; it is forbidden.")
    print("  4. When cyclic stabilizer is trivial, the binary word still quotients by")
    print("     reflection, usually halving lonely slots into t,-t pairs.")
    print("  5. The next labelled-event audit should test which owner/core labels remain")
    print("     reflection-fixed and which labels break only after forgetting ownership.")


if __name__ == "__main__":
    main()
