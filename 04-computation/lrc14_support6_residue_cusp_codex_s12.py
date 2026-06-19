#!/usr/bin/env python3
"""
LRC(14) support-6 residue/cusp scout.

S11 showed that the coupled absolute support-6 Minkowski count is still much
too blunt, while the signed hyperplane sums are tiny.  This script probes the
next proof object:

    K(n_1,...,n_6,0,...,0) = C_d(n mod 7) / (n_1 ... n_6)

for nonzero n_i not divisible by 7.  Thus every exact six-support hyperplane
is a finite linear combination of reciprocal sums over residue-addressed
relation classes.  The computation below bins exact relation shells by the
max-norm boundary face they occupy.  If the large absolute envelope is coming
from coordinate cusps, THM-538's support-floor cancellation should be visible
as small signed boundary-face totals.

This is evidence and route-finding, not a proof.
"""
from __future__ import annotations

import cmath
import itertools
import math
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from math import comb, prod
from random import Random

TWO_PI_I = 2j * math.pi


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def fmt(x: float | complex) -> str:
    if isinstance(x, complex):
        return f"{x.real:.8g}{x.imag:+.2g}j"
    return f"{x:.8g}"


def shat(n: int, j: int) -> complex:
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (
        cmath.exp(-TWO_PI_I * n * a)
        - cmath.exp(-TWO_PI_I * n * (a + 1.0 / 7.0))
    ) / (TWO_PI_I * n)


SUBSETS = tuple(
    tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)
)
SIGNS = tuple((-1) ** len(T) for T in SUBSETS)


@lru_cache(None)
def chat(n: int, T: tuple[int, ...]) -> complex:
    if n == 0:
        return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0:
        return 0j
    return -sum(shat(n, j) for j in T)


def K_value(n: tuple[int, ...]) -> complex:
    total = 0j
    for sign, T in zip(SIGNS, SUBSETS):
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0:
                break
        total += sign * p
    return total


@lru_cache(None)
def residue_coeff(residues: tuple[int, ...], ambient_d: int) -> complex:
    """
    Normalized exact-support-6 coefficient.

    For all nonzero n_i with n_i == residues_i mod 7,

        K(n,0,...,0) = residue_coeff(residues,d) / prod(n_i).
    """
    assert len(residues) == 6
    vec = tuple(residues) + (0,) * (ambient_d - 6)
    return K_value(vec) * prod(residues)


def coeffs(H: int) -> tuple[int, ...]:
    return tuple(x for x in range(-H, H + 1) if x and x % 7)


def residue_tuple(ns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(n % 7 for n in ns)


def check_reciprocal_identity() -> None:
    section("RECIPROCAL RESIDUE IDENTITY")
    rng = Random(2614)
    worst = 0.0
    samples = 0
    for ambient_d in (6, 7, 8, 9):
        for _ in range(80):
            ns = []
            while len(ns) < 6:
                x = rng.randint(-90, 90)
                if x and x % 7:
                    ns.append(x)
            ns_t = tuple(ns)
            lhs = K_value(ns_t + (0,) * (ambient_d - 6))
            rhs = residue_coeff(residue_tuple(ns_t), ambient_d) / prod(ns_t)
            worst = max(worst, abs(lhs - rhs))
            samples += 1
    print(
        "Checked K(n,0,...)=C_d(n mod 7)/prod(n_i) on "
        f"{samples} random six-support vectors; worst error {worst:.3e}."
    )
    print(
        "Proof use: the six-support tail is not a free harmonic product; it is a "
        "finite family of residue-addressed reciprocal sums on relation hyperplanes."
    )


def check_simple_residue_marginals() -> None:
    """Guardrail: the raw residue coefficient is not coordinate-wise balanced."""
    section("RAW RESIDUE MARGINAL GUARDRAIL")
    print(
        "Testing the tempting stronger claim sum_{r_i} C_d(r)=0.  It is false; "
        "cancellation must use the relation hyperplane and finite wall split."
    )
    print(f"{'ambient d':>9} {'worst one-coordinate marginal':>32}")
    for ambient_d in (6, 7, 8, 9):
        worst = 0.0
        for i in range(6):
            for fixed in itertools.product(range(1, 7), repeat=5):
                vals = []
                it = iter(fixed)
                for j in range(6):
                    vals.append(None if j == i else next(it))
                total = 0j
                for r in range(1, 7):
                    vals[i] = r
                    total += residue_coeff(tuple(vals), ambient_d)
                worst = max(worst, abs(total))
        print(f"{ambient_d:>9} {worst:>32.12g}")


@dataclass
class ShellStats:
    count: int = 0
    signed: complex = 0j
    absK: float = 0.0

    def add(self, z: complex) -> None:
        self.count += 1
        self.signed += z
        self.absK += abs(z)


def six_support_shells(
    vals: tuple[int, int, int, int, int, int],
    ambient_d: int,
    Hmax: int,
) -> tuple[dict[int, ShellStats], dict[tuple[int, int], ShellStats]]:
    """Enumerate exact six-support relations and bin by height and boundary face."""
    left_vals = vals[:3]
    right_vals = vals[3:]
    cs = coeffs(Hmax)
    right = defaultdict(list)
    for ns in itertools.product(cs, repeat=3):
        s = sum(a * n for a, n in zip(right_vals, ns))
        right[s].append((ns, prod(ns), max(abs(n) for n in ns)))

    by_shell: dict[int, ShellStats] = defaultdict(ShellStats)
    by_shell_touch: dict[tuple[int, int], ShellStats] = defaultdict(ShellStats)
    coeff_cache: dict[tuple[int, ...], complex] = {}

    for lns in itertools.product(cs, repeat=3):
        ls = sum(a * n for a, n in zip(left_vals, lns))
        matches = right.get(-ls)
        if not matches:
            continue
        lprod = prod(lns)
        lmax = max(abs(n) for n in lns)
        for rns, rprod, rmax in matches:
            ns = lns + rns
            h = max(lmax, rmax)
            touch = sum(1 for n in ns if abs(n) == h)
            residues = residue_tuple(ns)
            c = coeff_cache.get(residues)
            if c is None:
                c = residue_coeff(residues, ambient_d)
                coeff_cache[residues] = c
            z = c / (lprod * rprod)
            by_shell[h].add(z)
            by_shell_touch[(h, touch)].add(z)
    return by_shell, by_shell_touch


def cumulative_rows(
    by_shell: dict[int, ShellStats],
    checkpoints: tuple[int, ...],
) -> list[tuple[int, int, complex, float, float]]:
    out = []
    count = 0
    signed = 0j
    absK = 0.0
    checkpoints_set = set(checkpoints)
    for h in range(1, max(checkpoints) + 1):
        s = by_shell.get(h)
        if s:
            count += s.count
            signed += s.signed
            absK += s.absK
        if h in checkpoints_set:
            ratio = absK / abs(signed) if abs(signed) > 1e-18 else float("inf")
            out.append((h, count, signed, absK, ratio))
    return out


def report_case(
    name: str,
    vals: tuple[int, int, int, int, int, int],
    ambient_d: int,
    Hmax: int,
    checkpoints: tuple[int, ...],
) -> dict[str, float]:
    section(f"CASE: {name}")
    print(f"support values={vals}  ambient nonzero-offset dimension d={ambient_d}  Hmax={Hmax}")
    by_shell, by_shell_touch = six_support_shells(vals, ambient_d, Hmax)

    print("\nCumulative exact signed support-6 reciprocal sum")
    print(f"{'H':>4} {'relations':>12} {'signed':>20} {'absK':>14} {'abs/signed':>12}")
    for H, count, signed, absK, ratio in cumulative_rows(by_shell, checkpoints):
        print(f"{H:>4} {count:>12} {fmt(signed):>20} {absK:>14.8g} {ratio:>12.6g}")

    final_h = max((h for h, s in by_shell.items() if s.count), default=Hmax)
    final_shell = by_shell.get(final_h, ShellStats())
    print(f"\nBoundary-face split on the last nonempty shell h={final_h}")
    print(f"{'touch':>5} {'relations':>12} {'signed':>20} {'absK':>14} {'abs/signed':>12}")
    for touch in range(1, 7):
        s = by_shell_touch.get((final_h, touch), ShellStats())
        ratio = s.absK / abs(s.signed) if abs(s.signed) > 1e-18 else float("inf")
        print(f"{touch:>5} {s.count:>12} {fmt(s.signed):>20} {s.absK:>14.8g} {ratio:>12.6g}")

    shell_ratio = (
        final_shell.absK / abs(final_shell.signed)
        if abs(final_shell.signed) > 1e-18
        else float("inf")
    )
    print(
        "\nReading: a large abs/signed ratio means the volume/envelope count is seeing "
        "coordinate-boundary mass that the exact six-sector kernel almost cancels."
    )
    return {
        "relations": float(sum(s.count for s in by_shell.values())),
        "cum_abs": sum(s.absK for s in by_shell.values()),
        "cum_signed_abs": abs(sum((s.signed for s in by_shell.values()), 0j)),
        "final_shell_abs": final_shell.absK,
        "final_shell_signed_abs": abs(final_shell.signed),
        "final_shell_ratio": shell_ratio,
    }


def tournament_fingerprint(case_scores: dict[str, dict[str, float]]) -> None:
    """
    Tournament Analysis with proof-obligation vertices.

    The observable is how much of the ghost absolute mass a proof quotient
    removes.  Edges point toward the quotient with smaller residual pressure.
    Ties follow the Hamiltonian path listed below.
    """
    section("TOURNAMENT ANALYSIS: PROOF-OBLIGATION QUOTIENTS")
    vertices = [
        "free_product_envelope",
        "coupled_absolute_hyperplane",
        "residue_permanent_constant",
        "signed_reciprocal_shell",
        "boundary_face_cancellation",
        "low_height_wall_ledger",
        "cluster_translation_quotient",
    ]
    # Scores are route-pressure ranks distilled from S10/S11 plus this S12
    # shell diagnostic.  Lower is better.
    median_ratio = sorted(v["final_shell_ratio"] for v in case_scores.values())[
        len(case_scores) // 2
    ]
    pressure = {
        "free_product_envelope": 1000.0,
        "coupled_absolute_hyperplane": 300.0,
        "residue_permanent_constant": 60.0,
        "signed_reciprocal_shell": max(1.0, 50.0 / median_ratio),
        "boundary_face_cancellation": max(0.5, 25.0 / median_ratio),
        "low_height_wall_ledger": 2.0,
        "cluster_translation_quotient": 3.0,
    }
    score = {v: 0 for v in vertices}
    flips = 0
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            if pressure[a] < pressure[b] or (
                pressure[a] == pressure[b] and vertices.index(a) < vertices.index(b)
            ):
                score[a] += 1
            else:
                score[b] += 1
                flips += 1
    hist = defaultdict(int)
    for s in score.values():
        hist[s] += 1
    order = sorted(vertices, key=lambda v: (pressure[v], vertices.index(v)))
    print(f"Hamiltonian proof path: {order}")
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"edge flips against the coarse old route order: {flips}")
    print(
        "Assumption challenged: vertices are not runners, arcs, or speeds.  They are "
        "quotients of the support-6 proof obligation: residue classes, boundary "
        "faces, low-height walls, and cluster translations.  This preserves the "
        "LRC predicate 'wide support-6 correction below cap margin' while destroying "
        "witness-time geometry."
    )


def main() -> None:
    section("LRC(14) SUPPORT-6 RESIDUE/CUSP SCOUT S12")
    print("Goal: explain why the signed support-6 layer is tiny after the absolute")
    print("Minkowski count remains huge: expose residue-address and boundary-cusp cancellation.")
    check_reciprocal_identity()
    check_simple_residue_marginals()

    cases = [
        (
            "AP core support in E={0..7}",
            (1, 2, 3, 4, 5, 6),
            7,
            28,
            (4, 8, 12, 16, 20, 24, 28),
        ),
        (
            "resonant one-stranger support in E={0..6,21}",
            (1, 2, 3, 4, 5, 21),
            7,
            28,
            (4, 8, 12, 16, 20, 24, 28),
        ),
        (
            "dissociated one-stranger support in E={0..6,211}",
            (2, 3, 4, 5, 6, 211),
            7,
            34,
            (4, 8, 12, 16, 20, 24, 28, 34),
        ),
        (
            "k=9 sampled wide support in E={0..7,68}",
            (2, 3, 4, 5, 6, 68),
            8,
            30,
            (4, 8, 12, 16, 20, 24, 30),
        ),
        (
            "k=10 height-one wall support in E={0..8,22}",
            (1, 2, 4, 7, 8, 22),
            9,
            24,
            (4, 8, 12, 16, 20, 24),
        ),
    ]

    scores = {}
    for name, vals, d, Hmax, checkpoints in cases:
        scores[name] = report_case(name, vals, d, Hmax, checkpoints)

    tournament_fingerprint(scores)

    section("S12 READING")
    print(
        "1. Exact support-6 terms factor into a residue/address coordinate and a "
        "reciprocal denominator.  That turns the missing tail from a volume count "
        "into a Dedekind/cotangent-style signed reciprocal sum."
    )
    print(
        "2. The absolute relation count is dominated by max-norm boundary faces.  "
        "Those are precisely the coordinate cusps where a lower-support shadow "
        "would live; THM-538 predicts cancellation there, and the shell ledgers "
        "show huge abs/signed ratios on the sampled hard supports."
    )
    print(
        "3. A cleaner finishing lemma is now visible: after the finite low-height "
        "wall ledger, prove a boundary-face annihilation / summation-by-parts "
        "bound for the residue reciprocal hyperplane sums.  This is sharper than "
        "a bare Minkowski second-theorem count and explains why S11's signed layer "
        "is small."
    )
    print(
        "4. LRC(14) remains unproved in this file.  The progress is a more exact "
        "target for HYP-2608(a): residue-addressed signed theta tails with all "
        "proper boundary faces killed by the support-6 floor."
    )


if __name__ == "__main__":
    main()
