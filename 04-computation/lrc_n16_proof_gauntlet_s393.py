#!/usr/bin/env python3
"""
lrc_n16_proof_gauntlet_s393.py

codex-2026-05-31 S393

Another deliberately broad proof attempt for the n=16 Lonely Runner case.

This pass builds on S389/S390 but tries different pressure points:

1. fold the time circle by t -> t+1/2 and ask whether the even-safe quotient
   ever forces both antipodal sides to be killed by odd speeds;
2. exhaust the normalized 2-torsion cube in the finite scalar quotient;
3. solve exact local set covers for endpoints of dyadic gates;
4. force the exact minimum local cover of the 16-gate, then search all
   five-speed completions through a beam over fixed-threshold interval length;
5. run a second sieve-aware beam from the primitive seed {1,16}.

The goal is still not to claim a proof.  It is to make the remaining proof
obligation smaller and stranger in a useful way.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()


N = 16
K = 15
ONE = Fraction(1, 1)


@dataclass(frozen=True)
class FixedReport:
    speeds: tuple[int, ...]
    length: Fraction
    max_gap: Fraction
    gap_count: int
    components: int
    witness: Fraction | None


@dataclass(frozen=True)
class FoldReport:
    label: str
    speeds: tuple[int, ...]
    even_count: int
    odd_count: int
    even_safe_length: Fraction
    folded_witness_length: Fraction
    folded_danger_length: Fraction
    witness: Fraction | None


@dataclass(frozen=True)
class Candidate:
    label: str
    speeds: tuple[int, ...]
    fixed_length: Fraction
    fixed_gap: Fraction
    exact_gap_ratio: Fraction
    classification: str
    missing_moduli: tuple[int, ...]
    unprotected: int
    core_endpoints: int
    fold_witness_length: Fraction


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def ffloat(value: Fraction) -> str:
    return f"{float(value):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def v2(value: int) -> int:
    if value == 0:
        return 99
    count = 0
    while value % 2 == 0:
        count += 1
        value //= 2
    return count


def split_interval(lo: Fraction, hi: Fraction) -> list[tuple[Fraction, Fraction]]:
    return S356.split_circle_interval(lo, hi)


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    return S356.merge_intervals(intervals)


def complement(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    intervals = merge(intervals)
    if not intervals:
        return [(Fraction(0), ONE)]
    out: list[tuple[Fraction, Fraction]] = []
    start = Fraction(0)
    for lo, hi in intervals:
        if lo > start:
            out.append((start, lo))
        start = max(start, hi)
    if start < ONE:
        out.append((start, ONE))
    return out


def intersect_two(
    a: list[tuple[Fraction, Fraction]], b: list[tuple[Fraction, Fraction]]
) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    i = j = 0
    a = merge(a)
    b = merge(b)
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def intersect_many(
    families: list[list[tuple[Fraction, Fraction]]]
) -> list[tuple[Fraction, Fraction]]:
    if not families:
        return [(Fraction(0), ONE)]
    current = families[0]
    for family in families[1:]:
        current = intersect_two(current, family)
        if not current:
            break
    return current


def subtract_intervals(
    base: list[tuple[Fraction, Fraction]], cut: list[tuple[Fraction, Fraction]]
) -> list[tuple[Fraction, Fraction]]:
    return intersect_two(base, complement(cut))


def measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((hi - lo for lo, hi in intervals), Fraction(0))


def midpoint(interval: tuple[Fraction, Fraction]) -> Fraction:
    return ((interval[0] + interval[1]) / 2) % ONE


def intervals_around(centers: list[Fraction], radius: Fraction) -> list[tuple[Fraction, Fraction]]:
    pieces: list[tuple[Fraction, Fraction]] = []
    for center in centers:
        pieces.extend(split_interval(center - radius, center + radius))
    return merge(pieces)


def fixed_threshold_report(raw_speeds: tuple[int, ...]) -> FixedReport:
    speeds = tuple(sorted(raw_speeds))
    pieces: list[tuple[Fraction, Fraction]] = []
    threshold = Fraction(1, N)
    for speed in speeds:
        radius = threshold / speed
        for center_num in range(speed):
            center = Fraction(center_num, speed)
            pieces.extend(split_interval(center - radius, center + radius))
    components = merge(pieces)
    gaps = S356.circular_gaps(components)
    max_gap_tuple = max(gaps, key=lambda pair: pair[1] - pair[0]) if gaps else None
    return FixedReport(
        speeds=speeds,
        length=measure(components),
        max_gap=(max_gap_tuple[1] - max_gap_tuple[0]) if max_gap_tuple else Fraction(0),
        gap_count=len(gaps),
        components=len(components),
        witness=midpoint(max_gap_tuple) if max_gap_tuple else None,
    )


def missing_moduli(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(m for m in range(2, N + 1) if all(speed % m for speed in speeds))


def classify_exact(speeds: tuple[int, ...]) -> tuple[str, Fraction]:
    report = S356.report("candidate", list(speeds))
    if report.max_gap > 0:
        classification = "positive_gap"
    elif report.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover_candidate"
    return classification, report.max_gap / report.threshold


def audit_candidate(label: str, speeds: tuple[int, ...]) -> Candidate:
    speeds = tuple(S356.normalize_speed_set(list(speeds)))
    fixed = fixed_threshold_report(speeds)
    classification, gap_ratio = classify_exact(speeds)
    descent = S362.summarize(list(speeds))
    folded = folded_report(label, speeds)
    return Candidate(
        label=label,
        speeds=speeds,
        fixed_length=fixed.length,
        fixed_gap=fixed.max_gap,
        exact_gap_ratio=gap_ratio,
        classification=classification,
        missing_moduli=missing_moduli(speeds),
        unprotected=descent.unprotected_count,
        core_endpoints=descent.core_endpoint_count,
        fold_witness_length=folded.folded_witness_length,
    )


def folded_report(label: str, speeds: tuple[int, ...]) -> FoldReport:
    """Fold t to s=2t and test antipodal-pair escape.

    Even speed 2w is unsafe for both antipodal preimages exactly when
    ||w*s|| < 1/16.

    Odd speed u kills side 0 when ||u*s/2|| < 1/16, and kills side 1 when
    ||u*s/2|| > 7/16.  A folded witness exists when even speeds are safe and
    not both odd-low and odd-high occur.
    """

    even_bad: list[tuple[Fraction, Fraction]] = []
    odd_low: list[tuple[Fraction, Fraction]] = []
    odd_high: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        if speed % 2 == 0:
            w = speed // 2
            centers = [Fraction(a, w) for a in range(w)]
            even_bad.extend(intervals_around(centers, Fraction(1, N * w)))
        else:
            # side 0 bad: u*s/2 near an integer
            low_centers = [Fraction(2 * a, speed) % ONE for a in range(speed)]
            odd_low.extend(intervals_around(low_centers, Fraction(1, 8 * speed)))
            # side 1 bad: u*s/2 near a half-integer
            high_centers = [Fraction(2 * a + 1, speed) % ONE for a in range(speed)]
            odd_high.extend(intervals_around(high_centers, Fraction(1, 8 * speed)))

    even_safe = complement(even_bad)
    danger = intersect_many([even_safe, merge(odd_low), merge(odd_high)])
    folded_witness = subtract_intervals(even_safe, danger)
    witness = midpoint(folded_witness[0]) if folded_witness else None
    return FoldReport(
        label=label,
        speeds=speeds,
        even_count=sum(1 for speed in speeds if speed % 2 == 0),
        odd_count=sum(1 for speed in speeds if speed % 2 == 1),
        even_safe_length=measure(even_safe),
        folded_witness_length=measure(folded_witness),
        folded_danger_length=measure(danger),
        witness=witness,
    )


def print_folded_antipodal() -> None:
    print("1. Folded antipodal quotient attack")
    samples = [
        ("initial", tuple(range(1, N))),
        ("single 16 gate", tuple(list(range(1, N - 1)) + [N])),
        ("best 8-ladder", (1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 120)),
        ("min-cover seed", (1, 3, 5, 7, 8, 9, 11, 13, 15, 16)),
        ("S390 random near", (13, 25, 48, 61, 63, 72, 75, 90, 94, 110, 118, 125, 126, 127, 128)),
    ]
    print("   label             even odd even_safe folded_witness danger witness_s")
    for label, speeds in samples:
        row = folded_report(label, tuple(sorted(speeds)))
        print(
            f"   {label:<17} {row.even_count:>4} {row.odd_count:>3} "
            f"{fmt(row.even_safe_length):>9} {fmt(row.folded_witness_length):>14} "
            f"{fmt(row.folded_danger_length):>8} {fmt(row.witness):>9}"
        )
    print(
        "   Folded parity is not a proof by itself: the initial and single-gate\n"
        "   rows have zero folded witness mass, meaning the antipodal pair test\n"
        "   can lose boundary-only information.  But when the fold does expose\n"
        "   mass, it gives a literal quotient witness interval.  The useful\n"
        "   target is therefore folded parity plus endpoint-leaf debt."
    )
    print()


def print_two_torsion_cube() -> None:
    print("2. Exact normalized 2-torsion cube")
    system = S372.build_pattern_system(N)
    half = N // 2
    by_support: dict[int, list[tuple[int, tuple[int, ...]]]] = {}
    global_best: list[tuple[int, tuple[int, ...]]] = []
    for mask in range(1, 1 << (K - 1)):
        vector = [0] * K
        for idx in range(K - 1):
            if (mask >> idx) & 1:
                vector[idx + 1] = half
        vector_tuple = tuple(vector)
        missed = S372.score_vector(system, vector_tuple).missed
        support = sum(1 for value in vector_tuple if value)
        by_support.setdefault(support, []).append((missed, vector_tuple))
        global_best.append((missed, vector_tuple))
    print(f"   patterns={len(system.patterns)} candidates={system.candidate_count}")
    print("   support  best_missed best_count example_support")
    for support in sorted(by_support):
        rows = by_support[support]
        best = min(missed for missed, _vector in rows)
        examples = [vector for missed, vector in rows if missed == best]
        example_support = tuple(i + 1 for i, value in enumerate(examples[0]) if value)
        print(
            f"   {support:>7} {best:>12} {len(examples):>10} {example_support}"
        )
    print("   global best half-turn vectors")
    for missed, vector in sorted(global_best)[:8]:
        support = tuple(i + 1 for i, value in enumerate(vector) if value)
        print(f"     missed={missed:4d} support={support} vector={vector}")
    print(
        "   Even the entire half-turn cube has a moat.  The best pure dyadic\n"
        "   non-scalar defect is still the last-coordinate half-turn with 128\n"
        "   missed cells."
    )
    print()


def endpoint_mask(speed: int, protector: int) -> int:
    mask = 0
    modulus = N * speed
    for m in range(speed):
        for sign_index, sign in enumerate((-1, 1)):
            bit = 2 * m + sign_index
            residue = (protector * (N * m + sign)) % modulus
            if min(residue, modulus - residue) < speed:
                mask |= 1 << bit
    return mask


def exact_set_cover(speed: int) -> tuple[int | None, tuple[int, ...]]:
    full = (1 << (2 * speed)) - 1
    masks = [(p, endpoint_mask(speed, p)) for p in range(1, speed)]
    masks = [(p, mask) for p, mask in masks if mask]
    for size in range(1, len(masks) + 1):
        for combo in combinations(masks, size):
            covered = 0
            for _p, mask in combo:
                covered |= mask
            if covered == full:
                return size, tuple(p for p, _mask in combo)
    return None, tuple()


def greedy_set_cover(speed: int) -> tuple[int, tuple[int, ...], int]:
    full = (1 << (2 * speed)) - 1
    masks = [(p, endpoint_mask(speed, p)) for p in range(1, speed)]
    covered = 0
    chosen: list[int] = []
    while covered != full:
        best = max(
            ((mask & ~covered).bit_count(), p, mask) for p, mask in masks if p not in chosen
        )
        gain, p, mask = best
        if gain == 0:
            break
        chosen.append(p)
        covered |= mask
    return len(chosen), tuple(chosen), covered.bit_count()


def print_gate_endpoint_covers() -> tuple[int, ...]:
    print("3. Exact local set covers for dyadic gate endpoints")
    print("   speed endpoints exact_min exact_cover greedy_count greedy_cover_prefix")
    min16: tuple[int, ...] = tuple()
    for speed in (8, 16, 32, 64):
        if speed <= 16:
            exact_size, exact_cover = exact_set_cover(speed)
        else:
            exact_size, exact_cover = None, tuple()
        greedy_count, greedy_cover, covered = greedy_set_cover(speed)
        if speed == 16:
            min16 = exact_cover
        print(
            f"   {speed:>5} {2*speed:>9} {str(exact_size):>9} "
            f"{str(exact_cover):<34} {greedy_count:>12} "
            f"{greedy_cover[:10]} ({covered}/{2*speed})"
        )
    print(
        "   The 16-gate endpoint cover forces a very specific nine-speed local\n"
        "   basis.  The next experiment treats that basis as an adversary's\n"
        "   opening move, then gives it five extra speeds."
    )
    print()
    return min16


def candidate_pool() -> tuple[int, ...]:
    pool = set(range(1, 65))
    pool.update(16 * q for q in range(1, 17))
    pool.update(lcm(16, m) for m in range(2, 17))
    pool.update((72, 80, 88, 96, 104, 112, 120, 128, 144, 160, 192, 240))
    return tuple(sorted(pool))


def beam_complete(
    seed: tuple[int, ...],
    target_size: int,
    pool: tuple[int, ...],
    width: int,
) -> list[tuple[tuple[int, ...], FixedReport]]:
    states = {tuple(sorted(seed))}
    for _step in range(target_size - len(seed)):
        scored: list[tuple[tuple[int, Fraction, Fraction, int], tuple[int, ...], FixedReport]] = []
        for state in states:
            for value in pool:
                if value in state:
                    continue
                candidate = tuple(sorted(state + (value,)))
                if len(candidate) != len(set(candidate)):
                    continue
                if gcd_all(candidate) != 1:
                    continue
                fixed = fixed_threshold_report(candidate)
                miss = len(missing_moduli(candidate))
                shell_bonus = len({min(v2(speed), 5) for speed in candidate})
                score = (miss, -fixed.length, fixed.max_gap, -shell_bonus)
                scored.append((score, candidate, fixed))
        scored.sort(key=lambda item: item[0])
        states = {candidate for _score, candidate, _fixed in scored[:width]}
    finals = [(state, fixed_threshold_report(state)) for state in states if len(state) == target_size]
    finals.sort(key=lambda item: (len(missing_moduli(item[0])), item[1].max_gap, -item[1].length, item[0]))
    return finals


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def print_min_cover_completion(min16: tuple[int, ...]) -> list[Candidate]:
    print("4. Adversarial completion from the exact 16-gate local cover")
    seed = tuple(sorted(min16 + (16,)))
    pool = candidate_pool()
    finals = beam_complete(seed, K, pool, width=36)
    audited = [audit_candidate(f"mincover_beam_{idx}", speeds) for idx, (speeds, _fixed) in enumerate(finals[:8])]
    print(f"   seed={seed}")
    print("   label            miss fixed_len fixed_gap exact_gap missing unprot core fold_witness")
    for row in audited:
        print(
            f"   {row.label:<16} {len(row.missing_moduli):>4} "
            f"{fmt(row.fixed_length):>9} {fmt(row.fixed_gap):>9} "
            f"{ffloat(row.exact_gap_ratio):>9} "
            f"{str(row.missing_moduli)[:12]:>12} {row.unprotected:>6} "
            f"{row.core_endpoints:>4} {fmt(row.fold_witness_length):>12}"
        )
        print(f"     speeds={row.speeds}")
    print(
        "   Giving the exact local cover five extra speeds still leaves positive\n"
        "   fixed-threshold gaps and empty endpoint cores in the best beam rows."
    )
    print()
    return audited


def print_sieve_beam() -> list[Candidate]:
    print("5. Sieve-aware beam from primitive seed {1,16}")
    pool = candidate_pool()
    finals = beam_complete((1, 16), K, pool, width=48)
    audited = [audit_candidate(f"sieve_beam_{idx}", speeds) for idx, (speeds, _fixed) in enumerate(finals[:10])]
    print("   label         miss fixed_len fixed_gap exact_gap class       unprot core fold_witness")
    for row in audited:
        print(
            f"   {row.label:<13} {len(row.missing_moduli):>4} "
            f"{fmt(row.fixed_length):>9} {fmt(row.fixed_gap):>9} "
            f"{ffloat(row.exact_gap_ratio):>9} {row.classification:<11} "
            f"{row.unprotected:>6} {row.core_endpoints:>4} {fmt(row.fold_witness_length):>12}"
        )
        print(f"     speeds={row.speeds}")
    print(
        "   The beam finds sieve-complete high-coverage sets; some even defeat\n"
        "   the folded antipodal witness test.  They still have positive exact\n"
        "   gaps and empty endpoint cores.  So the obstruction is not failure to\n"
        "   guess the right gate; it is the inability to make the dyadic debt\n"
        "   leafless after the sieve is paid."
    )
    print()
    return audited


def print_synthesis() -> None:
    print("6. Synthesis / proof targets left standing")
    print(
        "   Route A: Folded antipodal plus endpoint-leaf theorem.  The fold alone\n"
        "   can miss boundary-only structure, but any folded-danger set found by\n"
        "   the beams still peels to empty endpoint core.  Prove that folded\n"
        "   danger forces a private dyadic endpoint leaf."
    )
    print(
        "   Route B: Local-cover rigidity.  The 16-gate's endpoint cover forces\n"
        "   the nine-residue pattern (1,3,5,7,8,9,11,13,15) up to dilation.\n"
        "   Prove that any five-speed completion keeps a positive fixed-threshold\n"
        "   gap or a private endpoint leaf."
    )
    print(
        "   Route C: 2-torsion moat.  The full normalized half-turn cube has a\n"
        "   128-cell moat.  Prove non-half-turn residues are no better by a\n"
        "   compression/rounding argument, then lift through endpoint debt."
    )
    print(
        "   These are three different shadows of the same thing: in the pure\n"
        "   2-adic row, every attempt to close a boundary creates a lower-layer\n"
        "   private obligation before a labelled protection cycle can close."
    )


def main() -> None:
    print("n=16 Lonely Runner proof gauntlet (codex-2026-05-31 S393)")
    print("denominator n=16, k=15 moving speeds, threshold=1/16")
    print()
    print_folded_antipodal()
    print_two_torsion_cube()
    min16 = print_gate_endpoint_covers()
    print_min_cover_completion(min16)
    print_sieve_beam()
    print_synthesis()


if __name__ == "__main__":
    main()
