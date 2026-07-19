#!/usr/bin/env python3
"""Exact positive-control referee for THM-1263's mod-23 rung.

The canonical twelve-set {1,...,12} has exact loneliness margin 1/13,
strictly below 2/25 and 2/23.  At every scale b/23 it supplies a runner
within 2/23 of one witnessed integer, and its antipodal residue counts are
1,...,1,2.  The blocker control {1,...,11,23} shows why the hypothesis
23 ∤ v_i cannot be erased: it is also below 2/23 but has a zero residue and
no doubled nonzero antipodal pair.

Tournament analysis uses antipodal pair obligations as vertices.  Orienting
by (multiplicity, canonical label) gives a transitive loss-audit tournament;
it forgets runner identity and residue sign, so it is not the proof object.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


P = 23
CANONICAL = tuple(range(1, 13))
BLOCKER = tuple(range(1, 12)) + (23,)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1263 referee failed: {message}")


def optimization_safe_require_probe() -> None:
    caught = False
    try:
        require(False, "deliberate optimization-safety probe")
    except RuntimeError as error:
        caught = "deliberate optimization-safety probe" in str(error)
    require(caught, "require probe did not fire")


def candidate_denominators(speeds: tuple[int, ...]) -> list[int]:
    values = {2 * speed for speed in speeds}
    for left, right in combinations(speeds, 2):
        values.add(left + right)
        values.add(abs(left - right))
    values.discard(0)
    return sorted(values)


def exact_margin(speeds: tuple[int, ...]) -> tuple[Fraction, int, int]:
    """Exact piecewise-linear maximum, returned with denominator/numerator."""
    best = Fraction(0)
    witness = (1, 0)
    for denominator in candidate_denominators(speeds):
        for numerator in range(1, denominator // 2 + 1):
            margin = min(
                Fraction(min((speed * numerator) % denominator,
                             denominator - (speed * numerator) % denominator),
                         denominator)
                for speed in speeds
            )
            if margin > best:
                best = margin
                witness = (denominator, numerator)
    return best, *witness


def pair_label(residue: int) -> int | None:
    residue %= P
    if residue == 0:
        return None
    return min(residue, P - residue)


def pair_counts(speeds: tuple[int, ...]) -> tuple[dict[int, int], int]:
    counts = {label: 0 for label in range(1, 12)}
    zeros = 0
    for speed in speeds:
        label = pair_label(speed)
        if label is None:
            zeros += 1
        else:
            counts[label] += 1
    return counts, zeros


def close_witness(speed: int, scale: int) -> tuple[int, Fraction]:
    quotient, residue = divmod(speed * scale, P)
    if residue > P // 2:
        quotient += 1
    error = abs(Fraction(speed * scale, P) - quotient)
    return quotient, error


def all_scale_witnesses(speeds: tuple[int, ...]) -> list[tuple[int, int, int, Fraction]]:
    witnesses = []
    for scale in range(P):
        choices = []
        for speed in speeds:
            integer, error = close_witness(speed, scale)
            if error < Fraction(2, P):
                choices.append((speed, integer, error))
        require(choices, f"missing corrected existential witness at scale {scale}")
        speed, integer, error = min(choices, key=lambda item: (item[2], item[0]))
        witnesses.append((scale, speed, integer, error))
    return witnesses


def positive_compositions(total: int, parts: int):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in positive_compositions(total - first, parts - 1):
            yield (first,) + rest


def tournament_loss_audit(counts: dict[int, int]) -> None:
    order = sorted(counts, key=lambda label: (counts[label], label))
    require(len(order) == 11, "tournament vertex count")
    scores = list(range(11))
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=eleven antipodal pair obligations {+-1},...,{+-11}")
    print("observable=(pair_multiplicity,canonical_label); gauge=canonical_label")
    print("score_histogram=" + ",".join(map(str, scores)))
    print("directed_3_cycles=0 scc_sizes=" + ",".join("1" for _ in order))
    print("hamiltonian_path_count=1")
    print("tie_hamiltonian_path=" + ",".join(map(str, order)))
    print("preserves=missing-pair and multiplicity predicates")
    print("destroys=runner identity and residue sign; tournament is not proof carrier")


def main() -> None:
    optimization_safe_require_probe()
    canonical_margin, canonical_q, canonical_a = exact_margin(CANONICAL)
    blocker_margin, blocker_q, blocker_a = exact_margin(BLOCKER)
    require(canonical_margin == Fraction(1, 13), "canonical exact margin")
    require(canonical_margin < Fraction(2, 25) < Fraction(2, 23), "threshold order")
    canonical_counts, canonical_zeros = pair_counts(CANONICAL)
    require(canonical_zeros == 0, "canonical nonzero residues")
    require(canonical_counts[11] == 2, "canonical double pair")
    require(all(canonical_counts[label] == 1 for label in range(1, 11)),
            "canonical singleton pairs")
    witnesses = all_scale_witnesses(CANONICAL)
    require(all(error < Fraction(2, P) for _, _, _, error in witnesses),
            "strict existential closeness")

    occupancy = list(positive_compositions(12, 11))
    require(len(occupancy) == 11, "positive occupancy census size")
    require(all(sorted(row) == [1] * 10 + [2] for row in occupancy),
            "surjection occupancy law")

    blocker_counts, blocker_zeros = pair_counts(BLOCKER)
    require(blocker_margin < Fraction(2, P), "blocker is genuinely below rung")
    require(blocker_zeros == 1, "blocker zero residue")
    require(all(value == 1 for value in blocker_counts.values()),
            "blocker has no doubled nonzero pair")

    print("THM-1263 MOD-23 NEAR-BIJECTION EXACT REFEREE")
    print("optimization_safe_require_probe=PASS")
    print(
        f"positive_control={CANONICAL} exact_margin={canonical_margin} "
        f"witness_time={canonical_a}/{canonical_q}"
    )
    print("threshold_chain=1/13<2/25<2/23")
    print(
        "pair_counts="
        + ",".join(f"{label}:{canonical_counts[label]}" for label in canonical_counts)
        + f" zero_residues={canonical_zeros}"
    )
    print("corrected_closeness=PASS scales=0..22 quantifiers=forall_b_exists_i_exists_m")
    print(
        "scale_witnesses="
        + ";".join(
            f"b{scale}:v{speed},m{integer},err={error}"
            for scale, speed, integer, error in witnesses
        )
    )
    print("surjective_occupancy_census=11 compositions; all pattern 2,1x10")
    print(
        f"divisible_control={BLOCKER} exact_margin={blocker_margin} "
        f"witness_time={blocker_a}/{blocker_q} zero_residues={blocker_zeros} "
        "nonzero_pair_counts=1x11 conclusion_without_hunit=FALSE"
    )
    tournament_loss_audit(canonical_counts)
    print("RESULT=positive control passes; divisible-by-23 exclusion is necessary")


if __name__ == "__main__":
    main()
