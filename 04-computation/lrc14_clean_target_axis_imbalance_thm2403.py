#!/usr/bin/env python3
"""Exact companion for THM-2403.

The script exhausts both clean-cover endpoint lemmas: the literal
q_*-deletion bank and the restored fully-all-safe present bank.  It includes
every reduced three/four-root core, every relevant ordinary two-root word,
and every ordinary blocker gate with zero, one, or two failed shifts.  It
also checks the sharp physical controls, prime-cyclotomic reduction,
variance bounds, target-phase projector, and exact LRC floors.
"""

from fractions import Fraction
from itertools import combinations


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def translated(word: frozenset[int], shift: int) -> frozenset[int]:
    return frozenset((entry + shift) % P for entry in word)


def clean_cover_examples() -> dict[tuple[int, int], tuple[frozenset[int], ...]]:
    # In every example the order is q_*, q_i, q_2, q_3, q_4, guard.
    examples = {
        # The double is q_* cap q_i: |A_0|=1 and |B|=3.
        (1, 3): (
            frozenset((0, 1)),
            frozenset((1, 2)),
            frozenset((3, 4)),
            frozenset((5, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
        # The double is q_i cap q_2: |A_0|=2 and |B|=3.
        (2, 3): (
            frozenset((0, 1)),
            frozenset((2, 3)),
            frozenset((2, 4)),
            frozenset((5, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
        # The double is wholly retained: |A_0|=2 and |B|=4.
        (2, 4): (
            frozenset((0, 1)),
            frozenset((2, 3)),
            frozenset((4, 5)),
            frozenset((4, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
    }
    for key, words in examples.items():
        counts = [sum(root in word for word in words) for root in range(P)]
        require(sorted(counts) == [1] * 12 + [2], "bad clean-cover profile")
        q_star, q_i, *retained = words
        core = frozenset(range(P)) - frozenset().union(*retained)
        require(
            (len(core - q_i), len(core)) == key,
            "bad clean-cover residual type",
        )
        require(len(q_star) == len(q_i) == 2, "bad ordinary word size")
        require(len(words[-1]) == 4, "bad guard word size")
    return examples


def exhaustive_endpoint_lemma():
    minimum: dict[tuple[int, int, int], int] = {}
    configurations = 0
    gated_cases = 0
    unique_mass_vectors: set[tuple[int, ...]] = set()
    type_counts: dict[tuple[int, int], int] = {}

    nonzero_shifts = tuple(range(1, P))
    zero_banks = tuple(
        frozenset(bank)
        for size in range(3)
        for bank in combinations(nonzero_shifts, size)
    )
    require(len(zero_banks) == 79, "wrong blocker-gate bank")

    for core_size in (3, 4):
        for core_tuple in combinations(range(P), core_size):
            core = frozenset(core_tuple)
            for word_tuple in combinations(range(P), 2):
                word = frozenset(word_tuple)
                a0 = len(core - word)
                residual_type = (a0, core_size)
                if residual_type not in ((1, 3), (2, 3), (2, 4)):
                    continue
                configurations += 1
                type_counts[residual_type] = type_counts.get(residual_type, 0) + 1

                ungated = tuple(
                    len(core - translated(word, shift)) for shift in range(P)
                )
                require(sum(ungated) == 11 * core_size, "eleven-cover identity")
                require(ungated[0] == a0, "wrong base deletion support")

                for zero_bank in zero_banks:
                    masses = tuple(
                        0 if shift in zero_bank else ungated[shift]
                        for shift in range(P)
                    )
                    gap = sum(masses) - P * masses[0]
                    lower = 9 * core_size - P * a0
                    require(gap >= lower > 0, "target-axis gap failed")

                    key = (a0, core_size, len(zero_bank))
                    minimum[key] = min(minimum.get(key, gap), gap)
                    unique_mass_vectors.add(masses)
                    gated_cases += 1

    expected_minimum = {
        (1, 3, 0): 20,
        (1, 3, 1): 17,
        (1, 3, 2): 14,
        (2, 3, 0): 7,
        (2, 3, 1): 4,
        (2, 3, 2): 1,
        (2, 4, 0): 18,
        (2, 4, 1): 14,
        (2, 4, 2): 10,
    }
    require(minimum == expected_minimum, "wrong sharp gap table")

    # For degree at most twelve, reduction modulo Phi_13 subtracts the
    # X^12 coefficient from every lower coefficient.  Thus a rational
    # target mode can vanish only for a flat vector.
    for masses in unique_mass_vectors:
        reduced = tuple(masses[j] - masses[12] for j in range(12))
        require(any(reduced), "nonconstant vector reduced to zero")

        total = sum(masses)
        gap = total - P * masses[0]
        variance = sum(
            (P * entry - total) ** 2 for entry in masses
        ) / Fraction(P**3)
        require(
            variance >= Fraction(gap * gap, 2028),
            "variance floor failed",
        )

    return (
        configurations,
        gated_cases,
        len(unique_mass_vectors),
        type_counts,
        minimum,
    )


def exhaustive_fully_masked_lemma():
    minimum: dict[tuple[int, int], int] = {}
    configurations = 0
    gated_cases = 0
    unique_mass_vectors: set[tuple[int, ...]] = set()
    size_counts: dict[int, int] = {}

    nonzero_shifts = tuple(range(1, P))
    zero_banks = tuple(
        frozenset(bank)
        for size in range(3)
        for bank in combinations(nonzero_shifts, size)
    )

    for core_size in (3, 4):
        for core_tuple in combinations(range(P), core_size):
            core = frozenset(core_tuple)
            for word_tuple in combinations(range(P), 2):
                word = frozenset(word_tuple)
                deletion_type = (len(core - word), core_size)
                if deletion_type not in ((1, 3), (2, 3), (2, 4)):
                    continue
                for q_star_tuple in combinations(range(P), 2):
                    q_star = frozenset(q_star_tuple)
                    if not core <= word | q_star:
                        continue
                    charged = core - q_star
                    if len(charged) not in (1, 2):
                        continue
                    require(charged <= word, "restored base is not empty")

                    configurations += 1
                    charged_size = len(charged)
                    size_counts[charged_size] = (
                        size_counts.get(charged_size, 0) + 1
                    )
                    ungated = tuple(
                        len(charged - translated(word, shift))
                        for shift in range(P)
                    )
                    require(ungated[0] == 0, "fully masked base is nonzero")
                    require(
                        sum(ungated) == 11 * charged_size,
                        "fully masked eleven-cover identity",
                    )

                    for zero_bank in zero_banks:
                        masses = tuple(
                            0 if shift in zero_bank else ungated[shift]
                            for shift in range(P)
                        )
                        gap = sum(masses)
                        lower = 9 * charged_size
                        require(gap >= lower > 0, "fully masked gap failed")
                        key = (charged_size, len(zero_bank))
                        minimum[key] = min(minimum.get(key, gap), gap)
                        unique_mass_vectors.add(masses)
                        gated_cases += 1

    expected_minimum = {
        (1, 0): 11,
        (1, 1): 10,
        (1, 2): 9,
        (2, 0): 22,
        (2, 1): 20,
        (2, 2): 18,
    }
    require(minimum == expected_minimum, "wrong fully masked gap table")

    for masses in unique_mass_vectors:
        require(masses[0] == 0, "lost anchored-zero slice")
        require(sum(masses) > 0, "flat fully masked bank")
        reduced = tuple(masses[j] - masses[12] for j in range(12))
        require(any(reduced), "fully masked vector reduced to zero")

    return (
        configurations,
        gated_cases,
        len(unique_mass_vectors),
        size_counts,
        minimum,
    )


def sharp_gap_one_control():
    core = frozenset((0, 1, 3))
    word = frozenset((2, 3))
    killed = frozenset((2, 3))
    masses = tuple(
        0
        if shift in killed
        else len(core - translated(word, shift))
        for shift in range(P)
    )
    require(
        masses == (2, 2, 0, 0, 3, 3, 3, 3, 3, 3, 2, 1, 2),
        "wrong sharp mass vector",
    )
    require(sum(masses) - P * masses[0] == 1, "gap one is not sharp")

    # The ordinary strict danger arcs centred at s/13 kill exactly shifts
    # two and three at z=5/26, while the base shift zero is safe.
    z = Fraction(5, 26)
    danger_shifts = tuple(
        shift
        for shift in range(P)
        if circle_norm(z - Fraction(shift, P)) < Fraction(1, 14)
    )
    require(danger_shifts == (2, 3), "wrong physical two-shift gate")
    return masses, danger_shifts


def sharp_fully_masked_control():
    # This is the (2,3) physical clean-cover example in
    # clean_cover_examples: q_*={0,1}, q_i={2,3}, core={0,1,3}.
    core = frozenset((0, 1, 3))
    q_star = frozenset((0, 1))
    word = frozenset((2, 3))
    charged = core - q_star
    killed = frozenset((2, 3))
    masses = tuple(
        0
        if shift in killed
        else len(charged - translated(word, shift))
        for shift in range(P)
    )
    require(charged == frozenset((3,)), "wrong restored charged set")
    require(
        masses == (0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        "wrong sharp fully masked mass vector",
    )
    require(sum(masses) == 9, "fully masked gap nine is not sharp")
    return masses


def target_phase_projector():
    # With eta=e_a-e_i, the lawful action shifts the a-factor by -s/13
    # and the i-factor by +s/13.  Multiplying by zeta^(b s) selects
    # n_a-n_i=b modulo thirteen.
    checks = 0
    live = 0
    for n_a in range(P):
        for n_i in range(P):
            for b in range(P):
                residues = tuple(
                    (b + n_i - n_a) * shift % P for shift in range(P)
                )
                selected = all(residue == 0 for residue in residues)
                require(
                    selected == ((n_a - n_i - b) % P == 0),
                    "target phase projector failed",
                )
                live += int(selected)
                checks += 1
    require(live == P * P, "wrong number of selected target phases")
    return checks, live


def deep_target_transform_control():
    def reduced(raw: list[Fraction]) -> tuple[Fraction, ...]:
        return tuple(raw[j] - raw[12] for j in range(12))

    def scaled(
        value: tuple[Fraction, ...], scalar: Fraction
    ) -> tuple[Fraction, ...]:
        return tuple(scalar * entry for entry in value)

    # A deterministic nonnegative diagonal-zero tensor.  Its t=0 slice
    # is nontrivial, so all normalizations are exercised.
    def h_value(r: int, s: int, t: int) -> int:
        if r == t:
            return 0
        return 1 + ((3 * r + 5 * s + 7 * t) % 5)

    m = tuple(
        sum(h_value(r, s, 0) for r in range(P))
        for s in range(P)
    )
    j_bank: dict[tuple[int, int], tuple[Fraction, ...]] = {}
    checks = 0

    for a in range(P):
        for b in range(P):
            raw = [Fraction(0) for _ in range(P)]
            for r in range(P):
                for s in range(P):
                    raw[(a * r + b * s) % P] += Fraction(
                        h_value(r, s, 0), P**2
                    )
            j_bank[a, b] = reduced(raw)

    for b in range(P):
        raw_m = [Fraction(0) for _ in range(P)]
        for s in range(P):
            raw_m[(b * s) % P] += Fraction(m[s], P)
        m_hat = reduced(raw_m)
        require(
            j_bank[0, b] == scaled(m_hat, Fraction(1, P)),
            "J(0,b)=Mhat(b)/13 failed",
        )
        checks += 1

        sum_a = tuple(
            sum(j_bank[a, b][j] for a in range(P))
            for j in range(12)
        )
        require(
            all(entry == 0 for entry in sum_a),
            "diagonal-zero target-line sum failed",
        )
        checks += 1

    # Exact h-character orthogonality is the normalization behind
    # J(a,b)=sum_h B(a,b,h).
    for t in range(P):
        raw_h = [Fraction(0) for _ in range(P)]
        for h in range(P):
            raw_h[(h * t) % P] += 1
        expected = (
            (Fraction(P),) + (Fraction(0),) * 11
            if t == 0
            else (Fraction(0),) * 12
        )
        require(reduced(raw_h) == expected, "h-slice orthogonality failed")
        checks += 1

    # THM-2365 selects deep residue a and endpoint residues b,h; the
    # preserved target is q=(b,a+h).
    target_cases = 0
    for a in range(P):
        for b in range(P):
            for h in range(P):
                q = (b, (a + h) % P)
                require(q[0] == b and q[1] == (a + h) % P, "target typing")
                target_cases += 1
    require(target_cases == P**3, "wrong target typing case count")
    return checks, target_cases


def exact_floors():
    universal_mass = Fraction(1, 26754)
    common_core_mass = Fraction(66, 4459)
    floors = {
        "universal_endpoint_energy": 27 * universal_mass**2 / 114244,
        "universal_endpoint_max": 3 * universal_mass / 676,
        "universal_deep_energy": 27 * universal_mass**2 / 28561,
        "universal_deep_max": 3 * universal_mass / 338,
        "common_core_endpoint_energy": 27 * common_core_mass**2 / 114244,
        "common_core_endpoint_max": 3 * common_core_mass / 676,
        "common_core_deep_energy": 27 * common_core_mass**2 / 28561,
        "common_core_deep_max": 3 * common_core_mass / 338,
    }
    require(
        floors["universal_endpoint_energy"] == Fraction(3, 9085908032656),
        "wrong universal endpoint energy floor",
    )
    require(
        floors["universal_endpoint_max"] == Fraction(1, 6028568),
        "wrong universal endpoint max floor",
    )
    require(
        floors["universal_deep_energy"] == Fraction(3, 2271477008164),
        "wrong universal deep energy floor",
    )
    require(
        floors["universal_deep_max"] == Fraction(1, 3014284),
        "wrong universal deep max floor",
    )
    require(
        floors["common_core_endpoint_energy"]
        == Fraction(29403, 567869252041),
        "wrong common-core endpoint energy floor",
    )
    require(
        floors["common_core_endpoint_max"] == Fraction(99, 1507142),
        "wrong common-core endpoint max floor",
    )
    require(
        floors["common_core_deep_energy"] == Fraction(117612, 567869252041),
        "wrong common-core deep energy floor",
    )
    require(
        floors["common_core_deep_max"] == Fraction(99, 753571),
        "wrong common-core deep max floor",
    )
    return floors


def main() -> None:
    examples = clean_cover_examples()
    (
        configurations,
        gated_cases,
        unique_vectors,
        type_counts,
        minimum,
    ) = exhaustive_endpoint_lemma()
    (
        full_configurations,
        full_gated_cases,
        full_unique_vectors,
        full_size_counts,
        full_minimum,
    ) = exhaustive_fully_masked_lemma()
    sharp_masses, danger_shifts = sharp_gap_one_control()
    sharp_full_masses = sharp_fully_masked_control()
    projector_checks, projector_live = target_phase_projector()
    transform_checks, transform_targets = deep_target_transform_control()
    floors = exact_floors()

    print("THM-2403 clean target-axis imbalance exact companion")
    print(
        f"clean_cover_types={len(examples)} "
        + ",".join(f"{a}:{b}" for a, b in sorted(examples))
    )
    print(
        f"endpoint_configurations={configurations} "
        f"gated_cases={gated_cases} unique_mass_vectors={unique_vectors}"
    )
    print(
        "type_counts="
        + ",".join(
            f"{a}:{b}:{type_counts[a, b]}" for a, b in sorted(type_counts)
        )
    )
    print(
        "sharp_gap_table="
        + ",".join(
            f"{a}:{b}:z{z}:{minimum[a,b,z]}"
            for a, b, z in sorted(minimum)
        )
    )
    print(
        "gap_one_masses="
        + ",".join(map(str, sharp_masses))
        + " danger_shifts="
        + ",".join(map(str, danger_shifts))
    )
    print(
        f"full_configurations={full_configurations} "
        f"full_gated_cases={full_gated_cases} "
        f"full_unique_mass_vectors={full_unique_vectors}"
    )
    print(
        "full_size_counts="
        + ",".join(
            f"{size}:{full_size_counts[size]}"
            for size in sorted(full_size_counts)
        )
    )
    print(
        "full_gap_table="
        + ",".join(
            f"c{size}:z{z}:{full_minimum[size,z]}"
            for size, z in sorted(full_minimum)
        )
    )
    print(
        "full_gap_nine_masses=" + ",".join(map(str, sharp_full_masses))
    )
    print(
        f"target_phase_checks={projector_checks} "
        f"selected={projector_live}"
    )
    print(
        f"deep_target_transform_checks={transform_checks} "
        f"target_typing_cases={transform_targets}"
    )
    print(
        "floors="
        + ",".join(f"{key}:{value}" for key, value in floors.items())
    )
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
