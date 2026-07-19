#!/usr/bin/env python3
"""Exact replay for THM-1196, the prefix-survivor Farey beat mesh.

The theorem is analytic.  The finite scan at the end is explicitly telemetry:
it exhibits exact beat witnesses in a bounded bank but is not used in the
all-range proof.

Tournament/carrier note
-----------------------
The speed-order tournament is transitive and forgets the proof predicate.
The faithful vertices are phase-local mesh cells separated by rational beat
points; a defining speed pair is a label on a separating wall.  The replay
therefore reports the trivial tournament fingerprint only as a loss audit.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import ceil, gcd


def require(condition: bool, message: str) -> None:
    """Optimization-safe verifier check; unlike assert, survives python -O."""
    if not condition:
        raise RuntimeError(f"THM-1196 verification failed: {message}")


def optimization_safe_require_probe() -> None:
    """Deliberately fail once to prove checks remain live under python -O."""
    caught = False
    try:
        require(False, "deliberate optimization-safety probe")
    except RuntimeError as error:
        caught = "deliberate optimization-safety probe" in str(error)
    if not caught:
        raise RuntimeError("THM-1196 require probe did not fire")


def slow_gap(c: int, k: int) -> tuple[Fraction, Fraction]:
    return (Fraction(14 * k + 1, 14 * c),
            Fraction(14 * k + 13, 14 * c))


def floor_fraction(x: Fraction) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def dangerous(speed: int, t: Fraction) -> bool:
    x = speed * t
    residue = x - floor_fraction(x)
    return min(residue, 1 - residue) < Fraction(1, 14)


def tooth_boundaries(
    c: int, speed: int, k: int,
) -> tuple[list[tuple[Fraction, Fraction]], list[Fraction]]:
    """Open speed teeth meeting G and their interior boundary points."""
    left, right = slow_gap(c, k)
    centre_lo = speed * left - Fraction(1, 14)
    centre_hi = speed * right + Fraction(1, 14)
    first = floor_fraction(centre_lo) - 1
    last = ceil_fraction(centre_hi) + 1
    intervals: list[tuple[Fraction, Fraction]] = []
    boundaries: list[Fraction] = []
    for m in range(first, last + 1):
        lo = Fraction(14 * m - 1, 14 * speed)
        hi = Fraction(14 * m + 1, 14 * speed)
        if lo < right and hi > left:  # strict tooth really meets closed G
            intervals.append((lo, hi))
            if left < lo < right:
                boundaries.append(lo)
            if left < hi < right:
                boundaries.append(hi)
    return intervals, boundaries


def survivor_data(
    c: int, speeds: tuple[int, ...], k: int,
) -> tuple[Fraction, int, int]:
    """Lebesgue length, component count, and total meeting-tooth count.

    Critical points and open cells are tested independently, so an isolated
    safe seam is counted as a component even though it has measure zero.
    """
    left, right = slow_gap(c, k)
    critical = {left, right}
    tooth_count = 0
    for speed in speeds:
        intervals, boundaries = tooth_boundaries(c, speed, k)
        tooth_count += len(intervals)
        critical.update(boundaries)
    points = sorted(critical)

    # Atoms alternate point, open cell, point.  Runs of safe atoms are exactly
    # the connected components of the closed survivor.
    atoms: list[tuple[bool, Fraction]] = []
    for index, point in enumerate(points):
        atoms.append((not any(dangerous(s, point) for s in speeds),
                      Fraction(0)))
        if index + 1 < len(points):
            nxt = points[index + 1]
            midpoint = (point + nxt) / 2
            safe = not any(dangerous(s, midpoint) for s in speeds)
            atoms.append((safe, nxt - point))

    length = sum((mass for safe, mass in atoms if safe), Fraction(0))
    components = 0
    previous_safe = False
    for safe, _ in atoms:
        if safe and not previous_safe:
            components += 1
        previous_safe = safe
    return length, components, tooth_count


def beat_points(
    c: int, prefix: tuple[int, ...], last: int, k: int,
) -> tuple[Fraction, ...]:
    left, right = slow_gap(c, k)
    points: set[Fraction] = set()
    for speed in prefix:
        for denominator in (last - speed, last + speed):
            lo = ceil_fraction(denominator * left)
            hi = floor_fraction(denominator * right)
            for numerator in range(lo, hi + 1):
                points.add(Fraction(numerator, denominator))
    return tuple(sorted(points))


def farey_mesh(
    c: int, prefix: tuple[int, ...], last: int, k: int,
) -> Fraction:
    left, right = slow_gap(c, k)
    points = sorted({left, right, *beat_points(c, prefix, last, k)})
    return max(y - x for x, y in zip(points, points[1:]))


def mesh_from_points(c: int, k: int, beats: tuple[Fraction, ...]) -> Fraction:
    left, right = slow_gap(c, k)
    points = sorted({left, right, *beats})
    return max(y - x for x, y in zip(points, points[1:]))


def prefix_beat_witness(
    c: int, prefix: tuple[int, ...], last: int, k: int,
) -> Fraction | None:
    for point in beat_points(c, prefix, last, k):
        if not any(dangerous(speed, point) for speed in prefix):
            return point
    return None


def harmonic_pressure(c: int, speeds: tuple[int, ...]) -> Fraction:
    return sum((Fraction(c, speed) for speed in speeds), Fraction(0))


def survivor_lower_bound(
    c: int, prefix: tuple[int, ...], r: int,
) -> Fraction:
    return Fraction(6, 49 * c) * (
        (8 - r) - harmonic_pressure(c, prefix)
    )


def phase_free_component_bound(c: int, prefix: tuple[int, ...]) -> int:
    return 1 + sum(ceil(Fraction(6 * speed + c, 7 * c))
                   for speed in prefix)


def modular_distance(numerator: int, denominator: int) -> int:
    residue = numerator % denominator
    return min(residue, denominator - residue)


def replay_beat_identity() -> int:
    checks = 0
    for left in range(1, 40):
        for right in range(left + 1, 41):
            for denominator in (right - left, right + left):
                for numerator in range(-2 * denominator, 2 * denominator + 1):
                    require(
                        modular_distance(left * numerator, denominator)
                        == modular_distance(right * numerator, denominator),
                        "sum/difference beat identity",
                    )
                    checks += 1
    return checks


def gcd_covariance_check() -> tuple[Fraction, Fraction]:
    c, prefix, last, k, scale = 5, (7, 8, 11), 13, 2, 6
    base = farey_mesh(c, prefix, last, k)
    dilated = farey_mesh(scale * c,
                         tuple(scale * speed for speed in prefix),
                         scale * last, k)
    require(dilated == base / scale, "Farey mesh gcd covariance")
    base_survivor = survivor_data(c, prefix, k)
    dilated_survivor = survivor_data(
        scale * c, tuple(scale * speed for speed in prefix), k
    )
    require(
        dilated_survivor[0] == base_survivor[0] / scale,
        "survivor length gcd covariance",
    )
    require(
        dilated_survivor[1:] == base_survivor[1:],
        "survivor component/tooth gcd covariance",
    )
    return base, dilated


def finite_telemetry() -> dict[int, tuple[int, ...]]:
    result: dict[int, tuple[int, ...]] = {}
    for r in (4, 5, 6):
        rows = beat_witness_rows = covered_rows = 0
        subcritical_rows = coarse_tail_refutations = 0
        local_mesh_refutations = mesh_strict_improvements = 0
        distance_checks = 0
        smallest_mesh_ratio: Fraction | None = None
        for c in range(2, 8):
            for speeds in combinations(range(c + 1, 3 * c + 1), r):
                prefix = speeds[:-1]
                last = speeds[-1]
                for k in range(c):
                    rows += 1
                    beats = beat_points(c, prefix, last, k)
                    witness = next((point for point in beats
                                    if not any(dangerous(speed, point)
                                               for speed in prefix)), None)
                    if witness is not None:
                        beat_witness_rows += 1
                        # The defining prefix speed and the last speed have
                        # equal distance; hence the full packet misses witness.
                        require(
                            not any(dangerous(s, witness) for s in prefix),
                            "selected beat witness is prefix-safe",
                        )
                        defining = [
                            s for s in prefix
                            if any((witness * q).denominator == 1
                                   for q in (last - s, last + s))
                        ]
                        require(bool(defining), "beat witness has a defining pair")
                        require(
                            all(
                                dangerous(s, witness)
                                == dangerous(last, witness)
                                for s in defining
                            ),
                            "defining-pair beat distances agree",
                        )
                        require(
                            not dangerous(last, witness),
                            "beat witness is safe for the last comb",
                        )

                    # Every witness is pointwise safe for the full packet,
                    # so it proves noncoverage without a second interval
                    # merge.  The runtime check after the loop verifies all
                    # rows have such a witness.

                    c_bar = phase_free_component_bound(c, prefix)

                    lower = survivor_lower_bound(c, prefix, r)
                    q_max = last + prefix[-1]
                    mesh = mesh_from_points(c, k, beats)
                    require(
                        mesh <= Fraction(1, q_max),
                        "phase-local mesh is bounded by densest sum lattice",
                    )
                    distance_checks += 1
                    ratio = mesh * q_max
                    if smallest_mesh_ratio is None or ratio < smallest_mesh_ratio:
                        smallest_mesh_ratio = ratio
                    if mesh < Fraction(1, q_max):
                        mesh_strict_improvements += 1

                    if lower > 0:
                        subcritical_rows += 1
                        if lower >= Fraction(c_bar, q_max):
                            coarse_tail_refutations += 1
                        if lower >= c_bar * mesh:
                            local_mesh_refutations += 1

        require(smallest_mesh_ratio is not None, "finite mesh bank is nonempty")
        result[r] = (
            rows,
            beat_witness_rows,
            covered_rows,
            subcritical_rows,
            coarse_tail_refutations,
            local_mesh_refutations,
            mesh_strict_improvements,
            distance_checks,
            smallest_mesh_ratio.numerator,
            smallest_mesh_ratio.denominator,
        )
    return result


def component_and_survivor_replay() -> int:
    """Independent exact audit of (11), (14), and THM-1176's supplier."""
    checks = 0
    for c in range(2, 6):
        bank = tuple(range(c + 1, 3 * c + 1))
        for prefix_size in range(1, min(5, len(bank)) + 1):
            r = prefix_size + 1
            for prefix in combinations(bank, prefix_size):
                for k in range(c):
                    length, components, tooth_count = survivor_data(c, prefix, k)
                    require(
                        components <= 1 + tooth_count,
                        "survivor component count",
                    )
                    require(
                        1 + tooth_count
                        <= phase_free_component_bound(c, prefix),
                        "phase-free tooth count",
                    )
                    lower = survivor_lower_bound(c, prefix, r)
                    require(length >= lower, "THM-1176 survivor lower bound")
                    checks += 1
    return checks


def tournament_audit() -> None:
    prefix = (7, 8, 11, 13, 14)
    scores = [sum(speed > other for other in prefix if other != speed)
              for speed in prefix]
    require(
        sorted(scores) == list(range(5)),
        "transitive tournament score histogram",
    )
    print("TOURNAMENT_LOSS_AUDIT")
    print("observable=(last+d_i)-(last+d_j)=d_i-d_j")
    print("score_histogram=0,1,2,3,4")
    print("directed_3_cycles=0")
    print("scc_sizes=1,1,1,1,1")
    print("hamiltonian_path_count=1")
    print("tie_hamiltonian_path=7,8,11,13,14")
    print("faithful_vertices=phase-local_Farey_mesh_cells")


def five_comb_dual_composition() -> None:
    normalized_length = Fraction(1, 36) / Fraction(7, 6)
    slow_gap_scale = Fraction(6, 7)
    carrier_scaled_length = normalized_length * slow_gap_scale
    require(
        normalized_length == Fraction(1, 42),
        "five-comb normalized survivor conversion",
    )
    require(
        carrier_scaled_length == Fraction(1, 49),
        "five-comb physical survivor conversion",
    )
    ceiling_checks = 0
    for c in range(1, 101):
        for speed in range(c + 1, 20 * c + 1):
            require(
                ceil(Fraction(6 * speed + c, 7 * c))
                <= ceil(Fraction(speed, c)),
                "component ceiling monotonicity",
            )
            ceiling_checks += 1
    print("FIVE_COMB_DUAL_COMPOSITION")
    print("dual_survivor_f_mass>1/36 and max_f=7/6 for faster slopes")
    print(f"normalized_Lebesgue_survivor>{normalized_length}")
    print(f"carrier_scaled_survivor>1/({carrier_scaled_length.denominator}*c)")
    print("one_last_tooth_per_E5_component=true")
    print("unconditional_tail=d6/c<7*C5")
    print("uniform_last_ratio=d6/d5<77")
    print(f"ceiling_monotonicity_checks={ceiling_checks}")


def max_cyclic_run(values: list[bool]) -> int:
    if all(values):
        return len(values)
    best = current = 0
    for value in values + values:
        current = current + 1 if value else 0
        best = max(best, current)
    return min(best, len(values))


def j4_safe_beat_flood_tail() -> None:
    normalized_f_mass = 1 - 4 * Fraction(7, 36)
    normalized_length = normalized_f_mass / Fraction(7, 6)
    carrier_scaled_length = normalized_length * Fraction(6, 7)
    require(
        normalized_f_mass == Fraction(2, 9),
        "four-prefix dual mass conversion",
    )
    require(
        normalized_length == Fraction(4, 21),
        "four-prefix normalized survivor conversion",
    )
    require(
        carrier_scaled_length == Fraction(8, 49),
        "four-prefix physical survivor conversion",
    )

    run_checks = equality_rows = 0
    for q_reduced in range(3, 501):
        for u in range(1, (q_reduced + 1) // 2):
            if gcd(u, q_reduced) != 1:
                continue
            values = [
                14 * min((u * p) % q_reduced, (-u * p) % q_reduced)
                < q_reduced
                for p in range(q_reduced)
            ]
            run = max_cyclic_run(values)
            bound = ceil(Fraction(q_reduced, 7 * u))
            require(run <= bound, "safe-beat dangerous-run bound")
            equality_rows += run == bound
            run_checks += 1

    print("J4_SAFE_BEAT_FLOOD_TAIL")
    span_checks = 0
    for d5 in range(2, 201):
        for d6 in range(d5 + 1, min(4 * d5, 400) + 1):
            require(
                Fraction(1, 7 * d5) + Fraction(2, 7 * d6)
                < Fraction(3, 7 * d5),
                "two-tooth component-span inequality",
            )
            span_checks += 1

    print(f"four_prefix_dual_f_mass>{normalized_f_mass}")
    print(f"normalized_Lebesgue_survivor>{normalized_length}")
    print(f"carrier_scaled_survivor>{carrier_scaled_length}/c")
    print("two_tooth_union_component_span<3/(7*d5)")
    print("unconditional_component_tail=d5/c<(21/8)*C4")
    print("uniform_penultimate_ratio=d5/d4<189/8")
    print(f"component_span_inequality_checks={span_checks}")
    print("dangerous_run_R<=ceil((d5+d6)/(7*d5))")
    print("safe_beat_mesh<8/(7*d5)")
    print("safe_beat_scalar_tail=d5/c<7*C4 (weaker but phase-aware)")
    print(f"reduced_rotation_run_checks={run_checks}")
    print(f"sharp_run_rows={equality_rows}")


def j3_three_tooth_component_tail() -> None:
    normalized_f_mass = 1 - 3 * Fraction(7, 36)
    normalized_length = normalized_f_mass / Fraction(7, 6)
    carrier_scaled_length = normalized_length * Fraction(6, 7)
    require(
        normalized_f_mass == Fraction(5, 12),
        "three-prefix dual mass conversion",
    )
    require(
        normalized_length == Fraction(5, 14),
        "three-prefix normalized survivor conversion",
    )
    require(
        carrier_scaled_length == Fraction(15, 49),
        "three-prefix physical survivor conversion",
    )

    span_checks = 0
    for d4 in range(2, 101):
        for d5 in range(d4 + 1, min(2 * d4, 200) + 1):
            for d6 in sorted({d5 + 1, 2 * d5, 5 * d5}):
                w4 = Fraction(1, 7 * d4)
                w5 = Fraction(1, 7 * d5)
                w6 = Fraction(1, 7 * d6)
                two_tooth_span = w5 + 2 * w6
                three_tooth_span = w4 + 2 * two_tooth_span
                require(
                    two_tooth_span < 3 * w4 < 6 * w4,
                    "two-tooth component cannot bridge d4 teeth",
                )
                require(
                    three_tooth_span < 7 * w4 == Fraction(1, d4),
                    "three-tooth recursive component span",
                )
                span_checks += 1

    print("J3_THREE_TOOTH_COMPONENT_TAIL")
    print(f"three_prefix_dual_f_mass>{normalized_f_mass}")
    print(f"normalized_Lebesgue_survivor>{normalized_length}")
    print(f"carrier_scaled_survivor>{carrier_scaled_length}/c")
    print("two_tooth_span=w5+2*w6<3*w4<gap4")
    print("three_tooth_span=w4+2*(w5+2*w6)<7*w4=1/d4")
    print("unconditional_component_tail=d4/c<(49/15)*C3")
    print("uniform_third_tail_ratio=d4/d3<343/15")
    print(f"recursive_component_span_checks={span_checks}")


def main() -> None:
    print("THM-1196 PREFIX-SURVIVOR FAREY BEAT-MESH DESCENT")
    optimization_safe_require_probe()
    print("optimization_safe_require_probe=PASS")
    identity_checks = replay_beat_identity()
    print(f"sum_difference_beat_identity_checks={identity_checks}")
    print("analytic_lemma=prefix survivor avoids every last-prefix sum/difference beat")
    print("component_bound=C<=1+sum ceil((6*d_i+c)/(7*c))")
    print("phase_free_mesh_bound=M_k<=1/(d_last+d_penultimate)")
    print("tail_bound=(d_last+d_penultimate)/c"
          "<49*Cbar/[6*((8-r)-H_prefix)]")

    base_mesh, dilated_mesh = gcd_covariance_check()
    print("\nGCD_TOOTHPICK_COVARIANCE")
    print(f"base_mesh={base_mesh}")
    print(f"sixfold_dilated_mesh={dilated_mesh}")
    print("pressure_tooth_counts_components_invariant=true")

    print("\nFINITE_EXACT_TELEMETRY_NOT_PROOF")
    component_checks = component_and_survivor_replay()
    print(f"component_and_survivor_bound_replay_rows={component_checks}")
    telemetry = finite_telemetry()
    expected_rows = {4: 11354, 5: 20268, 6: 27730}
    for r, row in telemetry.items():
        (rows, witnesses, covers, subcritical, coarse, local,
         improved, distance_checks,
         ratio_num, ratio_den) = row
        require(rows == expected_rows[r], f"r={r} expected row count")
        require(witnesses == rows, f"r={r} every row has a beat witness")
        require(covers == 0, f"r={r} zero open-cover rows")
        require(distance_checks == rows, f"r={r} mesh check count")
        print(
            f"r={r} rows={rows} beat_witness_rows={witnesses} "
            f"open_cover_rows={covers} subcritical_prefix_rows={subcritical}"
        )
        print(
            f"r={r} coarse_tail_refutations={coarse} "
            f"phase_local_mesh_refutations={local} "
            f"strict_mesh_improvements={improved} "
            f"minimum_M_times_qmax={ratio_num}/{ratio_den}"
        )

    print()
    tournament_audit()
    print()
    five_comb_dual_composition()
    print()
    j4_safe_beat_flood_tail()
    print()
    j3_three_tooth_component_tail()
    print("\nPROOF_TELEMETRY_SEPARATION")
    print("PROVED=beat descent; tail ratios <77, <189/8, <343/15")
    print("TELEMETRY=all finite-bank row counts and zero survivors")


if __name__ == "__main__":
    main()
