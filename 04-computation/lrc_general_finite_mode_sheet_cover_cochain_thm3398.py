#!/usr/bin/env python3
"""Exact companion for THM-3398: general finite-mode sheet covers.

The analytic theorem is uniform in the fibre degree.  This standard-library
companion attacks its two delicate points independently:

* every actual phase-grid hit pattern through q=29, including boundaries and
  wraparound, is the unique maximal consecutive block predicted by the mode
  formula; and
* mode-lattice pair overlap and bounded full-cover existence agree with
  independent rational-lattice and event-cell calculations.

All arithmetic is integral or Fraction-exact.  Runtime gates use RuntimeError,
so ``python -O`` retains every decision.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY_PINS = {
    "01-canon/theorems/THM-3395-small-sheet-typed-cover-star-cochain.md":
        "a64f472f167cf739d727da2386f38aa3f1c4175dd97cec6e87bdbed6ab256300",
    "04-computation/lrc14_small_sheet_typed_cover_star_cochain_thm3395.py":
        "0953401d98d62fd3118bd4a7bbeb50bd43a459a04dc120578ca6af355cada067",
    "05-knowledge/results/lrc14_small_sheet_typed_cover_star_cochain_thm3395.out":
        "f7ed05e16fdd3660741aa8a79600cf9920bbebd8087c8d25a252ecca0dbc1ce5",
    "04-computation/lrc14_q8_domino_mode_clutter_probe_20260814.py":
        "3e523a2ff8cbd6329782347c56fae2d8519a161c3d127697ca452f3891890b9c",
    "05-knowledge/results/lrc14_q8_domino_mode_clutter_probe_20260814.out":
        "0f5a421205bc559c8f12dce8462b4d570fcffba0e602740d1ea66c52cd84d045",
}

EXPECTED_LOCAL_SUMMARY = (
    868,
    42_328,
    2_444,
    47_200,
    "646f9f66df6447a55739e27cf57d4a19b56397786882b5d1517a8c1be676b833",
)
EXPECTED_PAIR_COMPARISONS = 40_187
EXPECTED_GLOBAL_PROFILE_SHA256 = "c4985ead4c085fc4100143dbf26988c8e83f9dc5598915094614453d8d3c04e4"
EXPECTED_SEMANTIC_SHA256 = "b6d3663a25b45b732d4acc597dc2772a5adfd9860934dab7a26c4174f268f6ff"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_repr(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def circle_distance(left: Fraction, right: Fraction) -> Fraction:
    forward = (left - right) % 1
    return min(forward, 1 - forward)


def phase_danger(y: Fraction, phase: int, phase_count: int) -> bool:
    value = (y + Fraction(phase, phase_count)) % 1
    return 14 * min(value, 1 - value) < 1


def source_danger(q: int, speed: int, t: Fraction, sheet: int) -> bool:
    value = (speed * (t + Fraction(sheet, q))) % 1
    return 14 * min(value, 1 - value) < 1


def owner_modes(q: int, speed: int) -> tuple[tuple[object, ...], ...]:
    """Return (block,h,w,r,s,g,m) selected-block modes.

    The semantic is monotone: every sheet in ``block`` is required to fire;
    sheets outside it may also fire.  This is deliberately not an exact-state
    partition.
    """
    common = math.gcd(q, speed)
    phase_count = q // common
    phase_step = (speed // common) % phase_count
    modes = []
    for size in range(1, (phase_count + 6) // 7 + 1):
        width = q - 7 * common * (size - 1)
        require(width > 0, ("positive mode width", q, speed, size, width))
        for start in range(phase_count):
            phase_block = frozenset((start + offset) % phase_count for offset in range(size))
            sheet_block = frozenset(
                sheet
                for sheet in range(q)
                if (phase_step * sheet) % phase_count in phase_block
            )
            centre = (-common * (2 * start + size - 1)) % (2 * q)
            modes.append((sheet_block, centre, width, start, size, common, phase_count))
    return tuple(modes)


def gap_values(
    q: int,
    left_speed: int,
    right_speed: int,
    left_mode: tuple[object, ...],
    right_mode: tuple[object, ...],
) -> tuple[int, ...]:
    """Exact p=2q*u*v*(x_u-x_v) values for strict mode overlap."""
    left_centre = int(left_mode[1])
    right_centre = int(right_mode[1])
    left_width = int(left_mode[2])
    right_width = int(right_mode[2])
    modulus = 2 * q * math.gcd(left_speed, right_speed)
    residue = (left_centre * right_speed - right_centre * left_speed) % modulus
    strict_numerator = left_width * right_speed + right_width * left_speed
    bound = (strict_numerator - 1) // 7
    first = residue % modulus
    lower = (-bound - first + modulus - 1) // modulus
    upper = (bound - first) // modulus
    return tuple(first + modulus * index for index in range(lower, upper + 1))


def independent_pair_overlap(
    q: int,
    left_speed: int,
    right_speed: int,
    left_mode: tuple[object, ...],
    right_mode: tuple[object, ...],
) -> bool:
    """Overlap from rational lattice distance, without using the p fibre."""
    left_offset = Fraction(int(left_mode[1]), 2 * q * left_speed)
    right_offset = Fraction(int(right_mode[1]), 2 * q * right_speed)
    step = Fraction(math.gcd(left_speed, right_speed), left_speed * right_speed)
    ratio = (left_offset - right_offset) / step
    floor = ratio.numerator // ratio.denominator
    fractional = ratio - floor
    distance = min(fractional, 1 - fractional) * step
    radius_sum = Fraction(int(left_mode[2]), 14 * q * left_speed) + Fraction(
        int(right_mode[2]), 14 * q * right_speed
    )
    return distance < radius_sum


def complete_cochain(
    q: int,
    speeds: Sequence[int],
    modes: Sequence[tuple[object, ...]],
) -> Optional[tuple[tuple[int, int, int], ...]]:
    if len(speeds) <= 1:
        return ()
    stars = tuple(
        gap_values(q, speeds[0], speeds[index], modes[0], modes[index])
        for index in range(1, len(speeds))
    )
    if any(not fibre for fibre in stars):
        return None
    for star in itertools.product(*stars):
        edges = {(0, index): star[index - 1] for index in range(1, len(speeds))}
        compatible = True
        for left in range(1, len(speeds)):
            for right in range(left + 1, len(speeds)):
                numerator = (
                    speeds[left] * star[right - 1]
                    - speeds[right] * star[left - 1]
                )
                if numerator % speeds[0] != 0:
                    compatible = False
                    break
                value = numerator // speeds[0]
                if value not in gap_values(
                    q, speeds[left], speeds[right], modes[left], modes[right]
                ):
                    compatible = False
                    break
                edges[(left, right)] = value
            if not compatible:
                break
        if compatible:
            return tuple(
                (left, right, edges[(left, right)])
                for left in range(len(speeds))
                for right in range(left + 1, len(speeds))
            )
    return None


def combine_congruences(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    a, modulus_a = left
    b, modulus_b = right
    common = math.gcd(modulus_a, modulus_b)
    require((b - a) % common == 0, ("incompatible CRT", left, right))
    reduced_b = modulus_b // common
    if reduced_b == 1:
        multiplier = 0
    else:
        multiplier = (
            ((b - a) // common)
            * pow(modulus_a // common, -1, reduced_b)
        ) % reduced_b
    modulus = math.lcm(modulus_a, modulus_b)
    return ((a + modulus_a * multiplier) % modulus, modulus)


def realize_cochain(
    q: int,
    speeds: Sequence[int],
    modes: Sequence[tuple[object, ...]],
    cochain: Sequence[tuple[int, int, int]],
) -> Fraction:
    """Effectivize the proof: generalized CRT, then interval Helly."""
    edge_map = {(left, right): value for left, right, value in cochain}
    potentials = [Fraction(0)]
    for index in range(1, len(speeds)):
        value = edge_map[(0, index)]
        potentials.append(-Fraction(value, 2 * q * speeds[0] * speeds[index]))

    offsets = tuple(
        Fraction(int(mode[1]), 2 * q * speed) - potential
        for speed, mode, potential in zip(speeds, modes, potentials)
    )
    scale = math.lcm(
        *(speed for speed in speeds),
        *(offset.denominator for offset in offsets),
    )
    congruences = tuple(
        (int(scale * offset), scale // speed)
        for speed, offset in zip(speeds, offsets)
    )
    combined = congruences[0]
    for congruence in congruences[1:]:
        combined = combine_congruences(combined, congruence)
    shift = Fraction(combined[0], scale)
    centres = tuple(potential + shift for potential in potentials)

    for speed, mode, centre in zip(speeds, modes, centres):
        lattice_coordinate = speed * centre - Fraction(int(mode[1]), 2 * q)
        require(lattice_coordinate.denominator == 1, ("centre lattice", lattice_coordinate))

    lowers = tuple(
        centre - Fraction(int(mode[2]), 14 * q * speed)
        for speed, mode, centre in zip(speeds, modes, centres)
    )
    uppers = tuple(
        centre + Fraction(int(mode[2]), 14 * q * speed)
        for speed, mode, centre in zip(speeds, modes, centres)
    )
    lower = max(lowers)
    upper = min(uppers)
    require(lower < upper, ("Helly interval", q, speeds, modes, cochain, lower, upper))
    witness = (lower + upper) / 2
    for speed, mode in zip(speeds, modes):
        require(
            all(source_danger(q, speed, witness, int(sheet)) for sheet in mode[0]),
            ("realized mode", q, speed, mode, witness),
        )
    return witness


def mode_cover_all_firing(
    q: int, speeds: Sequence[int]
) -> Optional[tuple[tuple[int, ...], tuple[tuple[object, ...], ...], tuple[tuple[int, int, int], ...], Fraction]]:
    ordered = tuple(sorted(speeds, key=lambda speed: (len(owner_modes(q, speed)), speed)))
    banks = tuple(owner_modes(q, speed) for speed in ordered)
    maximum_sizes = tuple(max(len(mode[0]) for mode in bank) for bank in banks)
    universe = frozenset(range(q))

    def visit(index: int, assigned: list[tuple[object, ...]], covered: frozenset[int]):
        if len(covered) + sum(maximum_sizes[index:]) < q:
            return None
        if index == len(ordered):
            if covered != universe:
                return None
            cochain = complete_cochain(q, ordered, tuple(assigned))
            if cochain is None:
                return None
            witness = realize_cochain(q, ordered, tuple(assigned), cochain)
            return ordered, tuple(assigned), cochain, witness
        speed = ordered[index]
        for mode in banks[index]:
            if any(
                not gap_values(q, ordered[prior], speed, assigned[prior], mode)
                for prior in range(index)
            ):
                continue
            result = visit(index + 1, assigned + [mode], covered | mode[0])
            if result is not None:
                return result
        return None

    return visit(0, [], frozenset())


def event_samples(q: int, speeds: Sequence[int]) -> tuple[int, tuple[int, ...]]:
    scale = 14 * q * math.lcm(*speeds)
    events = set()
    for speed in speeds:
        for sheet in range(q):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // q
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    ordered = tuple(sorted(events))
    samples = [2 * event for event in ordered]
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += scale
        samples.append((left + right) % (2 * scale))
    return scale, tuple(samples)


def integer_sample_danger(
    q: int, speed: int, sample: int, scale: int, sheet: int
) -> bool:
    denominator = 2 * scale * q
    numerator = speed * (q * sample + 2 * scale * sheet)
    residue = numerator % denominator
    return 14 * min(residue, denominator - residue) < denominator


def event_cover(q: int, speeds: Sequence[int]) -> Optional[Fraction]:
    scale, samples = event_samples(q, speeds)
    for sample in samples:
        if all(
            any(integer_sample_danger(q, speed, sample, scale, sheet) for speed in speeds)
            for sheet in range(q)
        ):
            return Fraction(sample, 2 * scale)
    return None


def local_mode_audit() -> tuple[object, ...]:
    owner_count = 0
    mode_count = 0
    wrap_count = 0
    pattern_samples = 0
    pattern_digest = hashlib.sha256()
    for q in range(2, 30):
        for speed in range(1, 2 * q + 1):
            owner_count += 1
            common = math.gcd(q, speed)
            phase_count = q // common
            modes = owner_modes(q, speed)
            expected_count = phase_count * ((phase_count + 6) // 7)
            require(len(modes) == expected_count, ("mode count", q, speed, len(modes)))
            require(len({mode[0] for mode in modes}) == len(modes), ("duplicate blocks", q, speed))
            mode_count += len(modes)

            for mode in modes:
                block, centre, width, start, size, _, _ = mode
                require(len(block) == common * int(size), ("block size", q, speed, mode))
                require(7 * (int(size) - 1) < phase_count, ("span", q, speed, mode))
                require(
                    (int(centre) + common * (2 * int(start) + int(size) - 1)) % (2 * q) == 0,
                    ("centre residue", q, speed, mode),
                )
                if int(start) + int(size) > phase_count:
                    wrap_count += 1
                centre_phase = Fraction(int(centre), 2 * q)
                radius_phase = Fraction(int(width), 14 * q)
                selected_phases = tuple((int(start) + offset) % phase_count for offset in range(int(size)))
                require(
                    all(phase_danger(centre_phase, phase, phase_count) for phase in selected_phases),
                    ("mode centre", q, speed, mode),
                )
                for sign in (-1, 1):
                    endpoint = centre_phase + sign * radius_phase
                    require(
                        not all(phase_danger(endpoint, phase, phase_count) for phase in selected_phases),
                        ("strict endpoint", q, speed, mode, sign),
                    )

            events = {
                (-Fraction(phase, phase_count) + Fraction(sign, 14)) % 1
                for phase in range(phase_count)
                for sign in (-1, 1)
            }
            ordered = tuple(sorted(events))
            samples = set(events)
            for index, left in enumerate(ordered):
                right = ordered[(index + 1) % len(ordered)]
                if index + 1 == len(ordered):
                    right += 1
                samples.add(((left + right) / 2) % 1)
            for y in sorted(samples):
                pattern_samples += 1
                actual = frozenset(
                    phase for phase in range(phase_count) if phase_danger(y, phase, phase_count)
                )
                if not actual:
                    pattern_digest.update(repr((q, speed, y, ())).encode("ascii") + b"\n")
                    continue
                candidates = []
                for mode in modes:
                    start = int(mode[3])
                    size = int(mode[4])
                    phase_block = frozenset((start + offset) % phase_count for offset in range(size))
                    if phase_block == actual:
                        candidates.append(mode)
                require(len(candidates) == 1, ("unique maximal block", q, speed, y, actual, candidates))
                mode = candidates[0]
                require(
                    circle_distance(y, Fraction(int(mode[1]), 2 * q))
                    < Fraction(int(mode[2]), 14 * q),
                    ("mode interval identity", q, speed, y, actual, mode),
                )
                pattern_digest.update(
                    repr((q, speed, y, tuple(sorted(actual)), int(mode[3]), int(mode[4]))).encode("ascii")
                    + b"\n"
                )
    return owner_count, mode_count, wrap_count, pattern_samples, pattern_digest.hexdigest()


def pair_lattice_audit() -> int:
    comparisons = 0
    for q in range(2, 12):
        speeds = tuple(range(1, min(q + 2, 10) + 1))
        for left_speed, right_speed in itertools.combinations(speeds, 2):
            for left_mode in owner_modes(q, left_speed):
                for right_mode in owner_modes(q, right_speed):
                    values = gap_values(q, left_speed, right_speed, left_mode, right_mode)
                    independent = independent_pair_overlap(
                        q, left_speed, right_speed, left_mode, right_mode
                    )
                    require(bool(values) == independent, (
                        "pair lattice mismatch", q, left_speed, right_speed,
                        left_mode, right_mode, values, independent,
                    ))
                    modulus = 2 * q * math.gcd(left_speed, right_speed)
                    residue = (
                        int(left_mode[1]) * right_speed
                        - int(right_mode[1]) * left_speed
                    ) % modulus
                    for value in values:
                        require(value % modulus == residue, ("gap congruence", value, modulus, residue))
                        require(
                            7 * abs(value)
                            < int(left_mode[2]) * right_speed + int(right_mode[2]) * left_speed,
                            ("gap strictness", q, left_speed, right_speed, value),
                        )
                    comparisons += 1
    return comparisons


def bounded_global_audit() -> tuple[tuple[object, ...], ...]:
    profile = []
    for q in range(2, 12):
        pool_limit = min(q + 3, 10) if q <= 8 else 8
        pool = tuple(speed for speed in range(1, pool_limit + 1) if speed % q != 0)
        max_rank = min(5 if q <= 8 else 4, len(pool))
        event_cache: dict[frozenset[int], bool] = {}
        mode_cache: dict[frozenset[int], bool] = {}
        compared = 0
        event_cells = 0
        positive = 0
        minimal = []
        for rank in range(1, max_rank + 1):
            for subset in itertools.combinations(pool, rank):
                key = frozenset(subset)
                inherited = any(
                    event_cache.get(frozenset(proper), False)
                    for size in range(1, rank)
                    for proper in itertools.combinations(subset, size)
                )
                if inherited:
                    event_result = True
                    mode_result = True
                else:
                    maximum_capacity = sum(
                        max(len(mode[0]) for mode in owner_modes(q, speed))
                        for speed in subset
                    )
                    if maximum_capacity < q:
                        event_result = False
                        mode_result = False
                    else:
                        event_cells += 1
                        event_witness = event_cover(q, subset)
                        mode_witness = mode_cover_all_firing(q, subset)
                        event_result = event_witness is not None
                        mode_result = mode_witness is not None
                        if mode_witness is not None:
                            _, modes, _, witness = mode_witness
                            require(
                                frozenset().union(*(mode[0] for mode in modes)) == frozenset(range(q)),
                                ("mode block cover", q, subset, mode_witness),
                            )
                            require(
                                all(
                                    any(source_danger(q, speed, witness, sheet) for speed in mode_witness[0])
                                    for sheet in range(q)
                                ),
                                ("effectivized full cover", q, subset, witness),
                            )
                        require(event_result == mode_result, (
                            "global event/mode mismatch", q, subset,
                            event_witness, mode_witness,
                        ))
                        if event_result:
                            minimal.append(subset)
                event_cache[key] = event_result
                mode_cache[key] = mode_result
                require(event_result == mode_result, ("cache mismatch", q, subset))
                compared += 1
                positive += int(event_result)
        profile.append((q, pool, max_rank, compared, event_cells, positive, tuple(minimal)))
    return tuple(profile)


def control_audit() -> tuple[object, ...]:
    # q=7 equality gives zero width and is excluded; q=8 has the first positive domino.
    require((7 + 6) // 7 == 1, "q7 maximum mode")
    require((8 + 6) // 7 == 2, "q8 maximum mode")
    require(7 - 7 * (2 - 1) == 0, "q7 zero boundary")
    require(8 - 7 * (2 - 1) == 1, "q8 positive boundary")

    tie_speeds = (1, 3, 5, 7)
    tie_blocks = (
        frozenset((0, 1)),
        frozenset((3, 6)),
        frozenset((2, 7)),
        frozenset((4, 5)),
    )
    tie_modes = tuple(
        next(mode for mode in owner_modes(8, speed) if mode[0] == block)
        for speed, block in zip(tie_speeds, tie_blocks)
    )
    tie_cochain = complete_cochain(8, tie_speeds, tie_modes)
    require(tie_cochain == tuple((i, j, 0) for i in range(4) for j in range(i + 1, 4)), tie_cochain)
    tie_time = Fraction(-1, 16)
    require(
        all(
            (speed * tie_time - Fraction(int(mode[1]), 16)).denominator == 1
            for speed, mode in zip(tie_speeds, tie_modes)
        ),
        ("tie centre lattices", tie_modes),
    )
    require(frozenset().union(*tie_blocks) == frozenset(range(8)), "tie partition")

    hostile_speeds = (12, 10, 14)
    hostile_blocks = (
        frozenset((0, 2, 4, 6)),
        frozenset((1, 5)),
        frozenset((3, 7)),
    )
    hostile_modes = tuple(
        next(mode for mode in owner_modes(8, speed) if mode[0] == block)
        for speed, block in zip(hostile_speeds, hostile_blocks)
    )
    hostile_gaps = tuple(
        gap_values(8, hostile_speeds[left], hostile_speeds[right], hostile_modes[left], hostile_modes[right])
        for left, right in itertools.combinations(range(3), 2)
    )
    require(hostile_gaps == ((-16, 16), (-16, 16), (-16, 16)), hostile_gaps)
    require(complete_cochain(8, hostile_speeds, hostile_modes) is None, "hostile closure")
    require(frozenset().union(*hostile_blocks) == frozenset(range(8)), "hostile partition")

    # The first mode beyond the q=8 domino: one q=9 triple plus three dominoes
    # form a zero-cochain partition at t=5/6.
    q9_speeds = (6, 1, 5, 7)
    q9_blocks = (
        frozenset((0, 3, 6)),
        frozenset((1, 2)),
        frozenset((5, 7)),
        frozenset((4, 8)),
    )
    q9_modes = tuple(
        next(mode for mode in owner_modes(9, speed) if mode[0] == block)
        for speed, block in zip(q9_speeds, q9_blocks)
    )
    q9_cochain = complete_cochain(9, q9_speeds, q9_modes)
    require(q9_cochain == tuple((i, j, 0) for i in range(4) for j in range(i + 1, 4)), q9_cochain)
    q9_time = Fraction(5, 6)
    require(
        all(
            (speed * q9_time - Fraction(int(mode[1]), 18)).denominator == 1
            for speed, mode in zip(q9_speeds, q9_modes)
        ),
        ("q9 centre lattices", q9_modes),
    )
    require(frozenset().union(*q9_blocks) == frozenset(range(9)), "q9 partition")

    core_mode = owner_modes(11, 11)
    require(
        len(core_mode) == 1
        and core_mode[0][0] == frozenset(range(11))
        and core_mode[0][1:3] == (0, 11),
        ("m=1 mode", core_mode),
    )
    return (
        ("q7_q8_widths", 0, 1),
        ("q8_zero_tie", tie_speeds, tie_blocks, tie_cochain, tie_time),
        ("q8_pairwise_nonclosure", hostile_speeds, hostile_blocks, hostile_gaps),
        ("q9_triple_domino_zero_tie", q9_speeds, q9_blocks, q9_cochain, q9_time),
        ("m1_core", core_mode),
    )


def main() -> None:
    for relative, expected in DEPENDENCY_PINS.items():
        actual = sha256_lf(ROOT / relative)
        require(actual == expected, ("dependency changed", relative, actual, expected))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    local = local_mode_audit()
    pairs = pair_lattice_audit()
    global_profile = bounded_global_audit()
    controls = control_audit()
    if EXPECTED_LOCAL_SUMMARY is not None:
        require(local == EXPECTED_LOCAL_SUMMARY, ("local summary", local))
    if EXPECTED_PAIR_COMPARISONS is not None:
        require(pairs == EXPECTED_PAIR_COMPARISONS, ("pair comparisons", pairs))
    global_profile_sha256 = digest_repr(global_profile)
    require(
        global_profile_sha256 == EXPECTED_GLOBAL_PROFILE_SHA256,
        ("global profile digest", global_profile_sha256),
    )

    semantic = digest_repr((local, pairs, global_profile, controls))
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    print("THM-3398 GENERAL FINITE-MODE SHEET-COVER COCHAIN")
    print(f"source_sha256_lf={sha256_lf(source)}")
    print(f"dependency_sha256_lf={tuple(DEPENDENCY_PINS.items())}")
    print("status=PROVED analytic all-q selected-block mode iff;FINITE-EXACT local q<=29 and pair/global q<=11")
    print("mode=at_least_block_not_exact_state;(g,m,a)=((u,q),q/g,u/g mod m);1<=s<=ceil(m/7)")
    print("block={ell:a*ell mod m in cyclic[r,r+s)};h=-g(2r+s-1) mod 2q;w=q-7g(s-1)>0")
    print("source_centres=(Z+h/(2q))/u;radius=w/(14qu)")
    print("cochain=p_ij=2q*u_i*u_j*(x_i-x_j);congruence_mod_2qgcd;7|p|<w_i*u_j+w_j*u_i;triangle_closure")
    print(f"local_mode_wrap_q2_q29={local}")
    print(f"pair_lattice_comparisons={pairs}")
    print(f"bounded_global_profile={global_profile}")
    print(f"controls={controls}")
    print("wrap_resolution=unique_short_unwrapping_up_to_integer;h_mod_2q_invariant;extra_owner_hits_allowed")
    print("effectivity=pairwise_generalized_CRT_plus_complete_cochain_gives_real_lifts;open_interval_Helly_gives_source_time")
    print("scope=sheet_cover_carrier_only;no_refined_ledger_decrement;no_LRC14")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
