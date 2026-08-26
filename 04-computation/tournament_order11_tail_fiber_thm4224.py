#!/usr/bin/env python3
"""Exact audit for the order-eleven THM-4224 tail-fiber extension.

This packet deliberately avoids the full order-eleven tournament bank.  It
checks the algebraic consequences of the externally proved forced-tail
bijection, checks the family by an independent subset DP, and exhausts the
small, presentation-dependent sphere obtained by reversing exactly two arcs
of one fixed transitive order.

Only the Python standard library is used.  All comparisons are exact integer
cross multiplications, and all checks remain active under ``python -O``.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
from typing import Iterable


def check(condition: bool, message: str) -> None:
    """Optimization-stable replacement for assert."""

    if not condition:
        raise RuntimeError(message)


def bits(mask: int) -> Iterable[int]:
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


@dataclass(frozen=True)
class Tournament:
    n: int
    out: tuple[int, ...]


@dataclass(frozen=True)
class PathData:
    finish: tuple[int, ...]
    paths: tuple[int, ...]
    full: int


@dataclass(frozen=True)
class Profile:
    h: int
    c: int
    end: tuple[int, ...]
    defect: tuple[int, ...]
    w: int
    energy: int
    target: int
    gap: int


@dataclass(frozen=True)
class PairData:
    i: int
    j: int
    n_cover: int
    b: int
    product: int
    minimum: int
    rho: int


def all_arcs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def from_reversals(n: int, reversals: Iterable[tuple[int, int]]) -> Tournament:
    reverse = frozenset(reversals)
    arcs = all_arcs(n)
    check(reverse <= frozenset(arcs), "invalid reversal")
    out = [0] * n
    for i, j in arcs:
        if (i, j) in reverse:
            out[j] |= 1 << i
        else:
            out[i] |= 1 << j
    return Tournament(n, tuple(out))


def encode(tournament: Tournament) -> str:
    return "".join(
        "1" if tournament.out[i] & (1 << j) else "0"
        for i, j in all_arcs(tournament.n)
    )


def reversal_set(tournament: Tournament) -> tuple[tuple[int, int], ...]:
    return tuple(
        (i, j)
        for i, j in all_arcs(tournament.n)
        if tournament.out[j] & (1 << i)
    )


def incoming_masks(tournament: Tournament) -> tuple[int, ...]:
    full = (1 << tournament.n) - 1
    return tuple(full ^ (1 << i) ^ tournament.out[i] for i in range(tournament.n))


def is_strong(tournament: Tournament) -> bool:
    full = (1 << tournament.n) - 1
    incoming = incoming_masks(tournament)
    for neighborhoods in (tournament.out, incoming):
        seen = 1
        while True:
            expanded = seen
            for vertex in bits(seen):
                expanded |= neighborhoods[vertex]
            if expanded == seen:
                break
            seen = expanded
        if seen != full:
            return False
    return True


def relabel(tournament: Tournament, image: tuple[int, ...]) -> Tournament:
    n = tournament.n
    check(tuple(sorted(image)) == tuple(range(n)), "relabeling is not a permutation")
    out = [0] * n
    for old_source in range(n):
        for old_target in bits(tournament.out[old_source]):
            out[image[old_source]] |= 1 << image[old_target]
    return Tournament(n, tuple(out))


def induced_tournament(tournament: Tournament, vertices: tuple[int, ...]) -> Tournament:
    """Induce and relabel in the supplied vertex order."""

    check(len(set(vertices)) == len(vertices), "repeated induced vertex")
    check(all(0 <= vertex < tournament.n for vertex in vertices), "bad induced vertex")
    new_index = {old: new for new, old in enumerate(vertices)}
    out = [0] * len(vertices)
    for old_source in vertices:
        for old_target in bits(tournament.out[old_source]):
            if old_target in new_index:
                out[new_index[old_source]] |= 1 << new_index[old_target]
    return Tournament(len(vertices), tuple(out))


def path_data(tournament: Tournament) -> PathData:
    n = tournament.n
    states = 1 << n
    full = states - 1
    incoming = incoming_masks(tournament)
    finish = [0] * (states * n)
    paths = [0] * states
    paths[0] = 1

    for mask in range(1, states):
        base = mask * n
        for last in bits(mask):
            rest = mask ^ (1 << last)
            if rest == 0:
                ways = 1
            else:
                previous_mask = rest & incoming[last]
                previous_base = rest * n
                ways = sum(finish[previous_base + previous] for previous in bits(previous_mask))
            finish[base + last] = ways
            paths[mask] += ways
    return PathData(tuple(finish), tuple(paths), full)


def cycle_count(tournament: Tournament) -> int:
    n = tournament.n
    states = 1 << n
    full = states - 1
    start = [0] * (states * n)
    start[n] = 1  # mask={0}, last=0
    for mask in range(1, states):
        if not (mask & 1):
            continue
        base = mask * n
        for last in bits(mask):
            ways = start[base + last]
            if ways == 0:
                continue
            available = tournament.out[last] & (full ^ mask)
            for nxt in bits(available):
                start[(mask | (1 << nxt)) * n + nxt] += ways
    return sum(
        start[full * n + last]
        for last in range(1, n)
        if tournament.out[last] & 1
    )


def tournament_profile(tournament: Tournament, data: PathData | None = None) -> Profile:
    if data is None:
        data = path_data(tournament)
    n = tournament.n
    full = data.full
    h = data.paths[full]
    end = tuple(data.finish[full * n + i] for i in range(n))
    defect = []
    for i in range(n):
        value = 0
        bit = 1 << i
        for mask in range(1, full + 1):
            if mask & bit:
                value += data.finish[mask * n + i] * data.paths[full ^ mask]
        check(value >= h, "negative one-defect coordinate")
        defect.append(value - h)
    defect_tuple = tuple(defect)
    defect_sum = sum(defect_tuple)
    w = (n - 1) * h + defect_sum
    energy = sum((h + value) ** 2 for value in defect_tuple)
    target = (w + h) * (w + 3 * h)
    gap = (
        defect_sum**2
        + 2 * (n - 4) * h * defect_sum
        + n * (n - 3) * h * h
        - 5 * sum(value * value for value in defect_tuple)
    )
    check(gap == target - 5 * energy, "five-copy identities disagree")
    return Profile(h, cycle_count(tournament), end, defect_tuple, w, energy, target, gap)


def two_path_cover(tournament: Tournament, data: PathData, i: int, j: int) -> int:
    n = tournament.n
    full = data.full
    i_bit = 1 << i
    j_bit = 1 << j
    total = 0
    for mask in range(1, full):
        if mask & i_bit and not mask & j_bit:
            total += data.finish[mask * n + i] * data.finish[(full ^ mask) * n + j]
    return total


def pair_data(
    tournament: Tournament,
    data: PathData,
    profile: Profile,
    i: int,
    j: int,
) -> PairData:
    n_cover = two_path_cover(tournament, data, i, j)
    b_value = n_cover + profile.c - profile.end[i] - profile.end[j]
    check(b_value >= 0, "negative exact two-owner count")
    product = profile.end[i] * profile.end[j]
    minimum = min(profile.end[i], profile.end[j])
    incoming = incoming_masks(tournament)
    indegree_i = incoming[i].bit_count()
    indegree_j = incoming[j].bit_count()
    common = (incoming[i] & incoming[j]).bit_count()
    rho = indegree_i * indegree_j - common
    return PairData(i, j, n_cover, b_value, product, minimum, rho)


def all_pair_data(
    tournament: Tournament, data: PathData, profile: Profile
) -> tuple[PairData, ...]:
    return tuple(
        pair_data(tournament, data, profile, i, j)
        for i in range(tournament.n)
        for j in range(i + 1, tournament.n)
    )


def reduced(numerator: int, denominator: int) -> tuple[int, int]:
    divisor = gcd(numerator, denominator)
    return numerator // divisor, denominator // divisor


def fraction_text(numerator: int, denominator: int) -> str:
    top, bottom = reduced(numerator, denominator)
    return f"{top}/{bottom}"


def x_family(m: int) -> Tournament:
    n = m + 5
    return from_reversals(n, ((0, 3), (3, n - 1)))


def core_constants() -> tuple[int, int, int]:
    """Derive THM-4224's constants 5, 32, and 7 by literal subset DP."""

    core = from_reversals(4, ((0, 3),))
    hamilton = [0] * 16
    end_three = [0] * 16
    for mask in range(16):
        vertices = tuple(bits(mask))
        induced = induced_tournament(core, vertices)
        data = path_data(induced)
        hamilton[mask] = data.paths[data.full]
        if mask & (1 << 3):
            local_three = vertices.index(3)
            end_three[mask] = data.finish[data.full * induced.n + local_three]
    full = 15
    partition_sum = sum(hamilton[mask] * hamilton[full ^ mask] for mask in range(16))
    terminal_exclusions = sum(
        hamilton[mask] * end_three[full ^ mask] for mask in range(16)
    )
    check((hamilton[full], partition_sum, terminal_exclusions) == (5, 32, 7), "core constants")
    return hamilton[full], partition_sum, terminal_exclusions


def formal_identity_checks() -> tuple[int, int, int]:
    # Affine pairs represent a*u+b with formal u=2^(t-1).
    def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return left[0] + right[0], left[1] + right[1]

    def scale(factor: int, value: tuple[int, int]) -> tuple[int, int]:
        return factor * value[0], factor * value[1]

    core_h, partition_sum, terminal_exclusions = core_constants()
    end_x = (core_h, 0)
    end_z = (0, core_h)
    n_cover = (partition_sum, 0)
    cycle = (0, 1)
    b_value = add(add(n_cover, scale(-1, end_x)), add(scale(-1, end_z), cycle))
    product = (core_h * core_h, 0)
    check(b_value == (27, -4), "formal B identity")
    check(add(b_value, scale(-1, product)) == (2, -4), "formal product excess")
    check(add(add(b_value, scale(-1, product)), (0, -5)) == (2, -9), "formal +min")
    check(add(scale(27, product), scale(-25, b_value)) == (0, 100), "formal 27/25")

    # Terminal affine pairs use v=2^(m-1).
    terminal_n = (partition_sum, -terminal_exclusions)
    terminal_end = (core_h, 0)
    terminal_b = add(add(terminal_n, scale(-1, terminal_end)), ((0, -5 + 1)))
    check(terminal_b == (27, -11), "formal terminal B")
    check(add(add(terminal_b, scale(-1, product)), (0, -5)) == (2, -16), "terminal +min")
    check(add(scale(27, product), scale(-25, terminal_b)) == (0, 275), "terminal 27/25")

    # Exact integer sign boundaries, plus an induction certificate beyond m=7.
    check(all((2**t - 4 > 0) == (t >= 3) for t in range(1, 65)), "product sign")
    check(all((2**t - 9 > 0) == (t >= 4) for t in range(1, 65)), "+min sign")
    check(2**7 > 2 * 7 + 6, "rook monotonicity base")
    # If 2^m>2m+6, doubling gives 2^(m+1)>4m+12; the next
    # threshold is 2m+8, and their formal affine difference is 2m+4>0.
    check((4 - 2, 12 - 8) == (2, 4), "rook induction affine step")
    check(2 * 7 + 4 > 0, "rook induction sign")
    return core_h, partition_sum, terminal_exclusions


def family_checks(digest: "sha256") -> tuple[Profile, tuple[PairData, ...], Profile, PairData]:
    checked = 0
    x6_profile: Profile | None = None
    x6_pairs: tuple[PairData, ...] | None = None
    x7_profile: Profile | None = None
    x7_terminal: PairData | None = None

    for m in range(2, 8):
        tournament = x_family(m)
        check(is_strong(tournament), "X_m is not strong")
        data = path_data(tournament)
        profile = tournament_profile(tournament, data)
        z = tournament.n - 1
        r = z - 1
        terminal = pair_data(tournament, data, profile, r, z)
        terminal_scale = 1 << (m - 1)
        check(profile.c == 1 and profile.end[z] == 5, "X_m base counts")
        check(profile.end[r] == 5 * terminal_scale, "terminal End formula")
        check(terminal.n_cover == 32 * terminal_scale - 7, "terminal N formula")
        check(terminal.b == 27 * terminal_scale - 11, "terminal B formula")
        check(terminal.rho == m * m + 5 * m + 7, "terminal rho formula")

        for t in range(1, m):
            x_t = 3 + t
            pair = pair_data(tournament, data, profile, x_t, z)
            scale = 1 << (t - 1)
            expected_rho = (t + 3) * (m + 3) - (t + 2)
            check(profile.end[x_t] == 5 * scale, "tail End formula")
            check(pair.n_cover == 32 * scale, "tail N formula")
            check(pair.b == 27 * scale - 4, "tail B formula")
            check(pair.rho == expected_rho, "tail rho formula")
            check(27 * pair.product - 25 * pair.b == 100, "tail 27/25 slack")
            check(pair.b - pair.product - pair.minimum == 2**t - 9, "tail +min formula")
            digest.update(
                f"family:{m}:{t}:{profile.h}:{profile.c}:{profile.end[x_t]}:"
                f"{pair.n_cover}:{pair.b}:{pair.rho}\n".encode()
            )
            checked += 1

        if m == 6:
            x6_profile = profile
            x6_pairs = all_pair_data(tournament, data, profile)
        if m == 7:
            x7_profile = profile
            x7_terminal = terminal

    check(checked == 21, "family check count")
    check(x6_profile is not None and x6_pairs is not None, "missing X6")
    check(x7_profile is not None and x7_terminal is not None, "missing X7")
    return x6_profile, x6_pairs, x7_profile, x7_terminal


def classified_two_reversal_shell(n: int) -> frozenset[tuple[tuple[int, int], ...]]:
    z = n - 1
    long_arc = (0, z)
    classified: set[tuple[tuple[int, int], ...]] = set()
    for other in all_arcs(n):
        if other not in (long_arc, (0, 1), (z - 1, z)):
            classified.add(tuple(sorted((long_arc, other))))
    for b in range(1, z):
        for c in range(1, b + 1):
            classified.add(tuple(sorted(((0, b), (c, z)))))
    return frozenset(classified)


def exact_two_reversal_order11(digest: "sha256") -> dict[str, object]:
    n = 11
    arcs = all_arcs(n)
    brute_strong: list[tuple[tuple[int, int], ...]] = []
    for reversal_pair in combinations(arcs, 2):
        tournament = from_reversals(n, reversal_pair)
        if is_strong(tournament):
            brute_strong.append(tuple(sorted(reversal_pair)))

    classified = classified_two_reversal_shell(n)
    check(frozenset(brute_strong) == classified, "two-reversal classification mismatch")
    check(len(brute_strong) == n * n - 2 * n - 2 == 97, "radius-two count")

    pair_checks = 0
    product_failures = 0
    plus_min_failures = 0
    candidate_failures = 0
    rho_repair_failures = 0
    five_copy_failures = 0
    five_copy_equalities = 0

    max_ratio_num = -1
    max_ratio_den = 1
    max_ratio_rows: list[tuple[str, int, int, PairData]] = []
    max_rho_num = -1
    max_rho_den = 1
    max_rho_rows: list[tuple[str, int, int, PairData]] = []
    min_gap: int | None = None
    min_gap_rows: list[tuple[str, Profile]] = []
    min_energy_num: int | None = None
    min_energy_den = 1
    min_energy_rows: list[tuple[str, Profile]] = []

    for reversals in sorted(brute_strong):
        tournament = from_reversals(n, reversals)
        word = encode(tournament)
        data = path_data(tournament)
        profile = tournament_profile(tournament, data)
        pairs = all_pair_data(tournament, data, profile)

        if profile.gap < 0:
            five_copy_failures += 1
        if profile.gap == 0:
            five_copy_equalities += 1
        if min_gap is None or profile.gap < min_gap:
            min_gap = profile.gap
            min_gap_rows = [(word, profile)]
        elif profile.gap == min_gap:
            min_gap_rows.append((word, profile))

        if min_energy_num is None or profile.target * min_energy_den < min_energy_num * profile.energy:
            min_energy_num = profile.target
            min_energy_den = profile.energy
            min_energy_rows = [(word, profile)]
        elif profile.target * min_energy_den == min_energy_num * profile.energy:
            min_energy_rows.append((word, profile))

        digest.update(
            f"shell:{word}:{profile.h}:{profile.c}:{profile.end}:{profile.defect}:"
            f"{profile.gap}:{profile.target}:{profile.energy}\n".encode()
        )
        for pair in pairs:
            pair_checks += 1
            if pair.b > pair.product:
                product_failures += 1
            if pair.b > pair.product + pair.minimum:
                plus_min_failures += 1
            if 25 * pair.b > 27 * pair.product:
                candidate_failures += 1
            if pair.b > pair.product + pair.rho:
                rho_repair_failures += 1

            if pair.b * max_ratio_den > max_ratio_num * pair.product:
                max_ratio_num = pair.b
                max_ratio_den = pair.product
                max_ratio_rows = [(word, pair.i, pair.j, pair)]
            elif pair.b * max_ratio_den == max_ratio_num * pair.product:
                max_ratio_rows.append((word, pair.i, pair.j, pair))

            excess = pair.b - pair.product
            if excess > 0:
                if excess * max_rho_den > max_rho_num * pair.rho:
                    max_rho_num = excess
                    max_rho_den = pair.rho
                    max_rho_rows = [(word, pair.i, pair.j, pair)]
                elif excess * max_rho_den == max_rho_num * pair.rho:
                    max_rho_rows.append((word, pair.i, pair.j, pair))

            digest.update(
                f"pair:{word}:{pair.i}:{pair.j}:{pair.n_cover}:{pair.b}:"
                f"{pair.product}:{pair.minimum}:{pair.rho}\n".encode()
            )

    check(pair_checks == 97 * 55 == 5335, "pair check count")
    check(five_copy_failures == 0 and five_copy_equalities == 0, "five-copy shell sign")
    check(candidate_failures == 0 and rho_repair_failures == 0, "pair candidate failure")
    check((product_failures, plus_min_failures) == (4, 3), "local failure counts")
    check(min_gap is not None and min_energy_num is not None, "missing shell extrema")

    return {
        "presentations": len(brute_strong),
        "candidate_pairs": len(tuple(combinations(arcs, 2))),
        "pair_checks": pair_checks,
        "product_failures": product_failures,
        "plus_min_failures": plus_min_failures,
        "candidate_failures": candidate_failures,
        "rho_repair_failures": rho_repair_failures,
        "five_copy_failures": five_copy_failures,
        "five_copy_equalities": five_copy_equalities,
        "min_gap": min_gap,
        "min_gap_rows": tuple(min_gap_rows),
        "min_energy_num": min_energy_num,
        "min_energy_den": min_energy_den,
        "min_energy_rows": tuple(min_energy_rows),
        "max_ratio_num": max_ratio_num,
        "max_ratio_den": max_ratio_den,
        "max_ratio_rows": tuple(max_ratio_rows),
        "max_rho_num": max_rho_num,
        "max_rho_den": max_rho_den,
        "max_rho_rows": tuple(max_rho_rows),
    }


def t_n_1(n: int) -> Tournament:
    z = n - 1
    out = [0] * n
    for i in range(z):
        for j in range(i + 1, z):
            out[i] |= 1 << j
    special = {0, n - 3}
    for i in range(z):
        if i in special:
            out[z] |= 1 << i
        else:
            out[i] |= 1 << z
    return Tournament(n, tuple(out))


def mapped_pair_dictionary(
    pairs: tuple[PairData, ...], image: tuple[int, ...]
) -> dict[tuple[int, int], tuple[int, int, int, int, int]]:
    result: dict[tuple[int, int], tuple[int, int, int, int, int]] = {}
    for pair in pairs:
        mapped = tuple(sorted((image[pair.i], image[pair.j])))
        result[mapped] = (pair.n_cover, pair.b, pair.product, pair.minimum, pair.rho)
    return result


def gauge_and_controls(
    shell: dict[str, object],
    x6_profile: Profile,
    x6_pairs: tuple[PairData, ...],
) -> dict[str, object]:
    n = 11
    x6 = x_family(6)
    x6_word = encode(x6)
    expected_x6_word = encode(from_reversals(n, ((0, 3), (3, 10))))
    check(x6_word == expected_x6_word, "X6 word")

    max_ratio_rows = shell["max_ratio_rows"]
    max_rho_rows = shell["max_rho_rows"]
    check(isinstance(max_ratio_rows, tuple) and len(max_ratio_rows) == 1, "local ratio tie")
    check(isinstance(max_rho_rows, tuple) and len(max_rho_rows) == 1, "rho ratio tie")
    max_word, max_i, max_j, max_pair = max_ratio_rows[0]
    rho_word, rho_i, rho_j, rho_pair = max_rho_rows[0]
    check((max_word, max_i, max_j) == (x6_word, 8, 10), "wrong local ratio hostile")
    check((max_pair.b, max_pair.product) == (428, 400), "wrong local ratio value")
    check((rho_word, rho_i, rho_j) == (x6_word, 9, 10), "wrong rho hostile")
    check((rho_pair.b - rho_pair.product, rho_pair.rho) == (53, 73), "wrong rho value")
    check((x6_profile.gap, x6_profile.target, x6_profile.energy) == (
        65_005_874,
        833_255_799,
        153_649_985,
    ), "X6 global profile")

    minimum_tournament = from_reversals(n, ((0, 8), (8, 10)))
    minimum_word = encode(minimum_tournament)
    min_gap_rows = shell["min_gap_rows"]
    min_energy_rows = shell["min_energy_rows"]
    check(isinstance(min_gap_rows, tuple) and len(min_gap_rows) == 3, "gap minimum orbit")
    check(isinstance(min_energy_rows, tuple) and len(min_energy_rows) == 3, "ratio minimum orbit")
    check(
        tuple(word for word, _profile in min_gap_rows)
        == tuple(word for word, _profile in min_energy_rows),
        "gap and ratio minimum presentations disagree",
    )
    check(minimum_word in {word for word, _profile in min_gap_rows}, "wrong global hostile")
    minimum_profile = next(profile for word, profile in min_gap_rows if word == minimum_word)
    check((minimum_profile.gap, minimum_profile.target, minimum_profile.energy) == (
        27_894_642,
        1_134_462_087,
        221_313_489,
    ), "T(11,1) profile")
    minimum_data = path_data(minimum_tournament)
    minimum_pairs = all_pair_data(minimum_tournament, minimum_data, minimum_profile)
    local_ratio = minimum_pairs[0]
    for pair in minimum_pairs[1:]:
        if pair.b * local_ratio.product > local_ratio.b * pair.product:
            local_ratio = pair
    check(local_ratio.b == local_ratio.product, "T(11,1) local maximum")

    # All three labelled minimum presentations are one T(11,1) orbit.  The
    # prefix 0..7 is fixed; an exact six-permutation check resolves the tail.
    tower = t_n_1(n)
    tower_maps: list[tuple[int, ...]] = []
    for word, _profile in min_gap_rows:
        reversals = []
        position = 0
        for arc in all_arcs(n):
            if word[position] == "0":
                reversals.append(arc)
            position += 1
        presentation = from_reversals(n, reversals)
        matching = []
        for tail_image in permutations((8, 9, 10)):
            image = tuple(range(8)) + tail_image
            if relabel(presentation, image) == tower:
                matching.append(image)
        check(len(matching) >= 1, "T(11,1) isomorphism")
        tower_maps.append(matching[0])

    # Swapping labels 0 and 1 changes the presentation radius but not the object.
    swap_zero_one = (1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    swapped_x6 = relabel(x6, swap_zero_one)
    swapped_reversals = reversal_set(swapped_x6)
    check(reversal_set(x6) == ((0, 3), (3, 10)), "X6 radius-two gauge")
    check(swapped_reversals == ((0, 1), (1, 3), (3, 10)), "X6 radius-three gauge")
    swapped_data = path_data(swapped_x6)
    swapped_profile = tournament_profile(swapped_x6, swapped_data)
    swapped_pairs = all_pair_data(swapped_x6, swapped_data, swapped_profile)
    expected_end = [0] * n
    expected_defect = [0] * n
    for old in range(n):
        expected_end[swap_zero_one[old]] = x6_profile.end[old]
        expected_defect[swap_zero_one[old]] = x6_profile.defect[old]
    check(swapped_profile.end == tuple(expected_end), "gauge End transport")
    check(swapped_profile.defect == tuple(expected_defect), "gauge defect transport")
    check(swapped_profile.gap == x6_profile.gap, "gauge gap transport")
    check(
        mapped_pair_dictionary(x6_pairs, swap_zero_one)
        == mapped_pair_dictionary(swapped_pairs, tuple(range(n))),
        "gauge pair transport",
    )

    c3 = from_reversals(3, ((0, 2),))
    c3_profile = tournament_profile(c3)
    check(is_strong(c3) and c3_profile.gap == 0, "C3 equality control")
    transitive = from_reversals(11, ())
    transitive_profile = tournament_profile(transitive)
    check(not is_strong(transitive) and transitive_profile.gap < 0, "transitive hostile control")

    return {
        "x6_word": x6_word,
        "x6_local": max_pair,
        "x6_rho": rho_pair,
        "minimum_word": minimum_word,
        "minimum_profile": minimum_profile,
        "minimum_local": local_ratio,
        "tower_maps": tuple(tower_maps),
        "swapped_word": encode(swapped_x6),
        "swapped_reversals": swapped_reversals,
        "c3_profile": c3_profile,
        "transitive_profile": transitive_profile,
    }


def main() -> None:
    digest = sha256()
    core_h, partition_sum, terminal_exclusions = formal_identity_checks()
    x6_profile, x6_pairs, x7_profile, x7_terminal = family_checks(digest)

    terminal_margins = tuple(
        (m, m * m + 5 * m + 18 - 2**m, m * m + 5 * m + 23 - 2**m)
        for m in range(1, 13)
    )
    check(all(row[1] >= 0 and row[2] >= 0 for row in terminal_margins[:6]), "rook m<=6")
    check(all(row[1] < 0 and row[2] < 0 for row in terminal_margins[6:]), "rook m>=7")
    check(
        (x7_terminal.b, x7_terminal.product, x7_terminal.minimum, x7_terminal.rho)
        == (1717, 1600, 5, 91),
        "X7 rook boundary",
    )

    shell = exact_two_reversal_order11(digest)
    controls = gauge_and_controls(shell, x6_profile, x6_pairs)

    print("tournament_order11_tail_fiber_thm4224 exact audit")
    print("status=SYMBOLIC_ALL_M_ALGEBRA_CHECKS,VERIFIED_EXACT_FAMILY,FINITE_EXACT_TWO_REVERSAL_SPHERE")
    print()
    print("[formal_tail_fiber]")
    print("scope=X_m,1<=t<m,u=2^(t-1)")
    print(
        f"literal_core_subset_DP=H:{core_h},partition_sum:{partition_sum},"
        f"terminal_exclusions:{terminal_exclusions}"
    )
    print("c=1 End_z=5 End_x_t=5u N_x_t,z=32u B_x_t,z=27u-4")
    print("product_excess=B-End_x_t*End_z=2u-4=2^t-4;positive_iff=t>=3")
    print("plus_min_excess=2u-9=2^t-9;positive_iff=t>=4")
    print("multiplicative_slack=(27/25)product-B=4;cross_multiplied_slack=100")
    print("terminal_boundary=N_r,z=32*2^(m-1)-7 B_r,z=27*2^(m-1)-11")
    print("terminal_cross_multiplied_slack=275")
    print("mechanism=nonempty_later_tail_forces_z_component_and_removes_seven_terminal_exclusions")
    print("symbolic_factor_and_sign_checks=PASS")
    print()
    print("[family_subset_dp]")
    print("universe=X_m,m=2..7,all_nonterminal_tail_pairs")
    print("checked_tournaments=6 checked_tail_pairs=21")
    print("checks=c,End,N,B,rho,plus_min,27_over_25,terminal_formulas")
    print("result=PASS")
    print()
    print("[successor_rook_boundary]")
    print("rho_r,z=m^2+5m+7 product_plus_rho_margin=m^2+5m+18-2^m")
    print("product_plus_min_plus_rho_margin=m^2+5m+23-2^m")
    print(
        "X6="
        f"B=853 product=800 min=5 rho=73 margins="
        f"{800 + 73 - 853},{800 + 5 + 73 - 853}"
    )
    print(
        "X7="
        f"B={x7_terminal.b} product={x7_terminal.product} min={x7_terminal.minimum} "
        f"rho={x7_terminal.rho} margins="
        f"{x7_terminal.product + x7_terminal.rho - x7_terminal.b},"
        f"{x7_terminal.product + x7_terminal.minimum + x7_terminal.rho - x7_terminal.b}"
    )
    print("threshold=both_repairs_hold_iff_m<=6;first_failure_m=7_order=12")
    print("all_m_sign_certificate=base_m7_plus_strictly_decreasing_forward_difference")
    print()
    print("[exact_two_reversal_order11]")
    print("gauge=fixed_labelled_transitive_order_0<...<10")
    print(f"candidate_reversal_pairs={shell['candidate_pairs']}")
    print(
        f"strong_presentations={shell['presentations']} formula=n^2-2n-2=97 "
        "classification_vs_reachability=PASS"
    )
    print(f"pair_checks={shell['pair_checks']}")
    print(
        "failures="
        f"five_copy:{shell['five_copy_failures']},"
        f"27_over_25:{shell['candidate_failures']},"
        f"product_plus_rho:{shell['rho_repair_failures']}"
    )
    print(
        "local_counts="
        f"product_failures:{shell['product_failures']},"
        f"plus_min_failures:{shell['plus_min_failures']}"
    )
    print(
        "max_B_over_product="
        f"{fraction_text(shell['max_ratio_num'], shell['max_ratio_den'])} "
        f"occurrences={len(shell['max_ratio_rows'])} word={controls['x6_word']} pair=8,10"
    )
    print(
        "max_positive_product_excess_over_rho="
        f"{fraction_text(shell['max_rho_num'], shell['max_rho_den'])} "
        f"occurrences={len(shell['max_rho_rows'])} word={controls['x6_word']} pair=9,10"
    )
    print(
        "min_five_copy_gap="
        f"{shell['min_gap']} occurrences={len(shell['min_gap_rows'])} "
        f"word={controls['minimum_word']}"
    )
    print(
        "min_five_copy_ratio="
        f"target:{shell['min_energy_num']},energy:{shell['min_energy_den']},"
        f"reduced:{fraction_text(shell['min_energy_num'], shell['min_energy_den'])} "
        f"occurrences={len(shell['min_energy_rows'])}"
    )
    print()
    print("[hostiles_and_gauge_controls]")
    print(
        "X6_local_hostile="
        f"B={controls['x6_local'].b} product={controls['x6_local'].product} "
        f"ratio={fraction_text(controls['x6_local'].b, controls['x6_local'].product)} "
        f"G={x6_profile.gap} target={x6_profile.target} energy={x6_profile.energy} "
        f"global_ratio={fraction_text(x6_profile.target, x6_profile.energy)}"
    )
    minimum_profile = controls["minimum_profile"]
    print(
        "T11_1_global_hostile="
        f"G={minimum_profile.gap} target={minimum_profile.target} energy={minimum_profile.energy} "
        f"ratio={fraction_text(minimum_profile.target, minimum_profile.energy)} "
        f"max_local_ratio={fraction_text(controls['minimum_local'].b, controls['minimum_local'].product)}"
    )
    print(f"T11_1_minimum_presentations={len(shell['min_gap_rows'])}")
    print(f"T11_1_isomorphism_maps_old_to_tower={controls['tower_maps']} check=PASS")
    print("X6_fixed_gauge_two_reversals=((0,3),(3,10))")
    print(f"X6_swap_0_1_fixed_gauge_three_reversals={controls['swapped_reversals']}")
    print(f"X6_swap_0_1_word={controls['swapped_word']}")
    print("gauge_transport=End,D,G,all_pair_N_B_product_min_rho:PASS")
    print(
        "controls="
        f"C3_strong_gap:{controls['c3_profile'].gap},"
        f"P11_strong:no,P11_gap:{controls['transitive_profile'].gap}"
    )
    print()
    print("[firewall]")
    print("exact_two_reversal_sphere_is_presentation_dependent_not_an_isomorphism_class_census")
    print("arbitrary_order11_strong_tournaments,27_over_25,HYP_9081,OS_plus=OPEN")
    print("raw_upper_zeta_and_normalized_rho_do_not_encode_weighted_successor_completions")
    print()
    print(f"audit_fingerprint_sha256={digest.hexdigest()}")
    print(f"script_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
    print("checks=PASS all_symbolic,family_DP,rook_threshold,shell_exact,hostiles,gauge,controls")


if __name__ == "__main__":
    main()
