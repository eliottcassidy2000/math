#!/usr/bin/env python3
"""Exact THM-3710 endpoint / THM-3713 offset comparison.

On the thirty typed THM-3672 charts, reconstruct the two signed C_13
successor-endpoint fibres and the C_13 deep-offset defect.  Test every
affine reindexing of C_13 against the most general two-input convolution,
with and without a chart-independent constant term.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
P = 13
SOURCE = ROOT / "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py"

PINS = {
    "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py":
        "a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8",
    "05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out":
        "0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6",
    "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py":
        "0e4a9e181263647e13d2a6738b6996c45df901d9d2b37d4d589dfddfbdd91480",
}

EXPECTED_LEDGER_DIGEST = "6d740186d329b552298d6e8e38d319b61472c9d56f93de29c1edbf65e94d0715"
PRIMES = (1000003, 1000033)


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":")).encode("ascii")
    return sha256(encoded).hexdigest()


def load_source():
    spec = importlib.util.spec_from_file_location("thm3672_source", SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM-3672")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def periodic_arc_prefix(x, period, start, length):
    """Measure of [0,x) in one periodic half-open arc."""
    quotient, residue = divmod(x, period)
    result = quotient * length
    end = start + length
    pieces = ((start, end),) if end <= period else ((start, period), (0, end - period))
    for left, right in pieces:
        result += max(0, min(residue, right) - left)
    return result


def endpoint_counter(intervals, period):
    counter = defaultdict(int)
    for left, right in intervals:
        counter[left % period] -= 1
        counter[right % period] += 1
    return dict(counter)


def add_scaled(target, source, scalar):
    for point, value in source.items():
        target[point] += scalar * value


def endpoint_fibre(counter, period, sign):
    """Coefficients on {(14*j+sign)/182:j in F_13}."""
    require(period % 182 == 0, period)
    scale = period // 182
    return tuple(
        counter.get(((14 * j + sign) % 182) * scale, 0)
        for j in range(P)
    )


def reconstruct_ledgers():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    source = load_source()
    ref = source.load_referee()
    require(ref.W[ref.TARGET_B] == 742586, ref.W[ref.TARGET_B])
    require(ref.NN == 50334435734703120, ref.NN)

    zero = (0,) * 9
    word = ref.build_boolean_set(ref.PATTERN_QA, zero)

    def marked(shift):
        return source.marked_intervals(
            ref,
            ref.build_boolean_set(ref.PATTERN_E, shift),
            word,
        )

    intervals = {"zero": marked(zero)}
    for k in source.UNITS:
        intervals[f"a{k}"] = marked(source.negative_target_dipole(ref.TARGET_A, k))
        intervals[f"b{k}"] = marked(source.negative_target_dipole(ref.TARGET_B, k))

    deep_speed = ref.W[ref.TARGET_B]
    deep_period = ref.NN // deep_speed
    half_width = ref.NN // (14 * deep_speed)
    phase_unit = ref.NN // (P * deep_speed)
    require(deep_period == P * phase_unit, (deep_period, phase_unit))

    def colour_mass(marked_intervals, r):
        center = (r * phase_unit) % deep_period
        start = (center - half_width) % deep_period
        return sum(
            periodic_arc_prefix(right, deep_period, start, 2 * half_width)
            - periodic_arc_prefix(left, deep_period, start, 2 * half_width)
            for left, right in marked_intervals
        )

    colours = {
        label: tuple(colour_mass(marked_intervals, r) for r in range(P))
        for label, marked_intervals in intervals.items()
    }

    successor_period = ref.NN // (P * deep_speed)
    counters = {
        label: endpoint_counter(marked_intervals, successor_period)
        for label, marked_intervals in intervals.items()
    }

    rows = []
    for k in source.UNITS:
        for ell in source.UNITS:
            if k == ell:
                continue
            # THM-3713: u=r-t, so the t=1 term samples raw colour u+1.
            defect = tuple(
                colours["zero"][u]
                + colours[f"a{k}"][u]
                - 2 * colours[f"b{ell}"][(u + 1) % P]
                for u in range(P)
            )
            counter = defaultdict(int)
            add_scaled(counter, counters["zero"], 1)
            add_scaled(counter, counters[f"a{k}"], 1)
            add_scaled(counter, counters[f"b{ell}"], -2)
            plus = endpoint_fibre(counter, successor_period, 1)
            minus = endpoint_fibre(counter, successor_period, -1)
            rows.append((k, ell, defect, plus, minus))

    require(len(rows) == 30, len(rows))
    canonical = next(row for row in rows if row[:2] == (1, 2))
    require(
        tuple(value % P for value in canonical[3])
        == (6, 0, 0, 0, 4, 10, 10, 12, 5, 12, 12, 12, 12),
        canonical[3],
    )
    require(
        tuple(value % P for value in canonical[4])
        == (7, 6, 0, 0, 9, 3, 3, 1, 8, 1, 1, 1, 7),
        canonical[4],
    )

    ledger_digest = digest(rows)
    require(ledger_digest == EXPECTED_LEDGER_DIGEST, ("ledger digest", ledger_digest))
    return rows, ledger_digest


def rank_mod(matrix, prime):
    if not matrix:
        return 0
    a = [[value % prime for value in row] for row in matrix]
    row_count, column_count = len(a), len(a[0])
    rank = 0
    for column in range(column_count):
        pivot = next((r for r in range(rank, row_count) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][column], -1, prime)
        a[rank] = [(value * inverse) % prime for value in a[rank]]
        for r in range(row_count):
            if r != rank and a[r][column]:
                scalar = a[r][column]
                a[r] = [
                    (value - scalar * pivot_value) % prime
                    for value, pivot_value in zip(a[r], a[rank])
                ]
        rank += 1
        if rank == row_count:
            break
    return rank


def rank_q(matrix):
    if not matrix:
        return 0
    a = [[Fraction(value) for value in row] for row in matrix]
    row_count, column_count = len(a), len(a[0])
    rank = 0
    for column in range(column_count):
        pivot = next((r for r in range(rank, row_count) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pivot_value = a[rank][column]
        for r in range(rank + 1, row_count):
            if not a[r][column]:
                continue
            scalar = a[r][column] / pivot_value
            for j in range(column, column_count):
                a[r][j] -= scalar * a[rank][j]
        rank += 1
        if rank == row_count:
            break
    return rank


def is_prime(value):
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def convolution_ranks(rows, prime, alpha):
    """Ranks for a twisted two-fibre convolution, then with an intercept."""
    design = []
    target = []
    for _, _, defect, plus, minus in rows:
        for u in range(P):
            design.append(
                [plus[(alpha * u - j) % P] for j in range(P)]
                + [minus[(alpha * u - j) % P] for j in range(P)]
            )
            target.append(defect[u])

    linear_rank = rank_mod(design, prime)
    linear_augmented_rank = rank_mod(
        [row + [value] for row, value in zip(design, target)], prime
    )
    affine_design = [row + [1] for row in design]
    affine_rank = rank_mod(affine_design, prime)
    affine_augmented_rank = rank_mod(
        [row + [1, value] for row, value in zip(design, target)], prime
    )
    return linear_rank, linear_augmented_rank, affine_rank, affine_augmented_rank


def variation_ranks(rows):
    base = rows[0]
    endpoint_differences = []
    colour_differences = []
    for row in rows[1:]:
        endpoint = row[3] + row[4]
        base_endpoint = base[3] + base[4]
        endpoint_differences.append(
            [value - base_value for value, base_value in zip(endpoint, base_endpoint)]
        )
        colour_differences.append(
            [value - base_value for value, base_value in zip(row[2], base[2])]
        )
    return (
        rank_q(endpoint_differences),
        rank_q(colour_differences),
        rank_q(
            [endpoint + colour for endpoint, colour in zip(endpoint_differences, colour_differences)]
        ),
    )


def main():
    require(all(is_prime(prime) for prime in PRIMES), PRIMES)
    rows, ledger_digest = reconstruct_ledgers()

    colour_mod_13 = tuple(tuple(value % P for value in row[2]) for row in rows)
    endpoint_mod_13 = tuple(
        tuple(value % P for value in row[3] + row[4]) for row in rows
    )
    require(set(colour_mod_13) == {(0,) * P}, set(colour_mod_13))
    require(len(set(endpoint_mod_13)) == 13, len(set(endpoint_mod_13)))

    keys = tuple((row[0], row[1]) for row in rows)
    partitions_equal = all(
        (colour_mod_13[i] == colour_mod_13[j])
        == (endpoint_mod_13[i] == endpoint_mod_13[j])
        for i in range(len(keys))
        for j in range(len(keys))
    )
    require(not partitions_equal, partitions_equal)

    exact_variation_ranks = variation_ranks(rows)
    require(exact_variation_ranks == (6, 5, 6), exact_variation_ranks)

    rank_tables = {}
    for prime in PRIMES:
        table = tuple(
            convolution_ranks(rows, prime, alpha)
            for alpha in range(1, P)
        )
        require(table == ((26, 27, 27, 28),) * 12, (prime, table))
        rank_tables[prime] = table

    print(f"charts={len(rows)}")
    print(f"ledger_sha256={ledger_digest}")
    print("raw_offset_profiles_mod13=1")
    print("raw_offset_profile_is_zero_mod13=True")
    print("signed_endpoint_profiles_mod13=13")
    print("endpoint_offset_partitions_equal_mod13=False")
    print(f"exact_chart_variation_ranks_endpoint_offset_joint={exact_variation_ranks}")
    print("affine_reindexings_tested=12")
    for prime in PRIMES:
        common = rank_tables[prime][0]
        print(f"mod_{prime}_all_alpha_linear_rank_pair={common[:2]}")
        print(f"mod_{prime}_all_alpha_affine_rank_pair={common[2:]}")
    print("rational_twisted_convolution_exists=False")
    print("rational_affine_twisted_convolution_exists=False")


if __name__ == "__main__":
    main()
