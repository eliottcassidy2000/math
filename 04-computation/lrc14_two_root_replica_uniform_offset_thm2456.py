#!/usr/bin/env python3
"""Exact finite companion for THM-2456.

All truth-bearing arithmetic uses integers or fractions.Fraction.  There are
no Python asserts, floating-point approximations, or external packages.
"""

from __future__ import annotations

from fractions import Fraction as F


P = 13
SOURCE_ROWS = 7
Vector = tuple[F, ...]
Kernel = tuple[F, ...]
Chart = tuple[tuple[int, ...], ...]
Table = tuple[Vector, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def basis(index: int) -> Vector:
    require(0 <= index < P, "basis index")
    return tuple(F(int(position == index)) for position in range(P))


ZERO = tuple(F(0) for _ in range(P))
ONE = tuple(F(1) for _ in range(P))


def add(left: Vector, right: Vector) -> Vector:
    return tuple(left[index] + right[index] for index in range(P))


def scale(scalar: F, vector: Vector) -> Vector:
    return tuple(scalar * entry for entry in vector)


def apply_kernel(kernel: Kernel, vector: Vector) -> Vector:
    """(Kf)(s)=sum_r h(r-s)f(r) on C_13."""

    require(len(kernel) == P and len(vector) == P, "cyclic dimensions")
    return tuple(
        sum(
            (
                kernel[(root - shift) % P] * vector[root]
                for root in range(P)
            ),
            F(0),
        )
        for shift in range(P)
    )


def kernel_matrix(kernel: Kernel) -> tuple[Vector, ...]:
    return tuple(
        tuple(kernel[(column - row) % P] for column in range(P))
        for row in range(P)
    )


def transpose(matrix: tuple[Vector, ...]) -> tuple[Vector, ...]:
    return tuple(
        tuple(matrix[row][column] for row in range(P))
        for column in range(P)
    )


def multiply(
    left: tuple[Vector, ...],
    right: tuple[Vector, ...],
) -> tuple[Vector, ...]:
    return tuple(
        tuple(
            sum(
                (
                    left[row][index] * right[index][column]
                    for index in range(P)
                ),
                F(0),
            )
            for column in range(P)
        )
        for row in range(P)
    )


def matrix_rank(matrix: tuple[Vector, ...]) -> int:
    work = [list(row) for row in matrix]
    rank = 0
    for column in range(P):
        pivot = next(
            (row for row in range(rank, P) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [entry / pivot_value for entry in work[rank]]
        for row in range(P):
            if row == rank or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(P)
            ]
        rank += 1
    return rank


def danger_kernel(step: int) -> Kernel:
    require(1 <= step < P, "nonzero danger step")
    return tuple(F(int(index in (0, step))) for index in range(P))


def safe_kernel(step: int) -> Kernel:
    danger = danger_kernel(step)
    return tuple(F(1) - entry for entry in danger)


def full_support_kernel() -> Kernel:
    # J+I: positive on every root and invertible.
    return tuple(F(2 if index == 0 else 1) for index in range(P))


def table_from_density(kernel: Kernel, replica: Vector, offset: F) -> Table:
    require(offset > 0, "positive owner offset")
    owner = add(replica, scale(offset, ONE))
    replica_row = apply_kernel(kernel, replica)
    owner_row = apply_kernel(kernel, owner)
    return (owner_row,) + (replica_row,) * (SOURCE_ROWS - 1)


def anchored_offset(table: Table) -> F:
    require(len(table) == SOURCE_ROWS, "source row count")
    require(all(row[0] == 0 for row in table[1:]), "replica anchor")
    return table[0][0]


def check_delta_replicas(table: Table) -> tuple[F, Vector]:
    a = anchored_offset(table)
    require(a > 0, "positive owner anchor")
    replica = table[1]
    require(
        all(table[row] == replica for row in range(1, SOURCE_ROWS)),
        "six replicas differ",
    )
    require(
        all(table[0][shift] - replica[shift] == a for shift in range(P)),
        "owner-minus-replica is not uniform",
    )
    require(
        all(
            table[0][shift] - table[row][shift] - a == 0
            for row in range(1, SOURCE_ROWS)
            for shift in range(1, P)
        ),
        "anchored rectangle survived",
    )
    return a, replica


def chart() -> list[list[int]]:
    return [[0 for _ in range(P)] for _ in range(SOURCE_ROWS)]


def freeze_chart(mutable: list[list[int]]) -> Chart:
    frozen = tuple(tuple(row) for row in mutable)
    require(
        all(entry in (0, 1) for row in frozen for entry in row),
        "chart stopped being Boolean",
    )
    return frozen


def average_rows(charts: tuple[Chart, ...]) -> tuple[Vector, ...]:
    require(charts, "nonempty chart bank")
    denominator = F(len(charts))
    return tuple(
        tuple(
            sum(
                (F(item[row][root]) for item in charts),
                F(0),
            )
            / denominator
            for root in range(P)
        )
        for row in range(SOURCE_ROWS)
    )


def require_source_exclusive(charts: tuple[Chart, ...]) -> None:
    require(
        all(
            sum(item[row][root] for row in range(SOURCE_ROWS)) <= 1
            for item in charts
            for root in range(P)
        ),
        "source rows overlap on one chart/root",
    )


def two_chart_bank(root: int) -> tuple[Chart, ...]:
    first = chart()
    second = chart()
    first[0] = [1] * P
    second[0][root] = 1
    for row in range(1, SOURCE_ROWS):
        first[row][root] = 1
    return freeze_chart(first), freeze_chart(second)


def eight_chart_bank(root: int) -> tuple[Chart, ...]:
    charts = [chart() for _ in range(8)]
    charts[0][0] = [1] * P
    charts[1][0][root] = 1
    for row in range(1, SOURCE_ROWS):
        charts[row + 1][row][root] = 1
    frozen = tuple(freeze_chart(item) for item in charts)
    require_source_exclusive(frozen)
    return frozen


def singleton_chart_bank(root: int) -> tuple[Chart, ...]:
    charts: list[Chart] = []
    for position in range(P):
        item = chart()
        item[0][position] = 1
        charts.append(freeze_chart(item))
    extra_owner = chart()
    extra_owner[0][root] = 1
    charts.append(freeze_chart(extra_owner))
    for row in range(1, SOURCE_ROWS):
        item = chart()
        item[row][root] = 1
        charts.append(freeze_chart(item))
    frozen = tuple(charts)
    require(len(frozen) == 20, "singleton chart count")
    require_source_exclusive(frozen)
    require(
        all(sum(sum(row) for row in item) <= 1 for item in frozen),
        "singleton chart gained a second occupied pair",
    )
    return frozen


def table_from_chart_average(kernel: Kernel, charts: tuple[Chart, ...]) -> Table:
    rows = average_rows(charts)
    return tuple(apply_kernel(kernel, row) for row in rows)


def paired_safe_table(kernel: Kernel, root: int) -> tuple[Table, int]:
    entries: list[tuple[int, Chart]] = []
    hostile = eight_chart_bank(root)
    for target_offset in range(P):
        if target_offset == 0:
            bank = hostile
        else:
            mutable = [chart() for _ in range(8)]
            mutable[0][0] = [1] * P
            bank = tuple(freeze_chart(item) for item in mutable)
            require_source_exclusive(bank)
        entries.extend((target_offset, item) for item in bank)

    require(len(entries) == 104, "paired-safe chart count")
    base_safe = tuple(int(shift not in (5, 6)) for shift in range(P))
    table: list[Vector] = []
    for row in range(SOURCE_ROWS):
        response: list[F] = []
        for shift in range(P):
            total = F(0)
            for target_offset, item in entries:
                gate = base_safe[(shift - target_offset) % P]
                total += F(gate) * apply_kernel(
                    kernel,
                    tuple(F(value) for value in item[row]),
                )[shift]
            response.append(total / len(entries))
        table.append(tuple(response))
    return tuple(table), len(entries)


def row_sum(vector: Vector) -> F:
    return sum(vector, F(0))


def table_sum(table: Table) -> F:
    return sum((row_sum(row) for row in table), F(0))


# All ordinary nonzero steps have exact full rank for both truth values.
danger_ranks = tuple(
    matrix_rank(kernel_matrix(danger_kernel(step)))
    for step in range(1, P)
)
safe_ranks = tuple(
    matrix_rank(kernel_matrix(safe_kernel(step)))
    for step in range(1, P)
)
require(danger_ranks == (P,) * (P - 1), "danger kernel lost invertibility")
require(safe_ranks == (P,) * (P - 1), "safe kernel lost invertibility")

full_kernel = full_support_kernel()
require(all(entry > 0 for entry in full_kernel), "full-support control")
require(
    matrix_rank(kernel_matrix(full_kernel)) == P,
    "full-support control lost invertibility",
)

# The only positive constant difference of two Boolean values is one.
positive_boolean_differences = {
    F(owner - replica)
    for owner in (0, 1)
    for replica in (0, 1)
    if owner - replica > 0
}
require(
    positive_boolean_differences == {F(1)},
    "one-chart Boolean difference changed",
)
one_chart_nontrivial = sum(
    1
    for mask in range(1, 1 << P)
    if all(
        ((mask >> root) & 1) + 1 in (0, 1)
        for root in range(P)
    )
)
require(one_chart_nontrivial == 0, "nontrivial one-chart profile appeared")

danger = danger_kernel(1)
safe = safe_kernel(1)
root = 2

# Two charts: exact averaged hostile, deliberately not source-exclusive.
two_charts = two_chart_bank(root)
two_rows = average_rows(two_charts)
require(two_rows[1] == scale(F(1, 2), basis(root)), "two-chart replica density")
require(
    two_rows[0] == add(two_rows[1], scale(F(1, 2), ONE)),
    "two-chart owner density",
)
two_table = tuple(apply_kernel(danger, row) for row in two_rows)
two_a, two_w = check_delta_replicas(two_table)
require(two_a == 1, "two-chart anchor")
require(two_w[1] == F(1, 2) and two_w[2] == F(1, 2), "two-chart active shifts")
require(sum(value != 0 for value in two_w) == 2, "two-chart support")
require(row_sum(two_table[0]) == 14, "two-chart owner row sum")
require(row_sum(two_w) == 1, "two-chart replica row sum")
require(table_sum(two_table) == 20, "two-chart total mass")

# The safe truth value has the same exact hostile mechanism.
safe_two_root = 0
safe_rows = average_rows(two_chart_bank(safe_two_root))
safe_table = tuple(apply_kernel(safe, row) for row in safe_rows)
safe_a, safe_w = check_delta_replicas(safe_table)
require(safe_a == F(11, 2), "safe two-chart anchor")
require(sum(value != 0 for value in safe_w) == 11, "safe hostile support")

# Sharp source-exclusive eight-chart hostile.
eight_charts = eight_chart_bank(root)
eight_rows = average_rows(eight_charts)
require(
    eight_rows[1] == scale(F(1, 8), basis(root)),
    "eight-chart replica density",
)
require(
    eight_rows[0] == add(eight_rows[1], scale(F(1, 8), ONE)),
    "eight-chart owner density",
)
require(
    all(
        7 * eight_rows[1][position] + F(1, 8) <= 1
        for position in range(P)
    ),
    "7f+c exclusivity constraint",
)
eight_table = tuple(apply_kernel(danger, row) for row in eight_rows)
eight_a, eight_w = check_delta_replicas(eight_table)
require(eight_a == F(1, 4), "eight-chart anchor")
require(eight_w[1] == F(1, 8) and eight_w[2] == F(1, 8), "eight-chart shifts")
require(table_sum(eight_table) == 5, "eight-chart total mass")
minimum_equal_charts = 8
require(7 * F(1, minimum_equal_charts) + F(1, minimum_equal_charts) == 1, "sharp N=8")
require(
    7 * F(1, minimum_equal_charts - 1)
    + F(1, minimum_equal_charts - 1)
    > 1,
    "N=7 should fail",
)

# Singleton-only exact control.
singleton_charts = singleton_chart_bank(root)
singleton_rows = average_rows(singleton_charts)
require(
    singleton_rows[1] == scale(F(1, 20), basis(root)),
    "singleton replica density",
)
require(
    singleton_rows[0]
    == add(singleton_rows[1], scale(F(1, 20), ONE)),
    "singleton owner density",
)
singleton_table = tuple(apply_kernel(danger, row) for row in singleton_rows)
singleton_a, singleton_w = check_delta_replicas(singleton_table)
require(singleton_a == F(1, 10), "singleton anchor")
require(
    singleton_w[1] == F(1, 20) and singleton_w[2] == F(1, 20),
    "singleton active shifts",
)
require(table_sum(singleton_table) == 2, "singleton total mass")

# A translated two-hole safe partner still admits an exact hostile.
paired_table, paired_count = paired_safe_table(danger, root)
paired_a, paired_w = check_delta_replicas(paired_table)
require(paired_count == 104, "paired chart count")
require(paired_a == F(11, 52), "paired anchor")
require(
    paired_w[1] == F(1, 104) and paired_w[2] == F(1, 104),
    "paired active shifts",
)
require(sum(value != 0 for value in paired_w) == 2, "paired replica support")
require(paired_table[0][1] == F(23, 104), "paired owner active value")
require(
    all(
        paired_table[0][shift] == F(11, 52)
        for shift in range(P)
        if shift not in (1, 2)
    ),
    "paired owner background",
)
require(table_sum(paired_table) == F(75, 26), "paired total mass")

# Exact Gram identities underlying the analytic singular-value statement.
danger_matrix = kernel_matrix(danger)
safe_matrix = kernel_matrix(safe)
danger_gram = multiply(transpose(danger_matrix), danger_matrix)
safe_gram = multiply(transpose(safe_matrix), safe_matrix)
expected_danger_gram = tuple(
    tuple(
        F(
            2 * int(row == column)
            + int((column - row) % P == 1)
            + int((row - column) % P == 1)
        )
        for column in range(P)
    )
    for row in range(P)
)
require(danger_gram == expected_danger_gram, "danger Gram identity")
require(
    all(
        safe_gram[row][column] - danger_gram[row][column] == 9
        for row in range(P)
        for column in range(P)
    ),
    "safe Gram should differ by 9J",
)

print("THM-2456 exact companion")
print(f"group_order={P}")
print(f"source_rows={SOURCE_ROWS}")
print("danger_kernel_ranks=12_of_12_full")
print("safe_kernel_ranks=12_of_12_full")
print("full_support_control_rank=13")
print(f"one_chart_nontrivial_profiles={one_chart_nontrivial}")
print(f"two_chart_a={two_a}")
print("two_chart_w=1:1/2,2:1/2")
print(f"two_chart_row_sums={row_sum(two_table[0])},{row_sum(two_w)}")
print(f"two_chart_total_mass={table_sum(two_table)}")
print(f"safe_two_chart_a={safe_a}")
print(f"safe_two_chart_nonzero_w={sum(value != 0 for value in safe_w)}")
print("exclusive_constraint=7f+c<=1")
print(f"eight_chart_count={len(eight_charts)}")
print(f"eight_chart_a={eight_a}")
print("eight_chart_w=1:1/8,2:1/8")
print(f"eight_chart_total_mass={table_sum(eight_table)}")
print(f"eight_chart_minimal_equal_bank={minimum_equal_charts}")
print(f"singleton_chart_count={len(singleton_charts)}")
print(f"singleton_a={singleton_a}")
print(f"singleton_total_mass={table_sum(singleton_table)}")
print(f"paired_safe_chart_count={paired_count}")
print(f"paired_safe_a={paired_a}")
print("paired_safe_w=1:1/104,2:1/104")
print(f"paired_safe_total_mass={table_sum(paired_table)}")
print("danger_gram=2I+S+S^-1")
print("safe_gram_minus_danger_gram=9J")
print("ALL CHECKS PASSED")
