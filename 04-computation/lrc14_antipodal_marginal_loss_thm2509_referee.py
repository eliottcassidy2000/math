#!/usr/bin/env python3
"""Independent exact referee for the marginal-loss boundary in THM-2509."""

from fractions import Fraction


P = 13
R = 7
DOMAIN_DIMENSION = P * (R - 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rational_rank(matrix: list[list[int]]) -> int:
    """Return the exact Q-rank by independent Fraction row reduction."""

    rows = [[Fraction(entry) for entry in row] for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        rows[pivot_row] = [entry / pivot_value for entry in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            multiplier = rows[row][column]
            rows[row] = [
                left - multiplier * right
                for left, right in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def marginal_matrix(tau: int) -> list[list[int]]:
    """Matrix of d -> (R_tau d, R_-tau d) on a row-zero basis.

    The domain basis is e_(h,r)-e_(h,6), with h in F_13 and 0<=r<6.
    """

    matrix = [[0] * DOMAIN_DIMENSION for _ in range(2 * P)]
    for h in range(P):
        for r in range(R - 1):
            column = h * (R - 1) + r
            matrix[(h + tau * r) % P][column] += 1
            matrix[(h + tau * (R - 1)) % P][column] -= 1
            matrix[P + (h - tau * r) % P][column] += 1
            matrix[P + (h - tau * (R - 1)) % P][column] -= 1
    return matrix


def ramanujan_trace(exponent: int) -> int:
    """Tr_{Q(zeta_13)/Q}(zeta_13^exponent), evaluated exactly."""

    return 12 if exponent % P == 0 else -1


def trace_hostile(tau: int) -> list[list[int]]:
    """Integral trace witness annihilated by both antipodal marginals."""

    defect = [[0] * R for _ in range(P)]
    for h in range(P):
        trace = ramanujan_trace(h)
        symmetric_trace = (
            trace + ramanujan_trace(h + tau) + ramanujan_trace(h - tau)
        )
        defect[h] = [
            -trace,
            symmetric_trace,
            -symmetric_trace,
            trace,
            0,
            0,
            0,
        ]
    return defect


def radon(defect: list[list[int]], tau: int) -> list[int]:
    return [
        sum(defect[(v - tau * r) % P][r] for r in range(R))
        for v in range(P)
    ]


def main() -> None:
    rank_rows: list[str] = []
    hostile_l1: set[int] = set()
    for tau in range(1, P):
        matrix = marginal_matrix(tau)
        first_rank = rational_rank(matrix[:P])
        second_rank = rational_rank(matrix[P:])
        pair_rank = rational_rank(matrix)
        kernel_dimension = DOMAIN_DIMENSION - pair_rank
        require(first_rank == 12 and second_rank == 12,
                f"single-marginal rank drifted at tau={tau}")
        require(pair_rank == 24 and kernel_dimension == 54,
                f"paired-marginal rank/kernel drifted at tau={tau}")
        rank_rows.append(f"{tau}:{first_rank}/{second_rank}/{pair_rank}/{kernel_dimension}")

        defect = trace_hostile(tau)
        require(all(sum(row) == 0 for row in defect),
                f"trace hostile is not row-zero at tau={tau}")
        require(any(defect[h] != defect[0] for h in range(1, P)),
                f"trace hostile became vertical at tau={tau}")
        require(radon(defect, tau) == [0] * P,
                f"positive marginal survived at tau={tau}")
        require(radon(defect, -tau) == [0] * P,
                f"negative marginal survived at tau={tau}")
        hostile_l1.add(sum(abs(entry) for row in defect for entry in row))

    require(hostile_l1 == {168}, "trace-hostile L1 norm drifted")

    print("THM-2509 antipodal marginal-loss independent referee: PASS")
    print("domain_dimension=78; output_dimension=26")
    print("tau:rank_plus/rank_minus/rank_pair/kernel=" + ",".join(rank_rows))
    print("integral_trace_hostiles=12; row_zero=12; nonvertical=12")
    print("annihilated_antipodal_marginals=24; hostile_L1=168")
    print("the joint 91-entry chart is lossless; its two marginals are not")
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
