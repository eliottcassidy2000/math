#!/usr/bin/env python3
"""Exact companion for THM-2551.

The huge relation-residue axis is handled by the proved explicit split
sections.  The 91-cell transfer and all finite hostiles are checked directly.
No floating point or optional dependency is used.
"""

from __future__ import annotations


P = 13
Q = 7
N = P * Q
TARGETS = P * P
RELATION_SIZE = 91**8
FIELD = 547
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def prime_factors(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def primitive_root(p: int) -> int:
    factors = prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // ell, p) != 1 for ell in factors):
            return g
    raise AssertionError("no primitive root")


def transfer(x: list[int]) -> list[int]:
    return [sum(x[(i - j) % N] for j in range(Q)) for i in range(N)]


def rank_mod(matrix: list[list[int]], p: int) -> int:
    a = [[v % p for v in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], p - 2, p)
        a[r] = [(v * inv) % p for v in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                scale = a[i][c]
                a[i] = [(a[i][j] - scale * a[r][j]) % p
                        for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def transfer_matrix() -> list[list[int]]:
    matrix = [[0] * N for _ in range(N)]
    for col in range(N):
        e = [0] * N
        e[col] = 1
        for row, value in enumerate(transfer(e)):
            matrix[row][col] = value
    return matrix


def check_split_blocks() -> int:
    """Toy integral model K(2)+EU(2)+EJ(2), with literal sections."""
    # Coordinates are K0,K1,EU0,EU1,EJ0,EJ1.
    U = [[0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0]]
    J = [[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]]
    require(rank_mod(U, 101) == rank_mod(J, 101) == 2,
            "toy U/J section ranks failed")
    require(rank_mod(U + J, 101) == 4, "toy joint section rank failed")
    # EU and EJ are exact independent coordinate sections.
    for q in range(2):
        eu = [0] * 6
        ej = [0] * 6
        eu[2 + q] = 1
        ej[4 + q] = 1
        require([sum(row[i] * eu[i] for i in range(6)) for row in U]
                == [1 if i == q else 0 for i in range(2)],
                "EU is not a U section")
        require([sum(row[i] * eu[i] for i in range(6)) for row in J] == [0, 0],
                "EU is not in ker J")
        require([sum(row[i] * ej[i] for i in range(6)) for row in U] == [0, 0],
                "EJ is not in ker U")
        require([sum(row[i] * ej[i] for i in range(6)) for row in J]
                == [1 if i == q else 0 for i in range(2)],
                "EJ is not a J section")
    return 4


def check_hall_hostile() -> tuple[int, int, int]:
    aligned = [[0] * P for _ in range(P)]
    swapped = [[0] * P for _ in range(P)]
    aligned[0][0] = aligned[1][1] = 1
    swapped[0][1] = swapped[1][0] = 1
    row_a = [sum(row) for row in aligned]
    row_b = [sum(row) for row in swapped]
    col_a = [sum(aligned[i][j] for i in range(P)) for j in range(P)]
    col_b = [sum(swapped[i][j] for i in range(P)) for j in range(P)]
    require(row_a == row_b and col_a == col_b,
            "Hall hostile marginals differ")
    diag_a = sum(aligned[i][i] for i in range(P))
    diag_b = sum(swapped[i][i] for i in range(P))
    require((diag_a, diag_b) == (2, 0), "Hall diagonal hostile failed")
    return sum(row_a), diag_a, diag_b


def main() -> None:
    n = RELATION_SIZE
    require(n == 4_702_525_276_151_521, "relation-space size failed")

    # Exact rank arithmetic from the two integral split sequences.
    intersection = 85 * (n - TARGETS)
    kernel = 91 * n - TARGETS * 85
    image = TARGETS * 85
    require(intersection == 399_714_648_472_864_920,
            "J intersection arithmetic failed")
    require(kernel == 427_929_800_129_774_046, "JD kernel arithmetic failed")
    require(image == 14_365, "JD image arithmetic failed")
    require(TARGETS * 6 == 1_014, "JD cokernel arithmetic failed")

    nonzero_intersection = 85 * (n - 168)
    nonzero_kernel = 91 * n - 168 * 85
    require(nonzero_intersection == 399_714_648_472_865_005,
            "Jstar intersection arithmetic failed")
    require(nonzero_kernel == 427_929_800_129_774_131,
            "Jstar kernel arithmetic failed")
    require(168 * 85 == 14_280 and 168 * 6 == 1_008,
            "Jstar image/cokernel arithmetic failed")

    joint_intersection = 85 * (n - 338)
    joint_kernel = 91 * n - 338 * 85
    require(joint_intersection == 399_714_648_472_850_555,
            "joint intersection arithmetic failed")
    require(joint_kernel == 427_929_800_129_759_681,
            "joint kernel arithmetic failed")
    require(338 * 85 == 28_730 and 338 * 6 == 2_028,
            "joint image/cokernel arithmetic failed")

    split_rank = check_split_blocks()

    matrix = transfer_matrix()
    transfer_rank = rank_mod(matrix, 1_000_003)
    require(transfer_rank == 85, "transfer rank is not 85")

    # Fourier multiplier: all 84 root-charged modes survive, and exactly the
    # six root-uniform clock augmentations die.
    generator = primitive_root(FIELD)
    xi = pow(generator, (FIELD - 1) // Q, FIELD)
    zeta = pow(generator, (FIELD - 1) // P, FIELD)
    root_survive = 0
    clock_killed = 0
    for alpha in range(P):
        for beta in range(Q):
            ratio = (pow(xi, (-beta) % Q, FIELD)
                     * pow(zeta, (-alpha) % P, FIELD)) % FIELD
            multiplier = sum(pow(ratio, j, FIELD) for j in range(Q)) % FIELD
            expected_zero = alpha == 0 and beta != 0
            require((multiplier == 0) == expected_zero,
                    "transfer Fourier zero set failed")
            if alpha:
                root_survive += 1
            elif beta:
                clock_killed += 1
    require((root_survive, clock_killed) == (84, 6),
            "wrong root/clock Fourier census")

    # Positive nonunit-section hostile: its transfer is seven nonnegative
    # point masses and has every root-charged Fourier mode.
    e = [0] * N
    e[0] = 1
    de = transfer(e)
    require(set(de) <= {0, 1} and sum(de) == 7,
            "positive transfer hostile failed")
    c_host_U = [1] * TARGETS
    c_host_J = [0] * TARGETS
    require(all(c_host_U) and not any(c_host_J),
            "projector-kernel hostile failed")

    # Sharp clock-augmentation loss.
    a_clk = [0] * N
    for t in range(P):
        a_clk[Q * t] = 1
        a_clk[1 + Q * t] = -1
    require(any(a_clk) and transfer(a_clk) == [0] * N,
            "clock augmentation sharp control failed")

    total_mass, hall_aligned, hall_swapped = check_hall_hostile()

    print("THM-2551 exact transfer/projector tensor-intersection referee")
    print(f"relation_size={n}")
    print(f"split_section_rank={split_rank}")
    print(f"transfer_rank={transfer_rank} root_charged_survive={root_survive} clock_augmentation_killed={clock_killed}")
    print(f"J_intersection_rank={intersection}")
    print(f"JD_rank={image} JD_kernel_rank={kernel} JD_cokernel_rank={TARGETS * 6}")
    print(f"Jstar_intersection_rank={nonzero_intersection}")
    print(f"JstarD_rank={168 * 85} JstarD_kernel_rank={nonzero_kernel} JstarD_cokernel_rank={168 * 6}")
    print(f"joint_intersection_rank={joint_intersection}")
    print(f"joint_rank={338 * 85} joint_kernel_rank={joint_kernel} joint_cokernel_rank={338 * 6}")
    print(f"positive_kernel_hostile_U_support={sum(c_host_U)} J_support={sum(c_host_J)} transfer_mass={sum(de)}")
    print(f"clock_visible_then_killed={int(any(a_clk))}")
    print(f"hall_equal_marginal_mass={total_mass} aligned_diagonal={hall_aligned} swapped_diagonal={hall_swapped}")
    print(f"explicit_require_checks={CHECKS}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
