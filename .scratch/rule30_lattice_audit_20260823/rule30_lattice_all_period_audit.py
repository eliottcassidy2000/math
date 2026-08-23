#!/usr/bin/env python3
"""Independent exact audit of the Rule 30 amplitude lattice.

For n >= 1 let

    (T_n x)_j = x_(2j mod n) + x_(2j+1 mod n).

The companion report proves the image and Smith form of every power T_n^r,
including even n.  This program supplies independent exact controls.  It does
not import the earlier scout and deliberately uses explicit RuntimeError gates
instead of Python assertions, so every check remains active under ``-O``.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def two_adic_order(n: int) -> int:
    require(n != 0, "two_adic_order is undefined at zero")
    n = abs(n)
    out = 0
    while n % 2 == 0:
        out += 1
        n //= 2
    return out


def period_parameters(n: int, r: int) -> tuple[int, int, int]:
    """Return (a,d,s), with n=2^a*m, d=n/2^min(a,r), s=max(0,r-a)."""
    require(n >= 1 and r >= 0, "invalid period or iterate")
    odd_part = n
    a = 0
    while odd_part % 2 == 0:
        a += 1
        odd_part //= 2
    peeled = min(a, r)
    return a, n >> peeled, r - peeled


def matrix_t(n: int) -> Matrix:
    require(n >= 1, "matrix_t requires positive length")
    rows = [[0 for _ in range(n)] for _ in range(n)]
    for j in range(n):
        rows[j][(2 * j) % n] += 1
        rows[j][(2 * j + 1) % n] += 1
    return Matrix(rows)


def forward_step(x: tuple[int, ...]) -> tuple[int, ...]:
    n = len(x)
    require(n >= 1, "empty vector")
    return tuple(x[(2 * j) % n] + x[(2 * j + 1) % n] for j in range(n))


def forward_iter(x: tuple[int, ...], r: int) -> tuple[int, ...]:
    require(r >= 0, "negative iterate")
    out = x
    for _ in range(r):
        out = forward_step(out)
    return out


def repeat_block(block: tuple[int, ...], n: int) -> tuple[int, ...]:
    require(block and n % len(block) == 0, "invalid repetition")
    return tuple(block[j % len(block)] for j in range(n))


def predicted_member(target: tuple[int, ...], r: int) -> bool:
    n = len(target)
    _, d, s = period_parameters(n, r)
    if any(target[j] != target[j % d] for j in range(n)):
        return False
    return sum(target[:d]) % (1 << s) == 0


def inverse_odd_step(target: tuple[int, ...]) -> tuple[int, ...] | None:
    """Unique rational inverse of T_n; integral exactly at even total sum."""
    n = len(target)
    require(n % 2 == 1, "inverse_odd_step requires odd length")
    inverse_two = pow(2, -1, n) if n > 1 else 0
    reordered = tuple(target[(inverse_two * k) % n] for k in range(n))
    numerators = tuple(
        sum((1 if t % 2 == 0 else -1) * reordered[(k + t) % n] for t in range(n))
        for k in range(n)
    )
    parity = sum(target) & 1
    require(all((value & 1) == parity for value in numerators), "inverse parity identity")
    if parity:
        return None
    source = tuple(value // 2 for value in numerators)
    require(forward_step(source) == target, "odd inverse formula")
    return source


def constructive_preimage(target: tuple[int, ...], r: int) -> tuple[int, ...] | None:
    """Construct a T_n^r-preimage exactly on the predicted all-period image."""
    require(r >= 0 and target, "invalid constructive_preimage input")
    if r == 0:
        return target
    n = len(target)
    if n % 2 == 1:
        prior = inverse_odd_step(target)
        if prior is None:
            return None
        return constructive_preimage(prior, r - 1)

    half = n // 2
    if any(target[j + half] != target[j] for j in range(half)):
        return None
    reduced_source = constructive_preimage(target[:half], r - 1)
    if reduced_source is None:
        return None
    source = [0 for _ in range(n)]
    for j, value in enumerate(reduced_source):
        source[2 * j] = value
    out = tuple(source)
    require(forward_iter(out, r) == target, "even recursive section")
    return out


def predicted_smith(n: int, r: int) -> tuple[int, ...]:
    _, d, s = period_parameters(n, r)
    return (1,) * (d - 1) + (1 << s,) + (0,) * (n - d)


def exact_smith(n: int, r: int) -> tuple[int, ...]:
    diagonal = smith_normal_form(matrix_t(n) ** r, domain=ZZ)
    return tuple(abs(int(diagonal[j, j])) for j in range(n))


def permutation_matrix(permutation: tuple[int, ...]) -> Matrix:
    n = len(permutation)
    require(sorted(permutation) == list(range(n)), "not a permutation")
    return Matrix([[int(permutation[i] == j) for j in range(n)] for i in range(n)])


def cyclic_shift_matrix(n: int) -> Matrix:
    return Matrix([[int(j == (i + 1) % n) for j in range(n)] for i in range(n)])


def projective_step(profile: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    n = len(profile)
    require(n >= 1, "empty projective profile")
    out = []
    for j in range(n):
        g0 = profile[(2 * j) % n]
        g1 = profile[(2 * j + 1) % n]
        g2 = profile[(2 * j + 2) % n]
        require(g0 != 1, "projective denominator boundary")
        out.append(-g0 * g1 * (1 - g2) / (1 - g0))
    return tuple(out)


def projective_step_mod(profile: tuple[int, ...], prime: int) -> tuple[int, ...] | None:
    out = []
    n = len(profile)
    for j in range(n):
        g0 = profile[(2 * j) % n]
        if g0 == 1:
            return None
        g1 = profile[(2 * j + 1) % n]
        g2 = profile[(2 * j + 2) % n]
        out.append((-g0 * g1 * (1 - g2) * pow(1 - g0, -1, prime)) % prime)
    return tuple(out)


def is_periodic(profile: tuple[object, ...], period: int) -> bool:
    return all(profile[j] == profile[j % period] for j in range(len(profile)))


def proportional(x: tuple[int, ...], y: tuple[int, ...]) -> bool:
    require(len(x) == len(y) and x and y, "invalid proportionality test")
    return all(x[0] * y[j] == y[0] * x[j] for j in range(len(x)))


def main() -> None:
    smith_checks = 0
    power_factorization_checks = 0
    kernel_checks = 0
    lift_checks = 0
    image_hostiles = 0
    factor_checks = 0
    permutation_checks = 0
    projective_checks = 0
    finite_field_checks = 0
    physical_checks = 0
    semantic_records: list[dict[str, object]] = []

    # Exact all-period Smith controls, using SymPy's independent integer SNF.
    for n in range(1, 25):
        for r in range(0, 9):
            _, d, s = period_parameters(n, r)
            block_size = n // d
            power = matrix_t(n) ** r
            got = exact_smith(n, r)
            want = predicted_smith(n, r)
            require(got == want, f"Smith mismatch n={n}, r={r}: {got} != {want}")
            block_sum = Matrix(
                [[int(block_size * i <= j < block_size * (i + 1)) for j in range(n)]
                 for i in range(d)]
            )
            repeat = Matrix([[int(j == i % d) for j in range(d)] for i in range(n)])
            factor = repeat * (matrix_t(d) ** s) * block_sum
            require(power == factor, f"power factorization n={n}, r={r}")
            require(power.rank() == d, f"power rank n={n}, r={r}")
            power_factorization_checks += 2
            for block in range(d):
                start = block_size * block
                for offset in range(1, block_size):
                    vector = Matrix(
                        [int(j == start + offset) - int(j == start) for j in range(n)]
                    )
                    require(power * vector == Matrix.zeros(n, 1), "block-difference kernel")
                    kernel_checks += 1
            smith_checks += 1

    # Constructive membership in both directions, including free and torsion
    # hostiles on the even-period boundary.
    for n in range(1, 33):
        for r in range(0, 11):
            _, d, s = period_parameters(n, r)
            for seed in range(3):
                block = [((seed + 2) * (j + 1) ** 3 - 5 * j + n + r) for j in range(d)]
                block[-1] += (-sum(block)) % (1 << s)
                target = repeat_block(tuple(block), n)
                require(predicted_member(target, r), "constructed target rejected")
                source = constructive_preimage(target, r)
                require(source is not None, "constructive image target did not lift")
                require(forward_iter(source, r) == target, "constructive round trip")
                lift_checks += 3

            if d < n:
                nonperiodic = (1,) + (0,) * (n - 1)
                require(not predicted_member(nonperiodic, r), "free hostile accepted")
                require(constructive_preimage(nonperiodic, r) is None, "free hostile lifted")
                image_hostiles += 2
            if s > 0:
                bad_block = (1,) + (0,) * (d - 1)
                bad_torsion = repeat_block(bad_block, n)
                require(not predicted_member(bad_torsion, r), "torsion hostile accepted")
                require(constructive_preimage(bad_torsion, r) is None, "torsion hostile lifted")
                image_hostiles += 2

            semantic_records.append(
                {"n": n, "r": r, "period_rank": d, "torsion_exponent": 1 << s}
            )

    # Verify T_(2m)=E T_m Q and the full rank-m kernel generated by pair
    # differences, rather than only the single alternating hostile.
    for n in range(2, 33, 2):
        half = n // 2
        t_n = matrix_t(n)
        t_half = matrix_t(half)
        repeat = Matrix([[int(j == i % half) for j in range(half)] for i in range(n)])
        pair_sum = Matrix([[int(j in (2 * i, 2 * i + 1)) for j in range(n)] for i in range(half)])
        require(t_n == repeat * pair_sum, "even T=E Q factorization")
        require(t_n * repeat == repeat * t_half, "even intertwining")
        require(t_n.rank() == half, "even rank")
        for j in range(half):
            kernel_vector = Matrix([int(k == 2 * j) - int(k == 2 * j + 1) for k in range(n)])
            require(t_n * kernel_vector == Matrix.zeros(n, 1), "pair kernel")
            factor_checks += 1
        factor_checks += 3

    # Odd-period robustness under arbitrary deterministic coordinate
    # permutations between (I+S) factors.
    for n in range(1, 18, 2):
        shift = cyclic_shift_matrix(n)
        product = Matrix.eye(n)
        for r in range(1, 7):
            permutation = tuple((3 * j + r) % n for j in range(n))
            if len(set(permutation)) != n:
                permutation = tuple((j + r) % n for j in range(n))
            product = permutation_matrix(permutation) * (Matrix.eye(n) + shift) * product
            diagonal = smith_normal_form(product, domain=ZZ)
            got = tuple(abs(int(diagonal[j, j])) for j in range(n))
            want = (1,) * (n - 1) + (1 << r,)
            require(got == want, "permuted odd Smith control")
            require(all(sum(int(product[i, j]) for i in range(n)) == (1 << r) for j in range(n)),
                    "permuted column-sum law")
            permutation_checks += 2

    # Direct projective period-halving controls over Q.
    for n in range(2, 25, 2):
        profile = tuple(Fraction(j + 2, j + n + 3) for j in range(n))
        image = projective_step(profile)
        require(is_periodic(image, n // 2), "even projective period did not halve")
        projective_checks += 1

    # Exhaustive finite-field cycle hostile: any defined cycle represented
    # with an even declared period already repeats on its odd core.
    for prime in (5, 7):
        for n in (2, 4, 6):
            odd_part = n >> two_adic_order(n)
            values = tuple(range(2, prime))
            for profile in itertools.product(values, repeat=n):
                current = tuple(profile)
                for _ in range(3):
                    nxt = projective_step_mod(current, prime)
                    if nxt is None or any(value in (0, 1) for value in nxt):
                        break
                    current = nxt
                    finite_field_checks += 1
                    if current == profile:
                        require(is_periodic(tuple(profile), odd_part), "genuine even cycle survived")

    # The all-odd physical constant-mode obstruction: for every odd unit c,
    # the next constant ratio -c^2 has gap exactly one.  Repeating once gives
    # two consecutive unit gaps.
    for precision in range(3, 13):
        for c in range(1, 1 << precision, 2):
            require(two_adic_order(1 + c * c) == 1, "odd-square gap law")
            physical_checks += 1

    # Positive and hostile controls delimiting what the carry sees.
    one_step_source = (1, 5, 9)
    one_step_raw = forward_step(one_step_source)
    require(tuple(two_adic_order(value) for value in one_step_raw) == (1, 1, 1),
            "uniform one-step physical normalization")
    one_step_parent = tuple(value // 2 for value in one_step_raw)
    require(all(value & 1 for value in one_step_parent), "normalized parent units")
    require(sum(one_step_source) == sum(one_step_parent), "one-step sum invoice")
    require(not proportional(one_step_source, one_step_parent), "nonclosing control closed")
    plane_hostile = (1, 1, -2)
    require(forward_iter(plane_hostile, 2) == plane_hostile, "period-three plane hostile")
    require(sum(plane_hostile) == 0 and any(value % 2 == 0 for value in plane_hostile),
            "plane hostile typing")
    physical_checks += 6

    digest = hashlib.sha256(
        "\n".join(json.dumps(record, sort_keys=True) for record in semantic_records).encode("ascii")
    ).hexdigest()

    print("RULE30 ALL-PERIOD AMPLITUDE-LATTICE AUDIT")
    print("THEOREM n=2^a*m, d=n/2^min(a,r), s=r-min(a,r)")
    print("IMAGE period-d vectors whose primitive-block sum is 0 mod 2^s")
    print("SMITH diag(1^(d-1),2^s,0^(n-d))")
    print("COKERNEL Z^(n-d) direct-sum Z/(2^s), with Z/(1)=0")
    print("SMITH_ACTIVE_CHECKS", smith_checks)
    print("POWER_FACTORIZATION_ACTIVE_CHECKS", power_factorization_checks)
    print("KERNEL_BASIS_ACTIVE_CHECKS", kernel_checks)
    print("CONSTRUCTIVE_LIFT_ACTIVE_CHECKS", lift_checks)
    print("IMAGE_HOSTILE_ACTIVE_CHECKS", image_hostiles)
    print("EVEN_FACTORIZATION_ACTIVE_CHECKS", factor_checks)
    print("PERMUTATION_COVARIANCE_ACTIVE_CHECKS", permutation_checks)
    print("PROJECTIVE_HALVING_ACTIVE_CHECKS", projective_checks)
    print("FINITE_FIELD_TRANSITION_ACTIVE_CHECKS", finite_field_checks)
    print("PHYSICAL_GAP_ACTIVE_CHECKS", physical_checks)
    print("SEMANTIC_RECORDS", len(semantic_records))
    print("SEMANTIC_SHA256", digest)
    print("SMITH_UNIVERSE", "1<=n<=24, 0<=r<=8")
    print("LIFT_UNIVERSE", "1<=n<=32, 0<=r<=10, three targets each")
    print("EVEN_BOUNDARY", "rank halves until odd core; later one cyclic 2-carry")
    print("NONCLOSING_CONTROL", "A=(1,5,9) normalizes to (3,5,7), same sum, no projective return")
    print("CARRY_ONLY_HOSTILE", "n=3,r=2,A=(1,1,-2) is fixed but not all odd")
    print("DYADIC_PROFILE_CONSEQUENCE", "a physical finite spatial period cannot divide a power of two")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
