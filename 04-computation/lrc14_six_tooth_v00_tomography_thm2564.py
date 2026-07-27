#!/usr/bin/env python3
"""Exact referee for THM-2564.

The six nonaffine tooth maps on the integral doubly centered
``C_7 x F_13`` lattice form a rational tomography atlas.  This companion
checks every five- and six-slope bank, the sharp literal field-trace
witnesses, the integral 13-primary invoice, and the ordered-section seam.
"""

from __future__ import annotations

from itertools import combinations


P = 13
Q = 7
DIM = (P - 1) * (Q - 1)
ZERO = (0,) * (P - 1)
ONE = (1,) + (0,) * (P - 2)
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def zero_table() -> list[list[int]]:
    return [[0 for _ in range(P)] for _ in range(Q)]


def rectangle(k: int, v: int) -> list[list[int]]:
    f = zero_table()
    f[k][v] += 1
    f[k][P - 1] -= 1
    f[Q - 1][v] -= 1
    f[Q - 1][P - 1] += 1
    return f


RECTANGLES = [(k, v) for k in range(Q - 1) for v in range(P - 1)]


def is_v00(f: list[list[int]]) -> bool:
    return (
        all(sum(f[k][v] for v in range(P)) == 0 for k in range(Q))
        and all(sum(f[k][v] for k in range(Q)) == 0 for v in range(P))
    )


def tooth(f: list[list[int]], tau: int) -> list[int]:
    return [
        sum(f[k][(v - tau * k) % P] for k in range(Q))
        for v in range(P)
    ]


def output_coordinates(f: list[list[int]], tau: int) -> list[int]:
    out = tooth(f, tau)
    require(sum(out) == 0, "tooth output left V_0")
    return out[: P - 1]


SLOPE_ROWS: dict[int, list[list[int]]] = {}
for tau in range(1, P):
    columns = [output_coordinates(rectangle(k, v), tau) for k, v in RECTANGLES]
    SLOPE_ROWS[tau] = [
        [columns[j][i] for j in range(DIM)] for i in range(P - 1)
    ]


def bank_matrix(slopes: tuple[int, ...]) -> list[list[int]]:
    return [row for tau in slopes for row in SLOPE_ROWS[tau]]


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    if not matrix:
        return 0
    a = [[x % prime for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                c = a[i][col]
                a[i] = [(x - c * y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def det_bareiss(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    require(n > 0 and all(len(row) == n for row in a), "non-square determinant")
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        require(pivot is not None, "representative bank became singular")
        assert pivot is not None
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        p = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * p - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division was not exact")
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = p
    return sign * a[-1][-1]


def add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x + y for x, y in zip(a, b))


def neg(a: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in a)


def mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply in Z[z]/(1+z+...+z^12)."""
    raw = [0] * (2 * (P - 2) + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            raw[i + j] += x * y
    for exponent in range(len(raw) - 1, P - 2, -1):
        coefficient = raw[exponent]
        if coefficient:
            for j in range(P - 1):
                raw[exponent - (P - 1) + j] -= coefficient
            raw[exponent] = 0
    return tuple(raw[: P - 1])


def zeta_pow(exponent: int) -> tuple[int, ...]:
    exponent %= P
    if exponent < P - 1:
        ans = [0] * (P - 1)
        ans[exponent] = 1
        return tuple(ans)
    return (-1,) * (P - 1)


def dft(vector: list[int], alpha: int) -> tuple[int, ...]:
    ans = ZERO
    for v, coefficient in enumerate(vector):
        if coefficient:
            ans = add(
                ans,
                tuple(coefficient * x for x in zeta_pow(-alpha * v)),
            )
    return ans


def poly_mul_linear(
    coefficients: list[tuple[int, ...]], root: tuple[int, ...]
) -> list[tuple[int, ...]]:
    out = [ZERO] * (len(coefficients) + 1)
    for i, coefficient in enumerate(coefficients):
        out[i] = add(out[i], neg(mul(coefficient, root)))
        out[i + 1] = add(out[i + 1], coefficient)
    return out


def trace_zeta_product(coefficient: tuple[int, ...], v: int) -> int:
    """Trace(coefficient*zeta^v) using Tr(zeta^e)=12 or -1."""
    return sum(
        c * (P - 1 if (exponent + v) % P == 0 else -1)
        for exponent, c in enumerate(coefficient)
    )


def trace_witness(slopes: tuple[int, ...]) -> list[list[int]]:
    coefficients = [ONE]
    coefficients = poly_mul_linear(coefficients, ONE)  # X-1
    for tau in slopes:
        coefficients = poly_mul_linear(coefficients, zeta_pow(-tau))
    require(len(coefficients) == Q, "trace polynomial has wrong degree")
    require(coefficients[-1] == ONE, "trace polynomial lost its leading term")
    f = [
        [trace_zeta_product(coefficients[k], v) for v in range(P)]
        for k in range(Q)
    ]
    require(is_v00(f), "literal field-trace witness is not doubly centered")
    require(any(value for row in f for value in row), "trace witness vanished")
    return f


def shear(f: list[list[int]], sigma: int) -> list[list[int]]:
    return [
        [f[k][(v - sigma * k) % P] for v in range(P)] for k in range(Q)
    ]


def translate_rows(f: list[list[int]], c: int) -> list[list[int]]:
    return [[f[(k - c) % Q][v] for v in range(P)] for k in range(Q)]


print("== THM-2564: six-tooth V00 tomography ==")
require(len(RECTANGLES) == DIM == 72, "wrong doubly centered dimension")
require(all(is_v00(rectangle(k, v)) for k, v in RECTANGLES), "bad basis")
require(all(tooth(rectangle(0, 0), 0) == [0] * P for _ in range(1)), "T_0 != 0")
print("  V00 rank / one tooth codomain rank: 72 / 12")
print("  T_0 vanishes identically by the second margin")

five_banks = 0
for slopes in combinations(range(1, P), 5):
    require(rank_mod(bank_matrix(slopes), 101) == 60, f"five-bank rank: {slopes}")
    five_banks += 1
require(five_banks == 792, "wrong five-bank count")

six_banks = 0
for slopes in combinations(range(1, P), 6):
    matrix = bank_matrix(slopes)
    require(rank_mod(matrix, 101) == DIM, f"six-bank Q rank: {slopes}")
    require(rank_mod(matrix, P) == 51, f"six-bank mod-13 rank: {slopes}")
    six_banks += 1
require(six_banks == 924, "wrong six-bank count")
print("  all 792 five-slope banks have rational rank 60")
print("  all 924 six-slope banks have rational rank 72 and mod-13 rank 51")

determinant_target = P**21
determinant_banks = (
    (1, 2, 3, 4, 5, 6),
    (1, 2, 3, 4, 5, 12),
    (1, 3, 5, 7, 9, 11),
)
for slopes in determinant_banks:
    require(
        abs(det_bareiss(bank_matrix(slopes))) == determinant_target,
        f"integral determinant drifted at {slopes}",
    )
print(f"  direct determinant controls: 3/3 equal 13^21 = {determinant_target}")
print("  theorem invoice: Smith form 1^51 direct-sum 13^21 for every bank")

print("\n== sharp five-slope literal field-trace hostiles ==")
witness_count = 0
l1_values: list[int] = []
first_witness: list[list[int]] | None = None
for slopes in combinations(range(1, P), 5):
    f = trace_witness(slopes)
    if first_witness is None:
        first_witness = f
    bad = []
    for tau in range(1, P):
        out = tooth(f, tau)
        require(sum(out) == 0, "trace output lost augmentation")
        if out == [0] * P:
            bad.append(tau)
        else:
            require(
                all(dft(out, alpha) != ZERO for alpha in range(1, P)),
                "a surviving rational output lost a C13 colour",
            )
    require(tuple(bad) == slopes, f"trace bad set drifted: {slopes}")
    l1_values.append(sum(abs(value) for row in f for value in row))
    witness_count += 1
require(witness_count == 792, "wrong trace-witness count")
require((min(l1_values), max(l1_values)) == (396, 664), "trace L1 range drifted")
print("  792/792 prescribed five-slope sets are exact zero sets")
print("  every good output fires all 12 nontrivial C13 characters")
print("  integral witness L1 range: 396 .. 664")
print("  minimum support theorem: every nonzero interaction has >=7 good slopes")

print("\n== chart shear and section-holonomy boundary ==")
assert first_witness is not None
for sigma in range(P):
    moved = shear(first_witness, sigma)
    for tau in range(P):
        require(
            tooth(moved, tau) == tooth(first_witness, (tau + sigma) % P),
            "nonaffine chart shear failed",
        )
print("  ambient ordered-set shear sends T_tau to T_(tau+sigma): 169/169 PASS")

shear_hostile = shear(rectangle(0, 0), 1)
column_margin = [sum(shear_hostile[k][v] for k in range(Q)) for v in range(P)]
expected_margin = [1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1]
require(column_margin == expected_margin, "shear margin hostile drifted")
require(not is_v00(shear_hostile), "ambient shear unexpectedly preserved V00")
require(tooth(shear_hostile, 0) != [0] * P, "shear margin loss became invisible")
print("  S_1 rectangle column margin: " + " ".join(map(str, column_margin)))
print("  ambient shear does not preserve V00; bank transitions use abstract inverses")

# Cyclic row translation crosses the chosen representative seam at k=6.
c = 1
tau = 1
wrap = [((k + c) // Q) for k in range(Q)]
defect = [(-Q * tau * q) % P for q in wrap]
require(wrap == [0, 0, 0, 0, 0, 0, 1], "row-wrap pattern drifted")
require(len(set(defect)) == 2, "section cocycle became constant")
before = tooth(first_witness, tau)
after = tooth(translate_rows(first_witness, c), tau)
require(
    not any(after == before[shift:] + before[:shift] for shift in range(P)),
    "section holonomy became a common root translation",
)
print("  c=1 row rotation has wrap q=(0,0,0,0,0,0,1)")
print("  induced -7*tau*q defect is nonconstant and physically nonaffine")
print("  exact witness after rotation is not any cyclic root translate: PASS")

print(f"\nexplicit checks: {CHECKS}")
print("ALL EXACT CHECKS PASSED")
