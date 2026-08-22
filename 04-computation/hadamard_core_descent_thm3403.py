#!/usr/bin/env python3
"""Exact controls for THM-3403 (Hadamard core descent).

Standard library only.  The program constructs its small matrices itself and
uses integer arithmetic throughout.  It also pins, but does not import or run,
the inert THM-3394 certificate bundle used by the twelve finite corollaries.
The universal statements of THM-3403 are proved analytically in the theorem;
these are positive, boundary, and hostile finite controls.
"""

from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY_PINS = (
    (
        "01-canon/theorems/THM-3394-twelve-formerly-missing-hadamard-orders-through-2000.md",
        "63e3e841b609f2e9a38bbaa13dbbd665d6c05d7594fdb92351653207c3cadf86",
    ),
    (
        "04-computation/hadamard_twelve_order_bank_thm3394.py",
        "7ae931b3cf268550287bd0621b9b85b8ea167126fadfb90d57b5106d0f82fb2d",
    ),
    (
        "04-computation/hadamard_twelve_order_signword_thm3394.b85",
        "68f7ceebb67005bf1b968171f7e6897cc33bde68adbd63f14bd45edfeb7b3f06",
    ),
    (
        "05-knowledge/results/hadamard_twelve_order_bank_thm3394.out",
        "d8efee90947015a7e6fc28a1685cc3d378357a85e1d4814953b32b17c5cd76a9",
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def scalar_eye(n, scalar):
    return [[scalar if i == j else 0 for j in range(n)] for i in range(n)]


def gram_is(a, scalar):
    return matmul(a, transpose(a)) == scalar_eye(len(a), scalar)


def determinant_bareiss(a):
    """Fraction-free exact determinant."""
    a = [row[:] for row in a]
    n = len(a)
    require(all(len(row) == n for row in a), "determinant input is not square")
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division lost exactness")
                a[i][j] = numerator // previous
        previous = pivot
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[-1][-1]


def rank_mod(a, p):
    a = [[x % p for x in row] for row in a]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][col], -1, p)
        a[rank] = [(inverse * x) % p for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                multiplier = a[i][col]
                a[i] = [
                    (x - multiplier * y) % p for x, y in zip(a[i], a[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def swap_columns(a, i, j):
    if i != j:
        for row in a:
            row[i], row[j] = row[j], row[i]


def smith_diagonal(a):
    """Smith invariant factors via exact Euclidean row and column moves."""
    a = [row[:] for row in a]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    require(all(len(row) == cols for row in a), "ragged Smith input")
    k = 0
    while k < rows and k < cols:
        choices = [
            (abs(a[i][j]), i, j)
            for i in range(k, rows)
            for j in range(k, cols)
            if a[i][j]
        ]
        if not choices:
            break
        _, i0, j0 = min(choices)
        a[k], a[i0] = a[i0], a[k]
        swap_columns(a, k, j0)

        while True:
            if a[k][k] < 0:
                a[k] = [-x for x in a[k]]

            reduced = False
            for i in range(k + 1, rows):
                if a[i][k]:
                    quotient = a[i][k] // a[k][k]
                    a[i] = [x - quotient * y for x, y in zip(a[i], a[k])]
                    if a[i][k]:
                        a[k], a[i] = a[i], a[k]
                    reduced = True
                    break
            if reduced:
                continue

            for j in range(k + 1, cols):
                if a[k][j]:
                    quotient = a[k][j] // a[k][k]
                    for i in range(rows):
                        a[i][j] -= quotient * a[i][k]
                    if a[k][j]:
                        swap_columns(a, k, j)
                    reduced = True
                    break
            if reduced:
                continue

            bad = next(
                (
                    (i, j)
                    for i in range(k + 1, rows)
                    for j in range(k + 1, cols)
                    if a[i][j] % a[k][k]
                ),
                None,
            )
            if bad is None:
                break
            i, _ = bad
            a[k] = [x + y for x, y in zip(a[k], a[i])]

        if a[k][k] < 0:
            a[k] = [-x for x in a[k]]
        k += 1

    diagonal = [abs(a[i][i]) for i in range(min(rows, cols))]
    require(
        all(a[i][j] == 0 for i in range(rows) for j in range(cols) if i != j),
        "Smith reducer did not diagonalize",
    )
    nonzero = [x for x in diagonal if x]
    require(
        all(y % x == 0 for x, y in zip(nonzero, nonzero[1:])),
        "Smith divisibility failed",
    )
    return tuple(nonzero + [0] * (len(diagonal) - len(nonzero)))


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for divisor in range(3, isqrt(n) + 1, 2):
        if n % divisor == 0:
            return False
    return True


def is_squarefree(n):
    for prime in range(2, isqrt(n) + 1):
        if n % (prime * prime) == 0:
            return False
    return True


def paley_type_i(q):
    require(is_prime(q) and q % 4 == 3, "Paley-I control requires a prime q=3 mod 4")

    def character(x):
        x %= q
        if x == 0:
            return 0
        return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

    jacobsthal = [[character(i - j) for j in range(q)] for i in range(q)]
    core = [
        [jacobsthal[i][j] - (1 if i == j else 0) for j in range(q)]
        for i in range(q)
    ]
    return [[1] * (q + 1)] + [[1] + row for row in core]


def kronecker(a, b):
    return [
        [aij * bij for aij in arow for bij in brow]
        for arow in a
        for brow in b
    ]


def sylvester(order):
    require(order >= 1 and order & (order - 1) == 0, "Sylvester order is not a power of two")
    h = [[1]]
    h2 = [[1, 1], [1, -1]]
    while len(h) < order:
        h = kronecker(h, h2)
    return h


def descended_design(h):
    n = len(h)
    require(all(len(row) == n for row in h), "Hadamard control is not square")
    require(all(x == 1 for x in h[0]), "first row is not normalized")
    require(all(row[0] == 1 for row in h), "first column is not normalized")
    return [[(1 - h[i][j]) // 2 for j in range(1, n)] for i in range(1, n)]


def expected_standard_snf(m):
    return (1,) + (2,) * (2 * m - 1) + (2 * m,) * (2 * m - 1) + (4 * m,)


def expected_design_snf(m):
    return (1,) * (2 * m - 1) + (m,) * (2 * m - 1) + (2 * m,)


def negative_support_bits(h, rows=True):
    if not rows:
        h = transpose(h)
    n = len(h)
    return [sum(((1 - row[j]) // 2) << j for j in range(n)) for row in h]


def closed_quadruple_count(h, rows=True):
    bits = negative_support_bits(h, rows=rows)
    all_ones = (1 << len(h)) - 1
    count = 0
    for i, j, k, ell in combinations(range(len(h)), 4):
        if bits[i] ^ bits[j] ^ bits[k] ^ bits[ell] in (0, all_ones):
            count += 1
    return count


def compressed_profile(values):
    pieces = []
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[end] == values[start]:
            end += 1
        exponent = end - start
        pieces.append(str(values[start]) if exponent == 1 else "%d^%d" % (values[start], exponent))
        start = end
    return "(" + ",".join(pieces) + ")"


def check_paley_control(q):
    h = paley_type_i(q)
    n = q + 1
    m = n // 4
    v = n - 1
    b = descended_design(h)
    require(gram_is(h, n), "H_%d Gram failure" % n)

    target = [
        [m * ((1 if i == j else 0) + 1) for j in range(v)] for i in range(v)
    ]
    require(matmul(b, transpose(b)) == target, "B_%d row Gram failure" % v)
    require(matmul(transpose(b), b) == target, "B_%d column Gram failure" % v)
    require(all(sum(row) == 2 * m for row in b), "B_%d row sums" % v)
    require(all(sum(row) == 2 * m for row in transpose(b)), "B_%d column sums" % v)

    det_b = determinant_bareiss(b)
    det_h = determinant_bareiss(h)
    require(abs(det_b) == 2 * m ** (2 * m), "B_%d determinant" % v)
    require(abs(det_h) == n ** (n // 2), "H_%d determinant" % n)
    require(det_h == (-2) ** v * det_b, "border determinant identity")
    polar = [[2 * b[j][i] - 1 for j in range(v)] for i in range(v)]
    require(matmul(b, polar) == scalar_eye(v, 2 * m), "B_%d polar inverse" % v)
    adjugate = [[(det_b // (2 * m)) * x for x in row] for row in polar]
    require(matmul(b, adjugate) == scalar_eye(v, det_b), "B_%d adjugate" % v)
    require(
        {abs(x) for row in adjugate for x in row} == {m ** (2 * m - 1)},
        "B_%d cofactor shell" % v,
    )

    toggle_abs = set()
    for i in range(v):
        for j in range(v):
            changed = [row[:] for row in b]
            changed[i][j] ^= 1
            toggle_abs.add(abs(determinant_bareiss(changed)))
    expected_toggle = (2 * m - 1) * m ** (2 * m - 1)
    require(toggle_abs == {expected_toggle}, "B_%d nonuniform bit-toggle shell" % v)

    require(is_squarefree(m), "chosen Paley control unexpectedly nonsquarefree")
    snf_b = smith_diagonal(b)
    snf_h = smith_diagonal(h)
    require(snf_b == expected_design_snf(m), "B_%d Smith form" % v)
    require(snf_h == expected_standard_snf(m), "H_%d Smith form" % n)

    rank2 = rank_mod(b, 2)
    row_closed = closed_quadruple_count(h, rows=True)
    col_closed = closed_quadruple_count(h, rows=False)
    if m % 2:
        require(rank2 == v - 1, "B_%d odd-m rank" % v)
        if m > 1:
            require((row_closed, col_closed) == (0, 0), "H_%d closed quadruple" % n)
    if m > 1 and m % 2 and is_prime(m):
        require(rank_mod(h, m) == 2 * m, "H_%d p-rank" % n)
        require(
            all(x % m == 0 for row in matmul(h, transpose(h)) for x in row),
            "H_%d code self-orthogonality" % n,
        )
    return (
        n,
        m,
        rank2,
        row_closed,
        col_closed,
        abs(det_b),
        expected_toggle,
        snf_b,
        snf_h,
    )


def verify_dependency_pins():
    for relative, expected in DEPENDENCY_PINS:
        path = REPO_ROOT / relative
        require(path.is_file(), "missing pinned dependency: " + relative)
        actual = sha256(path.read_bytes()).hexdigest()
        require(actual == expected, "dependency hash drift: " + relative)
        print("PIN %s %s" % (relative, actual))


def main():
    print("HADAMARD CORE DESCENT THM-3403: EXACT CONTROLS")
    verify_dependency_pins()

    for q in (3, 7, 11, 19):
        data = check_paley_control(q)
        n, m, rank2, row_closed, col_closed, det_b, toggle, snf_b, snf_h = data
        print(
            "Paley-I N=%d m=%d rank2(B)=%d closed=(%d,%d) detB=%d toggle=%d "
            "SNF(B)=%s SNF(H)=%s"
            % (
                n,
                m,
                rank2,
                row_closed,
                col_closed,
                det_b,
                toggle,
                compressed_profile(snf_b),
                compressed_profile(snf_h),
            )
        )

    h8 = paley_type_i(7)
    b8 = descended_design(h8)
    require(rank_mod(b8, 2) == 3, "order-8 even hostile rank changed")
    require(
        (closed_quadruple_count(h8), closed_quadruple_count(h8, rows=False)) == (14, 14),
        "order-8 closed-quadruple hostile changed",
    )
    h4 = paley_type_i(3)
    require(
        (closed_quadruple_count(h4), closed_quadruple_count(h4, rows=False)) == (1, 1),
        "order-4 sharp exception changed",
    )

    h16 = sylvester(16)
    snf16 = smith_diagonal(h16)
    expected_sylvester16 = (1,) + (2,) * 4 + (4,) * 6 + (8,) * 4 + (16,)
    require(snf16 == expected_sylvester16, "Sylvester-16 Smith hostile changed")
    require(snf16 != expected_standard_snf(4), "nonsquarefree hostile became standard")
    print("HOSTILES N=4 closed=(1,1); N=8 rank2(B)=3 closed=(14,14)")
    print("HOSTILE Sylvester N=16 SNF=(1,2^4,4^6,8^4,16), not standard")

    orders = (668, 716, 892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948, 1964)
    quarters = tuple(n // 4 for n in orders)
    require(all(n % 4 == 0 for n in orders), "THM-3394 order not divisible by four")
    require(all(is_prime(m) for m in quarters), "THM-3394 quarter-order primality failure")
    for n, m in zip(orders, quarters):
        print(
            "THM-3394 N=%d m=%d: D01(%d)=2*%d^%d; Smith/code/full-circuit tiers APPLY"
            % (n, m, n - 1, m, 2 * m)
        )
    print("FINITE CONTROLS ONLY; UNIVERSAL CLAIMS USE THE SELF-CONTAINED PROOF")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
