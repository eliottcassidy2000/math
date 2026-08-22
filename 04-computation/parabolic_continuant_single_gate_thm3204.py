#!/usr/bin/env python3
"""Exact controls for THM-3204's parabolic continuant single gate.

No floating point, no randomness, no imported executable.  Every gate is an
explicit ``require`` so that ordinary and ``-O`` replay are byte-identical.
"""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------- Smith form


def smith_invariant_factors(matrix):
    """Invariant factors of a square integer matrix, by elementary moves."""

    a = [list(row) for row in matrix]
    rows, cols = len(a), len(a[0])
    factors = []
    top = 0
    while top < rows and top < cols:
        best, pivot = None, None
        for r in range(top, rows):
            for c in range(top, cols):
                if a[r][c] and (best is None or abs(a[r][c]) < best):
                    best, pivot = abs(a[r][c]), (r, c)
        if pivot is None:
            break
        while True:
            r0, c0 = pivot
            a[top], a[r0] = a[r0], a[top]
            for row in a:
                row[top], row[c0] = row[c0], row[top]
            moved = False
            for r in range(top + 1, rows):
                if a[r][top]:
                    q = a[r][top] // a[top][top]
                    for c in range(top, cols):
                        a[r][c] -= q * a[top][c]
                    if a[r][top]:
                        moved = True
            for c in range(top + 1, cols):
                if a[top][c]:
                    q = a[top][c] // a[top][top]
                    for r in range(top, rows):
                        a[r][c] -= q * a[r][top]
                    if a[top][c]:
                        moved = True
            if not moved:
                break
            best, pivot = None, None
            for r in range(top, rows):
                for c in range(top, cols):
                    if a[r][c] and (best is None or abs(a[r][c]) < best):
                        best, pivot = abs(a[r][c]), (r, c)
        factors.append(abs(a[top][top]))
        top += 1
    for i in range(1, len(factors)):
        g = sp.igcd(factors[i - 1], factors[i])
        product = factors[i - 1] * factors[i]
        factors[i - 1], factors[i] = g, product // g
    return factors


def jacobi(sub, diag, sup):
    """Tridiagonal matrix with the given sub/diagonal/super entries."""

    n = len(diag)
    require(len(sub) == n - 1 and len(sup) == n - 1, "jacobi band lengths")
    return [[diag[i] if i == j else
             sub[j] if i == j + 1 else
             sup[i] if j == i + 1 else 0
             for j in range(n)] for i in range(n)]


def deleted_minor_det(matrix):
    """Determinant of the submatrix on rows 2..n and columns 1..n-1."""

    n = len(matrix)
    block = [[matrix[i + 1][j] for j in range(n - 1)] for i in range(n - 1)]
    return int(sp.Matrix(block).det()) if n > 1 else 1


# --------------------------------------- 1.  Unital Jacobi => cyclic cokernel

# A deterministic bank of unital Jacobi matrices: every sub-diagonal entry is
# +-1, the diagonal and super-diagonal range over a fixed integer window.
CYCLIC_CASES = 0
for size in range(1, 8):
    for shift in range(-3, 4):
        for spread in range(-2, 3):
            sub = [1 if (i + shift) % 2 == 0 else -1 for i in range(size - 1)]
            diag = [shift + spread * i for i in range(size)]
            sup = [1 + spread + i for i in range(size - 1)]
            matrix = jacobi(sub, diag, sup)
            determinant = int(sp.Matrix(matrix).det())
            factors = smith_invariant_factors(matrix)
            require(len(factors) == size or determinant == 0,
                    "full rank invariant factor count")
            if determinant == 0:
                continue
            require(factors[:-1] == [1] * (size - 1),
                    "unital Jacobi cokernel is cyclic")
            require(factors[-1] == abs(determinant),
                    "cyclic order equals |det|")
            require(abs(deleted_minor_det(matrix)) == 1,
                    "deleted minor is the sub-diagonal product")
            CYCLIC_CASES += 1

# Hostile control: break unitality of one sub-diagonal entry and the cokernel
# stops being cyclic.
HOSTILE_NONCYCLIC = jacobi([3], [3, 0], [3])
HOSTILE_FACTORS = smith_invariant_factors(HOSTILE_NONCYCLIC)
require(HOSTILE_FACTORS == [3, 3],
        "non-unital hostile is not cyclic")

# THM-3182/3183's rank-three reset really is obstructed.
SMITH_OBSTRUCTION = []


def scalar_frame(index, d_value, v_value):
    a = 2 * (index + 1) * (2 * index + 1) * v_value
    b = index * (index + 1) * (1 - 4 * d_value * v_value)
    c = d_value - index - 1
    return [[a, b, c], [1, 0, 0], [0, 0, d_value]]


def weighted_frame(index, d_value, v_value):
    return [[-(index + 1), 2 * v_value * (index + 1), d_value],
            [-(index + 1) * (2 * index + 3 + 2 * d_value),
             (index + 1) * (1 + 2 * v_value * (2 * index + 3)),
             d_value * (2 * index + 3)],
            [0, 0, d_value]]


def p_valuation(value, prime):
    require(value != 0, "valuation of zero")
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


for prime in (5, 7, 11, 13):
    for s_value in (1, 2, 3):
        for v_value in (1, 3, 5):
            d_value = prime + s_value
            if (1 - 4 * d_value * v_value) % prime == 0:
                continue
            if (s_value * v_value) % prime == 0:
                continue
            scalar = [p_valuation(x, prime) for x in
                      smith_invariant_factors(scalar_frame(prime - 1,
                                                           d_value, v_value))]
            weighted = [p_valuation(x, prime) for x in
                        smith_invariant_factors(weighted_frame(prime - 1,
                                                               d_value,
                                                               v_value))]
            require(scalar == [0, 0, 1], "scalar reset Smith type (1,1,p)")
            require(weighted == [0, 1, 1], "weighted reset Smith type (1,p,p)")
            SMITH_OBSTRUCTION.append((prime, s_value, v_value))
require(SMITH_OBSTRUCTION, "non-vacuous Smith obstruction bank")


# ------------------------------------ 2.  Chebyshev form and parabolic locus

t_sym, alpha_sym, c_sym, n_sym = sp.symbols("t alpha c n")

CHEBYSHEV_CHECKS = 0
for size in range(1, 12):
    matrix = jacobi([1] * (size - 1), [alpha_sym] * size, [1] * (size - 1))
    determinant = sp.expand(sp.Matrix(matrix).det())
    chebyshev = sp.expand(sp.chebyshevu(size, alpha_sym / 2))
    require(sp.simplify(determinant - chebyshev) == 0,
            "leading minor equals U_n(alpha/2)")
    at_minus_two = sp.expand(determinant.subs(alpha_sym, -2))
    require(at_minus_two == (-1) ** size * (size + 1),
            "parabolic value (-1)^n (n+1)")
    CHEBYSHEV_CHECKS += 1

# Closed forms of the constant continuant, both branches.
PARABOLIC_CHECKS = 0
for alpha_value in (-6, -4, -2, 2, 4, 6):
    c_value = sp.Rational(alpha_value ** 2, 4)
    previous, current = sp.Integer(0), sp.Integer(1)
    for step in range(1, 14):
        previous, current = current, alpha_value * current - c_value * previous
        require(current == (step + 1) * sp.Rational(alpha_value, 2) ** step,
                "parabolic closed form (n+1)(alpha/2)^n")
        PARABOLIC_CHECKS += 1


def continuant_zero_period(prime, alpha_value, c_value):
    """Least q>0 with D_(q-1)=0 mod prime, or 0 when D never vanishes."""

    previous, current = 0, 1
    seen = {}
    step = 0
    while True:
        state = (previous % prime, current % prime)
        if state in seen:
            return 0
        seen[state] = step
        previous, current = current, (alpha_value * current
                                      - c_value * previous) % prime
        step += 1
        if current % prime == 0:
            return step + 1


PERIOD_ROWS = []
for prime in (3, 5, 7, 11, 13):
    parabolic_periods = set()
    other_periods = set()
    for alpha_value in range(prime):
        for c_value in range(1, prime):
            period = continuant_zero_period(prime, alpha_value, c_value)
            require(period > 0,
                    "a unital constant continuant always vanishes somewhere")
            if (alpha_value ** 2 - 4 * c_value) % prime == 0:
                parabolic_periods.add(period)
            else:
                other_periods.add(period)
                require(period != prime,
                        "non-parabolic period is never p")
                require((prime ** 2 - 1) % period == 0,
                        "non-parabolic period divides p^2-1")
    require(parabolic_periods == {prime},
            "parabolic period is exactly p")
    PERIOD_ROWS.append((prime, sorted(parabolic_periods),
                        sorted(other_periods)))


# ---------------------------- 3.  Factorial exterior transfer parabolic locus

i_sym, v_sym, d_sym = sp.symbols("i v d")
a_i = 2 * (i_sym + 1) * (2 * i_sym + 1) * v_sym
b_i = i_sym * (i_sym + 1) * (1 - 4 * d_sym * v_sym)
STEP_DISCRIMINANT = sp.factor(sp.expand((a_i * d_sym) ** 2
                                        + 4 * d_sym * (b_i * d_sym)))
PARABOLIC_POLY = sp.expand(sp.cancel(STEP_DISCRIMINANT
                                     / (4 * d_sym ** 2 * (i_sym + 1))))
EXPECTED_POLY = sp.expand((i_sym + 1) * (2 * i_sym + 1) ** 2 * v_sym ** 2
                          - 4 * i_sym * d_sym * v_sym + i_sym)
require(sp.expand(PARABOLIC_POLY - EXPECTED_POLY) == 0,
        "reduced parabolic condition")

# Two consecutive parabolic steps force v = 0.
CONSECUTIVE_OBSTRUCTION = sp.expand(
    (i_sym + 1) ** 2 * (2 * i_sym + 1) ** 2
    - i_sym * (i_sym + 2) * (2 * i_sym + 3) ** 2)
require(sp.expand(CONSECUTIVE_OBSTRUCTION
                  + (8 * i_sym ** 3 + 20 * i_sym ** 2 + 12 * i_sym - 1)) == 0,
        "consecutive obstruction polynomial")

ELIMINATION = sp.expand(
    (i_sym + 1) * (PARABOLIC_POLY)
    - i_sym * (PARABOLIC_POLY.subs(i_sym, i_sym + 1)))
require(sp.expand(ELIMINATION - CONSECUTIVE_OBSTRUCTION * v_sym ** 2) == 0,
        "elimination of d between consecutive parabolic conditions")

CONSECUTIVE_ROWS = []
for index in range(1, 13):
    value = int(CONSECUTIVE_OBSTRUCTION.subs(i_sym, index))
    require(value < 0, "consecutive obstruction is nonzero and negative")
    solutions = sp.solve([PARABOLIC_POLY.subs(i_sym, index),
                          PARABOLIC_POLY.subs(i_sym, index + 1)],
                         [v_sym, d_sym], dict=True)
    require(all(solution.get(v_sym, None) == 0 for solution in solutions),
            "every simultaneous solution has v = 0")
    CONSECUTIVE_ROWS.append((index, value, len(solutions)))

# Positive control: a single step really can be made parabolic.
SINGLE_STEP_WITNESSES = []
for index in (1, 2, 3, 4):
    v_value = sp.Integer(1)
    d_value = sp.Rational(((index + 1) * (2 * index + 1) ** 2 * v_value ** 2
                           + index), 4 * index * v_value)
    check = PARABOLIC_POLY.subs({i_sym: index, v_sym: v_value,
                                 d_sym: d_value})
    require(sp.simplify(check) == 0, "single-step parabolic witness")
    neighbour = PARABOLIC_POLY.subs({i_sym: index + 1, v_sym: v_value,
                                     d_sym: d_value})
    require(sp.simplify(neighbour) != 0, "neighbour step stays non-parabolic")
    SINGLE_STEP_WITNESSES.append((index, sp.nsimplify(d_value)))


print("THM-3204 PARABOLIC CONTINUANT SINGLE GATE EXACT CONTROL")
print("unital_jacobi_cyclic_cokernel_cases=" + str(CYCLIC_CASES))
print("deleted_minor_identity=det(rows_2..n,cols_1..n-1)=prod(subdiagonal)")
print("hostile_non_unital_invariant_factors=" + repr(HOSTILE_FACTORS))
print("chebyshev_identity_checks=" + str(CHEBYSHEV_CHECKS))
print("chebyshev_form=det tridiag(1,alpha,1)_n = U_n(alpha/2)")
print("parabolic_locus=alpha=+-2  (Chebyshev endpoints)")
print("parabolic_closed_form_checks=" + str(PARABOLIC_CHECKS))
for prime, parabolic, other in PERIOD_ROWS:
    print("zero_period p=" + str(prime)
          + " parabolic=" + repr(parabolic)
          + " nonparabolic=" + repr(other))
print("smith_obstruction_bank=" + str(len(SMITH_OBSTRUCTION)))
print("scalar_reset_smith=(1,1,p)_admits_jacobi_form")
print("weighted_reset_smith=(1,p,p)_forbids_jacobi_form")
print("factorial_parabolic_condition="
      + str(sp.factor(PARABOLIC_POLY)))
print("consecutive_obstruction=-(8i^3+20i^2+12i-1)")
print("consecutive_rows=" + repr(CONSECUTIVE_ROWS))
print("single_step_parabolic_witnesses="
      + repr([(index, str(value)) for index, value in SINGLE_STEP_WITNESSES]))
print("scope=continuant_gate_structure_not_NC2_GMC2_JC2_or_LRC14")
print("ALL EXACT CHECKS PASSED")
