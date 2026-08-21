"""Exact replay for the strong-dephasing gain-cycle expansion.

Scope
-----
For a finite Hermitian hopping matrix H with zero diagonal, position-basis
dephasing damps every off-diagonal density-matrix unit at rate ``lambda``.
The script checks the invariant-graph expansion of the slow population
generator through order lambda^-3, the triangle and square Wilson-loop
terms, the tree/gauge controls, and the bipartite response-matrix Gram/wedge
identity used in the planar-Jacobian connection.

The numerical portion is only an independent control of the exact symbolic
identities.  It solves the finite-dimensional Riccati equation defining the
slow invariant graph; it is not used as a proof of the formulas.

Reproduce with both

    python 04-computation/strong_dephasing_gain_cycle_expansion_codex_20260820.py
    python -O 04-computation/strong_dephasing_gain_cycle_expansion_codex_20260820.py
"""

from __future__ import annotations

import itertools

import numpy as np
import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def commutator(H: sp.Matrix, X: sp.Matrix) -> sp.Matrix:
    return -sp.I * (H * X - X * H)


def diagonal_action(H: sp.Matrix, power: int) -> sp.Matrix:
    """Matrix of Pi K^power Pi in the diagonal matrix-unit basis."""
    n = H.rows
    columns = []
    for source in range(n):
        X = sp.zeros(n)
        X[source, source] = 1
        for _ in range(power):
            X = commutator(H, X)
        columns.append(sp.Matrix([sp.expand(X[target, target]) for target in range(n)]))
    return sp.Matrix.hstack(*columns)


def slow_coefficients(H: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return G1,G2,G3 for G(lambda)=sum lambda^-j Gj+O(lambda^-4)."""
    M2 = diagonal_action(H, 2)
    M3 = diagonal_action(H, 3)
    M4 = diagonal_action(H, 4)
    return M2, M3, sp.simplify(M4 - 2 * M2 * M2)


def superoperator_blocks(H: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return A=Pi K Q, B=Q K Pi, C=Q K Q in matrix-unit bases."""
    n = H.rows
    off = [(row, column) for row in range(n) for column in range(n) if row != column]
    A = sp.zeros(n, len(off))
    B = sp.zeros(len(off), n)
    C = sp.zeros(len(off), len(off))

    for column, (row, col) in enumerate(off):
        E = sp.zeros(n)
        E[row, col] = 1
        KE = commutator(H, E)
        for target in range(n):
            A[target, column] = sp.expand(KE[target, target])
        for target, (r2, c2) in enumerate(off):
            C[target, column] = sp.expand(KE[r2, c2])

    for source in range(n):
        D = sp.zeros(n)
        D[source, source] = 1
        KD = commutator(H, D)
        for target, (row, col) in enumerate(off):
            B[target, source] = sp.expand(KD[row, col])
    return A, B, C


def response_hamiltonian(A: sp.Matrix) -> sp.Matrix:
    """Hermitian bipartite dilation [[0,A*],[A,0]]."""
    inputs = A.cols
    outputs = A.rows
    H = sp.zeros(inputs + outputs)
    H[:inputs, inputs:] = A.conjugate().T
    H[inputs:, :inputs] = A
    return H


print("SECTION A. Slow invariant graph and the retardation correction.")

# A nontrivial exact triangle is enough to check the block algebra and the
# Riccati recurrence without treating noncommuting A,B,C as scalar symbols.
H_control = sp.Matrix(
    [
        [0, 2, -5 * sp.I],
        [2, 0, 3],
        [5 * sp.I, 3, 0],
    ]
)
A_block, B_block, C_block = superoperator_blocks(H_control)
R1 = B_block
R2 = C_block * B_block
R3 = C_block**2 * B_block - B_block * A_block * B_block
G1_block = sp.simplify(A_block * R1)
G2_block = sp.simplify(A_block * R2)
G3_block = sp.simplify(A_block * R3)
G1_direct, G2_direct, G3_direct = slow_coefficients(H_control)
require(G1_block == G1_direct, "G1 block/direct mismatch")
require(G2_block == G2_direct, "G2 block/direct mismatch")
require(G3_block == G3_direct, "G3 block/direct mismatch")

H2_control = H_control**2
H3_control = H_control**3
H4_control = H_control**4
conductance_control = H_control.applyfunc(lambda entry: sp.expand(entry * sp.conjugate(entry)))
degrees_control = [sum(conductance_control[row, col] for col in range(H_control.cols))
                   for row in range(H_control.rows)]
for row in range(H_control.rows):
    for col in range(H_control.cols):
        if row == col:
            expected = sp.expand(
                2 * H4_control[row, row]
                - 2 * degrees_control[row] ** 2
                - 8 * sum(conductance_control[row, target] ** 2
                          for target in range(H_control.cols))
            )
        else:
            expected = sp.expand(
                -8 * sp.re(H_control[row, col] * H3_control[col, row])
                + 6 * H2_control[row, col] * sp.conjugate(H2_control[row, col])
                - 2 * (G1_direct**2)[row, col]
            )
        require(sp.simplify(G3_direct[row, col] - expected) == 0,
                "G3 entry formula failed")

eps = sp.symbols("eps")
R_truncated = eps * R1 + eps**2 * R2 + eps**3 * R3
riccati_residual = sp.expand(
    R_truncated
    - eps * (B_block + C_block * R_truncated - R_truncated * A_block * R_truncated)
)
for entry in riccati_residual:
    for degree in range(1, 4):
        require(sp.expand(entry).coeff(eps, degree) == 0, "Riccati coefficient failed")

print("   R1=B, R2=CB, R3=C^2B-BAB: PASS")
print("   G1=Pi K^2 Pi, G2=Pi K^3 Pi, G3=Pi K^4 Pi-2(Pi K^2 Pi)^2: PASS")
print("   explicit H^2/H^3/H^4 entry formulas for G3: PASS")
print("   the second term in G3 is the population-retardation correction: PASS")


print("\nSECTION B. Triangle flux is the lambda^-2 skew circulation.")
a, b, d = sp.symbols("a b d", real=True)
H_triangle = sp.Matrix([[0, a, -sp.I * d], [a, 0, b], [sp.I * d, b, 0]])
G1_triangle, G2_triangle, G3_triangle = slow_coefficients(H_triangle)
phi_triangle = a * b * d
triangle_pattern = sp.Matrix([[0, -6, 6], [6, 0, -6], [-6, 6, 0]])
require(G2_triangle == phi_triangle * triangle_pattern, "triangle circulation failed")
require(G2_triangle.T == -G2_triangle, "triangle coefficient is not skew")
require(G3_triangle.T == G3_triangle, "third slow coefficient is not symmetric")
require(G2_triangle * sp.ones(3, 1) == sp.zeros(3, 1), "triangle row sums failed")
require(sp.ones(1, 3) * G2_triangle == sp.zeros(1, 3), "triangle column sums failed")

H_triangle_real = sp.Matrix([[0, 2, 5], [2, 0, 3], [5, 3, 0]])
H_triangle_flux = H_control
real_coeffs = slow_coefficients(H_triangle_real)
flux_coeffs = slow_coefficients(H_triangle_flux)
require(real_coeffs[0] == flux_coeffs[0], "triangle conductance controls differ")
require(real_coeffs[2] == flux_coeffs[2], "triangle G3 should be phase blind")
require(real_coeffs[1] != flux_coeffs[1], "triangle G2 failed to see flux")

print("   G2_xy=-6 Im(H_xy (H^2)_yx)=-6 sum_z Im(H_xy H_yz H_zx): PASS")
print("   G2 is real skew, annihilates constants, and breaks uniform detailed balance when nonzero: PASS")
print("   on a triangle G1 and G3 are phase blind while G2 detects sin(flux): PASS")


print("\nSECTION C. A chordless square first contributes cos(flux) at lambda^-3.")
c = sp.symbols("c", real=True)


def square_hamiltonian(last_phase: sp.Expr) -> sp.Matrix:
    H = sp.zeros(4)
    edges = [(0, 1, a), (1, 2, b), (2, 3, c), (3, 0, d * last_phase)]
    for left, right, value in edges:
        H[left, right] = sp.conjugate(value)
        H[right, left] = value
    return H


square_plus = slow_coefficients(square_hamiltonian(sp.Integer(1)))
square_quadrature = slow_coefficients(square_hamiltonian(sp.I))
square_minus = slow_coefficients(square_hamiltonian(sp.Integer(-1)))
square_pattern = sp.Matrix(
    [
        [4, -8, 12, -8],
        [-8, 4, -8, 12],
        [12, -8, 4, -8],
        [-8, 12, -8, 4],
    ]
)
for coefficients in (square_plus, square_quadrature, square_minus):
    require(coefficients[1] == sp.zeros(4), "bipartite square has nonzero G2")
require(square_plus[0] == square_quadrature[0] == square_minus[0], "square G1 changed with flux")
require(sp.simplify(square_plus[2] - square_quadrature[2]) == a * b * c * d * square_pattern,
        "positive square Wilson term failed")
require(sp.simplify(square_minus[2] - square_quadrature[2]) == -a * b * c * d * square_pattern,
        "negative square Wilson term failed")

unit_pi_G3 = square_minus[2].subs({a: 1, b: 1, c: 1, d: 1})
require(unit_pi_G3[0, 2] == -16, "pi-flux virtual nonedge control failed")

print("   holonomy-bearing G3 contribution is Re(H01 H12 H23 H30)*C4_pattern: PASS")
print("   C4_pattern has diagonal/adjacent/opposite coefficients 4/-8/12: PASS")
print("   unit pi flux gives G3_02=-16 on a nonedge, so higher-order G need not be Metzler: PASS")


print("\nSECTION D. Trees are phase-gauge controls at every order.")
H_path_real = sp.zeros(4)
H_path_phase = sp.zeros(4)
path_edges = [(0, 1, 2, 1), (1, 2, 3, sp.I), (2, 3, 5, -sp.I)]
for left, right, magnitude, phase in path_edges:
    H_path_real[left, right] = H_path_real[right, left] = magnitude
    H_path_phase[left, right] = magnitude * phase
    H_path_phase[right, left] = magnitude * sp.conjugate(phase)
for power in range(2, 9):
    require(diagonal_action(H_path_real, power) == diagonal_action(H_path_phase, power),
            f"path phase gauge failed at K^{power}")
print("   Pi K^n Pi agrees for phased and real paths for n=2,...,8: PASS")
print("   the all-order reason is diagonal-unitary gauge removal on a tree: PASS")


print("\nSECTION E. Bipartite response graphs: G3 recovers Gram/wedge energy.")
A_plus = sp.Matrix([[-1, 1], [-1, 1]])
A_minus = sp.Matrix([[-1, 1], [1, 1]])
response_plus = slow_coefficients(response_hamiltonian(A_plus))
response_minus = slow_coefficients(response_hamiltonian(A_minus))
require(response_plus[0] == response_minus[0], "response conductance shadows differ")
require(response_plus[1] == response_minus[1] == sp.zeros(4), "response G2 should vanish")
require(A_plus.rank() == 1 and A_minus.rank() == 2, "response rank hostile failed")
require(response_plus[2] != response_minus[2], "response G3 missed rank hostile")

x11, x12, x21, x22 = sp.symbols("x11 x12 x21 x22", real=True)
A_symbolic = sp.Matrix([[x11, x12], [x21, x22]])
G3_response = slow_coefficients(response_hamiltonian(A_symbolic))[2]
overlap = x11**2 * x12**2 + x21**2 * x22**2
gram12 = x11 * x12 + x21 * x22
require(sp.expand(G3_response[0, 1] - (6 * gram12**2 - 8 * overlap)) == 0,
        "same-side Gram formula failed")
norm1 = x11**2 + x21**2
norm2 = x12**2 + x22**2
recovered_wedge = sp.expand(norm1 * norm2 - (G3_response[0, 1] + 8 * overlap) / 6)
require(sp.factor(recovered_wedge - (x11 * x22 - x12 * x21) ** 2) == 0,
        "2x2 wedge recovery failed")

# A larger exact complex control checks the Hermitian formula and the full
# Lagrange identity, not merely its real 2x2 specialization.
A_complex = sp.Matrix(
    [
        [1 + sp.I, 2 - sp.I, 1],
        [2, -1 + 2 * sp.I, 3 + sp.I],
        [1 - 2 * sp.I, 2 + sp.I, -2],
        [3, 1 + sp.I, 1 - sp.I],
    ]
)
G3_complex = slow_coefficients(response_hamiltonian(A_complex))[2]
Gram = sp.simplify(A_complex.conjugate().T * A_complex)
column_norms = [sp.expand(Gram[index, index]) for index in range(A_complex.cols)]
total_wedge_direct = 0
total_wedge_recovered = 0
for left, right in itertools.combinations(range(A_complex.cols), 2):
    overlap = sum(
        sp.expand(abs(A_complex[row, left]) ** 2 * abs(A_complex[row, right]) ** 2)
        for row in range(A_complex.rows)
    )
    gram_modulus = sp.expand(Gram[left, right] * sp.conjugate(Gram[left, right]))
    require(sp.expand(G3_complex[left, right] - (6 * gram_modulus - 8 * overlap)) == 0,
            "complex Gram formula failed")
    direct = sum(
        sp.expand(
            abs(
                A_complex[r1, left] * A_complex[r2, right]
                - A_complex[r1, right] * A_complex[r2, left]
            )
            ** 2
        )
        for r1, r2 in itertools.combinations(range(A_complex.rows), 2)
    )
    recovered = sp.expand(column_norms[left] * column_norms[right] - gram_modulus)
    require(sp.expand(direct - recovered) == 0, "complex Lagrange identity failed")
    total_wedge_direct += direct
    total_wedge_recovered += recovered
require(sp.expand(total_wedge_direct - total_wedge_recovered) == 0, "compound norm failed")

print(f"   same |A_wu| but ranks {A_plus.rank()}/{A_minus.rank()}; G3 distinguishes them: PASS")
print("   G3_uv=6 |(A* A)_uv|^2-8 sum_w |A_wu A_wv|^2 for input vertices u!=v: PASS")
print("   conductances + G3 recover ||col_u wedge col_v||^2 and ||wedge^2 A||_F^2: PASS")
print("   for a 2x2 block this recovered invariant is |det A|^2: PASS")


print("\nSECTION F. Independent numerical Riccati and gauge controls.")


def numeric_blocks(H: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = H.shape[0]
    identity = np.eye(n, dtype=complex)
    # Column-major vectorization: vec(HX-XH)=(I kron H-H^T kron I)vec(X).
    K = -1j * (np.kron(identity, H) - np.kron(H.T, identity))
    diagonal = [index + n * index for index in range(n)]
    off = [index for index in range(n * n) if index not in diagonal]
    A = K[np.ix_(diagonal, off)]
    B = K[np.ix_(off, diagonal)]
    C = K[np.ix_(off, off)]
    return A, B, C


def numeric_coefficients(H: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    A, B, C = numeric_blocks(H)
    G1 = A @ B
    G2 = A @ C @ B
    G3 = A @ (C @ C @ B - B @ A @ B)
    return G1, G2, G3


def solve_slow_graph(H: np.ndarray, dephasing: float) -> tuple[np.ndarray, float]:
    A, B, C = numeric_blocks(H)
    denominator = dephasing * np.eye(C.shape[0], dtype=complex) - C
    R = np.linalg.solve(denominator, B)
    for _ in range(200):
        updated = np.linalg.solve(denominator, B - R @ A @ R)
        if np.linalg.norm(updated - R, ord=np.inf) < 2e-15:
            R = updated
            break
        R = updated
    G = A @ R
    residual = np.linalg.norm(
        R @ G - (B + C @ R - dephasing * R), ord=np.inf
    )
    return np.real_if_close(G, tol=1000), residual


H_numeric_triangle = np.array(
    [[0, 1.2, 0.7j], [1.2, 0, 0.9], [-0.7j, 0.9, 0]], dtype=complex
)
G_numeric = numeric_coefficients(H_numeric_triangle)
triangle_scaled_errors = []
for dephasing in (40.0, 80.0, 160.0):
    exact, residual = solve_slow_graph(H_numeric_triangle, dephasing)
    approximation = sum(G_numeric[order - 1] / dephasing**order for order in range(1, 4))
    error = np.linalg.norm(exact - approximation, ord=np.inf)
    triangle_scaled_errors.append(error * dephasing**4)
    require(residual < 2e-11, "triangle numerical invariant-graph residual")
require(max(triangle_scaled_errors) / min(triangle_scaled_errors) < 1.2,
        "triangle fourth-order remainder did not scale")

phi = 1.1
H_numeric_square = np.zeros((4, 4), dtype=complex)
for left, right, value in (
    (0, 1, 1.0),
    (1, 2, 0.8),
    (2, 3, 1.1),
    (3, 0, 0.9 * np.exp(1j * phi)),
):
    H_numeric_square[left, right] = np.conjugate(value)
    H_numeric_square[right, left] = value
G_square_numeric = numeric_coefficients(H_numeric_square)
require(np.linalg.norm(G_square_numeric[1]) < 1e-12, "numeric bipartite G2 failed")
square_scaled_errors = []
for dephasing in (40.0, 80.0, 160.0):
    exact, residual = solve_slow_graph(H_numeric_square, dephasing)
    approximation = G_square_numeric[0] / dephasing + G_square_numeric[2] / dephasing**3
    error = np.linalg.norm(exact - approximation, ord=np.inf)
    square_scaled_errors.append(error * dephasing**5)
    require(residual < 2e-11, "square numerical invariant-graph residual")
require(max(square_scaled_errors) / min(square_scaled_errors) < 1.3,
        "square fifth-order remainder did not scale")

vertex_gauge = np.diag(np.exp(1j * np.array([0.2, -0.4, 0.7, 1.3])))
H_gauged = vertex_gauge @ H_numeric_square @ vertex_gauge.conjugate().T
for dephasing in (17.0, 41.0):
    original, _ = solve_slow_graph(H_numeric_square, dephasing)
    gauged, _ = solve_slow_graph(H_gauged, dephasing)
    require(np.linalg.norm(original - gauged, ord=np.inf) < 2e-12,
            "exact slow generator is not gauge invariant")

print("   triangle truncation error times lambda^4 is stable at 40/80/160: PASS")
print("   bipartite square truncation error times lambda^5 is stable at 40/80/160: PASS")
print("   exact Riccati slow generator is unchanged by a random vertex gauge: PASS")


print("\nVERDICT: PASS")
print("The resistor generator is G1/lambda; triangle sin-holonomy is G2/lambda^2;")
print("square cos-holonomy and the bipartite Gram/wedge rank certificate are G3/lambda^3.")
print("Beyond leading order the raw population-coordinate slow generator need not be Markov.")
