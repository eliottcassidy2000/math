#!/usr/bin/env python3
"""Exact tangent/second-jet audit for THM-3558.

The script studies every affine plane through a pair of the three collision
points of the fixed THM-1300 Keller map.  It proves that the two image germs
cannot lie in one target hypersurface smooth at the collision value.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    """An active truth gate, including under ``python -O``."""
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


def matrix_equal(left: sp.Matrix, right: sp.Matrix) -> bool:
    return left.shape == right.shape and all(
        sp.simplify(a - b) == 0 for a, b in zip(left, right)
    )


def zero_mod_13(expr: sp.Expr, a: sp.Symbol) -> bool:
    """Test a rational-function identity in Q[a]/(a^2-13)."""
    numerator = sp.cancel(expr).as_numer_denom()[0]
    remainder = sp.rem(sp.Poly(numerator, a), sp.Poly(a**2 - 13, a))
    return remainder.as_expr() == 0


def matrix_equal_mod_13(left: sp.Matrix, right: sp.Matrix, a: sp.Symbol) -> bool:
    return left.shape == right.shape and all(
        zero_mod_13(u - v, a) for u, v in zip(left, right)
    )


x, y, z = sp.symbols("x y z")
u = 1 + x * y
F = sp.Matrix(
    [
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    ]
)
J = F.jacobian((x, y, z))

p0 = sp.Matrix([0, 0, sp.Rational(-1, 4)])
pp = sp.Matrix([1, sp.Rational(-3, 2), sp.Rational(13, 2)])
pm = sp.Matrix([-1, sp.Rational(3, 2), sp.Rational(13, 2)])
q = sp.Matrix([sp.Rational(-1, 4), 0, 0])
points = (p0, pp, pm)


def at(expr: sp.Expr | sp.Matrix, point: sp.Matrix):
    substitution = dict(zip((x, y, z), point))
    return expr.subs(substitution)


require(sp.factor(J.det()) == -2, "ambient Keller determinant")
for index, point in enumerate(points):
    require(matrix_equal(at(F, point), q), f"collision point {index}")

# A row normal n to a source plane is transported to the target branch normal
# m=n J_F(p)^(-1).  Parallel target normals are therefore necessary for two
# image germs to lie in one target hypersurface smooth at q.
B, C = sp.symbols("B C")

# Pair p0,p+.
n_0p = sp.Matrix([[sp.Rational(6, 4) * B - sp.Rational(27, 4) * C, B, C]])
m_0 = sp.simplify(n_0p * at(J, p0).inv())
m_p = sp.simplify(n_0p * at(J, pp).inv())
cross_0p = sp.Matrix(m_0).cross(sp.Matrix(m_p)).applyfunc(sp.factor)
expected_0p = sp.Matrix(
    [[
        -sp.Rational(27, 32) * (B**2 - 6 * B * C + 3 * C**2),
        sp.Rational(27, 16) * (4 * B**2 - 38 * B * C + 77 * C**2),
        -sp.Rational(3, 4) * (8 * B**2 - 53 * B * C + C**2),
    ]]
)
require(matrix_equal(cross_0p, expected_0p), "p0,p+ normal cross product")
res_0p = sp.factor(
    sp.resultant(
        B**2 - 6 * B * C + 3 * C**2,
        4 * B**2 - 38 * B * C + 77 * C**2,
        B,
    )
)
require(res_0p == -647 * C**4, "p0,p+ tangent resultant")

# Pair p0,p-.  The sign-reversed formulas are an independent symmetry check.
n_0m = sp.Matrix([[sp.Rational(6, 4) * B + sp.Rational(27, 4) * C, B, C]])
m_0m = sp.simplify(n_0m * at(J, p0).inv())
m_m = sp.simplify(n_0m * at(J, pm).inv())
cross_0m = sp.Matrix(m_0m).cross(sp.Matrix(m_m)).applyfunc(sp.factor)
expected_0m = sp.Matrix(
    [[
        -sp.Rational(27, 32) * (B**2 + 6 * B * C + 3 * C**2),
        -sp.Rational(27, 16) * (4 * B**2 + 38 * B * C + 77 * C**2),
        sp.Rational(3, 4) * (8 * B**2 + 53 * B * C + C**2),
    ]]
)
require(matrix_equal(cross_0m, expected_0m), "p0,p- normal cross product")
res_0m = sp.factor(
    sp.resultant(
        B**2 + 6 * B * C + 3 * C**2,
        4 * B**2 + 38 * B * C + 77 * C**2,
        B,
    )
)
require(res_0m == -647 * C**4, "p0,p- tangent resultant")

# If C=0, the first component in either cross product is a nonzero multiple
# of B^2.  The resultants handle C!=0.  Hence neither pair has a nonzero
# projective normal with matching target tangent planes.

# Pair p+,p-.  This time a quadratic exceptional locus survives first order.
r, c = sp.symbols("r c")
n_pm = sp.Matrix([[3 * r, 2 * r, 2 * c]])
m_plus = sp.simplify(n_pm * at(J, pp).inv())
m_minus = sp.simplify(n_pm * at(J, pm).inv())
cross_pm = sp.Matrix(m_plus).cross(sp.Matrix(m_minus)).applyfunc(sp.factor)
expected_pm = sp.Matrix([[0, 27 * (r**2 - 13 * c**2), 12 * (r**2 - 13 * c**2)]])
require(matrix_equal(cross_pm, expected_pm), "p+,p- normal cross product")

# On the exceptional planes normalize c=1 and write a=r, a^2=13:
#   3 a x + 2 a y + 2 z - 13 = 0.
# Compute the target graph F3=g(F1,F2) to second order on each source branch.
a = sp.symbols("a")
z_plane = (13 - a * (3 * x + 2 * y)) / 2
H = [sp.expand(component.subs(z, z_plane)) for component in F]


def target_graph_jet(source_point: tuple[sp.Expr, sp.Expr]):
    substitution = {x: source_point[0], y: source_point[1]}
    A = sp.Matrix(
        [[sp.diff(H[i], variable).subs(substitution) for variable in (x, y)] for i in (0, 1)]
    )
    w_row = sp.Matrix(
        [[sp.diff(H[2], variable).subs(substitution) for variable in (x, y)]]
    )
    gradient = sp.simplify(w_row * A.inv())
    source_hessians = [sp.hessian(component, (x, y)).subs(substitution) for component in H]
    residual = source_hessians[2] - gradient[0] * source_hessians[0] - gradient[1] * source_hessians[1]
    target_hessian = sp.simplify(A.inv().T * residual * A.inv())
    return A, gradient, target_hessian


A_plus, grad_plus, hess_plus = target_graph_jet((1, sp.Rational(-3, 2)))
A_minus, grad_minus, hess_minus = target_graph_jet((-1, sp.Rational(3, 2)))

expected_A_plus = sp.Matrix(
    [
        [3 * (a - 3) / 16, (a - 3) / 8],
        [-3 * (3 * a - 1) / 8, -(3 * a - 25) / 4],
    ]
)
expected_A_minus = sp.Matrix(
    [
        [3 * (a + 3) / 16, (a + 3) / 8],
        [3 * (3 * a + 1) / 8, (3 * a + 25) / 4],
    ]
)
require(matrix_equal(A_plus, expected_A_plus), "plus branch graph-coordinate Jacobian")
require(matrix_equal(A_minus, expected_A_minus), "minus branch graph-coordinate Jacobian")
require(zero_mod_13(A_plus.det() - 9 * (a - 3) / 8, a), "plus graph determinant")
require(zero_mod_13(A_minus.det() - 9 * (a + 3) / 8, a), "minus graph determinant")

expected_gradient = sp.Matrix([[-32 * a / 9, sp.Rational(4, 9)]])
require(matrix_equal_mod_13(grad_plus, expected_gradient, a), "plus graph gradient")
require(matrix_equal_mod_13(grad_minus, expected_gradient, a), "minus graph gradient")

expected_hess_plus = sp.Matrix(
    [
        [320 * a / 9 + sp.Rational(10816, 27), 16 * a / 9 + sp.Rational(224, 9)],
        [16 * a / 9 + sp.Rational(224, 9), -sp.Rational(16, 27)],
    ]
)
expected_hess_minus = sp.Matrix(
    [
        [320 * a / 9 - sp.Rational(10816, 27), sp.Rational(224, 9) - 16 * a / 9],
        [sp.Rational(224, 9) - 16 * a / 9, sp.Rational(16, 27)],
    ]
)
require(matrix_equal_mod_13(hess_plus, expected_hess_plus, a), "plus graph Hessian")
require(matrix_equal_mod_13(hess_minus, expected_hess_minus, a), "minus graph Hessian")

hessian_difference = expected_hess_plus - expected_hess_minus
expected_difference = sp.Matrix(
    [
        [sp.Rational(21632, 27), 32 * a / 9],
        [32 * a / 9, -sp.Rational(32, 27)],
    ]
)
require(matrix_equal(hessian_difference, expected_difference), "branch Hessian mismatch")
require(hessian_difference[1, 1] != 0, "nonzero second-jet obstruction")

print("THM-3558 exact tangent/second-jet audit")
print("ambient determinant: -2")
print("collision fibre: p0,p+,p- -> (-1/4,0,0)")
print("p0,p+ tangent resultant:", res_0p)
print("p0,p- tangent resultant:", res_0m)
print("p+,p- tangent cross factor: r^2-13*c^2")
print("exceptional plane graph determinants: 9*(a-3)/8, 9*(a+3)/8")
print("common target gradient: (-32*a/9,4/9), a^2=13")
print("Hessian difference:")
print(expected_difference)
print("all active truth gates passed")
