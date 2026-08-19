#!/usr/bin/env python3
"""Exact low-degree invariant-graph search for the fixed sporadic Keller map.

For h in Qbar[x,y] of total degree at most D, form

    E_h(x,y) = F3(x,y,h(x,y)) - h(F1(x,y,h), F2(x,y,h)).

The graph z=h(x,y) is F-invariant iff E_h is the zero polynomial.  For each
D=0,1,2,3 this script builds the exact coefficient ideal over QQ and asks for
a Groebner unit certificate.  It also checks the normal-multiplier identity
which any invariant graph would inherit from det(JF)=-2.

This is a finite coefficient-space experiment.  It does not classify graphs
of degree above three, nongraph hypersurfaces, or polynomial semiconjugacies.
"""

from __future__ import annotations

import hashlib
from collections import defaultdict

import sympy as sp


x, y, z = sp.symbols("x y z")
u = 1 + x * y
F1 = u**3 * z + y**2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_expression(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Return a deterministic primitive QQ polynomial representative."""

    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    _, primitive = poly.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive.as_expr()


def coefficient_system(D: int):
    indices = [(i, j) for total in range(D + 1) for i in range(total + 1) for j in [total - i]]
    coeffs = sp.symbols(" ".join(f"c_{i}_{j}" for i, j in indices), seq=True)
    h = sum(c * x**i * y**j for c, (i, j) in zip(coeffs, indices))
    f1h = sp.expand(F1.subs(z, h))
    f2h = sp.expand(F2.subs(z, h))
    f3h = sp.expand(F3.subs(z, h))
    evaluated_h = sp.expand(h.subs({x: f1h, y: f2h}, simultaneous=True))
    defect = sp.Poly(sp.expand(f3h - evaluated_h), x, y)

    by_degree: dict[int, list[sp.Expr]] = defaultdict(list)
    seen: set[str] = set()
    for (i, j), value in defect.terms():
        eq = primitive_expression(value, coeffs)
        if eq == 0:
            continue
        key = sp.srepr(eq)
        if key not in seen:
            seen.add(key)
            by_degree[i + j].append(eq)

    equations = [eq for degree in sorted(by_degree) for eq in by_degree[degree]]
    return indices, coeffs, h, defect, by_degree, equations


def unit_basis(equations, coeffs):
    basis = sp.groebner(equations, *coeffs, order="grevlex", domain=sp.QQ, method="f5b")
    is_unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    return basis, is_unit


def modular_unit(equations, coeffs, prime: int) -> bool:
    basis = sp.groebner(
        equations, *coeffs, order="grevlex", modulus=prime, method="f5b"
    )
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def modular_selector(equations, coeffs, prime: int):
    """Use a finite field only to select a small candidate QQ certificate."""

    core = list(equations)
    require(modular_unit(core, coeffs, prime), "full modular system is not a unit ideal")
    partitions = 2
    while len(core) >= 2:
        width = (len(core) + partitions - 1) // partitions
        removed = False
        for start in range(0, len(core), width):
            trial = core[:start] + core[start + width :]
            if trial and modular_unit(trial, coeffs, prime):
                core = trial
                partitions = max(2, partitions - 1)
                removed = True
                break
        if not removed:
            if partitions >= len(core):
                break
            partitions = min(len(core), 2 * partitions)

    # Make the returned modular core deletion-minimal.  This is selection,
    # never the characteristic-zero proof; the caller recomputes over QQ.
    changed = True
    while changed:
        changed = False
        for index in range(len(core) - 1, -1, -1):
            trial = core[:index] + core[index + 1 :]
            if trial and modular_unit(trial, coeffs, prime):
                core = trial
                changed = True
    return core


def greedy_unit_core(equations, coeffs):
    """Find a short prefix certificate, then delete redundant equations."""

    # Unit-ideal generation is monotone under adjoining equations.  A binary
    # search for the first unit prefix avoids dozens of expensive deletions
    # from the full high-degree coefficient list.
    low, high = 1, len(equations)
    while low < high:
        middle = (low + high) // 2
        _, is_unit = unit_basis(equations[:middle], coeffs)
        if is_unit:
            high = middle
        else:
            low = middle + 1
    core = list(equations[:low])
    changed = True
    while changed:
        changed = False
        for index in range(len(core) - 1, -1, -1):
            trial = core[:index] + core[index + 1 :]
            if not trial:
                continue
            _, is_unit = unit_basis(trial, coeffs)
            if is_unit:
                core = trial
                changed = True
    basis, is_unit = unit_basis(core, coeffs)
    require(is_unit, "greedy core lost the unit certificate")
    return core, basis


def normal_multiplier(h: sp.Expr, f1h: sp.Expr, f2h: sp.Expr) -> sp.Expr:
    hx = sp.diff(h, x).subs({x: f1h, y: f2h}, simultaneous=True)
    hy = sp.diff(h, y).subs({x: f1h, y: f2h}, simultaneous=True)
    return sp.expand(-x**3 - u**3 * hx - 3 * x * u**2 * hy)


def constant_nonzero_multiplier_system(q0: sp.Expr, coeffs):
    """Equations saying q0 is a nonzero constant, with one inverse variable."""

    inverse = sp.Symbol("q_inverse")
    variables = tuple(coeffs) + (inverse,)
    qpoly = sp.Poly(q0, x, y)
    equations = []
    for monomial, value in qpoly.terms():
        if monomial != (0, 0):
            equations.append(primitive_expression(value, variables))
    constant = qpoly.coeff_monomial(1)
    equations.append(primitive_expression(inverse * constant - 1, variables))
    # Preserve order while removing duplicate primitive equations.
    unique = []
    seen = set()
    for equation in equations:
        key = sp.srepr(equation)
        if key not in seen:
            seen.add(key)
            unique.append(equation)
    return variables, unique


def main() -> None:
    jac = sp.expand(sp.det(sp.Matrix([F1, F2, F3]).jacobian([x, y, z])))
    require(jac == -2, "fixed map failed the Keller control")
    print("fixed map control: det JF = -2")
    Xx, Xy, Xp, Yx, Yy, Yp, Ex, Ey, Q = sp.symbols(
        "Xx Xy Xp Yx Yy Yp Ex Ey Q"
    )
    generic_det = sp.det(sp.Matrix([[Xx, Xy, Xp], [Yx, Yy, Yp], [Ex, Ey, Q]]))
    block_formula = (
        Q * (Xx * Yy - Xy * Yx)
        - Yp * (Xx * Ey - Xy * Ex)
        + Xp * (Yx * Ey - Yy * Ex)
    )
    require(sp.expand(generic_det - block_formula) == 0,
            "generic adapted-coordinate determinant formula failed")
    print("generic graph-coordinate determinant formula verified exactly")
    c_orbit = 2 - 3 * x * y - x**2 * z
    require(sp.expand(F3 - x * c_orbit) == 0, "orbit-cylinder factorization failed")
    plane_f1 = sp.expand(F1.subs(x, 0))
    plane_f2 = sp.expand(F2.subs(x, 0))
    plane_jac = sp.det(sp.Matrix([plane_f1, plane_f2]).jacobian([y, z]))
    require((plane_f1, plane_f2, sp.expand(F3.subs(x, 0))) == (z + 4 * y**2, y, 0),
            "x=0 plane restriction failed")
    require(plane_jac == -1, "x=0 plane restriction is not Keller")
    require(sp.expand(x * (3 * y + x * z) - 2 + c_orbit) == 0,
            "orbit-cylinder unit identity failed")
    gradient = [sp.diff(c_orbit, variable) for variable in (x, y, z)]
    require(all(value.subs({x: 0, y: 0}) == 0 for value in gradient),
            "orbit-cylinder gradient hostile failed")
    collision_plus = {x: 1, y: sp.Rational(-3, 2), z: sp.Rational(13, 2)}
    collision_minus = {x: -1, y: sp.Rational(3, 2), z: sp.Rational(13, 2)}
    require(c_orbit.subs(collision_plus) == c_orbit.subs(collision_minus) == 0,
            "orbit collision does not lie on C=0")
    require(
        tuple(component.subs(collision_plus) for component in (F1, F2, F3))
        == tuple(component.subs(collision_minus) for component in (F1, F2, F3))
        == (sp.Rational(-1, 4), 0, 0),
        "orbit-cylinder collision control failed",
    )
    print("near-miss control: F3=x*C, C=2-3xy-x^2z")
    print("  x=0 -> z=0 gives (y,z)->(z+4y^2,y), planar Jacobian -1")
    print("  C=0 contains the colliding orbit pair and maps to z=0")
    print("  on C=0, x*(3y+xz)=2, so x is a nonconstant unit: C=0 is not A2")
    print("  grad(C) vanishes on x=y=0, so C is not a coordinate polynomial")
    print("universe: h in Qbar[x,y], total degree <= D; graph invariance E_h == 0")
    print("coefficient ideals and Groebner calculations are over QQ")

    digest_parts: list[str] = []
    for D in range(4):
        indices, coeffs, h, defect, by_degree, equations = coefficient_system(D)
        max_xy_degree = max((sum(monomial) for monomial, _ in defect.terms()), default=-1)
        print()
        print(f"D={D}: unknowns={len(coeffs)}, nonzero coefficient equations={len(equations)}, "
              f"max_xy_degree={max_xy_degree}")
        print("  equations by x,y total degree: " +
              ", ".join(f"{degree}:{len(by_degree[degree])}" for degree in sorted(by_degree)))

        f1h = sp.expand(F1.subs(z, h))
        f2h = sp.expand(F2.subs(z, h))
        q0 = normal_multiplier(h, f1h, f2h)

        if D < 3:
            certificate_equations = equations
            selection_note = "full coefficient ideal"
            certificate_variables = coeffs
        else:
            certificate_variables, certificate_equations = constant_nonzero_multiplier_system(q0, coeffs)
            selection_note = (
                "normal-multiplier necessary ideal "
                f"({len(certificate_equations)} equations, including q0 inverse)"
            )
        basis, is_unit = unit_basis(certificate_equations, certificate_variables)
        require(is_unit, f"D={D} coefficient ideal did not reduce to the unit ideal")
        print(f"  {selection_note}")
        print("  exact characteristic-zero Groebner basis: [1]")

        if D <= 1:
            core, core_basis = greedy_unit_core(equations, coeffs)
            print(f"  greedy exact unit core size={len(core)}")
            for number, equation in enumerate(core, 1):
                print(f"    g{number} = {sp.sstr(equation)}")
            require(len(core_basis.polys) == 1 and core_basis.polys[0].as_expr() == 1,
                    f"D={D} core basis is not [1]")
            digest_parts.extend([f"D={D}", *(sp.srepr(eq) for eq in core)])
        else:
            # The full characteristic-zero unit basis is the certificate.  A
            # deletion-minimal core is deliberately not sought: repeated
            # Groebner runs are much more expensive and add no logical force.
            equation_digest = hashlib.sha256(
                "\n".join(sp.srepr(eq) for eq in equations).encode("utf-8")
            ).hexdigest()
            print(f"  coefficient-system sha256={equation_digest}")
            digest_parts.extend([f"D={D}", equation_digest])

        # If the coefficient equations make E identically zero, then Ex=Ey=0
        # in the generic determinant formula checked above.  Consequently
        # -2=Jac(F1|_Gamma,F2|_Gamma)*q0, forcing both factors to be units.
        print("  graph consequence: invariance forces q0 and planar Jacobian to be nonzero constants")

    digest = hashlib.sha256("\n".join(digest_parts).encode("utf-8")).hexdigest()
    print()
    print("result: no invariant polynomial graph z=h(x,y) exists in degrees D=0,1,2,3")
    print("scope: finite degree <=3 only; nongraph and semiconjugate surfaces remain open")
    print(f"semantic_core_sha256={digest}")


if __name__ == "__main__":
    main()
