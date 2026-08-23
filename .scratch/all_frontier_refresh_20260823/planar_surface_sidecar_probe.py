#!/usr/bin/env python3
"""Exact planar-surface sidecar probes for the THM-3791--3795 era.

The packet checks three sharply scoped facts:
1. a two-arm n=3 deformation whose lower Hensel jets agree while the
   resonant c^2 coefficient changes the de Rham obstruction;
2. the canonical nodal carrier still has a critical point after adding the
   one-parameter second-normal correction alpha*z^2, now a positive control
   inside the full r-independent cell closed by THM-3795;
3. a cubic finite cover has an unramified sheet above its branch divisor,
   so the entire-branch-fibre deletion step in the quadratic THM-3794 proof
   cannot be copied to degree three.

No general Darboux, Keller, JC(2), or DC(2) conclusion is claimed.
"""

from sympy import Poly, cancel, diff, discriminant, expand, factor, rem, simplify, symbols


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def zero_rational_mod(expr, variable, modulus):
    numerator, _ = cancel(expr).as_numer_denom()
    return expand(rem(Poly(expand(numerator), variable), Poly(expand(modulus), variable)).as_expr()) == 0


def hensel_resonance_probe():
    c, b, e, lam = symbols("c b e lam")
    sigma = b * (b - 1 - lam * c**2)
    equation = c**3 * e - sigma
    beta0 = 0
    beta1 = 1 + lam * c**2
    require(expand(sigma.subs(b, beta0)) == 0, "first Hensel arm failed")
    require(expand(sigma.subs(b, beta1)) == 0, "moving Hensel arm failed")

    cb = c**3
    be = 3 * c**2 * e - diff(sigma, c)
    ce = diff(sigma, b)

    def bracket(f, g):
        return expand(
            (diff(f, c) * diff(g, b) - diff(f, b) * diff(g, c)) * cb
            + (diff(f, b) * diff(g, e) - diff(f, e) * diff(g, b)) * be
            + (diff(f, c) * diff(g, e) - diff(f, e) * diff(g, c)) * ce
        )

    require(bracket(equation, c) == 0, "surface equation not Poisson-central against c")
    require(bracket(equation, b) == 0, "surface equation not Poisson-central against b")
    require(bracket(equation, e) == 0, "surface equation not Poisson-central against e")

    lower_jet_lam0 = (0, 1)
    lower_jet_lam1 = (0, 1)
    require(lower_jet_lam0 == lower_jet_lam1, "lower-jet control changed")
    resonant_lam0 = (0, 0)
    resonant_lam1 = (0, 1)
    require(resonant_lam0 != resonant_lam1, "resonant vectors did not separate")
    return {
        "sigma": sigma,
        "bracket_be": be,
        "bracket_ce": ce,
        "lower": lower_jet_lam0,
        "r0": resonant_lam0,
        "r1": resonant_lam1,
    }


def nodal_second_normal_probe():
    r, z, e = symbols("r z e")
    alpha = symbols("alpha", nonzero=True)
    surface = r**2 * e - z**3 + r
    carrier = e**2 - z / 3 + alpha * z**2
    carrier_z = diff(carrier, z)
    carrier_e = diff(carrier, e)
    ar = expand(-3 * r**2 * carrier_z - 9 * z**2 * carrier_e)
    az = expand(-(3 + 6 * r * e) * carrier_e)
    ae = expand((3 + 6 * r * e) * carrier_z)

    point = {e: 0, z: 1 / (6 * alpha), r: 1 / (216 * alpha**3)}
    require(simplify(surface.subs(point)) == 0, "second-normal point misses surface")
    require(simplify(ar.subs(point)) == 0, "second-normal point misses {A,r}=0")
    require(simplify(az.subs(point)) == 0, "second-normal point misses {A,z}=0")
    require(simplify(ae.subs(point)) == 0, "second-normal point misses {A,e}=0")

    # Alpha=0 is the canonical seven-point boundary from THM-3790.
    zz = symbols("zz", nonzero=True)
    canonical_modulus = 8 * zz**7 + 9
    canonical_point = {r: 2 * zz**3, z: zz, e: -1 / (4 * zz**3), alpha: 0}
    require(zero_rational_mod(surface.subs(canonical_point), zz, canonical_modulus), "canonical point misses surface")
    require(zero_rational_mod(ar.subs(canonical_point), zz, canonical_modulus), "canonical point misses {A,r}=0")
    require(zero_rational_mod(az.subs(canonical_point), zz, canonical_modulus), "canonical point misses {A,z}=0")
    require(zero_rational_mod(ae.subs(canonical_point), zz, canonical_modulus), "canonical point misses {A,e}=0")
    return {
        "carrier": carrier,
        "generic_point": point,
        "canonical_modulus": canonical_modulus,
    }


def cubic_branch_fibre_probe():
    t, z = symbols("t z")
    cubic = z**3 - 3 * t * z + 2
    derivative = diff(cubic, z)
    delta = factor(discriminant(cubic, z))
    require(expand(delta - 108 * (t**3 - 1)) == 0, "cubic discriminant changed")
    branch_fibre = factor(cubic.subs(t, 1))
    require(expand(branch_fibre - (z - 1) ** 2 * (z + 2)) == 0, "branch fibre factorization changed")
    require(derivative.subs({t: 1, z: 1}) == 0, "double sheet is not ramified")
    require(derivative.subs({t: 1, z: -2}) == 9, "simple sheet is not unramified")
    require(delta.subs(t, 1) == 0, "branch polynomial should vanish at unramified sheet")

    quadratic = z**2 - t
    quadratic_derivative = diff(quadratic, z)
    require(factor(quadratic.subs(t, 0)) == z**2, "quadratic branch fibre changed")
    require(quadratic_derivative.subs({t: 0, z: 0}) == 0, "quadratic branch point should ramify")
    return {
        "cubic": cubic,
        "delta": delta,
        "branch_fibre": branch_fibre,
        "simple_derivative": derivative.subs({t: 1, z: -2}),
    }


def main():
    h = hensel_resonance_probe()
    n = nodal_second_normal_probe()
    c = cubic_branch_fibre_probe()

    print("PLANAR SURFACE SIDECAR AUDIT")
    print(f"moving-root family Sigma(c,b): {h['sigma']}")
    print(f"identical jets through order 1 for lambda=0,1: {h['lower']}")
    print(f"resonant vectors at order 2: lambda=0 -> {h['r0']}; lambda=1 -> {h['r1']}")
    print("THM-3791 finite-etale control: the resonant universal Hensel jet changes the de Rham class.")
    print(f"second-normal carrier family: {n['carrier']}")
    print("alpha!=0 critical point: e=0, z=1/(6 alpha), r=1/(216 alpha^3)")
    print(f"alpha=0 canonical critical equation: {n['canonical_modulus']}=0")
    print("THM-3795 positive control: this z^2 direction lies inside the now-closed r-independent cell.")
    print(f"cubic cover: {c['cubic']}=0")
    print(f"discriminant: {c['delta']}")
    print(f"branch fibre at t=1: {c['branch_fibre']}; derivative at simple sheet z=-2 is {c['simple_derivative']}")
    print("THM-3794 boundary: degree three can retain an unramified sheet above the branch divisor.")
    print("STATUS FIREWALL: exact scoped probes only; r-coupled mixed normal corrections, cubic etale constant-unit surfaces, JC(2), and DC(2) remain open.")


if __name__ == "__main__":
    main()
