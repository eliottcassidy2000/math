#!/usr/bin/env python3
"""Exact scratch referee for two-scale Rule 30 projective profiles at n=3.

The algebraic theorem is carried by the displayed twisted-amplitude proof.
This script checks every matrix, parametrization, cycle identity, saturation
factor, and small 2-adic hostile with explicit exceptions (no assert gates).
"""

import hashlib
import json

import sympy as sp


x, h, lam, t = sp.symbols("x h lambda t")


def require(condition, label):
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def recurrence(g):
    n = len(g)
    return tuple(sp.factor(
        -g[(2*j) % n] * g[(2*j+1) % n]
        * (1-g[(2*j+2) % n]) / (1-g[(2*j) % n])
    ) for j in range(n))


def twisted_operator(q):
    # A_(j+3)=q*A_j and B_j=A_(2j)+A_(2j+1), j=0,1,2.
    return sp.Matrix([[1, 1, 0], [q, 0, 1], [0, q, q]])


def ratios(A, q):
    return (
        sp.factor(-A[1] / A[0]),
        sp.factor(-A[2] / A[1]),
        sp.factor(-q * A[0] / A[2]),
    )


def reduce_cyclotomic(expr):
    numerator, denominator = sp.fraction(sp.cancel(expr))
    coefficient_field = sp.QQ.frac_field(lam)
    mod = sp.Poly(h**2 + h + 1, h, domain=coefficient_field)
    num = sp.rem(sp.Poly(numerator, h, domain=coefficient_field), mod).as_expr()
    den = sp.rem(sp.Poly(denominator, h, domain=coefficient_field), mod).as_expr()
    return sp.cancel(num / den)


def matrix_reduce_cyclotomic(M):
    return M.applyfunc(reduce_cyclotomic)


def v2_integer(n):
    n = abs(int(n))
    if n == 0:
        return 99
    out = 0
    while n % 2 == 0:
        n //= 2
        out += 1
    return out


def main():
    # Product/holonomy law for every odd period, checked symbolically at n=3.
    g0, g1, g2 = sp.symbols("g0 g1 g2")
    g = (g0, g1, g2)
    rg = recurrence(g)
    p = sp.prod(g)
    require(sp.factor(sp.prod(rg) + p**2) == 0, "product maps p -> -p^2")
    # A two-cycle has p -> -p^2 -> -p^4 = p, hence p^3=-1.

    T_h = twisted_operator(h)
    T_h2 = twisted_operator(h**2)
    M_h = sp.simplify(T_h2 * T_h)
    expected_M = sp.Matrix([
        [h + 1, 1, 1],
        [h**2, h**2 + h, h],
        [h**3, h**3, h**3 + h**2],
    ])
    require((M_h-expected_M).applyfunc(sp.expand) == sp.zeros(3),
            "twisted two-scale matrix")
    char_general = sp.factor(M_h.charpoly(lam).as_expr())
    expected_char = -(h-lam)*(h**2-lam)*(1+h+h**2+h**3-lam)
    require(sp.expand(char_general-expected_char) == 0, "general characteristic")

    # Ordinary holonomy h=1.
    T1 = twisted_operator(1)
    M1 = T1**2
    require(sp.factor(M1.charpoly(lam).as_expr()) == (lam-4)*(lam-1)**2,
            "h=1 characteristic")
    A = sp.Matrix([t+1, -2, 1-t])
    B = T1*A
    require(B == sp.Matrix([t-1, 2, -t-1]), "period-three first scale")
    require(M1*A == A, "lambda-one plane parametrization")
    gt = ratios(A, 1)
    ht = ratios(B, 1)
    expected_g = (2/(t+1), (1-t)/2, (t+1)/(t-1))
    expected_h = (2/(1-t), (t+1)/2, (t-1)/(t+1))
    require(all(sp.factor(a-b) == 0 for a, b in zip(gt, expected_g)), "g(t)")
    require(all(sp.factor(a-b) == 0 for a, b in zip(ht, expected_h)), "h(t)")
    rgt = recurrence(gt)
    require(all(sp.factor(a-b) == 0 for a, b in zip(rgt, ht)), "R(g(t))=g(-t)")
    rrgt = recurrence(rgt)
    require(all(sp.factor(a-b) == 0 for a, b in zip(rrgt, gt)), "R^2(g(t))=g(t)")
    require(sp.factor(sp.prod(gt) + 1) == 0, "ordinary holonomy product")
    saturation = sp.factor(sp.prod([
        item for value in gt + ht for item in (value, 1-value)
    ]))
    require(saturation != 0, "saturation product nonzero rational function")
    saturation_num = sp.factor(sp.fraction(sp.cancel(saturation))[0])
    require(sp.factor(saturation_num).subs(t, 2) != 0, "t=2 saturated control")
    require(all(sp.factor(a-b) == 0 for a, b in zip(
        tuple(item.subs(t, -t) for item in gt), ht)), "involution t -> -t")
    # t=0 is the only saturated projectively fixed point of this plane;
    # beta=0 is the other T eigenline but has a zero amplitude.
    require(all(sp.factor(a.subs(t, 0)-b.subs(t, 0)) == 0
                for a, b in zip(gt, ht)), "t=0 fixed boundary")

    # Independent direct elimination on the ordinary-holonomy product p=-1.
    # Put g2=-1/(g0*g1); after clearing the saturated denominators, all three
    # two-cycle equations share one curve factor.  Off it they force the
    # isolated constant point.
    aa, bb = sp.symbols("a b")
    direct_g = (aa, bb, -1/(aa*bb))
    direct_r2 = recurrence(recurrence(direct_g))
    direct_differences = tuple(sp.factor(sp.together(direct_r2[j]-direct_g[j]))
                               for j in range(3))
    curve_factor = aa*bb-aa+1
    direct_numerators = tuple(sp.factor(sp.fraction(item)[0])
                              for item in direct_differences)
    require(all(sp.rem(sp.Poly(num, aa, bb), sp.Poly(curve_factor, aa, bb)) == 0
                for num in direct_numerators), "common direct curve factor")
    direct_quotients = tuple(sp.factor(num/curve_factor) for num in direct_numerators)
    require(direct_quotients == (-(aa+1), -(bb+1), -(aa*bb-1)),
            "direct residual factors")

    # Nontrivial cubic holonomy: h^2+h+1=0.  The three eigenvalues are
    # h,h^2,1, and every eigenline has a zero amplitude coordinate.
    M_cyc = matrix_reduce_cyclotomic(M_h)
    char_cyc = reduce_cyclotomic(char_general)
    require(sp.factor(char_cyc) == (lam-1)*(lam**2+lam+1),
            "primitive holonomy characteristic")
    eigen_controls = (
        (h, sp.Matrix([-1, 1, 0])),
        (h**2, sp.Matrix([0, -1, 1])),
        (1, sp.Matrix([1, 0, -h])),
    )
    for eigenvalue, vector in eigen_controls:
        residual = matrix_reduce_cyclotomic(M_cyc*vector-eigenvalue*vector)
        require(residual == sp.zeros(3, 1), ("primitive eigenvector", eigenvalue))
        require(any(entry == 0 for entry in vector), ("zero amplitude", eigenvalue))

    # Period one: genuine algebraic two-cycles solve g^2-g+1=0; this
    # polynomial has no F_2 root, hence no Q_2 odd-unit realization.
    require(all((a*a+a+1) % 2 != 0 for a in (0, 1)), "Phi_6 has no F2 root")

    # Bounded residue control for the physical unit obstruction on the
    # period-three P^1: nu_2(t+1)=nu_2(1-t)=1 never occurs.
    residue_records = []
    for residue in range(16):
        vp = v2_integer(residue+1)
        vm = v2_integer(1-residue)
        residue_records.append((residue, min(vp, 5), min(vm, 5)))
        require(not (vp == 1 and vm == 1), ("unit obstruction", residue))

    example_t = sp.Integer(2)
    example_g = tuple(sp.factor(value.subs(t, example_t)) for value in gt)
    example_h = tuple(sp.factor(value.subs(t, example_t)) for value in ht)
    require(recurrence(example_g) == example_h, "rational genuine cycle step one")
    require(recurrence(example_h) == example_g, "rational genuine cycle step two")

    semantic = {
        "M": [[str(item) for item in row] for row in M_h.tolist()],
        "char": str(char_general),
        "char_h1": str(sp.factor(M1.charpoly(lam).as_expr())),
        "char_primitive": str(sp.factor(char_cyc)),
        "g": tuple(map(str, gt)),
        "hprofile": tuple(map(str, ht)),
        "example": (tuple(map(str, example_g)), tuple(map(str, example_h))),
        "direct_quotients": tuple(map(str, direct_quotients)),
        "residues": tuple(residue_records),
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_TWO_SCALE_PERIOD3_TWISTED_CLASSIFICATION")
    print("status=FINITE-EXACT_SCRATCH_ALGEBRA;RULE30_OPEN")
    print("product_law=p_to_minus_p_squared;two_cycle_holonomy=h_cubed_1")
    print("twisted_matrix=" + str(M_h.tolist()).replace(" ", ""))
    print("characteristic=" + str(char_general).replace(" ", ""))
    print("ordinary_holonomy_characteristic=(lambda-4)*(lambda-1)^2")
    print("ordinary_genuine_family_parameter=t_in_P1_minus_{infinity,-1,0,1}")
    print("ordinary_cycle=R(g(t))=g(-t);R^2(g(t))=g(t)")
    print("ordinary_direct_elimination=(a*b-a+1)*((a+1),(b+1),(a*b-1))")
    print("rational_control_t2_g=" + str(example_g).replace(" ", ""))
    print("rational_control_t2_Rg=" + str(example_h).replace(" ", ""))
    print("primitive_holonomy_eigenvalues=(h,h^2,1);all_eigenlines_zero_coordinate")
    print("physical_gate=no_t_has_v2(t+1)=v2(1-t)=1;no_all_odd_amplitudes")
    print("constant_line=g=(-1,-1,-1);gaps=(1,1);forbidden_consecutive_one")
    print("mealy_gate=not_reached_after_valuation_failure")
    print("semantic_sha256=" + digest)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()

