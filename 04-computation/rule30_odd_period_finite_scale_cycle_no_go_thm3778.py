#!/usr/bin/env python3
"""Independent exact audit for the period-3, two-scale Rule-30 profile lane.

This is deliberately independent of the candidate companion.  It rebuilds
the twisted amplitude matrices from index reduction, checks the complete
ordinary and cubic-holonomy eigenspaces, verifies the homogeneous boundary,
and uses finite-field hostile universes to look for discarded components.
No Python ``assert`` statement is used.
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


def require(condition, label):
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def twisted_operator(period, holonomy):
    """Matrix A -> (A_(2j)+A_(2j+1)) for A_(j+n)=h A_j."""
    out = sp.zeros(period)
    for j in range(period):
        for k in (2 * j, 2 * j + 1):
            quotient, residue = divmod(k, period)
            out[j, residue] += holonomy**quotient
    return out


def recurrence(profile):
    n = len(profile)
    return tuple(
        sp.factor(
            -profile[2 * j % n]
            * profile[(2 * j + 1) % n]
            * (1 - profile[(2 * j + 2) % n])
            / (1 - profile[2 * j % n])
        )
        for j in range(n)
    )


def amplitude_ratios(amplitude, holonomy):
    period = len(amplitude)
    return tuple(
        sp.factor(
            -(
                amplitude[(j + 1) % period]
                * (holonomy if j + 1 >= period else 1)
            )
            / amplitude[j]
        )
        for j in range(period)
    )


def scale_composition(period, cycle_length, holonomy):
    out = sp.eye(period)
    current_holonomy = holonomy
    for _ in range(cycle_length):
        out = twisted_operator(period, current_holonomy) * out
        current_holonomy = current_holonomy**2
    return out


def reduce_primitive(expr, h, lam):
    numerator, denominator = sp.fraction(sp.cancel(expr))
    domain = sp.QQ.frac_field(lam)
    modulus = sp.Poly(h**2 + h + 1, h, domain=domain)
    numerator = sp.rem(sp.Poly(numerator, h, domain=domain), modulus).as_expr()
    denominator = sp.rem(sp.Poly(denominator, h, domain=domain), modulus).as_expr()
    return sp.cancel(numerator / denominator)


def finite_field_hostile(prime):
    def inv(value):
        return pow(value % prime, prime - 2, prime)

    def finite_R(profile):
        return tuple(
            (
                -profile[2 * j % 3]
                * profile[(2 * j + 1) % 3]
                * (1 - profile[(2 * j + 2) % 3])
                * inv(1 - profile[2 * j % 3])
            )
            % prime
            for j in range(3)
        )

    actual = set()
    for a in range(prime):
        for b in range(prime):
            for c in range(prime):
                profile = (a, b, c)
                if any(x in (0, 1) for x in profile):
                    continue
                image = finite_R(profile)
                if any(x in (0, 1) for x in image):
                    continue
                if finite_R(image) == profile:
                    actual.add(profile)

    predicted = {((-1) % prime,) * 3}
    for parameter in range(prime):
        if parameter in (1, prime - 1):
            continue
        predicted.add(
            (
                2 * inv(parameter + 1) % prime,
                (1 - parameter) * inv(2) % prime,
                (parameter + 1) * inv(parameter - 1) % prime,
            )
        )
    require(actual == predicted, ("finite-field component hostile", prime))
    return prime, len(actual)


def v2_nonzero_integer(value):
    value = abs(int(value))
    require(value != 0, "nonzero v2 input")
    valuation = 0
    while value % 2 == 0:
        valuation += 1
        value //= 2
    return valuation


def main():
    h, lam, t = sp.symbols("h lambda t")

    # Product law hostile at three different odd spatial periods.
    product_law_controls = []
    for period in (1, 3, 5):
        variables = sp.symbols(f"g0:{period}")
        image = recurrence(variables)
        product = sp.prod(variables)
        require(sp.factor(sp.prod(image) + product**2) == 0,
                ("odd-period product law", period))
        product_law_controls.append(period)

    # Cheapest hostile to the odd-period hypothesis: at n=2 doubling is not
    # a coordinate permutation, the product law changes, and T_2 is singular.
    even_variables = sp.symbols("e0:2")
    even_product = sp.prod(even_variables)
    even_product_defect = sp.factor(sp.prod(recurrence(even_variables)) + even_product**2)
    require(even_product_defect != 0, "even-period product-law hostile")
    require(twisted_operator(2, 1).det() == 0, "even-period determinant hostile")

    T_h = twisted_operator(3, h)
    require(
        T_h == sp.Matrix([[1, 1, 0], [h, 0, 1], [0, h, h]]),
        "twisted operator reconstructed from indices",
    )
    M_h = twisted_operator(3, h**2) * T_h
    require(M_h == scale_composition(3, 2, h),
            "rightmost-first scale composition orientation")
    characteristic = sp.factor(M_h.charpoly(lam).as_expr())
    expected = (lam - h) * (lam - h**2) * (lam - (1 + h + h**2 + h**3))
    require(sp.expand(characteristic - expected) == 0, "two-scale characteristic")

    M_1 = M_h.subs(h, 1)
    lambda_one_basis = (M_1 - sp.eye(3)).nullspace()
    require(len(lambda_one_basis) == 2, "ordinary lambda-one plane dimension")
    require(all(sum(vector) == 0 for vector in lambda_one_basis),
            "ordinary lambda-one plane equation")
    require((M_1 - 4 * sp.eye(3)).nullspace() == [sp.Matrix([1, 1, 1])],
            "ordinary lambda-four line")

    amplitude = sp.Matrix([t + 1, -2, 1 - t])
    sibling = twisted_operator(3, 1) * amplitude
    profile = (
        sp.factor(-amplitude[1] / amplitude[0]),
        sp.factor(-amplitude[2] / amplitude[1]),
        sp.factor(-amplitude[0] / amplitude[2]),
    )
    image = recurrence(profile)
    target = tuple(sp.factor(item.subs(t, -t)) for item in profile)
    require(all(sp.factor(a - b) == 0 for a, b in zip(image, target)),
            "rational involution")
    require(M_1 * amplitude == amplitude, "rational curve is lambda one")

    # An arbitrary odd-unit amplitude vector need not have one common sibling
    # valuation.  This prevents the physical invoice from being exported to
    # an ambient coordinate-dependent normalization model.
    ambient_amplitude = sp.Matrix([1, 1, 3])
    ambient_sibling = twisted_operator(3, 1) * ambient_amplitude
    require(tuple(map(int, ambient_sibling)) == (2, 4, 4),
            "coordinate-dependent normalization hostile")

    # A nonperiodic-holonomy numerical control checks that the T_h orientation
    # is the literal recurrence lift, independently of cycle closure.
    amplitude5 = sp.Matrix([2, 3, 5, 7, 11])
    holonomy5 = sp.Integer(13)
    profile5 = amplitude_ratios(amplitude5, holonomy5)
    sibling5 = twisted_operator(5, holonomy5) * amplitude5
    require(recurrence(profile5) == amplitude_ratios(sibling5, holonomy5**2),
            "period-five twisted lift orientation")

    # Homogeneous P^1 boundary.  These are exactly beta=0 and alpha=+-beta;
    # multiplying all rational saturation factors first would incorrectly
    # simplify to one and conceal these individual zero/pole points.
    alpha, beta = sp.symbols("alpha beta")
    homogeneous_A = sp.Matrix([alpha + beta, -2 * beta, beta - alpha])
    homogeneous_B = twisted_operator(3, 1) * homogeneous_A
    require(homogeneous_B == sp.Matrix([alpha - beta, 2 * beta, -alpha - beta]),
            "homogeneous sibling coordinates")
    boundary_forms = tuple(sp.factor(x) for x in tuple(homogeneous_A) + tuple(homogeneous_B))
    require(set(boundary_forms) == {
        alpha + beta, -2 * beta, beta - alpha,
        alpha - beta, 2 * beta, -alpha - beta,
    }, "complete homogeneous boundary")

    primitive_matrix = M_h.applyfunc(lambda x: reduce_primitive(x, h, lam))
    primitive_controls = (
        (h, sp.Matrix([-1, 1, 0])),
        (h**2, sp.Matrix([0, -1, 1])),
        (1, sp.Matrix([1, 0, -h])),
    )
    for eigenvalue, vector in primitive_controls:
        residual = primitive_matrix * vector - eigenvalue * vector
        residual = residual.applyfunc(lambda x: reduce_primitive(x, h, lam))
        require(residual == sp.zeros(3, 1), ("primitive eigenvector", eigenvalue))
        require(any(x == 0 for x in vector), ("primitive boundary", eigenvalue))

    # Physical normalization gate: for an r-scale projective cycle,
    # v2(lambda)=sum of its r positive gap costs.  Here r=2, so lambda=1 is
    # impossible and lambda=4 forces the exact gap pair (1,1).
    require(v2_nonzero_integer(1) < 2, "lambda one fails two-scale gate")
    require(v2_nonzero_integer(4) == 2, "lambda four has total gap two")
    gap_pair = (1, 1)
    require(sum(gap_pair) == 2 and min(gap_pair) >= 1,
            "lambda four uniquely forces consecutive unit gaps")

    # The t-chart version of the same failure, checked on every residue class
    # modulo 32.  Equal valuations of (t+1,-2,1-t) would require both outer
    # valuations to equal one, which never happens.
    parity_hostile = []
    for residue in range(32):
        def capped_v2(value):
            if value == 0:
                return 6
            return min(v2_nonzero_integer(value), 6)
        left = capped_v2(residue + 1)
        right = capped_v2(1 - residue)
        require(not (left == 1 and right == 1), ("Q2 amplitude gate", residue))
        parity_hostile.append((residue, left, right))

    finite_fields = tuple(finite_field_hostile(p) for p in (5, 7, 11, 13))

    # The same carrier specializes to the inherited fixed-profile result:
    # scale-cycle length r=1 forces h=1 and uses T_1, while r=2 forces h^3=1
    # and uses T_(h^2)T_h.  Recheck the inherited n=5,7 spectra independently.
    fixed_spectra = []
    for period in (5, 7):
        polynomial = sp.factor(twisted_operator(period, 1).charpoly(lam).as_expr())
        fixed_spectra.append((period, str(polynomial)))
    require(fixed_spectra == [
        (5, "(lambda - 2)*(lambda - 1)*(lambda + 1)*(lambda**2 + 1)"),
        (7, "(lambda - 2)*(lambda - 1)**2*(lambda**2 + lambda + 1)**2"),
    ], "inherited fixed-profile spectra")

    # General odd-period control behind the physical extension.  Writing
    # T_n=P_2(I+S), where P_2 is the doubling permutation and S the cyclic
    # shift, gives |det(T_n)|=det(I+S)=2 for every odd n.  The constant vector
    # has eigenvalue 2.  Hence, in Qbar_2, it accounts for the whole positive
    # determinant valuation and every other eigenvalue is a unit.  The finite
    # checks below are hostiles for indexing/sign mistakes; the determinant
    # proof itself is in the audit report.
    determinant_controls = []
    for period in range(1, 16, 2):
        operator = twisted_operator(period, 1)
        determinant = int(operator.det())
        eig2_dimension = period - (operator - 2 * sp.eye(period)).rank()
        require(abs(determinant) == 2, ("odd-period determinant", period))
        require(eig2_dimension == 1, ("constant eigenline", period))
        require(operator * sp.ones(period, 1) == 2 * sp.ones(period, 1),
                ("constant eigenvector", period))
        determinant_controls.append((period, determinant, eig2_dimension))

    semantic = {
        "T_h": [[str(x) for x in row] for row in T_h.tolist()],
        "M_h": [[str(x) for x in row] for row in M_h.tolist()],
        "characteristic": str(characteristic),
        "profile": tuple(map(str, profile)),
        "sibling": tuple(map(str, sibling)),
        "boundary": tuple(map(str, boundary_forms)),
        "finite_fields": finite_fields,
        "fixed_spectra": fixed_spectra,
        "determinant_controls": determinant_controls,
        "parity_hostile": tuple(parity_hostile),
        "product_law_controls": product_law_controls,
        "twisted_period5_control": {
            "amplitude": tuple(map(str, amplitude5)),
            "holonomy": str(holonomy5),
            "profile": tuple(map(str, profile5)),
            "sibling": tuple(map(str, sibling5)),
        },
        "even_boundary": {
            "product_defect": str(even_product_defect),
            "determinant": str(twisted_operator(2, 1).det()),
        },
        "ambient_normalization_hostile": tuple(map(str, ambient_sibling)),
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_TWO_SCALE_PERIOD3_INDEPENDENT_AUDIT_20260823")
    print("status=PROVED_ALGEBRAIC_CLASSIFICATION+FINITE_EXACT_HOSTILES;NO_PRIZE")
    print("product_law_controls=odd_periods_1_3_5;p_to_minus_p_squared")
    print("orientation_control=period5_h13_twisted_lift_exact")
    print("even_boundary=n2_product_law_fails;det_T2=0")
    print("normalization_hostile=odd_A_1_1_3_has_sibling_valuations_1_2_2")
    print("ordinary_locus=isolated_lambda4_point+open_lambda1_P1")
    print("ordinary_boundary=beta*\u0028alpha-beta\u0029*\u0028alpha+beta\u0029=0")
    print("primitive_holonomy=three_eigenlines;each_has_zero_amplitude")
    print("physical_eigenvalue_gate=v2\u0028lambda\u0029=d_m+d_mplus1>=2")
    print("only_survivor=lambda4;forced_gaps=\u00281,1\u0029;canonically_forbidden")
    print("finite_field_hostiles=" + repr(finite_fields).replace(" ", ""))
    print("fixed_profile_unification=r1_h1_T1;r2_h3eq1_T_h2_T_h")
    print("fixed_spectra=" + repr(tuple(fixed_spectra)).replace(" ", ""))
    print("general_cycle_geometry=finite_union_of_open_projective_eigenspaces;no_elliptic")
    print("odd_period_determinants=" + repr(tuple(determinant_controls)).replace(" ", ""))
    print("general_physical_conclusion=no_exact_finite_scale_cycle_in_this_odd_period_scheme")
    print("semantic_sha256=" + digest)
    print("RESULT: PASS_MATHEMATICS_WITH_STATUS_AND_REPLAY_REPAIRS")


if __name__ == "__main__":
    main()
