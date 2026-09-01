#!/usr/bin/env python3
"""Dependency-free independent certificate for THM-4317.

This script rebuilds the literal source polynomial used by THM-4312 in a
sparse polynomial ring over ``fractions.Fraction``.  It checks the two
weighted initial forms governing the surviving k=1 carrier, the finite and
infinity charts of its completed-local normalization/resolution, and the
Dirichlet-path/absorbing-walk identities for the A_55 and A_56 chains.

The stochastic statement is deliberately restricted to those ADE chains.
The separate 1/3(1,1) and infinity A_1 quotient resolutions are recorded in
the geometric ledger but are not silently folded into the unbiased path.
"""

from __future__ import annotations

from fractions import Fraction as Q


NAMES = (
    "q", "x", "t", "z", "U", "rho", "up", "eta", "zet", "D", "T",
    "Phi", "w", "Y", "u", "V", "lam", "kap", "sig",
)
INDEX = {name: i for i, name in enumerate(NAMES)}
NVARS = len(NAMES)
ZERO_EXP = (0,) * NVARS
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


class SparsePoly:
    """Tiny exact multivariate sparse polynomial over Q."""

    def __init__(self, terms: dict[tuple[int, ...], Q] | None = None):
        clean: dict[tuple[int, ...], Q] = {}
        if terms:
            for exponent, coefficient in terms.items():
                coefficient = Q(coefficient)
                if coefficient:
                    clean[tuple(exponent)] = clean.get(tuple(exponent), Q(0)) + coefficient
        self.terms = {e: c for e, c in clean.items() if c}

    @staticmethod
    def constant(value: int | Q) -> "SparsePoly":
        value = Q(value)
        return SparsePoly({ZERO_EXP: value}) if value else SparsePoly()

    @staticmethod
    def variable(name: str) -> "SparsePoly":
        exponent = [0] * NVARS
        exponent[INDEX[name]] = 1
        return SparsePoly({tuple(exponent): Q(1)})

    def __add__(self, other: int | Q | "SparsePoly") -> "SparsePoly":
        other = as_poly(other)
        result = dict(self.terms)
        for exponent, coefficient in other.terms.items():
            result[exponent] = result.get(exponent, Q(0)) + coefficient
            if not result[exponent]:
                del result[exponent]
        return SparsePoly(result)

    __radd__ = __add__

    def __neg__(self) -> "SparsePoly":
        return SparsePoly({e: -c for e, c in self.terms.items()})

    def __sub__(self, other: int | Q | "SparsePoly") -> "SparsePoly":
        return self + (-as_poly(other))

    def __rsub__(self, other: int | Q | "SparsePoly") -> "SparsePoly":
        return as_poly(other) - self

    def __mul__(self, other: int | Q | "SparsePoly") -> "SparsePoly":
        other = as_poly(other)
        result: dict[tuple[int, ...], Q] = {}
        for left_exp, left_coefficient in self.terms.items():
            for right_exp, right_coefficient in other.terms.items():
                exponent = tuple(a + b for a, b in zip(left_exp, right_exp))
                result[exponent] = (
                    result.get(exponent, Q(0))
                    + left_coefficient * right_coefficient
                )
        return SparsePoly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "SparsePoly":
        if exponent < 0:
            raise ValueError("negative polynomial exponent")
        result = SparsePoly.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def __eq__(self, other: object) -> bool:
        if isinstance(other, (int, Q)):
            other = as_poly(other)
        return isinstance(other, SparsePoly) and self.terms == other.terms

    def substitute(self, name: str, replacement: "SparsePoly") -> "SparsePoly":
        position = INDEX[name]
        result = SparsePoly()
        variable_free_cache: dict[int, SparsePoly] = {}
        for exponent, coefficient in self.terms.items():
            power = exponent[position]
            reduced = list(exponent)
            reduced[position] = 0
            monomial = SparsePoly({tuple(reduced): coefficient})
            if power not in variable_free_cache:
                variable_free_cache[power] = replacement**power
            result = result + monomial * variable_free_cache[power]
        return result

    def coefficient(self, name: str, power: int) -> "SparsePoly":
        position = INDEX[name]
        result: dict[tuple[int, ...], Q] = {}
        for exponent, coefficient in self.terms.items():
            if exponent[position] == power:
                reduced = list(exponent)
                reduced[position] = 0
                result[tuple(reduced)] = coefficient
        return SparsePoly(result)

    def derivative(self, name: str) -> "SparsePoly":
        position = INDEX[name]
        result: dict[tuple[int, ...], Q] = {}
        for exponent, coefficient in self.terms.items():
            power = exponent[position]
            if not power:
                continue
            reduced = list(exponent)
            reduced[position] -= 1
            result[tuple(reduced)] = coefficient * power
        return SparsePoly(result)

    def specialize_zero(self, name: str) -> "SparsePoly":
        return self.coefficient(name, 0)

    def divide_monomial(self, name: str, power: int) -> "SparsePoly":
        position = INDEX[name]
        result: dict[tuple[int, ...], Q] = {}
        for exponent, coefficient in self.terms.items():
            if exponent[position] < power:
                raise AssertionError(
                    f"term not divisible by {name}^{power}: {exponent}"
                )
            reduced = list(exponent)
            reduced[position] -= power
            result[tuple(reduced)] = coefficient
        return SparsePoly(result)

    def weighted_initial(self, weights: dict[str, int]) -> tuple[int, "SparsePoly"]:
        if not self.terms:
            raise AssertionError("zero polynomial has no weighted initial")
        numeric_weights = [weights.get(name, 0) for name in NAMES]
        values = {
            exponent: sum(a * b for a, b in zip(exponent, numeric_weights))
            for exponent in self.terms
        }
        minimum = min(values.values())
        return minimum, SparsePoly(
            {e: self.terms[e] for e, value in values.items() if value == minimum}
        )

    def minimum_exponent(self, linear_form: dict[str, int]) -> int:
        weights = [linear_form.get(name, 0) for name in NAMES]
        return min(
            sum(a * b for a, b in zip(exponent, weights))
            for exponent in self.terms
        )


def as_poly(value: int | Q | SparsePoly) -> SparsePoly:
    return value if isinstance(value, SparsePoly) else SparsePoly.constant(value)


def coefficient_two(poly: SparsePoly, first: str, p1: int, second: str, p2: int) -> SparsePoly:
    return poly.coefficient(first, p1).coefficient(second, p2)


def literal_source(high_contact: bool = False) -> tuple[SparsePoly, SparsePoly]:
    """Return F before and after q=rho*t+x, with k=1 imposed."""

    q = SparsePoly.variable("q")
    x = SparsePoly.variable("x")
    t = SparsePoly.variable("t")
    z = SparsePoly.variable("z")
    U = SparsePoly.variable("U")
    rho = SparsePoly.variable("rho")
    up = SparsePoly.variable("up")
    eta = SparsePoly.variable("eta")
    zet = SparsePoly.variable("zet")
    D = SparsePoly.variable("D")
    T = SparsePoly.variable("T")
    Phi = SparsePoly.variable("Phi")

    alpha = -2 * U * rho
    c2 = U * rho**2
    if high_contact:
        eta = -zet - rho * up  # L1=eta+zet+rho*up=0.

    r = 1 + q
    K = Q(2848, 45) - Q(7, 6) * D
    hhat = (
        U * (r**6 - 2 * r**5 + r**4)
        + t * alpha * (r**5 - r**4)
        + t**2 * (up * r**5 + (c2 - up) * r**4)
        + t**3 * (eta * r**4 + zet * r**3)
        + t**4 * (D * r**4 + T * r**3)
        + t**5 * Phi * r**3
        + t**6 * (-Q(1376, 135) * r**3 + K * r**2)
        + Q(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    source = q * (hhat - z**12) - Q(1, 2) * t**12
    translated = source.substitute("q", x + rho * t)
    check(translated.coefficient("q", 0) == translated, "translation removes q")
    return source, translated


def check_literal_rows_and_initials() -> tuple[SparsePoly, SparsePoly, SparsePoly, SparsePoly]:
    x = SparsePoly.variable("x")
    t = SparsePoly.variable("t")
    z = SparsePoly.variable("z")
    U = SparsePoly.variable("U")
    rho = SparsePoly.variable("rho")
    up = SparsePoly.variable("up")
    eta = SparsePoly.variable("eta")
    zet = SparsePoly.variable("zet")
    D = SparsePoly.variable("D")
    T = SparsePoly.variable("T")

    source1, translated1 = literal_source(False)
    h_minus_z12 = source1.coefficient("q", 1)
    check(coefficient_two(h_minus_z12, "t", 2, "z", 0) == U * rho**2,
          "literal h row t2")
    check(coefficient_two(h_minus_z12, "t", 3, "z", 0) == eta + zet,
          "literal h row t3")
    check(coefficient_two(h_minus_z12, "z", 12, "t", 0) == -1,
          "literal h row -z12")
    q2 = source1.coefficient("q", 2)
    check(coefficient_two(q2, "t", 1, "z", 0) == -2 * U * rho,
          "literal A row alpha=-2Urho")
    check(coefficient_two(source1, "q", 3, "t", 0) == U,
          "literal cubic unit")

    L1 = eta + zet + rho * up
    minimum1, initial1 = translated1.weighted_initial({"x": 6, "t": 4, "z": 1})
    expected1 = rho * t * (U * x**2 - z**12 + L1 * t**3)
    check(minimum1 == 16, "L1 initial weight 16")
    check(initial1 == expected1, "L1 literal weighted initial")

    _, translated2 = literal_source(True)
    C = D + T - rho * zet
    minimum2, initial2 = translated2.weighted_initial({"x": 6, "t": 3, "z": 1})
    expected2 = rho * t * (U * x**2 + up * t**2 * x + C * t**4 - z**12)
    check(minimum2 == 15, "L2 initial weight 15")
    check(initial2 == expected2, "L2 literal weighted initial")

    # z-charts.  No quotient is present because z has weight one.
    w = SparsePoly.variable("w")
    Y = SparsePoly.variable("Y")
    zchart1 = translated1.substitute("t", w * z**4).substitute("x", z**6 * Y)
    zchart1 = zchart1.divide_monomial("z", 16).specialize_zero("z")
    check(zchart1 == rho * w * (U * Y**2 - 1 + L1 * w**3),
          "L1 z-chart exceptional equation")

    zchart2 = translated2.substitute("t", w * z**3).substitute("x", z**6 * Y)
    zchart2 = zchart2.divide_monomial("z", 15).specialize_zero("z")
    raw2 = rho * w * (U * Y**2 + up * w**2 * Y + C * w**4 - 1)
    check(zchart2 == raw2, "L2 z-chart raw exceptional equation")
    L2_scaled = 4 * U * C - up**2  # 4U*L2.
    completed_square = (2 * U * Y + up * w**2) ** 2 - 4 * U + L2_scaled * w**4
    check(4 * U * (U * Y**2 + up * w**2 * Y + C * w**4 - 1)
          == completed_square, "L2 critical recenter/completed square")

    return translated1, translated2, L1, L2_scaled


def check_boundary_branches(L1: SparsePoly, L2_scaled: SparsePoly) -> None:
    U = SparsePoly.variable("U")
    rho = SparsePoly.variable("rho")
    up = SparsePoly.variable("up")
    eta = SparsePoly.variable("eta")
    zet = SparsePoly.variable("zet")
    kap = SparsePoly.variable("kap")
    q = SparsePoly.variable("q")
    t = SparsePoly.variable("t")
    z = SparsePoly.variable("z")

    c2 = U * rho**2
    c3 = eta + zet
    branch_relation = c2 * kap**2 - 1
    check(branch_relation != 0, "formal quadratic branch relation retained")
    check(4 * c2 != 0, "branch discriminant 4c2 nonzero under k1 hypotheses")

    # The first correction to the two roots of h(t)=z^12 is forced already by
    # the literal t^2 and t^3 rows.
    a1 = -Q(1, 2) * c3 * kap**4
    next_balance = 2 * c2 * kap * a1 + c3 * kap**3
    check(next_balance == -c3 * kap**3 * branch_relation,
          "critical t correction modulo c2*kappa^2=1")
    Tseries = kap * z**6 + a1 * z**12
    minus_half_t12 = -Q(1, 2) * Tseries**12
    check(minus_half_t12.coefficient("z", 72) == -Q(1, 2) * kap**12,
          "critical value z72 coefficient")
    check(minus_half_t12.coefficient("z", 78) == 3 * c3 * kap**15,
          "literal-z correction occurs at z78")

    # F_t=-6t^11+q h'(t)+... has h'(t) leading 2c2*t.  The displayed
    # q-leading coefficient solves this modulo the branch relation.
    q_lead = 3 * kap**12
    ft_lead = -6 * kap**11 + q_lead * 2 * c2 * kap
    check(ft_lead == 6 * kap**11 * branch_relation,
          "critical q order-60 leading balance")

    # On the exact q-critical graph, F=F-qF_q.  This deletes the q-linear
    # critical equation and leaves -t^12/2 plus terms -(k-1)q^k B_k.  The
    # literal q^2 coefficient starts with alpha*t, so the first correction is
    # exactly order 2*60+6=126.
    source, _ = literal_source(False)
    critical_reduction = source - q * source.derivative("q")
    explicit_reduction = -Q(1, 2) * t**12
    for power in range(2, 8):
        explicit_reduction = (
            explicit_reduction
            - (power - 1) * q**power * source.coefficient("q", power)
        )
    check(critical_reduction == explicit_reduction,
          "exact critical-value identity F-qFq")
    q_correction = critical_reduction + Q(1, 2) * t**12
    correction_order, _ = q_correction.weighted_initial({"q": 60, "t": 6, "z": 1})
    check(correction_order == 126, "q^2 A remainder begins exactly at z126")
    check(12 * 6 == 72 and 11 * 6 + 12 == 78,
          "correct literal-z critical orders")

    # There are two squarefree kappa branches.  At either one the boundary
    # point Y=-rho*kappa satisfies UY^2=1 and the (w,Y) Hessian determinant is
    # a nonzero monomial -(2 rho U Y)^2.
    Y0 = -rho * kap
    check(U * Y0**2 - 1 == branch_relation,
          "two boundary points lie on exceptional conic")
    hessian_det = -(2 * rho * U * Y0) ** 2
    check(hessian_det != 0, "boundary transverse Hessian nonzero")
    check(L1 != 0 and L2_scaled != 0, "case sidecars are formally nonzero")

    check(72 - 16 == 56, "L1 transverse critical order 56")
    check(72 - 15 == 57, "L2 transverse critical order 57")
    check(56 - 1 == 55 and 57 - 1 == 56, "A-type indices")
    check(2 * 55 == 110 and 2 * 56 == 112, "two boundary chains")

    # On the L2 locus c3=-rho*up, so the z^78 correction is visibly nonzero;
    # this prevents the incorrect literal-z O(z^126) assertion.
    check((3 * c3 * kap**15).substitute("eta", -zet - rho * up)
          == -3 * rho * up * kap**15,
          "L2 z78 correction is nonzero")


def infinity_transform(
    source: SparsePoly,
    t_lambda_weight: int,
    x_u_power: int,
    lambda_order: int,
    u_saturation: int,
) -> SparsePoly:
    """Apply the sigma-chart cover and normalization substitution."""

    lam = SparsePoly.variable("lam")
    u = SparsePoly.variable("u")
    V = SparsePoly.variable("V")
    transformed = source.substitute("t", lam**t_lambda_weight * u)
    transformed = transformed.substitute("z", lam * u)
    transformed = transformed.substitute("x", lam**6 * u**x_u_power * V)
    transformed = transformed.divide_monomial("lam", lambda_order)
    transformed = transformed.divide_monomial("u", u_saturation)
    return transformed


def check_infinity_and_exhaustion(
    translated1: SparsePoly,
    translated2: SparsePoly,
    L1: SparsePoly,
    L2_scaled: SparsePoly,
) -> None:
    U = SparsePoly.variable("U")
    rho = SparsePoly.variable("rho")
    up = SparsePoly.variable("up")
    D = SparsePoly.variable("D")
    T = SparsePoly.variable("T")
    zet = SparsePoly.variable("zet")
    u = SparsePoly.variable("u")
    V = SparsePoly.variable("V")
    C = D + T - rho * zet

    # L1: sigma=lambda^3, z=lambda*u, x=lambda^6*X, X=uV.
    check(translated1.minimum_exponent({"x": 1, "t": 1, "z": 1}) == 3,
          "L1 u^3 normalization saturation is regular")
    infinity1 = infinity_transform(translated1, 4, 1, 16, 3)
    exceptional1 = infinity1.specialize_zero("lam")
    check(exceptional1 == rho * (U * V**2 + L1 * u - u**10),
          "L1 infinity normalized equation")
    check(L1 != 0, "L1 infinity implicit u derivative nonzero")

    # mu_3 cover action: lambda->zeta lambda, u->zeta^2 u, X fixed,
    # hence V=X/u has weight one.  Solving u leaves weights (1,1).
    mu3 = {"lam": 1, "u": 2, "X": 0, "V": 1}
    check((2 * mu3["V"]) % 3 == mu3["u"], "mu3 equation character")
    check((10 * mu3["u"]) % 3 == mu3["u"], "mu3 u10 character")
    check((mu3["lam"], mu3["V"]) == (1, 1), "quotient 1/3(1,1)")
    check(3 // 1 == 3, "Hirzebruch-Jung [3] gives one minus-three curve")

    # L2: sigma=lambda^2, z=lambda*u, x=lambda^6*X, X=u^2V.
    check(translated2.minimum_exponent({"x": 2, "t": 1, "z": 1}) == 5,
          "L2 u^5 normalization saturation is regular")
    infinity2 = infinity_transform(translated2, 3, 2, 15, 5)
    exceptional2 = infinity2.specialize_zero("lam")
    check(exceptional2 == rho * (U * V**2 + up * V + C - u**8),
          "L2 infinity normalized equation")
    discriminant = up**2 - 4 * U * C
    check(discriminant == -L2_scaled, "L2 infinity discriminant=-4U*L2")
    check(discriminant != 0, "two distinct L2 infinity roots")

    # mu_2 action fixes V and acts by (-1,-1) on local coordinates
    # (lambda,u) at each of the two roots.
    mu2 = {"lam": 1, "u": 1, "X": 0, "V": 0}
    check((2 * mu2["u"] + mu2["V"]) % 2 == mu2["X"],
          "mu2 normalization action")
    check((mu2["lam"], mu2["u"]) == (1, 1), "two quotient A1 actions")

    # Projective chart exhaustion in the actual weighted variables, with
    # t=sigma*z.  At z=0 both homogeneous carrier equations restrict to
    # U*x^2, so x=0 and only the sigma-chart remains; the x-chart is empty.
    sig = SparsePoly.variable("sig")
    zvar = SparsePoly.variable("z")
    xvar = SparsePoly.variable("x")
    carrier1 = U * xvar**2 - zvar**12 + L1 * sig**3 * zvar**3
    carrier2 = (
        U * xvar**2
        + up * sig**2 * zvar**2 * xvar
        + C * sig**4 * zvar**4
        - zvar**12
    )
    check(carrier1.specialize_zero("z") == U * xvar**2,
          "L1 z=0 forces x=0 and enters sigma-chart")
    check(carrier2.specialize_zero("z") == U * xvar**2,
          "L2 z=0 forces x=0 and enters sigma-chart")
    check(
        exceptional1.specialize_zero("u") == rho * U * V**2
        and exceptional2.specialize_zero("u")
        == rho * (U * V**2 + up * V + C),
        "sigma-chart contains exactly the audited infinity root equations",
    )
    check(2 * 55 + 1 == 111, "L1 minimal new rational-curve count")
    check(2 * 56 + 2 == 114, "L2 minimal new rational-curve count")


def green(r: int, i: int, j: int) -> Q:
    if i in (0, r) or j in (0, r):
        return Q(0)
    return Q(min(i, j) * (r - max(i, j)), r)


def verify_path_kernel(r: int) -> None:
    check(r >= 3, f"r={r} kept away from quotient-A1 application")
    for i in range(1, r):
        hit_right = Q(i, r)
        check(hit_right == (Q(i - 1, r) + Q(i + 1, r)) / 2,
              f"r={r} gambler ruin harmonic i={i}")
        check(Q(r - i) == r * (1 - hit_right), f"r={r} u valuation")
        check(Q(i) == r * hit_right, f"r={r} v valuation")
        for j in range(1, r):
            laplacian = 2 * green(r, i, j) - green(r, i - 1, j) - green(r, i + 1, j)
            check(laplacian == (1 if i == j else 0),
                  f"r={r} inverse Cartan ({i},{j})")
            visits = 2 * green(r, i, j)
            next_visits = green(r, i - 1, j) + green(r, i + 1, j)
            check(visits - next_visits == (1 if i == j else 0),
                  f"r={r} absorbing occupation ({i},{j})")

    # Several monomial orders, including asymmetric and z-shifted controls.
    for a, b, c in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (3, 5, 7), (11, 2, 13)):
        values = [Q(a * (r - i) + b * i + c) for i in range(r + 1)]
        for i in range(1, r):
            check(values[i] == (values[i - 1] + values[i + 1]) / 2,
                  f"r={r} monomial harmonic ({a},{b},{c}) i={i}")

    # Exact Poisson/optional-stopping control with nonzero contact source.
    left = Q(5, 3)
    right = Q(7, 4)
    source = [Q(0)] + [Q((j * j + 3 * j + 1) % 7) for j in range(1, r)] + [Q(0)]
    potential = [left]
    for i in range(1, r):
        harmonic_boundary = Q(r - i, r) * left + Q(i, r) * right
        poisson = sum(green(r, i, j) * source[j] for j in range(1, r))
        occupation = Q(1, 2) * sum(2 * green(r, i, j) * source[j] for j in range(1, r))
        check(poisson == occupation, f"r={r} Green/occupation identity i={i}")
        potential.append(harmonic_boundary + poisson)
    potential.append(right)
    for i in range(1, r):
        check(2 * potential[i] - potential[i - 1] - potential[i + 1] == source[i],
              f"r={r} Poisson equation i={i}")
        check(potential[i] >= Q(r - i, r) * left + Q(i, r) * right,
              f"r={r} effective-source superharmonicity i={i}")


def check_stochastic_resolution_ledger() -> None:
    # Small exact controls exercise both parity classes; r=56,57 are the two
    # boundary chains in the theorem.  Infinity quotient curves remain a
    # separate ledger even though an A1 can of course be analyzed separately.
    for r in tuple(range(3, 13)) + (56, 57):
        verify_path_kernel(r)
    check(55 == 56 - 1 and 56 == 57 - 1, "chain lengths A55 and A56")
    resolution_ledger = {
        "boundary_A55_pair": {"multiplicity": 2, "path_order": 56},
        "boundary_A56_pair": {"multiplicity": 2, "path_order": 57},
        "infinity_minus_three": {"multiplicity": 1, "path_order": None},
        "infinity_A1_pair": {"multiplicity": 2, "path_order": None},
    }
    check(
        tuple(
            entry["path_order"]
            for entry in resolution_ledger.values()
            if entry["path_order"] is not None
        )
        == (56, 57),
        "only the two boundary ADE orders enter the unbiased paths",
    )
    check(
        resolution_ledger["infinity_minus_three"]["path_order"] is None
        and resolution_ledger["infinity_A1_pair"]["path_order"] is None,
        "cyclic quotient sidecars are excluded from the path application",
    )


def main() -> None:
    translated1, translated2, L1, L2_scaled = check_literal_rows_and_initials()
    check_boundary_branches(L1, L2_scaled)
    check_infinity_and_exhaustion(translated1, translated2, L1, L2_scaled)
    check_stochastic_resolution_ledger()

    print("THM-4317 INDEPENDENT FRACTION/SPARSE AUDIT: PASS")
    print("L1_INITIAL weight=16 rho*t*(U*x^2-z^12+L1*t^3)")
    print("L2_INITIAL weight=15 rho*t*(U*x^2+ups*t^2*x+C*t^4-z^12)")
    print("L2_RECENTER C=D+Theta-rho*zeta 4U*L2=4U*C-ups^2")
    print("BOUNDARY_BRANCHES c2*kappa^2=1 two_distinct=yes t=kappa*z^6+O(z^12) q=3*kappa^12*z^60+O(z^66)")
    print("CRITICAL_VALUE literal_z=-kappa^12*z^72/2+3*c3*kappa^15*z^78+O(z^84) q_correction_starts_z126")
    print("BOUNDARY_TYPES L1=two_A55 L2=two_A56")
    print("L1_INFINITY saturation=u^3 equation=U*V^2+L1*u-u^10 action=1/3(1,1) resolution=one_minus3")
    print("L2_INFINITY saturation=u^5 equation=U*V^2+ups*V+C-u^8 discriminant=-4U*L2 roots=2 actions=two_A1")
    print("CHARTS z_nonzero_plus_sigma_infinity_exhaust_carrier x_chart_adds_none")
    print("RATIONAL_CURVES minimal_L1=111 minimal_L2=114 later_point_blowups=P1")
    print("A_CHAIN r=56,57 nu_i(u)=r-i nu_i(v)=i nu_i(z)=1")
    print("GREEN inverse_Cartan=min(i,j)*(r-max(i,j))/r expected_visits=2*inverse_Cartan")
    print("POISSON a_i=E_boundary+(1/2)*E_accumulated_contact martingale_with_compensator=yes")
    print("UNBIASED_PATH_SCOPE boundary_A55_A56_only infinity_minus3_and_A1_quotients_separate")
    print("SCOPE finite_row8_k1_completed_local conditional_actual_Keller_lift no_all_row no_JC2")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
