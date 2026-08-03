#!/usr/bin/env python3
"""Exact exceptional-deck pullback of the JC cofactor sidecars.

On the fixed ``C=c+x``, ``d=k=1`` slice in the two THM-3212 accessory
fields, reconstruct the degree-36 linear-subresultant base field ``A`` and
the degree-two exceptional field

    B=A[t]/(P2*t**2-Q2*t+R2).

The degree-36 polynomial ``linear_a`` and the degree-32 quadratic pivot
``P2`` are kept distinct.  The script evaluates the localized gradients and
the cubic-to-linear elimination cofactors at ``y=-t``, computes their deck
trace/norm packets, identifies the quadratic inverse different, and checks
compatibility with the already derived first-normal velocity.

This is a finite-exact fixed-slice scout.  The inverse-different generator is
not asserted to be a primitive Keller cofactor, and the critical deck itself
is an exact obstruction to a polynomial Keller mate on this slice.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md":
        "514b5c07a70ea0dd0020857af3278008c3df2a03af6826784d96d059a3b26111",
    "04-computation/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py":
        "6f050a583004172f812c3f7729427079d5df45c3a985c2e470b2a0d34ad8f337",
    "05-knowledge/results/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.out":
        "0603517a0e97c1eb3d6b60051cfb02c2a8074ac3468028a36b97e7b8398f5670",
    "04-computation/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py":
        "a719b2582b93a0a6d110b1f13b65e9d54800e8669914da9f21a9371545bbae31",
    "05-knowledge/results/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.out":
        "67067d9448caa6a809520b190208a561ab4cc14517455d6da0eef9210ccce1ff",
}
EXPECTED_DELTA_DIGESTS = {
    "4111": "6b66fda17d31412a565a23649197ca94e9372c7368f01a31f60da61714b328fc",
    "3211": "fa7ce86608c083d3e8e722cf5de425a9a3d3ae3ddbdbc5e829beb6e10fb7187a",
}
EXPECTED_CASE_PACKETS = (
    (
        "83bf3f87c6604485781ed8b44b4b2ba65d8fe8fb45781a5e86593790b72c6e10",
        "a716987b7b8096ba9aacd0b31e520f6d504611976222b17589ee377674323cb8",
        "39a8b14847a5c8981a7dd8c799a8c7deeeb30a0f569daa95fc14db402d422692",
        "262bf1b4077b33afcf4c98eab5a4f9a0ef7f2d2652a8e408d5a600495ef4d47c",
        "5cdc30990101b93ad8f0e1e3f4c66950bd8a4b6df05cfb8a6f701a29d71000ac",
        "837d9ed3890d4a20349c9bb4f67329b4d1648414ccf86a982adc991f32deebf0",
    ),
    (
        "0a81253a6b16e7e46e6f601b96d9561c20aa374faefbd421dd777dc50ab38c2f",
        "5a33fa1b60027e2a971f4d50bcf6b0510eec23ffc7453897f53f604769009b91",
        "7ce4c1cee9ffef1d1804385d24a03e82367d312e11e63b9ec083905e67626d89",
        "8068e4edc0ef4bc20326cfa291caecc0673ced8488f6fd5ecd8baab8e738e842",
        "3248875cd752b42e2154811bdd2fc1c2d0a1ddff89f62ace4ed703245955bbd3",
        "5bf8ba6c20e9792def6e622d4900e8ff658852e61b285b45b2ced51b801b2cf3",
    ),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def element_text(element) -> str:
    return "|".join(
        f"{monomial}:{coefficient}"
        for monomial, coefficient in sorted(element.to_dict().items())
    )


def element_digest(element) -> str:
    return sha256(element_text(element).encode("utf-8")).hexdigest()


def pair_text(element) -> str:
    return element_text(element[0]) + "\n" + element_text(element[1])


def pair_digest(element) -> str:
    return sha256(pair_text(element).encode("utf-8")).hexdigest()


def packet_digest(elements) -> str:
    return sha256("\n--\n".join(pair_text(item) for item in elements).encode("utf-8")).hexdigest()


def inverse_mod(element, modulus):
    inverse, _, common = element.gcdex(modulus)
    require(common.degree() == 0, "requested A-element is a unit")
    return (inverse / common.LC) % modulus


def response_case(name: str):
    u, x = sp.symbols("u x")
    if name == "4111":
        accessory = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    elif name == "3211":
        accessory = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2
    else:
        raise RuntimeError(f"unknown accessory case: {name}")

    field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = field.ext
    ring = field.poly_ring(x)
    X = ring.gens[0]
    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
        extras = (9, 0)
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)

    gamma = -7 * A0
    quadratic = X**2 - alpha * X + accessory_v
    Dsource = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    Esource = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * Dsource * T**2 / gamma**2
    Aresponse = 2 * S * Esource * T / gamma
    boundary = S**3 * T**8 * X**extras[0] * (X - 1) ** extras[1]
    require((V.degree(), Aresponse.degree(), boundary.degree()) == (16, 8, 44), name)
    require(
        2 * V * Aresponse.diff(X) - Aresponse * V.diff(X) == 2 * V,
        f"{name} response identity",
    )
    return {
        "name": name,
        "field": field,
        "ring": ring,
        "X": X,
        "V": V,
        "Aresponse": Aresponse,
        "boundary": boundary,
    }


def linear_row_at_parameter(case, parameter):
    """Return the literal closed-form linear pseudo-remainder at c=parameter."""

    field = case["field"]
    X = case["X"]
    V = case["V"]
    Aresponse = case["Aresponse"]
    derivative_V = V.diff(X)
    C = X + field.convert(parameter)
    P2 = 2 * derivative_V + 8 * V**2
    Q2 = 2 * derivative_V + 12 * V**2
    R2 = V * (4 * V**2 + 8 * C * V**2 + derivative_V * (2 * C + Aresponse))
    first_q = 2 + 4 * C * V
    first_r = V * (2 * C + Aresponse)
    D = 6 * P2 - 4 * Q2
    ell1 = P2 * (P2 * first_q - 4 * R2) - D * Q2
    ell0 = P2**2 * first_r - D * R2
    return ell1, ell0


class QuadraticDeck:
    """Arithmetic in A[t]/(P*t^2-Q*t+R), with A=K[x]/(linear_a)."""

    def __init__(self, modulus, P, Q, R):
        self.modulus = modulus
        self.ring = modulus.ring
        self.P = P % modulus
        self.Q = Q % modulus
        self.R = R % modulus
        self.inverse_P = inverse_mod(self.P, modulus)
        self.trace_t = (self.Q * self.inverse_P) % modulus
        self.norm_t = (self.R * self.inverse_P) % modulus

    def A(self, value):
        return (self.ring.zero, value % self.modulus)

    def t(self):
        return (self.ring.one, self.ring.zero)

    def reduce(self, element):
        return (element[0] % self.modulus, element[1] % self.modulus)

    def add(self, left, right):
        return self.reduce((left[0] + right[0], left[1] + right[1]))

    def neg(self, element):
        return self.reduce((-element[0], -element[1]))

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def scale(self, scalar, element):
        return self.reduce((scalar * element[0], scalar * element[1]))

    def mul(self, left, right):
        leading = left[0] * right[0]
        t_coefficient = (
            leading * self.trace_t + left[0] * right[1] + left[1] * right[0]
        )
        constant = left[1] * right[1] - leading * self.norm_t
        return self.reduce((t_coefficient, constant))

    def conjugate(self, element):
        return self.reduce((-element[0], element[1] + element[0] * self.trace_t))

    def trace(self, element):
        return (2 * element[1] + element[0] * self.trace_t) % self.modulus

    def norm(self, element):
        return self.mul(element, self.conjugate(element))[1]

    def inverse(self, element):
        norm = self.norm(element)
        require(norm != self.ring.zero, "requested B-element is nonzero")
        return self.scale(inverse_mod(norm, self.modulus), self.conjugate(element))

    def equal(self, left, right):
        return self.sub(left, right) == (self.ring.zero, self.ring.zero)


def audit_case(case):
    name = case["name"]
    field = case["field"]
    ring = case["ring"]
    X = case["X"]
    V = case["V"]
    Aresponse = case["Aresponse"]
    boundary = case["boundary"]
    derivative_V = V.diff(X)
    second_derivative_V = derivative_V.diff(X)

    raw_a0, raw_b0 = linear_row_at_parameter(case, 0)
    raw_a1, raw_b1 = linear_row_at_parameter(case, 1)
    raw_a0 = raw_a0.exquo(boundary)
    raw_b0 = raw_b0.exquo(boundary)
    raw_a1 = raw_a1.exquo(boundary)
    raw_b1 = raw_b1.exquo(boundary)
    require(raw_a0 == raw_a1, f"{name} linear coefficient independent of c")
    common_unit = raw_a0.LC
    linear_a = raw_a0 / common_unit
    b0 = raw_b0 / common_unit
    b1 = (raw_b1 - raw_b0) / common_unit
    require(
        (linear_a.degree(), b0.degree(), b1.degree()) == (36, 44, 36),
        f"{name} degree-36 linear pivot packet",
    )
    require(linear_a.factor_list()[1][0][0].degree() == 36, f"{name} irreducible linear_a")
    cstar = (-b0 * inverse_mod(b1, linear_a)) % linear_a
    Cstar = X + cstar

    P2_raw = 2 * derivative_V + 8 * V**2
    Q2_raw = 2 * derivative_V + 12 * V**2
    R2_raw = V * (
        4 * V**2 + 8 * Cstar * V**2 + derivative_V * (2 * Cstar + Aresponse)
    )
    require(
        (P2_raw.degree(), Q2_raw.degree()) == (32, 32),
        f"{name} degree-32 quadratic pivot",
    )
    P2, Q2, R2 = (P2_raw % linear_a, Q2_raw % linear_a, R2_raw % linear_a)
    delta = (Q2**2 - 4 * P2 * R2) % linear_a
    require(
        all(linear_a.gcd(item).degree() == 0 for item in (V, derivative_V, P2, delta)),
        f"{name} deck units",
    )
    require(
        element_digest(delta) == EXPECTED_DELTA_DIGESTS[name],
        f"{name} inherited nonsquare discriminant drift",
    )
    deck = QuadraticDeck(linear_a, P2, Q2, R2)
    zero = (ring.zero, ring.zero)
    one = deck.A(ring.one)
    t = deck.t()
    y = deck.neg(t)

    # The defining exceptional relation is literal, not only a discriminant
    # statement.  It makes y=-t a common root of both gradient cubics.
    F0 = deck.add(
        deck.add(deck.scale(P2, deck.mul(t, t)), deck.scale(-Q2, t)),
        deck.A(R2),
    )
    require(F0 == zero, f"{name} exceptional relation")
    L = deck.add(deck.add(deck.mul(y, y), y), deck.A((Cstar * V) % linear_a))
    R1 = deck.add(
        deck.scale(2, deck.mul(L, deck.add(deck.scale(2, y), one))),
        deck.A((V * Aresponse) % linear_a),
    )
    R2_value = deck.add(
        deck.add(deck.A(V**3 % linear_a), deck.scale(V**2, y)),
        deck.mul(
            L,
            deck.add(deck.scale(-derivative_V, y), deck.A(2 * V**2 % linear_a)),
        ),
    )
    require(R1 == zero and R2_value == zero, f"{name} common exceptional root")

    # Independently evaluate the physical gradient at z=y/V, then verify the
    # triangular change of gradient coordinates
    #   R1=V*P_z,  R2=V^3*P_x-(V'*y/2)*R1.
    # Since V is a unit in A, the actual gradient vanishes on B.  This is the
    # first Keller-mate gate.
    inverse_V = inverse_mod(V % linear_a, linear_a)
    z = deck.scale(inverse_V, y)
    square_base = deck.add(
        deck.add(deck.scale(V, deck.mul(z, z)), z),
        deck.A(Cstar),
    )
    direct_Pz = deck.add(
        deck.scale(
            2,
            deck.mul(square_base, deck.add(deck.scale(2 * V, z), one)),
        ),
        deck.A(Aresponse),
    )
    direct_Px = deck.add(
        deck.add(
            deck.scale(
                2,
                deck.mul(
                    square_base,
                    deck.add(deck.scale(derivative_V, deck.mul(z, z)), one),
                ),
            ),
            deck.scale(Aresponse.diff(X), z),
        ),
        one,
    )
    Pz = deck.scale(inverse_V, R1)
    Px = deck.add(
        deck.scale(inverse_mod(V**3 % linear_a, linear_a), R2_value),
        deck.scale(
            inverse_mod(2 * V**3 % linear_a, linear_a),
            deck.mul(deck.scale(derivative_V, y), R1),
        ),
    )
    require(
        deck.equal(R1, deck.scale(V, direct_Pz))
        and deck.equal(
            R2_value,
            deck.sub(
                deck.scale(V**3, direct_Px),
                deck.scale(derivative_V / 2, deck.mul(y, R1)),
            ),
        ),
        f"{name} localized/physical gradient change",
    )
    require(
        Px == zero and Pz == zero and direct_Px == zero and direct_Pz == zero,
        f"{name} physical gradient pullback",
    )

    # Pull back the universal elimination-cofactor pair.  Here ``Qelim`` is
    # the linear pseudo-division quotient, not the polynomial-map mate.
    D = (6 * P2 - 4 * Q2) % linear_a
    Qelim = deck.add(deck.scale(4 * P2, y), deck.A(D))
    U = deck.sub(deck.A(P2**2 % linear_a), deck.scale(derivative_V, Qelim))
    W = deck.scale(-4, Qelim)
    cofactor_identity = deck.sub(U, deck.scale(derivative_V / 4, W))
    require(
        deck.equal(cofactor_identity, deck.A(P2**2 % linear_a)),
        f"{name} elimination cofactor identity on deck",
    )

    Fprime = deck.add(deck.scale(2 * P2, t), deck.A(-Q2))
    require(
        deck.equal(deck.mul(Fprime, Fprime), deck.A(delta)),
        f"{name} quadratic derivative square",
    )
    require(
        deck.equal(Qelim, deck.add(deck.A(-24 * V**2 % linear_a), deck.scale(-2, Fprime))),
        f"{name} invariant/anti-invariant Qelim split",
    )
    require(deck.equal(deck.conjugate(Fprime), deck.neg(Fprime)), f"{name} deck sign")
    require(
        deck.equal(deck.add(t, deck.conjugate(t)), deck.A(deck.trace_t))
        and deck.equal(deck.mul(t, deck.conjugate(t)), deck.A(deck.norm_t)),
        f"{name} tautological and conjugate root factorization",
    )
    require(
        deck.equal(W, deck.add(deck.A(96 * V**2 % linear_a), deck.scale(8, Fprime)))
        and deck.equal(
            U,
            deck.add(
                deck.A((P2**2 + 24 * derivative_V * V**2) % linear_a),
                deck.scale(2 * derivative_V, Fprime),
            ),
        ),
        f"{name} cofactor invariant/anti-invariant split",
    )

    # The two cofactor values are individually nonzero in the field B.  Their
    # invariant combination descends, but their projective ratio does not.
    norms = tuple(deck.norm(item) for item in (Qelim, U, W))
    require(all(item != ring.zero for item in norms), f"{name} cofactor units in B")
    traces = tuple(deck.trace(item) for item in (Qelim, U, W))
    require(
        traces == (
            (-48 * V**2) % linear_a,
            (2 * P2**2 + 48 * derivative_V * V**2) % linear_a,
            (192 * V**2) % linear_a,
        ),
        f"{name} cofactor traces",
    )
    projective_deck_defect = deck.sub(
        deck.mul(U, deck.conjugate(W)),
        deck.mul(deck.conjugate(U), W),
    )
    require(
        deck.equal(projective_deck_defect, deck.scale(-16 * P2**2, Fprime)),
        f"{name} cofactor ratio separates deck branches",
    )

    # For the monic exceptional polynomial, f'=Fprime/P2.  Its reciprocal
    # generates the inverse different and is a unit because the cover is
    # etale.  It is anti-invariant, trace zero, and is not a supplied Keller
    # cofactor for P.
    monic_derivative = deck.scale(inverse_mod(P2, linear_a), Fprime)
    inverse_different = deck.inverse(monic_derivative)
    require(
        deck.trace(monic_derivative) == ring.zero
        and deck.norm(monic_derivative) == (-delta * inverse_mod(P2**2, linear_a)) % linear_a,
        f"{name} monic different trace and norm",
    )
    require(
        deck.equal(deck.mul(monic_derivative, inverse_different), one)
        and deck.equal(deck.conjugate(inverse_different), deck.neg(inverse_different)),
        f"{name} inverse different unit and deck parity",
    )
    require(
        deck.norm(inverse_different) == (-P2**2 * inverse_mod(delta, linear_a)) % linear_a,
        f"{name} inverse different norm",
    )

    # Reconstruct the first normal packet and show that the quadratic deck
    # involution lifts compatibly to first order.  If T(u)=Q(u)/P(u), then
    # sigma_u(t(u))=T(u)-t(u); hence sigma(dot t)=dot T-dot t.
    a_prime = linear_a.diff(X) % linear_a
    b_x = (b0.diff(X) + cstar * b1.diff(X)) % linear_a
    inverse_a_prime = inverse_mod(a_prime, linear_a)
    inverse_b1 = inverse_mod(b1, linear_a)
    P2_x = (2 * second_derivative_V + 16 * V * derivative_V) % linear_a
    Q2_x = (2 * second_derivative_V + 24 * V * derivative_V) % linear_a
    bracket = (
        4 * V**2 + 8 * Cstar * V**2
        + derivative_V * (2 * Cstar + Aresponse)
    )
    R2_x = (
        derivative_V * bracket
        + V * (
            8 * V * derivative_V + 8 * V**2
            + 16 * Cstar * V * derivative_V
            + second_derivative_V * (2 * Cstar + Aresponse)
            + derivative_V * (2 + Aresponse.diff(X))
        )
    ) % linear_a
    R2_c = (V * (8 * V**2 + 2 * derivative_V)) % linear_a
    F1_2 = (P2_x * inverse_a_prime) % linear_a
    F1_1 = (-Q2_x * inverse_a_prime + R2_c * inverse_b1) % linear_a
    F1_0 = (
        R2_x * inverse_a_prime
        - R2_c * b_x * inverse_a_prime * inverse_b1
    ) % linear_a
    F1_value = deck.add(
        deck.add(deck.scale(F1_2, deck.mul(t, t)), deck.scale(F1_1, t)),
        deck.A(F1_0),
    )
    velocity = deck.neg(deck.mul(F1_value, deck.inverse(Fprime)))
    trace_derivative = (
        (-F1_1 * P2 - Q2 * F1_2) * inverse_mod(P2**2, linear_a)
    ) % linear_a
    require(
        deck.trace(velocity) == trace_derivative,
        f"{name} first-normal trace derivative",
    )
    require(
        deck.equal(deck.conjugate(velocity), deck.sub(deck.A(trace_derivative), velocity)),
        f"{name} first-normal deck compatibility",
    )

    print(
        f"case={name};linear_a_degree=36;P2_degree_before_mod=32;"
        "symbols_distinct=1;B_over_A_degree=2;B_over_K_degree=72;"
        "relative_geometric_fibre_points=2;total_geometric_points_over_K=72"
    )
    print(
        f"case={name};tautological_root=y=-t;"
        f"gradient_pullback_digest={packet_digest((Px, Pz))};"
        "R1=R2=Px=Pz=0_in_B;gradient_ideal_proper=1;"
        "Jac(P,arbitrary_Q)=0_in_B;mate_gate=FAILS_BEFORE_mu(P)"
    )
    print(
        f"case={name};Qelim_digest={pair_digest(Qelim)};"
        f"U_W_digest={packet_digest((U, W))};"
        "U-(Vprime/4)W=P2^2;Qelim=-24V^2-2F0prime;"
        "W=96V^2+8F0prime;U=P2^2+24VprimeV^2+2VprimeF0prime"
    )
    print(
        f"case={name};cofactor_trace_digest="
        f"{sha256('|'.join(element_text(item) for item in traces).encode('utf-8')).hexdigest()};"
        f"cofactor_norm_digest="
        f"{sha256('|'.join(element_text(item) for item in norms).encode('utf-8')).hexdigest()};"
        "Qelim_U_W_units_in_B=1;displayed_descents=(trace,norm,U-(Vprime/4)W);"
        "projective_pair_deck_invariant=0;"
        f"projective_deck_defect_digest={pair_digest(projective_deck_defect)}"
    )
    print(
        f"case={name};monic_Fprime_digest={pair_digest(monic_derivative)};"
        f"inverse_different_digest={pair_digest(inverse_different)};"
        "trace(Fprime_monic)=0;norm(Fprime_monic)=-delta/P2^2;"
        "inverse_different_unit=1;deck_parity=anti_invariant;"
        "inverse_different_is_not_a_supplied_Keller_cofactor=1"
    )
    print(
        f"case={name};velocity_digest={pair_digest(velocity)};"
        f"trace_velocity_digest={element_digest(trace_derivative)};"
        "sigma(velocity)=trace_derivative-velocity;"
        "first_normal_quadratic_deck_compatible=1"
    )
    return (
        pair_digest(Qelim),
        packet_digest((U, W)),
        pair_digest(projective_deck_defect),
        pair_digest(inverse_different),
        pair_digest(velocity),
        element_digest(trace_derivative),
    )


def source_audit() -> None:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imported_roots.update(
        node.module.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    require(
        not imported_roots.intersection({"importlib", "runpy", "subprocess"}),
        "no imported execution path",
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals},"
        "imported_execution_paths=0)"
    )


def main() -> None:
    print("finite-exact exceptional quadratic deck cofactor/integrability scout")
    for relative_path, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / relative_path) == expected_hash,
            f"dependency drift: {relative_path}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    packets = tuple(audit_case(response_case(name)) for name in ("4111", "3211"))
    require(packets == EXPECTED_CASE_PACKETS, "case digest packet drift")
    print(f"case_digest_packets={packets}")
    print(
        "consequence=B_restores_tautological_relative_root_label_and_etale_inverse_different;"
        "elimination_cofactor_pair_separates_conjugate_branches_but_does_not_descend;"
        "critical_gradient_pullback_kills_every_Jac(P,Q)_before_divergence_class_gate"
    )
    print(
        "scope=FINITE-EXACT_FIXED_SLICE_two_accessory_fields;"
        "degree_2_over_A_degree_72_over_K_per_MISTAKE_362;"
        "not_Keller_cofactor_not_mate_not_inverse_map_not_JC2_not_DC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
