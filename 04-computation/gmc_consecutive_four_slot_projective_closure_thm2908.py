#!/usr/bin/env python3
"""Exact companion for THM-2908.

The finite affine chart is reconstructed from the hash-pinned THM-2879
companion.  An exact FLINT resultant and factorization certify that the
cubic branch eliminant and selector-cleared endpoint holonomy have no
common root at any integer depth.  A separate SymPy calculation closes
the projective line at infinity.

This script requires ``python-flint``.  Every truth-bearing gate uses
``require`` rather than ``assert``.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path

import sympy as sp

try:
    from flint import fmpz_mpoly_ctx
except ModuleNotFoundError as error:
    raise RuntimeError(
        "THM-2908 exact replay requires python-flint"
    ) from error


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canonical_digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(exponent) for exponent in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


def newton_values(polynomial: sp.Poly, base: int) -> list[sp.Integer]:
    values = [
        sp.Integer(polynomial.eval(base + offset))
        for offset in range(polynomial.degree() + 1)
    ]
    answer: list[sp.Integer] = []
    for _ in range(polynomial.degree() + 1):
        answer.append(values[0])
        values = [
            values[index + 1] - values[index]
            for index in range(len(values) - 1)
        ]
    return answer


def has_no_residue_root(polynomial: sp.Poly, prime: int) -> bool:
    coefficients = [
        int(coefficient) % prime
        for coefficient in polynomial.all_coeffs()
    ]
    for residue in range(prime):
        value = 0
        for coefficient in coefficients:
            value = (value * residue + coefficient) % prime
        if value == 0:
            return False
    return True


dependency = Path(__file__).with_name(
    "gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py"
)
dependency_bytes = dependency.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(dependency_bytes).hexdigest()
    == "44012d84c88a22f246ef350f7f9a364116ac1fc839347361dee64c0a9c4a6e27",
    "THM-2879 exact dependency hash changed",
)
spec = importlib.util.spec_from_file_location("thm2879_exact", dependency)
require(spec is not None and spec.loader is not None, "cannot load THM-2879")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
n, x, y = source.n, source.x, source.y

# Finite high chart: cubic eliminant and primitive linear selector.
invariant_pairs = tuple(
    sp.together(value).as_numer_denom()
    for value in (source.cubic_one, source.cubic_two)
)
invariants = tuple(numerator for numerator, _ in invariant_pairs)
expected_cubic_denominators = (
    72
    * (2 * n + 1) ** 2
    * (2 * n + 3) ** 2
    * (3 * n + 1)
    * (3 * n + 2)
    * (3 * n + 4)
    * (3 * n + 5),
    36
    * (2 * n + 1)
    * (2 * n + 3) ** 2
    * (3 * n + 1)
    * (3 * n + 2)
    * (3 * n + 4)
    * (3 * n + 5),
)
require(
    all(
        sp.cancel(denominator / expected) == 1
        for (_, denominator), expected
        in zip(invariant_pairs, expected_cubic_denominators)
    ),
    "cleared cubic denominator changed",
)
subresultants = sp.subresultants(*invariants, x)
require(
    [sp.degree(value, x) for value in subresultants] == [4, 3, 2, 1, 0],
    "cubic subresultant profile changed",
)
eliminant = next(
    factor
    for factor, _ in sp.factor_list(subresultants[-1])[1]
    if sp.degree(factor, y) == 15
)
if sp.Poly(eliminant, y).LC().subs(n, 1) < 0:
    eliminant = -eliminant
eliminant_poly = sp.Poly(eliminant, n, y, domain=sp.QQ)

linear = sp.Poly(subresultants[-2], x)
linear_content = sp.gcd(linear.nth(1), linear.nth(0))
selector_a = sp.cancel(linear.nth(1) / linear_content)
selector_n = sp.cancel(-linear.nth(0) / linear_content)
if selector_a.subs({n: 1, y: 1}) < 0:
    selector_a, selector_n = -selector_a, -selector_n
selector_a_poly = sp.Poly(selector_a, n, y, domain=sp.QQ)
selector_n_poly = sp.Poly(selector_n, n, y, domain=sp.QQ)

content_quadratic = (
    (4 * n + 6) * y**2 + (4 * n + 6) * y + n + 2
)
require(
    sp.expand(
        content_quadratic
        - ((4 * n + 6) * (y + sp.Rational(1, 2)) ** 2 + sp.Rational(1, 2))
    )
    == 0
    and sp.Poly(content_quadratic, y).LC() == 4 * n + 6,
    "positive completed-square content identity changed",
)
expected_content = (
    864
    * (n + 1)
    * (n + 2) ** 3
    * (2 * n + 1) ** 4
    * (3 * n + 1) ** 3
    * (3 * n + 2) ** 2
    * content_quadratic
)
content_ratio = sp.cancel(linear_content / expected_content)
require(
    content_ratio.is_Rational and content_ratio != 0,
    "linear-subresultant content changed",
)
require(
    sp.discriminant(content_quadratic, y) == -4 * (2 * n + 3),
    "content discriminant changed",
)

# Homogeneously clear the degree-five endpoint numerator at A X=N.
endpoint_numerator = sp.together(
    source.endpoint_holonomy
).as_numer_denom()[0]
endpoint_denominator = sp.together(
    source.endpoint_holonomy
).as_numer_denom()[1]
expected_endpoint_denominator = (
    1024
    * (n + 3)
    * (2 * n + 1) ** 2
    * (2 * n + 3) ** 4
    * (4 * n + 1)
    * (4 * n + 3)
    * (4 * n + 5)
    * (4 * n + 7)
)
require(
    sp.cancel(endpoint_denominator / expected_endpoint_denominator) == 1,
    "cleared endpoint denominator changed",
)
endpoint_poly = sp.Poly(endpoint_numerator, x)
endpoint_degree = endpoint_poly.degree()
require(endpoint_degree == 5, "endpoint degree changed")
powers_a = [sp.Poly(1, n, y, domain=sp.QQ)]
powers_n = [sp.Poly(1, n, y, domain=sp.QQ)]
for _ in range(endpoint_degree):
    powers_a.append(powers_a[-1] * selector_a_poly)
    powers_n.append(powers_n[-1] * selector_n_poly)
selector_cleared = sum(
    (
        sp.Poly(endpoint_poly.nth(index), n, y, domain=sp.QQ)
        * powers_n[index]
        * powers_a[endpoint_degree - index]
        for index in range(endpoint_degree + 1)
    ),
    sp.Poly(0, n, y, domain=sp.QQ),
)
require(
    (
        eliminant_poly.degree(n),
        eliminant_poly.degree(y),
        len(eliminant_poly.terms()),
        selector_a_poly.degree(n),
        selector_a_poly.degree(y),
        len(selector_a_poly.terms()),
        selector_n_poly.degree(n),
        selector_n_poly.degree(y),
        len(selector_n_poly.terms()),
        selector_cleared.degree(n),
        selector_cleared.degree(y),
        len(selector_cleared.terms()),
    )
    == (21, 15, 352, 22, 10, 251, 22, 10, 253, 121, 55, 6820),
    "finite-chart source profile changed",
)

# Verify the full cubic resultant factorization used to infer that every
# real cubic-divisible plane lies on P=0.  In particular, the n-only
# degree-ten cofactor is positive and G has no real zero.
cubic_resultant = sp.Poly(subresultants[-1], n, y, domain=sp.QQ)
cubic_unit, cubic_factors = sp.factor_list(
    cubic_resultant.as_expr(),
    n,
    y,
)
require(
    sorted(
        (
            sp.degree(factor, n),
            sp.degree(factor, y),
            exponent,
        )
        for factor, exponent in cubic_factors
    )
    == [
        (1, 0, 2),
        (1, 0, 2),
        (1, 0, 3),
        (1, 0, 4),
        (1, 0, 5),
        (1, 2, 2),
        (10, 0, 1),
        (21, 15, 1),
    ],
    "full cubic resultant factor profile changed",
)
cubic_g_factor = next(
    factor
    for factor, _ in cubic_factors
    if sp.degree(factor, y) == 2
)
require(
    sp.cancel(cubic_g_factor / content_quadratic).is_Rational
    and sp.cancel(cubic_g_factor / content_quadratic) != 0,
    "cubic resultant lost the exact G factor",
)
cubic_nonlinear_cofactor = sp.Poly(
    next(
        factor
        for factor, _ in cubic_factors
        if sp.degree(factor, n) == 10 and sp.degree(factor, y) == 0
    ),
    n,
)
require(
    all(
        coefficient > 0
        for coefficient in cubic_nonlinear_cofactor.all_coeffs()
    ),
    "degree-ten cubic cofactor lost coefficient positivity",
)
require(
    all(
        (
            sp.degree(factor, n) == 10
            or sp.degree(factor, y) > 0
            or (
                sp.Poly(factor, n).eval(0) > 0
                and sp.Poly(factor, n).LC() > 0
            )
        )
        for factor, _ in cubic_factors
    ),
    "an n-only cubic resultant factor is not positive at n>=0",
)

# Exact direct resultant over Z[n].  FLINT returns a certified
# factorization; reconstructing the product checks the API boundary.
context = fmpz_mpoly_ctx.get(["n", "y"])
_, eliminant_integer = eliminant_poly.clear_denoms(convert=True)
_, endpoint_integer = selector_cleared.clear_denoms(convert=True)
eliminant_flint = context.from_dict(
    {
        monomial: int(coefficient)
        for monomial, coefficient in eliminant_integer.as_dict().items()
    }
)
endpoint_flint = context.from_dict(
    {
        monomial: int(coefficient)
        for monomial, coefficient in endpoint_integer.as_dict().items()
    }
)
direct_resultant = eliminant_flint.resultant(endpoint_flint, 1)
resultant_unit, resultant_factors = direct_resultant.factor()
reconstructed_resultant = context.from_dict(
    {(0, 0): int(resultant_unit)}
)
for factor, exponent in resultant_factors:
    reconstructed_resultant *= factor**exponent
require(
    reconstructed_resultant == direct_resultant,
    "FLINT factor product does not reconstruct resultant",
)
require(
    (direct_resultant.degrees()[0], len(direct_resultant))
    == (2804, 2805),
    "direct resultant profile changed",
)

factor_records: list[
    tuple[int, int, str, str, int | None, sp.Poly]
] = []
for factor, exponent in resultant_factors:
    require(factor.degrees()[1] == 0, "factor retained y")
    factor_poly = sp.Poly.from_dict(
        {
            (monomial[0],): int(coefficient)
            for monomial, coefficient in factor.to_dict().items()
        },
        (n,),
        domain=sp.ZZ,
    )
    if factor_poly.LC() < 0:
        factor_poly = -factor_poly
    degree = factor_poly.degree()
    if degree == 1:
        require(
            factor_poly.eval(0) > 0 and factor_poly.LC() > 0,
            "linear resultant factor is not positive at n>=0",
        )
        certificate = "linear-positive"
        root_free_prime = None
    else:
        differences = newton_values(factor_poly, 1)
        if factor_poly.eval(0) > 0 and all(value > 0 for value in differences):
            certificate = "GN-positive"
            root_free_prime = None
        else:
            certificate = "root-free-mod-p"
            root_free_prime = {65: 43, 203: 83}.get(degree)
            require(
                root_free_prime is not None
                and has_no_residue_root(factor_poly, root_free_prime),
                f"fixed root-free-prime certificate failed for degree {degree}",
            )
    factor_records.append(
        (
            degree,
            exponent,
            certificate,
            canonical_digest(factor_poly),
            root_free_prime,
            factor_poly,
        )
    )

require(
    sorted((degree, exponent) for degree, exponent, *_ in factor_records)
    == [
        (1, 15),
        (1, 24),
        (1, 37),
        (1, 48),
        (1, 105),
        (1, 360),
        (5, 161),
        (10, 50),
        (19, 3),
        (65, 10),
        (203, 1),
    ],
    "direct resultant factor profile changed",
)
expected_digests = {
    5: "bb1a42b92837273e8d421ae22300fbf785edf0418b4837aec6807a9e8e388fce",
    10: "cd924e68200c7b69719b57aee6a712e6f74eba50fa1e2f30d858b687ed693da8",
    19: "c725b7cf35b3efcd4304ffe50a34e5f2d61069497a4caf8fad3958256621f8a5",
    65: "2182971d85ea1aed55dfc0b4211544bc76b502896ccb7fd003cd75172706cda6",
    203: "eac4e0d0fdb5a6a787b2614bd3704d8d777e7bf41fae9ddf9f5e1ce1cd89ab22",
}
require(
    {
        degree: digest
        for degree, _, _, digest, _, _ in factor_records
        if degree > 1
    }
    == expected_digests,
    "nonlinear resultant factor digest changed",
)
require(
    {
        degree: certificate
        for degree, _, certificate, _, _, _ in factor_records
        if degree > 1
    }
    == {
        5: "GN-positive",
        10: "GN-positive",
        19: "GN-positive",
        65: "root-free-mod-p",
        203: "root-free-mod-p",
    },
    "affine factor certificate classes changed",
)
require(
    {
        degree: root_free_prime
        for degree, _, _, _, root_free_prime, _ in factor_records
        if degree in {65, 203}
    }
    == {65: 43, 203: 83},
    "root-free-prime certificates changed",
)
require(
    all(
        certificate in {"linear-positive", "GN-positive"}
        or root_free_prime is not None
        for _, _, certificate, _, root_free_prime, _ in factor_records
    ),
    "an affine factor lacks an integer-nonvanishing certificate",
)

# Projective line at infinity.  Finite t represents
# span(d_n-t*d_(n+1), d_(n+2)).
t = sp.symbols("t")
u_infinity = {0: sp.Integer(1), 1: -t}
v_infinity = {2: sp.Integer(1)}
g0 = source.multilinear((u_infinity, u_infinity))
g1 = source.multilinear((u_infinity, v_infinity))
g2 = source.multilinear((v_infinity, v_infinity))
t0 = source.multilinear((u_infinity, u_infinity, u_infinity))
t1 = source.multilinear((u_infinity, u_infinity, v_infinity))
t2 = source.multilinear((u_infinity, v_infinity, v_infinity))
t3 = source.multilinear((v_infinity, v_infinity, v_infinity))
infinity_i1 = sp.factor(
    3 * t1 * g0 * g2 - t3 * g0**2 - 2 * t0 * g1 * g2
)
infinity_i2 = sp.factor(
    3 * t2 * g0 * g2 - 2 * t3 * g1 * g0 - t0 * g2**2
)
infinity_i1_num = sp.Poly(
    sp.together(infinity_i1).as_numer_denom()[0],
    t,
    n,
    domain=sp.QQ,
)
infinity_i2_num = sp.Poly(
    sp.together(infinity_i2).as_numer_denom()[0],
    t,
    n,
    domain=sp.QQ,
)
infinity_resultant = sp.Poly(
    sp.resultant(
        infinity_i1_num.as_expr(),
        infinity_i2_num.as_expr(),
        t,
    ),
    n,
)
infinity_unit, infinity_factors = sp.factor_list(
    infinity_resultant.as_expr(),
    n,
)
require(
    [
        (sp.degree(factor, n), exponent)
        for factor, exponent in infinity_factors
    ]
    == [(1, 2), (1, 3), (1, 3), (1, 3), (1, 4), (1, 10), (5, 1), (19, 1)],
    "infinity resultant factor profile changed",
)
require(
    all(
        sp.Poly(factor, n).eval(0) > 0
        and sp.Poly(factor, n).LC() > 0
        and (
            sp.degree(factor, n) == 1
            or all(
                coefficient > 0
                for coefficient in sp.Poly(factor, n).all_coeffs()
            )
        )
        for factor, _ in infinity_factors
    ),
    "infinity factor lost positivity at n>=0",
)

infinity_u = {1: sp.Integer(1)}
infinity_v = {2: sp.Integer(1)}
ig0 = source.multilinear((infinity_u, infinity_u))
ig1 = source.multilinear((infinity_u, infinity_v))
ig2 = source.multilinear((infinity_v, infinity_v))
it0 = source.multilinear((infinity_u, infinity_u, infinity_u))
it1 = source.multilinear((infinity_u, infinity_u, infinity_v))
it2 = source.multilinear((infinity_u, infinity_v, infinity_v))
it3 = source.multilinear((infinity_v, infinity_v, infinity_v))
t_infinity_i1 = sp.factor(
    sp.together(
        3 * it1 * ig0 * ig2 - it3 * ig0**2 - 2 * it0 * ig1 * ig2
    ).as_numer_denom()[0]
)
t_infinity_i2 = sp.factor(
    sp.together(
        3 * it2 * ig0 * ig2 - 2 * it3 * ig1 * ig0 - it0 * ig2**2
    ).as_numer_denom()[0]
)
require(
    t_infinity_i1 == -(n + 2) ** 2 * (28 * n**2 + 87 * n + 66),
    "t=infinity first cubic invariant changed",
)
require(
    t_infinity_i2 == -2 * (n + 2) * (4 * n + 5),
    "t=infinity second cubic invariant changed",
)

infinity_records = "\n".join(
    f"{monomial[0]}:{coefficient}"
    for monomial, coefficient in infinity_resultant.terms()
)
print("THM-2908 CONSECUTIVE FOUR-SLOT PROJECTIVE CLOSURE")
print(
    "finite_profiles="
    "P:21,15,352;A:22,10,251;N:22,10,253;F:121,55,6820"
)
print(
    "cubic_resultant_profile="
    + ",".join(
        f"{sp.degree(factor, n)},{sp.degree(factor, y)}^{exponent}"
        for factor, exponent in cubic_factors
    )
)
print(
    "cubic_degree10_all_coefficients_positive=true;"
    f"constant={cubic_nonlinear_cofactor.eval(0)}"
)
print(f"selector_content_ratio={content_ratio}")
print("selector_content_discriminant=-4*(2*n+3)")
print(
    "direct_resultant="
    f"degree:{direct_resultant.degrees()[0]};terms:{len(direct_resultant)}"
)
print(
    "direct_factor_profile="
    + ",".join(
        f"{degree}^{exponent}:{certificate}"
        + (
            f":p={root_free_prime}"
            if root_free_prime is not None
            else ""
        )
        for degree, exponent, certificate, _, root_free_prime, _
        in sorted(factor_records)
    )
)
print(
    "direct_nonlinear_digests="
    + ",".join(
        f"{degree}:{digest}"
        for degree, _, _, digest, _, _ in sorted(factor_records)
        if degree > 1
    )
)
print(
    "infinity_profiles="
    f"I1:{infinity_i1_num.degree(t)},{infinity_i1_num.degree(n)},"
    f"{len(infinity_i1_num.terms())};"
    f"I2:{infinity_i2_num.degree(t)},{infinity_i2_num.degree(n)},"
    f"{len(infinity_i2_num.terms())}"
)
print(
    f"infinity_resultant=degree:{infinity_resultant.degree()};"
    "factor_profile:"
    + ",".join(
        f"{sp.degree(factor, n)}^{exponent}"
        for factor, exponent in infinity_factors
    )
)
print(
    "infinity_resultant_digest="
    f"{sha256((infinity_records + chr(10)).encode()).hexdigest()}"
)
print(
    "t_infinity_invariants="
    "-(n+2)^2*(28*n^2+87*n+66);-2*(n+2)*(4*n+5)"
)
print("scope=consecutive first-window four-slot SFC only")
print("all_exact_checks=PASS")
