#!/usr/bin/env python3
"""Exact referee for every degree-22 support-three coefficient stratum.

The already audited BCD companion contains the universal normalized fluxes,
the fixed quintic, and two independent finite-field-free Hensel engines (one
root and one unordered root pair).  This script pins that dependency by hash,
reconstructs the five-coefficient eliminant before choosing a scale, and runs
the two terminal-rank tests on the eight triples not already closed by
THM-2617 (BDW) or THM-2636 (BCD).  It also derives the two uniform y=0
boundary eliminants, for B nonzero and B zero respectively.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py"
BASE_SHA256 = "0866a29f665aedc6d2c226f35943852e56907ff821e705a0dbca2651e71fa15c"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited THM-2636 dependency changed",
)

spec = importlib.util.spec_from_file_location("thm2636_base", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load THM-2636")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)

t, v, zeta = base.t, base.v, base.zeta
c, d = base.c, base.d
B0, C0, D0, E0, W0 = s.symbols("B0 C0 D0 E0 W0")

# These are the normalized fluxes before choosing which nonzero coefficient
# supplies rho.  Substitution B0=t^2, C0=c*t^3, ..., W0=w*t^6 recovers the
# universal eliminant used by THM-2636 exactly.
F1 = (
    819896 * B0 * zeta
    - 1449459 * v * zeta
    + 83853 * zeta
    - 2981440 * B0 * v
    + 24640 * B0
    + 9370240 * C0 * v
    - 232320 * C0
    + 2044416 * D0
    - 14992384 * E0
    + 3689532 * v**2
    - 101640 * v
    + 252
)
F2 = (
    15944049 * zeta**2
    + 65591680 * B0 * zeta
    - 206145280 * C0 * zeta
    - 162339408 * v * zeta
    + 2236080 * zeta
    + 1443016960 * B0 * v**2
    - 71554560 * B0 * v
    + 98560 * B0
    + 449771520 * C0 * v
    - 1239040 * C0
    - 1978994688 * D0 * v
    + 16355328 * D0
    - 239878144 * E0
    - 1319329792 * W0
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)

general_raw = s.resultant(F1, F2, zeta)
general_content, general_primitive = s.Poly(
    general_raw, B0, C0, D0, E0, W0, v, domain=s.QQ
).primitive()
R_general = general_primitive.as_expr()
require(general_content == 28344976, "general resultant content changed")
require(len(general_primitive.terms()) == 60, "general eliminant support changed")
scaled_universal = s.expand(
    R_general.subs(
        {
            B0: t**2,
            C0: base.c * t**3,
            D0: base.d * t**4,
            E0: base.e * t**5,
            W0: base.w * t**6,
        }
    )
)
require(
    s.expand(scaled_universal - base.R_universal) == 0,
    "pre-scale eliminant does not recover THM-2636",
)


class TerminalHenselSystem(base.HenselSystem):
    """THM-2636 Hensel engine with a variable terminal order/support."""

    def terminal_certificate(
        self,
        r_polys: list[base.VP],
        terminal: int,
        expected_support: set[tuple[int, int]],
    ) -> dict[str, object]:
        require(
            self.v_equal(self.v_mul(self.q0, self.s0), r_polys[0]),
            f"{self.name}: fixed factor/cofactor product failed",
        )
        qs: list[base.VP] = [self.q0]
        ss: list[base.VP] = [self.s0]
        controls = 0
        zero = self.v_from_field([0])
        for n in range(1, terminal + 1):
            rn = r_polys[n] if n < len(r_polys) else zero
            convolution = zero
            for index in range(1, n):
                convolution = self.v_add(
                    convolution, self.v_mul(qs[index], ss[n - index])
                )
            residual = self.v_add(rn, self.v_neg(convolution))
            qn = self.v_remainder(
                self.v_mul(residual, self.inv_s0_mod_q0), self.q0
            )
            numerator = self.v_add(residual, self.v_neg(self.v_mul(qn, self.s0)))
            sn, remainder = self.v_divmod_monic(numerator, self.q0)
            require(
                all(not coefficient for coefficient in remainder),
                f"{self.name}: Hensel division failed at order {n}",
            )
            qs.append(qn)
            ss.append(sn)
            reconstructed = zero
            for index in range(n + 1):
                reconstructed = self.v_add(
                    reconstructed, self.v_mul(qs[index], ss[n - index])
                )
            require(
                self.v_equal(reconstructed, rn),
                f"{self.name}: product control failed at order {n}",
            )
            controls += 1

        terminal_equations = qs[terminal] + ss[terminal]
        require(
            len(terminal_equations) == 5,
            f"{self.name}: terminal equation count changed",
        )
        actual_support = set().union(
            *(set(coefficient) for coefficient in terminal_equations)
        )
        require(
            actual_support == expected_support,
            f"{self.name}: terminal support {actual_support} != {expected_support}",
        )
        monomial_basis = sorted(expected_support)
        matrix = [
            [coefficient.get(monomial, s.Integer(0)) for monomial in monomial_basis]
            for coefficient in terminal_equations
        ]
        width = len(monomial_basis)
        minors: list[s.Expr] = []
        for rows in combinations(range(5), width):
            minors.append(self.determinant([matrix[index] for index in rows]))
        nonzero = [minor for minor in minors if not self.field.zero(minor)]
        require(nonzero, f"{self.name}: terminal matrix lost full column rank")
        require(
            all(self.field.coprime_numerator(minor) for minor in nonzero),
            f"{self.name}: terminal minor is not uniform on conjugates",
        )
        first = self.field.numerator(nonzero[0])
        return {
            "equations": len(terminal_equations),
            "rank": width,
            "nonzero_minors": len(nonzero),
            "minor_degree": first.degree(),
            "minor_terms": len(first.terms()),
            "controls": controls,
            "support": tuple(monomial_basis),
        }


def scaled_eliminant(
    substitutions: dict[s.Symbol, s.Expr | int],
) -> tuple[s.Expr, s.Poly, s.Rational, s.Expr, list[s.Expr]]:
    expression = s.expand(R_general.subs(substitutions))
    raw_poly = s.Poly(expression, t, v, c, d, domain=s.QQ)
    content, primitive = raw_poly.primitive()
    leading = s.Poly(primitive.as_expr(), v).coeff_monomial(v**5)
    require(leading != 0 and not leading.free_symbols, "v-leading term is not constant")
    monic = s.expand(primitive.as_expr() / leading)
    require(s.Poly(monic, v).LC() == 1, "monic normalization failed")
    degree_t = s.Poly(monic, t).degree()
    coefficients = [
        s.Poly(monic, t).coeff_monomial(t**n) for n in range(degree_t + 1)
    ]
    require(
        s.expand(coefficients[0] - base.P5_expr) == 0,
        "fixed quintic changed under specialization",
    )
    return monic, primitive, content, leading, coefficients


TRIPLES = (
    # name, active parameter labels, scale substitution, terms, t-degree,
    # terminal order, terminal monomial support in the generic (c,d) chart.
    (
        "BCE",
        ("C", "E"),
        {B0: t**2, C0: c * t**3, D0: 0, E0: d * t**5, W0: 0},
        41,
        10,
        11,
        {(0, 1), (1, 0), (2, 1), (3, 0)},
    ),
    (
        "BCW",
        ("C", "W"),
        {B0: t**2, C0: c * t**3, D0: 0, E0: 0, W0: d * t**6},
        38,
        10,
        11,
        {(1, 0), (1, 1), (3, 0)},
    ),
    (
        "BDE",
        ("D", "E"),
        {B0: t**2, C0: 0, D0: c * t**4, E0: d * t**5, W0: 0},
        36,
        10,
        11,
        {(0, 1), (1, 1)},
    ),
    (
        "BEW",
        ("E", "W"),
        {B0: t**2, C0: 0, D0: 0, E0: c * t**5, W0: d * t**6},
        31,
        10,
        11,
        {(1, 0), (1, 1)},
    ),
    (
        "CDE",
        ("D", "E"),
        {B0: 0, C0: t**3, D0: c * t**4, E0: d * t**5, W0: 0},
        25,
        10,
        11,
        {(0, 1), (2, 0)},
    ),
    (
        "CDW",
        ("D", "W"),
        {B0: 0, C0: t**3, D0: c * t**4, E0: 0, W0: d * t**6},
        22,
        8,
        9,
        {(0, 0), (0, 1)},
    ),
    (
        "CEW",
        ("E", "W"),
        {B0: 0, C0: t**3, D0: 0, E0: c * t**5, W0: d * t**6},
        21,
        10,
        11,
        {(1, 0), (1, 1)},
    ),
    (
        "DEW",
        ("E", "W"),
        {B0: 0, C0: 0, D0: t**4, E0: c * t**5, W0: d * t**6},
        19,
        10,
        11,
        {(1, 1)},
    ),
)

certificates: list[dict[str, object]] = []
for name, labels, substitution, terms, degree_t, terminal, support in TRIPLES:
    monic, primitive, content, leading, coefficients = scaled_eliminant(substitution)
    require(len(primitive.terms()) == terms, f"{name}: term count changed")
    require(s.Poly(monic, t).degree() == degree_t, f"{name}: t-degree changed")
    require(s.Poly(monic, v).degree() == 5, f"{name}: v-degree changed")

    root_system = TerminalHenselSystem(
        f"{name}_root", base.root_field, base.root_q0, base.root_s0
    )
    pair_system = TerminalHenselSystem(
        f"{name}_pair", base.pair_field, base.pair_q0, base.pair_s0
    )
    root_polys = [root_system.expression_to_vpoly(value) for value in coefficients]
    pair_polys = [pair_system.expression_to_vpoly(value) for value in coefficients]
    root_result = root_system.terminal_certificate(root_polys, terminal, support)
    pair_result = pair_system.terminal_certificate(pair_polys, terminal, support)
    require(
        root_result["rank"] == len(support) == pair_result["rank"],
        f"{name}: root/pair rank disagreement",
    )
    certificates.append(
        {
            "name": name,
            "labels": labels,
            "terms": terms,
            "degree_t": degree_t,
            "terminal": terminal,
            "support": tuple(sorted(support)),
            "content": content,
            "leading": leading,
            "root": root_result,
            "pair": pair_result,
        }
    )

# The boundary y=0 is independent of the scale chart.  First normalize B=1.
u, Z, Cb, Db, Eb, Wb = s.symbols("u Z Cb Db Eb Wb")
N1_B = 1331 * (616 - 1089 * u) * Z + 9370240 * Cb * u - 14992384 * Eb
N2_B = (
    15944049 * Z**2
    - 206145280 * Cb * Z
    + 1443016960 * u**2
    - 1978994688 * Db * u
    - 1319329792 * Wb
    - 1190488992 * u**3
)
B_boundary_content, B_boundary_primitive = s.Poly(
    s.resultant(N1_B, N2_B, Z), u, Cb, Db, Eb, Wb, domain=s.QQ
).primitive()
B_boundary_expected = (
    -2264031 * u**5
    + 5305608 * u**4
    - (3763584 * Db + 3829056) * u**3
    + (-1267200 * Cb**2 + 4257792 * Db - 2509056 * Wb + 878080) * u**2
    + (1433600 * Cb**2 - 1204224 * Db + 2838528 * Wb) * u
    - 2293760 * Cb * Eb
    + 3244032 * Eb**2
    - 802816 * Wb
)
require(
    B_boundary_content == 1104726788605792
    and s.expand(B_boundary_primitive.as_expr() - B_boundary_expected) == 0,
    "B-nonzero y=0 eliminant content or primitive part changed",
)
require(
    s.Poly(B_boundary_expected, u).degree() == 5
    and s.Poly(B_boundary_expected, u).LC() == -2264031,
    "B-nonzero boundary lost constant quintic leading term",
)

# On B=0 the open first-flux chart is -1089u != 0.  The same elimination
# gives a five-term quintic, valid uniformly for all four B-free triples.
N1_0 = -1449459 * u * Z + 9370240 * Cb * u - 14992384 * Eb
N2_0 = (
    15944049 * Z**2
    - 206145280 * Cb * Z
    - 1978994688 * Db * u
    - 1319329792 * Wb
    - 1190488992 * u**3
)
zero_boundary_content, zero_boundary_primitive = s.Poly(
    s.resultant(N1_0, N2_0, Z), u, Cb, Db, Eb, Wb, domain=s.QQ
).primitive()
zero_boundary_expected = (
    -22869 * u**5
    - 38016 * Db * u**3
    - (12800 * Cb**2 + 25344 * Wb) * u**2
    + 32768 * Eb**2
)
require(
    s.expand(zero_boundary_primitive.as_expr() - zero_boundary_expected) == 0,
    "B-zero y=0 eliminant changed",
)
require(
    zero_boundary_content == 109367952071973408
    and s.Poly(zero_boundary_expected, u).LC() == -22869,
    "B-zero boundary content/leading term changed",
)

print("degree-22 complete support-three weighted Hensel closure")
print(f"base_thm2636_sha256={BASE_SHA256}")
print(f"general_resultant_content={general_content}")
print(f"general_terms={len(general_primitive.terms())}")
print("fixed_quintic_degree=5")
print("root_field_degree=5")
print("pair_field_degree=10")
for certificate in certificates:
    name = certificate["name"]
    labels = ",".join(certificate["labels"])
    support = ";".join(f"{i},{j}" for i, j in certificate["support"])
    print(
        f"{name}:params={labels}:terms={certificate['terms']}:"
        f"degrees={certificate['degree_t']},5:terminal={certificate['terminal']}:"
        f"support={support}:content={certificate['content']}:"
        f"v5={certificate['leading']}"
    )
    for field_name in ("root", "pair"):
        result = certificate[field_name]
        print(
            f"{name}_{field_name}:equations={result['equations']}:"
            f"rank={result['rank']}:nonzero_minors={result['nonzero_minors']}:"
            f"first_minor_degree={result['minor_degree']}:"
            f"first_minor_terms={result['minor_terms']}:"
            f"product_controls={result['controls']}:uniform=True"
        )
print(f"B_nonzero_y0_content={B_boundary_content}")
print("B_nonzero_y0_degree=5:v5=-2264031")
print(f"B_zero_y0_content={zero_boundary_content}")
print("B_zero_y0_degree=5:v5=-22869:open_chart_u_nonzero=True")
print("inherited_closed_controls=BCD,BDW")
print("newly_closed_triples=BCE,BCW,BDE,BEW,CDE,CDW,CEW,DEW")
print("all_ten_support_three_strata=EMPTY")
print("ALL CHECKS PASSED")
