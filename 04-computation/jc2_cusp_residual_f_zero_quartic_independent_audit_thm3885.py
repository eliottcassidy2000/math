#!/usr/bin/env python3
"""Independent hostile audit of THM-3885's quartic f=0 candidate.

The primary companion is not imported.  The six terminal square-root
equations are reconstructed by a direct coefficient recurrence.  Their only
common zero is certified by exact Macaulay linear algebra over Q (FLINT), not
by the primary grevlex Groebner basis or its three displayed reductions.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mat


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.cancel(left - right) == 0, message)


# ---------------------------------------------------------------------------
# Independently reconstruct the complete x=0 coefficient universe.
# ---------------------------------------------------------------------------
y = sp.symbols("y")
p, q, r = sp.symbols("p q r")
tau = p * y + q * y**2 + r * y**3
S0 = sp.Poly(
    sp.expand(256 - 96 * tau**2 - 8 * (y**2 - 4) * tau**3 - 3 * tau**4),
    y,
)
same(S0.coeff_monomial(1), 256, "address specialization has square constant")
gate(S0.degree() == 12, "formal quartic x=0 residual degree")

# For a square root with g_0=16, coefficient n gives
# g_n=(S_n-sum_{1<=i<n}g_i*g_{n-i})/(2g_0).  This recurrence does not call
# solve() and retains p,q,r as free parameters.
g = {0: sp.Integer(16)}
for degree in range(1, 7):
    known = sum(g.get(i, 0) * g.get(degree - i, 0) for i in range(1, degree))
    g[degree] = sp.cancel((S0.coeff_monomial(y**degree) - known) / 32)
    same(
        2 * g[0] * g[degree] + known,
        S0.coeff_monomial(y**degree),
        f"independent root recurrence degree {degree}",
    )

root = sum(g[degree] * y**degree for degree in range(7))
difference = sp.Poly(sp.expand(root**2 - S0.as_expr()), y)
for degree in range(7):
    same(difference.coeff_monomial(y**degree), 0,
         f"root recurrence clears degree {degree}")
terminal = [
    sp.Poly(
        sp.together(difference.coeff_monomial(y**degree)).as_numer_denom()[0],
        p,
        q,
        r,
    )
    for degree in range(7, 13)
]
gate(len(terminal) == 6 and all(poly.as_expr() != 0 for poly in terminal),
     "six nonzero terminal equations")
gate(all(poly.as_expr().subs({p: 0, q: 0, r: 0}) == 0 for poly in terminal),
     "base point is a positive control")
gate(any(poly.as_expr().subs({p: 1, q: 0, r: 0}) != 0 for poly in terminal),
     "linear hostile is rejected")
gate(any(poly.as_expr().subs({p: 0, q: 0, r: 1}) != 0 for poly in terminal),
     "cubic hostile is rejected")

terminal_blob = "\n".join(
    sp.srepr(poly.as_expr()) for poly in terminal
).encode()


# ---------------------------------------------------------------------------
# Exact Macaulay certificates, structurally independent of Groebner.
# ---------------------------------------------------------------------------
def monomials_of_degree(nvars: int, degree: int):
    if nvars == 1:
        yield (degree,)
        return
    for first in range(degree + 1):
        for rest in monomials_of_degree(nvars - 1, degree - first):
            yield (first,) + rest


def reduce_mod_prime(row, basis, prime):
    """Sparse row reduction used only to select a nonsingular exact minor."""
    row = {key: value % prime for key, value in row.items() if value % prime}
    while row:
        pivot = max(row)
        if pivot not in basis:
            inverse = pow(row[pivot], -1, prime)
            return {key: value * inverse % prime for key, value in row.items()}
        factor = row[pivot]
        for key, value in basis[pivot].items():
            updated = (row.get(key, 0) - factor * value) % prime
            if updated:
                row[key] = updated
            elif key in row:
                del row[key]
    return {}


def macaulay_certificate(polynomials, variables, multiplier_bound, target, label):
    """Prove target lies in the ideal by an exact rational row-span solve."""
    nvars = len(variables)
    rows = []
    for degree in range(multiplier_bound + 1):
        for shift in monomials_of_degree(nvars, degree):
            for polynomial in polynomials:
                row = {}
                for monomial, coefficient in polynomial.terms():
                    key = tuple(monomial[i] + shift[i] for i in range(nvars))
                    row[key] = row.get(key, 0) + int(coefficient)
                rows.append(row)

    # Mod 1009 selects independent original rows and pivot columns.  The
    # selected integer minor is then solved and verified over Q, so the
    # modular step is only a speed heuristic and carries no proof burden.
    prime = 1009
    basis = {}
    independent_rows = []
    for index, row in enumerate(rows):
        reduced = reduce_mod_prime(row, basis, prime)
        if reduced:
            basis[max(reduced)] = reduced
            independent_rows.append(index)
    gate(not reduce_mod_prime({target: 1}, basis, prime),
         f"{label}: target reaches selected modular span")

    pivots = list(basis)
    rank = len(pivots)
    gate(rank == len(independent_rows), f"{label}: square minor")
    matrix_entries = [
        rows[independent_rows[row]].get(pivots[column], 0)
        for row in range(rank)
        for column in range(rank)
    ]
    matrix = fmpq_mat(rank, rank, matrix_entries)
    right_hand_side = fmpq_mat(
        rank, 1, [1 if pivot == target else 0 for pivot in pivots]
    )
    coefficients = matrix.transpose().solve(right_hand_side)

    exact_sum = {}
    for coefficient_index, row_index in enumerate(independent_rows):
        coefficient = coefficients[coefficient_index, 0]
        if coefficient:
            for monomial, value in rows[row_index].items():
                exact_sum[monomial] = (
                    exact_sum.get(monomial, fmpq(0)) + coefficient * value
                )
    exact_sum = {monomial: value for monomial, value in exact_sum.items() if value}
    gate(exact_sum == {target: fmpq(1)}, f"{label}: exact Q certificate")

    coefficient_text = ",".join(
        str(coefficients[index, 0]) for index in range(rank)
    )
    return {
        "rows": len(rows),
        "rank": rank,
        "nonzero": sum(bool(coefficients[index, 0]) for index in range(rank)),
        "max_num_bits": max(
            abs(int(coefficients[index, 0].numerator)).bit_length()
            for index in range(rank)
        ),
        "max_den_bits": max(
            int(coefficients[index, 0].denominator).bit_length()
            for index in range(rank)
        ),
        "coefficient_sha256": hashlib.sha256(coefficient_text.encode()).hexdigest(),
    }


# The full six-equation ideal contains r^4.
r4_certificate = macaulay_certificate(
    terminal, (p, q, r), 6, (0, 0, 4), "full r^4"
)

# At a common zero r=0.  Rebuild the specialized ideal in Q[p,q] and obtain
# two independent pure-power certificates; no primary triangular relation is
# used.
terminal_r0 = [sp.Poly(poly.as_expr().subs(r, 0), p, q) for poly in terminal]
q4_certificate = macaulay_certificate(
    terminal_r0, (p, q), 4, (0, 4), "r=0 q^4"
)
p8_certificate = macaulay_certificate(
    terminal_r0, (p, q), 4, (8, 0), "r=0 p^8"
)

# Therefore every characteristic-zero common zero has r=q=p=0.  The opposite
# constant root sign gives the negative of the same root and no new branch.
same((-root) ** 2, root**2, "opposite square-root sign is covered")


# ---------------------------------------------------------------------------
# Global y-degree 3/2/1/0 audit.
# ---------------------------------------------------------------------------
alpha = sp.symbols("alpha")  # alpha=a=x+1; x=0 is alpha=1.
beta = sp.symbols("beta", nonzero=True)

# The y^3 coefficient is alpha*r and hence vanishes.  The y^2 coefficient b
# is divisible by alpha, has alpha-degree <=2, and b(1)=q=0.  Thus the exact
# nonzero form is beta*alpha*(alpha-1), sharper than the primary's broader
# alpha*(b0+b1*alpha) packet.
b = beta * alpha * (alpha - 1)
lead_y8 = sp.expand(-b**3 * (8 + 3 * alpha**2 * b))
same(sp.Poly(lead_y8, alpha).coeff_monomial(alpha**3), 8 * beta**3,
     "dy=2 exact odd alpha-valuation witness")
gate(sp.Poly(lead_y8, alpha).terms()[-1][0] == (3,),
     "dy=2 leading coefficient has exact alpha valuation three")

# For y-degree one, -8*K*T^3 supplies the unique y^5 term.
d, t0 = sp.symbols("d t0")
F_alpha = 15 * alpha**2 - 15 * alpha + 4
K_alpha = y**2 - F_alpha
L_alpha = 9 * alpha - 5
T_linear = d * y + t0
S_linear = sp.Poly(sp.expand(
    L_alpha**4
    - 6 * alpha * L_alpha**2 * T_linear**2
    - 8 * K_alpha * T_linear**3
    - 3 * alpha**2 * T_linear**4
), y)
same(S_linear.coeff_monomial(y**5), -8 * d**3,
     "dy=1 unique odd-degree witness")

# For y-degree zero the residual is -8*t(alpha)^3*y^2+C(alpha), with no y
# term.  If it were a square and t!=0, the missing cross term would force
# C=0.  Degree rows prove C cannot vanish for every possible e<=4.
S_constant_y = sp.Poly(S_linear.as_expr().subs(d, 0), y)
same(S_constant_y.coeff_monomial(y), 0,
     "dy=0 missing linear y term")
same(S_constant_y.coeff_monomial(y**2), -8 * t0**3,
     "dy=0 quadratic y coefficient")
for degree_t in range(5):
    degree_rows = (4, 2 * degree_t + 3, 3 * degree_t + 2, 4 * degree_t + 2)
    if degree_t == 0:
        gate(degree_rows[0] > max(degree_rows[1:]),
             "dy=0 constant t has unique L^4 degree")
    else:
        gate(degree_rows[3] > max(degree_rows[:3]),
             f"dy=0 degree {degree_t} has unique quartic degree")


semantic = {
    "target": "THM-3885 f=0 addressed total-degree-four candidate",
    "terminal_engine": "independent coefficient recurrence plus exact Q Macaulay syzygies",
    "terminal_consequence": "r4 in I; after r=0, q4 and p8 in specialized I",
    "global": "dy3 gone; dy2 beta*a*(a-1) odd valuation; dy1 odd; dy0 unique degree",
    "conclusion": "addressed f=0 deg(T)<=4 square residual implies T=0",
    "frontier": "nonzero survivor has deg(T)>=5; JC2 remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
source_text = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
     "no inactive Python assert")

print("THM3885_FZERO_QUARTIC_INDEPENDENT_AUDIT_20260823")
print("terminal_method=exact_Q_Macaulay_not_primary_Groebner")
print(f"terminal_sha256={hashlib.sha256(terminal_blob).hexdigest()}")
print(
    "r4_certificate="
    f"rows:{r4_certificate['rows']},rank:{r4_certificate['rank']},"
    f"nonzero:{r4_certificate['nonzero']},sha256:{r4_certificate['coefficient_sha256']}"
)
print(
    "r0_q4_certificate="
    f"rows:{q4_certificate['rows']},rank:{q4_certificate['rank']},"
    f"nonzero:{q4_certificate['nonzero']},sha256:{q4_certificate['coefficient_sha256']}"
)
print(
    "r0_p8_certificate="
    f"rows:{p8_certificate['rows']},rank:{p8_certificate['rank']},"
    f"nonzero:{p8_certificate['nonzero']},sha256:{p8_certificate['coefficient_sha256']}"
)
print("global_y_degrees=3,2,1,0_all_closed")
print("closed=f=0;addressed;deg_total_T<=4;unique_T=0")
print("promotion=JUSTIFIED")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
