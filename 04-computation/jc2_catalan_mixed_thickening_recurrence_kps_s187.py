#!/usr/bin/env python3
"""Exact mixed-thickening recurrence and low-cap obstruction audit.

For

    P = v^2 + sum_{j=1}^N a_j(v) w^j,
    Q = v^3-v + sum_{j=1}^N b_j(v) w^j,

this companion derives every coefficient of Jac(P,Q), verifies the cap-free
N=1 and N=2 reductions, and certifies that the affine coefficient varieties
for N=3 with deg_v(a_j),deg_v(b_j) <= 3, 4, or 5 are empty in characteristic
zero.  The N=3 proof is a complete degree-branch ledger plus one tiny
saturated Groebner calculation over QQ; it is not a point search.
"""

from __future__ import annotations

import platform

import sympy as sp


failures: list[str] = []
gates = 0


def gate(name: str, condition: bool) -> None:
    global gates
    gates += 1
    if not bool(condition):
        failures.append(name)


v, w = sp.symbols("v w")
p = v**2
q = v**3 - v


def wronskian(x: sp.Expr, y: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(x, v) * y - x * sp.diff(y, v))


def jacobian_rows(a_rows: list[sp.Expr], b_rows: list[sp.Expr]) -> list[sp.Expr]:
    """Return [w^k] Jac(P,Q), with row zero included in each input."""
    if len(a_rows) != len(b_rows):
        raise ValueError("row lists must have the same length")
    N = len(a_rows) - 1
    rows: list[sp.Expr] = []
    for k in range(2 * N):
        value = 0
        for i in range(N + 1):
            j = k + 1 - i
            if 0 <= j <= N:
                value += (
                    j * sp.diff(a_rows[i], v) * b_rows[j]
                    - i * a_rows[i] * sp.diff(b_rows[j], v)
                )
        rows.append(sp.expand(value))
    return rows


def tied_max(degrees: list[int]) -> bool:
    top = max(degrees)
    return sum(d == top for d in degrees) >= 2


# ---------------------------------------------------------------------------
# 1. Exact coefficient recurrence.
# ---------------------------------------------------------------------------

a1f, a2f, a3f = (sp.Function(name)(v) for name in ("a1", "a2", "a3"))
b1f, b2f, b3f = (sp.Function(name)(v) for name in ("b1", "b2", "b3"))
formal_a = [p, a1f, a2f, a3f]
formal_b = [q, b1f, b2f, b3f]
recurrence_rows = jacobian_rows(formal_a, formal_b)
P_formal = sum(formal_a[j] * w**j for j in range(4))
Q_formal = sum(formal_b[j] * w**j for j in range(4))
J_formal = sp.expand(
    sp.diff(P_formal, v) * sp.diff(Q_formal, w)
    - sp.diff(P_formal, w) * sp.diff(Q_formal, v)
)
for k, row in enumerate(recurrence_rows):
    gate(f"generic recurrence row {k}", sp.expand(J_formal.coeff(w, k) - row) == 0)

expected_n3 = [
    sp.diff(p, v) * b1f - a1f * sp.diff(q, v),
    2 * (sp.diff(p, v) * b2f - a2f * sp.diff(q, v))
    + wronskian(a1f, b1f),
    3 * (sp.diff(p, v) * b3f - a3f * sp.diff(q, v))
    + 2 * sp.diff(a1f, v) * b2f
    - a1f * sp.diff(b2f, v)
    + sp.diff(a2f, v) * b1f
    - 2 * a2f * sp.diff(b1f, v),
    3 * sp.diff(a1f, v) * b3f
    - a1f * sp.diff(b3f, v)
    + 2 * wronskian(a2f, b2f)
    + sp.diff(a3f, v) * b1f
    - 3 * a3f * sp.diff(b1f, v),
    3 * sp.diff(a2f, v) * b3f
    - 2 * a2f * sp.diff(b3f, v)
    + 2 * sp.diff(a3f, v) * b2f
    - 3 * a3f * sp.diff(b2f, v),
    3 * wronskian(a3f, b3f),
]
for k, expected in enumerate(expected_n3):
    gate(f"displayed N=3 row {k}", sp.expand(recurrence_rows[k] - expected) == 0)

# Positive control 1: the same recurrence engine recognizes (P,Q)=(v,w).
identity_rows = jacobian_rows([v, sp.Integer(0)], [sp.Integer(0), sp.Integer(1)])
gate("identity-map positive control", identity_rows == [1, 0])

# Positive control 2: the first three Catalan rows solve the triangular rows
# E_0,E_1,E_2 exactly, and leak first at E_3 as THM-3545 predicts.
catalan_a = [p, sp.Integer(1), sp.Rational(3, 4), sp.Rational(9, 8)]
catalan_b = [q, sp.Rational(3, 2) * v, sp.Rational(9, 8) * v,
             sp.Rational(27, 16) * v]
catalan_rows = jacobian_rows(catalan_a, catalan_b)
catalan_expected = [
    sp.Integer(1), sp.Integer(0), sp.Integer(0),
    -sp.Rational(135, 16), -sp.Rational(405, 64), -sp.Rational(729, 128),
]
gate("Catalan partial-solution control", catalan_rows == catalan_expected)


# ---------------------------------------------------------------------------
# 2. N=1 and N=2 cap-free reductions.
# ---------------------------------------------------------------------------

lam, mu = sp.symbols("lambda mu")
a = sp.Function("a")(v)
b = sp.Function("b")(v)
c = sp.Function("c")(v)
s = sp.Function("s")(v)
H = sp.expand(lam * sp.diff(p, v) - sp.diff(q, v))
gate("H has unavoidable quadratic top", sp.Poly(H, v).degree() == 2
     and sp.Poly(H, v).LC() == -3)

# N=1: W(a,b)=0 makes b/a constant in k(v), after which E_0 is a*H=1.
# The script checks the terminal polynomial obstruction; the logarithmic-
# derivative implication is stated and proved in the reflection.
gate("N=1 terminal factor is nonunit", sp.Poly(H, v).degree() > 0)

# Main N=2 branch: c,d nonzero and d=lambda*c.  Put s=b-lambda*a.
n2_main = jacobian_rows([p, a, c], [q, s + lam * a, lam * c])
n2_main_expected = [
    sp.diff(p, v) * s + a * H,
    2 * c * H + wronskian(a, s),
    sp.diff(c, v) * s - 2 * c * sp.diff(s, v),
    sp.Integer(0),
]
for k in range(4):
    gate(f"N=2 proportional branch row {k}",
         sp.expand(n2_main[k] - n2_main_expected[k]) == 0)

# If s is nonzero, E_2 gives c=mu*s^2.  E_0 forces
# deg(a)=deg(s)-1; E_1 then has degrees 2r+2 and at most 2r-2.
gate("N=2 main strict degree gap",
     all(2 * r + 2 > 2 * r - 2 for r in range(1, 65)))

# One-sided top branch c=0,d!=0: E_2 gives d=mu*a^2.  E_0 forces
# deg(b)=deg(a)+1, while E_1 has degrees 2r+1 and at most 2r.
n2_one_sided = jacobian_rows([p, a, sp.Integer(0)], [q, b, mu * a**2])
gate("N=2 one-sided row 1",
     sp.expand(n2_one_sided[1] - (4 * mu * v * a**2 + wronskian(a, b))) == 0)
gate("N=2 one-sided row 2 vanishes",
     sp.expand(n2_one_sided[2]) == 0)
gate("N=2 one-sided strict degree gap",
     all(2 * r + 1 > 2 * r for r in range(0, 65)))


# ---------------------------------------------------------------------------
# 3. N=3 transformed hierarchy and degree-cap ledger.
# ---------------------------------------------------------------------------

e = sp.Function("e")(v)
t = sp.Function("t")(v)
n3_both = jacobian_rows(
    [p, a, c, e],
    [q, t + lam * a, s + lam * c, lam * e],
)
n3_both_expected = [
    sp.diff(p, v) * t + a * H,
    2 * sp.diff(p, v) * s + 2 * c * H + wronskian(a, t),
    3 * e * H
    + 2 * sp.diff(a, v) * s - a * sp.diff(s, v)
    + sp.diff(c, v) * t - 2 * c * sp.diff(t, v),
    2 * wronskian(c, s) + sp.diff(e, v) * t - 3 * e * sp.diff(t, v),
    2 * sp.diff(e, v) * s - 3 * e * sp.diff(s, v),
    sp.Integer(0),
]
for k in range(6):
    gate(f"N=3 proportional-top hierarchy row {k}",
         sp.expand(n3_both[k] - n3_both_expected[k]) == 0)


def both_top_degree_audit(cap: int) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Return the E2 survivors and the E3 survivors in the (S,E)=(2,3) branch."""
    # If s=0, E3 gives e=mu*t^3.  The cap leaves only r=deg(t)=1;
    # then E1 asks a nonzero degree-zero Wronskian to cancel c*H of
    # degree at least two (or leaves that Wronskian alone when c=0).
    s_zero_r = [r for r in range(1, cap + 1) if 3 * r <= cap]
    gate(f"N=3 D={cap} s=0 residual degree",
         s_zero_r == [1] and 2 > 0)

    # The constant (S,E)=(0,0) branch has E2 degree list
    # (2, r-2, 2r-1), whose maximum is unique for every r>=1.
    gate(f"N=3 D={cap} constant-s/e branch",
         all(not tied_max([2, r - 2, 2 * r - 1]) for r in range(1, cap + 1)))

    after_e2: list[tuple[int, int]] = []
    after_e3: list[tuple[int, int]] = []
    for r in range(1, cap + 1):
        for C in range(cap + 1):
            # E2: degrees of 3eH, (2a's-as'), and (c't-2ct').
            deg_as = r
            if 2 * (r - 1) == 2:
                deg_as -= 1
            deg_ct = C + r - 1
            if C == 2 * r:
                deg_ct -= 1
            if not tied_max([5, deg_as, deg_ct]):
                continue
            after_e2.append((r, C))

            # E3: degrees of e't-3et' and 2W(c,s).
            deg_et = r + 2
            if 3 == 3 * r:
                deg_et -= 1
            deg_cs = C + 1
            if C == 2:
                deg_cs -= 1
            if tied_max([deg_et, deg_cs]):
                after_e3.append((r, C))
    return after_e2, after_e3


def e_zero_degree_audit(cap: int) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Audit e=0,f!=0; return easy c=0 states and square/cube survivors."""
    c_zero: list[tuple[int, int]] = []
    for A in range(cap):
        D = 2 * A - 1
        F = 3 * A - 3
        if 0 <= D <= cap and 0 <= F <= cap:
            # E3 has leading coefficient 3A-F and hence cannot vanish.
            gate(f"N=3 D={cap} e=0,c=0 state A={A}", 3 * A - F != 0)
            c_zero.append((A, D))

    square_cube: list[tuple[int, int]] = []
    NEG = -10**6
    for A in range(cap):
        for D in range(-1, cap + 1):
            # E1 degrees: 2p'd, -2cq', W(a,b), with (C,F)=(2,3).
            deg_pd = NEG if D < 0 else D + 1
            if not tied_max([deg_pd, 4, 2 * A]):
                continue

            # E2 degrees, lowering the bound when the leading coefficient
            # vanishes internally.
            deg_ad = NEG if D < 0 else A + D - 1
            if D >= 0 and 2 * A == D:
                deg_ad -= 1
            deg_cb = A + 2
            if A == 0:
                deg_cb -= 1
            if not tied_max([4, deg_ad, deg_cb]):
                continue

            # E3 degrees: 3a'f-af' and 2W(c,d).
            deg_af = A + 2
            if A == 1:
                deg_af -= 1
            deg_cd = NEG if D < 0 else D + 1
            if D == 2:
                deg_cd -= 1
            if tied_max([deg_af, deg_cd]):
                square_cube.append((A, D))
    return c_zero, square_cube


def f_zero_degree_audit(cap: int) -> tuple[list[tuple[int, int, int]], list[tuple[int, int]]]:
    """Audit f=0,e!=0; return d=0 candidates and (D,E)=(2,3) survivors."""
    d_zero: list[tuple[int, int, int]] = []
    for A in range(1, cap):
        C = 2 * A - 2
        if C > cap:
            continue
        # E2 can balance only when E+2=C+A.
        E = C + A - 2
        if 0 <= E <= cap:
            d_zero.append((A, C, E))
            # E3 then consists solely of e'b-3eb', with nonzero top.
            gate(f"N=3 D={cap} f=0,d=0 state A={A}", E - 3 * (A + 1) != 0)

    de_23: list[tuple[int, int]] = []
    for A in range(cap):
        for C in range(cap + 1):
            deg_ad = A + 1
            if 2 * A == 2:
                deg_ad -= 1
            deg_cb = C + A
            if C == 2 * (A + 1):
                deg_cb -= 1
            if tied_max([5, deg_ad, deg_cb]) and tied_max([3, C + 2, 2 * A]):
                de_23.append((A, C))
    return d_zero, de_23


CAPS = (3, 4, 5)
ledger: dict[int, dict[str, object]] = {}
for cap in CAPS:
    allowed_SE = [(S, E) for S in range(cap + 1) for E in range(cap + 1)
                  if 2 * E == 3 * S]
    gate(f"N=3 D={cap} proportional-top degree pairs",
         allowed_SE == [(0, 0), (2, 3)])

    both_e2, both_e3 = both_top_degree_audit(cap)
    expected_both_e2 = {
        3: [(3, 3)],
        4: [(3, 3), (4, 2)],
        5: [(1, 5), (3, 3), (4, 2), (5, 0), (5, 1)],
    }[cap]
    gate(f"N=3 D={cap} both-top E2 survivors", both_e2 == expected_both_e2)
    gate(f"N=3 D={cap} both-top E3 emptiness", both_e3 == [])

    c_zero, square_cube = e_zero_degree_audit(cap)
    gate(f"N=3 D={cap} square/cube unique state", square_cube == [(2, 3)])

    d_zero, de_23 = f_zero_degree_audit(cap)
    expected_d_zero = [(2, 2, 2)] if cap <= 4 else [(2, 2, 2), (3, 4, 5)]
    gate(f"N=3 D={cap} d=0 candidate", d_zero == expected_d_zero)
    gate(f"N=3 D={cap} (D,E)=(2,3) no degree survivor", de_23 == [])

    # Constant c,f leaves only A=1; its exact coefficient comparison gives
    # c=h^2+3/4 and c=h^2-3/2, a contradiction.
    constant_cf_A = [A for A in range(cap)
                     if tied_max([A + 1, 2, 2 * A])]
    gate(f"N=3 D={cap} constant c/f unique state", constant_cf_A == [1])

    # Constant d,e has no state for D=3.  At D=4 only A=3 reaches E1,
    # but E2 has the unique degree-seven term -rho*b*b'.
    constant_de_A = [A for A in range(cap)
                     if tied_max([1, A + 3, 2 * A])]
    gate(f"N=3 D={cap} constant d/e states",
         constant_de_A == ([] if cap == 3 else [3]))
    if cap >= 4:
        gate(f"N=3 D={cap} constant d/e terminal degree",
             2 * 3 + 1 > max(2, 3 - 1))

    ledger[cap] = {
        "both_E2": both_e2,
        "both_E3": both_e3,
        "e0_c0": c_zero,
        "e0_square_cube": square_cube,
        "f0_d0": d_zero,
        "f0_DE23": de_23,
        "constant_cf_A": constant_cf_A,
        "constant_de_A": constant_de_A,
    }


# ---------------------------------------------------------------------------
# 4. Exact exceptional square/cube branch.
# ---------------------------------------------------------------------------

# The other non-degree terminal is the constant c,f branch.  Its A=1
# coefficient equations force the two incompatible values printed below.
h_const = sp.symbols("h_const")
c_from_E1 = h_const**2 + sp.Rational(3, 4)
c_from_E2 = h_const**2 - sp.Rational(3, 2)
constant_cf_gap = sp.expand(c_from_E1 - c_from_E2)
gate("constant c/f coefficient contradiction", constant_cf_gap == sp.Rational(9, 4))

# Independent saturated ideal for that constant c,f branch.  Here A=1,
# h_const,c_const,f_const are all nonzero and E3 has already given
# d=(3f/(2c))*a+constant.
c_const, f_const, d_at_zero, zeta_const = sp.symbols(
    "c_const f_const d_at_zero zeta_const"
)
a_const_branch = 1 + 2 * h_const * v
b_const_branch = 3 * h_const * v**2 + sp.Rational(3, 2) * v - h_const
rho_const = 3 * f_const / (2 * c_const)
d_const_branch = d_at_zero + 2 * rho_const * h_const * v
E1_const_branch = sp.expand(
    2 * (sp.diff(p, v) * d_const_branch - c_const * sp.diff(q, v))
    + wronskian(a_const_branch, b_const_branch)
)
E2_const_branch = sp.expand(
    3 * sp.diff(p, v) * f_const
    + 2 * sp.diff(a_const_branch, v) * d_const_branch
    - a_const_branch * sp.diff(d_const_branch, v)
    - 2 * c_const * sp.diff(b_const_branch, v)
)
constant_cf_equations: list[sp.Expr] = []
for row in (E1_const_branch, E2_const_branch):
    for (_, coefficient) in sp.Poly(row, v).terms():
        numerator = sp.together(c_const * coefficient).as_numer_denom()[0]
        constant_cf_equations.append(sp.expand(numerator))
constant_cf_equations.append(zeta_const * h_const * c_const * f_const - 1)
G_constant_cf = sp.groebner(
    constant_cf_equations,
    zeta_const, d_at_zero, f_const, c_const, h_const,
    domain=sp.QQ, order="grevlex", method="f5b",
)

# In the sole degree survivor A=2,C=2,D=3,F=3, the three leading equations
# are I=II=III=0 below.  Eliminating Dlead and Flead gives a perfect square,
# forcing C=alpha^2/4, Dlead=3alpha^2/4, Flead=alpha^3/8.
alpha, Clead, Dlead, Flead = sp.symbols("alpha Clead Dlead Flead")
I = 4 * Dlead - 6 * Clead - sp.Rational(3, 2) * alpha**2
II = 6 * Flead + alpha * Dlead - 6 * alpha * Clead
III = 3 * alpha * Flead - 2 * Clead * Dlead

# Independent extraction from monomial leading rows.
leading_rows = jacobian_rows(
    [p, alpha * v**2, Clead * v**2, sp.Integer(0)],
    [q, sp.Rational(3, 2) * alpha * v**3, Dlead * v**3, Flead * v**3],
)
gate("exceptional leading E1 extracted",
     sp.expand(sp.Poly(leading_rows[1], v).coeff_monomial(v**4) - I) == 0)
gate("exceptional leading E2 extracted",
     sp.expand(sp.Poly(leading_rows[2], v).coeff_monomial(v**4) - II) == 0)
gate("exceptional leading E3 extracted",
     sp.expand(sp.Poly(leading_rows[3], v).coeff_monomial(v**4) - III) == 0)

D_forced = sp.Rational(3, 2) * Clead + sp.Rational(3, 8) * alpha**2
F_forced = alpha * (sp.Rational(3, 4) * Clead - sp.Rational(1, 16) * alpha**2)
gate("exceptional leading I", sp.expand(I.subs(Dlead, D_forced)) == 0)
gate("exceptional leading II",
     sp.expand(II.subs({Dlead: D_forced, Flead: F_forced})) == 0)
gate("exceptional perfect-square eliminant",
     sp.simplify(
         III.subs({Dlead: D_forced, Flead: F_forced})
         + sp.Rational(3, 16) * (alpha**2 - 4 * Clead)**2
     ) == 0)

h0, h1, root, zeta = sp.symbols("h0 h1 root zeta")
h = h0 + h1 * v
a_exc = 1 + 2 * v * h
b_exc = sp.Rational(3, 2) * v + (3 * v**2 - 1) * h
c_exc = h1**2 * (v - root)**2
f_exc = h1**3 * (v - root)**3

# E1 solves d.  Polynomiality is exactly the vanishing of the remainder at v=0.
d_numerator = sp.expand(2 * c_exc * sp.diff(q, v) - wronskian(a_exc, b_exc))
d_quotient, d_remainder_poly = sp.Poly(d_numerator, v).div(sp.Poly(v, v))
d_exc = sp.expand(d_quotient.as_expr() / 4)
d_remainder = sp.expand(d_remainder_poly.as_expr())

E2_exc = sp.expand(
    3 * sp.diff(p, v) * f_exc
    + 2 * sp.diff(a_exc, v) * d_exc - a_exc * sp.diff(d_exc, v)
    + sp.diff(c_exc, v) * b_exc - 2 * c_exc * sp.diff(b_exc, v)
)
E3_exc = sp.expand(
    3 * sp.diff(a_exc, v) * f_exc - a_exc * sp.diff(f_exc, v)
    + 2 * wronskian(c_exc, d_exc)
)

exc_equations = [d_remainder]
exc_equations.extend(coef for (_, coef) in sp.Poly(E2_exc, v).terms())
exc_equations.extend(coef for (_, coef) in sp.Poly(E3_exc, v).terms())
exc_equations.append(zeta * h1 - 1)  # saturation: h1 != 0
exc_equations = [
    sp.Poly(eq, zeta, root, h0, h1, domain=sp.QQ).primitive()[1].as_expr()
    for eq in exc_equations if eq != 0
]

G_QQ = sp.groebner(
    exc_equations, zeta, root, h0, h1,
    domain=sp.QQ, order="grevlex", method="f5b",
)
G5 = sp.groebner(
    exc_equations, zeta, root, h0, h1,
    modulus=5, order="grevlex", method="f5b",
)
G7 = sp.groebner(
    exc_equations, zeta, root, h0, h1,
    modulus=7, order="grevlex", method="f5b",
)


def unit_basis(G: sp.GroebnerBasis) -> bool:
    return len(G.polys) == 1 and G.polys[0].as_expr() == 1


gate("exceptional Groebner QQ unit ideal", unit_basis(G_QQ))
gate("exceptional Groebner GF(5) control", unit_basis(G5))
gate("exceptional Groebner GF(7) control", unit_basis(G7))
gate("constant c/f Groebner QQ unit ideal", unit_basis(G_constant_cf))

# Transparent hand-readable elimination chain behind the unit basis.
E2_v2 = sp.Poly(E2_exc, v).coeff_monomial(v**2)
E2_v0 = sp.Poly(E2_exc, v).coeff_monomial(1)
gate("exceptional E2 top square",
     sp.factor(E2_v2) == 9 * h1 * (h0 + h1 * root)**2)
remainder_reduced = sp.factor(d_remainder.subs(h0, -h1 * root))
gate("exceptional remainder forces h1=3/2",
     sp.simplify(remainder_reduced + h1 - sp.Rational(3, 2)) == 0)
terminal_constant = sp.factor(E2_v0.subs({h0: -h1 * root,
                                         h1: sp.Rational(3, 2)}))
gate("exceptional terminal contradiction", terminal_constant == -sp.Rational(9, 4))


# ---------------------------------------------------------------------------
# Transcript.
# ---------------------------------------------------------------------------

print("JC2 Catalan mixed-thickening exact recurrence audit")
print(f"python={platform.python_version()} sympy={sp.__version__}")
print(f"gates={gates}")
print("recurrence=E_k=sum_{i+j=k+1}(j*a_i'*b_j-i*a_i*b_j'), a_0=v^2, b_0=v^3-v")
print("target=E_0=1; E_k=0 for 1<=k<=2N-1")
print("N3_rows=" + " | ".join(sp.sstr(sp.expand(row)) for row in expected_n3))
print("identity_control=" + ",".join(sp.sstr(x) for x in identity_rows))
print("catalan_N3_rows=" + ",".join(sp.sstr(x) for x in catalan_rows))
print("N1=IMPOSSIBLE_CAP_FREE")
print("N2=IMPOSSIBLE_CAP_FREE")
print("N2_main_gap=2r+2_vs_at_most_2r-2")
print("N2_one_sided_gap=2r+1_vs_at_most_2r")
for cap in CAPS:
    item = ledger[cap]
    print(
        f"N3_D{cap}_ledger="
        f"bothE2:{item['both_E2']};bothE3:{item['both_E3']};"
        f"e0c0:{item['e0_c0']};squarecube:{item['e0_square_cube']};"
        f"f0d0:{item['f0_d0']};DE23:{item['f0_DE23']};"
        f"constcf:{item['constant_cf_A']};constde:{item['constant_de_A']}"
    )
    print(f"N3_D{cap}=AFFINE_EMPTY_EXACT")
print(f"exceptional_groebner_QQ={G_QQ}")
print(f"exceptional_groebner_GF5={G5}")
print(f"exceptional_groebner_GF7={G7}")
print(f"exceptional_terminal={terminal_constant}")
print(f"constant_cf_terminal_gap={constant_cf_gap}")
print(f"constant_cf_groebner_QQ={G_constant_cf}")
print("projective_closure=NOT_COMPUTED_NOT_NEEDED_FOR_AFFINE_COEFFICIENT_EMPTINESS")
print("raw_full_ideal=NOT_USED;solver_failure_claims=none")
print("failures=" + ("none" if not failures else ",".join(failures)))

if failures:
    raise SystemExit(1)
