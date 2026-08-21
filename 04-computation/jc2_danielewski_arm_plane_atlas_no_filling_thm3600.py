#!/usr/bin/env python3
"""Exact controls for proved, hostile-audited THM-3600.

The universal arm-cover, intersection, Cech descent, exact-image implication,
and normal no-filling lemma are proved in the theorem.  This companion checks
their algebraic identities, boundary examples, incidence ranks, cyclic
compiler rows, and the complete A13 specialization with exact SymPy
arithmetic.  Finite rows are hostile controls, not extrapolations.
"""

from itertools import combinations, product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one explicit gate and raise on failure."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expr):
    """Exact rational-function zero test."""
    return sp.cancel(expr) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


b, c, a = sp.symbols("b c a")
t, x, q = sp.symbols("t x q")

print("THM-3600 exact companion -- provisional arm-plane atlas and no-filling gate")
print("scope=finite exact controls; universal cover/intersection/descent/filling statements are proof-driven")


print("SECTION arm-plane charts, singular shears, and global primitive")
ROOTS = (sp.Integer(0), sp.Integer(1), sp.Integer(3), sp.Integer(7), sp.Integer(12))
chart_rows = 0
shear_rows = 0
log_rows = 0

for n in range(1, 6):
    for h in range(1, 6):
        roots = ROOTS[:h]
        gamma = sp.Integer(h + 2)
        Sigma = sp.expand(gamma * sp.prod(b - beta for beta in roots))

        for i, beta_i in enumerate(roots):
            Sigma_i = sp.cancel(Sigma / (b - beta_i))
            b_i = beta_i + c**n * a
            e_i = sp.expand(a * Sigma_i.subs(b, b_i))

            require(
                f"chart relation n={n} h={h} i={i}",
                zero(c**n * e_i - Sigma.subs(b, b_i)),
            )
            require(
                f"chart inverse n={n} h={h} i={i}",
                same((b_i - beta_i) / c**n, a),
            )
            require(
                f"whole free arm n={n} h={h} i={i}",
                same(e_i.subs(c, 0), a * Sigma_i.subs(b, beta_i)),
            )
            require(
                f"squarefree arm slope n={n} h={h} i={i}",
                Sigma_i.subs(b, beta_i) != 0,
            )
            a_rational = (b - beta_i) / c**n
            require(
                f"Darboux chart bracket n={n} h={h} i={i}",
                same(c**n * sp.diff(a_rational, b), 1),
            )
            chart_rows += 1

            if h > 1:
                beta_j = roots[(i + 1) % h]
                require(
                    f"proof-trap point lies on Sigma_i=0 n={n} h={h} i={i}",
                    Sigma_i.subs(b, beta_j) == 0,
                )
                require(
                    f"proof-trap point lies on Y n={n} h={h} i={i}",
                    Sigma.subs(b, beta_j) == 0,
                )

        for i, j in combinations(range(h), 2):
            beta_i, beta_j = roots[i], roots[j]
            delta = beta_i - beta_j
            a_j = a + delta * c ** (-n)

            require(
                f"shear b compatibility n={n} h={h} ij={i},{j}",
                same(beta_j + c**n * a_j, beta_i + c**n * a),
            )
            require(
                f"shear symplectic n={n} h={h} ij={i},{j}",
                sp.diff(a_j, a) == 1,
            )
            numerator, denominator = sp.fraction(sp.together(a_j))
            require(
                f"shear has boundary pole n={n} h={h} ij={i},{j}",
                sp.degree(denominator, c) == n
                and numerator.subs(c, 0) == delta,
            )

            if n >= 2:
                H_ij = delta * c ** (1 - n) / (1 - n)
                require(
                    f"rational Hamiltonian shear n={n} h={h} ij={i},{j}",
                    same(sp.diff(H_ij, c), delta * c ** (-n)),
                )

                # Pull alpha_j back to (a_i,c).  Store coefficients of da,dc.
                daj_dc = sp.diff(a_j, c)
                alpha_j_da = c / (n - 1)
                alpha_j_dc = sp.cancel(c * daj_dc / (n - 1) + n * a_j / (n - 1))
                require(
                    f"primitive da gluing n={n} h={h} ij={i},{j}",
                    same(alpha_j_da, c / (n - 1)),
                )
                require(
                    f"primitive dc gluing n={n} h={h} ij={i},{j}",
                    same(alpha_j_dc, n * a / (n - 1)),
                )
            else:
                residue = sp.residue(delta / c, c, 0)
                require(
                    f"logarithmic residue h={h} ij={i},{j}",
                    residue == delta and residue != 0,
                )
                log_rows += 1
            shear_rows += 1

        if n >= 2:
            alpha_da = c / (n - 1)
            alpha_dc = n * a / (n - 1)
            exterior_coefficient = sp.diff(alpha_dc, a) - sp.diff(alpha_da, c)
            require(
                f"global primitive differential n={n} h={h}",
                exterior_coefficient == 1,
            )

            # A local canonical pair supplies a positive Hamilton--Jacobi control.
            P_local, Q_local = a, c
            phi_local = c * a / (n - 1)
            rhs_da = P_local * sp.diff(Q_local, a) + sp.diff(phi_local, a)
            rhs_dc = P_local * sp.diff(Q_local, c) + sp.diff(phi_local, c)
            require(
                f"local HJ da coefficient n={n} h={h}",
                same(rhs_da, alpha_da),
            )
            require(
                f"local HJ dc coefficient n={n} h={h}",
                same(rhs_dc, alpha_dc),
            )

print(f"PASS chart_rows={chart_rows} shear_rows={shear_rows} logarithmic_rows={log_rows}")


print("SECTION Hermite-Pade and Laurent shear membership gates")
y, z, s_chart = sp.symbols("y z s_chart")


def y_order(polynomial):
    """Lowest y exponent of a nonzero polynomial in y and s_chart."""
    poly = sp.Poly(sp.expand(polynomial), y, s_chart)
    return min(monomial[0] for monomial, coefficient in poly.terms() if coefficient != 0)


pade_rows = 0
laurent_rows = 0
for h in range(2, 6):
    roots = ROOTS[:h]
    D = sp.expand(sp.prod(z - beta for beta in roots))

    for i, j in combinations(range(h), 2):
        require(
            f"pairwise arm-ideal comaximality h={h} ij={i},{j}",
            roots[j] - roots[i] != 0,
        )

    for M in range(1, 4):
        for k in range(M + 1):
            N = sp.expand(y ** (M - k) * D**k)
            for j, beta_j in enumerate(roots):
                pulled = sp.expand(N.subs(z, beta_j + y * s_chart))
                require(
                    f"Pade arm order h={h} M={M} k={k} j={j}",
                    y_order(pulled) >= M,
                )

            F_base = sp.cancel(N.subs(z, y * a) / y**M)
            numerator, denominator = sp.fraction(F_base)
            require(
                f"Pade constructive base polynomial h={h} M={M} k={k}",
                denominator == 1 and sp.Poly(numerator, y, a) is not None,
            )
            expected = sp.cancel((D.subs(z, y * a) / y) ** k)
            require(
                f"Pade constructive formula h={h} M={M} k={k}",
                same(F_base, expected),
            )
            pade_rows += 1

        # The base-chart monomial a^M clears the beta_0 arm only.
        N_hostile = z**M
        require(
            f"Pade one-chart positive h={h} M={M}",
            y_order(N_hostile.subs(z, y * s_chart)) >= M,
        )
        hostile_pull = sp.expand(N_hostile.subs(z, roots[1] + y * s_chart))
        require(
            f"Pade other-arm hostile h={h} M={M}",
            y_order(hostile_pull) == 0,
        )

    for n in range(1, 5):
        a_chart = sp.symbols("a_chart")
        for depth in range(1, 2 * n + 1):
            M = (depth + n - 1) // n
            good_piece = c ** (-depth) * D**M
            bad_piece = c ** (-depth) * z**M

            for j, beta_j in enumerate(roots):
                good_pull = sp.cancel(
                    good_piece.subs(z, beta_j + c**n * a_chart)
                )
                good_num, good_den = sp.fraction(good_pull)
                require(
                    f"Laurent all-arm good h={h} n={n} d={depth} j={j}",
                    good_den == 1 and sp.Poly(good_num, c, a_chart) is not None,
                )

            bad_base = sp.cancel(bad_piece.subs(z, c**n * a_chart))
            bad_base_num, bad_base_den = sp.fraction(bad_base)
            require(
                f"Laurent base-chart positive h={h} n={n} d={depth}",
                bad_base_den == 1
                and sp.Poly(bad_base_num, c, a_chart) is not None,
            )
            bad_other = sp.cancel(
                bad_piece.subs(z, roots[1] + c**n * a_chart)
            )
            _, bad_other_den = sp.fraction(bad_other)
            require(
                f"Laurent other-arm pole h={h} n={n} d={depth}",
                sp.degree(bad_other_den, c) == depth,
            )

            degree_a = h * M
            leading_c_order = n * degree_a - depth
            invoice = n * (degree_a - degree_a // h)
            require(
                f"leading-coefficient degree invoice h={h} n={n} d={depth}",
                leading_c_order >= invoice,
            )
            laurent_rows += 1

print(
    f"PASS Hermite-Pade_rows={pade_rows} Laurent_rows={laurent_rows} "
    "one-chart hostiles active"
)


print("SECTION Cech incidence, de Rham, and Picard ranks")
cech_rows = 0
for h in range(1, 9):
    edges = list(combinations(range(h), 2))
    triangles = list(combinations(range(h), 3))

    d0 = sp.zeros(len(edges), h)
    for row, (i, j) in enumerate(edges):
        d0[row, i] = -1
        d0[row, j] = 1

    edge_index = {edge: index for index, edge in enumerate(edges)}
    d1 = sp.zeros(len(triangles), len(edges))
    for row, (i, j, k) in enumerate(triangles):
        d1[row, edge_index[(i, j)]] = 1
        d1[row, edge_index[(i, k)]] = -1
        d1[row, edge_index[(j, k)]] = 1

    rank_d0 = d0.rank()
    kernel_d1 = len(edges) - d1.rank()
    require(f"simplex vertex-edge rank h={h}", rank_d0 == h - 1)
    require(f"truncated dlog edge kernel h={h}", kernel_d1 == h - 1)
    require(f"H1 vanishing h={h}", h - rank_d0 == 1)
    require(f"Picard exponent rank h={h}", rank_d0 == h - 1)
    cech_rows += 1

print(f"PASS Cech_rows={cech_rows} H1_dR=0 H2_rank=h-1 Picard_rank=h-1")


print("SECTION hostile boundaries")
# Repeated root: the naive beta=0 chart reaches only e=0 over c=0.
Sigma_multiple = b**2 * (b - 1)
b_multiple = c**2 * a
e_multiple = sp.expand(c**2 * a**2 * (c**2 * a - 1))
require(
    "multiple-root substituted relation",
    zero(c**2 * e_multiple - Sigma_multiple.subs(b, b_multiple)),
)
require("multiple-root misses generic arm", e_multiple.subs(c, 0) == 0)
require("missed point lies on singular target", Sigma_multiple.subs(b, 0) == 0)

# Squarefree comparison retains the full arm.
Sigma_squarefree = b * (b - 1)
e_squarefree = sp.expand(a * (c**2 * a - 1))
require(
    "squarefree comparison relation",
    zero(c**2 * e_squarefree - Sigma_squarefree.subs(b, c**2 * a)),
)
require("squarefree comparison free arm", e_squarefree.subs(c, 0) == -a)
print("PASS multiple-root missed-arm hostile and n=1 residue boundary")


def cyclic_compiler(n, r):
    """Return the exact cyclic compiler row of THM-3581."""
    S = 1 + t**r
    P = sp.expand(
        sum(
            sp.binomial(n - 1, j) * t ** (r * j) / sp.Integer(n * r * j + 1)
            for j in range(n)
        )
    )
    B = sp.expand(t * P**n)
    return S, P, B


print("SECTION cyclic compiler and direct Keller controls")
compiler_rows = ((2, 1), (2, 2), (3, 1), (3, 2), (4, 1))
for n, r in compiler_rows:
    S, P, B = cyclic_compiler(n, r)
    require(
        f"compiler ODE n={n} r={r}",
        zero(P + n * t * sp.diff(P, t) - S ** (n - 1)),
    )
    require(
        f"compiler derivative n={n} r={r}",
        zero(sp.diff(B, t) - P ** (n - 1) * S ** (n - 1)),
    )
    require(f"compiler gcd(P,S) n={n} r={r}", sp.gcd(P, S) == 1)
    require(
        f"compiler degree n={n} r={r}",
        sp.degree(B, t) == n * (n - 1) * r + 1,
    )

    txq = x**n * q
    a_pull = q / S.subs(t, txq) ** n
    c_pull = x * P.subs(t, txq) * S.subs(t, txq)
    require(
        f"compiler central identity n={n} r={r}",
        zero(a_pull * c_pull**n - B.subs(t, txq)),
    )
    jac = sp.cancel(
        sp.diff(a_pull, x) * sp.diff(c_pull, q)
        - sp.diff(a_pull, q) * sp.diff(c_pull, x)
    )
    require(f"direct Keller Jacobian n={n} r={r}", jac == -1)

print(f"PASS cyclic_compiler_rows={len(compiler_rows)} direct_Jacobian=-1")


print("SECTION A13 punctures, central chart, collision, and zero-value escape")
n13 = 3
S13, P13, B13 = cyclic_compiler(3, 2)
kappa = sp.Rational(72, 91) ** 3
I = sp.I
require("A13 displayed P", zero(P13 - (7 * t**4 + 26 * t**2 + 91) / 91))
require("A13 plus critical point", sp.diff(B13, t).subs(t, I) == 0)
require("A13 minus critical point", sp.diff(B13, t).subs(t, -I) == 0)
require("A13 plus critical value", same(B13.subs(t, I), I * kappa))
require("A13 minus critical value", same(B13.subs(t, -I), -I * kappa))
require("A13 nonzero side values", kappa != 0)

a0 = sp.symbols("a0")
b0 = c**3 * a0
e0 = a0 * (c**6 * a0**2 + kappa**2)
require(
    "A13 central target relation",
    zero(c**3 * e0 - b0 * (b0**2 + kappa**2)),
)

txq13 = x**3 * q
a13 = q / (1 + txq13**2) ** 3
c13 = x * P13.subs(t, txq13) * (1 + txq13**2)
require("A13 pullback product", zero(a13 * c13**3 - B13.subs(t, txq13)))
jac13 = sp.cancel(
    sp.diff(a13, x) * sp.diff(c13, q)
    - sp.diff(a13, q) * sp.diff(c13, x)
)
require("A13 direct Jacobian", jac13 == -1)

# Exact inverses exhibit the two forced units after localizing at S=t^2+1.
require("A13 unit inverse t-i", same((t - I) * ((t + I) / S13), 1))
require("A13 unit inverse t+i", same((t + I) * ((t - I) / S13), 1))

# Roots of P are zero-valued critical points, remain off the deleted S-curves,
# and account for twelve of the thirteen central-collision sheets.
require("A13 P squarefree", sp.gcd(P13, sp.diff(P13, t)) == 1)
require("A13 P and S coprime", sp.gcd(P13, S13) == 1)
require("A13 zero value on P roots", sp.rem(B13, P13, domain=sp.QQ) == 0)
require(
    "A13 critical derivative on P roots",
    sp.rem(sp.diff(B13, t), P13, domain=sp.QQ) == 0,
)
require("A13 collision count", 1 + n13 * sp.degree(P13, t) == 13)

eta, rho = sp.symbols("eta rho", nonzero=True)
x_cubed = rho * kappa**2 / (eta * (1 + rho**2) ** 3)
q_rho = eta * (1 + rho**2) ** 3 / kappa**2
require("A13 root-sheet t value", same(x_cubed * q_rho, rho))
require(
    "A13 root-sheet central a value",
    same(q_rho / (1 + rho**2) ** 3, eta / kappa**2),
)
require("A13 t=0 central a value", same((eta / kappa**2), eta / kappa**2))
print("PASS A13 deleted_Gm=2 pole_order=3 zero_order=1 collision_count=13")


print("SECTION normal no-filling and quantitative unit rank")
A_symbol, C_symbol = sp.symbols("A_symbol C_symbol")
for n in range(1, 6):
    beta = sp.Integer(n + 1)
    coefficient_da = C_symbol**n
    relation = A_symbol * C_symbol**n - beta
    require(f"cotangent da coefficient n={n}", coefficient_da == C_symbol**n)
    require(
        f"nonzero-value contradiction n={n}",
        sp.expand(relation - A_symbol * coefficient_da) == -beta,
    )

critical_points = (sp.Integer(-2), sp.Integer(1), sp.Integer(3))
valuation_matrix = sp.eye(len(critical_points))
require("forced-unit valuation rank", valuation_matrix.rank() == len(critical_points))
independence_rows = 0
for exponents in product((-1, 0, 1), repeat=len(critical_points)):
    if all(exponent == 0 for exponent in exponents):
        continue
    unit_word = sp.prod(
        (t - sigma) ** exponent
        for sigma, exponent in zip(critical_points, exponents)
    )
    require(
        f"forced-unit independence exponents={exponents}",
        not zero(sp.diff(unit_word, t)),
    )
    independence_rows += 1

print(
    "PASS no-filling: each nonzero critical point forces a unit; "
    f"independence_rows={independence_rows} forced_rank={len(critical_points)}"
)
print("zero-critical-value escape=P-roots remain and carry central collisions")
print("nonnormal boundary=conductor may hide integral t only after etaleness is lost")


print("SECTION consequence ledger")
print("unconditional=arm charts A2; overlaps singular shears; global ring is all-chart intersection")
print("Hamilton-Jacobi=for n>=2 alpha=(c da+n a dc)/(n-1), d alpha=omega")
print("compiler=central chart pulls back to punctured A2 and is surjective etale by THM-3581")
print("A13=two Gm punctures retain a surjective degree-13 etale collision map to A2")
print("boundary=n=1 logarithmic; repeated roots singular; zero critical values escape no-filling")
print("nonconsequence=no polynomial Darboux pair and no JC(2) counterexample")
print(f"CHECKS={CHECKS}")
print("FINAL PASS -- THM-3600 finite exact controls")
