"""Exact symbolic certificate for THM-2576.

The script computes the closure of F(V(L)) for the fixed sporadic Keller map.
It derives a small inverse-coordinate norm problem, removes its three chart
artifacts, and verifies that the residual factor H is primitive irreducible.
It also checks the independent normalized parametrization on hostile controls.

All arithmetic is over QQ.  No ``assert`` is used, so ``python3 -O`` is a
substantive replay.  The 361-term H is specified canonically by the resultant
identity and by a SHA-256 hash of its lexicographic coefficient ledger.
"""

import hashlib

import sympy as sp


X, s, a, b, c = sp.symbols("X s a b c")

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
S = 27 * a * c**2 - 9 * b * c + 8
E = L * X**3 + T * X - 2 * c


def require(condition, message):
    """Raise in ordinary and optimized Python if an exact check fails."""
    if not condition:
        raise RuntimeError(message)


def zero(expr, message):
    require(sp.expand(expr) == 0, message)


print("=" * 78)
print("[1] Inverse conic pair and its degree-one subresultant")
print("=" * 78)

C1 = 3 * a * X**2 - (1 + s) * (b * X - s)
C2 = a * X**3 + c * (1 + s) ** 3 - X * (1 + s) * (2 + s)

D = (
    -3 * X**2 * a * c
    + X**2 * b**2 * c
    - X**2 * b
    + 2 * X * b * c
    - 2 * X
    + c
)
N = (
    -3 * X**3 * a * b * c
    + 4 * X**3 * a
    - 6 * X**2 * a * c
    + X**2 * b**2 * c
    - X**2 * b
    + 2 * X * b * c
    - 2 * X
    + c
)

subresultants = sp.subresultants(C1, C2, s)
require(len(subresultants) == 4, "unexpected conic-pair subresultant length")
zero(subresultants[-2] - (D * s + N), "degree-one subresultant changed")
zero(subresultants[-1] - a * X**3 * E, "fiber cubic changed")

print("  Sres_1(C1,C2;s)=D(X)s+N(X).  [PASS]")
print("  Res_s(C1,C2)=a X^3[L X^3+(4-3bc)X-2c].  [PASS]")
print("  Thus generically s=-N/D on the inverse X-sheet.")
print()

print("=" * 78)
print("[2] Pull the source Jelonek equation into the inverse X-sheet")
print("=" * 78)

source_y = s / X
source_z = (X * (2 - 3 * s) - c) / X**3
source_L = (
    27 * X**2 * source_z**2
    - 18 * X * source_y * source_z
    + 16 * X
    + source_y**3 * source_z
    - source_y**2
)
P, P_den = sp.cancel(source_L).as_numer_denom()
require(P_den == X**6, "unexpected source-L denominator")
require(sp.degree(P, s) == 4 and len(sp.Poly(P, X, s, c).terms()) == 10,
        "unexpected source-L numerator")

# Q=D^4 P(X,-N/D), expanded without a rational-expression blowup.
P_s = sp.Poly(P, s)
Q = sp.expand(
    sum(coeff * (-N) ** power * D ** (4 - power)
        for (power,), coeff in P_s.terms())
)
require(sp.degree(Q, X) == 15, "unexpected inverse norm degree")
require(len(sp.Poly(Q, X, a, b, c).terms()) == 398,
        "unexpected inverse norm support")

print("  X^6 L_source=P(X,s), with deg_s(P)=4 and 10 terms.  [PASS]")
print("  Q(X)=D^4 P(X,-N/D) has X-degree 15 and 398 terms.  [PASS]")
print()

print("=" * 78)
print("[3] Saturated resultant and the irreducible image divisor H")
print("=" * 78)

# SymPy's PRS implementation normalizes the polynomial order when the first
# degree is smaller, but does not restore the odd*odd swap sign.  Here
# deg(E)=3<15=deg(Q), so ``resultant(E,Q)`` has the sign of Res(Q,E), not the
# standard Sylvester/root-product Res(E,Q).  Compute the higher-degree-first
# value and restore Res(E,Q)=(-1)^(3*15) Res(Q,E) explicitly.
prs_low_first = sp.resultant(E, Q, X)
raw_resultant = -sp.resultant(Q, E, X)
require(sp.expand(raw_resultant + prs_low_first) == 0,
        "PRS/Sylvester resultant sign correction changed")
H_fraction = sp.cancel(raw_resultant / (a**8 * c**18 * S**8))
H_num, H_den = H_fraction.as_numer_denom()
require(H_den == 1, "artifact quotient is not polynomial")
H = sp.expand(H_num)
H_poly = sp.Poly(H, a, b, c)

factor_content, factors = sp.factor_list(H)
require(factor_content == 1 and len(factors) == 1, "H is reducible")
require(factors[0][1] == 1 and sp.expand(factors[0][0] - H) == 0,
        "H factor normalization changed")
require(sp.gcd(sp.Poly(L, a, b, c), H_poly).total_degree() == 0,
        "H and L share a component")

multidegree = tuple(H_poly.degree(v) for v in (a, b, c))
term_count = len(H_poly.terms())
ledger = "\n".join(f"{monomial}:{coefficient}"
                   for monomial, coefficient in H_poly.terms())
ledger_hash = hashlib.sha256(ledger.encode("ascii")).hexdigest()

require(multidegree == (14, 21, 12), "H multidegree changed")
require(H_poly.total_degree() == 25, "H total degree changed")
require(term_count == 361, "H support changed")
require(
    ledger_hash == "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2",
    "H coefficient ledger changed",
)

print("  Res_X(E,Q)=+a^8 c^18 S^8 H (standard Sylvester convention).  [PASS]")
print("  a, c, S are inverse-chart artifacts, removed before geometry.")
print("  H is primitive and irreducible over QQ; gcd(H,L)=1.  [PASS]")
print("  multidegree(H)=(14,21,12), total degree=25, terms=361.")
print("  H coefficient-ledger SHA-256:")
print("   ", ledger_hash)
print()

print("=" * 78)
print("[4] Independent normalized-image controls")
print("=" * 78)

tau, lam = sp.symbols("tau lambda")
source_x = lam**2 * (3 - tau * lam) / 27
source_y = lam * (4 - tau * lam) / 3
source_z = tau
u = 1 + source_x * source_y

image_a = sp.expand(u**3 * source_z + source_y**2 * u * (4 + 3 * source_x * source_y))
image_b = sp.expand(
    source_y
    + 3 * source_x * u**2 * source_z
    + 3 * source_x * source_y**2 * (4 + 3 * source_x * source_y)
)
image_c = sp.expand(2 * source_x - 3 * source_x**2 * source_y - source_x**3 * source_z)

# Exact hostile grid: this is an independent pointwise path, while the theorem
# proves the full identity by the resultant relation on a dense open set.
sample_count = 0
for tau_value in (-2, -1, 0, 1, 2):
    for lambda_value in (-2, -1, 0, 1, 2):
        substitution = {tau: tau_value, lam: lambda_value}
        target = {
            a: image_a.subs(substitution),
            b: image_b.subs(substitution),
            c: image_c.subs(substitution),
        }
        require(H_poly.eval(target) == 0, "normalized image misses H")
        sample_count += 1

image_jacobian = sp.Matrix((image_a, image_b, image_c)).jacobian((tau, lam))
rank_minor = sp.factor(
    image_jacobian.extract((0, 1), (0, 1)).det().subs({tau: 1, lam: 0})
)
require(rank_minor == sp.Rational(4, 3), "normalized image is not a surface")

generic_substitution = {tau: 1, lam: 1}
generic_source_x = source_x.subs(generic_substitution)
generic_source_u = u.subs(generic_substitution)
generic_target = {
    a: image_a.subs(generic_substitution),
    b: image_b.subs(generic_substitution),
    c: image_c.subs(generic_substitution),
}
generic_D = D.subs({X: generic_source_x, **generic_target})
generic_artifact = (a * c * S).subs(generic_target)
require(generic_source_x * generic_source_u * generic_D * generic_artifact != 0,
        "dense inverse chart lacks a positive control")

axis_target = {
    a: image_a.subs({tau: 1, lam: 0}),
    b: image_b.subs({tau: 1, lam: 0}),
    c: image_c.subs({tau: 1, lam: 0}),
}
require(axis_target == {a: 1, b: 0, c: 0}, "axis hostile changed")
require(H_poly.eval(axis_target) == 0 and L.subs(axis_target) == 16,
        "H and L were not separated by the axis hostile")

off_target = None
off_value = None
for av in range(-2, 3):
    for bv in range(-2, 3):
        for cv in range(-2, 3):
            value = H_poly.eval({a: av, b: bv, c: cv})
            if value != 0:
                off_target = (av, bv, cv)
                off_value = value
                break
        if off_target is not None:
            break
    if off_target is not None:
        break
require(off_target is not None, "H vanished on the hostile search box")

print(f"  H(F(nu(tau,lambda)))=0 on {sample_count} exact hostile samples.  [PASS]")
print("  (tau,lambda)=(1,1) lies in the honest X*u*a*c*S*D chart.  [PASS]")
print("  rank d(F o nu)=2 at (tau,lambda)=(1,0); minor=4/3.  [PASS]")
print("  image hostile (1,0,0): H=0 but L=16, so H != L.  [PASS]")
print(f"  off-image control {off_target}: H={off_value} !=0.  [PASS]")
print()

print("ALL EXACT THM-2576 COMPOSITE JELONEK CHECKS PASSED")
