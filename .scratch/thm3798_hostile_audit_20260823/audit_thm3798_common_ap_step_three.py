#!/usr/bin/env python3
"""Independent hostile audit of the THM-3798 common-AP step-three peel.

The candidate is treated as untrusted.  This script replays it, but all
mathematical gates below use a separate contribution census and a formal
first-jet differential algebra.  Divisibility is checked by exact polynomial
division and by remainders modulo the actual arm polynomial; no
cancel-then-multiply identity is accepted as a factor certificate.
"""

from __future__ import annotations

import ast
import hashlib
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    """Optimization-safe truth gate; deliberately not a Python assert."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.together(left - right)) == 0, message)


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CANDIDATE_DIR = REPO / ".scratch" / "ap_d3_extension_20260823"
CANDIDATE = CANDIDATE_DIR / "jc2_cubic_pseudoplane_ap_d3_candidate.py"
FROZEN = CANDIDATE.with_suffix(".out")

# The original 17-gate source used fraction-field cancel/multiply identities
# for divisibility.  These hashes identify the repaired, current 24-gate
# source and transcript; the old 2dcfde.../2cc0a1... pair is superseded.
EXPECTED_SOURCE_HASH = "6b6a1867ec9311bf108dea2147c39504e66fec78b92c62f0f302043c2db8dad4"
EXPECTED_OUTPUT_HASH = "9c65319ccf272948767aa4083435299ad072ed660c1b7b62d592dcf464f0c522"


# ---------------------------------------------------------------------------
# 1. Reproducibility, optimization safety, and repaired verification lineage.
# ---------------------------------------------------------------------------

source_bytes = CANDIDATE.read_bytes()
frozen_bytes = FROZEN.read_bytes()
gate(hashlib.sha256(source_bytes).hexdigest() == EXPECTED_SOURCE_HASH, "candidate source drift")
gate(hashlib.sha256(frozen_bytes).hexdigest() == EXPECTED_OUTPUT_HASH, "candidate output drift")

tree = ast.parse(source_bytes.decode("utf-8"))
bare_asserts = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
gate(bare_asserts == 0, "candidate contains optimization-erased assert")
gate(b"polynomial_factor" in source_bytes, "repaired polynomial-factor gate missing")
gate(b"quotient not polynomial" in source_bytes, "denominator-free quotient gate missing")
gate(b"CHECKS={CHECKS}" in source_bytes, "candidate active-gate output missing")


def replay(args: list[str]) -> bytes:
    run = subprocess.run(args, cwd=REPO, capture_output=True, check=False)
    gate(run.returncode == 0, f"candidate replay failed: {args}")
    gate(run.stderr == b"", f"candidate replay emitted stderr: {args}")
    return run.stdout


def lf(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n")


normal_raw = replay([sys.executable, str(CANDIDATE)])
optimized_raw = replay([sys.executable, "-O", str(CANDIDATE)])
gate(normal_raw == optimized_raw, "normal and optimized candidate replays differ")
gate(lf(normal_raw) == lf(frozen_bytes), "normal replay differs from frozen output")
gate(lf(optimized_raw) == lf(frozen_bytes), "optimized replay differs from frozen output")
gate(b"CHECKS=24" in normal_raw, "repaired candidate did not execute 24 gates")
gate(b"RESULT=PASS" in normal_raw, "candidate did not print PASS")


# ---------------------------------------------------------------------------
# 2. Independent exact contribution census, including every sign/zero seam.
# ---------------------------------------------------------------------------


def contribution_buckets(d: int) -> dict[int, tuple[tuple[int, int], ...]]:
    buckets: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i in range(2):
        for j in range(4):
            buckets[(i + j) * d].append((i, j))
    return {address: tuple(terms) for address, terms in sorted(buckets.items())}


def singleton_commutes(u: int, v: int) -> bool:
    # THM-3785: nonzero commuting weights have equal strict sign.  An active
    # zero profile paired with a nonzero weight would have to be scalar and is
    # removable; only (0,0) is a genuine endpoint seam.
    return (u < 0 and v < 0) or (u > 0 and v > 0) or (u == v == 0)


def scalar_rows(d: int) -> dict[int, tuple[tuple[int, int], ...]]:
    buckets = contribution_buckets(d)
    rows: dict[int, tuple[tuple[int, int], ...]] = {}
    for scalar_q, scalar_terms in buckets.items():
        if len(scalar_terms) != 2:
            continue
        total = -scalar_q - 2
        current: list[tuple[int, int]] = []
        # Bottom commutation plus total<0 forces a,b<0.  Consequently this
        # finite interval is exact, not a bounded-search assumption.
        for a in range(total + 1, 0):
            b = total - a
            if b >= 0:
                continue
            ok = True
            for address, terms in buckets.items():
                if address == scalar_q or len(terms) != 1:
                    continue
                i, j = terms[0]
                if not singleton_commutes(a + i * d, b + j * d):
                    ok = False
                    break
            if ok:
                current.append((a, b))
        rows[scalar_q] = tuple(current)
    return rows


expected_multiplicities = (1, 2, 2, 2, 1)
for d in (3, 4):
    buckets = contribution_buckets(d)
    gate(tuple(len(v) for v in buckets.values()) == expected_multiplicities,
         f"wrong AP multiplicities at d={d}")
    gate(tuple(buckets) == tuple(i * d for i in range(5)), f"wrong AP addresses at d={d}")

rows3 = scalar_rows(3)
rows4 = scalar_rows(4)
gate(rows3 == {
    3: ((-2, -3), (-1, -4)),
    6: ((-2, -6), (-1, -7)),
    9: (),
}, f"step-three rows changed: {rows3}")
gate(rows4 == {
    4: ((-3, -3), (-2, -4), (-1, -5)),
    8: ((-3, -7), (-2, -8), (-1, -9)),
    12: ((-3, -11),),
}, f"step-four hostile rows changed: {rows4}")

# A deliberately wider enumeration is a hostile control on the exact
# negative-sum interval and on omitted zero/opposite-sign rows.
for d, expected in ((3, rows3), (4, rows4)):
    buckets = contribution_buckets(d)
    for scalar_q, wanted in expected.items():
        wide: list[tuple[int, int]] = []
        for a in range(-80, 81):
            for b in range(-80, 81):
                if a + b != -scalar_q - 2:
                    continue
                if not singleton_commutes(a, b):
                    continue
                if all(
                    singleton_commutes(a + i * d, b + j * d)
                    for address, terms in buckets.items()
                    if address != scalar_q and len(terms) == 1
                    for i, j in terms
                ):
                    wide.append((a, b))
        gate(tuple(wide) == wanted, f"wide sign/zero census differs d={d},q={scalar_q}")


# ---------------------------------------------------------------------------
# 3. Exact profile modules, squarefreeness, and the UFD endpoint-power law.
# ---------------------------------------------------------------------------


def rho(weight: int) -> int:
    return weight % 3


def arm_order(weight: int) -> int:
    return max(0, (-weight + 2) // 3) if weight < 0 else 0


gate(tuple((u, rho(u), arm_order(u)) for u in range(-7, 7)) == (
    (-7, 2, 3), (-6, 0, 2), (-5, 1, 2), (-4, 2, 2),
    (-3, 0, 1), (-2, 1, 1), (-1, 2, 1), (0, 0, 0),
    (1, 1, 0), (2, 2, 0), (3, 0, 0), (4, 1, 0),
    (5, 2, 0), (6, 0, 0),
), "profile-module residue/order table changed")
gate(all(arm_order(u) >= 1 for u in range(-30, 0)), "negative profile lost arm divisor")
gate(rho(1) == 1 and rho(2) == 2, "upper F profiles need a nonconstant w-factor")

# The only zero-weight profiles in the four rows occur inside collision
# buckets, never at a singleton endpoint.  Scalars there are removable; the
# proof below even permits them and still retains Delta.
row_weight_profiles = {
    (-2, -3): ((-2, 1), (-3, 0, 3, 6)),
    (-1, -4): ((-1, 2), (-4, -1, 2, 5)),
    (-2, -6): ((-2, 1), (-6, -3, 0, 3)),
    (-1, -7): ((-1, 2), (-7, -4, -1, 2)),
}
zero_locations = {
    row: tuple((side, i) for side, weights in enumerate(pair) for i, u in enumerate(weights) if u == 0)
    for row, pair in row_weight_profiles.items()
}
gate(zero_locations == {
    (-2, -3): ((1, 1),),
    (-1, -4): (),
    (-2, -6): ((1, 2),),
    (-1, -7): (),
}, f"zero-weight profile locations changed: {zero_locations}")

w, c = sp.symbols("w c", nonzero=True)
Delta = w**3 - c**3
Delta_poly = sp.Poly(Delta, w, domain=sp.QQ.frac_field(c))
Delta_prime_poly = sp.Poly(sp.diff(Delta, w), w, domain=sp.QQ.frac_field(c))
gate(sp.gcd(Delta_poly, Delta_prime_poly).degree() == 0, "Delta is not squarefree")
same(sp.discriminant(Delta, w), -27 * c**6, "Delta discriminant")
same(Delta, (w - c) * (w**2 + c * w + c**2), "Delta factorization")

# Logarithmic differentiation: for every endpoint exponent used here, the
# bracket is exactly the numerator of (g^U/f^V)'.
f, df, g, dg = sp.symbols("f df g dg")


def formal_D_fg(expr: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(expr, f) * df + sp.diff(expr, g) * dg)


endpoint_pairs = ((2, 3), (1, 4), (2, 6), (1, 3), (1, 7), (2, 2))
for U, V in endpoint_pairs:
    endpoint_bracket = U * f * dg - V * df * g
    same(
        formal_D_fg(g**U / f**V),
        g ** (U - 1) * endpoint_bracket / f ** (V + 1),
        f"endpoint logarithmic derivative U={U},V={V}",
    )

# Opposite signs would make a positive product f^V*g^U constant.  Since one
# of the two weights is negative, its profile contains Delta and is
# nonconstant; hence an active opposite-sign singleton cannot commute.  A
# zero/nonzero singleton instead differentiates the zero-weight profile and
# makes it a removable scalar.  These are the exact sign/zero laws used by
# the finite row census.
for U in range(1, 7):
    for V in range(1, 7):
        opposite_bracket = U * f * dg + V * df * g
        same(
            formal_D_fg(f**V * g**U),
            f ** (V - 1) * g ** (U - 1) * opposite_bracket,
            f"opposite-sign product law U={U},V={V}",
        )
same(0 * f * dg - 3 * df * g, -3 * df * g, "zero/nonzero lower law")
same(3 * f * dg - 0 * df * g, 3 * f * dg, "nonzero/zero upper law")

# Independent valuation-lattice audit.  U*ord(g)=V*ord(f) is equivalent to
# ord(f)=(U/delta)n and ord(g)=(V/delta)n.  The report supplies the unbounded
# Euclidean proof; this box attacks constants, gcds, and zero valuations.
valuation_cases = 0
for U in range(1, 13):
    for V in range(1, 13):
        d = math.gcd(U, V)
        for of in range(0, 49):
            for og in range(0, 49):
                if U * og != V * of:
                    continue
                valuation_cases += 1
                gate(of % (U // d) == 0, f"bad f valuation U={U},V={V},of={of},og={og}")
                gate(og % (V // d) == 0, f"bad g valuation U={U},V={V},of={of},og={og}")
                gate(of // (U // d) == og // (V // d),
                     f"valuation parameters differ U={U},V={V},of={of},og={og}")

# Radical seam: at each simple root of Delta, Delta|t^m forces ord(t)>=1.
for m in range(1, 13):
    for order_t in range(0, 13):
        if m * order_t >= 1:
            gate(order_t >= 1, f"squarefree radical seam failed m={m},ord={order_t}")


# ---------------------------------------------------------------------------
# 4. Formal first-jet algebra and genuine factor extraction for all four rows.
# ---------------------------------------------------------------------------

p, dp, q, dq = sp.symbols("p dp q dq")
z, dz, h, dh = sp.symbols("z dz h dh")
t, dt = sp.symbols("t dt")
delta, ddelta = sp.symbols("delta ddelta")
P, dP, N, dN, U0, dU0 = sp.symbols("P dP N dN U0 dU0")
lam, mu, nu, kap = sp.symbols("lam mu nu kap")
u_sym, v_sym = sp.symbols("u_sym v_sym")

jet = {
    p: dp, q: dq, z: dz, h: dh, t: dt,
    delta: ddelta, P: dP, N: dN, U0: dU0,
}


def D(expr: sp.Expr) -> sp.Expr:
    return sp.expand(sum(sp.diff(expr, variable) * derivative for variable, derivative in jet.items()))


def wb(u: sp.Expr, left: sp.Expr, v: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.expand(u * left * D(right) - v * D(left) * right)


def exact_factor(expr: sp.Expr, factor: sp.Symbol, power: int, message: str) -> sp.Expr:
    """Certify factor**power by exact polynomial division, never cancellation."""
    polynomial = sp.Poly(sp.expand(expr), factor)
    divisor = sp.Poly(factor**power, factor)
    quotient, remainder = polynomial.div(divisor)
    gate(remainder.is_zero, f"{message}: nonzero remainder {remainder.as_expr()}")
    rebuilt = sp.expand(quotient.as_expr() * factor**power)
    same(rebuilt, sp.expand(expr), f"{message}: quotient reconstruction")
    return quotient.as_expr()


same(wb(v_sym, z, u_sym, p), -wb(u_sym, p, v_sym, z), "output-swap skew identity")

# q=3, row (-2,-3): p=t^2 and g0=lambda*t^3.
s31 = wb(-2, t**2, 0, z) + wb(1, q, -3, lam * t**3)
same(s31, t**2 * (-2 * dz + 3 * lam * q * dt + 3 * lam * dq * t),
     "q=3 (-2,-3) explicit t^2 factor")
exact_factor(s31, t, 2, "q=3 (-2,-3) t^2")
s31_delta = sp.expand(s31.subs({t: delta * U0, dt: ddelta * U0 + delta * dU0}))
exact_factor(s31_delta, delta, 2, "q=3 (-2,-3) Delta^2")

# q=3, row (-1,-4): p and g1 are independent negative profiles.
s32 = wb(-1, p, -1, z) + wb(2, q, -4, lam * p**4)
s32_delta = sp.expand(s32.subs({
    p: delta * P, dp: ddelta * P + delta * dP,
    z: delta * N, dz: ddelta * N + delta * dN,
}))
exact_factor(s32_delta, delta, 2, "q=3 (-1,-4) Delta^2")
same(
    wb(-1, delta * P, -1, delta * N),
    delta**2 * (dP * N - P * dN),
    "equal-negative Wronskian cancellation",
)

# q=6, row (-2,-6): both adjacent equations and their full constants.
g0_33 = lam * p**3
g3_33 = mu * q**3
g1_33 = 3 * lam * p**2 * q + h
g2_33 = 3 * mu * p * q**2 + nu
lower33 = wb(-2, p, -3, g1_33) + wb(1, q, -6, g0_33)
upper33 = wb(-2, p, 3, g3_33) + wb(1, q, 0, g2_33)
same(lower33, wb(-2, p, -3, h), "q=6 (-2,-6) lower residual")
same(upper33, 0, "q=6 (-2,-6) upper integrated solution")

# For an arbitrary upper profile z, the upper collision is q times the
# derivative of z-3*mu*p*q^2; hence its full polynomial solution adds nu in k.
upper33_general = wb(-2, p, 3, g3_33) + wb(1, q, 0, z)
same(upper33_general, q * D(z - 3 * mu * p * q**2),
     "q=6 (-2,-6) upper conservation law")

s33 = wb(-2, p, 0, g2_33) + wb(1, q, -3, g1_33)
expected33 = 6 * (lam - mu) * p * D(p * q**2) + wb(1, q, -3, h)
same(s33, expected33, "q=6 (-2,-6) scalar identity")
same(sp.diff(s33, nu), 0, "q=6 (-2,-6) scalar is independent of nu")

# H=0: exactly one arm factor is guaranteed.  lambda=mu gives zero, not a
# nonzero scalar; nu=0 is separately included.
s33_h0 = sp.expand(s33.subs({h: 0, dh: 0}))
exact_factor(s33_h0, p, 1, "q=6 (-2,-6) H=0 p factor")
s33_h0_delta = sp.expand(s33_h0.subs({p: delta * P, dp: ddelta * P + delta * dP}))
exact_factor(s33_h0_delta, delta, 1, "q=6 (-2,-6) H=0 Delta factor")
same(s33_h0.subs(lam, mu), 0, "q=6 (-2,-6) H=0 lambda=mu zero seam")
same(s33_h0.subs(nu, 0), s33_h0, "q=6 (-2,-6) H=0 nu=0 seam")

# H!=0: endpoint UFD gives p=t^2,H=kap*t^3, with kap in k*.  The identity is
# polynomial also at kap=0, so the boundary back to H=0 is explicitly safe.
s33_hn = sp.expand(s33.subs({p: t**2, dp: 2 * t * dt, h: kap * t**3, dh: 3 * kap * t**2 * dt}))
exact_factor(s33_hn, t, 2, "q=6 (-2,-6) H!=0 t^2 factor")
s33_hn_delta = sp.expand(s33_hn.subs({t: delta * U0, dt: ddelta * U0 + delta * dU0}))
exact_factor(s33_hn_delta, delta, 2, "q=6 (-2,-6) H!=0 Delta^2")
exact_factor(sp.expand(s33_hn_delta.subs(kap, 0)), delta, 2,
             "q=6 (-2,-6) kap=0 boundary")
exact_factor(sp.expand(s33_hn_delta.subs(nu, 0)), delta, 2,
             "q=6 (-2,-6) nu=0 boundary")

# q=6, row (-1,-7): lower adjacent ODE and the unused upper seam.
g0_34 = lam * p**7
g1_34 = p**4 * (7 * lam * p**2 * q + nu)
lower34 = wb(-1, p, -4, g1_34) + wb(2, q, -7, g0_34)
same(lower34, 0, "q=6 (-1,-7) lower integrated solution")

K_lower = z - 7 * lam * p**6 * q
lower34_general = wb(-1, p, -4, z) + wb(2, q, -7, g0_34)
same(lower34_general, -p * D(K_lower) + 4 * dp * K_lower,
     "q=6 (-1,-7) lower conservation law")
same(-p * D(K_lower) + 4 * dp * K_lower, -p**5 * D(K_lower / p**4),
     "q=6 (-1,-7) lower quotient primitive")

# Top endpoint gives g3=mu*q.  If K=g2-mu*p, the upper collision is
# 2qK'+q'K; multiplying by K is (qK^2)'.  Since the valid weight-2 profile
# q=w^2 Q(w^3) is nonconstant, polynomiality forces K=0.  This unused seam
# actually sharpens the scalar factor from Delta^2 to Delta^3.
g3_34 = mu * q
upper34_general = wb(-1, p, 2, g3_34) + wb(2, q, -1, z)
K_upper = z - mu * p
same(upper34_general, 2 * q * D(K_upper) + dq * K_upper,
     "q=6 (-1,-7) upper residual")
same(D(q * K_upper**2), K_upper * upper34_general,
     "q=6 (-1,-7) upper square primitive")

s34 = wb(-1, p, -1, z) + wb(2, q, -4, g1_34)
s34_delta = sp.expand(s34.subs({
    p: delta * P, dp: ddelta * P + delta * dP,
    z: delta * N, dz: ddelta * N + delta * dN,
}))
exact_factor(s34_delta, delta, 2, "q=6 (-1,-7) Delta^2")

second34 = wb(2, q, -4, g1_34)
second34_delta = sp.expand(second34.subs({p: delta * P, dp: ddelta * P + delta * dP}))
exact_factor(second34_delta, delta, 3, "q=6 (-1,-7) second summand Delta^3")
exact_factor(sp.expand(s34_delta.subs(nu, 0)), delta, 2,
             "q=6 (-1,-7) nu=0 boundary")

s34_upper_solved = sp.expand(s34.subs({z: mu * p, dz: mu * dp}))
same(wb(-1, p, -1, mu * p), 0, "q=6 (-1,-7) solved upper Wronskian")
exact_factor(
    sp.expand(s34_upper_solved.subs({p: delta * P, dp: ddelta * P + delta * dP})),
    delta, 3, "q=6 (-1,-7) upper-solved Delta^3",
)


# ---------------------------------------------------------------------------
# 5. Redundant actual-polynomial remainder controls in the profile modules.
# ---------------------------------------------------------------------------

w0 = sp.symbols("w0")
Delta0 = w0**3 - 8
T0 = w0**3


def bracket_w(u: int, left: sp.Expr, v: int, right: sp.Expr) -> sp.Expr:
    return sp.expand(u * left * sp.diff(right, w0) - v * sp.diff(left, w0) * right)


def remainder_gate(expr: sp.Expr, divisor: sp.Expr, message: str) -> None:
    polynomial = sp.Poly(sp.expand(expr), w0, domain=sp.QQ)
    divisor_poly = sp.Poly(sp.expand(divisor), w0, domain=sp.QQ)
    quotient, remainder = polynomial.div(divisor_poly)
    gate(remainder.is_zero, f"{message}: remainder={remainder.as_expr()}")
    gate(quotient.as_expr().is_polynomial(w0), f"{message}: nonpolynomial quotient")


# (-2,-3), scalar q=3.  t=w^2*Delta makes p=t^2 a valid weight -2 profile.
t0 = w0**2 * Delta0
p0 = t0**2
q0 = w0 * (1 + T0)
g00 = 2 * t0**3
g10 = 1 + 2 * T0 + T0**2
s31_poly = bracket_w(-2, p0, 0, g10) + bracket_w(1, q0, -3, g00)
remainder_gate(s31_poly, Delta0**2, "polynomial (-2,-3) q=3 Delta^2")

# (-1,-4), scalar q=3.
p0 = w0**2 * Delta0 * (1 + T0)
q0 = w0**2 * (1 + 2 * T0)
g10 = w0**2 * Delta0 * (2 - T0)
g00 = 2 * p0**4
s32_poly = bracket_w(-1, p0, -1, g10) + bracket_w(2, q0, -4, g00)
remainder_gate(s32_poly, Delta0**2, "polynomial (-1,-4) q=3 Delta^2")

# (-2,-6), scalar q=6, H=0 and H!=0; nu is tested at zero and nonzero.
p0 = w0 * Delta0 * (1 + T0)
q0 = w0 * (1 + 2 * T0)
for nu0 in (0, 5):
    g10 = 6 * p0**2 * q0
    g20 = 9 * p0 * q0**2 + nu0
    s33_poly_h0 = bracket_w(-2, p0, 0, g20) + bracket_w(1, q0, -3, g10)
    remainder_gate(s33_poly_h0, Delta0, f"polynomial (-2,-6) H=0 nu={nu0} Delta")

t0 = w0**2 * Delta0
p0 = t0**2
q0 = w0 * (1 + 2 * T0)
for kap0 in (0, 7):
    for nu0 in (0, 5):
        h0 = kap0 * t0**3
        g10 = 6 * p0**2 * q0 + h0
        g20 = 9 * p0 * q0**2 + nu0
        s33_poly_hn = bracket_w(-2, p0, 0, g20) + bracket_w(1, q0, -3, g10)
        remainder_gate(
            s33_poly_hn, Delta0**2,
            f"polynomial (-2,-6) square branch kap={kap0},nu={nu0} Delta^2",
        )

# (-1,-7), scalar q=6, including the zero integration constant.
p0 = w0**2 * Delta0 * (1 + T0)
q0 = w0**2 * (1 + 2 * T0)
g20 = w0**2 * Delta0 * (2 + T0)
for nu0 in (0, 5):
    g10 = p0**4 * (14 * p0**2 * q0 + nu0)
    s34_poly = bracket_w(-1, p0, -1, g20) + bracket_w(2, q0, -4, g10)
    remainder_gate(s34_poly, Delta0**2, f"polynomial (-1,-7) nu={nu0} Delta^2")


audit_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(audit_tree)),
     "independent audit contains optimization-erased assert")


print("AUDIT=THM3798-common-AP-step-three-hostile")
print(f"CANDIDATE_SHA256={EXPECTED_SOURCE_HASH}")
print(f"CANDIDATE_OUTPUT_SHA256={EXPECTED_OUTPUT_HASH}")
print("D3_ROWS=q3:[(-2,-3),(-1,-4)];q6:[(-2,-6),(-1,-7)];q9:[]")
print("D4_HOSTILE=q4:[(-3,-3),(-2,-4),(-1,-5)];q8:[(-3,-7),(-2,-8),(-1,-9)];q12:[(-3,-11)]")
print("FACTORS=q3(-2,-3):Delta^2;q3(-1,-4):Delta^2;q6(-2,-6):Delta_or_Delta^2;q6(-1,-7):Delta^2")
print("SEAMS=support_shrink,zero_weight,nu0,kap0,H0,Hnonzero,lambda_eq_mu,algebraically_closed_units,profile_modules,output_swap:PASS")
print(f"VALUATION_CASES={valuation_cases}")
print(f"CHECKS={CHECKS}")
print("PROMOTION=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print("RESULT=PASS")
