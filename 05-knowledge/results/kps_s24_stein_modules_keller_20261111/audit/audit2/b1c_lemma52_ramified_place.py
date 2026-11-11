"""b1c_lemma52_ramified_place.py - second audit, Lemma 5.2: a place with e = 2, f = 2 and a non-constant c.

Surface X = {(u, eta, zeta) : zeta^2 = u eta^2 + 1, eta zeta != 0}, map (u, v) = (u, eta^2): etale of degree 4.
Dual fields:  d_u u = 1, d_u eta = 0, d_u zeta = eta^2/(2 zeta);   d_v u = 0, d_v eta = 1/(2 eta), d_v zeta = u/(2 zeta).
Over the line at infinity (tau = 1/u, s = v/u) there is ONE place: with sigma^2 = s and pi = 1/eta,
     pi^2 = tau / s  (c = 1/s, not constant),   residue field C(sigma)  (f = 2),   e = 2,
     u = sigma^-2 pi^-2,  eta = pi^-1,  zeta = sigma^-1 pi^-2 sqrt(1 + sigma^2 pi^4).
w = (1/2) ord_pi.  Tested: w(d_u phi) >= w(phi) + 1, w(d_v phi) >= w(phi) + 1 on Laurent monomials and random
combinations, and the formulas of the proof on the series:
     D~_1 = (1/2) pi d/dpi,   D~_2 (c_j pi^j) = (d c_j/ds - (j/2) s^-1 c_j) pi^j   (d sigma/ds = 1/(2 sigma)),
     d_v = tau D~_2 = sigma^2 pi^2 D~_2,   d_u = - tau D~_1 - s tau D~_2.
"""
import random
from fractions import Fraction as Fr
import sympy as sp
from b1b_lemma52_slice_places import LP, Ser

random.seed(523)
NT = 40                                      # number of terms in pi


def binom_half(k):
    out = Fr(1)
    for i in range(k):
        out *= (Fr(1, 2) - i) / (i + 1)
    return out


sigma = LP.var(1)
u_ser = Ser(-2, [LP.var(-2)], exact=True)
eta_ser = Ser(-1, [LP.const(1)], exact=True)
# zeta = sigma^-1 pi^-2 sum_k binom(1/2, k) sigma^(2k) pi^(4k)
zc = [LP() for _ in range(NT)]
for k in range(NT // 4 + 1):
    if 4 * k < NT:
        zc[4 * k] = LP.var(2 * k - 1, binom_half(k))
zeta_ser = Ser(-2, zc)
# check the relation zeta^2 = u eta^2 + 1 on the series
rel = zeta_ser * zeta_ser - u_ser * eta_ser * eta_ser - Ser.const(1)
assert not rel.c, "relation fails on the series"
eta_inv = Ser(1, [LP.const(1)], exact=True)
zeta_inv = zeta_ser.inverse()

U, E, Z = sp.symbols("u eta zeta")


def du(h):
    return sp.expand(sp.diff(h, U) + sp.diff(h, Z) * E ** 2 / (2 * Z))


def dv(h):
    return sp.expand(sp.diff(h, E) / (2 * E) + sp.diff(h, Z) * U / (2 * Z))


def ev(h):
    """series of a Laurent polynomial in (u, eta, zeta)."""
    h = sp.expand(h)
    tot = Ser(0, [], exact=True)
    if h == 0:
        return tot
    for term in sp.Add.make_args(h):
        c, rest = term.as_independent(U, E, Z)
        pw = rest.as_powers_dict()
        s = Ser.const(1)
        for base, e in pw.items():
            e = int(e)
            if base == 1:
                continue
            if base == U:
                b = u_ser
            elif base == E:
                b = eta_ser if e > 0 else eta_inv
            elif base == Z:
                b = zeta_ser if e > 0 else zeta_inv
            else:
                raise ValueError(base)
            for _ in range(abs(e)):
                s = s * b
        tot = tot + s.scale(Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)))
    return tot


def ds(a):                                   # d/ds on Laurent polynomials in sigma: (1/(2 sigma)) d/dsigma
    return a.deriv() * LP.var(-1, Fr(1, 2))


if __name__ == "__main__":
    # sanity: the dual fields on u, v = eta^2 and compatibility with the relation
    assert du(U) == 1 and du(E ** 2) == 0 and dv(E ** 2) == 1 and dv(U) == 0
    assert sp.simplify(du(Z ** 2 - U * E ** 2 - 1)) == 0 and sp.simplify(dv(Z ** 2 - U * E ** 2 - 1).subs(Z ** 2, U * E ** 2 + 1)) == 0
    tests = [U ** i * E ** j * Z ** k for i in range(0, 3) for j in range(-2, 3) for k in range(-2, 3)]
    for _ in range(40):
        tests.append(sum(random.randint(-4, 4) * U ** random.randint(0, 2) * E ** random.randint(-2, 2) * Z ** random.randint(-2, 2) for _ in range(3)))
    tests += [Z ** 2 - U * E ** 2 + Z ** -1, (Z ** 2 - U * E ** 2) * E + U, (Z ** 2 - U * E ** 2 - 1) + E ** -3, Z ** 3 - U * E ** 2 * Z]
    tot = bad = eq = und = fbad = fchk = 0
    seen_vals = set()
    for h in tests:
        if sp.expand(h) == 0:
            continue
        S0, Su, Sv = ev(h), ev(du(h)), ev(dv(h))
        v0 = S0.valuation()
        if v0 is None or v0 == "zero":
            und += 1
            continue
        w0 = Fr(v0, 2)
        seen_vals.add(w0)
        for S1, hz in ((Su, sp.expand(du(h)) == 0), (Sv, sp.expand(dv(h)) == 0)):
            tot += 1
            if hz:
                continue
            v1 = S1.valuation()
            if v1 is None or v1 == "zero":
                if S1.prec() < 2 * (w0 + 1):
                    und += 1
                continue
            w1 = Fr(v1, 2)
            if w1 < w0 + 1:
                bad += 1
            elif w1 == w0 + 1:
                eq += 1
        # formulas of the proof: tau = sigma^2 pi^2, s = sigma^2
        D1 = Ser(S0.val, [a.scale(Fr(S0.val + j, 2)) for j, a in enumerate(S0.c)], S0.exact)
        D2 = Ser(S0.val, [ds(a) - (a * LP.var(-2)).scale(Fr(S0.val + j, 2)) for j, a in enumerate(S0.c)], S0.exact)
        tau = Ser(2, [LP.var(2)], exact=True)
        pred_v = tau * D2
        pred_u = -(tau * D1) - (tau * D2).scale(LP.var(2))
        for pred, act in ((pred_u, Su), (pred_v, Sv)):
            top = min(pred.prec(), act.prec())
            top = int(top) if top != float("inf") else max(pred.val + len(pred.c), act.val + len(act.c))
            lo = min(pred.val, act.val)
            fchk += 1
            if any((pred.coeff(e) - act.coeff(e)) for e in range(lo, top)):
                fbad += 1
    print("place with e = 2, f = 2, c = 1/s on zeta^2 = u eta^2 + 1 over (u, eta^2):")
    print("   %d inequalities tested: %d violations, %d equalities, %d undecided at this precision; half-integer valuations occur: %s"
          % (tot, bad, eq, und, any(v.denominator == 2 for v in seen_vals)))
    print("   formulas of the proof (d_v = tau D~_2, d_u = -tau D~_1 - s tau D~_2 with the term (j/e)(c'/c)): %d mismatches of %d" % (fbad, fchk))
    S = ev(Z ** 2 - U * E ** 2)
    print("   cancellation example: w(zeta^2 - u eta^2) = %s although w(zeta^2) = w(u eta^2) = -2 (the valuation is not monomial here)" % Fr(S.valuation(), 2))
