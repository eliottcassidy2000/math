"""Audit check 3 (THM-4708 part 3), by a method different from the author's resultants / Hilbert polynomials.

Logic.  m_N is a Hauptmodul of Gamma_0(N)+ with rational q-expansion; t = (eta1/etaN)^{24/(N-1)} is a Hauptmodul of
Gamma_0(N) with rational q-expansion and t + b_N/t = m_N - a_N.  Since j(tau) and j(N tau) are in C(X_0(N)) = C(t) with
rational q-expansions, j(tau), j(N tau) lie in Q(t).  If m_N(tau) is rational then [Q(t):Q] <= 2, hence
[Q(j(tau)):Q] = h(disc(tau)) <= 2 (and likewise for N tau).  So it suffices to run over all orders O with h(O) <= 2.

For each such discriminant D' and each reduced primitive form Q of discriminant D', the N+1 points of Gamma_0(N)\\H over
j(tau_Q) are tau_Q and -1/(tau_Q + k), k = 0..N-1.  We evaluate m_N there by eta products (mpmath, 300 digits), at the
Gamma_0(N)+-equivalent point (P or W_N P = -1/(N P)) of larger imaginary part.  The multiset of these (N+1) h(D') values
is Galois-stable, so prod (X - m_N(P)) lies in Q[X]; we round its coefficients (reporting the worst distance to Z) and
factor it with sympy.  Rational roots = rational CM values with disc(tau) = D'.

Also: the class-number <= 2 list is derived here by counting reduced forms (fundamental |D| <= 20000; all orders
f^2 D_K with h(D_K) <= 2 and f <= 30, which is more than enough since h(f^2 D_K) >= h(D_K) phi(f)/3).
"""
import sys, os, json
from math import gcd, isqrt
import mpmath as mp
import sympy as sp

mp.mp.dps = 300
HERE = os.path.dirname(os.path.abspath(__file__))


def reduced_forms(D):
    """Primitive reduced positive definite forms (a,b,c), b^2-4ac = D < 0."""
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            out.append((a, b, c))
        a += 1
    return out


def is_fundamental(D):
    if D % 4 == 1:
        return sp.factorint(-D) and all(e == 1 for e in sp.factorint(-D).values())
    if D % 4 == 0:
        d = D // 4
        if d % 4 in (2, 3):
            return all(e == 1 for e in sp.factorint(-d).values())
    return False


def class_number_le2_list():
    fund = [D for D in range(-3, -20001, -1) if D % 4 in (0, 1) and is_fundamental(D)
            and len(reduced_forms(D)) <= 2]
    allD = set()
    for DK in fund:
        for f in range(1, 31):
            D = f * f * DK
            if len(reduced_forms(D)) <= 2:
                allD.add(D)
    return sorted(fund, reverse=True), sorted(allD, reverse=True)


def eta_ratio_pow(tau, N, e):
    """(eta(tau)/eta(N tau))^e with e (N-1) = 24  ->  q^{-1} prod (1-q^n)^e / (1-q^{Nn})^e."""
    q = mp.exp(2j * mp.pi * tau)
    eps = mp.mpf(10) ** (-mp.mp.dps - 20)
    num = mp.mpc(1)
    qn = q
    while abs(qn) > eps:
        num *= (1 - qn)
        qn *= q
    den = mp.mpc(1)
    qN = q ** N
    qn = qN
    while abs(qn) > eps:
        den *= (1 - qn)
        qn *= qN
    return (num / den) ** e / q


PARAM = {5: (6, 22, 125), 7: (4, 13, 49)}


def m_value(tau, N):
    e, a, b = PARAM[N]
    t = eta_ratio_pow(tau, N, e)
    return t + a + b / t


def best_rep(P, N):
    W = -1 / (N * P)
    return P if P.imag >= W.imag else W


def disc_of_point(A, B, C):
    """tau root of A x^2 + B x + C (A>0), returns primitive discriminant."""
    g = gcd(gcd(A, abs(B)), abs(C))
    A, B, C = A // g, B // g, C // g
    return B * B - 4 * A * C, (A, B, C)


def form_of_image(A, B, C, mat):
    """If tau is a root of A x^2 + B x + C and z = (a tau + b)/(c tau + d), return the form of z (up to sign)."""
    a, b, c, d = mat
    # tau = (d z - b)/(-c z + a);  substitute
    # A (d z - b)^2 + B (d z - b)(-c z + a) + C (-c z + a)^2
    z2 = A * d * d - B * d * c + C * c * c
    z1 = -2 * A * d * b + B * (d * a + b * c) - 2 * C * c * a
    z0 = A * b * b - B * b * a + C * a * a
    if z2 < 0:
        z2, z1, z0 = -z2, -z1, -z0
    return z2, z1, z0


def main():
    fund, allD = class_number_le2_list()
    h1 = [D for D in allD if len(reduced_forms(D)) == 1]
    h2 = [D for D in allD if len(reduced_forms(D)) == 2]
    print(f'fundamental D with h<=2, |D|<=20000: {len(fund)} (max |D| = {-min(fund)})')
    print(f'all orders with h(D)=1: {len(h1)}: {h1}')
    print(f'all orders with h(D)=2: {len(h2)}: {h2}')
    results = {}
    for N in (5, 7):
        found = {}
        worst = mp.mpf(0)
        for Dp in allD:
            forms = reduced_forms(Dp)
            pts = []
            for (a, b, c) in forms:
                tauQ = (-b + mp.sqrt(mp.mpf(Dp))) / (2 * a)  # sqrt of negative -> imaginary
                tauQ = mp.mpc(-mp.mpf(b) / (2 * a), mp.sqrt(-mp.mpf(Dp)) / (2 * a))
                # identity coset: P = tau_Q ; N P = N tau_Q
                cosets = [((1, 0, 0, 1), 'id')] + [((0, -1, 1, k), f'S T^{k}') for k in range(N)]
                for mat, lab in cosets:
                    aa, bb, cc, dd = mat
                    P = (aa * tauQ + bb) / (cc * tauQ + dd)
                    fA, fB, fC = form_of_image(a, b, c, mat)
                    DP, _ = disc_of_point(fA, fB, fC)
                    # N P is a root of fA (x/N)^2 + fB (x/N) + fC = 0  ->  fA x^2 + N fB x + N^2 fC
                    DNP, _ = disc_of_point(fA, N * fB, N * N * fC)
                    val = m_value(best_rep(P, N), N)
                    pts.append((val, (a, b, c), lab, DP, DNP))
            # polynomial
            X = sp.Symbol('X')
            coeffs = [mp.mpc(1)]
            for (val, *_r) in pts:
                new = [mp.mpc(0)] * (len(coeffs) + 1)
                for i, cf in enumerate(coeffs):
                    new[i] += cf
                    new[i + 1] -= cf * val
                coeffs = new
            icoeffs = []
            for cf in coeffs:
                r = mp.nint(cf.real)
                err = abs(cf - r) / max(1, abs(r)) if abs(r) > 0 else abs(cf - r)
                worst = max(worst, abs(cf - r))
                icoeffs.append(int(r))
            poly = sp.Poly(icoeffs, X)
            fl = sp.factor_list(poly)
            lin = [(-f.as_poly(X).all_coeffs()[1]) / f.as_poly(X).all_coeffs()[0] for f, e in fl[1] if f.as_poly(X).degree() == 1]
            for r in lin:
                r = sp.Rational(r)
                # which points realise it
                where = [(Q, lab, DP, DNP) for (val, Q, lab, DP, DNP) in pts if abs(val - mp.mpf(r.p) / r.q) < mp.mpf(10) ** -100]
                found.setdefault(str(r), []).extend([dict(Dtau=DP, DNtau=DNP, form=Q, coset=lab, Dlist=Dp) for (Q, lab, DP, DNP) in where])
        results[N] = found
        print(f'\nN = {N}: worst |coefficient - nearest integer| over all {len(allD)} polynomials = {mp.nstr(worst, 5)}')
        for r in sorted(found, key=lambda s: sp.Rational(s)):
            ds = sorted({(w['Dtau'], w['DNtau']) for w in found[r]})
            print(f'  m_{N} = {r:>8}: (D(tau), D(N tau)) over its points = {ds}; #points = {len(found[r])}')
    with open(os.path.join(HERE, 'a3_cm_classification.json'), 'w') as f:
        json.dump({str(N): results[N] for N in results}, f, indent=1, default=str)
    claimed = {7: {-10648: -427, -64: -91, -8: -35, -1: -7, 0: -3, 2: -12, 8: -19, 24: -27, 27: -28, 125: -112},
               5: {-15228: -235, -828: -115, -28: -35, -3: -15, 0: -4, 4: -11, 18: -16, 36: -19, 72: -40, 147: -60}}
    for N in (5, 7):
        got = {int(sp.Rational(r)): sorted({w['DNtau'] for w in results[N][r]} | {w['Dtau'] for w in results[N][r]})
               for r in results[N]}
        ok_vals = set(got) == set(claimed[N])
        ok_D = all(got[m] == [claimed[N][m]] for m in claimed[N] if m in got)
        print(f'N={N}: value list matches THM-4708: {ok_vals}; D labels match (single D at both Fricke points): {ok_D}')
        if not ok_vals:
            print('   extra:', set(got) - set(claimed[N]), ' missing:', set(claimed[N]) - set(got))


if __name__ == '__main__':
    main()
