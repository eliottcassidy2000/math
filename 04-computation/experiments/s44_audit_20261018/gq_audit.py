#!/usr/bin/env python3
"""S44 audit, section 6: G_q = Cay(F_(q^2), mu_(q+1)), q an odd prime.
(i)  lambda(a) = sum_{u in mu_(q+1)} psi(Tr(a u)) = -Kl_q(N(a)) checked for ALL a != 0 (q <= 31), exhaustively;
(ii) Hoffman chromatic bound 1 + (q+1)/max Kl at q = 23, 29, 31, 59 in 50-digit arithmetic;
(iii) arc-transitivity of x -> zeta x + b checked exhaustively for small q;
(iv) the Delsarte LP (= Schrijver theta' for this orbital scheme) for alpha(G_q): does it equal Hoffman's bound
     n(-lambda_min)/(lambda_max - lambda_min)?  (Lovasz Thm 9 gives theta = Hoffman; theta' <= theta may be smaller.)
"""
import math, cmath, json
import mpmath
from scipy.optimize import linprog


def nonsquare(q):
    for g in range(2, q):
        if pow(g, (q - 1) // 2, q) == q - 1:
            return g


def field(q):
    g = nonsquare(q)
    mul = lambda x, y: ((x[0] * y[0] + g * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
    norm = lambda x: (x[0] * x[0] - g * x[1] * x[1]) % q
    tr = lambda x: (2 * x[0]) % q
    return mul, norm, tr


def kloosterman(q, prec=None):
    inv = [0] + [pow(x, q - 2, q) for x in range(1, q)]
    if prec:
        mpmath.mp.dps = prec
        return {c: mpmath.fsum(mpmath.cos(2 * mpmath.pi * ((x + c * inv[x]) % q) / q) for x in range(1, q)) for c in range(1, q)}
    return {c: sum(math.cos(2 * math.pi * ((x + c * inv[x]) % q) / q) for x in range(1, q)) for c in range(1, q)}


def check_eigen(q):
    mul, norm, tr = field(q)
    circle = [(a, b) for a in range(q) for b in range(q) if norm((a, b)) == 1]
    assert len(circle) == q + 1
    e = [cmath.exp(2j * math.pi * t / q) for t in range(q)]
    Kl = kloosterman(q)
    worst = 0.0
    for a0 in range(q):
        for a1 in range(q):
            if (a0, a1) == (0, 0):
                continue
            lam = sum(e[tr(mul((a0, a1), u))] for u in circle)
            worst = max(worst, abs(lam - (-Kl[norm((a0, a1))])))
    return worst


def check_arc_transitive(q):
    mul, norm, tr = field(q)
    circle = [(a, b) for a in range(q) for b in range(q) if norm((a, b)) == 1]
    # stabilizer of 0 in the group {x -> zeta x + b} acts on the neighbours of 0 (= circle) by multiplication: transitive?
    orbit = {mul(z, circle[0]) for z in circle}
    return len(orbit) == len(circle)


def delsarte_lp(q):
    Kl = kloosterman(q)
    cs = [c for c in range(2, q)]  # classes c != 0, 1 (class 1 = edges, forced a_1 = 0)
    # maximize sum a_c  s.t.  -sum_c a_c (-Kl(c t)) <= 1  for all t in F_q^x, a_c >= 0
    # inner distribution a_c (a_0 = 1, a_1 = 0); Delsarte: 1 + sum_c a_c * (-Kl(c t))/(q+1) >= 0 for every t,
    # i.e. sum_c a_c Kl(ct) <= q + 1 (eigenvalue of class c at characters of norm t, divided by its valency q+1)
    A = [[Kl[(c * t) % q] for c in cs] for t in range(1, q)]
    b = [float(q + 1)] * (q - 1)
    res = linprog([-1.0] * len(cs), A_ub=A, b_ub=b, bounds=[(0, None)] * len(cs), method='highs')
    lp = 1 + (-res.fun)
    kmax = max(Kl.values())
    hoff_alpha = q * q * kmax / (q + 1 + kmax)
    return lp, hoff_alpha


def main():
    out = {}
    for q in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        err = check_eigen(q)
        at = check_arc_transitive(q)
        out[q] = dict(eig_err=err, arc_transitive_stabilizer=at)
        print('q=%2d: max |lambda(a) + Kl(N(a))| over all %d nonzero a = %.2e; stabilizer transitive on neighbours: %s' % (q, q * q - 1, err, at))
    for q in [23, 29, 31, 59]:
        Kl = kloosterman(q, prec=50)
        km = max(Kl.values())
        hb = 1 + (q + 1) / km
        out[q] = dict(out.get(q, {}), maxKl=mpmath.nstr(km, 20), hoffman=mpmath.nstr(hb, 20))
        print('q=%2d: max Kl = %s, Hoffman = %s, (q+1)/3 - maxKl = %s' % (q, mpmath.nstr(km, 15), mpmath.nstr(hb, 15), mpmath.nstr(mpmath.mpf(q + 1) / 3 - km, 5)))
    for q in [13, 17, 19, 23, 29, 31, 37]:
        lp, ha = delsarte_lp(q)
        print('q=%2d: Delsarte LP bound on alpha = %.6f ; Hoffman alpha bound = %.6f ; equal: %s' % (q, lp, ha, abs(lp - ha) < 1e-6))
        out.setdefault(q, {}).update(delsarte_lp=lp, hoffman_alpha=ha)
    json.dump({str(k): v for k, v in out.items()}, open('gq_audit.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
