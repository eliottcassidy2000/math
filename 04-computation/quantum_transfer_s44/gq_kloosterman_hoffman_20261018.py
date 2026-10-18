#!/usr/bin/env python3
"""Spectrum of the finite unit-distance graphs G_q = Cay(F_(q^2), mu_(q+1)) via Kloosterman sums (S44; builds on
mac-mini-2026-10-18-sixcolour's kappa table and Hoffman margins, merge_18_19_11_22_collatz_hn_20261018.md 2.2-2.2b).

Claims checked (q an odd prime):
 (i)  the eigenvalue at the additive character psi_a (a in F_(q^2)^x) is lambda(a) = sum_{u in mu_(q+1)} psi(Tr(a u))
      = -Kl_q(N(a)),  Kl_q(c) = sum_{x in F_q^x} e((x + c/x)/q)   [circle sum = minus Kloosterman sum];
 (ii) hence lambda_min = -max_c Kl_q(c), and Hoffman's ratio bound for the chromatic number is
      kappa(q) >= 1 + (q+1) / max_{c in F_q^x} Kl_q(c);
 (iii) G_q is arc-transitive (u -> zeta u, x -> x + b), so the Lovasz theta function equals the Hoffman ratio bound:
      no LP/SDP of Delsarte-Lovasz type can improve it (consistent with mac-mini's Delsarte/Schrijver LP at q = 29).
Prints the table of max Kl, Hoffman bound and its margin above the next-lower integer, for primes 5 <= q <= QMAX.
Usage: python gq_kloosterman_hoffman_20261018.py [QMAX]
Session opus-2026-10-18-S44.
"""
import sys, math, cmath, json, os


def primes_upto(n):
    return [p for p in range(3, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def nonsquare(q):
    for g in range(2, q):
        if pow(g, (q - 1) // 2, q) == q - 1:
            return g


def circle_eigen_check(q):
    """F_(q^2) = F_q[s]/(s^2 - g); elements (a, b) = a + b s.  Check lambda(a) = -Kl(N(a)) for a sample of a."""
    g = nonsquare(q)
    def mul(x, y):
        return ((x[0] * y[0] + g * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
    def norm(x):
        return (x[0] * x[0] - g * x[1] * x[1]) % q
    def tr(x):
        return (2 * x[0]) % q
    circle = [(a, b) for a in range(q) for b in range(q) if norm((a, b)) == 1]
    assert len(circle) == q + 1
    e = [cmath.exp(2j * math.pi * t / q) for t in range(q)]
    inv = [0] + [pow(x, q - 2, q) for x in range(1, q)]
    Kl = {c: sum(e[(x + c * inv[x]) % q] for x in range(1, q)).real for c in range(1, q)}
    worst = 0.0
    for a in [(1, 0), (2, 1), (0, 1), (3, 5 % q), (q - 1, 2)]:
        lam = sum(e[tr(mul(a, u))] for u in circle)
        worst = max(worst, abs(lam - (-Kl[norm(a)])))
    return worst, Kl


def main():
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rows = []
    maxerr = 0.0
    for q in primes_upto(QMAX):
        if q < 5:
            continue
        err, Kl = circle_eigen_check(q) if q <= 60 else (0.0, None)
        if Kl is None:
            e = [cmath.exp(2j * math.pi * t / q) for t in range(q)]
            inv = [0] + [pow(x, q - 2, q) for x in range(1, q)]
            Kl = {c: sum(e[(x + c * inv[x]) % q] for x in range(1, q)).real for c in range(1, q)}
        maxerr = max(maxerr, err)
        mk = max(Kl.values())
        hb = 1 + (q + 1) / mk
        rows.append(dict(q=q, maxKl=mk, weil=2 * math.sqrt(q), hoffman=hb, margin=hb - math.floor(hb)))
    print('(i) lambda(a) = -Kl_q(N(a)) checked for primes q <= 60 at 5 sample a each: max |error| = %.2e' % maxerr)
    print('(ii) Hoffman bound kappa(q) >= 1 + (q+1)/max Kl_q   [q, max Kl, 2 sqrt q, Hoffman, margin above floor]')
    for r in rows:
        flag = '  <-- near miss' if r['margin'] < 0.01 else ''
        print('  q=%4d  maxKl=%9.5f  2sqrt(q)=%8.4f  Hoffman=%.7f  margin=%.7f%s' % (r['q'], r['maxKl'], r['weil'], r['hoffman'], r['margin'], flag))
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gq_kloosterman_hoffman.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
