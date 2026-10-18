#!/usr/bin/env python3
"""S43 referee: independent checks of section 4 of factorization_vertex_supersymmetric_bridge_20261018.md (2026-10-18).

D. SU(2) N_f = 4 Schur index by numerical contour integration (trapezoidal rule on |z| = 1, mpmath q-Pochhammers) at a
   tiny q, coefficients peeled off; then the MLDE D^(2) f = alpha E_4 f, f = q^(7/12) I, alpha = 35/144 (exact).
E. Lee-Yang characters from the Rogers-Ramanujan SUM sides sum q^(n^2+n)/(q)_n, sum q^(n^2)/(q)_n; MLDE alpha = 11/3600.
F. Beukers' level-5 identity: sum b_n t^n with t = q prod (1-q^n)^(5 (n/5)) vs 1 + sum Re[(3 -+ i) sum_{d|n} chi(d)] q^n,
   both choices chi(2) = +-i; and t = q H^5/G^5 = (chi_0/chi_{-1/5})^5.
Usage: python audit_vertex_4d.py [D E F]
"""
import sys
from fractions import Fraction as Fr
from math import comb

def sigma(n, k=1):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)

def mlde_residuals(c, h, alpha, N):
    E2 = [Fr(1)] + [Fr(-24 * sigma(n)) for n in range(1, N + 1)]
    E4 = [Fr(1)] + [Fr(240 * sigma(n, 3)) for n in range(1, N + 1)]
    res = []
    for n in range(N + 1):
        t = (h + n) ** 2 * c[n]
        t -= Fr(1, 6) * sum(E2[k] * (h + n - k) * c[n - k] for k in range(n + 1))
        t -= alpha * sum(E4[k] * c[n - k] for k in range(n + 1))
        res.append(t)
    return res

def partD(NC=12):
    import mpmath as mp
    mp.mp.dps = 170
    q = mp.mpf(10) ** -12
    sq = mp.sqrt(q)
    M = 40
    tot = mp.mpf(0)
    qq = mp.qp(q, q)
    for j in range(M):
        z = mp.exp(2j * mp.pi * j / M)
        num = (1 - z ** 2) * (1 - z ** -2) / 2 * qq ** 2 * mp.qp(q * z ** 2, q) ** 2 * mp.qp(q * z ** -2, q) ** 2
        den = mp.qp(sq * z, q) ** 8 * mp.qp(sq / z, q) ** 8
        tot += num / den
    I = tot / M
    coeffs = []
    rem = I
    halfint = []
    for n in range(NC):
        cn = rem / q ** n
        k = int(mp.nint(mp.re(cn)))
        coeffs.append(k)
        rem -= k * q ** n
        # size of what is left at order q^(n+1/2): should be ~ q^(n+1) (no half-integer powers)
        halfint.append(mp.nstr(abs(rem) / q ** (n + mp.mpf(1) / 2), 3))
    print('[D] SU(2) N_f=4 Schur index (contour integral at q = 1e-12): I = %s' % coeffs, flush=True)
    print('[D] |remainder|/q^(n+1/2) after peeling q^n (half-integer powers would give O(1)): %s' % halfint[:8], flush=True)
    h = Fr(7, 12)
    alpha = h * (h - Fr(1, 6))
    res = mlde_residuals([Fr(x) for x in coeffs], h, alpha, NC - 1)
    print('[D] MLDE with alpha = %s (h = 7/12, c = -14): residuals zero for n <= %d: %s' % (alpha, NC - 1, all(r == 0 for r in res)), flush=True)

def qpoch_series(n, N):
    """1/(q;q)_n as a series to q^N"""
    c = [Fr(0)] * (N + 1)
    c[0] = Fr(1)
    for k in range(1, n + 1):
        for m in range(k, N + 1):
            c[m] += c[m - k]
    return c

def partE(N=60):
    H = [Fr(0)] * (N + 1)
    G = [Fr(0)] * (N + 1)
    n = 0
    while n * n <= N:
        s = qpoch_series(n, N)
        for m in range(N + 1):
            if m + n * n <= N:
                G[m + n * n] += s[m]
            if m + n * n + n <= N:
                H[m + n * n + n] += s[m]
        n += 1
    # product sides
    def prod(res):
        c = [Fr(0)] * (N + 1)
        c[0] = Fr(1)
        for k in range(1, N + 1):
            if k % 5 in res:
                for m in range(k, N + 1):
                    c[m] += c[m - k]
        return c
    rr_ok = H == prod((2, 3)) and G == prod((1, 4))
    a = Fr(11, 3600)
    okH = all(r == 0 for r in mlde_residuals(H, Fr(11, 60), a, N))
    okG = all(r == 0 for r in mlde_residuals(G, Fr(-1, 60), a, N))
    print('[E] Rogers-Ramanujan sum = product to q^%d: %s; H = %s...; vacuum q^(11/60)H solves MLDE(11/3600): %s; q^(-1/60)G: %s; '
          'H has no q^1 term: %s' % (N, rr_ok, [int(x) for x in H[:10]], okH, okG, H[1] == 0), flush=True)
    return H, G

def partF(N=60):
    H, G = partE(N)
    # t = q prod (1 - q^n)^(5 (n/5))
    def leg5(n):
        r = n % 5
        return 0 if r == 0 else (1 if r in (1, 4) else -1)
    base = [Fr(0)] * (N + 1)
    base[0] = Fr(1)
    for n in range(1, N + 1):
        e = 5 * leg5(n)
        for _ in range(abs(e)):
            if e > 0:
                for m in range(N, n - 1, -1):
                    base[m] -= base[m - n]
            else:
                for m in range(n, N + 1):
                    base[m] += base[m - n]
    t = [Fr(0)] + base[:N]
    # (chi_0/chi_{-1/5})^5 = q (H/G)^5
    def mul(a, b):
        out = [Fr(0)] * (N + 1)
        for i in range(N + 1):
            if a[i]:
                for j in range(N + 1 - i):
                    out[i + j] += a[i] * b[j]
        return out
    def inv(a):
        r = [Fr(0)] * (N + 1)
        r[0] = 1 / a[0]
        for n in range(1, N + 1):
            r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
        return r
    ratio = mul(H, inv(G))
    r5 = [Fr(1)] + [Fr(0)] * N
    for _ in range(5):
        r5 = mul(r5, ratio)
    t2 = [Fr(0)] + r5[:N]
    b = [sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1)) for n in range(N + 1)]
    F = [Fr(0)] * (N + 1)
    tp = [Fr(1)] + [Fr(0)] * N
    for n in range(N + 1):
        for i in range(N + 1):
            F[i] += b[n] * tp[i]
        tp = mul(tp, t)
    out = {}
    for lab, chi2 in (('chi(2)=+i', 1j), ('chi(2)=-i', -1j)):
        chi = {0: 0, 1: 1, 2: chi2, 3: chi2.conjugate(), 4: -1}
        ok = all(abs(((3 - 1j) * sum(chi[d % 5] for d in range(1, n + 1) if n % d == 0)).real - float(F[n])) < 1e-9 for n in range(1, N + 1))
        out[lab] = ok
    print('[F] t = q prod(1-q^n)^(5(n/5)) equals q (H/G)^5 to q^%d: %s; sum b_n t^n = %s ...; matches 1 + sum Re[(3-i) sum chi(d)] q^n: %s'
          % (N, t == t2, [int(x) for x in F[:12]], out), flush=True)

if __name__ == '__main__':
    for qq in (sys.argv[1:] or ['D', 'F']):
        {'D': partD, 'E': partE, 'F': partF}[qq]()
