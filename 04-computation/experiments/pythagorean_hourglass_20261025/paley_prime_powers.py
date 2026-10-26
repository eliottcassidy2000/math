# Prime-power test of the five- and six-vertex Paley census formulas.  For q = p^r = 3 (mod 4) the Paley tournament
# lives on F_q; the formulas are tested with Frobenius power traces:
#   A_q = alpha^r + beta^r            (alpha + beta = a_p(f8), alpha beta = p^3; f8 = eta(2z)^4 eta(4z)^4)
#   H_q = (u^r + v^r)^2               (u + v = h(p), uv = p^2, b(p) = i h(p) for LMFDB 32.3.c.a, so H_p = |b(p)|^2)
#   TT5 = q(q-1)(q^3 - 24q^2 + 168q - 319 - A_q)/1024
#   TT6 = q(q-1)(q^4 - 49q^3 + 828q^2 - 5313q + 10193 - 7(q-15)A_q - 8H_q)/32768
# Also re-checks both formulas at the small primes 3, 7, 11 (where TT6 = 0).  mac-mini-2026-10-25-pyth.
import itertools, sys
import numpy as np

def irreducible(p, r):
    """monic f (coefficients low -> high, length r+1) irreducible over F_p, r prime: no root and x^(p^r) = x mod f"""
    for tail in itertools.product(range(p), repeat=r):
        f = list(tail) + [1]
        if f[0] == 0: continue
        if any(sum(c * pow(x, i, p) for i, c in enumerate(f)) % p == 0 for x in range(p)): continue
        y = [0, 1] + [0] * (r - 2)                      # x
        for _ in range(r): y = polypowmod(y, p, f, p)
        if y == [0, 1] + [0] * (r - 2): return f
    raise ValueError

def polypowmod(a, e, f, p):
    res = [1] + [0] * (len(f) - 2)
    while e:
        if e & 1: res = polymulmod(res, a, f, p)
        a = polymulmod(a, a, f, p); e >>= 1
    return res

def polymulmod(a, b, f, p):
    r = len(f) - 1
    c = [0] * (2 * r - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % p
    c = c + [0] * max(0, r - len(c))
    for k in range(len(c) - 1, r - 1, -1):          # reduce by monic f
        t = c[k]
        if t:
            for i in range(r + 1):
                c[k - r + i] = (c[k - r + i] - t * f[i]) % p
    return c[:r]

def paley_qr_adjacency(p, r):
    q = p ** r
    if r == 1:
        sq = sorted({(x * x) % p for x in range(1, p)})
        Q = np.array(sq)
        return q, np.isin((Q[None, :] - Q[:, None]) % p, sq).astype(np.int64)
    f = irreducible(p, r)
    digits = np.array([[(x // p ** i) % p for i in range(r)] for x in range(q)])
    def toint(c): return sum(int(ci) * p ** i for i, ci in enumerate(c))
    sq = sorted({toint(polymulmod(list(digits[x]), list(digits[x]), f, p)) for x in range(1, q)})
    assert len(sq) == (q - 1) // 2
    isq = np.zeros(q, dtype=bool); isq[sq] = True
    Qd = digits[sq]                                         # m x r
    diff = (Qd[None, :, :] - Qd[:, None, :]) % p            # m x m x r : Q[j] - Q[i]
    idx = (diff * (p ** np.arange(r))[None, None, :]).sum(axis=2)
    A = isq[idx].astype(np.int64)
    assert (A + A.T + np.eye(len(sq), dtype=np.int64) == 1).all(), 'not a tournament: q must be 3 mod 4'
    return q, A

def census(A):
    """(TT4, TT5) of the tournament with adjacency A"""
    m = len(A)
    M = A @ A.T
    tt4 = 0
    for w in range(m):
        Y = np.nonzero(A[w])[0]
        d = M[Y, w]; tt4 += int((d * (d - 1) // 2).sum())
    tt5 = 0
    for w in range(m):
        Y = np.nonzero(A[w])[0]
        AY = A[np.ix_(Y, Y)]; MY = AY @ AY.T
        for u in range(len(Y)):
            Z = np.nonzero(AY[u])[0]
            d = MY[Z, u]; tt5 += int((d * (d - 1) // 2).sum())
    return tt4, tt5

a_f8 = {3: -4, 7: 24, 11: -44}
h_b = {3: 4, 7: -8}
def power_sum(s1, prod, r):
    s = [2, s1]
    for k in range(2, r + 1): s.append(s1 * s[-1] - prod * s[-2])
    return s[r]

# h(11) from the eta-theta expansion
sys.argv = [sys.argv[0]]
M = 40
def ser_mul(a, b):
    c = [0] * (M + 1)
    for i, x in enumerate(a):
        for j in range(M + 1 - i): c[i + j] += x * b[j]
    return c
def eta_prod(spec):
    shift = sum(s * e for s, e in spec) // 24
    c = [0] * (M + 1); c[0] = 1
    for s, e in spec:
        for _ in range(e):
            for n in range(1, M // s + 1):
                for k in range(M, s * n - 1, -1): c[k] -= c[k - s * n]
    return [0] * shift + c[:M + 1 - shift]
def theta(sc):
    c = [0] * (M + 1); n = 0
    while sc * n * n <= M: c[sc * n * n] += 1 if n == 0 else 2; n += 1
    return c
e48 = eta_prod([(4, 2), (8, 2)]); t2 = theta(2); t4 = theta(4)
imb = ser_mul(e48, [x - y for x, y in zip(ser_mul(t2, t2), ser_mul(t4, t4))])
f8 = eta_prod([(2, 4), (4, 4)])
for p in (3, 7, 11): assert f8[p] == a_f8.get(p, f8[p]); a_f8[p] = f8[p]; h_b[p] = imb[p]
print("a_p(f8):", a_f8, " h(p) (b(p) = i h):", h_b)

for (p, r) in [(3, 1), (7, 1), (11, 1), (3, 3), (3, 5), (7, 3), (11, 3), (3, 7)]:
    q, A = paley_qr_adjacency(p, r)
    tt5_qr, tt6_qr = census(A)
    TT5, TT6 = q * tt5_qr, q * tt6_qr
    Aq = power_sum(a_f8[p], p ** 3, r)
    Hq = power_sum(h_b[p], p ** 2, r) ** 2
    F5 = q * (q - 1) * (q ** 3 - 24 * q ** 2 + 168 * q - 319 - Aq)
    F6 = q * (q - 1) * (q ** 4 - 49 * q ** 3 + 828 * q ** 2 - 5313 * q + 10193 - 7 * (q - 15) * Aq - 8 * Hq)
    ok5 = F5 == 1024 * TT5; ok6 = F6 == 32768 * TT6
    print(f"q = {p}^{r} = {q}: TT5 = {TT5} ({'OK' if ok5 else 'FAIL, formula ' + str(F5 / 1024)}), "
          f"TT6 = {TT6} ({'OK' if ok6 else 'FAIL, formula ' + str(F6 / 32768)}); A_q = {Aq}, H_q = {Hq}", flush=True)
