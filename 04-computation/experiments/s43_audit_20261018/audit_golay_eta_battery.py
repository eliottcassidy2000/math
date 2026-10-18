#!/usr/bin/env python3
"""S43 referee: independent checks of the Golay/M24 construction, the eta-product table, the Apery-like battery
definitions and the fractional-period battery (2026-10-18). No session code is imported.

Parts:
  G. extended QR code from translates of {0} u QR (mod 23) plus the all-ones word: dimension, weight enumerator;
     Conway's generators preserve it; group order by Schreier-Sims (sympy); Aut(P_23) cycle shapes.
  H. weight-2 eta-products vs elliptic curves 11a1, 14a1, 15a1, 20a1, 24a1, 36a1 (a_p = p + 1 - #E(F_p)).
  I. sequence definitions of the battery: recurrences, Zagier D vs Apery b, Cooper s18 (OSS 2016 formula) vs the session's
     formula, Yang-Zudilin; half-period tests (mod p) of the correctly defined zeta(3)-type sequences against +-a_4A(p).
  J. fractional period n == -1/3 (mod p): s7(n) vs a(p) of eta(3tau)^8 modulo p, primes 7 <= p < PMAX; and the
     'mechanism' congruences (mod p, mod p^2).
Usage: python audit_golay_eta_battery.py [G H I J]
"""
import sys, time
from math import comb, factorial
from fractions import Fraction as Fr

def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]

# ---------------- G ----------------
def partG():
    P = 23
    INF = 23
    QR = sorted({x * x % P for x in range(1, P)})
    words = []
    for a in range(P):
        w = 0
        for r in [0] + QR:
            w |= 1 << ((a + r) % P)
        words.append(w)
    words.append((1 << 24) - 1)
    # Gaussian elimination over F2 with pivots
    basis = []
    for w in words:
        for b in basis:
            if w ^ b < w:
                w ^= b
        if w:
            basis.append(w)
            basis.sort(reverse=True)
    dim = len(basis)
    code = [0]
    for b in basis:
        code += [c ^ b for c in code]
    we = {}
    for c in code:
        k = bin(c).count('1')
        we[k] = we.get(k, 0) + 1
    print('[G] dimension %d; weight enumerator %s' % (dim, dict(sorted(we.items()))), flush=True)
    codeset = set(code)
    inv = lambda x: pow(x, P - 2, P)
    def perm(f):
        return [f(x) for x in range(24)]
    alpha = perm(lambda x: INF if x == INF else (x + 1) % P)
    beta = perm(lambda x: INF if x == INF else 2 * x % P)
    gamma = perm(lambda x: 0 if x == INF else (INF if x == 0 else (-inv(x)) % P))
    delta = perm(lambda x: x if x in (0, INF) else (pow(x, 3, P) * inv(9) % P if x in QR else 9 * pow(x, 3, P) % P))
    def apply(pm, c):
        out = 0
        for i in range(24):
            if c >> i & 1:
                out |= 1 << pm[i]
        return out
    pres = {nm: all(apply(pm, b) in codeset for b in basis) for nm, pm in (('alpha', alpha), ('beta', beta), ('gamma', gamma), ('delta', delta))}
    print('[G] generators map a basis into the code (hence preserve it): %s' % pres, flush=True)
    try:
        from sympy.combinatorics import Permutation, PermutationGroup
        G = PermutationGroup([Permutation(alpha), Permutation(beta), Permutation(gamma), Permutation(delta)])
        print('[G] |<alpha,beta,gamma,delta>| = %d (|M24| = 244823040); transitivity degree check: %s' % (G.order(), G.is_transitive()), flush=True)
    except Exception as e:
        print('[G] sympy group order failed:', e, flush=True)
    def shape(pm):
        seen = [0] * 24
        sh = {}
        for i in range(24):
            if not seen[i]:
                L = 0
                j = i
                while not seen[j]:
                    seen[j] = 1
                    j = pm[j]
                    L += 1
                sh[L] = sh.get(L, 0) + 1
        return tuple(sorted(sh.items()))
    cnt = {}
    for a in QR:
        for b in range(P):
            pm = perm(lambda x, a=a, b=b: INF if x == INF else (a * x + b) % P)
            assert all(apply(pm, w) in codeset for w in basis)
            s = shape(pm)
            cnt[s] = cnt.get(s, 0) + 1
    print('[G] Aut(P_23) = {x -> a x + b : a in QR}: cycle shapes %s' % cnt, flush=True)
    # Paley tournament automorphism group check: maps preserving the arc relation j - i in QR
    QRs = set(QR)
    auts = 0
    for a in range(1, P):
        for b in range(P):
            f = lambda x: (a * x + b) % P
            if all(((f(j) - f(i)) % P in QRs) == ((j - i) % P in QRs) for i in range(P) for j in range(P) if i != j):
                auts += 1
    print('[G] affine maps x -> ax+b preserving P_23: %d (all automorphisms of P_p are affine: x -> a x + b, a in QR)' % auts, flush=True)

# ---------------- H ----------------
def eta_product(shape, N):
    """coefficients of prod_a prod_n (1 - q^(a n))^(m_a), up to q^N"""
    c = [0] * (N + 1)
    c[0] = 1
    for a, m in shape.items():
        for _ in range(m):
            for n in range(1, N // a + 1):
                st = a * n
                for k in range(N, st - 1, -1):
                    c[k] -= c[k - st]
    return c

def count_points(coeffs, p):
    """#E(F_p) for y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6 (projective)"""
    a1, a2, a3, a4, a6 = coeffs
    cnt = 1
    for x in range(p):
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0:
                cnt += 1
    return cnt

def partH():
    curves = {  # Cremona optimal curves [a1,a2,a3,a4,a6]
        '11a1': ((1, 0, 0, 0, 0, 0), {1: 2, 11: 2}, [0, -1, 1, -10, -20], 11),
        '14a1': (None, {1: 1, 2: 1, 7: 1, 14: 1}, [1, 0, 1, 4, -6], 14),
        '15a1': (None, {1: 1, 3: 1, 5: 1, 15: 1}, [1, 1, 1, -10, -10], 15),
        '20a1': (None, {2: 2, 10: 2}, [0, 1, 0, 4, 4], 20),
        '24a1': (None, {2: 1, 4: 1, 6: 1, 12: 1}, [0, -1, 0, -4, 4], 24),
        '36a1': (None, {6: 4}, [0, 0, 0, 0, 1], 36),
    }
    N = 200
    for name, (_, shape, ainv, cond) in curves.items():
        c = eta_product(shape, N)
        off = sum(a * m for a, m in shape.items()) // 24
        coef = lambda n: c[n - off] if n - off >= 0 else 0
        bad = []
        for p in primes_upto(N):
            if cond % p == 0:
                continue
            ap = p + 1 - count_points(ainv, p)
            if ap != coef(p):
                bad.append((p, ap, coef(p)))
        print('[H] %s vs eta-product %s: mismatches at good primes < %d: %s' % (name, shape, N, bad), flush=True)

# ---------------- I ----------------
SEQS = {
    'Apery A': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)),
    'Apery b (Zagier D)': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1)),
    'session Zagier D = sum C(n,k)C(2k,k)C(2n-2k,n-k)': lambda n: sum(comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
    'session Zagier E = sum C(n,2k)C(2k,k)^2 4^(n-2k)': lambda n: sum(comb(n, 2 * k) * comb(2 * k, k) ** 2 * 4 ** (n - 2 * k) for k in range(n // 2 + 1)),
    'Cooper s18 (OSS 2016 / Cooper 2012)': lambda n: 1 if n == 0 else sum((-1) ** k * comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) * ((comb(2 * n - 3 * k - 1, n) if 2 * n - 3 * k - 1 >= 0 else 0) + comb(2 * n - 3 * k, n)) for k in range(n // 3 + 1)),
    'session s18 = sum (-1)^k C(n,k)C(2k,k)C(2n-2k,n-k)(C(2k,n)+C(2n-2k,n))': lambda n: sum((-1) ** k * comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) * (comb(2 * k, n) + comb(2 * n - 2 * k, n)) for k in range(n + 1)),
    'Cooper s10 = Yang-Zudilin = sum C(n,k)^4': lambda n: sum(comb(n, k) ** 4 for k in range(n + 1)),
    'session Yang-Zudilin = sum C(n,k)^2 C(2k,k)^2': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) ** 2 for k in range(n + 1)),
    'Cooper s7': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range(n + 1)),
    'Domb': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
    'Almkvist-Zudilin (delta)': lambda n: sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * comb(n + k, n) * factorial(3 * k) // factorial(k) ** 3 for k in range(n // 3 + 1)),
    '(epsilon) sum C(n,k)^2 C(2k,n)^2': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, n) ** 2 for k in range(n + 1)),
    'Franel': lambda n: sum(comb(n, k) ** 3 for k in range(n + 1)),
    'Zagier C': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) for k in range(n + 1)),
    'Zagier B': lambda n: sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * factorial(3 * k) // factorial(k) ** 3 for k in range(n // 3 + 1)),
    'Zagier F': lambda n: sum((-1) ** k * 8 ** (n - k) * comb(n, k) * sum(comb(k, j) ** 3 for j in range(k + 1)) for k in range(n + 1)),
}

def fit_recurrence(u, order2=True):
    """Find (a,b,c[,d]) with (n+1)^r u(n+1) = poly u(n) - poly u(n-1); try Zagier 2nd-order and Cooper 3rd-order forms."""
    import sympy as sp
    a, b, c, d = sp.symbols('a b c d')
    eqs2, eqs3 = [], []
    for n in range(1, 12):
        eqs2.append(sp.Eq((n + 1) ** 2 * u[n + 1], (a * n * n + a * n + b) * u[n] - c * n * n * u[n - 1]))
        eqs3.append(sp.Eq((n + 1) ** 3 * u[n + 1], (2 * n + 1) * (a * n * n + a * n + b) * u[n] - n * (c * n * n + d) * u[n - 1]))
    s2 = sp.solve(eqs2[:3], [a, b, c], dict=True)
    s3 = sp.solve(eqs3[:4], [a, b, c, d], dict=True)
    out = []
    for s in s2:
        if all(sp.simplify(e.lhs - e.rhs).subs(s) == 0 for e in eqs2):
            out.append(('Zagier (a,b,c)', tuple(s[x] for x in (a, b, c))))
    for s in s3:
        if all(sp.simplify(e.lhs - e.rhs).subs(s) == 0 for e in eqs3):
            out.append(('Cooper (a,b,c,d)', tuple(s[x] for x in (a, b, c, d))))
    return out

def partI(PMAX=600):
    for name, f in SEQS.items():
        u = [f(n) for n in range(14)]
        print('[I] %-72s first terms %s; recurrence %s' % (name, u[:7], fit_recurrence(u) or 'none of the two shapes'), flush=True)
    same = all(SEQS['session Zagier D = sum C(n,k)C(2k,k)C(2n-2k,n-k)'](n) == SEQS['session Zagier E = sum C(n,2k)C(2k,k)^2 4^(n-2k)'](n) for n in range(40))
    s18same = all(SEQS['Cooper s18 (OSS 2016 / Cooper 2012)'](n) == SEQS['session s18 = sum (-1)^k C(n,k)C(2k,k)C(2n-2k,n-k)(C(2k,n)+C(2n-2k,n))'](n) for n in range(40))
    print('[I] session "Zagier D" == session "Zagier E" for n < 40: %s; session s18 == OSS s18 for n < 40: %s' % (same, s18same), flush=True)
    # half-period test mod p for the zeta(3)-type sequences against +-a_4A(p)
    N = PMAX + 2
    c4A = eta_product({2: 4, 4: 4}, N)
    a4A = lambda m: c4A[m - 1] if m >= 1 else 0
    ps = [p for p in primes_upto(PMAX) if p >= 5]
    for name in ('Apery A', 'Cooper s7', 'Domb', 'Cooper s10 = Yang-Zudilin = sum C(n,k)^4', 'Almkvist-Zudilin (delta)', '(epsilon) sum C(n,k)^2 C(2k,n)^2', 'Cooper s18 (OSS 2016 / Cooper 2012)', 'session s18 = sum (-1)^k C(n,k)C(2k,k)C(2n-2k,n-k)(C(2k,n)+C(2n-2k,n))'):
        f = SEQS[name]
        vals = {p: f((p - 1) // 2) for p in ps}
        res = []
        for eps_name, eps in (('+1', lambda p: 1), ('(-1/p)', lambda p: 1 if p % 4 == 1 else -1)):
            m1 = all((vals[p] - eps(p) * a4A(p)) % p == 0 for p in ps)
            m2 = all((vals[p] - eps(p) * a4A(p)) % (p * p) == 0 for p in ps)
            res.append('%s: mod p %s, mod p^2 %s' % (eps_name, m1, m2))
        print('[I] %-60s at (p-1)/2 vs eps*a_4A(p), primes 5..%d: %s' % (name, PMAX, '; '.join(res)), flush=True)

# ---------------- J ----------------
def partJ(PMAX=1500):
    t = time.time()
    N = PMAX + 2
    c3B = eta_product({3: 8}, N)       # eta(3tau)^8 = q prod (1 - q^(3n))^8
    a3B = lambda m: c3B[m - 1] if m >= 1 else 0
    ps = [p for p in primes_upto(PMAX) if p >= 7]
    s7 = lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range((n + 1) // 2, n + 1))
    bad1, bad2 = [], []
    zero2 = []
    for p in ps:
        n = (-pow(3, -1, p)) % p     # least n >= 0 with n == -1/3 mod p
        v = s7(n)
        if (v - a3B(p)) % p:
            bad1.append(p)
        if (v - a3B(p)) % (p * p) == 0:
            zero2.append(p)
    print('[J] s7(n), n == -1/3 mod p (least n), vs a(p) of eta(3tau)^8 mod p, primes 7 <= p < %d: %d tested, failures %s; '
          'number also holding mod p^2: %d (%.1fs)' % (PMAX, len(ps), bad1, len(zero2), time.time() - t), flush=True)
    # mechanism checks at s = 1/3 for p < 300
    s = Fr(1, 3)
    res = {}
    for test in ('A vs sum b^2', 's7 vs conv'):
        for mod in (1, 2):
            ok = True
            for p in [q for q in ps if q < 300]:
                n = (-pow(3, -1, p)) % p
                m = p ** mod
                bk = [Fr(1)]
                for k in range(1, n + 1):
                    bk.append(bk[-1] * (s + k - 1) * (1 - s + k - 1) / (k * k))
                if test == 'A vs sum b^2':
                    lhs = sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))
                    r = sum(x * x for x in bk)
                else:
                    lhs = s7(n)
                    r = (-1) ** n * sum(bk[k] * bk[n - k] for k in range(n + 1))
                if r.denominator % p == 0:
                    ok = None
                    break
                if (lhs - r.numerator * pow(r.denominator, -1, m)) % m:
                    ok = False
                    break
            res['%s mod p^%d' % (test, mod)] = ok
    print('[J] mechanism at s = 1/3 (primes 7..300):', res, flush=True)

if __name__ == '__main__':
    for q in (sys.argv[1:] or ['G', 'H', 'I', 'J']):
        {'G': partG, 'H': partH, 'I': partI, 'J': partJ}[q]()
