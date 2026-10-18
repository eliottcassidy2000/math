#!/usr/bin/env python3
"""Audit F2, check (6): Fibonacci started at (5,1) and (7,1) -- independent recomputation of every claim of section 3
of first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md (3.1 - 3.3)."""
from fractions import Fraction
from decimal import Decimal, getcontext
from math import isqrt, comb
from sympy import isprime

getcontext().prec = 80


def gib(a, b, lo, hi):
    """dict n -> G_n for lo <= n <= hi (negative indices by running the recurrence backwards)."""
    G = {0: a, 1: b}
    for n in range(2, hi + 1): G[n] = G[n - 1] + G[n - 2]
    for n in range(-1, lo - 1, -1): G[n] = G[n + 2] - G[n + 1]
    return G


F = gib(0, 1, -40, 3000)
L = gib(2, 1, -40, 3000)
NOTE = {(5, 1): [5, 1, 6, 7, 13, 20, 33, 53, 86, 139, 225, 364, 589, 953, 1542, 2495, 4037, 6532, 10569, 17101, 27670,
                 44771, 72441, 117212],
        (7, 1): [7, 1, 8, 9, 17, 26, 43, 69, 112, 181, 293, 474, 767, 1241, 2008, 3249, 5257, 8506, 13763, 22269, 36032,
                 58301, 94333, 152634]}
# OEIS data (read from oeis.org/A022096 and /A022098, offset 0)
OEIS = {(5, 1): [1, 6, 7, 13, 20, 33, 53, 86, 139, 225, 364, 589, 953, 1542, 2495, 4037, 6532, 10569, 17101, 27670, 44771,
                 72441, 117212, 189653, 306865, 496518, 803383, 1299901, 2103284, 3403185],
        (7, 1): [1, 8, 9, 17, 26, 43, 69, 112, 181, 293, 474, 767, 1241, 2008, 3249, 5257, 8506, 13763, 22269, 36032, 58301,
                 94333, 152634, 246967, 399601, 646568, 1046169, 1692737, 2738906, 4431643]}


# Z[phi] arithmetic: (x, y) = x + y phi, phi^2 = phi + 1
def zm(p, q): return (p[0] * q[0] + p[1] * q[1], p[0] * q[1] + p[1] * q[0] + p[1] * q[1])
def zpow(p, k):
    if k < 0:
        p = zinv(p); k = -k
    r = (1, 0)
    for _ in range(k): r = zm(r, p)
    return r
def znorm(p): return p[0] * p[0] + p[0] * p[1] - p[1] * p[1]
def zconj(p): return (p[0] + p[1], -p[1])      # phi -> 1 - phi
def zinv(p):
    n = znorm(p); c = zconj(p)
    assert c[0] % n == 0 and c[1] % n == 0, 'not a unit'
    return (c[0] // n, c[1] // n)
PHI = (0, 1); SQ5 = (-1, 2)                   # sqrt5 = 2 phi - 1

out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

for (a, b) in [(5, 1), (7, 1)]:
    G = gib(a, b, -40, 3000)
    Nab = b * b - a * b - a * a
    say('== (%d,%d): N = b^2 - ab - a^2 = %d' % (a, b, Nab))
    say('  listed 24 terms correct:', [G[n] for n in range(24)] == NOTE[(a, b)])
    say('  = %d followed by the OEIS terms (30 terms):' % a, [G[n] for n in range(1, 31)] == OEIS[(a, b)])
    say('  G_n = a F_(n-1) + b F_n for -30 <= n <= 2000:', all(G[n] == a * F[n - 1] + b * F[n] for n in range(-30, 2001)))
    say('  G_n = (a L_n + (2b - a) F_n)/2:', all(2 * G[n] == a * L[n] + (2 * b - a) * F[n] for n in range(-30, 2001)))
    if (a, b) == (5, 1):
        say('  G_n = (5 L_n - 3 F_n)/2:', all(2 * G[n] == 5 * L[n] - 3 * F[n] for n in range(-30, 2001)))
        say('  G_n = F_(n+3) - F_(n-4) = F_(n-5) + ... + F_(n+1):',
            all(G[n] == F[n + 3] - F[n - 4] == sum(F[k] for k in range(n - 5, n + 2)) for n in range(-30, 2001)))
        say('  OEIS indexing check: A022096(n) = G_(n+1) = F(n+4) - F(n-3) [McRae sum F(n-4..n+2)]:',
            all(G[n + 1] == F[n + 4] - F[n - 3] == sum(F[k] for k in range(n - 4, n + 3)) for n in range(0, 40)),
            '; the OEIS formula "a(n) = F(n+3) - F(n-4)" (Dresden-Neale) evaluated at offset 0 gives',
            [F[n + 3] - F[n - 4] for n in range(4)], '(= our G_n, i.e. a(n-1))')
        theta = (b - a, a)                       # b + a(phi - 1)
        rhs = zm(zpow(PHI, -4), (zpow(PHI, 7)[0] - 1, zpow(PHI, 7)[1]))
        say('  theta = 5 phi - 4 = phi^-4 (phi^7 - 1):', theta == (-4, 5) == rhs)
    else:
        say('  G_n = 4 L_n - L_(n+1) = L_(n+2) - L_(n-3) = L_(n-4) + ... + L_n:',
            all(G[n] == 4 * L[n] - L[n + 1] == L[n + 2] - L[n - 3] == sum(L[k] for k in range(n - 4, n + 1)) for n in range(-30, 2001)))
        say('  OEIS indexing check: A022098(n) = G_(n+1) = L(n+3) - L(n-2) [Dresden-Nyc], = L(-3+n) + ... + L(n+1) [del Arte]:',
            all(G[n + 1] == L[n + 3] - L[n - 2] == sum(L[k] for k in range(n - 3, n + 2)) for n in range(0, 40)))
        theta = (b - a, a)
        rhs = zm(zm(zpow(PHI, -3), SQ5), (zpow(PHI, 5)[0] - 1, zpow(PHI, 5)[1]))
        say('  theta = 7 phi - 6 = phi^-3 sqrt5 (phi^5 - 1):', theta == (-6, 7) == rhs)
    say('  N(theta) = %d = Cassini constant:' % znorm(theta), znorm(theta) == Nab)
    say('  orbit: theta phi^n = G_n phi + G_(n-1) for -20 <= n <= 60:',
        all(zm(theta, zpow(PHI, n)) == (G[n - 1], G[n]) for n in range(-20, 61)))
    say('  Cassini G_n^2 - G_(n-1) G_(n+1) = (-1)^n * %d for -30 <= n <= 2000:' % (-Nab),
        all(G[n] ** 2 - G[n - 1] * G[n + 1] == (-1) ** n * (-Nab) for n in range(-30, 2001)))
    # growth
    s5 = Decimal(5).sqrt(); phi = (1 + s5) / 2; psi = (1 - s5) / 2
    lim = b + a / phi; alpha = lim / s5; alphac = (b + a / psi) / s5 if False else None
    # conjugate weight: G_n = alpha phi^n + alpha' psi^n with alpha' = conj(theta)/(-sqrt5)
    th = Decimal(theta[0]) + Decimal(theta[1]) * phi; thc = Decimal(theta[0]) + Decimal(theta[1]) * psi
    al = th / s5; alc = -thc / s5
    say('  G_n/F_n -> b + a/phi = %.4f ; alpha = %.4f, conj alpha = %.4f, alpha*conj(alpha) = %s (-N/5 = %s)'
        % (float(lim), float(al), float(alc), round(float(al * alc), 12), Fraction(-Nab, 5)))
    say('  alpha closed form:', '(25 - 3 sqrt5)/10' if (a, b) == (5, 1) else '(7 - sqrt5)/2', '=',
        float((25 - 3 * s5) / 10) if (a, b) == (5, 1) else float((7 - s5) / 2), '; matches:',
        abs(al - ((25 - 3 * s5) / 10 if (a, b) == (5, 1) else (7 - s5) / 2)) < Decimal(10) ** -60)
    rok = [n for n in range(0, 200) if int((al * phi ** n).to_integral_value()) == G[n]]
    first_all = min(n for n in range(200) if all(m in rok for m in range(n, 200)))
    say('  G_n = round(alpha phi^n) holds for all n >= %d (checked to 200; error |alpha\'| |psi|^n -> 0); fails at n = %s'
        % (first_all, [n for n in range(first_all) if n not in rok]))
    # primes dividing N: eigen-sequences
    for p in sorted({q for q in range(2, abs(Nab) + 1) if Nab % q == 0 and isprime(q)}):
        lam = b * pow(a, -1, p) % p
        eig = all(G[n] % p == a * pow(lam, n, p) % p for n in range(0, 3000))
        order = next(k for k in range(1, p) if pow(lam, k, p) == 1)
        say('  p = %d | N: lambda = %d (root of x^2-x-1: %s), eigen-sequence: %s, order of lambda %d, any term divisible: %s'
            % (p, lam, (lam * lam - lam - 1) % p == 0, eig, order, any(G[n] % p == 0 for n in range(0, 3000))))
    # squares
    sq = [(n, G[n]) for n in range(0, 3000) if G[n] >= 0 and isqrt(G[n]) ** 2 == G[n]]
    say('  squares G_n, 0 <= n < 3000:', sq)
    # Siegel identity
    say('  (2 G_(n+1) - G_n)^2 = 5 G_n^2 + 4 (-1)^n N for -30 <= n <= 2000:',
        all((2 * G[n + 1] - G[n]) ** 2 == 5 * G[n] ** 2 + 4 * (-1) ** n * Nab for n in range(-30, 2001)))
    # prime terms
    pt = [(n, G[n]) for n in range(0, 120) if isprime(G[n])]
    say('  prime terms (n < 120):', pt)

# Fibonacci squares control (Cohn): 0, 1, 144
say('Fibonacci squares n < 3000:', [(n, F[n]) for n in range(0, 3000) if isqrt(F[n]) ** 2 == F[n]])
# N(phi^r - 1) = -L_r for odd r, 2 - L_r for even r
ok = all(znorm((zpow(PHI, r)[0] - 1, zpow(PHI, r)[1])) == (-L[r] if r % 2 else 2 - L[r]) for r in range(1, 60))
say('N(phi^r - 1) = -L_r (r odd), 2 - L_r (r even), 1 <= r < 60:', ok, '; L_5 =', L[5], ', L_7 =', L[7])


def period(a, b, m):
    x, y = a % m, b % m; s = (x, y); k = 0
    while True:
        x, y = y, (x + y) % m; k += 1
        if (x, y) == s: return k
M = 3000
d51 = [m for m in range(1, M) if period(5, 1, m) != period(0, 1, m)]
d71L = [m for m in range(1, M) if period(7, 1, m) != period(2, 1, m)]
d71F = [m for m in range(1, M) if period(7, 1, m) != period(0, 1, m)]
say('periods, m < %d: (5,1) differs from Pisano at %d moduli %s ; (7,1) differs from Lucas at %d ; (7,1) differs from Pisano at %d moduli, all multiples of 5: %s'
    % (M, len(d51), d51[:10], len(d71L), len(d71F), all(m % 5 == 0 for m in d71F)))
say('period of (5,1) mod 29 = %d = pi(29) = %d ; mod 841: %d = pi(841) = %d' % (period(5, 1, 29), period(0, 1, 29), period(5, 1, 841), period(0, 1, 841)))
