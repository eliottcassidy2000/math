"""Audit: do the 18 nonzero rational CM values of THM-4708 pass (B1') at every prime 7 <= p <= PMAX (p not | m)?
Sequences from their three-term recurrences (equal to the binomial definitions for n < 699, a2_hauptmoduln.py)."""
import sys
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
def rec(N, f):
    u = [1, f(0)]
    for n in range(1, N - 1):
        num = f(n, u[n], u[n - 1])
        q, r = divmod(num, (n + 1) ** 3); assert r == 0
        u.append(q)
    return u
def s7(N):
    u = [1, 4]
    for n in range(1, N - 1):
        num = (2*n+1)*(13*n*n+13*n+4)*u[n] + 3*n*(9*n*n-1)*u[n-1]
        q, r = divmod(num, (n+1)**3); assert r == 0; u.append(q)
    return u
def U5(N):
    u = [1, 6]
    for n in range(1, N - 1):
        num = 2*(2*n+1)*(11*n*n+11*n+3)*u[n] + 4*n*(4*n*n-1)*u[n-1]
        q, r = divmod(num, (n+1)**3); assert r == 0; u.append(q)
    return u
primes = [p for p in range(7, PMAX + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
for name, gen, ms in (('s7', s7, [-10648, -64, -8, -1, 2, 8, 24, 27, 125]),
                      ('A274786', U5, [-15228, -828, -28, -3, 4, 18, 36, 72, 147])):
    u = gen(2 * PMAX + 2)
    fails = {}
    for m in ms:
        for p in primes:
            if m % p == 0: continue
            P2 = p * p
            x = pow(m, -1, P2)
            F = 0
            for k in range(p - 1, -1, -1): F = (F * x + u[k]) % P2
            G = 0
            for k in range(2 * p - 1, p - 1, -1): G = (G * x + u[k]) % P2
            if (G * pow(x, p, P2) - u[1] * x * F) % P2:
                fails.setdefault(m, []).append(p)
    print(f'{name}: {len(primes)} primes 7..{PMAX}; failures of (B1\') at the nonzero CM values: {fails if fails else "none"}')
