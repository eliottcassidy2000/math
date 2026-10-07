"""Side checks for the S19 note (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, sections 1.8, 3, 5.
  (H) 3x+d has unboundedly many one-run cycles: d_R = lcm_{r<R} (2^(r+e_r) - 3^(r+1)) carries the cycles
      x_r = d_R (3^(r+1) - 2^(r+1)) / (2^(r+e_r) - 3^(r+1)), e_r least with 2^(r+e_r) > 3^(r+1).
  (S) The Mersenne run as an exact saddle passage at -1: g(x) = (3x+1)/2, g^k(x) + 1 = (3/2)^k (x+1);
      g^a(2^a - 1) = 3^a - 1 (2-adic depth a in, 3-adic depth a out); the product formula for 3/2.
  (P) Paley tournaments P_p, p = 3 mod 4 prime < 400: colour refinement after individualizing one vertex gives 3 classes;
      after individualizing one arc 0 -> 1 (two witnessed symmetric choices) the colouring is discrete.
Prints ALL CHECKS PASSED."""
from math import gcd
from fractions import Fraction
from sympy import isprime, factorint

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m)
    if not c:
        FAIL.append(m)


def T_d(x, d):
    return x // 2 if x % 2 == 0 else (3 * x + d) // 2


print('(H) one-run cycles of 3x+d')
for R in (3, 5, 7, 9):
    pairs = []
    for r in range(R):
        e = 2
        while 2 ** (r + e) <= 3 ** (r + 1):
            e += 1
        pairs.append((r, e))
    d = 1
    for r, e in pairs:
        Dd = 2 ** (r + e) - 3 ** (r + 1)
        d = d * Dd // gcd(d, Dd)
    cycles = []
    for r, e in pairs:
        Dd = 2 ** (r + e) - 3 ** (r + 1)
        num = d * (3 ** (r + 1) - 2 ** (r + 1))
        assert num % Dd == 0
        x = num // Dd
        y, orb = x, [x]
        for _ in range(2 * (r + e) + 2):
            y = T_d(y, d)
            if y == x:
                break
            orb.append(y)
        par = ''.join('1' if z % 2 else '0' for z in orb)
        cycles.append((y == x, len(orb), par))
    ok = all(c[0] for c in cycles) and len(set(c[1] for c in cycles)) == R and d % 2 == 1 and d % 3 != 0
    check(ok, f'R = {R}: d = {d} (gcd(d,6) = 1) has {R} distinct positive one-run cycles, periods {[c[1] for c in cycles]}')

print('(S) the Mersenne run is an exact saddle passage at the fixed point -1 of g(x) = (3x+1)/2')
okS = True
for a in range(2, 60):
    x = Fraction(2 ** a - 1)
    for k in range(a):
        x = (3 * x + 1) / 2
    okS &= (x == 3 ** a - 1)
check(okS, 'g^a(2^a - 1) = 3^a - 1 for 2 <= a < 60 (entry: 2-adic distance 2^-a from -1; exit: 3-adic distance 3^-a from -1)')
check(Fraction(3, 2) * 2 * Fraction(1, 3) == 1, '|3/2|_R |3/2|_2 |3/2|_3 = (3/2)(2)(1/3) = 1: the run is volume-neutral on R x Q_2 x Q_3; Dulac exponent log 3/log 2 = log_2 3')


def refine(n, adj, col):
    while True:
        sig = []
        for v in range(n):
            outs = tuple(sorted(col[u] for u in range(n) if adj[v][u]))
            ins = tuple(sorted(col[u] for u in range(n) if adj[u][v]))
            sig.append((col[v], outs, ins))
        keys = {s: i for i, s in enumerate(sorted(set(sig)))}
        new = [keys[s] for s in sig]
        if len(set(new)) == len(set(col)):
            return new
        col = new


print('(P) Paley tournaments: one witnessed arc choice + colour refinement')
res = []
for p in range(3, 400):
    if p % 4 == 3 and isprime(p):
        qr = {(k * k) % p for k in range(1, p)}
        adj = [[(j - i) % p in qr for j in range(p)] for i in range(p)]
        c1 = refine(p, adj, [1 if v == 0 else 0 for v in range(p)])
        c2 = refine(p, adj, [1 if v == 0 else (2 if v == 1 else 0) for v in range(p)])
        res.append((p, len(set(c1)), len(set(c2))))
check(all(c1 == 3 for _, c1, _ in res) and all(c2 == p for p, _, c2 in res),
      f'{len(res)} primes p = 3 mod 4 below 400: one individualized vertex leaves 3 classes; one individualized arc gives a discrete colouring')
print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}')
