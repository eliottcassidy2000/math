"""Side checks for the S19 note (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, sections 1.8, 3, 5.
  (H) Remark H: d_R = lcm_{r<R} (2^(r+e_r) - 3^(r+1)) gives 3x + d_R at least R one-run cycles, but they are NON-PRIMITIVE
      (gcd(x, d) > 1: scaled copies m*Gamma of cycles of 3x + Delta_r, since T_(md)(mx) = m T_d(x)); Lagarias (1990) proves far
      more (infinitely many k with >= k^(1-eps) primitive cycles).  For d = 35 we also list the primitive cycles found by search.
  (S) The Mersenne run as an exact saddle passage at -1: g(x) = (3x+1)/2, g^k(x) + 1 = (3/2)^k (x+1);
      g^a(2^a - 1) = 3^a - 1 (2-adic depth a in, 3-adic depth a out); the product formula for 3/2.
  (P) Paley tournaments P_p, p = 3 mod 4 prime < 2000: colour refinement after individualizing one vertex gives 3 classes;
      after individualizing the arc 0 -> 1 (two witnessed symmetric choices: a vertex, then an out-neighbour) it is discrete.
  (L) Round 2 of that refinement reads the Legendre family: 8 N_(e,d)(x) = p - ed - e - d + ed S(x) + O(1), where
      S(x) = sum_y chi(y(y-1)(y-x)) = -a_p(Y^2 = X(X-1)(X-x)) and the O(1) term (from y in {0,1,x}) is constant on round-1 classes.
Prints ALL CHECKS PASSED.  Runtime about 30 s."""
from math import gcd
import numpy as np
from sympy import isprime, primerange

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m)
    if not c:
        FAIL.append(m)


def T_d(x, d):
    return x // 2 if x % 2 == 0 else (3 * x + d) // 2


print('(H) one-run cycles of 3x+d from d_R = lcm(Delta_r): at least R, all non-primitive')
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
        cycles.append((y == x, len(orb), gcd(x, d)))
    ok = all(c[0] for c in cycles) and len(set(c[1] for c in cycles)) == R and d % 2 == 1 and d % 3 != 0
    nonprim = all(c[2] > 1 for c in cycles[1:])   # r = 0 has Delta = 1, so x = d itself
    check(ok and nonprim, f'R = {R}: d = {d} carries {R} distinct one-run cycles, periods {[c[1] for c in cycles]}, '
          f'gcd(x, d) = {[c[2] for c in cycles]} (non-primitive)')
# all positive cycles of 3x+35 met from starts < 20000 (orbits followed until they repeat): primitive = gcd(cycle, 35) = 1
d = 35
cyc_min = set()
for s in range(1, 20000):
    x, pos, k = s, {}, 0
    while x not in pos and k < 5000 and x < 10 ** 9:
        pos[x] = k
        x = T_d(x, d)
        k += 1
    if x in pos:
        cyc = [y for y, i in pos.items() if i >= pos[x]]
        cyc_min.add(min(cyc))
allc = sorted(cyc_min)
prim = [m for m in allc if gcd(m, d) == 1]
check(len(allc) >= 9 and prim == [13, 17],
      f'3x+35: positive cycles met from starts < 20000 have minima {allc} ({len(allc)} cycles); primitive (gcd 1 with 35): {prim}; '
      'the construction yields only non-primitive ones (minima 35, 25, 7 type)')

print('(S) the Mersenne run is an exact saddle passage at the fixed point -1 of g(x) = (3x+1)/2')
from fractions import Fraction
okS = True
for a in range(2, 60):
    x = Fraction(2 ** a - 1)
    for k in range(a):
        x = (3 * x + 1) / 2
    okS &= (x == 3 ** a - 1)
check(okS, 'g^a(2^a - 1) = 3^a - 1 for 2 <= a < 60 (entry: 2-adic distance 2^-a from -1; exit: 3-adic distance 3^-a from -1)')
check(Fraction(3, 2) * 2 * Fraction(1, 3) == 1, '|3/2|_R |3/2|_2 |3/2|_3 = (3/2)(2)(1/3) = 1; Dulac exponent log 3/log 2 = log_2 3')


def refine(A, col):
    n = A.shape[0]
    while True:
        k = col.max() + 1
        oh = np.zeros((n, k), dtype=np.float32)
        oh[np.arange(n), col] = 1
        cnt = (A @ oh).astype(np.int64)            # out-neighbour colour counts (in-counts follow for tournaments)
        sig = np.concatenate([col[:, None], cnt], axis=1)
        _, new = np.unique(sig, axis=0, return_inverse=True)
        new = new.ravel()
        if new.max() + 1 == k:
            return col
        col = new


print('(P) Paley tournaments: two witnessed choices (a vertex, an out-neighbour) + colour refinement')
res = []
for p in primerange(3, 2000):
    if p % 4 != 3:
        continue
    qr = np.zeros(p, dtype=bool)
    qr[(np.arange(1, p) ** 2) % p] = True
    A = qr[(np.arange(p)[None, :] - np.arange(p)[:, None]) % p].astype(np.float32)
    c1 = np.zeros(p, dtype=np.int64); c1[0] = 1
    c2 = np.zeros(p, dtype=np.int64); c2[0] = 1; c2[1] = 2
    res.append((p, refine(A, c1).max() + 1, refine(A, c2).max() + 1))
check(all(n1 == 3 for _, n1, _ in res) and all(n2 == p for p, _, n2 in res),
      f'{len(res)} primes p = 3 mod 4 below 2000: one individualized vertex leaves 3 classes; one individualized arc 0 -> 1 gives a discrete colouring')


def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


print('(L) round 2 reads the Legendre family Y^2 = X(X-1)(X-x)')
okL = True
for p in (7, 11, 19, 23, 43, 47):
    for x in range(2, p):
        S = sum(chi(y * (y - 1) * (y - x), p) for y in range(p))
        bds = set()
        for e in (1, -1):
            for dd in (1, -1):
                N = sum(1 for y in range(p) if y not in (0, 1, x) and chi(y, p) == e and chi(y - 1, p) == dd and chi(y - x, p) == 1)
                full = sum((1 + e * chi(y, p)) * (1 + dd * chi(y - 1, p)) * (1 + chi(y - x, p)) for y in range(p) if y not in (0, 1, x))
                okL &= (8 * N == full)
                bd = full - (p - e * dd - e - dd + e * dd * S)
                okL &= abs(bd) <= 12
check(okL, '8 N_(e,d)(x) = p - ed - e - d + ed S(x) + O(1) with |O(1)| <= 12, S(x) = -a_p(Legendre E_x), for p in {7, 11, 19, 23, 43, 47}')
print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}')
