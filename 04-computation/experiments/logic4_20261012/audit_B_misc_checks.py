"""Audit B (2026-10-12): small exact checks for THM-4622(6), THM-4623(2), THM-4624(3).
(1) Primitive Pythagorean hypotenuses c <= 20000: every prime factor is 1 mod 4 (THM-4622(6)).
(2) Gamma_Q contains (0, 1/5) and (0, 1/25): explicit integer combinations of rational unit vectors, so the
    components of UD(Q^2) are not contained in lines (what the 'density' step of THM-4622(6) needs).
(3) The only solutions of T^3 x = x in N are 1, 2, 4 (used implicitly in THM-4623(2): 'phi(c) is periodic
    of period dividing 3, i.e. in {1,2,4}').
(4) THM-4624(3): p^f - 1 | p^e + 1 with f | 2e, for p <= 97, e <= 12, happens iff (p,f) in {(2,1), (2,2) with
    e odd, (3,1)}; and N = (q+1)/(p^f-1) <= D + 1 = 2e/f + 1 exactly for the six table entries.
"""
from fractions import Fraction
from math import gcd

# (1)
bad = []
for c in range(1, 20001):
    for a in range(1, c):
        b2 = c * c - a * a
        b = int(round(b2 ** 0.5))
        if b * b == b2 and b > 0 and gcd(gcd(a, b), c) == 1:
            n = c; d = 2; ok = True
            while d * d <= n:
                if n % d == 0:
                    if d % 4 != 1: ok = False
                    while n % d == 0: n //= d
                d += 1
            if n > 1 and n % 4 != 1: ok = False
            if not ok: bad.append((a, b, c))
            break
print("(1) primitive hypotenuses c <= 20000 with a prime factor != 1 mod 4:", bad)

# (2)
F = Fraction
u1, u2, e2 = (F(3, 5), F(4, 5)), (F(-3, 5), F(4, 5)), (F(0), F(1))
def add(*vs): return (sum(v[0] for v in vs), sum(v[1] for v in vs))
def sc(k, v): return (k * v[0], k * v[1])
p5 = add(sc(2, add(u1, u2, sc(-1, e2))), sc(-1, e2))          # 2[(3/5,4/5)+(-3/5,4/5)-(0,1)]-(0,1)
w1, w2 = (F(7, 25), F(24, 25)), (F(-7, 25), F(24, 25))
p25 = add(sc(2, add(w1, w2, sc(-1, e2), sc(-4, p5))), sc(-1, p5))  # 2[(7/25,24/25)+(-7/25,24/25)-(0,1)-4(0,1/5)]-(0,1/5)
print("(2) (0,1/5) in Gamma_Q:", p5 == (F(0), F(1, 5)), "; (0,1/25) in Gamma_Q:", p25 == (F(0), F(1, 25)))

# (3)
def T(n): return n // 2 if n % 2 == 0 else 3 * n + 1
sol = [n for n in range(1, 10 ** 6) if T(T(T(n))) == n]
print("(3) solutions of T^3 n = n, n < 10^6:", sol)

# (4)
def is_prime(n): return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))
hits = []; small = []
for p in [p for p in range(2, 98) if is_prime(p)]:
    for e in range(1, 13):
        q = p ** e
        for f in range(1, 2 * e + 1):
            if (2 * e) % f: continue
            if (q + 1) % (p ** f - 1) == 0:
                hits.append((p, e, f))
                N, D = (q + 1) // (p ** f - 1), 2 * e // f
                if N <= D + 1: small.append((q, (p, f), (N, D)))
expected = all((p == 2 and f == 1) or (p == 2 and f == 2 and e % 2 == 1) or (p == 3 and f == 1) for p, e, f in hits)
missing = [(p, e, f) for p in (2, 3) for e in range(1, 13) for f in ((1, 2) if p == 2 else (1,))
           if (p, e, f) not in hits and not (p == 2 and f == 2 and e % 2 == 0)]
print("(4) all divisibility hits are of the three types:", expected, "; expected cases missing:", missing)
print("    N <= D+1 cases:", small)
