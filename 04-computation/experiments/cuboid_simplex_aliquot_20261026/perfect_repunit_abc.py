# (1) preimages of the perfect numbers 6, 28, 496 under s; (2) coincidences s(p^a) = s(q^b) between prime powers
# (s(p^a) = 1 + p + ... + p^(a-1), so these are repunits in two prime bases; Goormaghtigh: 31 = 11111_2 = 111_5);
# (3) the Gersonides/Catalan contacts of the power chain of 3 with the doubling chain of 1 (|3^j - 2^i| = 1) and the
# radicals of the resulting abc triples; (4) abc quality of the Mersenne triples 1 + M_p = 2^p behind even perfect numbers.
from math import log, gcd
import numpy as np
from sympy import isprime, primerange, factorint
def rad(n):
    r = 1
    for p in factorint(n): r *= p
    return r
for v in (6, 28, 496):
    N = (v - 1) ** 2 + 1
    sig = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1): sig[d::d] += d
    pre = [n for n in range(2, N + 1) if sig[n] - n == v and not isprime(n)]
    print(f"s^-1({v}) = {pre}")
# repunit coincidences in prime bases
vals = {}
for p in primerange(2, 2000):
    t, e = 1, 0
    acc = 1
    for a in range(2, 64):
        acc = acc * p + 1 if a > 2 else 1 + p        # 1 + p + ... + p^(a-1)
        if acc > 10 ** 12: break
        vals.setdefault(acc, []).append((p, a))
both = {v: L for v, L in vals.items() if len({p for p, a in L}) > 1}
print("values sigma(p^(a-1)) = 1 + p + ... + p^(a-1) shared by two prime bases (<= 1e12, p < 2000):",
      {v: [(p, a) for p, a in L] for v, L in sorted(both.items())})
# Gersonides contacts
cont = [(3 ** j, 2 ** i) for j in range(0, 40) for i in range(0, 64) if abs(3 ** j - 2 ** i) == 1]
print("|3^j - 2^i| = 1:", cont)
for x, y in cont:
    a, c = sorted((x, y)); b = 1
    A, B, C = sorted((1, a, c))[0], None, None
    lo, hi = min(x, y), max(x, y)
    trip = (1, lo, hi)
    r = rad(1 * lo * hi)
    q = log(hi) / log(r) if r > 1 else float('nan')
    print(f"   abc triple 1 + {lo} = {hi}: rad = {r}, quality log c / log rad = {q:.4f}")
for p in (2, 3, 5, 7, 13):
    M = 2 ** p - 1; n = 2 ** (p - 1) * M
    r = rad(M * 2 ** p)
    print(f"   perfect {n} = 2^{p-1} * {M}: Mersenne triple 1 + {M} = {2 ** p}, rad = {r}, quality {log(2 ** p) / log(r):.4f}")
