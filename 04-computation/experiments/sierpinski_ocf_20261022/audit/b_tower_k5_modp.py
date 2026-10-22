# Audit extension: H(S_5) (243 vertices) by the auditor's shape-sum law, computed modulo 64 primes < 2^28
# with numpy and reconstructed by CRT; compared with the session's 404-digit value.
import sys, re
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from b_tower_H import path_covers, lexprod, shape_sum, C3, walk_counts
from sympy import prevprime
from sympy.ntheory.modular import crt

c1 = path_covers(3, C3)
S2 = lexprod(C3, C3)
c2 = path_covers(9, S2)
c3 = shape_sum([c2, c2, c2], 27)          # exact ordered covers of S_3 (validated in b_tower_H.py)
b = 27
G = {}
for s in range(3):
    for L in range(1, 3 * b + 1):
        k = walk_counts(s, L)
        if max(k) <= b:
            G[k] = G.get(k, 0) + 1
primes = []
p = 2 ** 28
while len(primes) < 64:
    p = prevprime(p)
    primes.append(p)
res = []
for p in primes:
    cv = np.array([0] + [c3[m] % p for m in range(1, b + 1)], dtype=np.int64)
    cur = np.zeros((b + 1,) * 3, dtype=np.int64)
    cur[0, 0, 0] = 1
    cS4 = [0]
    for M in range(1, 3 * b + 1):
        nxt = np.zeros_like(cur)
        for (k0, k1, k2), mult in G.items():
            nxt[k0:, k1:, k2:] += (mult * cur[:b + 1 - k0, :b + 1 - k1, :b + 1 - k2]) % p
            nxt %= p
        cur = nxt
        t = (cur * cv[None, None, :]) % p
        t = t.sum(axis=2) % p
        t = (t * cv[None, :]) % p
        t = t.sum(axis=1) % p
        t = (t * cv) % p
        cS4.append(int(t.sum() % p))
    cc = cS4 + [0]
    h = 3 * sum(cc[a] * (cc[a] ** 2 + cc[a] * cc[a + 1] + cc[a + 1] ** 2) for a in range(1, 3 * b + 1)) % p
    res.append(h)
H5, mod = crt(primes, res)
H5 = int(H5)
s = open(__file__.rsplit('/', 2)[0] + '/ternary_fast_k6.out').read()
sess = int(re.search(r'k=5 n=243: H\(S_k\) = (\d+)', s).group(1))
v = 0; x = H5
while x % 3 == 0: x //= 3; v += 1
print("H(S_5) digits:", len(str(H5)), " v3:", v, " modulus digits:", len(str(int(mod))))
print("equals session's k=5 value:", H5 == sess)
