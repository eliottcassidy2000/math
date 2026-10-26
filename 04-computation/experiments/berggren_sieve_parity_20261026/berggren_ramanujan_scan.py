# Near-Ramanujan scan of the Berggren graphs G_q (vertex adjacency with multiplicities, see berggren_hashimoto2.py)
# and of the signed version (branch character: the A-continuation carries sign -1), primes 5 <= q < QMAX.
# Reports the primes where the bound 2 sqrt 3 is exceeded and the largest ratio.
import sys
import numpy as np
from collections import Counter
from sympy import primerange
INF = 'inf'
def inv(q, t): return 0 if t == INF else (INF if t == 0 else pow(t, -1, q))
def neg(q, t): return t if t == INF else (-t) % q
def Dm(q, t): return INF if t == INF else (t + 2) % q
def S(q, t): return neg(q, inv(q, t))
def U(q, t): return Dm(q, S(q, t))
def A_(q, t): return Dm(q, S(q, neg(q, t)))
def tau(q, t): return Dm(q, neg(q, t))
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
QMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 5
r = 2 * np.sqrt(3); bad = []; worst = (0, 0)
nq = 0
for q in primerange(QMIN, QMAX):
    nq += 1
    pts = list(range(q)) + [INF]
    M = {t: Counter([U(q, t), A_(q, t), Dm(q, t), tau(q, t)]) for t in pts}
    key = {t: frozenset(M[t].items()) for t in pts}
    vk = list(set(key.values())); vid = {k: i for i, k in enumerate(vk)}
    head = {t: vid[key[t]] for t in pts}
    nv = len(vk); Ad = np.zeros((nv, nv))
    for k in vk:
        for f, mult in dict(k).items(): Ad[vid[k], head[f]] += mult
    mu = np.linalg.eigvals(Ad)
    nontriv = [m for m in mu if abs(m - 4) > 1e-7]
    mx = max(abs(m) for m in nontriv)
    if mx / r > worst[0]: worst = (mx / r, q)
    if mx > r + 1e-9: bad.append((q, round(float(mx), 4)))
print(f"primes {QMIN} <= q < {QMAX} ({nq} primes): Ramanujan bound 2 sqrt 3 exceeded at {len(bad)} primes: {bad}")
print(f"fraction non-Ramanujan: {len(bad)/nq:.4f}; max nontrivial |mu| over the range: {max([b[1] for b in bad], default=0)}")
print(f"largest ratio max|mu| / (2 sqrt 3) = {worst[0]:.5f} at q = {worst[1]}")
