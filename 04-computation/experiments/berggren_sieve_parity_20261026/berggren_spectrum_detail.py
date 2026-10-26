# Detail of the Berggren-mod-q spectra: modulus classes of the nontrivial eigenvalues of T and T_chi on P1(F_q).
import numpy as np
from sympy import primerange
from collections import Counter
INF = 'inf'
def act(q, t, g):
    if g == 'D': return INF if t == INF else (t + 2) % q
    if t == INF: return 2 % q
    if t == 0: return INF
    inv = pow(t, -1, q)
    return (2 - inv) % q if g == 'U' else (2 + inv) % q
r3 = 1 / np.sqrt(3)
exc = []
for q in primerange(3, 400):
    pts = list(range(q)) + [INF]; idx = {p: i for i, p in enumerate(pts)}; N = len(pts)
    T = np.zeros((N, N)); Tc = np.zeros((N, N))
    for p in pts:
        for g, s in (('U', 1), ('A', -1), ('D', 1)):
            T[idx[p], idx[act(q, p, g)]] += 1 / 3; Tc[idx[p], idx[act(q, p, g)]] += s / 3
    ev = np.linalg.eigvals(T); evc = np.linalg.eigvals(Tc)
    mods = np.abs(ev); mc = np.abs(evc)
    nontriv = sorted(mods)[:-1]                      # drop the eigenvalue 1
    cls = Counter('one' if abs(x - 1) < 1e-6 else 'r3' if abs(x - r3) < 1e-6 else 'zero' if x < 1e-6 else
                  ('big' if x > r3 + 1e-6 else 'small') for x in nontriv)
    clsc = Counter('third' if abs(x - 1 / 3) < 1e-6 else 'r3' if abs(x - r3) < 1e-6 else 'zero' if x < 1e-6 else
                   ('big' if x > r3 + 1e-6 else 'small') for x in mc)
    if cls['big'] or cls['one'] or clsc['big'] or q < 30:
        big = sorted(x for x in nontriv if x > r3 + 1e-6)
        bigc = sorted(x for x in mc if x > r3 + 1e-6)
        exc.append((q, q % 8, dict(cls), dict(clsc), [round(x, 4) for x in big], [round(x, 4) for x in bigc]))
for e in exc: print(e)
