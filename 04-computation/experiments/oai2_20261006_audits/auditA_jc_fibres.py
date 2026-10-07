"""Audit A: fibre point counts of THM-1300's map over F_q, including non-prime q (GF(p^2), GF(p^3)), brute force over all (x,y,z).
F = (u^3 z + y^2 u (4 + 3xy), y + 3 x u^2 z + 3 x y^2 (4 + 3xy), 2x - 3x^2 y - x^3 z), u = 1 + xy."""
import sys, numpy as np
from auditA_gf import GF
expected = lambda q: {1: (2*q*q - 2*q + 1, q*q - q + 1), 2: (q*q - q + 1, q*q + 1), 3: (2*q*q - q, q*q - q)}
for (p, m) in [(5, 1), (7, 1), (11, 1), (13, 1), (5, 2), (7, 2), (11, 2), (5, 3)]:
    F = GF(p, m); q = F.q
    A = np.array(F.add, dtype=np.int32); M = np.array(F.mul, dtype=np.int32); N = np.array(F.neg, dtype=np.int32)
    c = lambda k: k % p                      # the field element k*1 (prime-field constants have index k mod p)
    ad = lambda a, b: A[a, b]; mu = lambda a, b: M[a, b]; sub = lambda a, b: A[a, N[b]]
    X, Y, Z = np.meshgrid(np.arange(q), np.arange(q), np.arange(q), indexing='ij')
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    xy = mu(X, Y); u = ad(c(1), xy)
    u2 = mu(u, u); u3 = mu(u2, u); y2 = mu(Y, Y); x2 = mu(X, X); x3 = mu(x2, X)
    w = ad(c(4), mu(c(3), xy))                # 4 + 3xy
    F1 = ad(mu(u3, Z), mu(mu(y2, u), w))
    F2 = ad(ad(Y, mu(mu(mu(c(3), X), u2), Z)), mu(mu(mu(c(3), X), y2), w))
    F3 = sub(sub(mu(c(2), X), mu(mu(c(3), x2), Y)), mu(x3, Z))
    ok = True
    out = []
    for j, Fj in ((1, F1), (2, F2), (3, F3)):
        cnt = np.bincount(Fj, minlength=q)
        z0, nz = int(cnt[0]), set(int(v) for v in cnt[1:])
        e0, e1 = expected(q)[j]
        good = (z0 == e0) and nz == {e1}
        ok &= good
        out.append(f"F{j}: zero {z0} (exp {e0}), nonzero {sorted(nz)} (exp {e1}) {'ok' if good else 'MISMATCH'}")
    print(f"q={q} (p={p}, m={m}): " + "; ".join(out), flush=True)
