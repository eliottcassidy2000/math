# Check the closed two-loop formula order by order in z (Shat -> z Shat) against the brute-force excess-1
# cluster sums of cluster_check.py (L = 4, 6).
import sys
import numpy as np
sys.path.insert(0, '.')
from cluster_check import tournament

def X1_series(Sh, LMAX):
    n = Sh.shape[0]
    P = [np.eye(n)]
    for _ in range(LMAX): P.append(P[-1] @ Sh)
    eps = lambda a: (-1) ** (a // 2)
    coef = {}
    for L in range(2, LMAX + 1):
        f8 = 0.0
        for a in range(2, L, 2):
            b = L - a
            if b < 2 or b % 2: continue
            da, db = np.diag(P[a]), np.diag(P[b])
            f8 += np.sum(2 * eps(a) * da * db - 1.5 * eps(a) * eps(b) * da * db)
        th = 0.0
        for a in range(1, L):
            for b in range(1, L - a):
                c = L - a - b
                if c < 1: continue
                M = P[a] * P[b] * P[c]
                if a % 2 == 0 and b % 2 == 0 and c % 2 == 0:
                    k = 4 * (eps(a) + eps(b) + eps(c)) - 8 * eps(a) * eps(b) * eps(c)
                elif (a % 2) + (b % 2) + (c % 2) == 2:
                    ev = [x for x in (a, b, c) if x % 2 == 0][0]
                    k = -4 * eps(ev)
                else: k = 0
                th += k * np.sum(M) / 12
        db = 0.0
        for a in range(2, L, 2):
            for b in range(2, L - a, 2):
                c = L - a - b
                if c < 2 or c % 2: continue
                if a % 4 == 2 and b % 4 == 2:
                    db += (-32 * eps(c)) * (np.diag(P[a]) @ P[c] @ np.diag(P[b])) / 8
        coef[L] = (f8, th, db, f8 + th + db)
    return coef

for kind, n, bf in (('rot', 15, {4: -1.991e-01, 6: -6.025e-02}), ('paley', 19, {4: -1.627e-01, 6: -6.814e-03})):
    A = tournament(kind, n); S = A - A.T
    c = X1_series(S / n, 6)
    for L in (4, 6):
        print(f"{kind} n={n} L={L}: formula = {c[L][3]:+.4e} (F8 {c[L][0]:+.3e}, theta {c[L][1]:+.3e}, dumbbell {c[L][2]:+.1e});  brute force = {bf[L]:+.4e}")
