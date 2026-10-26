# Six-vertex level of the Paley census tower: TT6(P_p) = p * TT5(T[QR]), TT5(X) = sum_w TT4(N+_X(w)),
# TT4(Y) = sum_u t3(N+_Y(u)), t3(Z) = sum C(d, 2).  Then test TT6 = A(p) + B(p) a_p(f8) exactly (A deg 6, B deg <= 3),
# f8 = eta(2z)^4 eta(4z)^4.  mac-mini-2026-10-25-pyth.
import sys
from fractions import Fraction
import numpy as np
from sympy import primerange
sys.path.insert(0, '.')
from tt5_modular_search import eta_product, solve

def tt6(p):
    sq = sorted({(x * x) % p for x in range(1, p)})
    Q = np.array(sq); m = len(Q)
    A = np.isin((Q[None, :] - Q[:, None]) % p, sq).astype(np.int64)        # tournament on QR
    tot = 0
    for w in range(m):                                                     # TT5(QR) = sum_w TT4(N+(w))
        Y = np.nonzero(A[w])[0]
        AY = A[np.ix_(Y, Y)]
        MY = AY @ AY.T                                                     # common out-neighbours within Y
        for u in range(len(Y)):
            Z = np.nonzero(AY[u])[0]
            d = MY[Z, u]                                                   # out-degree inside N+_Y(u) of each z
            tot += int((d * (d - 1) // 2).sum())
    return p * tot
if __name__ == '__main__':
    Pmax = int(sys.argv[1])
    rows = []
    for p in primerange(19, Pmax):
        if p % 4 != 3: continue
        rows.append((p, tt6(p)))
        print(f"p={p} TT6={rows[-1][1]}", flush=True)
    a = eta_product([(2, 4), (4, 4)], Pmax + 2)
    for degB in (2, 3, 4):
        nun = 7 + degB + 1
        if len(rows) < nun + 4: break
        train = rows[2:2 + nun]
        A_ = [[Fraction(p) ** k for k in range(7)] + [Fraction(p) ** k * a[p] for k in range(degB + 1)] for p, t in train]
        sol = solve(A_, [Fraction(t) for p, t in train])
        if sol is None: print('singular'); continue
        F = lambda p: sum(sol[k] * Fraction(p) ** k for k in range(7)) + sum(sol[7 + k] * Fraction(p) ** k for k in range(degB + 1)) * a[p]
        bad = [(p, t, float(F(p))) for p, t in rows if F(p) != t]
        print(f"TT6 = A(p) + B(p) a_p(f8), deg B = {degB}: exact on {len(rows) - len(bad)} of {len(rows)} primes", bad[:2])
