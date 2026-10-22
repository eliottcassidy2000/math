# Paley residual y_q = q (1 - kappa(P_q)/rho(P_q)), rho(P_q) = ((q+1)/(q-1))^((q-1)/2); exact H(P_q) from ie_ham64
# (q prime) and paley_extremal_20261021/P27.out (q = 27, F_27).  mac-mini-2026-10-22-sierp.
from mpmath import mp, mpf, factorial, matrix, lu_solve, e
mp.dps = 30
H = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349, 27: 418425945688463432031,
     31: 19876729650849384218015865}
def kap(q): return mpf(H[q]) * 2 ** (q - 1) / factorial(q)
def rho(q): return (mpf(q + 1) / (q - 1)) ** ((q - 1) // 2)
def extrap(ns, f):
    k = len(ns); A = matrix(k, k); b = matrix(k, 1)
    for i, n in enumerate(ns):
        for j in range(k): A[i, j] = mpf(1) / mpf(n) ** j
        b[i] = f(n)
    return lu_solve(A, b)[0]
for q in sorted(H):
    k, r = kap(q), rho(q)
    print(f"q={q:3d} kappa={mp.nstr(k, 12)} rho={mp.nstr(r, 12)} kappa/rho={mp.nstr(k / r, 10)} y_q={mp.nstr(q * (1 - k / r), 8)}  q(1-kappa/e)={mp.nstr(q * (1 - k / e), 8)}")
qs = sorted(H)
y = lambda q: q * (1 - kap(q) / rho(q))
for k in (2, 3, 4):
    print(f"extrapolated y_inf from last {k} (q = {qs[-k:]}): {mp.nstr(extrap(qs[-k:], y), 8)}")
for k in (2, 3, 4):
    print(f"extrapolated kappa_inf from last {k}: {mp.nstr(extrap(qs[-k:], kap), 10)}  (e = {mp.nstr(e, 10)})")
