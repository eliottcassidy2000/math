# Two-loop prediction of y = n(1 - kappa/rho) on ALL regular tournaments (not only circulants):
#   X_1 = sum_v [2 G_vv H_vv - (3/2) G_vv^2]                      (figure-eight)
#       + sum_(u,v) G_uv (H_uv^2 - O_uv^2 - (2/3) G_uv^2)          (theta)
#       - d^T G d,  d = diag(H - G)                                (dumbbell; 0 for vertex-transitive T)
#   y_pred = -n X_1 - (D2 + D1^2 - D1)/2.
# Usage: two_loop_regular.py n HK_BINARY [res/mod].  Classes from gentourng -d m -D m (optionally split res/mod),
# H by Held-Karp (hk_count reads gentourng codes on stdin).  Also compares with the empirical fit a + b M4.
import subprocess, sys, math
import numpy as np

def two_loop(S):
    n = S.shape[0]; Sh = S / n; I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2); O = Sh @ np.linalg.inv(I - S2)
    gd, hd = np.diag(G), np.diag(H)
    F8 = np.sum(2 * gd * hd - 1.5 * gd ** 2)
    TH = np.sum(G * (H * H - O * O - (2 / 3) * G * G))
    d = hd - gd
    DB = -d @ G @ d
    x = np.linalg.eigvals(Sh).imag
    D1 = np.sum(2 * x ** 2 / (1 - x ** 4)); D2 = np.sum(4 * x ** 2 * (1 + x ** 4) / (1 - x ** 4) ** 2)
    rho = float(np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2))))
    M4 = float(np.sum(x ** 4) / 2)
    return rho, -n * (F8 + TH + DB) - (D2 + D1 ** 2 - D1) / 2, n * F8, n * TH, n * DB, M4

if __name__ == '__main__':
    n = int(sys.argv[1]); HK = sys.argv[2]; m = (n - 1) // 2
    args = ['gentourng', '-q', f'-d{m}', f'-D{m}', str(n)] + (sys.argv[3:4] if len(sys.argv) > 3 else [])
    gen = subprocess.run(args, capture_output=True, text=True).stdout
    codes = [c for c in gen.split() if c and not c.startswith('>')]
    Hs = [int(x) for x in subprocess.run([HK], input=gen, capture_output=True, text=True).stdout.split()]
    assert len(Hs) == len(codes), (len(Hs), len(codes))
    iu = np.triu_indices(n, 1)
    Y, P, DBs, M4s = [], [], [], []
    for c, Hv in zip(codes, Hs):
        A = np.zeros((n, n)); bits = np.array([1 if ch == '1' else 0 for ch in c])
        A[iu[0], iu[1]] = bits; A[iu[1], iu[0]] = 1 - bits
        S = A - A.T
        rho, yp, f8, th, db, M4 = two_loop(S)
        kap = Hv * 2 ** (n - 1) / math.factorial(n)
        Y.append(n * (1 - kap / rho)); P.append(yp); DBs.append(db); M4s.append(M4)
    Y, P, DBs, M4s = map(np.array, (Y, P, DBs, M4s))
    rk = lambda v: np.argsort(np.argsort(v))
    b, a = np.polyfit(M4s, Y, 1)
    fitres = Y - (a + b * M4s)
    print(f"n={n}{' ' + sys.argv[3] if len(sys.argv) > 3 else ''}: {len(Y)} regular classes")
    print(f"  two-loop: corr(y, y_pred) = {np.corrcoef(Y, P)[0, 1]:+.4f}, spearman = {np.corrcoef(rk(Y), rk(P))[0, 1]:+.4f}, "
          f"mean(y - y_pred) = {np.mean(Y - P):+.4f}, rms scatter = {np.std(Y - P):.4f}; y in [{Y.min():.3f}, {Y.max():.3f}] (sd {np.std(Y):.4f})")
    print(f"  dumbbell share: n*DB in [{DBs.min():+.4f}, {DBs.max():+.4f}]; without dumbbells corr = {np.corrcoef(Y, P + DBs)[0, 1]:+.4f}")
    print(f"  empirical a + b M4 (fitted a = {a:.3f}, b = {b:.3f}): rms scatter = {np.std(fitres):.4f}, corr = {np.corrcoef(Y, a + b * M4s)[0, 1]:+.4f}")
    # does the two-loop corrected predictor kappa_2 = rho (1 - y_pred/n) rank kappa and find its maximiser?
    K = np.array([Hv * 2 ** (n - 1) / math.factorial(n) for Hv in Hs])
    R = K / (1 - Y / n)                      # rho recovered
    K2 = R * (1 - P / n)
    print(f"  ranking: spearman(kappa, rho) = {np.corrcoef(rk(K), rk(R))[0, 1]:+.4f}, spearman(kappa, rho(1 - y_pred/n)) = {np.corrcoef(rk(K), rk(K2))[0, 1]:+.4f}; "
          f"argmax kappa found: {int(np.argmax(K)) == int(np.argmax(K2))}; max rel. error |kappa_2/kappa - 1| = {np.max(np.abs(K2 / K - 1)):.2e}", flush=True)
