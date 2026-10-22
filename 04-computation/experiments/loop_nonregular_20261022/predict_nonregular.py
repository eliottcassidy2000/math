"""First- and second-order loop-expansion prediction of y = n(1 - kappa/rho) for NON-regular tournaments, tested on
almost-regular even-n tournaments (scores n/2, n/2 - 1).  Session opus-2026-10-22-S54.

Exact: kappa = sum_L phi_L g_n(L), Phi(z) = sum phi_L z^L = exp(sum_(e,k) X_(e,k)(z)), g_n(L) = n^L/(n)_L
= exp(a/n + b/n^2 + ...), a = C(L,2), b = L(L-1)(2L-1)/12.  X_(0,0) is the one-loop cycle sum (log rho at z = 1);
the leaf classes X_(e,k) come from leaf_engine.py (validated against brute force).  Measured orders on almost-regular
tournaments (class_scaling.out): first order X_(-1,2) (path trees, = -z^2 d^T (I + z^2 Shat^2)^-1 d), X_(1,0); second
order X_(0,1), X_(0,2), X_(2,0); third order X_(-1,3), X_(-1,4), X_(0,3), X_(1,1).

  kappa/rho - 1 = S1 + S2 + O(n^-3),
  S1 = X_(-1,2) + X_(1,0) + E0[a]/n,
  S2 = X_(0,1) + X_(0,2) + X_(2,0) + F^2/2 + F E0[a]/n + dE[a]/n + (E0[a^2]/2 + E0[b])/n^2,   F = X_(-1,2) + X_(1,0),
  dE[a] = (Delta_2 + (2 k_1 - 1) Delta_1)/2,  Delta_j = (z d/dz)^j F at z = 1,  k_j = (z d/dz)^j X_(0,0) at z = 1.
Predictions: y1 = -n S1, y2 = -n (S1 + S2), y2r = y2 - n (third-order classes) ("resummed").
Baselines: mac-mini's two-loop y_pred (= -n X_(1,0) - E0[a]) and y_pred + t, t = -n X_(-1,2) (almost_regular_tree.py).
Usage: python predict_nonregular.py n FILE [OUT]   (FILE lines: code H [aut], gentourng/ar_tool code convention)"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE

FIRST = [(-1, 2), (1, 0)]
SECOND = [(0, 1), (0, 2), (2, 0)]
THIRD = [(-1, 3), (-1, 4), (0, 3), (1, 1)]


def cumulants(f, kmax=4, r=0.05, K=32):
    """(z d/dz)^k f at z = 1, k <= kmax, by a Cauchy sum in theta = log z."""
    th = r * np.exp(2j * np.pi * np.arange(K) / K)
    vals = np.array([f(np.exp(t)) for t in th])
    return [(math.factorial(k) * np.mean(vals * th ** (-k))).real for k in range(kmax + 1)]


def X1_closed(Sh, z):
    n = Sh.shape[0]; I = np.eye(n); X = z * Sh; X2 = X @ X
    G = -X2 @ np.linalg.inv(I + X2); H = X2 @ np.linalg.inv(I - X2); O = X @ np.linalg.inv(I - X2)
    gd, hd = np.diag(G), np.diag(H); d = hd - gd
    return np.sum(2 * gd * hd - 1.5 * gd ** 2) + np.sum(G * (H * H - O * O - (2 / 3) * G * G)) - d @ G @ d


def tree_closed(Sh, z):
    n = Sh.shape[0]; d = Sh @ np.ones(n)
    return -z * z * d @ np.linalg.solve(np.eye(n) + z * z * Sh @ Sh, d)


def predict(S):
    n = S.shape[0]; Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    rho = float(np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2))))
    k = cumulants(lambda z: 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2))))
    k1, k2, k3, k4 = k[1], k[2], k[3], k[4]
    EL = k1; EL2 = k2 + k1 ** 2; EL3 = k3 + 3 * k2 * k1 + k1 ** 3
    EL4 = k4 + 4 * k3 * k1 + 3 * k2 ** 2 + 6 * k2 * k1 ** 2 + k1 ** 4
    Ea = (EL2 - EL) / 2; Ea2 = (EL4 - 2 * EL3 + EL2) / 4; Eb = (2 * EL3 - 3 * EL2 + EL) / 12
    P = LE.LEF.stack(Sh.astype(complex), 1.0)
    Xc = {c: LE.X(Sh, c[0], c[1], 1.0, P).real for c in FIRST + SECOND + THIRD}
    Dl = cumulants(lambda z: tree_closed(Sh.astype(complex), z) + X1_closed(Sh.astype(complex), z), kmax=2)
    dEa = (Dl[2] + (2 * k1 - 1) * Dl[1]) / 2
    F = Xc[(-1, 2)] + Xc[(1, 0)]
    S1 = F + Ea / n
    S2 = sum(Xc[c] for c in SECOND) + F * F / 2 + F * Ea / n + dEa / n + (Ea2 / 2 + Eb) / n ** 2
    S3c = sum(Xc[c] for c in THIRD)
    ypred_mm = -n * Xc[(1, 0)] - Ea                       # mac-mini's two-loop y_pred
    t = -n * Xc[(-1, 2)]
    return dict(rho=rho, y_mm=ypred_mm, y_mm_t=ypred_mm + t, y1=-n * S1, y2=-n * (S1 + S2), y2r=-n * (S1 + S2 + S3c),
                **{f'X{c}': v for c, v in Xc.items()}, dEa=dEa, Ea=Ea)


def load(n, path):
    iu = np.triu_indices(n, 1)
    for line in open(path):
        parts = line.split()
        if not parts:
            continue
        code, H = parts[0], int(parts[1])
        A = np.zeros((n, n)); b = np.array([1 if ch == '1' else 0 for ch in code])
        A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
        yield code, H, (A - A.T).astype(float)


if __name__ == '__main__':
    n = int(sys.argv[1]); path = sys.argv[2]; out = sys.argv[3] if len(sys.argv) > 3 else None
    rows = []
    for code, H, S in load(n, path):
        p = predict(S)
        kap = H * 2 ** (n - 1) / math.factorial(n)
        p['y'] = n * (1 - kap / p['rho'])
        p['code'] = code; p['H'] = H
        rows.append(p)
    Y = np.array([r['y'] for r in rows])
    rk = lambda v: np.argsort(np.argsort(v))
    print(f'n={n}: {len(rows)} tournaments from {os.path.basename(path)}; y in [{Y.min():.4f}, {Y.max():.4f}], sd {Y.std():.4f}')
    for key, label in (('y_mm', "two-loop y_pred (regular formula)"), ('y_mm_t', 'y_pred + path trees t (mac-mini)'),
                       ('y1', 'first order S1'), ('y2', 'second order S1 + S2'), ('y2r', 'second order + third-order classes')):
        P = np.array([r[key] for r in rows]); R = Y - P
        print(f'  {label:38s}: corr {np.corrcoef(Y, P)[0, 1]:+.4f}, spearman {np.corrcoef(rk(Y), rk(P))[0, 1]:+.4f}, '
              f'mean(y - pred) {R.mean():+.4f}, sd(y - pred) {R.std():.4f}, max|y - pred| {np.abs(R).max():.4f}', flush=True)
    means = {k: np.mean([r[k] for r in rows]) for k in rows[0] if k.startswith('X(') or k in ('dEa', 'Ea')}
    print('  class means (n * X): ' + ', '.join(f'{k} {n * v:+.4f}' for k, v in means.items() if k.startswith('X(')) +
          f'; mean dE[a] {means["dEa"]:+.4f}, mean E0[a] {means["Ea"]:+.4f}')
    if out:
        with open(out, 'w') as fh:
            keys = ['code', 'H', 'y', 'y_mm', 'y_mm_t', 'y1', 'y2', 'y2r'] + [k for k in rows[0] if k.startswith('X(')]
            fh.write(' '.join(keys) + '\n')
            for r in rows:
                fh.write(' '.join(str(r[k]) if k in ('code', 'H') else f'{r[k]:.10g}' for k in keys) + '\n')
