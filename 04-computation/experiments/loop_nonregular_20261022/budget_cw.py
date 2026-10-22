"""Exact error budget at z = 1 for the weight-graded loop expansion, and the coefficient-wise (Gamma-average) predictor.
Session opus-2026-10-22-S54.

Exact (forest numbers c_L from ar_tool forest):  log(kappa/rho) = A + B,  A = log Phi(1) - log rho (cluster part),
B = log(kappa/Phi(1)) (finite-size part).  The Session30 decorated-core engine (evaluate(M, w, z), imported unchanged)
gives the weight-w class sums W_w(z) (w = e + leaves).
  cluster truncation errors:  A - W1(1),  A - W1(1) - W2(1);
  formal finite-size factor:   B_f = p1/n + [D^2 W1 + (2 c1 - 1) D W1]/(2n) + (p2 - p1^2/2)/n^2  (log of the Session30
                               inverse-factorial operators, second order);
  coefficient-wise predictor: psi_L = [z^L] exp(X0 + W1 + W2) (Cauchy sums, K points on |z| = 1),
                               kappa_cw = sum_(L <= n - 1) psi_L n^L/(n)_L, i.e. the exact Gamma average
                               kappa = E[Phi(n/T)], T ~ Gamma(n + 1), applied to the truncated cluster series.
Writes one row per tournament and prints means and the y-units table (K = 0 skips the coefficient-wise predictor).
Usage: python budget_cw.py n FILE OUT [MAX] [K]"""
import sys, os, math, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'session30_nonregular_loops_20261022',
                                       'session30-decorated-engine'))
sys.path.insert(0, ENGINE)
from decorated import evaluate       # Session30 engine (unchanged)
from predict import predict          # Session30 predictor (unchanged)

n = int(sys.argv[1]); path = sys.argv[2]; out = sys.argv[3]
MAX = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
K = int(sys.argv[5]) if len(sys.argv) > 5 else 32
lines = [l for l in open(path) if l.strip()][:MAX]
codes = [l.split()[0] for l in lines]
cl = subprocess.run([os.path.join(HERE, 'ar_tool2.exe'), 'forest', str(n)], input='\n'.join(codes) + '\n',
                    capture_output=True, text=True).stdout.split('\n')
iu = np.triu_indices(n, 1)
g = np.array([n ** L / math.prod(range(n - L + 1, n + 1)) if L else 1.0 for L in range(n)])
zs = np.exp(2j * np.pi * np.arange(K) / K) if K else None
keys = ['y', 'y1', 'y2', 'ycw0', 'ycw1', 'ycw2', 'A', 'B', 'W1', 'W2', 'Bf', 'Bcw', 'logPhit']
with open(out, 'w') as fh:
    fh.write('code H ' + ' '.join(keys) + '\n')
    rows = []
    for line, fl in zip(lines, cl[:len(lines)]):
        code, H = line.split()[0], int(line.split()[1])
        A_ = np.zeros((n, n), dtype=int); b = np.array([1 if ch == '1' else 0 for ch in code])
        A_[iu[0], iu[1]] = b; A_[iu[1], iu[0]] = 1 - b
        S = A_ - A_.T; M = S / n
        c = [int(t) for t in fl.split()]
        assert sum(c) == H * 2 ** (n - 1)
        p = predict(S)
        rho = p['rho']; kap = H * 2 ** (n - 1) / math.factorial(n)
        Phi1 = sum(c[L] / (math.factorial(n - L) * n ** L) for L in range(n))
        A = math.log(Phi1 / rho); B = math.log(kap / Phi1)
        x2 = np.linalg.eigvalsh(1j * M) ** 2
        if K:
            X0 = np.array([0.5 * np.sum(np.log1p(z * z * x2) - np.log1p(-z * z * x2)) for z in zs])
            W1z = np.array([evaluate(M, 1, z) for z in zs]); W2z = np.array([evaluate(M, 2, z) for z in zs])
            cw = []
            for Xz in (X0, X0 + W1z, X0 + W1z + W2z):
                psi = (np.fft.fft(np.exp(Xz)) / K).real[:n]
                cw.append(float(psi @ g))
            psi2 = (np.fft.fft(np.exp(X0 + W1z + W2z)) / K).real[:n]
            Phit = float(psi2.sum())
        else:                                            # K = 0: skip the coefficient-wise predictor
            cw = [rho, rho, rho]; Phit = rho
        W1 = p['A1']; W2 = p['A2']
        r = dict(y=n * (1 - kap / rho), y1=p['y1'], y2=p['y2'], ycw0=n * (1 - cw[0] / rho), ycw1=n * (1 - cw[1] / rho),
                 ycw2=n * (1 - cw[2] / rho), A=A, B=B, W1=W1, W2=W2, Bf=0.0, Bcw=math.log(cw[2] / Phit), logPhit=math.log(Phit / rho))
        # formal finite-size factor in log form (Session30 operators): b1 = c1 etc. from predict's P1, P2 and D-terms
        d1 = sum(float(np.sum(x2 ** m)) * (2 * m) / m for m in range(1, 151, 2))
        r['Bf'] = p['P1'] / n + ((p['D2_A1'] + (2 * d1 - 1) * p['D_A1']) / 2) / n + (p['P2'] - p['P1'] ** 2 / 2) / n ** 2
        rows.append(r)
        fh.write(f'{code} {H} ' + ' '.join(f'{r[k]:.10g}' for k in keys) + '\n'); fh.flush()
R = {k: np.array([r[k] for r in rows]) for k in keys}
print(f'n={n}: {len(rows)} tournaments ({os.path.basename(path)}), K = {K}')
print(f'  exact cluster part A: mean {R["A"].mean():+.5f};  n (A - W1) {n * np.mean(R["A"] - R["W1"]):+.4f}, '
      f'n (A - W1 - W2) {n * np.mean(R["A"] - R["W1"] - R["W2"]):+.4f} (sd {n * np.std(R["A"] - R["W1"] - R["W2"]):.4f})')
print(f'  exact finite-size part B: mean {R["B"].mean():+.5f};  n (B - Bformal) {n * np.mean(R["B"] - R["Bf"]):+.4f} '
      f'(sd {n * np.std(R["B"] - R["Bf"]):.4f});  n (B - Bcw) {n * np.mean(R["B"] - R["Bcw"]):+.4f}')
print(f'  sum: n (log kappa/rho - [W1 + W2 + Bformal]) {n * np.mean(R["A"] + R["B"] - R["W1"] - R["W2"] - R["Bf"]):+.4f};  '
      + (f'truncated series at z = 1 restricted to L <= n - 1: n (A - logPhit) {n * np.mean(R["A"] - R["logPhit"]):+.4f}' if K else ''))
for k in ('y1', 'y2') + (('ycw0', 'ycw1', 'ycw2') if K else ()):
    res = R['y'] - R[k]
    print(f'  {k:5s}: corr {np.corrcoef(R["y"], R[k])[0, 1]:+.4f}, mean(y - pred) {res.mean():+.4f}, scatter {res.std():.4f}, RMSE {np.sqrt(np.mean(res ** 2)):.4f}')
