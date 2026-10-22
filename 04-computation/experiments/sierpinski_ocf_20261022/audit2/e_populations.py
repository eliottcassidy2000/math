# Audit2 (D): recompute the HYP-9365 population statistics independently.
#  - all circulant classes n = 9..25 listed in the repo (tournament_pseudoprimes_20261021/circ_*.out,
#    lattice_bridges_20261021/circ25_H.txt): H recomputed with hcount.c (inclusion-exclusion), compared with the repo,
#    then y = n(1 - kappa/rho) vs the two-loop y_pred (e_predictor.py, dense matrices; FFT cross-check);
#  - all regular tournaments n = 9 (15) and n = 11 (1223) from gentourng: H by Held-Karp (hcount.c), same statistics,
#    plus the ranking claims (Spearman(kappa, rho), Spearman(kappa, rho (1 - y_pred/n)), argmax).
# Spearman with average ranks (scipy).  Usage: e_populations.py HCOUNT_BINARY [circ|reg|all]
import sys, re, glob, subprocess, math
import numpy as np
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, '.')
from e_predictor import two_loop_dense, two_loop_circ, y_exact, rho_circ_exact
from common import circulant_A, skew, code_to_A, regular_classes

HC = sys.argv[1]
mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
BASE = '../../tournament_pseudoprimes_20261021/'
LB = '../../lattice_bridges_20261021/circ25_H.txt'

def stats(Y, P, label, extra=''):
    Y, P = np.array(Y), np.array(P)
    r = pearsonr(Y, P)[0]; s = spearmanr(Y, P)[0]
    print(f"  {label}: {len(Y)} classes  corr(y, y_pred) = {r:+.4f}  spearman = {s:+.4f}  mean(y - y_pred) = {np.mean(Y - P):+.4f}  "
          f"sd(y - y_pred) = {np.std(Y - P):.4f}  y in [{Y.min():.3f}, {Y.max():.3f}]{extra}", flush=True)

if mode in ('circ', 'all'):
    rows = {}
    for f in sorted(glob.glob(BASE + 'circ_*.out')):
        for line in open(f):
            m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
            if m: rows.setdefault(int(m.group(1)), []).append((tuple(int(x) for x in m.group(2).split(',')), int(m.group(4))))
    for line in open(LB):
        m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
        if m: rows.setdefault(25, []).append((tuple(int(x) for x in m.group(1).split(',')), int(m.group(2))))
    print('Circulant classes: H recomputed (hcount iecirc) vs repo; two-loop statistics')
    for n in sorted(rows):
        if n < 9: continue
        Y, P, nbad, maxdev = [], [], 0, 0.0
        for conn, Hrepo in rows[n]:
            Hm = int(subprocess.run([HC, 'iecirc', str(n)] + [str(c) for c in conn], capture_output=True, text=True).stdout)
            if Hm != Hrepo: nbad += 1; print(f"   H MISMATCH n={n} S={conn}: repo {Hrepo} vs {Hm}")
            d = two_loop_dense(skew(circulant_A(n, conn)))
            c = two_loop_circ(n, conn)
            maxdev = max(maxdev, abs(d['y_pred'] - c['y_pred']), abs(d['rho'] - c['rho']))
            rho = rho_circ_exact(n, conn)
            y, _ = y_exact(Hm, n, rho)
            Y.append(y); P.append(d['y_pred'])
        stats(Y, P, f"n={n}", f"  [H mismatches: {nbad}; dense-vs-FFT max dev {maxdev:.1e}]")

if mode in ('reg', 'all'):
    print('All regular tournaments (gentourng -d m -D m n): H by Held-Karp')
    for n in (9, 11):
        codes = regular_classes(n)
        out = subprocess.run([HC, 'hk'], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split()
        Hs = [int(x) for x in out]
        assert len(Hs) == len(codes)
        Y, P, K, R, DB = [], [], [], [], []
        for c, Hv in zip(codes, Hs):
            d = two_loop_dense(skew(code_to_A(c)))
            y, kap = y_exact(Hv, n, d['rho'])
            Y.append(y); P.append(d['y_pred']); K.append(float(kap)); R.append(d['rho']); DB.append(d['nDB'])
        Y, P, K, R = map(np.array, (Y, P, K, R))
        K2 = R * (1 - P / n)
        imax = int(np.argmax(K)); imax2 = int(np.argmax(K2)); imaxR = int(np.argmax(R))
        stats(Y, P, f"n={n}", f"  n*dumbbell in [{min(DB):+.2e}, {max(DB):+.2e}]")
        print(f"     ranking: spearman(kappa, rho) = {spearmanr(K, R)[0]:+.4f}, spearman(kappa, rho(1 - y_pred/n)) = {spearmanr(K, K2)[0]:+.4f}; "
              f"argmax kappa = class {imax} (H = {Hs[imax]}, code {codes[imax]}), argmax of corrected predictor = class {imax2}, argmax rho = class {imaxR}; "
              f"max |kappa_2/kappa - 1| = {np.max(np.abs(K2 / K - 1)):.2e}", flush=True)
        # how many distinct H values / ties
        print(f"     distinct H values: {len(set(Hs))} of {len(Hs)} classes")
