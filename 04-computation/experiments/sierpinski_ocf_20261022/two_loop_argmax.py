# Does the parameter-free two-loop predictor kappa_2 = rho (1 - y_pred/n) find the H-maximal circulant at each n
# (THM-4751: flat class for n <= 11, R_n for 13 <= n <= 25)?  Also the rank of the true maximiser under rho alone.
import re, glob, math, sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '.')
from two_loop_predictor import predict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m: rows[int(m.group(1))][tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(4))
for line in open('../lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m: rows[25][tuple(int(x) for x in m.group(1).split(','))] = int(m.group(2))
print(' n  classes  true argmax (H)          two-loop argmax      rho argmax       rank of true max: by rho / by kappa_2   margin true(1st/2nd)  pred(1st/2nd)')
for n in sorted(rows):
    if n < 7: continue
    items = []
    for S, H in rows[n].items():
        rho, yp, *_ = predict(n, S)
        items.append((S, H, rho, rho * (1 - yp / n)))
    tH = max(items, key=lambda t: t[1]); tR = max(items, key=lambda t: t[2]); tK = max(items, key=lambda t: t[3])
    rk_rho = sorted(items, key=lambda t: -t[2]).index(tH) + 1
    rk_k2 = sorted(items, key=lambda t: -t[3]).index(tH) + 1
    sH = sorted([t[1] for t in items], reverse=True); sK = sorted([t[3] for t in items], reverse=True)
    def lab(S):
        d = (n - 1) // 2
        if tuple(S) == tuple(range(1, d + 1)): return 'R_n'
        sq = {(x * x) % n for x in range(1, n)}
        if len(sq) == d and set(S) == sq: return 'Paley'
        return 'S=' + ','.join(map(str, S))
    print(f"{n:2d}  {len(items):4d}   {lab(tH[0]):22s}  {lab(tK[0]):18s}  {lab(tR[0]):12s}   {rk_rho:3d} / {rk_k2:3d}      "
          f"{sH[0] / sH[1] - 1 if len(sH) > 1 else 0:+.4f}        {sK[0] / sK[1] - 1 if len(sK) > 1 else 0:+.4f}")
