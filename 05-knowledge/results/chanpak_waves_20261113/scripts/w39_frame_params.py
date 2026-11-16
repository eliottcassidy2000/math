# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
import os
# Wave 39: parameters of the span-3 ribbon semigroup in the thin-semigroup frame.
#  (a) max per-letter growth rho(M_W)^(1/|W|) over all words of length <= 6 (lower bound for the joint spectral radius);
#  (b) Lyapunov (typical) growth of H for uniform random words, and the counting exponent alpha_lower = log 8 / log JSR;
#  (c) mixing of H mod l for l = 5, 7, 11, 13 over all words of length n (TV to uniform), n = 6..9.
import itertools, numpy as np
from collections import Counter
f = open(os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '04-computation', 'experiments', 'ribbon_spectrum_20261104'), 'span3_intrep.txt'))
r, L = map(int, f.readline().split())
y0 = np.array(list(map(int, f.readline().split())), dtype=np.int64)
hc = np.array(list(map(int, f.readline().split())), dtype=np.int64)
M = {}
for _ in range(L):
    mask = int(f.readline())
    M[mask] = np.array([list(map(int, f.readline().split())) for _ in range(r)], dtype=np.int64)
Mf = {b: M[b].astype(float) for b in M}
best = (0, None)
for n in range(1, 7):
    for w in itertools.product(range(8), repeat=n):
        P = np.identity(r)
        for b in w:
            P = Mf[b] @ P
        rad = max(abs(np.linalg.eigvals(P)))
        g = rad ** (1.0 / n)
        if g > best[0] + 1e-9:
            best = (g, w)
    print('words <= %d letters: max rho(M_W)^(1/|W|) = %.5f at %s' % (n, best[0], best[1]), flush=True)
jsr = best[0]
rng = np.random.default_rng(39)
logs = []
for t in range(400):
    v = y0.astype(float)
    acc = 0.0
    n = 400
    for _ in range(n):
        v = Mf[rng.integers(8)] @ v
        s = abs(hc.astype(float) @ v)
        nv = np.linalg.norm(v)
        acc += np.log(nv)
        v = v / nv
    logs.append((acc + np.log(abs(hc.astype(float) @ v) + 1e-300)) / n)
lyap = float(np.mean(logs))
print('Lyapunov growth of H (uniform random words, n = 400): %.4f per letter (sd of (log H)/n %.4f)' % (np.exp(lyap), np.std(logs)))
print('counting exponent lower bound alpha >= log 8 / log JSR = %.3f; with typical growth log 8 / log lyap = %.3f'
      % (np.log(8) / np.log(jsr), np.log(8) / lyap))
for n in range(6, 10):
    cnt = {l: Counter() for l in (5, 7, 11, 13)}
    stack = [((), y0.copy())]
    vals = []
    def rec(depth, v):
        if depth == n:
            vals.append(int(hc @ v))
            return
        for b in range(8):
            rec(depth + 1, M[b] @ v)
    rec(0, y0.copy())
    out = []
    for l in (5, 7, 11, 13):
        c = Counter(x % l for x in vals)
        tv = 0.5 * sum(abs(c.get(a, 0) / len(vals) - 1 / l) for a in range(l))
        out.append('l=%d TV %.4f' % (l, tv))
    print('all %d words of length %d: H mod l: %s' % (len(vals), n, '; '.join(out)), flush=True)
