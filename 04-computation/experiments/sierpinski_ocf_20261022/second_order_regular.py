# Second-order (three-loop) prediction on all regular 11-tournaments (non-vertex-transitive: bridged cores contribute)
# and on the ternary towers S_3, S_4 (n = 27, 81).  mac-mini-2026-10-22-sierp.
import subprocess, sys, math, re
import numpy as np
sys.path.insert(0, '.')
from second_order import predict2
from mpmath import mpf, factorial
HK = sys.argv[1]
n = 11; m = 5
gen = subprocess.run(['gentourng', '-q', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout
codes = [c for c in gen.split() if c]
Hs = [int(x) for x in subprocess.run([HK], input=gen, capture_output=True, text=True).stdout.split()]
iu = np.triu_indices(n, 1)
r1, r2, K, K1, K2 = [], [], [], [], []
for c, H in zip(codes, Hs):
    A = np.zeros((n, n)); b = np.array([1 if ch == '1' else 0 for ch in c]); A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
    rho, first, second, X2 = predict2(A - A.T)
    k = H * 2 ** (n - 1) / math.factorial(n)
    ex = k / rho - 1
    r1.append(ex - first); r2.append(ex - first - second); K.append(k); K1.append(rho * (1 + first)); K2.append(rho * (1 + first + second))
r1, r2, K, K1, K2 = map(np.array, (r1, r2, K, K1, K2))
rk = lambda v: np.argsort(np.argsort(v))
print(f"all regular n=11 ({len(K)} classes): residual 1st order mean {r1.mean():+.3e} sd {r1.std():.3e}; 2nd order mean {r2.mean():+.3e} sd {r2.std():.3e}")
print(f"  ranking of kappa: spearman with 1st-order predictor {np.corrcoef(rk(K), rk(K1))[0, 1]:+.4f}, with 2nd-order {np.corrcoef(rk(K), rk(K2))[0, 1]:+.4f}; "
      f"argmax found: 1st {int(np.argmax(K1)) == int(np.argmax(K))}, 2nd {int(np.argmax(K2)) == int(np.argmax(K))}", flush=True)
H = {}
for line in open('ternary_fast_k6.out'):
    mm = re.match(r'k=(\d+) n=(\d+): H\(S_k\) = (\d+)', line)
    if mm: H[int(mm.group(1))] = (int(mm.group(2)), int(mm.group(3)))
for k in (2, 3, 4):
    nn, h = H[k]
    S = np.zeros((nn, nn))
    for u in range(nn):
        for v in range(nn):
            if u == v: continue
            t = (v - u) % nn
            while t % 3 == 0: t //= 3
            S[u, v] = 1 if t % 3 == 1 else -1
    rho, first, second, X2 = predict2(S)
    ex = float(mpf(h) * mpf(2) ** (nn - 1) / factorial(nn)) / rho - 1
    print(f"tower S_{k} (n={nn}): exact {ex:+.6e}; resid 1st {ex - first:+.3e} (n^2 x {nn * nn * (ex - first):+.3f}); resid 2nd {ex - first - second:+.3e} (n^3 x {nn ** 3 * (ex - first - second):+.3f})", flush=True)
