# HYP-9479 out-of-sample runs (n = 26, 28; splits (9,5,9), (10,5,10), (9,7,9)).
import numpy as np, math, sys
def block(k_max):
    """yield (k, g on [2^k, 2^(k+1))) for k = 0..k_max, built by doubling: g(2u)=g(u)+u+(u&1), g(2u+1)=g(u)."""
    prev = np.array([1], dtype=np.int64); u = np.array([1], dtype=np.int64)   # block k=0: v=1, g(1)=1
    yield 0, prev
    for k in range(1, k_max + 1):
        cur = np.empty(2 * prev.size, dtype=np.int64)
        cur[0::2] = prev + u + (u & 1)
        cur[1::2] = prev
        u = np.arange(1 << k, 1 << (k + 1), dtype=np.int64)
        prev = cur
        yield k, cur
ths = {'phi': (math.sqrt(5) - 1) / 2, 'sqrt2-1': math.sqrt(2) - 1, 'pi-3': math.pi - 3, 'e-2': math.e - 2, 'sqrt3-1': math.sqrt(3) - 1}
for k, vals in block(27):
    n = k + 1
    if n not in (26, 28): continue
    out = []
    for name, t in ths.items():
        s = 0j
        for i in range(0, vals.size, 1 << 24):
            s += np.exp(2j * np.pi * t * vals[i:i + (1 << 24)].astype(np.float64)).sum()
        out.append(f'{name}={abs(s) / 2**((n - 1) / 2):.3f}')
    c = np.bincount(vals)
    E = int((c.astype(np.int64)**2).sum())
    print(f'n={n}: sqrt-ratios ' + ' '.join(out) + f' | E/|X| = {E / vals.size:.4f}, distinct/|X| = {np.count_nonzero(c) / vals.size:.4f}', flush=True)
# P3: diagonal dominance
def Hw_vec(L):
    h = np.zeros(1 << L, dtype=np.int64)
    for i in range(1 << L):
        H, A, p = 1, 0, 0
        for j in range(L - 1, -1, -1):
            b = (i >> j) & 1
            if b: A = 2 * A + 1; p = 1
            else: H += A + p; A = 2 * A; p = 0
        h[i] = H
    return h
def Cw_vec(L):  # C(w) = sum_j (1 - w_j) 2^(j-1), w_1 = most significant letter
    c = np.zeros(1 << L, dtype=np.int64)
    for i in range(1 << L):
        c[i] = sum((1 - ((i >> (L - 1 - j)) & 1)) << j for j in range(L))
    return c
from fractions import Fraction
thetas = [Fraction(3, 37), Fraction(5, 101), Fraction(11, 257), (math.sqrt(5) - 1) / 2]
res = []
for (a, L, l) in [(9, 5, 9), (10, 5, 10), (9, 7, 9)]:
    Hv = Hw_vec(L); Av = np.arange(1 << L, dtype=np.int64); Cv = Cw_vec(L)
    first = (Av >> (L - 1)) & 1; last = Av & 1
    epsv = (first == 0).astype(np.int64)
    Cw_all = Cw_vec(l)
    for th in thetas:
        t = float(th)
        for p in (0, 1):
            Au = (1 << (a - 1)) + 2 * np.arange(1 << (a - 2), dtype=np.int64) + p
            for eps in (0, 1):
                wfirst = np.array([(i >> (l - 1)) & 1 for i in range(1 << l)])
                Cw = Cw_all[(wfirst == (1 - eps))]
                Hhat = Hv + p * epsv + last * eps
                M1 = np.exp(2j * np.pi * t * np.outer(Au, Cv).astype(np.float64))       # Ku x 2^L
                D = np.exp(2j * np.pi * t * Hhat.astype(np.float64))
                M2 = np.exp(2j * np.pi * t * np.outer(Av, Cw).astype(np.float64))       # 2^L x Kw
                F = (M1 * D) @ M2
                Q = float((np.abs(F)**2).sum())
                diag = (1 << L) * Au.size * Cw.size
                res.append((a, L, l, str(th) if isinstance(th, Fraction) else 'phi', p, eps, Q / diag))
vals = [r[-1] for r in res]
print(f'P3: {len(res)} cases; Q/diag min {min(vals):.4f} max {max(vals):.4f}; outside [0.85, 1.2]: {sum(1 for x in vals if not 0.85 <= x <= 1.2)}')
for r in res:
    if not 0.85 <= r[-1] <= 1.2: print('  outside:', r)
