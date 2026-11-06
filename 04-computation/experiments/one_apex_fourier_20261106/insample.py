import numpy as np, math
def gtable(N):
    g = np.zeros(N, dtype=np.int64); g[0] = 1
    for v in range(1, N):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
NMAX = 24
G = gtable(1 << NMAX)
ths = {'phi': (math.sqrt(5) - 1) / 2, 'sqrt2-1': math.sqrt(2) - 1, 'pi-3': math.pi - 3, 'e-2': math.e - 2, 'sqrt3-1': math.sqrt(3) - 1}
print('square-root ratios |S_n(theta)| / 2^((n-1)/2):')
for n in range(14, NMAX + 1, 2):
    vals = G[1 << (n - 1): 1 << n].astype(np.float64)
    row = [abs(np.exp(2j * np.pi * t * vals).sum()) / 2**((n - 1) / 2) for t in ths.values()]
    print(f'  n={n}: ' + ' '.join(f'{k}={r:.3f}' for k, r in zip(ths, row)))
print('energy ratio E_n / |X_n| (pairs of words of length n, leading 1, with equal count):')
for n in range(12, NMAX + 1, 2):
    vals = G[1 << (n - 1): 1 << n]
    c = np.bincount(vals)
    E = int((c.astype(np.int64)**2).sum())
    print(f'  n={n}: E/|X| = {E / (1 << (n - 1)):.4f}  distinct values / |X| = {np.count_nonzero(c) / (1 << (n - 1)):.4f}')
