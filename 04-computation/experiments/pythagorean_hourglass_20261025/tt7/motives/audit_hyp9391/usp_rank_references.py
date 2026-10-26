# Audit of HYP-9391 (c): can the statistics distinguish USp(4) from higher-rank symplectic laws?  Trace laws of
# USp(2n), n = 2, 3, 4, by quadrature of the Weyl density prod_{i<j} (cos t_i - cos t_j)^2 prod sin^2 t_i on a
# midpoint grid; the N(0,1) limit; and KS distances of the auditor's new parts (part_c_analysis.py fit on [307, 9999],
# saved as mE.npy) to each law.  Usage: usp_rank_references.py mE.npy ; writes usp_rank_references.out
import sys, os, itertools
import numpy as np
from math import erf, sqrt
HERE = os.path.dirname(os.path.abspath(__file__))
data = np.load(sys.argv[1]); m = data[:, 1:]
out = []
def law(n, g):
    th = (np.arange(g) + 0.5) * np.pi / g
    grids = np.meshgrid(*([th] * n), indexing='ij')
    w = np.ones_like(grids[0])
    for i in range(n): w *= np.sin(grids[i]) ** 2
    for i, j in itertools.combinations(range(n), 2): w *= (np.cos(grids[i]) - np.cos(grids[j])) ** 2
    x = sum(2 * np.cos(G) for G in grids)
    w = (w / w.sum()).ravel(); x = x.ravel()
    return x, w
def ks(sample, x, wt):
    o = np.argsort(x); xs, cw = x[o], np.cumsum(wt[o])
    s = np.sort(sample); F = np.interp(s, xs, cw)
    hi = np.arange(1, len(s) + 1) / len(s); lo = np.arange(len(s)) / len(s)
    return max(np.max(np.abs(hi - F)), np.max(np.abs(lo - F)))
laws = {'USp(4)': law(2, 1500), 'USp(6)': law(3, 160), 'USp(8)': law(4, 56)}
g = np.linspace(-8, 8, 400001); gw = np.exp(-g ** 2 / 2); gw /= gw.sum()
laws['N(0,1)'] = (g, gw)
for nm, (x, w) in laws.items():
    out.append(f"{nm:7s}: E2={np.sum(w*x**2):.3f} E4={np.sum(w*x**4):.3f} E6={np.sum(w*x**6):.2f} E8={np.sum(w*x**8):.1f} "
               f"P(|x|<0.2)={np.sum(w*(np.abs(x)<0.2)):.4f} P(|x|>2)={np.sum(w*(np.abs(x)>2)):.4f} P(|x|>3)={np.sum(w*(np.abs(x)>3)):.4f}")
out.append(f"KS distances of m1..m4 ({len(m)} primes, fit on [307, 9999]); 5% critical ~ {1.36/np.sqrt(len(m)):.3f}:")
for j in range(4):
    out.append(f"  m{j+1}: " + ', '.join(f"{nm}: {ks(m[:, j], x, w):.3f}" for nm, (x, w) in laws.items()))
print('\n'.join(out))
open(os.path.join(HERE, 'usp_rank_references.out'), 'w').write('\n'.join(out) + '\n')
