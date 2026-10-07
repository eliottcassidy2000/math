import re
from fractions import Fraction
from mpmath import mp, mpf, log
mp.dps = 50
txt = open('zk_mitm_15_25.out').read()
Z = {0: 1}
for m in re.finditer(r'DFS Z_(\d+) = (\d+)', txt): Z[int(m.group(1))] = int(m.group(2))
for m in re.finditer(r'Z_(\d+) = (\d+)\s*$', txt, re.M): Z[int(m.group(1))] = int(m.group(2))
Bmax = {}; Bmin = {}
for m in re.finditer(r'L=(\d+) B_L\[0\]=Z_L=(\d+) maxB=(\d+) \(at s=(\d+)\) minB=(\d+) \(at s=(\d+)\)', txt):
    L = int(m.group(1)); Bmax[L] = (int(m.group(3)), int(m.group(4))); Bmin[L] = (int(m.group(5)), int(m.group(6)))
K = max(Z)
# OEIS A181610 b-file check
oeis = {}
for line in open('A181610_b.txt') if False else []: pass
print('K range 0..%d' % K)
D = {k: 2*Z[k+1] - 9*Z[k] for k in range(K)}
print('Delta_k = O_k - E_k = 2 Z_{k+1} - 9 Z_k:')
print([D[k] for k in range(K)])
# exact partial sums for c
c_part = Fraction(1)
rows = []
for k in range(K):
    c_part += Fraction(D[k] * 2**k, 9**(k+1))
    # check: equals Z_{k+1} (2/9)^{k+1}
    assert c_part == Fraction(Z[k+1] * 2**(k+1), 9**(k+1))
print('identity Z_{k+1}(2/9)^{k+1} = 1 + sum_{j<=k} Delta_j 2^j/9^{j+1} verified exactly for k<%d' % K)
for k in [10, 20, 25, 30, 35, 38, 39, 40]:
    print('k=%d  Z_k (2/9)^k = %s' % (k, mp.nstr(mpf(Z[k]) * (mpf(2)/9)**k, 30)))
# growth of |Delta|
import math
print('log|Delta_k|/k and |Delta_k|^(1/k):')
for k in range(5, K):
    if D[k] != 0: print(k, D[k], '%.4f' % (abs(D[k])**(1.0/k)), ' Delta_k/sqrt(Z_k)=%.3e' % (D[k]/math.sqrt(Z[k])))
# density sup/inf
print('density f_L = 2^L B_L / 9^L : sup, inf, sup/inf, rigorous growth bounds (minB)^(1/L), (maxB)^(1/L)')
for L in sorted(Bmax):
    fM = mpf(Bmax[L][0]) * 2**L / mpf(9)**L; fm = mpf(Bmin[L][0]) * 2**L / mpf(9)**L
    lo = mpf(Bmin[L][0])**(mpf(1)/L); hi = mpf(Bmax[L][0])**(mpf(1)/L)
    print(L, mp.nstr(fM, 10), mp.nstr(fm, 10), mp.nstr(fM/fm, 8), ' growth in [%s, %s]  dim_5 in [%s, %s]' % (mp.nstr(lo, 8), mp.nstr(hi, 8), mp.nstr(log(lo)/log(5), 8), mp.nstr(log(hi)/log(5), 8)))
print('log_5(9/2) =', mp.nstr(log(mpf(9)/2)/log(5), 20))
import json
json.dump({'Z': {str(k): str(v) for k, v in Z.items()}, 'Delta': {str(k): v for k, v in D.items()}}, open('zk_values.json', 'w'), indent=0)
