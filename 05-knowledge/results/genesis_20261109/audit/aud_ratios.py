"""AUDIT G1: exact ratios R_H, R_tau, and the THM-589 CV^2 identity (own code)."""
from fractions import Fraction as F
from math import factorial, e
import itertools
maxH = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661, 9: 3357, 10: 15745, 11: 95095}
S = {n: set(int(x) for x in open('S_%d.txt' % n)) for n in range(1, 12)}
assert all(max(S[n]) == maxH[n] for n in range(1, 12))
EH = lambda n: F(factorial(n), 2 ** (n - 1))
print('R_H(n) = max H(n-1)/E[H(n)] (exact maxima from census, n <= 12):')
print('  ', ', '.join('%d: %.4f' % (n, float(maxH[n - 1] / EH(n))) for n in range(5, 13)))
r = [float(maxH[n - 1] / EH(n)) for n in range(5, 13)]
print('   monotone decreasing?', all(r[i] > r[i + 1] for i in range(len(r) - 1)))
lb12, lb13 = 531205, 3719831          # THM-329 search values (not proved maximal)
print('R_H(13): with 545,835 (session) %.4f; with THM-329 lower bound 531,205: >= %.4f; elementary upper bound 12*95095=%d: <= %.4f'
      % (545835 / EH(13), lb12 / EH(13), 12 * 95095, 12 * 95095 / EH(13)))
print('R_H(14): circulant 3,711,175 -> %.5f; THM-329 3,719,831 -> %.5f; elementary upper bound 13*12*95095=%d -> %.4f'
      % (3711175 / EH(14), lb13 / EH(14), 13 * 12 * 95095, 13 * 12 * 95095 / EH(14)))
print('R_tau(m) = 2^(m-1)(m-1)^(m-3)/m^(m-2) vs 2^(m-1)/(e m):')
ok = True
for m in range(3, 60):
    R = F(2 ** (m - 1) * (m - 1) ** (m - 3), m ** (m - 2))
    ok &= float(R) >= 2 ** (m - 1) / (e * m)
    if m <= 9: print('   m=%d R_tau=%.4f lower bound %.4f' % (m, float(R), 2 ** (m - 1) / (e * m)))
print('   bound holds m=3..59:', ok)
def W(n):
    # c[r][k] = sum over compositions of r into k odd parts of 2^{#parts>=3}
    c = [[0] * (n + 1) for _ in range(n + 1)]
    c[0][0] = 1
    for r_ in range(1, n + 1):
        for k in range(1, r_ + 1):
            c[r_][k] = sum((1 if p == 1 else 2) * c[r_ - p][k - 1] for p in range(1, r_ + 1, 2))
    return sum(factorial(k) * c[n][k] for k in range(n + 1))
print('W(n), n=1..8:', [W(n) for n in range(1, 9)], '(THM-589: 1,2,8,32,158,928,6350,49752)')
for n in range(3, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    hs = []
    for msk in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if msk >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
        hs.append(sum(1 for p in itertools.permutations(range(n)) if all(A[p[k]][p[k + 1]] for k in range(n - 1))))
    m1 = F(sum(hs), len(hs)); m2 = F(sum(h * h for h in hs), len(hs))
    print('   n=%d: E H = %s (= n!/2^(n-1): %s), CV^2 = %s, W(n)/n! - 1 = %s' % (n, m1, m1 == EH(n), m2 / m1 ** 2 - 1, F(W(n), factorial(n)) - 1))
print('n*CV^2 = n(W(n)/n! - 1) for n = 10, 20, 40, 80, 160:', [round(float(n * (F(W(n), factorial(n)) - 1)), 4) for n in (10, 20, 40, 80, 160)])
