"""(i) closed form for the average-TV profile of P_K (auditor's derivation): avgTV(n) = P(B > n) + P(B = n)/2,
B ~ Bin(K-2, 1/2), checked exactly against the class computation (independent Phi) for K = 2..18;
(ii) the circulant of the note's convention (weight w(j) at offset j from px+1) has eigenvalues {1/(2 eta - 1)}."""
import math
from fractions import Fraction
import numpy as np


def v2(n):
    return (n & -n).bit_length() - 1


def phi(c, m, p=3):
    t = (p * c + 1) % (1 << m)
    if t == 0:
        return 1, 1
    v = v2(t)
    if v >= m:
        return 1, 1
    m2 = m - v
    return (t >> v) % (1 << m2), m2


for K in range(2, 19):
    n = 1 << (K - 1)
    cls = [(2 * i + 1, K) for i in range(n)]
    for s in range(K):
        avg = sum((1 - Fraction(1, 1 << (m - 1))) for c, m in cls) / n
        B = lambda k: Fraction(math.comb(K - 2, k), 1 << (K - 2)) if 0 <= k <= K - 2 else Fraction(0)
        pred = sum(B(k) for k in range(s + 1, K - 1)) + B(s) / 2
        assert avg == pred, (K, s, avg, pred)
        cls = [phi(c, m) for c, m in cls]
print('avg TV(n) = P(B > n) + P(B = n)/2, B ~ Bin(K-2, 1/2): exact for K = 2..18, all n; equals 1/2 at n = (K-2)/2')
for L in (2, 6, 18, 54):
    w = np.array([2.0 ** (-j) for j in range(1, L + 1)]) / (1 - 2.0 ** (-L))
    r = np.zeros(L)
    for j in range(1, L + 1):
        r[j % L] = w[j - 1]                       # weight w(j) at offset j
    C = np.array([np.roll(r, i) for i in range(L)])
    ev = np.linalg.eigvals(C)
    eta = np.exp(2j * np.pi * np.arange(L) / L)
    pred = 1 / (2 * eta - 1)
    d = max(min(abs(e - q) for q in pred) for e in ev)
    print(f'L = {L}: eigenvalues of Circ(w) (note convention) = {{1/(2 eta - 1)}} up to {d:.1e}; '
          f'max | |z - 1/3| - 2/3 | = {max(abs(abs(e - 1/3) - 2/3) for e in ev):.1e}; min |z| = {min(abs(ev)):.4f}')
