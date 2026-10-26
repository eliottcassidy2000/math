# Are any of the four new motives quadratic twists of each other?  For normalised new parts m_i, m_j and a quadratic
# character e, E[e(p) m_i m_j] is near +-1 if m_i = +-e m_j and near 0 otherwise.  Also robustness of the rank bound:
# max |m_i| over primes >= 2000 (where the poorly fitted low-order base terms contribute < 0.05).
import numpy as np
from sympy import jacobi_symbol
exec(open('moments_new.py').read().split('print(f"{len(qs)} primes in')[0])
qi = qs.astype(int)
chars = {'1': np.ones(len(qi)), '(2/p)': np.array([jacobi_symbol(2, int(p)) for p in qi], dtype=float),
         '(3/p)': np.array([jacobi_symbol(3, int(p)) for p in qi], dtype=float), '(-3/p)': np.array([jacobi_symbol(-3 % int(p), int(p)) for p in qi], dtype=float),
         '(5/p)': np.array([jacobi_symbol(5, int(p)) for p in qi], dtype=float), '(7/p)': np.array([jacobi_symbol(7, int(p)) for p in qi], dtype=float)}
print(f"{len(qs)} primes in [{int(qs.min())}, {int(qs.max())}]")
nm = ['m1', 'm2', 'm3', 'm4']
for e, ev in chars.items():
    M = np.array([[np.mean(ev * m[:, i] * m[:, j]) for j in range(4)] for i in range(4)])
    off = max(abs(M[i, j]) for i in range(4) for j in range(4) if i != j)
    print(f"  twist {e:6s}: largest off-diagonal |E[e m_i m_j]| = {off:.3f}; diagonal {np.round(np.diag(M), 3)}")
big = qs >= 2000
print("max |m_i| over primes >= 2000:", np.round(np.max(np.abs(m[big]), axis=0), 3), f"({big.sum()} primes)")
print("fraction of primes >= 2000 with |m_i| > 2:", np.round(np.mean(np.abs(m[big]) > 2, axis=0), 3))
