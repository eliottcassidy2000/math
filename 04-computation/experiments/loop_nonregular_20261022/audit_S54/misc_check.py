"""Auditor: small remaining HYP-9370 numbers: min 1/max|x_j| on the n = 8, 10 censuses (vs Session30's bound
x_j^2 <= (n-1)/(2n), i.e. 1/max|x_j| >= sqrt(2n/(n-1)) > sqrt 2), and the one-loop series of the first n = 8 class."""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn
HERE = os.path.dirname(os.path.abspath(__file__))
for fn, n in (('ar8_classes.txt', 8), ('ar10_classes.txt', 10)):
    m = []
    for l in open(os.path.join(HERE, '..', fn)):
        A = tourn.from_code(l.split()[0]); M = (A - A.T) / n
        x = np.abs(np.linalg.eigvals(M).imag)
        m.append(1 / x.max())
    print(f'{fn}: min 1/max|x_j| = {min(m):.4f}; bound sqrt(2n/(n-1)) = {math.sqrt(2 * n / (n - 1)):.4f}')
A = tourn.from_code(open(os.path.join(HERE, '..', 'ar8_classes.txt')).readline().split()[0]); n = 8
x = np.linalg.eigvals((A - A.T) / n).imag
p = lambda k: float(np.sum(x ** k))
# X0 = sum_{m odd} p_{2m} w^m / m ; rho series = exp(X0)
c1 = p(2); c3 = p(6) / 3
print(f'first n = 8 class: rho(z) = 1 + {c1:.4f} w + {c1 ** 2 / 2:.4f} w^2 + {c1 ** 3 / 6 + c3:.4f} w^3 + ...')
