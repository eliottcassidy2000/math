"""Independent exact audit of arXiv:2508.07806v1, section 5.1.
Source PDF SHA256: 964fdfc1046e0f34e652f97ce507e803bf50581e7a6042bb0c461d7a5bd3ad66
Requires SymPy. No numerical tolerance is used.
"""
import json
from math import isqrt
import sympy as S
p, q, r, s, m = S.symbols('p q r s m')
A = [4 * (p * r - q * s) ** 2, 4 * (p * r + q * s) ** 2, (p * p - q * q + r * r - s * s) ** 2, (m * m * p * p - m * m * q * q + r * r - s * s) ** 2 / m ** 2]
B = [S.Integer(0), -4 * (p * p - q * q) * (r * r - s * s), 4 * (p * p - s * s) * (q * q - r * r), 4 * (m * m * p * p - s * s) * (m * m * q * q - r * r) / m ** 2]
H = lambda x, y: x ** 4 - 6 * x * x * y * y + y ** 4
check = lambda f: S.cancel(f) == 0
if not check(A[2] + B[3] - A[3] - B[2] - (m * m - 1) / m ** 2 * (H(r, s) - m * m * H(p, q))):
    raise RuntimeError('verification failed at source line 14')
if not check(A[2] + B[2] - A[3] - B[3] - (m * m - 1) / m ** 2 * ((r * r + s * s) ** 2 - m * m * (p * p + q * q) ** 2)):
    raise RuntimeError('verification failed at source line 15')
R = [[2 * (p * r - q * s), 2 * (p * s - q * r), 2 * (p * q - r * s), 2 * (m * m * p * q - r * s) / m], [2 * (p * r + q * s), 2 * (p * s + q * r), 2 * (p * q + r * s), 2 * (m * m * p * q + r * s) / m], [p * p - q * q + r * r - s * s, p * p - q * q - r * r + s * s, p * p + q * q - r * r - s * s, None], [(m * m * p * p - m * m * q * q + r * r - s * s) / m, (m * m * p * p - m * m * q * q - r * r + s * s) / m, None, (m * m * p * p + m * m * q * q - r * r - s * s) / m]]
if not all((check(A[i] + B[j] - R[i][j] ** 2) for i in range(4) for j in range(4) if R[i][j] is not None)):
    raise RuntimeError('verification failed at source line 20')
control = {p: 1, q: 2, r: 3, s: 4, m: 5}
if not [25 * f.subs(control) for f in [A[2], A[3], B[2], B[3]]] == [2500, 6724, 7500, 3276]:
    raise RuntimeError('verification failed at source line 22')
u, v, x, y = S.symbols('u v x y')
D = x * x - 2 * x * y - y * y
N = x * x + 2 * x * y - y * y
t = N / D
sub = {p: u * x + v * y, q: u * y - v * x, r: t * (u * x - v * y), s: t * (u * y + v * x), m: t * t}
K = 8 * x * y * (x * x - y * y) * (x * x + y * y) / D ** 2
cross_roots = [K * (u * u - 2 * u * v - v * v), K * (u * u + 2 * u * v - v * v)]
if not check((A[2] + B[3]).subs(sub) - cross_roots[0] ** 2):
    raise RuntimeError('verification failed at source line 28')
if not check((A[3] + B[2]).subs(sub) - cross_roots[1] ** 2):
    raise RuntimeError('verification failed at source line 29')
from itertools import product
samples = 0
As = [S.cancel(z.subs(sub)) for z in A]
Bs = [S.cancel(z.subs(sub)) for z in B]
for vals in product(range(-2, 3), repeat=4):
    ss = dict(zip((u, v, x, y), vals))
    if D.subs(ss) == 0 or N.subs(ss) == 0:
        continue
    av = [z.subs(ss) for z in As]
    bv = [z.subs(ss) for z in Bs]
    for i, j in product(range(4), repeat=2):
        z = av[i] + bv[j]
        if not (z >= 0 and isqrt(int(S.numer(z))) ** 2 == S.numer(z) and (isqrt(int(S.denom(z))) ** 2 == S.denom(z))):
            raise RuntimeError('verification failed at source line 42')
    samples += 1
AA = [44782864, 340218025, 1738222864, 2777290000]
BB = [0, 25777136, 1222007600, 1719217136]
roots = [[isqrt(a + b) for b in BB] for a in AA]
if not all((roots[i][j] ** 2 == AA[i] + BB[j] for i, j in product(range(4), repeat=2))):
    raise RuntimeError('verification failed at source line 46')
print(json.dumps({'difference_identities': 'PASS', 'generic_square_identities': 14, 'repaired_cross_square_identities': 2, 'bounded_parameter_samples': samples, 'published_example_roots': roots}, indent=2))
