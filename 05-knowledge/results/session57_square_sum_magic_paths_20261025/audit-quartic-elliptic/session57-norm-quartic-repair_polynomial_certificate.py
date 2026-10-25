"""Direct audit of printed (5.1), arXiv:2508.07806v1, no rational parametrization."""
import json
import sympy as S
from pathlib import Path
u, v, x, y = S.symbols('u v x y')
V = x * x + y * y
W = x * x - y * y
H = x ** 4 - 6 * x * x * y * y + y ** 4
A = [H ** 2 * W ** 2 * (u * u + v * v) ** 2, H ** 2 * V ** 2 * (u * u - v * v) ** 2, W ** 2 * (V ** 2 * u * u - 16 * u * v * x * x * y * y - V ** 2 * v * v) ** 2, W ** 2 * (V ** 2 * u * u + 16 * u * v * x * x * y * y - V ** 2 * v * v) ** 2]
B = [S.Integer(0), H ** 2 * ((y - x) * u + (x + y) * v) * ((x + y) * u + (x - y) * v) * ((x - y) * u + (x + y) * v) * ((x + y) * u - (x - y) * v), -W ** 2 * (V * u + 4 * u * x * y - V * v) * (V * u - 4 * u * x * y - V * v) * (V * u + 4 * v * x * y + V * v) * (V * u - 4 * v * x * y + V * v), -W ** 2 * (V * u + 4 * u * x * y + V * v) * (V * u - 4 * u * x * y + V * v) * (V * u + 4 * v * x * y - V * v) * (V * u - 4 * v * x * y - V * v)]
E = u * u - v * v
F = u * u + v * v
T = x * x * y * y
compact = [[H * W * F, 2 * u * v * V * H, 2 * W * (4 * T * E + u * v * V ** 2), 2 * W * (4 * T * E - u * v * V ** 2)], [H * V * E, 2 * x * y * F * H, 2 * x * y * (V ** 2 * E + 4 * u * v * W ** 2), 2 * x * y * (V ** 2 * E - 4 * u * v * W ** 2)], [W * (V ** 2 * E - 16 * u * v * T), 4 * x * y * (W ** 2 * E - u * v * V ** 2), 4 * x * y * F * W * V, 4 * x * y * W * V * (E - 2 * u * v)], [W * (V ** 2 * E + 16 * u * v * T), 4 * x * y * (W ** 2 * E + u * v * V ** 2), 4 * x * y * W * V * (E + 2 * u * v), 4 * x * y * F * W * V]]
if not all((S.Poly(S.expand(A[i] + B[j] - compact[i][j] ** 2), u, v, x, y).is_zero for i in range(4) for j in range(4))):
    raise RuntimeError('verification failed at source line 20')
roots = []
for i in range(4):
    row = []
    for j in range(4):
        coefficient, factors = S.factor_list(A[i] + B[j], u, v, x, y)
        if not S.sqrt(coefficient).is_Integer:
            raise RuntimeError('verification failed at source line 26')
        if not all((e % 2 == 0 for f, e in factors)):
            raise RuntimeError('verification failed at source line 27')
        root = S.sqrt(coefficient) * S.prod((f ** (e // 2) for f, e in factors))
        if not S.Poly(S.expand(A[i] + B[j] - root ** 2), u, v, x, y).is_zero:
            raise RuntimeError('verification failed at source line 29')
        if not S.Poly(root, u, v, x, y).domain == S.ZZ:
            raise RuntimeError('verification failed at source line 30')
        row.append(str(root))
    roots.append(row)
print(json.dumps({'source': 'https://arxiv.org/pdf/2508.07806v1', 'equation': '5.1', 'parameter_map': {'u': 'u1', 'v': 'u2', 'x': 'v1', 'y': 'v2'}, 'identities_verified': 16, 'root_domain': 'ZZ[u,v,x,y]', 'all_integer_specializations_including_boundaries': True, 'roots': roots}, indent=2))
