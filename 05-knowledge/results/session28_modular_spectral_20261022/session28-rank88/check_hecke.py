# Exact checks remain active under python -O. See REPORT.md for proofs and provenance.
import sympy as s, json
from pathlib import Path
N = 48

def mul(a, b):
    return [sum((a[j] * b[i - j] for j in range(i + 1))) for i in range(N + 1)]

def power(a, n):
    r = [1] + [0] * N
    for _ in range(n):
        r = mul(r, a)
    return r

def sigma(n, k):
    return sum((d ** k for d in range(1, n + 1) if n % d == 0))
E4 = [1] + [240 * sigma(n, 3) for n in range(1, N + 1)]
E6 = [1] + [-504 * sigma(n, 5) for n in range(1, N + 1)]
D = [(a - b) // 1728 for a, b in zip(power(E4, 3), power(E6, 2))]
B = [mul(power(D, j), power(E4, 11 - 3 * j)) for j in range(4)]
A = s.Matrix([[b[i] for b in B] for i in range(4)])
v = A.inv() * s.Matrix([1, 0, 0, 0])
f = [sum((v[j] * B[j][i] for j in range(4))) for i in range(N + 1)]

def T(a, p, M):
    return [a[p * n] + (p ** 43 * a[n // p] if n % p == 0 else 0) for n in range(M + 1)]
T2 = s.Matrix.hstack(*[A.inv() * s.Matrix(T(b, 2, 3)) for b in B])
T3 = s.Matrix.hstack(*[A.inv() * s.Matrix(T(b, 3, 3)) for b in B])
if not T2 * T3 == T3 * T2:
    raise AssertionError('T2 * T3 == T3 * T2')
for p, M in [(2, T2), (3, T3)]:
    for j, b in enumerate(B):
        target = T(b, p, N // p)
        if not target == [sum((M[k, j] * B[k][i] for k in range(4))) for i in range(N // p + 1)]:
            raise AssertionError('target == [sum((M[k, j] * B[k][i] for k in range(4))) for i in range(N // p + 1)]')
x = s.Symbol('x')
poly = T2.charpoly(x).as_expr()
cp = s.div(poly, x - (1 + 2 ** 43))[0]
cyc = s.Matrix.hstack(*[T2 ** j * v for j in range(4)])
cusp = T2[1:, 1:]
if not s.factor(cusp.charpoly(x).as_expr() - cp) == 0:
    raise AssertionError('s.factor(cusp.charpoly(x).as_expr() - cp) == 0')
out = {'basis': 'Delta^j E4^(11-3j), j=0..3', 'theta_basis_coefficients': list(map(str, v)), 'theta_q0_to_q12': list(map(str, f[:13])), 'T2': [[str(x) for x in row] for row in T2.tolist()], 'T3': [[str(x) for x in row] for row in T3.tolist()], 'T2_characteristic_factorization': str(s.factor(poly)), 'cusp_cubic': str(cp), 'cusp_irreducible_over_Q': s.Poly(cp, x).is_irreducible, 'cusp_real_root_intervals': str(s.polys.polytools.intervals(cp, eps=s.Rational(1, 1000))), 'cyclic_determinant': str(cyc.det()), 'T2_theta_q0_to_q8': list(map(str, T(f, 2, 8))), 'checks': 'T2/T3 commute; each operator independently matches expansions through q24/q16; cyclic determinant nonzero'}
print(json.dumps(out, indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
