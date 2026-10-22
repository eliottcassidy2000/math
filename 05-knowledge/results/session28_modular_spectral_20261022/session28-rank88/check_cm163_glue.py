# Exact checks remain active under python -O. See REPORT.md for proofs and provenance.
import sympy as s, json
from sympy.matrices.normalforms import hermite_normal_form
from pathlib import Path
p = 163
A = s.Matrix([[2, 1], [1, 82]])
G = s.diag(A, A, A, A)
v = [s.Matrix([1, 0, 1, 18]), s.Matrix([0, 1, -18, 1])]
if not all((x.dot(y) % p == 0 for x in v for y in v)):
    raise AssertionError('all((x.dot(y) % p == 0 for x in v for y in v))')
g = [s.Matrix([t for a in x for t in (82 * a, -a)]) for x in v]
Y = hermite_normal_form(s.Matrix.hstack(p * s.eye(8), *g))
B = Y / p
H = B.T * G * B
if not all((t.q == 1 for t in H)):
    raise AssertionError('all((t.q == 1 for t in H))')
if not all((H[i, i] % 2 == 0 for i in range(8))):
    raise AssertionError('all((H[i, i] % 2 == 0 for i in range(8)))')
if not H.det() == 1:
    raise AssertionError('H.det() == 1')
if not B.det() == s.Rational(1, p * p):
    raise AssertionError('B.det() == s.Rational(1, p * p)')
M = s.diag(*[s.Matrix([[0, -41], [1, 1]])] * 4)
alpha = B.inv() * M * B
if not all((t.q == 1 for t in alpha)):
    raise AssertionError('all((t.q == 1 for t in alpha))')
if not alpha * alpha - alpha + 41 * s.eye(8) == s.zeros(8):
    raise AssertionError('alpha * alpha - alpha + 41 * s.eye(8) == s.zeros(8)')
if not alpha.T * H + H * alpha == H:
    raise AssertionError('alpha.T * H + H * alpha == H')
if not alpha.T * H * alpha == 41 * H:
    raise AssertionError('alpha.T * H * alpha == 41 * H')
out = {'binary_gram': A.tolist(), 'discriminant': p, 'isotropic_generators': [list(x) for x in v], 'isotropic_norms': [str(x.dot(x)) for x in v], 'mutual_pairing': str(v[0].dot(v[1])), 'index': p * p, 'overlattice_basis_numerator': Y.tolist(), 'basis_denominator': p, 'overlattice_gram': H.tolist(), 'determinant': 1, 'alpha_integral_matrix': alpha.tolist(), 'alpha_minpoly': 'x^2-x+41', 'alpha_adjoint_relation': 'alpha*+alpha=I', 'alpha_similarity': 'alpha*alpha=41I', 'identification': 'Even positive-definite unimodular rank8, hence E8 by classical classification; explicit CM endomorphism retained'}

def conv(x):
    return int(x) if isinstance(x, s.Integer) else x
print(json.dumps(out, indent=2, default=conv))
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2, default=conv) + '\n')
