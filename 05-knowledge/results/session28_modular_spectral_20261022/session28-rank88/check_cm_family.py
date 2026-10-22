# Exact checks remain active under python -O. See REPORT.md for proofs and provenance.
import sympy as s, json
from sympy.matrices.normalforms import hermite_normal_form
from pathlib import Path
rows = []
for p in s.primerange(3, 200):
    if p % 4 != 3:
        continue
    a, b = next(((a, b) for a in range(p) for b in range(p) if (a * a + b * b + 1) % p == 0))
    c = (p + 1) // 4
    A = s.Matrix([[2, 1], [1, 2 * c]])
    G = s.diag(A, A, A, A)
    v = [s.Matrix([1, 0, a, b]), s.Matrix([0, 1, -b, a])]
    g = [s.Matrix([t for x in w for t in (2 * c * x, -x)]) for w in v]
    B = hermite_normal_form(s.Matrix.hstack(p * s.eye(8), *g)) / p
    H = B.T * G * B
    M = s.diag(*[s.Matrix([[0, -c], [1, 1]])] * 4)
    alpha = B.inv() * M * B
    J = 2 * alpha - s.eye(8)
    if not (all((x.q == 1 for x in H)) and all((H[i, i] % 2 == 0 for i in range(8))) and (H.det() == 1)):
        raise AssertionError('all((x.q == 1 for x in H)) and all((H[i, i] % 2 == 0 for i in range(8))) and (H.det() == 1)')
    if not (all((x.q == 1 for x in alpha)) and J * J == -p * s.eye(8) and (J.T * H == -H * J)):
        raise AssertionError('all((x.q == 1 for x in alpha)) and J * J == -p * s.eye(8) and (J.T * H == -H * J)')
    rows.append({'p': p, 'a': a, 'b': b, 'norm_alpha': c, 'index': p * p, 'checks_passed': True})
for p in [3, 7, 11, 19, 23, 43, 67, 163]:
    if not not any(((1 + b * b) % p == 0 for b in range(p))):
        raise AssertionError('not any(((1 + b * b) % p == 0 for b in range(p)))')
a, b, d = s.symbols('a b d', real=True)
U = s.Matrix([[-2, 1], [-2, 0]])
W = s.Matrix([[a, b], [b, d]])
if not s.solve(list(W * U - U.T * W), (a,), dict=True) == [{a: -2 * b - 2 * d}]:
    raise AssertionError('s.solve(list(W * U - U.T * W), (a,), dict=True) == [{a: -2 * b - 2 * d}]')
if not s.expand(W.det().subs(a, -2 * b - 2 * d)) == s.expand(-(b + d) ** 2 - d ** 2):
    raise AssertionError('s.expand(W.det().subs(a, -2 * b - 2 * d)) == s.expand(-(b + d) ** 2 - d ** 2)')
out = {'all_primes_3_mod4_below200': rows, 'count': len(rows), 'two_copy_isotropy_negative_controls': 8, 'U2_symmetric_gram_condition': 'a=-2b-2d, determinant=-(b+d)^2-d^2; no positive definite Gram'}
print(json.dumps(out, indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
