# Exact checks remain active under python -O. See REPORT.md for proofs and provenance.
import sympy as s, json
from pathlib import Path
x = s.Symbol('x')
P = x ** 3 + 2209944 * x ** 2 - 15663522502656 * x - 19976984434430705664
p = s.Poly(P, x)
a, b, c = p.all_coeffs()[1:]
C = s.Matrix([[0, 0, -c], [1, 0, -b], [0, 1, -a]])
G = s.Matrix(3, 3, lambda i, j: s.trace(C ** (i + j)))
mins = [G[:i, :i].det() for i in [1, 2, 3]]
if not all((m > 0 for m in mins)):
    raise AssertionError('all((m > 0 for m in mins))')
if not G * C == C.T * G:
    raise AssertionError('G * C == C.T * G')
if not G.det() == s.discriminant(P, x) > 0:
    raise AssertionError('G.det() == s.discriminant(P, x) > 0')
m = -s.diff(P, x) / (3 * P)
if not s.cancel(m - s.trace((C - x * s.eye(3)).inv()) / 3) == 0:
    raise AssertionError('s.cancel(m - s.trace((C - x * s.eye(3)).inv()) / 3) == 0')
out = {'cusp_companion': [[str(v) for v in row] for row in C.tolist()], 'positive_trace_gram': [[str(v) for v in row] for row in G.tolist()], 'Sylvester_leading_minors': list(map(str, mins)), 'GC_equals_CtG': True, 'rational_resolvent': str(m), 'spectrum_scope': '3 real eigenvalues of T2 on S44; not zeros of zeta or singularities of theta in q'}
print(json.dumps(out, indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
