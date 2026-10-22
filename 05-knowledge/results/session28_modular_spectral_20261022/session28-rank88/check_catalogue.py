# Exact checks remain active under python -O. See REPORT.md for proofs and provenance.
from pathlib import Path
import re, json, hashlib
from fractions import Fraction as Q
p = Path(__file__).parent
G = json.loads((p / 'gram88.json').read_text())
source_hash = 'ab2f5fa56e2085534cbce64d5afbc95c422614bf7b12afec17c0d6fd1e239369'
if not (len(G) == 88 and all((len(row) == 88 for row in G))):
    raise AssertionError('len(G) == 88 and all((len(row) == 88 for row in G))')
if not all((G[i][j] == G[j][i] for i in range(88) for j in range(88))):
    raise AssertionError('all((G[i][j] == G[j][i] for i in range(88) for j in range(88)))')
if not all((G[i][i] % 2 == 0 for i in range(88))):
    raise AssertionError('all((G[i][i] % 2 == 0 for i in range(88)))')
A = [[Q(x) for x in r] for r in G]
piv = []
for i in range(88):
    d = A[i][i]
    if not d > 0:
        raise AssertionError('d > 0')
    piv.append(d)
    for j in range(i + 1, 88):
        for k in range(j, 88):
            A[j][k] -= A[i][j] * A[i][k] / d
            A[k][j] = A[j][k]
det = Q(1)
for d in piv:
    det *= d
if not det == 1:
    raise AssertionError('det == 1')
out = {'source': 'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88b.html', 'source_sha256': source_hash, 'rank': 88, 'symmetric': True, 'even': True, 'positive_LDL_pivots': list(map(str, piv)), 'determinant': str(det), 'minimum_status': '8 is CITED from Sakamoto/catalogue, not independently certified by these matrix checks'}
(p / 'gram88.json').write_text(json.dumps(G) + '\n')
(p / 'catalogue_check.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'positive_LDL_pivots'}, indent=2))
