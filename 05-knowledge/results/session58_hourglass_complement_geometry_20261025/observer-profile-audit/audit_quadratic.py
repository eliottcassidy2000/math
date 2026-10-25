def _publication_require(value, message):
    if not value:
        raise RuntimeError(message)
if '--negative-check' in __import__('sys').argv:
    _publication_require(False, 'injected publication guard failure')
import json
from pathlib import Path
from fractions import Fraction as F
P = Path(__file__).resolve().parent.parent / 'observer-profile-algebra'
O = Path(__file__).resolve().parent.parent / 'observer-profile-audit'
a = json.loads((P / 'direct_certificate.json').read_text())
q = json.loads((P / 'quadratic_certificate.json').read_text())
A = a['matrix']
cols = a['columns']
features = q['selected_features']
minor = []
for row in A:
    new = []
    for label in features:
        if label[0] == 'linear':
            new.append(row[cols.index(label[1])])
        elif label[0] == 'square':
            new.append(row[cols.index(label[1])] ** 2)
        else:
            raise RuntimeError('Unexpected feature')
    minor.append(new)
_publication_require(minor == q['selected_integer_minor'], 'independent audit assertion failed')
_publication_require(len(features) == 34 and sum((f[0] == 'linear' for f in features)) == 32, 'independent audit assertion failed')
_publication_require([f for f in features if f[0] == 'square'] == [['square', 3], ['square', 15]], 'independent audit assertion failed')
B = [list(map(F, row)) for row in minor]
D = F(1)
pivots = []
for j in range(len(B)):
    k = next((i for i in range(j, len(B)) if B[i][j]), None)
    _publication_require(k is not None, 'independent audit assertion failed')
    if j != k:
        B[j], B[k] = (B[k], B[j])
        D = -D
    pivot = B[j][j]
    D *= pivot
    pivots.append(str(pivot))
    for i in range(j + 1, len(B)):
        c = B[i][j] / pivot
        if c:
            B[i] = [x - c * y for x, y in zip(B[i], B[j])]
_publication_require(D.denominator == 1 and int(D) == int(q['determinant']), 'independent audit assertion failed')
_publication_require(int(D) % q['prime'] == q['determinant_mod_prime'], 'independent audit assertion failed')
weights = {int(k): v for k, v in a['signed_kernel'].items()}
M = a['positive_mass']
contrasts = []
for mask in [3, 15]:
    col = cols.index(mask)
    raw = sum((weights[r] * row[col] ** 2 for r, row in zip(a['representatives'], A)))
    saved = next((z for z in q['mixture_contrasts'] if z['feature'] == ['square', mask]))
    _publication_require(raw == saved['signed_mixture_numerator'], 'independent audit assertion failed')
    contrasts.append({'mask': mask, 'raw': raw, 'normalized_positive_minus_negative': str(F(raw, M * 720 ** 2))})
pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
labels = {str(mask): [list(e) for k, e in enumerate(pairs) if mask >> k & 1] for mask in [3, 15]}
result = {'rank': 34, 'selected_minor_matches_labeled_features': True, 'determinant': str(D), 'modulus': q['prime'], 'determinant_modulus': int(D) % q['prime'], 'normalization_denominator_power_of_720': 36, 'square_feature_edges': labels, 'mixture_contrasts': contrasts}
(O / 'independent_quadratic_certificate.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
