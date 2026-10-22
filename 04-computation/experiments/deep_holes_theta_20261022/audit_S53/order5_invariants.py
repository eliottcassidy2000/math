"""Audit check (S53): orbit invariants of the two order-5 strict-Fricke classes versus the true c''(A3^8), and the
'geometric form' of HYP-9368 (A): for a strict Fricke point v of order n whose unique minimal vector is v itself
(nearest Leech point 0), c = v + v/(n-1) = n v/(n-1) should be a deep hole of Coxeter number n - 1 (A3^8 for n = 5).
Invariant: the distribution of y.lambda over the 196560 minimal vectors lambda, y = 5v."""
import numpy as np
from collections import Counter
from cypari import pari
from audit_common import ROWS, GLAM_I, GLAM_P, h_plus_one_point, coset_profile, shell_vecs, GLAM_F

mv = pari.qfminim(GLAM_P, 4, 200000, 0)
L4 = np.array([[int(mv[2][i, j]) for i in range(24)] for j in range(196560 // 2)], dtype=np.int64)
L4 = np.vstack([L4, -L4])
GL4 = L4 @ GLAM_I


def normalise(Ym, e):
    """translate c = Ym/e so that its unique minimal vector is v itself; return y = e v (integer coords)."""
    vals, vecs = shell_vecs(GLAM_F, GLAM_I, np.ascontiguousarray(Ym, dtype=np.int64), e, 2.0 + 1e-6)
    mn = min(int(s) for s in vals)
    a = [np.asarray(v_, dtype=np.int64) for s, v_ in zip(vals, vecs) if int(s) == mn]
    assert len(a) == 1
    return e * a[0] + Ym          # e (a + Ym/e)


def invariant(y):
    return tuple(sorted(Counter((GL4 @ y).tolist()).items()))


row = [r for r in ROWS if r[0] == 'a3^8'][0]
Ym, m = h_plus_one_point(row[3], row[1])
y_c = normalise(Ym, m)
assert int(y_c @ GLAM_I @ y_c) == 40
EX = {
    "c''(A3^8)": y_c,
    'census class n2=24': np.array([-7, -4, 20, 5, -21, 8, 1, 1, 1, 1, -5, -8, -4, -3, -2, 3, 0, -4, 1, 2, -6, -2, -4, -1]),
    'census class n2=23': np.array([-9, -2, 23, 6, -22, 7, 0, 1, 1, 2, -6, -8, -5, -4, -3, 2, 0, -3, 1, 1, -6, -2, -4, 0]),
}
inv_c = invariant(y_c)
for name, y in EX.items():
    y = np.asarray(y, dtype=np.int64)
    inv = invariant(y)
    # geometric form: c = 5v/4 = y/4
    prof = coset_profile(-y, 4, 2.0)
    mn = min(prof)
    deep = (mn == 2 * 16)
    nvert = prof.get(32, 0)
    print(f'{name:20s}: y.lambda distribution = {dict(inv)}')
    print(f'{"":20s}  same invariant as c\'\'(A3^8): {inv == inv_c}; y/4: min norm {mn}/16, '
          f'{"deep hole with " + str(nvert) + " vertices (A3^8 has 32)" if deep else "NOT a deep hole"}')
