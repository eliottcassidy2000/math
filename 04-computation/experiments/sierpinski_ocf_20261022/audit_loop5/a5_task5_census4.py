# Audit (HYP-9371 section 8), task 5: the block-counting argument needs the <= 4-vertex census of a doubly regular
# tournament to depend on n only.  Direct count of the four 4-vertex types (by score sequence) on the 14 DRT(23)s.
import os, sys, json, itertools
from collections import Counter
import numpy as np
D = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, D)
from drt_switch import code_to_S
cls = json.load(open(os.path.join(D, 'drt23_classes.json'))); dbl = json.load(open(os.path.join(D, 'drt23_doubled.json')))
res = {}
for lab, code in [(c['label'], c['code']) for c in cls] + [('D23_dbl', dbl['code'])]:
    A = (code_to_S(code) == 1).astype(int)
    cnt = Counter()
    for q in itertools.combinations(range(23), 4):
        sub = A[np.ix_(q, q)]; cnt[tuple(sorted(sub.sum(1)))] += 1
    c3 = Counter(tuple(sorted(A[np.ix_(t, t)].sum(1))) for t in itertools.combinations(range(23), 3))
    res[lab] = (dict(sorted(cnt.items())), dict(sorted(c3.items())))
vals = set(repr(v) for v in res.values())
print(f"14 DRT(23) tournaments: distinct (4-census, 3-census) vectors: {len(vals)}; census = {next(iter(res.values()))}")
