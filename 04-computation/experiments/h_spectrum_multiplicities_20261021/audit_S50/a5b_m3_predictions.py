"""Audit S50 - a5b: the (M3) 'predictions' of h_multiplicities.out (n = 10..12) against the full sequence law.
h_multiplicities.py builds its tuple series with nmax = 9 (cores of total size <= 9), so cores needing 10+ vertices
(e.g. C_3 C_3 R_4 for 45, C_3 R_4 R_4 for 75, four C_3's for 81, C_3 + a 35-atom for 105, two 11-atoms for 121)
are dropped.  Session values are copied from h_multiplicities.out (data)."""
import json, os
from collections import defaultdict, Counter
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
inv = json.load(open(os.path.join(HERE, 'inventory_375.json')))
strong = defaultdict(Counter)
for k, rows in inv['orders'].items():
    for c, h, sc, aut in rows:
        if h < 375:
            strong[int(k)][h] += 1
def tuples(strong, hmax, nmax):
    B = defaultdict(int); B[(0, 0, 1)] = 1; fr = {(0, 0, 1): 1}
    while fr:
        nx = defaultdict(int)
        for (k, t, h), c in fr.items():
            for m, d in strong.items():
                if t + m <= nmax:
                    for h1, c1 in d.items():
                        if h * h1 <= hmax:
                            nx[(k + 1, t + m, h * h1)] += c * c1
        for key, c in nx.items():
            B[key] += c
        fr = nx
    return B
def N(B, n, h):
    return sum(c * comb(n - t + k, k) for (k, t, hh), c in B.items() if hh == h and n >= t)
full = tuples(strong, 374, 40); trunc = tuples(strong, 374, 9)
sess = {3: [8, 9, 10], 5: [7, 8, 9], 9: [27, 35, 44], 11: [6, 7, 8], 13: [6, 7, 8], 15: [42, 56, 72], 25: [20, 27, 35],
        27: [50, 77, 112], 35: [8, 10, 12], 45: [93, 145, 209], 49: [31, 39, 47], 63: [6, 8, 10], 75: [64, 100, 144],
        81: [96, 134, 176], 105: [55, 73, 91], 121: [36, 48, 60]}
bad = []
for h, v in sess.items():
    true = [N(full, n, h) for n in (10, 11, 12)]
    tr = [N(trunc, n, h) for n in (10, 11, 12)]
    assert tr == v, (h, tr, v)          # the session's numbers are exactly the truncated series
    if true != v:
        bad.append((h, v, true))
print('session (M3) predictions that are WRONG (h, session N(10..12), true N(10..12)):')
for row in bad:
    print('    ', row)
print('all other (M3) predictions correct:', [h for h in sess if h not in [b[0] for b in bad]])
