"""AUDIT G1: H of 6-tournaments with minimum out-degree >= 2 (the halves of the 3^6 8^6 escapers), and the check that
the 106 classes with an FC-FC 3-back-arc cut are exactly the classes with scores 3^6 8^6."""
from collections import Counter, defaultdict
from aud_struct import ham, decode, cut_census
n = 6; pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
hs = defaultdict(Counter)
for m in range(1 << 15):
    out = [0] * n
    for b, (i, j) in enumerate(pairs):
        if m >> b & 1: out[i] |= 1 << j
        else: out[j] |= 1 << i
    sc = tuple(sorted(bin(o).count('1') for o in out))
    if min(sc) >= 2:
        hs[sc][ham(out, list(range(6)))] += 1
for sc in sorted(hs): print('scores', sc, 'H values (labelled counts):', sorted(hs[sc].items()))
rows = [l.split() for l in open('classes193_props.txt')]
both = sc_only = fc_only = 0
for r in rows:
    out = decode(r[0], 12)
    a = any(c[2] and c[3] for c in cut_census(out, 12))
    s = (r[5] == 'sc=333333888888')
    both += a and s; sc_only += s and not a; fc_only += a and not s
print('FC-FC 3-cut and scores 3^6 8^6:', both, ' scores only:', sc_only, ' FC-FC cut only:', fc_only)
