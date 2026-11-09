"""Which first-creator pairs (|A|, H(A)) + (|B|, H(B)) appear as the two sides of a 3-back-arc cut of an order-12
escaper (A earlier, B later; 3 arcs B -> A).  Counts classes; a class may contribute several pairs.  (genesis lane)"""
import sys
from collections import Counter
from verify_escapers import adj_from, ham
from spectra import in_S
pairs = Counter(); per_class = Counter()
for line in open(sys.argv[1]):
    n, code = line.split()[:2]; n = int(n)
    out = adj_from(code, n); full = (1 << n) - 1
    seen = set()
    for S in range(1, full):
        k = bin(S).count('1')
        if k < 2 or k > n - 2: continue
        c = sum(bin(out[v] & ~S & full).count('1') for v in range(n) if S >> v & 1)
        if c != 3: continue
        A, B = full & ~S, S
        hA, hB = ham(out, A), ham(out, B)
        if in_S(n - k - 1, hA) is False and in_S(k - 1, hB) is False:
            seen.add(((n - k, hA), (k, hB)))
    for p in seen: pairs[p] += 1
    per_class[len(seen)] += 1
for p, v in pairs.most_common(): print(v, 'classes:  A = order %d, H %d   |   B = order %d, H %d' % (p[0][0], p[0][1], p[1][0], p[1][1]))
print('FC-FC cuts per class:', dict(per_class))
