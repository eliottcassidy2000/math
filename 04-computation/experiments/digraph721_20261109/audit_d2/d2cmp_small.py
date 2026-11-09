"""audit D2: compare d2small.out (independent exhaustive enumeration, n <= 7) with the census, per (n, d, parity)."""
import sys
from collections import Counter
mine = Counter()
for line in open(sys.argv[1]):
    p = line.split()
    if p and p[0] == 'C': mine[(int(p[1]), int(p[2]), p[3])] += int(p[4])
cen = Counter()
for line in open(sys.argv[2]):
    p = line.split()
    n = int(p[0])
    if n <= 7: cen[(n, int(p[2]), p[3])] += 1
keys = set(mine) | set(cen)
diff = [(k, mine[k], cen[k]) for k in sorted(keys) if mine[k] != cen[k]]
print('buckets compared:', len(keys), 'classes mine:', sum(mine.values()), 'census:', sum(cen.values()), 'differences:', len(diff))
for x in diff[:20]: print('DIFF', x)
