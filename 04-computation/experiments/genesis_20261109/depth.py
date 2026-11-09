"""Creator depth of tournaments (genesis lane): the least d such that deleting some d vertices leaves a first
creator (an m-tournament U with H(U) not in S(m - 1)).  For order-12 escapers d >= 2 by definition; we test d = 2
(all 66 pairs; exact S(9) from klein's census) and report the depth-2 status of every class given on stdin."""
import sys, itertools
from verify_escapers import adj_from, ham
from spectra import in_S

def depth2(code, n):
    out = adj_from(code, n)
    full = (1 << n) - 1
    hits = 0
    for u, w in itertools.combinations(range(n), 2):
        h = ham(out, full & ~(1 << u) & ~(1 << w))
        s = in_S(n - 3, h)
        if s is None:
            return None
        if not s:
            hits += 1
    return hits

if __name__ == '__main__':
    from collections import Counter
    c = Counter(); rows = []
    for line in sys.stdin:
        p = line.split()
        if not p: continue
        n, code = int(p[0]), p[1]
        d = depth2(code, n)
        c[d is not None and d > 0] += 1
        rows.append(d)
    print('classes with a first creator among their (n-2)-subtournaments:', c[True], ' without:', c[False])
    print('number of first-creator (n-2)-subtournaments per class: min', min(rows), 'max', max(rows))
