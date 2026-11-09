"""For each order-12 escaper: its minimum cuts (3 back arcs) with both sides of size >= 2, and whether both sides are
first creators (H(T[S]) not in S(|S|-1)); also whether some side is a first creator.  (genesis lane, EXPLORATORY
structure census of the escapers found; exact per class.)"""
import sys
from collections import Counter
from verify_escapers import adj_from, ham
from spectra import in_S

def fc(out, S, k):
    if k <= 1:
        return True          # a single vertex is the first creator of 1
    s = in_S(k - 1, ham(out, S))
    return s is False

stats = Counter(); examples = {}
for line in open(sys.argv[1]):
    n, code = line.split()[:2]; n = int(n)
    out = adj_from(code, n); full = (1 << n) - 1
    H = ham(out, full)
    kinds = set()
    for S in range(1, full):
        k = bin(S).count('1')
        if k < 2 or k > n - 2:
            continue
        c = sum(bin(out[v] & ~S & full).count('1') for v in range(n) if S >> v & 1)   # arcs S -> rest
        if c != 3:
            continue
        # S is the 'later' block (3 arcs go back from S to the rest); rest is the earlier block
        a, b = fc(out, full & ~S, n - k), fc(out, S, k)
        kinds.add((min(k, n - k), a and b, a or b))
    best = max(kinds, key=lambda t: (t[1], t[2], t[0])) if kinds else None
    key = 'no nontrivial 3-cut' if best is None else ('both sides first creators' if best[1] else ('one side first creator' if best[2] else 'neither side first creator'))
    stats[key] += 1
    examples.setdefault(key, (H, sorted(kinds)))
for k, v in stats.most_common():
    print(v, k, ' e.g. H=%d cuts(min side, both FC, some FC)=%s' % (examples[k][0], examples[k][1][:4]))
