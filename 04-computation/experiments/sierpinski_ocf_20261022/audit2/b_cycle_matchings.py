# Audit2 (B): the cycle-gluing count of the one-loop proof.  p labelled runs (segments with ends (j,0), (j,1)); a perfect
# matching of the 2p ends glues them into a 2-regular multigraph on the runs; count matchings giving ONE cycle, and
# compare with (p-1)! 2^(p-1).  Also count, among them, the orientation pattern relative to run 1 (each run forward or
# backward), which must be uniform: (p-1)! matchings per orientation pattern of runs 2..p.
from math import factorial
from collections import Counter
def matchings(items):
    if not items: yield []; return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        for m in matchings(items[1:i] + items[i + 1:]): yield [(a, b)] + m
for p in range(1, 8):
    ends = [(j, e) for j in range(p) for e in (0, 1)]
    single = 0; pattern = Counter()
    for m in matchings(ends):
        mate = {}
        for a, b in m: mate[a] = b; mate[b] = a
        # walk: start at run 0, leave through end (0,1)
        seen = {0}; cur = (0, 1); orient = {0: +1}; ok = True
        while True:
            nxt = mate[cur]
            j, e = nxt
            if j == 0:
                ok = (len(seen) == p) and e == 0
                break
            if j in seen: ok = False; break
            seen.add(j); orient[j] = +1 if e == 0 else -1     # entered at its start: traversed forward
            cur = (j, 1 - e)
        if ok:
            single += 1; pattern[tuple(orient[j] for j in range(1, p))] += 1
    uniform = len(set(pattern.values())) == 1 and len(pattern) == 2 ** (p - 1)
    print(f"p={p}: single-cycle matchings {single}, (p-1)! 2^(p-1) = {factorial(p - 1) * 2 ** (p - 1)}, "
          f"{'match' if single == factorial(p - 1) * 2 ** (p - 1) else 'MISMATCH'}; orientation patterns uniform: {uniform}")
