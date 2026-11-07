"""First-creator criticality (kind-pasteur S16): k(T) = #vertices whose deletion is a first creator; hypo analogues."""
from collections import Counter, defaultdict
NMAX = 9
cls = {}; least = {}
for n in range(1, NMAX + 1):
    L = []
    for line in open('cls_%d.txt' % n):
        head, _, tail = line.partition('|'); a = head.split()
        L.append((int(a[1]), int(a[4]), [tuple(map(int, x.split(':'))) for x in tail.split()]))
    cls[n] = L
    for h, _, _ in L: least.setdefault(h, n)
kdist = defaultdict(Counter); hypo = Counter(); orphans = []
for n in range(4, NMAX + 1):
    for h, sc, dels in cls[n]:
        k = sum(1 for _, hh in dels if least[hh] == n - 1)
        if least[h] == n:
            kdist[n][k] += 1
            if k == 0: orphans.append((n, h, sc))
        elif k == n:
            hypo[n] += 1
for n in range(4, NMAX + 1):
    tot = sum(kdist[n].values())
    print('order %d: first creators %d; k(T) (first-creator deletions) distribution: %s; all-n (every deletion a first creator): %d'
          % (n, tot, dict(sorted(kdist[n].items())), kdist[n][n]))
print('orphan first creators (k = 0):', orphans)
print('hypo analogue (NOT a first creator, every deletion a first creator):', dict(hypo))
