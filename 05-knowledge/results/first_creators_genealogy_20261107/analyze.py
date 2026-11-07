"""First-creator fibers of odd Hamiltonian-path counts (kind-pasteur S16), from classgen output cls_<n>.txt, n <= NMAX."""
import sys
from collections import defaultdict, Counter
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9
cls = {}          # n -> list of dicts
S = {}            # n -> set of H values
for n in range(1, NMAX + 1):
    L = []
    for line in open('cls_%d.txt' % n):
        head, _, tail = line.partition('|')
        a = head.split()
        dels = [tuple(map(int, x.split(':'))) for x in tail.split()]
        L.append(dict(code=int(a[0]), H=int(a[1]), c3=int(a[2]), strong=int(a[3]), SC=int(a[4]), aut=int(a[5]),
                      cc=int(a[6]), dels=dels))
    cls[n] = L; S[n] = set(d['H'] for d in L)
least = {}
for n in range(1, NMAX + 1):
    for h in S[n]:
        least.setdefault(h, n)
fib = defaultdict(list)
for n in range(1, NMAX + 1):
    for d in cls[n]:
        if least[d['H']] == n:
            fib[d['H']].append(d)
# (a) fiber structure
rows = []
for h in sorted(fib):
    F = fib[h]; s = sum(d['SC'] for d in F); ns = len(F) - s
    assert ns % 2 == 0
    codes = set(d['code'] for d in F)
    assert all(d['cc'] in codes for d in F)          # closed under converse
    rows.append((h, least[h], len(F), s, ns // 2))
kinds = Counter()
for h, n, size, s, p in rows:
    kinds[('SC-only' if p == 0 else 'NS-only' if s == 0 else 'mixed', 'unique' if (s, p) in [(1, 0), (0, 1)] else 'multi')] += 1
print('values with least order <= %d: %d' % (NMAX, len(rows)))
print('kinds:', dict(kinds))
print('(s,p) = (1,0) or (0,1) [the conjectured dichotomy]:', sum(1 for r in rows if (r[3], r[4]) in [(1, 0), (0, 1)]))
print('first 40 rows (h, n(h), fiber size, #SC, #NS pairs):')
print(rows[:40])
par = Counter((s % 2) for _, _, _, s, _ in rows)
print('parity of the number of SC first creators: even', par[0], 'odd', par[1])
# (b) genealogy: does a first creator have a vertex deletion that is a first creator of a smaller value?
has_fc_parent = Counter(); value_has = Counter()
for h, F in fib.items():
    n = least[h]
    if n < 4: continue
    anyv = False
    for d in F:
        ok = any(least.get(hh, 99) == n - 1 for _, hh in d['dels'])
        has_fc_parent[ok] += 1; anyv |= ok
    value_has[anyv] += 1
print('first creators with a first-creator parent (some T - v new at order n-1):', dict(has_fc_parent))
print('values with at least one such first creator:', dict(value_has))
missing = sorted(h for h, F in fib.items() if least[h] >= 4 and not any(any(least.get(hh, 99) == least[h] - 1 for _, hh in d['dels']) for d in F))
print('values with no first-creator lineage (first 40):', missing[:40], 'count', len(missing))
