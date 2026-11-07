from collections import Counter, defaultdict
NMAX = 9
least = {}; fib = defaultdict(list); maxH = {}
for n in range(1, NMAX + 1):
    rows = []
    for line in open('cls_%d.txt' % n):
        a = line.split('|')[0].split(); rows.append((int(a[1]), int(a[4]), int(a[5])))
    maxH[n] = max(r[0] for r in rows)
    for h, sc, aut in rows: least.setdefault(h, n)
    for h, sc, aut in rows:
        if least[h] == n: fib[h].append((sc, aut))
uniqSC = []; uniqNS = []; kind = {}
for h, F in fib.items():
    s = sum(x for x, _ in F); p = (len(F) - s) // 2
    kind[h] = 'SConly' if p == 0 else ('NSonly' if s == 0 else 'mixed')
    if (s, p) == (1, 0): uniqSC.append(h)
    if (s, p) == (0, 1): uniqNS.append(h)
print('dichotomy values, unique SC class  (%d):' % len(uniqSC), sorted(uniqSC))
print('dichotomy values, one converse pair (%d):' % len(uniqNS), sorted(uniqNS))
for m in (4, 8, 3):
    tab = defaultdict(Counter)
    for h, k in kind.items():
        if h > 1: tab[k][h % m] += 1
    print('h mod %d by kind:' % m, {k: dict(sorted(v.items())) for k, v in tab.items()})
# position within the order's new range
pos = defaultdict(list)
for h, k in kind.items():
    n = least[h]
    if n >= 6: pos[k].append(h / maxH[n])
for k, v in pos.items():
    v.sort(); print(k, 'count', len(v), 'median h/maxH(n(h)) %.3f' % v[len(v)//2], 'top-decile share %.2f' % (sum(1 for x in v if x > 0.9) / len(v)))
