"""S17: minimum in/out-degree (min score distance) against first-creator status, orders 3-10."""
import sys, struct
def decode_scores(code, n):
    out = [0] * n; b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> b & 1: out[i] += 1
            else: out[j] += 1
            b += 1
    return out
def run(n, records, Sprev):
    tab = {}
    for code, H in records:
        sc = decode_scores(code, n)
        md = min(min(s, n - 1 - s) for s in sc)
        fc = H not in Sprev
        t = tab.setdefault(md, [0, 0, None, None])
        t[0] += 1; t[1] += fc
        if t[2] is None or H < t[2]: t[2] = H
        if not fc and (t[3] is None or H > t[3]): t[3] = H
    print('order %d: min degree -> classes, first creators, least H, largest non-FC H' % n)
    for md in sorted(tab):
        print('   %d: %d %d %s %s' % (md, *tab[md]))
def readset(k): return set(int(x) for x in open('S%d.txt' % k).read().split())
for n in range(4, 10):
    recs = []
    for line in open('../s16/cls_%d.txt' % n):
        p = line.split()
        recs.append((int(p[0]), int(p[1])))
    run(n, recs, readset(n - 1))
