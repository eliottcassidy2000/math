"""S17: value-level genealogy ledger for orders 4-10.  For every value h first created at order n <= 10:
the number of first-creator classes, and the H values of their first-creator deletions (parent values, new at n - 1).
Orders <= 9 from S16's cls_n.txt (deletion codes and H), order 10 from classes10.bin (34-byte records)."""
import struct
S = {k: set(int(x) for x in open('S%d.txt' % k).read().split()) for k in range(1, 10)}
S[10] = set(int(x) for x in open('S10_klein.txt').read().split())
out = open('lineage_le10.txt', 'w')
out.write('# value-level genealogy: h  n(h)  #first-creator classes  #with a first-creator parent  parent values (H of first-creator deletions, new at n(h)-1)\n')
for n in range(4, 11):
    fcount = {}; withp = {}; parents = {}
    def handle(h, dels):
        if h in S[n - 1]: return
        fcount[h] = fcount.get(h, 0) + 1
        pv = {d for d in dels if d not in S[n - 2]}
        if pv: withp[h] = withp.get(h, 0) + 1
        parents.setdefault(h, set()).update(pv)
    if n <= 9:
        for line in open('../s16/cls_%d.txt' % n):
            left, right = line.split('|')
            h = int(left.split()[1])
            dels = [int(d.split(':')[1]) for d in right.split()]
            handle(h, dels)
    else:
        with open('classes10.bin', 'rb') as fh:
            while True:
                rec = fh.read(34)
                if len(rec) < 34: break
                h = struct.unpack_from('<I', rec, 8)[0]
                if h in S[9]: continue
                dels = struct.unpack_from('<10H', rec, 12)
                handle(h, dels)
    for h in sorted(fcount):
        pv = sorted(parents[h])
        ps = ','.join(map(str, pv)) if len(pv) <= 12 else '%s,...,%s (%d values)' % (','.join(map(str, pv[:6])), ','.join(map(str, pv[-3:])), len(pv))
        out.write('%d %d %d %d %s\n' % (h, n, fcount[h], withp.get(h, 0), ps if pv else 'ORPHAN'))
    print('order', n, 'new values', len(fcount), 'orphan values', [h for h in fcount if not parents[h]])
out.close()
