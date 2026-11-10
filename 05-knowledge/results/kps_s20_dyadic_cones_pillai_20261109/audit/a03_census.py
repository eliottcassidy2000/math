# a03: own census over the certified class lists of order <= 8 (computed directly with a_paths.exe,
# NOT from the S19 port-statistics file):
#  - strong classes with Q a power of two (claim: 11 classes, H in {5,19,217,223,599}, Q = 2^5,2^8,2^12,2^12,2^13);
#  - dyadic census (claims: 148 attained odd values <= 300; 69 with a dyadic port; 181 distinct dyadic triples;
#    odd parts of responses at H = 5: {1,3,5,7,9}; at H = 13: {7,9,13,15});
#  - writes all rooted triples to a03_triples_le8.txt for later use (n code v h e R).
import sys, collections
from a_lib import *
S19 = sys.argv[1]
rows = []
qpow = []
allstrong = collections.Counter()
for n in range(1, 9):
    codes = [int(l.split()[0]) for l in open('%s/cls_%d.txt' % (S19, n)) if l.strip()]
    tours = [from_code(c, n) for c in codes]
    res = run(tours)
    for c, A, r in zip(codes, tours, res):
        st = is_strong(A)
        if st: allstrong[n] += 1
        if st and ispow2(r['Q']):
            qpow.append((n, c, r['H'], r['Q']))
        for v in range(n):
            rows.append((n, c, v, r['H'], r['e'][v], r['R'][v]))
print('strong classes per order (is_strong, order 1 counted as strong):', dict(allstrong))
print('strong classes with Q a power of two:')
for n, c, H, Q in qpow:
    print('  n=%d code=%d H=%d Q=2^%d' % (n, c, H, Q.bit_length() - 1))
print('  total incl. order 1: %d; of order >= 3: %d' % (len(qpow), sum(1 for x in qpow if x[0] >= 3)))
with open('a03_triples_le8.txt', 'w') as f:
    for row in rows: f.write(' '.join(map(str, row)) + '\n')
vals = set(); dy = {}; triples = set(); odd_at = collections.defaultdict(set)
eviol = 0; Rlt = 0
for n, c, v, h, e, R in rows:
    vals.add(h); odd_at[h].add(oddpart(R))
    if (R - e) % 2: eviol += 1
    if R < h: Rlt += 1
    if ispow2(R):
        triples.add((h, e, R)); j = R.bit_length() - 1
        if h not in dy or j < dy[h]: dy[h] = j
small = sorted(x for x in vals if x <= 300)
print('rooted triples:', len(rows), ' parity R = e mod 2 violations:', eviol, ' R < h violations:', Rlt)
print('attained odd values <= 300 at order <= 8:', len(small), ' missing odd values <= 300:', [x for x in range(1, 301, 2) if x not in vals])
print('values <= 300 with a dyadic port:', sum(1 for x in small if x in dy))
print('distinct dyadic triples (h,e,R):', len(triples))
print('response odd parts at H = 5:', sorted(odd_at[5]), '  at H = 13:', sorted(odd_at[13]))
print('dyadic triples with h = 5:', sorted(t for t in triples if t[0] == 5), ' with h = 13:', sorted(t for t in triples if t[0] == 13))
print('least j per value (<= 300):', ' '.join('%d:%d' % (h, dy[h]) for h in sorted(dy) if h <= 300))
