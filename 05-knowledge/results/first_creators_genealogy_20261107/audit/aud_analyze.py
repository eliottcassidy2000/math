"""Auditor's analysis for THM-4922 / HYP-9512 (independent of analyze.py / hypo.py / split.py).

Reads the auditor's own class files  <dir>/aud_cls_<n>.txt  (n = 1..9; line: code H c3 strong SC aut hc | H(T-v)...)
and the auditor's order-10 per-H table  <dir>/aud_o10_byH.txt, and optionally the packet's order10.out and klein's
census files for comparison.  Usage:
    python3 aud_analyze.py RUNDIR PACKETDIR KLEINDATADIR
"""
import sys, os, gzip
from collections import defaultdict, Counter

RUN, PACK, KLEIN = sys.argv[1], sys.argv[2], sys.argv[3]
NMAX = 9

cls = {}
for n in range(1, NMAX + 1):
    rows = []
    for line in open(os.path.join(RUN, 'aud_cls_%d.txt' % n)):
        head, _, tail = line.partition('|')
        a = head.split()
        rows.append(dict(code=int(a[0]), H=int(a[1]), c3=int(a[2]), strong=int(a[3]), SC=int(a[4]), aut=int(a[5]),
                         hc=int(a[6]), dels=[int(x) for x in tail.split()]))
    cls[n] = rows

# order 10 per-H table
o10 = {}
for line in open(os.path.join(RUN, 'aud_o10_byH.txt')):
    if line.startswith('#'):
        continue
    a = list(map(int, line.split()))
    o10[a[0]] = dict(all=a[1], allsc=a[2], fib=a[3], s=a[4], p=a[5], full=a[6], nonstrong=a[7], kmin=a[8],
                     minc3=a[9], cntc3=a[10], minhc=a[11], cnthc=a[12])

S = {n: set(r['H'] for r in cls[n]) for n in cls}
S[10] = set(o10)
least = {}
for n in range(1, 11):
    for h in S[n]:
        least.setdefault(h, n)
maxH = {n: max(S[n]) for n in S}

out = []
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); out.append(s)

P('== class counts (A000568) and SC counts (A002785)')
P('classes', [len(cls[n]) for n in range(1, 10)] + [sum(v['all'] for v in o10.values())])
P('SC classes', [sum(r['SC'] for r in cls[n]) for n in range(1, 10)] + [sum(v['allsc'] for v in o10.values())])
P('|S(n)| =', [len(set().union(*[S[m] for m in range(1, n + 1)])) for n in range(1, 11)])
P('|N(n)| =', [sum(1 for h in least if least[h] == n) for n in range(1, 11)])
P('max H(n) =', [maxH[n] for n in range(1, 11)])

# ---------------- fibers
fib = defaultdict(list)
for n in range(1, NMAX + 1):
    for r in cls[n]:
        if least[r['H']] == n:
            fib[r['H']].append(r)
rows = {}
for h, F in fib.items():
    s = sum(r['SC'] for r in F); ns = len(F) - s
    assert ns % 2 == 0
    rows[h] = (least[h], len(F), s, ns // 2)
for h, v in o10.items():
    if v['fib']:
        rows[h] = (10, v['fib'], v['s'], v['p'])
P('== fibers')
le9 = {h: v for h, v in rows.items() if v[0] <= 9}
P('values with n(h) <= 9:', len(le9), '(including h = 1:', 1 in le9, ')')
d10 = sorted(h for h, v in le9.items() if (v[2], v[3]) == (1, 0))
d01 = sorted(h for h, v in le9.items() if (v[2], v[3]) == (0, 1))
P('(s,p)=(1,0):', len(d10), d10)
P('(s,p)=(0,1):', len(d01), d01)
P('dichotomy values with n(h) <= 8:', sum(1 for h, v in le9.items() if v[0] <= 8 and v[2] + v[3] == 1),
  'of', sum(1 for h, v in le9.items() if v[0] <= 8))
fails = sorted(h for h, v in le9.items() if v[2] + v[3] != 1)
P('smallest failures (h, n(h), fiber, s, p):', [(h,) + rows[h] for h in fails[:12]])
for h in (9, 15, 17, 23, 29, 33, 37):
    P('  h=%d: n=%d fiber=%d s=%d p=%d' % ((h,) + rows[h]))
P('max fiber n(h)<=9:', max((v[1], h) for h, v in le9.items()), ' max fiber at n(h)=10:',
  max((v[1], h) for h, v in rows.items() if v[0] == 10))
for n in range(3, 11):
    P('  max fiber at least order %d: %s' % (n, max(((v[1], h) for h, v in rows.items() if v[0] == n), default=None)))
kinds = Counter()
for h, (n, size, s, p) in le9.items():
    kinds['SC-only' if p == 0 else 'NS-only' if s == 0 else 'mixed'] += 1
P('kinds n(h) <= 9:', dict(kinds))
d10_10 = sorted(h for h, v in rows.items() if v[0] == 10 and (v[2], v[3]) in [(1, 0), (0, 1)])
P('order 10 dichotomy values:', len(d10_10), 'range', (min(d10_10), max(d10_10)))
P('  order-10 dichotomy (1,0):', [h for h in d10_10 if rows[h][2] == 1])
P('  max H(10) fiber:', rows[maxH[10]])

# ---------------- blue / black split
P('== blue/black split (blue: some SC first creator)')
tab = Counter()
for h, (n, size, s, p) in rows.items():
    if h == 1:
        continue
    tab[(n <= 9, 'even' if n % 2 == 0 else 'odd', h % 4, 'blue' if s > 0 else 'black')] += 1
for le in (True, False):
    P('  n(h) <= 9' if le else '  n(h) = 10')
    for par in ('even', 'odd'):
        for m in (1, 3):
            P('    %s %d: blue %d black %d' % (par, m, tab[(le, par, m, 'blue')], tab[(le, par, m, 'black')]))
blue9 = sum(1 for h, v in le9.items() if v[2] > 0)
P('  blue/black with n(h)<=9 incl. h=1:', blue9, len(le9) - blue9)
# per-order blue share among values not forced black
for n in range(4, 11):
    vals = [h for h, v in rows.items() if v[0] == n and not (n % 2 == 0 and h % 4 == 3)]
    b = sum(1 for h in vals if rows[h][2] > 0)
    P('  order %d: values not forced black %d, blue %d (%.3f)' % (n, len(vals), b, b / len(vals) if vals else 0))
ev = [h for h, v in le9.items() if v[0] % 2 == 0 and h % 4 == 1 and h > 1]
P('  even n(h)<=8, h=1 mod 4: blue share %.3f' % (sum(1 for h in ev if rows[h][2] > 0) / len(ev)))
od = [h for h, v in le9.items() if v[0] % 2 == 1 and h > 1]
P('  odd n(h)<=9: blue share %.3f' % (sum(1 for h in od if rows[h][2] > 0) / len(od)))

# ---------------- SC mod 4
P('== SC classes by H mod 4')
for n in range(1, 10):
    c = Counter(r['H'] % 4 for r in cls[n] if r['SC'])
    P('  n=%d: 1 mod 4: %d, 3 mod 4: %d' % (n, c[1], c[3]))

# ---------------- position of kinds
pos = defaultdict(list)
for h, (n, size, s, p) in le9.items():
    if n >= 6:
        k = 'SC-only' if p == 0 else 'NS-only' if s == 0 else 'mixed'
        pos[k].append(h / maxH[n])
for k, v in sorted(pos.items()):
    v.sort(); P('  %s: count %d, median h/maxH(n(h)) %.3f' % (k, len(v), v[len(v) // 2]))

# ---------------- k(T), orphans, hypo analogues
P('== k(T) = #{v : T - v is a first creator}')
for n in range(2, 10):
    kd = Counter(); hypo = 0; orph = []; nonstrong = 0
    for r in cls[n]:
        k = sum(1 for g in r['dels'] if least[g] == n - 1)
        if least[r['H']] == n:
            kd[k] += 1
            if k == 0: orph.append((r['H'], r['SC'], r['aut']))
            if not r['strong']: nonstrong += 1
        elif k == n:
            hypo += 1
    tot = sum(kd.values())
    P('  order %d: first creators %d, k-dist %s, k=n %d (%.3f), orphans %s, hypo analogues %d, non-strong first creators %d'
      % (n, tot, dict(sorted(kd.items())), kd[n], kd[n] / tot if tot else 0, orph, hypo, nonstrong))
P('  first creators at orders 4-9:', sum(1 for n in range(4, 10) for r in cls[n] if least[r['H']] == n))
o63 = [r for r in cls[8] if r['H'] == 63]
P('  order-8 classes with H=63:', [(r['SC'], r['aut'], sorted(r['dels']), r['strong'], r['c3']) for r in o63])
P('  least orders of their deletion values:', sorted(set(least[g] for r in o63 for g in r['dels'])))

# ---------------- thresholds of the proved genealogy
P('== genealogy threshold n*maxH(n-2)')
for n in range(4, 11):
    thr = n * maxH[n - 2]
    new = [h for h, v in rows.items() if v[0] == n]
    P('  n=%d: threshold %d; new values %d, above threshold %d' % (n, thr, len(new), sum(1 for h in new if h > thr)))

# ---------------- least-level refinement (c3) and Hamiltonian-cycle refinement
P('== refinement by least c3 level / least hc at the least order')
def refine_counts():
    resc3 = {}; reshc = {}
    for n in range(1, 10):
        byh = defaultdict(list)
        for r in cls[n]:
            if least[r['H']] == n:
                byh[r['H']].append(r)
        for h, F in byh.items():
            m = min(r['c3'] for r in F); resc3[h] = sum(1 for r in F if r['c3'] == m)
            m = min(r['hc'] for r in F); reshc[h] = sum(1 for r in F if r['hc'] == m)
    for h, v in o10.items():
        if v['fib']:
            resc3[h] = v['cntc3']; reshc[h] = v['cnthc']
    return resc3, reshc
rc3, rhc = refine_counts()
P('  values n(h) <= 10:', len(rc3), '; >2 classes at least order & least c3:', sum(1 for v in rc3.values() if v > 2),
  '; >2 classes at least order & least hc:', sum(1 for v in rhc.values() if v > 2))
# NB for order 10 the per-H minima are over ALL order-10 classes with that H, which equal the first creators
# because these H are new at order 10.

# ---------------- comparison with the packet's order10.out
P('== comparison with order10.out')
pv = {}
for line in open(os.path.join(PACK, 'order10.out')):
    if line.startswith('V '):
        a = list(map(int, line.split()[1:]))
        pv[a[0]] = tuple(a[1:])
mine = {h: (v['fib'], v['s'], v['p'], v['full']) for h, v in o10.items() if v['fib']}
P('  packet V lines:', len(pv), ' auditor values:', len(mine), ' identical tables:', pv == mine)
diff = [h for h in set(pv) | set(mine) if pv.get(h) != mine.get(h)]
P('  differing h:', sorted(diff)[:20])
m4 = Counter((h % 4, 'blue' if mine[h][1] > 0 else 'black') for h in mine)
P('  order-10 new values by h mod 4 and kind:', dict(m4))

# ---------------- comparison with klein's census (h, c3, strong) class counts, orders 3..10
P("== comparison with klein's census_<n>.txt (T h c3 strong count)")
def klein_T(n):
    fn = os.path.join(KLEIN, 'census_%d.txt' % n)
    op = open
    if not os.path.exists(fn):
        fn += '.gz'; op = lambda f: gzip.open(f, 'rt')
    d = Counter(); hcC = Counter()
    for line in op(fn):
        a = line.split()
        if a and a[0] == 'T':
            d[(int(a[1]), int(a[2]), int(a[3]))] += int(a[4])
        elif a and a[0] == 'C':
            hcC[(int(a[1]), int(a[2]))] += int(a[3])
    return d, hcC
for n in range(3, 10):
    kd, khc = klein_T(n)
    md = Counter((r['H'], r['c3'], r['strong']) for r in cls[n])
    mhc = Counter((r['H'], r['hc']) for r in cls[n])
    P('  order %d: (h,c3,strong) table identical: %s; (h,hc) table identical: %s' % (n, kd == md, khc == mhc))
kd10, khc10 = klein_T(10)
ok = True; okhc = True
byh = defaultdict(dict); byhhc = defaultdict(dict)
for (h, c3, st), c in kd10.items():
    byh[h][c3] = byh[h].get(c3, 0) + c
for (h, hc), c in khc10.items():
    byhhc[h][hc] = c
for h, v in o10.items():
    m = min(byh[h]);
    if (m, byh[h][m]) != (v['minc3'], v['cntc3']) or sum(byh[h].values()) != v['all']:
        ok = False
    m = min(byhhc[h])
    if (m, byhhc[h][m]) != (v['minhc'], v['cnthc']):
        okhc = False
P('  order 10: per-h total, least c3 and count at least c3 agree with klein:', ok and set(byh) == set(o10),
  '; least hc and count agree:', okhc)
nonstrong_new10 = sum(c for (h, c3, st), c in kd10.items() if st == 0 and h not in set().union(*[S[m] for m in range(1, 10)]))
P('  klein order 10: non-strong classes with a value new at order 10:', nonstrong_new10)
kd11, _ = klein_T(11)
S10all = set().union(*[S[m] for m in range(1, 11)])
ns11 = sum(c for (h, c3, st), c in kd11.items() if st == 0 and h not in S10all)
new11 = set(h for (h, c3, st) in kd11 if h not in S10all)
P('  klein order 11: values new at order 11: %d; non-strong classes with a new value: %d' % (len(new11), ns11))

with open(os.path.join(RUN, 'aud_analyze.out'), 'w') as f:
    f.write('\n'.join(out) + '\n')
# per-value fiber table for every h with n(h) <= 10
with open(os.path.join(RUN, 'aud_fibers_le10.txt'), 'w') as f:
    f.write('# h n(h) fiber_classes SC_classes converse_pairs  (first creators; auditor census)\n')
    for h in sorted(rows):
        n, size, s, p = rows[h]
        f.write('%d %d %d %d %d\n' % (h, n, size, s, p))
