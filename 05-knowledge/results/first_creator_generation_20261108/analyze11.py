"""S17: evaluate the order 10 -> 11 extension run (R1-R5).  usage: python3 analyze11.py"""
import glob, re
from array import array
HCAP = 100000
A = [0] * HCAP; F = [0] * HCAP
for p in sorted(glob.glob('o11_*.cnt')):
    a = array('Q')
    with open(p, 'rb') as fh:
        a.fromfile(fh, 2 * HCAP)
    for i in range(HCAP):
        if a[i]: A[i] += a[i]
        if a[HCAP + i]: F[i] += a[HCAP + i]
S10 = set(int(x) for x in open('S10_klein.txt').read().split())
S11 = set(int(x) for x in open('S11_klein.txt').read().split())
got = {h for h in range(HCAP) if A[h]}
print('S(11) from extensions: %d values, max %d; equals klein S(11): %s (missing %d, extra %d)'
      % (len(got), max(got), got == S11, len(S11 - got), len(got - S11)))
new = sorted(got - S10)
unc = [h for h in new if F[h] == 0]
only = [h for h in new if F[h] == A[h]]
print('new values at order 11: %d' % len(new))
print('R1: new values not reached from an order-10 first creator: %d %s' % (len(unc), unc[:50]))
print('R5: new values reached only from first-creator parents: %d = %.4f' % (len(only), len(only) / len(new)))
tot = {'classes': 0, 'fc': 0, 'sterile': 0, 'hypo': 0, 'chk': 0, 'fail': 0, 'low': 0}
hist = [0] * 12; minf = None
for p in sorted(glob.glob('o11_*.log')):
    s = open(p).read()
    m = re.search(r'classes (\d+), first creators (\d+), sterile first creators (\d+), hypo pairs (\d+), checked classes (\d+), check failures (\d+)', s)
    for k, v in zip(['classes', 'fc', 'sterile', 'hypo', 'chk', 'fail'], m.groups()): tot[k] += int(v)
    hl = [l for l in s.splitlines() if l.startswith('fertile-pattern histogram')][0]
    hist = [x + int(y) for x, y in zip(hist, hl.rsplit('):', 1)[1].split())]
    lm = re.search(r'fewer than half of the patterns fertile: (\d+); minimum fertile (\d+) \(code (\d+), H (\d+)\)', s)
    tot['low'] += int(lm.group(1))
    cand = (int(lm.group(2)), int(lm.group(3)), int(lm.group(4)))
    if minf is None or cand < minf: minf = cand
    for line in s.splitlines():
        if line.startswith('HYPO') or line.startswith('LOWFERT'): print('  ', line)
print('classes %(classes)d, first creators %(fc)d; brute-force checked classes %(chk)d, failures %(fail)d' % tot)
print('R2: hypo pairs %(hypo)d' % tot)
print('R3: sterile first creators %(sterile)d' % tot)
print('R4: first creators with fewer than 512 fertile patterns: %(low)d; minimum fertile %%s' % tot % (minf,))
print('fertile histogram (bucket b: [2^(b-1), 2^b)):', hist)
# multiplicity statistics for new values
import statistics
mfc = sorted(F[h] for h in new)
print('first-creator (class, pattern) realizations per new value: min %d, median %d, max %d' % (mfc[0], mfc[len(mfc) // 2], mfc[-1]))
low = sorted(new, key=lambda h: F[h])[:15]
print('new values with the fewest first-creator realizations:', [(h, F[h], A[h]) for h in low])
