"""Share of first-creator isomorphism classes at orders 3..11 from klein-2026-11-06's exact censuses.
A class of order n is a first creator iff H(T) is not in S(n-1).  census_n.txt lines: 'T h c3? strong? count'.
We use the LAST field as the class count of value h (checked against the per-order class totals)."""
import gzip, os
from spectra import in_S, KLEIN
tot_expected = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536, 10: 9733056, 11: 903753248}
for n in range(3, 12):
    fn = os.path.join(KLEIN, 'census_%d.txt' % n)
    op = open
    if not os.path.exists(fn):
        fn += '.gz'; op = gzip.open
    per = {}
    with op(fn, 'rt') as f:
        for line in f:
            if line.startswith('T '):
                p = line.split()
                per[int(p[1])] = per.get(int(p[1]), 0) + int(p[-1])
    tot = sum(per.values())
    fc = sum(c for h, c in per.items() if not in_S(n - 1, h))
    newv = sum(1 for h in per if not in_S(n - 1, h))
    print('n=%2d classes=%10d (expected %10d) first-creator classes=%10d share=%.4f  new values=%d' % (
        n, tot, tot_expected[n], fc, fc / tot, newv))
