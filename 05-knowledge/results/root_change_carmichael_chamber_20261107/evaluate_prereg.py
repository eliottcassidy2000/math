"""Evaluate the S14 pre-registered predictions P1-P4 (PREREGISTRATION.md) on the output of `rc_census 1000000000 2 512`.

usage: python3 evaluate_prereg.py rc_1e9.txt  > PREREG_OUTCOMES.txt
Every check raises SystemExit on malformed input; outcomes are printed, never asserted.
"""
import re
import sys
from pathlib import Path

KNOWN_ABOVE_1E5 = [100087, 111373, 116215, 125497, 136607, 138901, 205903, 229417, 485833, 532207, 586303,
                   797083, 1120831, 1847899, 3780967]
WINDOWS = [(67108864, 134217728), (134217728, 268435456), (268435456, 536870912), (536870912, 1000000000)]

text = Path(sys.argv[1]).read_bytes().decode('utf-8').replace('\r\n', '\n')
lines = text.split('\n')

m = re.search(r'^<=2 root changes: odd misses (\d+), largest (\d+)$', text, re.M)
if not m:
    raise SystemExit('no <=2 summary line')
miss2, largest2 = int(m.group(1)), int(m.group(2))
m = re.search(r'^misses \(<=2\) above 10\^5:(.*)$', text, re.M)
if not m:
    raise SystemExit('no misses line')
above = [int(x) for x in m.group(1).split()]
p1 = (miss2 == 140 and largest2 == 3780967 and above == KNOWN_ABOVE_1E5)
print('P1 coverage: <=2 misses %d, largest %d, misses above 1e5 = %s' % (miss2, largest2, above))
print('P1 %s (prediction: exactly the 140 misses below 1e8, largest 3780967, none in (1e8, 1e9])'
      % ('CONFIRMED' if p1 else 'FAILED'))

win = {}
pat = re.compile(r'^  \[(\d+),(\d+)\): n=(\d+) mean=([\d.]+) min=(\d+) min_mod3=\((\d+),(\d+),(\d+)\) sat255=(\d+) '
                 r'c0..8=([\d,]+)$')
in_mult = False
for ln in lines:
    if ln.startswith('multiplicity'):
        in_mult = True
        continue
    if in_mult:
        mm = pat.match(ln)
        if mm:
            lo, hi = int(mm.group(1)), int(mm.group(2))
            win[(lo, hi)] = dict(n=int(mm.group(3)), mean=float(mm.group(4)), mn=int(mm.group(5)),
                                 m3=tuple(int(mm.group(i)) for i in (6, 7, 8)), sat=int(mm.group(9)),
                                 c=[int(x) for x in mm.group(10).split(',')])
for w in WINDOWS:
    if w not in win:
        raise SystemExit('missing window %s' % (w,))

p2 = all(win[w]['mn'] >= 10 for w in WINDOWS)
print('P2 window minima: %s' % [(w, win[w]['mn']) for w in WINDOWS])
print('P2 %s (prediction: min M_2 >= 10 in each of the four windows)' % ('CONFIRMED' if p2 else 'FAILED'))

singles = [(lo, hi, v['c'][1]) for (lo, hi), v in sorted(win.items()) if v['c'][1] and lo > 1]
m = re.search(r'^odd N > 10\^5 with multiplicity <= 2:(.*)$', text, re.M)
low = [tuple(int(y) for y in x.split(':')) for x in m.group(1).split()] if m else None
if low is None:
    raise SystemExit('no low-multiplicity line')
ones = [n for n, c in low if c == 1]
p3 = (not singles) and (not ones) and win[(1, 2)]['c'][1] == 1
print('P3 windows with singletons above N = 1: %s; N > 1e5 with M_2 = 1: %s' % (singles, ones))
print('P3 %s (prediction: no odd 1 < N <= 1e9 with M_2(N) = 1)' % ('CONFIRMED' if p3 else 'FAILED'))

p4 = all(win[w]['m3'][1] < win[w]['m3'][0] and win[w]['m3'][1] < win[w]['m3'][2] for w in WINDOWS)
print('P4 window minima by residue mod 3 (0,1,2): %s' % [(w, win[w]['m3']) for w in WINDOWS])
print('P4 %s (prediction: the N = 1 (mod 3) minimum is strictly smallest in each window)'
      % ('CONFIRMED' if p4 else 'FAILED'))

print('low-multiplicity N > 1e5 (M_2 <= 2) above 1e8: %s' % [(n, c) for n, c in low if n > 10**8])
print('all windows (lo, n, mean, min, min_mod3, sat255, c0..8):')
for (lo, hi), v in sorted(win.items()):
    print('  [%d,%d) n=%d mean=%.2f min=%d m3=%s sat=%d c=%s' % (lo, hi, v['n'], v['mean'], v['mn'], v['m3'], v['sat'],
                                                               v['c']))
