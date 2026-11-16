#!/usr/bin/env python3
"""merge_hratio4.py: merge the outputs of the eight restricted runs `hratio4 mats3.txt frag.txt L b0` (b0 = 0..7).
usage: merge_hratio4.py out_0 ... out_7.  Per length: the exact global maxima of H(Wl)/H(W) and Z_l/H, the total number
of pairs (W, l) attaining each (summed over the runs whose maximum equals the global one), the first maximiser, and the
word count check (sum of the runs' counts = 8^n for n >= 1; length 0 is the single empty word)."""
import sys, re
from fractions import Fraction
pat = re.compile(r"length (\d+): words (\d+); max H\(Wl\)/H\(W\) = (\d+)/(\d+) = [\d.]+, attained by (\d+) pairs \(W,l\), first: letter (\d+), word([ \d]*); max Z_l/H = (\d+)/(\d+) = [\d.]+, attained by (\d+) pairs, first: letter (\d+), word([ \d]*)")
runs = []
for fn in sys.argv[1:]:
    d = {}
    for line in open(fn):
        m = pat.match(line.strip())
        if m:
            n = int(m.group(1))
            d[n] = dict(words=int(m.group(2)), H=Fraction(int(m.group(3)), int(m.group(4))), Ht=int(m.group(5)), Hl=int(m.group(6)), Hw=m.group(7).split(),
                        Z=Fraction(int(m.group(8)), int(m.group(9))), Zt=int(m.group(10)), Zl=int(m.group(11)), Zw=m.group(12).split())
    runs.append(d)
lengths = sorted(set.intersection(*[set(r) for r in runs]))
print(f"merged {len(runs)} restricted runs, lengths {lengths[0]}..{lengths[-1]}")
for n in lengths:
    words = 1 if n == 0 else sum(r[n]['words'] for r in runs)
    assert words == 8 ** n, (n, words)
    for key, tk, lk, wk, name in (('H', 'Ht', 'Hl', 'Hw', 'H(Wl)/H(W)'), ('Z', 'Zt', 'Zl', 'Zw', 'Z_l/H')):
        best = max(r[n][key] for r in runs)
        ties = sum(r[n][tk] for r in runs if r[n][key] == best)
        first = next(r[n] for r in runs if r[n][key] == best)
        print(f"length {n}: words {words}; max {name} = {best.numerator}/{best.denominator} = {float(best):.6f}, attained by {ties} pairs (W,l), first: letter {first[lk]}, word {' '.join(first[wk])}")
