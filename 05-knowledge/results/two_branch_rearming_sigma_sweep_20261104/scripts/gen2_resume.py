#!/usr/bin/env python3
"""gen2_resume.py (klein-2026-11-04) -- one-order restart of gen2_driver.py: parents = data/gen2_order{N-1}.txt + data/gen1_order{N-1}.txt (the same parent rule as the driver), used for the order-14 run (about 80 min; output data/gen2_order14.txt).  Usage: python3 gen2_resume.py 14 (run in the directory holding the gen2 binary and the data files)."""
import subprocess, sys
LABELG = '/opt/homebrew/bin/labelg'
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
N = int(sys.argv[1])  # target order; parents = gen2_order{N-1} + gen1_order{N-1}
cur = [l.split()[1] for l in open('gen2_order%d.txt' % (N-1))] + [l.split()[1] for l in open('gen1_order%d.txt' % (N-1))]
print('order %d: %d ready parents' % (N-1, len(cur)), flush=True)
targets = [4289, 4291, 12869, 12871, 12875, 12877, 38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635]
inp = '\n'.join(' '.join(map(str, rows_from_d6(c))) for c in cur) + '\n'
out = subprocess.run(['./gen2', str(N - 1)], input=inp, capture_output=True, text=True).stdout.split('\n')
recs = [l.split() for l in out if l.strip()]
can = subprocess.run([LABELG, '-q'], input='\n'.join(r[1] for r in recs) + '\n', capture_output=True, text=True).stdout.split()
cls = {}
for r, c in zip(recs, can): cls[c] = int(r[0])
vals = sorted(set(cls.values()))
print('order %d: %d digit-2 classes; R2(%d): %d values in [%d, %d]; tree values present: %s' % (N, len(cls), N, len(vals), vals[0], vals[-1], [t for t in targets if t in set(vals)]), flush=True)
with open('gen2_order%d.txt' % N, 'w') as fh:
    for c, h in cls.items(): fh.write('%d %s\n' % (h, c))
