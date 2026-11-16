#!/usr/bin/env python3
"""merge_slices.py out.txt slice1 slice2 ...: merge partial two-apex censuses (two_apex_slice.c) into one census file in the
format of two_apex_small.c: header, then 'V h count minpop orient' lines (count summed, minpop minimised, orient OR-ed)."""
import sys, re
out, slices = sys.argv[1], sys.argv[2:]
vals = {}; tapes = 0; maxh = 0; hmax = None; m = None; cmax = None
for fn in slices:
    for line in open(fn):
        if line.startswith('#'):
            mm = re.search(r"popcount\(b\) <= (\d+), spine (\d+): (\d+) tapes, max H = (\d+), distinct values below (\d+)", line)
            cmax, m, t, mh, hm = map(int, mm.groups()); tapes += t; maxh = max(maxh, mh); hmax = hm
        elif line.startswith('V '):
            _, h, c, mp, o = line.split(); h, c, mp, o = int(h), int(c), int(mp), int(o)
            if h in vals: vals[h] = (vals[h][0] + c, min(vals[h][1], mp), vals[h][2] | o)
            else: vals[h] = (c, mp, o)
with open(out, 'w') as f:
    f.write(f"# two-apex tapes with popcount(b) <= {cmax}, spine {m}: {tapes} tapes, max H = {maxh}, distinct values below {hmax}: {len(vals)}\n")
    for h in sorted(vals): f.write(f"V {h} {vals[h][0]} {vals[h][1]} {vals[h][2]}\n")
print(f"merged {len(slices)} slices: {tapes} tapes, max H {maxh}, {len(vals)} distinct values below {hmax}")
