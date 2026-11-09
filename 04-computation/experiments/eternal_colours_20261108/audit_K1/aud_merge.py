#!/usr/bin/env python3
# audit K1: merge per-range outputs of aud_sc.c / aud_sc12.c ('h count strongcount') by summing counts.
import sys, collections
cnt = collections.Counter(); sc = collections.Counter()
for fn in sys.argv[1:]:
    for line in open(fn):
        h, c, s = map(int, line.split()); cnt[h] += c; sc[h] += s
for h in sorted(cnt): print(h, cnt[h], sc[h])
