#!/usr/bin/env python3
"""merge_parts.py hist OUT part1 part2 ...   : sum "H count" histograms.
merge_parts.py bits OUT part1 part2 ...   : bitwise OR of equally long bitset files (pairs.c output)."""
import sys
mode, out, parts = sys.argv[1], sys.argv[2], sys.argv[3:]
if mode == "hist":
    d = {}
    for f in parts:
        for line in open(f):
            h, c = map(int, line.split()); d[h] = d.get(h, 0) + c
    open(out, "w").write("".join(f"{h} {c}\n" for h, c in sorted(d.items())))
else:
    acc = None
    for f in parts:
        b = open(f, "rb").read()
        acc = bytearray(b) if acc is None else bytearray(x | y for x, y in zip(acc, b))
    open(out, "wb").write(acc)
