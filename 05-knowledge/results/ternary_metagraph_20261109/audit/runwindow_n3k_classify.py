"""Classify the n = 3k run-window outcomes of runwindow_n3k.out by the lengths of the two 1-runs around each
isolated 0 (pattern 1^a 0 1^b with a, b maximal)."""
import re
import collections

stats = collections.Counter()
for line in open("runwindow_n3k.out"):
    m = re.match(r"\s+w=(\d+) isolated-0-between-1-runs=(\d) run-window realises=(\d)", line)
    if not m:
        continue
    w, iso, real = m.group(1), int(m.group(2)), int(m.group(3))
    if not iso:
        stats[("no isolated 0", real)] += 1
        continue
    L = len(w)
    worst = 0     # max over isolated 0s of min(left 1-run length, right 1-run length)
    for a in range(1, L - 1):
        if w[a] == "0" and w[a - 1] == "1" and w[a + 1] == "1":
            l = a - 1
            while l - 1 >= 0 and w[l - 1] == "1":
                l -= 1
            r = a + 1
            while r + 1 < L and w[r + 1] == "1":
                r += 1
            worst = max(worst, min(a - l, r - a))
    stats[("isolated 0, min neighbouring 1-run length >= 3" if worst >= 3 else
           "isolated 0, min neighbouring 1-run length " + str(worst), real)] += 1
for key in sorted(stats):
    print(f"{key[0]:55s} realised={key[1]}: {stats[key]}")
