"""Auditor: at n = 3k (6..18), for every tight colouring (x_0 = 1), the lengths of the 1-runs on both sides of each
isolated 0; joined with the run-window outcome from runwindow_n3k.out.  Tests the characterisation
'run-window fails iff some isolated 0 has both neighbouring 1-runs of length >= 2'."""
import re
import collections


def P(xs):
    c, last = 0, -10
    for a in xs:
        if a - last >= 3:
            c, last = c + 1, a
    return c


real = {}
for line in open("runwindow_n3k.out"):
    m = re.match(r"\s+w=(\d+) isolated-0-between-1-runs=(\d) run-window realises=(\d)", line)
    if m:
        real[m.group(1)] = int(m.group(3))
tab = collections.Counter()
for n in range(6, 19, 3):
    L, k = n - 2, n // 3
    for x in range(1 << L):
        if not x & 1:
            continue
        w = "".join("1" if x >> i & 1 else "0" for i in range(L))
        A = [i for i in range(L) if w[i] == "1"]
        B = [i for i in range(L) if w[i] == "0"]
        if P(A) + P(B) != k:
            continue
        r = real.get(w, 1)          # words not listed in runwindow_n3k.out are realised and have no isolated 0
        pairs = []
        for a in range(1, L - 1):
            if w[a] == "0" and w[a - 1] == "1" and w[a + 1] == "1":
                l = a - 1
                while l - 1 >= 0 and w[l - 1] == "1":
                    l -= 1
                rr = a + 1
                while rr + 1 < L and w[rr + 1] == "1":
                    rr += 1
                pairs.append((a - l, rr - a))
        both2 = any(min(p) >= 2 for p in pairs)
        has2 = any(2 in p for p in pairs)
        tab[(n, both2, r)] += 1
        if has2:
            tab[(n, "a neighbouring 1-run of length exactly 2")] += 1
for key in sorted(tab, key=str):
    print(key, tab[key])
