#!/usr/bin/env python3
"""auditC_ext_post.py -- audit C: post-processing of the brute-force extension spectra (ext<n>.out, lines
"E <tour> <h> <sorted distinct values of H(X+v) over all cuts>") for THM-4876 (5) and HYP-9469 (1), (2).

Per order n: k_min(h); increment-2 cases by h; longest window of consecutive odd values in E(X) (max over X, and
whether any window covers [2h+1, 3h]); whether any X has E(X) >= odd[h, 2h-1]; number of X with E(X) an odd interval;
cross-check of k_min against the deletion census del<n+1>.out (KMIN line).
"""
import os, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def longest_window(E):
    best = (0, None, None)
    start = prev = E[0]
    for x in E[1:] + [None]:
        if x is not None and x == prev + 2:
            prev = x
            continue
        L = (prev - start) // 2 + 1
        if L > best[0]:
            best = (L, start, prev)
        if x is not None:
            start = prev = x
    return best


def main():
    for n in range(3, 9):
        fn = os.path.join(HERE, "ext%d.out" % n)
        if not os.path.exists(fn):
            continue
        kmin = {}
        inc2 = Counter()
        nX = 0
        maxwin = (0, None, None, None)
        cover_2h3h = 0
        cover_h2h = 0
        intervals = 0
        min_ratio_q_p = None
        min_max_ratio = None
        by_h_best_L = defaultdict(int)
        with open(fn) as f:
            for line in f:
                if not line.startswith("E "):
                    continue
                parts = line.split()
                h = int(parts[2])
                E = [int(x) for x in parts[3:]]
                assert E[0] == h and E == sorted(set(E))
                nX += 1
                pos = [x for x in E if x > h]
                if pos:
                    k = (pos[0] - h) // 2
                    kmin[h] = min(kmin.get(h, 10 ** 9), k)
                    if k == 1:
                        inc2[h] += 1
                L, p, q = longest_window(E)
                if L > maxwin[0]:
                    maxwin = (L, p, q, h)
                if p is not None and p <= 2 * h + 1 and q >= 3 * h:
                    cover_2h3h += 1
                if all(x in set(E) for x in range(h, 2 * h, 2)):
                    cover_h2h += 1
                if E == list(range(h, E[-1] + 1, 2)):
                    intervals += 1
                r = (q - p) / h
                if min_ratio_q_p is None or r < min_ratio_q_p:
                    min_ratio_q_p = r
                mr = E[-1] / h
                if min_max_ratio is None or mr < min_max_ratio:
                    min_max_ratio = mr
                # largest L with every odd in [h, L] in E
                Lh = h
                S = set(E)
                while Lh + 2 in S:
                    Lh += 2
                by_h_best_L[h] = max(by_h_best_L[h], Lh)
        hs_good = [h for h in by_h_best_L if by_h_best_L[h] >= 2 * h - 1]
        print("order %d (X): %d tournaments; increment-2 cases by h: %s; k_min(h): %s" % (n, nX, dict(sorted(inc2.items())), {h: kmin[h] for h in sorted(kmin) if h <= 60}))
        print("   longest window of consecutive odd values: length %d, [%s, %s] at h=%s; X with a window covering [2h+1,3h]: %d; X with E(X) >= odd[h,2h-1]: %d; h with SOME realiser X having E(X) >= odd[h,2h-1]: %s; X with E(X) an odd interval: %d; min (q-p)/h = %.3f; min max(E)/h = %.3f" % (
            maxwin[0], maxwin[1], maxwin[2], maxwin[3], cover_2h3h, cover_h2h, hs_good, intervals, min_ratio_q_p, min_max_ratio))
        if 5 in kmin:
            print("   k_min(h)/h for h >= 5: min %.4f max %.4f; all k_min(h) >= h/4: %s" % (min(kmin[h] / h for h in kmin if h >= 5), max(kmin[h] / h for h in kmin if h >= 5), all(4 * kmin[h] >= h for h in kmin if h >= 5)))
        # cross-check against the deletion census of order n+1
        fd = os.path.join(HERE, "del%d.out" % (n + 1))
        if os.path.exists(fd):
            with open(fd) as f:
                for line in f:
                    if line.startswith("KMIN"):
                        kd = {int(a): int(b) for a, b in (t.split(":") for t in line.split()[1:])}
                        print("   k_min from the deletion census at order %d agrees: %s" % (n + 1, kd == kmin))
                        if kd != kmin:
                            print("      del:", kd, "\n      ext:", kmin)


if __name__ == "__main__":
    main()
