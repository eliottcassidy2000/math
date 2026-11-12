# cells.py  (lane energy, S94)
# Exact collision statistics of the one-apex counts eta(h) = 2h + 1 - 2 q(h) on n-bit words (top bit 1).
# For every unordered colliding pair h != h' (eta(h) = eta(h')): t' = common prefix length, s = common suffix
# length, w = n - t' - s, window pair (x, x').
# Output: (a) number of unordered colliding pairs in each cell (t', s);
#         (b) for each cell the number of distinct unordered window pairs {x, x'} that occur
#             (= odd window pairs admitting at least one solution (P, c) in that cell);
#         (c) for each width w: collisions, distinct window pairs over all cells of that width, and 4^(w-1)/2... totals;
#         (d) Coll(t) = ordered colliding pairs with common prefix >= t, per word.
# All arithmetic is exact (int64; eta < 2^(n+1)).
import sys, time
import numpy as np
from numba import njit

@njit(cache=True)
def eta_table(n):
    N = 1 << n
    e = np.empty(N, dtype=np.int64)
    e[0] = 1
    for h in range(1, N):
        p = h >> 1
        if h & 1:
            e[h] = e[p] + p + 2 - (p & 1)
        else:
            e[h] = e[p] + 2 * p
    return e

@njit(cache=True)
def count_pairs(vals, idx, n, cell, keys):
    # vals sorted, idx = words; returns number of pairs, fills cell[t'][s], keys[...] (packed t', s, x, x')
    m = vals.shape[0]
    npairs = 0
    i = 0
    while i < m:
        j = i + 1
        while j < m and vals[j] == vals[i]:
            j += 1
        if j - i > 1:
            for a in range(i, j):
                for b in range(a + 1, j):
                    h = idx[a]; hp = idx[b]
                    d = h ^ hp
                    # bit length of d
                    bl = 0
                    dd = d
                    while dd > 0:
                        bl += 1
                        dd >>= 1
                    tp = n - bl
                    s = 0
                    dd = d
                    while (dd & 1) == 0:
                        s += 1
                        dd >>= 1
                    w = n - tp - s
                    cell[tp, s] += 1
                    x = (h >> s) & ((1 << w) - 1)
                    xp = (hp >> s) & ((1 << w) - 1)
                    if x > xp:
                        x, xp = xp, x
                    if npairs < keys.shape[0]:
                        keys[npairs] = (np.int64(tp) << 58) | (np.int64(s) << 53) | (np.int64(x) << 26) | np.int64(xp)
                    npairs += 1
        i = j
    return npairs

def run(n, out):
    t0 = time.time()
    e = eta_table(n)
    lo = 1 << (n - 1)
    vals = e[lo:].copy()
    del e
    order = np.argsort(vals, kind='stable')
    vals = vals[order]
    idx = (order + lo).astype(np.int64)
    del order
    distinct = int(np.count_nonzero(np.diff(vals)) + 1)
    cap = int(3.0 * (1 << (n - 1))) + 1000
    keys = np.zeros(cap, dtype=np.int64)
    cell = np.zeros((n + 1, n + 1), dtype=np.int64)
    npairs = count_pairs(vals, idx, n, cell, keys)
    assert npairs <= cap, (npairs, cap)
    keys = keys[:npairs]
    words = 1 << (n - 1)
    print(f"n={n} words={words} distinct={distinct} ratio={distinct/words:.5f} unordered_pairs={npairs} pairs/word={npairs/words:.5f} time={time.time()-t0:.1f}s", file=out)
    # distinct window pairs per cell
    keys.sort()
    uniq, first = np.unique(keys, return_index=True)
    ucell = np.zeros((n + 1, n + 1), dtype=np.int64)
    tps = (uniq >> 58) & 31
    ss = (uniq >> 53) & 31
    np.add.at(ucell, (tps, ss), 1)
    # per width: distinct window pairs over all cells of that width
    wk = (n - tps - ss).astype(np.int64)
    key2 = (wk << 53) | (uniq & ((np.int64(1) << 53) - 1))
    u2 = np.unique(key2)
    wd = np.bincount((u2 >> 53).astype(np.int64), minlength=n + 1)
    print("cells: t' s w  pairs  distinct_window_pairs  pairs/2^(n-max(t',s))  distinct/4^(w-1)*2", file=out)
    for tp in range(1, n):
        for s in range(0, n):
            c = int(cell[tp, s])
            if c:
                w = n - tp - s
                tot_odd_unordered = (4 ** (w - 1)) // 2 if w >= 2 else 1
                print(f"  {tp:2d} {s:2d} {w:2d} {c:10d} {int(ucell[tp, s]):10d} {c / 2.0 ** (n - max(tp, s)):9.5f} {int(ucell[tp, s]) / tot_odd_unordered:10.3e}", file=out)
    print("by width: w  pairs  distinct_window_pairs_all_cells  pairs/2^n  distinct/odd_pairs  distinct/2^((n+w)/2)", file=out)
    for w in range(1, n):
        c = 0
        for tp in range(1, n):
            s = n - w - tp
            if s >= 0:
                c += int(cell[tp, s])
        if c:
            tot_odd_unordered = (4 ** (w - 1)) // 2 if w >= 2 else 1
            print(f"  {w:2d} {c:10d} {int(wd[w]):10d} {c / 2.0 ** n:9.5f} {int(wd[w]) / tot_odd_unordered:10.3e} {int(wd[w]) / 2.0 ** ((n + w) / 2):9.5f}", file=out)
    print("by prefix: t  ordered pairs with common prefix >= t, per word, times 2^t", file=out)
    cum = 0
    rows = []
    for tp in range(n - 1, 0, -1):
        cum += int(cell[tp].sum())
        rows.append((tp, cum))
    for tp, cum in reversed(rows):
        print(f"  {tp:2d} {2 * cum:10d} {2 * cum / words:10.6f} {2 * cum / words * 2.0 ** tp:9.4f}", file=out)
    print("by offset d = s - t': pairs", file=out)
    off = {}
    for tp in range(1, n):
        for s in range(0, n):
            c = int(cell[tp, s])
            if c:
                off[s - tp] = off.get(s - tp, 0) + c
    for d in sorted(off):
        print(f"  {d:3d} {off[d]:10d} {off[d] / npairs:9.5f}", file=out)
    out.flush()

if __name__ == "__main__":
    n0 = int(sys.argv[1]); n1 = int(sys.argv[2])
    for n in range(n0, n1 + 1):
        run(n, sys.stdout)
