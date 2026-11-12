# S94 density lane. One-apex collisions eta(h) = eta(h') among n-bit words, grouped by CELL (t', s) and by the
# signed-digit difference class eps of the windows (the class determines (X, Xi) and conversely, THM-5141 (A)2).
# For each class with at least one solution: sol = number of ordered colliding pairs (b, P, c) in the class.
# Test of the hypothesis H_class:  sol <= C * max(1, 2^(min(t',s) + z - w)),  z = number of zero digits of eps.
import sys, numpy as np
from numba import njit
@njit(cache=True)
def build_q(n):
    N = 1 << n
    q = np.zeros(N, np.int64)
    for h in range(1, N):
        q[h] = q[h >> 1] + ((h & 1) * (((h >> 1) + 1) >> 1))
    return q
@njit(cache=True)
def keys_of_pairs(n, order, eta, out):
    m = order.shape[0]; base = 1 << (n - 1); i = 0; k = 0
    while i < m:
        j = i; v = eta[order[i]]
        while j < m and eta[order[j]] == v: j += 1
        for a in range(i, j):
            ha = order[a] + base
            for b in range(i, j):
                if a == b: continue
                hb = order[b] + base
                d = ha ^ hb
                top = 0; dd = d
                while dd > 1: dd >>= 1; top += 1
                s = 0; dd = d
                while (dd & 1) == 0: dd >>= 1; s += 1
                w = top - s + 1
                mask = (1 << w) - 1
                plus = ((ha & ~hb) >> s) & mask
                minus = ((hb & ~ha) >> s) & mask
                tp = n - 1 - top
                out[k] = (((tp << 5) | s) << 52) | (plus << 26) | minus
                k += 1
        i = j
    return k
def popcount(a):
    a = a.copy(); c = np.zeros(a.shape, np.int64)
    while a.any():
        c += a & 1; a >>= 1
    return c
for n in [int(a) for a in sys.argv[1:]]:
    q = build_q(n); base = 1 << (n - 1)
    h = np.arange(base, 1 << n, dtype=np.int64)
    eta = 2 * h + 1 - 2 * q[base:]
    order = np.argsort(eta, kind='stable')
    out = np.zeros(3 * base, np.int64)
    k = keys_of_pairs(n, order, eta, out)
    keys, cnt = np.unique(out[:k], return_counts=True)
    tp = keys >> 57; s = (keys >> 52) & 31
    plus = (keys >> 26) & ((1 << 26) - 1); minus = keys & ((1 << 26) - 1)
    supp = plus | minus
    w = np.zeros(keys.shape, np.int64); t = supp.copy()
    while t.any():
        w += (t > 0); t >>= 1
    z = w - popcount(supp)
    m = np.minimum(tp, s)
    expo = m + z - w
    ratio = cnt / np.maximum(1.0, 2.0 ** expo)
    print("n=%d ordered colliding pairs=%d classes (cell, eps) with a solution=%d  max solutions in a class=%d  max of sol/max(1,2^(m+z-w))=%.2f"
          % (n, k, keys.size, int(cnt.max()), float(ratio.max())))
    print("  w : classes with a solution | all odd classes 4*3^(w-2) | pairs | mean sol per solved class | max sol | max ratio")
    for ww in range(1, n):
        sel = (w == ww)
        if sel.any():
            allc = 4 * 3 ** (ww - 2) if ww >= 2 else 2
            print("  %2d : %9d | %14d | %9d | %6.2f | %5d | %6.2f" % (ww, int(sel.sum()), allc, int(cnt[sel].sum()),
                  cnt[sel].mean(), int(cnt[sel].max()), float(ratio[sel].max())))
    # histogram of the ratio
    hh = np.bincount(np.minimum(np.floor(ratio).astype(np.int64), 20))
    print("  histogram of floor(ratio) (last bin >= 20):", hh.tolist())
    sys.stdout.flush()
