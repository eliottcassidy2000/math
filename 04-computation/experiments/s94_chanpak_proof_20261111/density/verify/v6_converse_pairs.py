# Verifier (lane density, S94). The reverse-complement involution h = 1|x  ->  h* = 1|~rev(x) (x = the n-1 low bits)
# and its share in the lane's collision statistics (D-07, D-08, D-09, D-11).
# (1) exact check of eta(h*) = eta(h) on all n-bit words; (2) the lane's statistics with the pairs (h, h*) removed.
import sys, numpy as np
from numba import njit
@njit(cache=True)
def eta_array(n):
    base = 1 << (n - 1)
    out = np.empty(base, np.int64)
    for k in range(base):
        h = base + k
        qq = 0
        for i in range(n):
            if (h >> i) & 1:
                qq += ((h >> (i + 1)) + 1) >> 1
        out[k] = 2 * h + 1 - 2 * qq
    return out
@njit(cache=True)
def conv_index(n):
    base = 1 << (n - 1); w = n - 1
    out = np.empty(base, np.int64)
    for k in range(base):
        r = 0
        for i in range(w):
            if (k >> i) & 1: r |= 1 << (w - 1 - i)
        out[k] = base - 1 - r
    return out
@njit(cache=True)
def scan(n, order, eta, conv, keys):
    base = 1 << (n - 1); m = order.shape[0]
    histT = np.zeros(n + 1, np.int64); histN = np.zeros(n + 1, np.int64)   # by width: trivial pairs / the others
    i = 0; tot = 0; triv = 0
    while i < m:
        j = i + 1
        while j < m and eta[order[j]] == eta[order[i]]: j += 1
        for a in range(i, j):
            ka = order[a]
            for b in range(i, j):
                if a == b: continue
                kb = order[b]
                d = ka ^ kb
                top = 0
                while (d >> (top + 1)) != 0: top += 1
                s = 0
                while ((d >> s) & 1) == 0: s += 1
                w = top - s + 1; tp = n - 1 - top
                if conv[ka] == kb:
                    triv += 1; histT[w] += 1
                else:
                    histN[w] += 1
                    xa = (ka >> s) & ((1 << w) - 1); xb = (kb >> s) & ((1 << w) - 1)
                    keys[tot] = (np.int64(tp) << 56) | (np.int64(s) << 50) | (np.int64(xa & ~xb) << 25) | np.int64(xb & ~xa)
                    tot += 1
        i = j
    return tot, triv, histT, histN
def popc(a):
    a = a.copy(); c = np.zeros(a.shape, np.int64)
    while a.any():
        c += a & 1; a >>= 1
    return c
for n in [int(a) for a in sys.argv[1:]]:
    eta = eta_array(n); conv = conv_index(n); base = 1 << (n - 1)
    viol = int((eta[conv] != eta).sum()); fixed = int((conv == np.arange(base)).sum())
    print("n=%d: eta(h*) = eta(h) on all %d words: violations %d; fixed points of the involution: %d" % (n, base, viol, fixed))
    order = np.argsort(eta, kind='mergesort')
    keys = np.zeros(2 * base + 1000000, np.int64)
    tot, triv, histT, histN = scan(n, order, eta, conv, keys)
    allp = tot + triv
    print("  ordered colliding pairs %d = %d reverse-complement pairs (h,h*) + %d others; per word: %.4f = %.4f + %.4f"
          % (allp, triv, tot, allp / base, triv / base, tot / base))
    line = []
    for w in range(1, n):
        if histT[w] + histN[w]:
            dT = (n + w) / 2 - np.log2(histT[w]) if histT[w] else float('nan')
            dN = (n + w) / 2 - np.log2(histN[w]) if histN[w] else float('nan')
            line.append("w=%d: %d+%d [d_triv=%.2f d_other=%.2f]" % (w, histT[w], histN[w], dT, dN))
    print("  by width, (h,h*) pairs + others, d = (n+w)/2 - log2 count:", "; ".join(line))
    k, cnt = np.unique(keys[:tot], return_counts=True)
    tp = k >> 56; s = (k >> 50) & 63; plus = (k >> 25) & ((1 << 25) - 1); minus = k & ((1 << 25) - 1)
    w = n - tp - s; z = w - popc(plus | minus); m = np.minimum(tp, s)
    ratio = cnt / np.maximum(1.0, 2.0 ** (m + z - w))
    print("  with the (h,h*) pairs removed: classes with a solution %d, max sol %d, max ratio sol/max(1,2^(m+z-w)) = %.1f" % (k.size, int(cnt.max()), float(ratio.max())))
    out = []
    for ww in range(max(3, n - 10), n):
        sel = w == ww
        if sel.any(): out.append("w=%d: %d classes, mean %.3f, max %d, max ratio %.1f" % (ww, int(sel.sum()), cnt[sel].mean(), int(cnt[sel].max()), float(ratio[sel].max())))
    print("  solved classes by width (others only):", "; ".join(out))
    sys.stdout.flush()
