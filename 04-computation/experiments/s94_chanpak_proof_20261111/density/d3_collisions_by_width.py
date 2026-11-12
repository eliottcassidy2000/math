# S94 density lane. One-apex collisions eta(h) = eta(h') among n-bit words (top bit 1), classified by the window
# width w (first to last differing bit) and by cell (t', s).  Compares the true number of colliding ordered pairs of
# width w with the bound 4*4^w per cell of THM-5141 (D) and with the heuristic 2^((n+w)/2).
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
def classify(n, order, eta):
    # order: indices (h - 2^(n-1)) sorted by eta
    hist = np.zeros((n + 1, n + 1), np.int64)   # [t', s]
    m = order.shape[0]
    base = 1 << (n - 1)
    i = 0
    npairs = 0
    maxfib = 0
    while i < m:
        j = i
        v = eta[order[i]]
        while j < m and eta[order[j]] == v:
            j += 1
        if j - i > maxfib: maxfib = j - i
        for a in range(i, j):
            ha = order[a] + base
            for b in range(i, j):
                if a == b: continue
                d = ha ^ (order[b] + base)
                top = 0
                dd = d
                while dd > 1:
                    dd >>= 1; top += 1
                s = 0
                dd = d
                while (dd & 1) == 0:
                    dd >>= 1; s += 1
                hist[n - 1 - top, s] += 1
                npairs += 1
        i = j
    return hist, npairs, maxfib
for n in [int(a) for a in sys.argv[1:]]:
    q = build_q(n)
    base = 1 << (n - 1)
    h = np.arange(base, 1 << n, dtype=np.int64)
    eta = 2 * h + 1 - 2 * q[base:]
    order = np.argsort(eta, kind='stable')
    hist, npairs, maxfib = classify(n, order, eta)
    print("n=%d words=%d ordered colliding pairs=%d (%.4f per word) distinct values=%d largest fibre=%d"
          % (n, base, npairs, npairs / base, len(np.unique(eta)), maxfib))
    print("  w : pairs(w)   log2 pairs   (n+w)/2   2w   | max cell (t',s) count, 4*4^w")
    for w in range(1, n):
        tot = 0; mx = 0; arg = (0, 0)
        for tp in range(1, n):
            s = n - w - tp
            if s < 0: continue
            c = hist[tp, s]; tot += c
            if c > mx: mx = c; arg = (tp, s)
        if tot:
            print("  %2d : %10d  %7.2f  %6.1f  %3d  | %d at %s, %d" % (w, tot, np.log2(tot), (n + w) / 2, 2 * w, mx, arg, 4 * 4 ** w))
    sys.stdout.flush()
