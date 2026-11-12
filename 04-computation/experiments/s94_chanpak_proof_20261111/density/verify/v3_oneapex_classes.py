# Verifier (lane density, S94). Independent recount of D-07, D-08, D-09, D-11.
# Independence: q(h) from the unrolled sum q(h) = sum_i h_i * ceil(floor(h/2^(i+1))/2) evaluated word by word (the lane
# used the array recursion); pair enumeration and class keys written from scratch.
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
                f = h >> (i + 1)
                qq += (f + 1) >> 1
        out[k] = 2 * h + 1 - 2 * qq
    return out
@njit(cache=True)
def scan(n, order, eta, keys):
    # returns total ordered colliding pairs, cell histogram, max fibre, number of distinct values; fills keys
    base = 1 << (n - 1); m = order.shape[0]
    hist = np.zeros((n + 1, n + 1), np.int64)
    i = 0; tot = 0; maxfib = 0; distinct = 0
    while i < m:
        j = i + 1
        while j < m and eta[order[j]] == eta[order[i]]: j += 1
        distinct += 1
        if j - i > maxfib: maxfib = j - i
        for a in range(i, j):
            ha = base + order[a]
            for b in range(i, j):
                if a == b: continue
                hb = base + order[b]
                d = ha ^ hb
                top = 0
                while (d >> (top + 1)) != 0: top += 1
                s = 0
                while ((d >> s) & 1) == 0: s += 1
                w = top - s + 1; tp = n - 1 - top
                hist[tp, s] += 1
                xa = (ha >> s) & ((1 << w) - 1); xb = (hb >> s) & ((1 << w) - 1)
                plus = xa & ~xb; minus = xb & ~xa
                if keys.shape[0] > tot:
                    keys[tot] = (np.int64(tp) << 56) | (np.int64(s) << 50) | (np.int64(plus) << 25) | np.int64(minus)
                tot += 1
        i = j
    return tot, hist, maxfib, distinct
def revw(a, w):
    # bit reversal of the low w bits, vectorised (w is an array)
    r = np.zeros(a.shape, np.int64)
    wm = int(w.max())
    for i in range(wm):
        bit = (a >> i) & 1
        sh = w - 1 - i
        r += np.where(sh >= 0, bit << np.maximum(sh, 0), 0)
    return r
def popc(a):
    a = a.copy(); c = np.zeros(a.shape, np.int64)
    while a.any():
        c += a & 1; a >>= 1
    return c
for n in [int(a) for a in sys.argv[1:]]:
    eta = eta_array(n); base = 1 << (n - 1)
    order = np.argsort(eta, kind='mergesort')
    keys = np.zeros(2 * base + 1000000, np.int64)
    tot, hist, maxfib, distinct = scan(n, order, eta, keys)
    assert tot <= keys.shape[0]
    print("n=%d words=%d ordered colliding pairs=%d (%.4f per word) distinct values=%d (%.4f of 2^(n-1)) largest fibre=%d"
          % (n, base, tot, tot / base, distinct, distinct / base, maxfib))
    # D-08: by width
    line = []
    for w in range(1, n):
        t = 0; best = 0; arg = None
        for tp in range(1, n):
            s = n - w - tp
            if s < 0: continue
            t += hist[tp, s]
            if hist[tp, s] > best: best = hist[tp, s]; arg = (tp, s)
        if t:
            line.append("w=%d:%d[d=%.2f,%s]" % (w, t, (n + w) / 2 - np.log2(t), arg))
    print("  pairs by width, d = (n+w)/2 - log2 pairs, dominant cell:", " ".join(line))
    # D-09 / D-11: classes
    k, cnt = np.unique(keys[:tot], return_counts=True)
    tp = k >> 56; s = (k >> 50) & 63; plus = (k >> 25) & ((1 << 25) - 1); minus = k & ((1 << 25) - 1)
    supp = plus | minus
    w = np.zeros(k.shape, np.int64); t = supp.copy()
    while t.any():
        w += (t > 0); t >>= 1
    assert np.all(w == n - tp - s)
    z = w - popc(supp)
    m = np.minimum(tp, s)
    X = plus - minus; Xi = revw(plus, w) - revw(minus, w)
    g = np.gcd(np.abs(X), np.abs(Xi)); mx = np.maximum(np.abs(X), np.abs(Xi))
    ratio = cnt / np.maximum(1.0, 2.0 ** (m + z - w))
    ratio2 = ratio / (1.0 + (2.0 ** m) * g / mx)                 # the 'corrected hypothesis' with the lattice factor
    ratio3 = ratio / (1.0 + (2.0 ** np.maximum(tp - 1, s)) * g / mx)  # with the factor of THM-5141 (C)
    print("  classes (cell, eps) with a solution: %d; max solutions in a class %d; max ratio sol/max(1,2^(m+z-w)) = %.1f;"
          " with lattice factor (1+2^m g/max): %.1f; with (1+2^max(t'-1,s) g/max): %.1f"
          % (k.size, int(cnt.max()), ratio.max(), ratio2.max(), ratio3.max()))
    for ww in (n - 1, n - 2):
        sel = w == ww
        eps_only = np.unique((plus[sel] << 25) | minus[sel]).size
        cells = np.unique((tp[sel] << 6) | s[sel]).size
        print("  width %d: (cell, eps) classes with a solution %d in %d cell(s); distinct eps %d; all odd classes per cell %d; pairs %d; mean sol %.2f"
              % (ww, int(sel.sum()), cells, eps_only, 4 * 3 ** (ww - 2), int(cnt[sel].sum()), cnt[sel].mean()))
    means = []
    for ww in range(1, n):
        sel = w == ww
        if sel.any(): means.append("w=%d:%.2f/%d" % (ww, cnt[sel].mean(), int(cnt[sel].max())))
    print("  mean/max solutions per solved class by width:", " ".join(means))
    # anatomy of the largest classes at the two top odd widths and at a small width
    def digits(p, mi, ww):
        return "".join("+" if (p >> i) & 1 else ("-" if (mi >> i) & 1 else "0") for i in range(ww - 1, -1, -1))
    kk = keys[:tot]
    for ww in sorted(set([n - 1, n - 3, 5, 3])):
        sel = np.where(w == ww)[0]
        if sel.size == 0: continue
        j = sel[np.argmax(cnt[sel])]
        pal = "X=Xi" if X[j] == Xi[j] else ("X=-Xi" if X[j] == -Xi[j] else "height>1" if mx[j] // g[j] > 1 else "?")
        print("  largest class at width %d: sol=%d cell (t',s)=(%d,%d) z=%d eps(top..bottom)=%s  X=%d Xi=%d gcd=%d height=%d (%s)  lattice bound per pair and parity class 1+g*min(2^(t'-1)/|X|,2^s/|Xi|) = %.3f"
              % (ww, int(cnt[j]), int(tp[j]), int(s[j]), int(z[j]), digits(int(plus[j]), int(minus[j]), ww), int(X[j]), int(Xi[j]), int(g[j]),
                 int(mx[j] // g[j]), pal, 1 + g[j] * min(2.0 ** (tp[j] - 1) / abs(X[j]), 2.0 ** s[j] / abs(Xi[j]))))
    # share of the solutions sitting in height-1 classes, and max sol among classes of height > 1
    h1 = (mx == g)
    print("  solutions in height-1 classes: %d of %d; max sol in a class of height > 1: %d (ratio %.1f)"
          % (int(cnt[h1].sum()), tot, int(cnt[~h1].max()), float(ratio[~h1].max())))
    hh = ~h1
    print("  classes of height > 1: max ratio %.1f; max ratio with lattice factor (1+2^m g/max) %.1f; with (1+2^max(t'-1,s) g/max) %.1f"
          % (float(ratio[hh].max()), float(ratio2[hh].max()), float(ratio3[hh].max())))
    # largest classes in the cell (1,0) (unique (P,c) = (1,2)): all their solutions are distinct common-bit vectors
    sel = np.where((tp == 1) & (s == 0))[0]
    if sel.size:
        j = sel[np.argmax(cnt[sel])]
        print("  cell (1,0) [P=1, c=2 forced]: largest class sol=%d, z=%d, eps=%s, height=%d" % (int(cnt[j]), int(z[j]), digits(int(plus[j]), int(minus[j]), n - 1), int(mx[j] // g[j])))
        top = sel[np.argsort(-cnt[sel])[:6]]
        print("     top classes in cell (1,0): " + "; ".join("sol=%d z=%d height=%d" % (int(cnt[u]), int(z[u]), int(mx[u] // g[u])) for u in top))
    sys.stdout.flush()
