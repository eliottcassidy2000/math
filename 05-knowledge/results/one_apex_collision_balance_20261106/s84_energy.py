"""S84: energy and collision structure of one-apex counts on reduced words Y_n (top bit 1, bottom bit 0).

g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2).  Y_n = even v in [2^(n-1), 2^n).
Outputs, per n:
  E^Y/|Y|, converse and non-converse parts, multiplicity histogram of values,
  centred-window statistics of non-converse collision pairs:
    a = common prefix length, s = common suffix length, w = window width,
    dX = w - bitlen|x-x'|, dXi = w - bitlen|C(x)-C(x')| (carry depths),
    checks of the proved lemma (|a-s| <= 2 + max(dX,dXi) or a <= 4+dXi and s <= 4+dX)
    and of the sharper empirical form |a-s| <= 1 + max(dX,dXi);
  prefix-level law: non-converse ordered pairs with prefix length a, per |Y|.
Usage: python s84_energy.py n_min n_max
"""
import sys, time
import numpy as np
from numba import njit


@njit(cache=True)
def gtable(N):
    g = np.empty(1 << N, dtype=np.int64)
    g[0] = 1
    size = 1
    for _ in range(N):
        for v in range(size - 1, -1, -1):
            gv = g[v]
            g[2 * v + 1] = gv
            g[2 * v] = gv + v + (v & 1)
        size *= 2
    g[0] = 1
    return g


@njit(cache=True)
def rc(x, n):
    r = 0
    for i in range(n):
        r |= (1 - ((x >> i) & 1)) << (n - 1 - i)
    return r


@njit(cache=True)
def bitlen(x):
    b = 0
    while x > 0:
        b += 1
        x >>= 1
    return b


@njit(cache=True)
def analyse(vs, gs, n, out_hist, out_lvl, out_dev, out_flags):
    # vs sorted by gs; gs sorted
    m = len(vs)
    i = 0
    conv = 0
    nonconv = 0
    viol_proved = 0
    viol_sharp = 0
    while i < m:
        j = i + 1
        while j < m and gs[j] == gs[i]:
            j += 1
        r = j - i
        if r < out_hist.shape[0]:
            out_hist[r] += 1
        else:
            out_hist[out_hist.shape[0] - 1] += 1
        if r >= 2:
            for p in range(i, j):
                for q in range(i, j):
                    if p == q:
                        continue
                    v = vs[p]
                    w2 = vs[q]
                    if rc(v, n) == w2:
                        conv += 1
                        continue
                    nonconv += 1
                    x = v ^ w2
                    t = bitlen(x) - 1
                    low = x & (-x)
                    s = bitlen(low) - 1
                    a = n - 1 - t
                    w = t - s + 1
                    mask = (1 << w) - 1
                    xv = (v >> s) & mask
                    xw = (w2 >> s) & mask
                    dX = w - bitlen(abs(xv - xw))
                    cv = rc(xv, w)
                    cw = rc(xw, w)
                    dXi = w - bitlen(abs(cv - cw))
                    dm = max(dX, dXi)
                    ok_proved = (abs(a - s) <= 2 + dm) or (a <= 4 + dXi and s <= 4 + dX)
                    if not ok_proved:
                        viol_proved += 1
                    if abs(a - s) > 1 + dm:
                        viol_sharp += 1
                    if a < out_lvl.shape[0]:
                        out_lvl[a] += 1
                    k = a - s + 20
                    if 0 <= k < out_dev.shape[0] and dm < out_dev.shape[1]:
                        out_dev[k, dm] += 1
        i = j
    out_flags[0] = conv
    out_flags[1] = nonconv
    out_flags[2] = viol_proved
    out_flags[3] = viol_sharp


def run(n, g):
    t0 = time.time()
    lo, hi = 1 << (n - 1), 1 << n
    v = np.arange(lo, hi, 2, dtype=np.int64)
    gv = g[v]
    order = np.argsort(gv, kind='stable')
    vs = v[order]
    gs = gv[order]
    hist = np.zeros(64, dtype=np.int64)
    lvl = np.zeros(64, dtype=np.int64)
    dev = np.zeros((41, 16), dtype=np.int64)
    flags = np.zeros(4, dtype=np.int64)
    analyse(vs, gs, n, hist, lvl, dev, flags)
    Y = len(v)
    E = sum(r * r * hist[r] for r in range(64))
    conv, nonconv, vp, vsh = flags
    # smoothed energies E^Y(K) = #{(v,v'): |g(v)-g(v')| <= K}
    smooth = {}
    for K in (2, 4, 8, 16):
        hi_ = np.searchsorted(gs, gs + K, side='right')
        lo_ = np.searchsorted(gs, gs - K, side='left')
        EK = int((hi_ - lo_).sum())
        smooth[K] = round((EK - E) / (K * Y), 5)
    out = {
        'n': n, '|Y|': Y, 'E/|Y|': round(E / Y, 5), 'conv/|Y|': round(conv / Y, 5),
        'nonconv/|Y|': round(nonconv / Y, 5),
        'viol_proved': int(vp), 'viol_sharp(+1)': int(vsh),
        'max_mult': int(max(r for r in range(64) if hist[r] > 0)),
        'values/|Y|': round(int(hist.sum()) / Y, 5),
        'mult_frac(even 2q)': {int(r): round(hist[r] / Y, 6) for r in range(2, 64, 2) if hist[r] > 0},
        'mult_hist': {int(r): int(hist[r]) for r in range(64) if hist[r] > 0},
        'level_per_|Y|': {int(a): round(lvl[a] / Y, 5) for a in range(64) if lvl[a] > 0},
        'slope (E(K)-E(0))/(K|Y|)': smooth,
        'secs': round(time.time() - t0, 1),
    }
    return out, dev


def layers(n, g):
    """X_n energy split into same-layer (Y_{n-k} with k trailing ones) and cross-layer parts."""
    vals = []
    lab = []
    same = 0
    for k in range(0, n - 1):
        m = n - k
        if m < 2:
            break
        v = np.arange(1 << (m - 1), 1 << m, 2, dtype=np.int64)
        gv = np.sort(g[v])
        _, c = np.unique(gv, return_counts=True)
        same += int((c.astype(np.int64) ** 2).sum())
        vals.append(gv)
    allv = np.concatenate(vals + [np.array([1], dtype=np.int64)])  # the all-ones word has g = 1
    _, c = np.unique(allv, return_counts=True)
    En = int((c.astype(np.int64) ** 2).sum())
    X = 1 << (n - 1)
    return {'n': n, 'E/|X|': round(En / X, 5), 'same-layer/|X|': round(same / X, 5),
            'cross-layer/|X|': round((En - same - 1) / X, 5)}


if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2])
    g = gtable(n1)
    pred = 14 - 16 * np.log(2)
    do_layers = len(sys.argv) > 3 and sys.argv[3] == 'layers'
    for n in range(n0, n1 + 1):
        out, dev = run(n, g)
        print(out)
        print('   Poisson-pair prediction 14-16ln2 = %.5f ; observed - predicted = %+.5f' % (pred, out['E/|Y|'] - pred))
        if do_layers:
            print('  ', layers(n, g), ' (independent-layer model 6-4ln2 = %.5f)' % (6 - 4 * np.log(2)))
        sys.stdout.flush()
