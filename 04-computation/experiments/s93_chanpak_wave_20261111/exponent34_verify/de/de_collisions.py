#!/usr/bin/env python3
"""Independent verifier (part D/E) - exact collision census of the one-apex counts x(h).

For every n in [n0, n1] enumerate ALL unordered pairs h != h' of n-bit words (top bit 1) with x(h) = x(h')
and record, for each pair:
  t' = length of the longest common prefix, s = length of the longest common suffix, w = n - t' - s,
  rr = max(longest run of equal bits in h, longest run in h').
Output: counts cnt[t'][s][rr] (unordered pairs), plus per-pair checks of
  * oddness of X, Xi; the collision equation c X - P Xi = gamma; |gamma| <= 2^w        (statement (C))
  * sign(X) = sign(Xi)                                                                 (THM-4890 (C) 1 transported)
  * |t' - s - 1| <= 1 + |dX - dXi|                                                     (THM-4890 (C) 4 transported)
  * |s - t'| <= rho + 2 with rho = max(dX, dXi) + 1                                     (Carries paragraph of (C))
  * max |s - t'| among pairs whose two words have no run longer than R, for every R     (|s - t'| <= R + 3)
Everything is exact integer arithmetic.  Own code; nothing imported from the coordinator.
"""
import sys, json, time
import numpy as np
from numba import njit


def q_table(nmax):
    """q(0)=0, q(2h)=q(h), q(2h+1)=q(h)+ceil(h/2), for h < 2^nmax."""
    N = 1 << nmax
    q = np.zeros(N, dtype=np.int64)
    for k in range(nmax):
        lo, hi = 1 << k, 1 << (k + 1)
        h = np.arange(lo, hi, dtype=np.int64)
        half = h >> 1
        q[lo:hi] = q[half] + (h & 1) * ((half + 1) >> 1)
    return q


def q_slow(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return q_slow(h // 2)
    return q_slow(h // 2) + (h // 2 + 1) // 2


def g_slow(v):
    """THM-4890 normalisation: g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2)."""
    if v == 0:
        return 1
    if v % 2 == 1:
        return g_slow(v // 2)
    u = v // 2
    return g_slow(u) + u + (u % 2)


@njit(cache=True)
def bitlen(v):
    b = 0
    while v > 0:
        b += 1
        v >>= 1
    return b


@njit(cache=True)
def rev_bits(v, w):
    r = 0
    for i in range(w):
        r = (r << 1) | ((v >> i) & 1)
    return r


@njit(cache=True)
def max_run(h, n):
    best = 1
    cur = 1
    prev = h & 1
    for i in range(1, n):
        b = (h >> i) & 1
        if b == prev:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
            prev = b
    return best


@njit(cache=True)
def census(n, q):
    lo = 1 << (n - 1)
    hi = 1 << n
    m = hi - lo
    xs = np.empty(m, dtype=np.int64)
    mr = np.empty(m, dtype=np.int64)
    for i in range(m):
        h = lo + i
        xs[i] = 2 * h + 1 - 2 * q[h]
        mr[i] = max_run(h, n)
    order = np.argsort(xs, kind='mergesort')
    cnt = np.zeros((n + 1, n + 1, n + 1), dtype=np.int64)
    stats = np.zeros(16, dtype=np.int64)
    # stats: 0 npairs, 1 viol_odd, 2 viol_eq, 3 viol_gamma, 4 viol_sign, 5 viol_bal4890, 6 viol_balhere,
    #        7 max(|s-t'| - max(dX,dXi)), 8 max |gamma| numerator at the pair maximising |gamma|/2^w (unused), 9 max group size
    #        10 number of pairs with |gamma| == 2^w, 11 max over pairs of (|t'-s-1| - |dX-dXi|)
    stats[7] = -1000
    stats[11] = -1000
    maxdiff_by_run = np.zeros(n + 1, dtype=np.int64) - 1
    maxgam_num = 0
    maxgam_w = 1
    i = 0
    while i < m:
        j = i + 1
        while j < m and xs[order[j]] == xs[order[i]]:
            j += 1
        if j - i > stats[9]:
            stats[9] = j - i
        if j - i >= 2:
            for a in range(i, j):
                h = lo + order[a]
                for b in range(a + 1, j):
                    h2 = lo + order[b]
                    d = h ^ h2
                    tp = n - bitlen(d)
                    s = 0
                    while ((d >> s) & 1) == 0:
                        s += 1
                    w = n - tp - s
                    rr = mr[order[a]]
                    if mr[order[b]] > rr:
                        rr = mr[order[b]]
                    cnt[tp, s, rr] += 1
                    stats[0] += 1
                    mask = (1 << w) - 1
                    x1 = (h >> s) & mask
                    x2 = (h2 >> s) & mask
                    P = h >> (w + s)
                    Z = h & ((1 << s) - 1)
                    c = (1 << (s + 1)) - rev_bits(Z, s)
                    X = x1 - x2
                    Xi = rev_bits(x1, w) - rev_bits(x2, w)
                    if (X & 1) == 0 or (Xi & 1) == 0:
                        stats[1] += 1
                    gam = 2 * (q[x1] - q[x2]) + (P & 1) * (((x1 >> (w - 1)) & 1) - ((x2 >> (w - 1)) & 1)) \
                        + (c & 1) * ((x1 & 1) - (x2 & 1))
                    if c * X - P * Xi != gam:
                        stats[2] += 1
                    ag = abs(gam)
                    if ag > (1 << w):
                        stats[3] += 1
                    if ag == (1 << w):
                        stats[10] += 1
                    # compare ag / 2^w with maxgam_num / 2^maxgam_w
                    if ag * (1 << maxgam_w) > maxgam_num * (1 << w):
                        maxgam_num = ag
                        maxgam_w = w
                    if (X > 0) != (Xi > 0):
                        stats[4] += 1
                    dX = w - bitlen(abs(X))
                    dXi = w - bitlen(abs(Xi))
                    dd = abs(dX - dXi)
                    if abs(tp - s - 1) > 1 + dd:
                        stats[5] += 1
                    if abs(tp - s - 1) - dd > stats[11]:
                        stats[11] = abs(tp - s - 1) - dd
                    mxd = dX if dX > dXi else dXi
                    if abs(s - tp) > mxd + 3:
                        stats[6] += 1
                    if abs(s - tp) - mxd > stats[7]:
                        stats[7] = abs(s - tp) - mxd
                    if abs(s - tp) > maxdiff_by_run[rr]:
                        maxdiff_by_run[rr] = abs(s - tp)
        i = j
    return cnt, stats, maxdiff_by_run, maxgam_num, maxgam_w


def main():
    n0 = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    n1 = int(sys.argv[2]) if len(sys.argv) > 2 else 22
    out_json = sys.argv[3] if len(sys.argv) > 3 else 'de_collisions.json'
    t0 = time.time()
    q = q_table(n1)
    # sanity of the recursion table and of x(h) = g(2h)
    for h in list(range(0, 3000)) + [12345, 54321, 99991, (1 << n1) - 1, (1 << (n1 - 1)) + 77]:
        assert q[h] == q_slow(h)
        assert 2 * h + 1 - 2 * q[h] == g_slow(2 * h), h
    hh = np.arange(1, 1 << n1, dtype=np.int64)
    assert np.all(q[1:] >= 0) and np.all(2 * q[1:] <= hh - 1)          # 0 <= q(h) <= (h-1)/2
    xfull = 2 * np.arange(1 << n1, dtype=np.int64) + 1 - 2 * q
    assert np.all(xfull[1:] >= hh + 2) and np.all(xfull[1:] <= 2 * hh + 1)
    assert not np.any((xfull == 7) | (xfull == 21)), "7 or 21 attained?!"
    print("table checks passed: q recursion, x(h) = g(2h) on samples, 0 <= q(h) <= (h-1)/2, h+2 <= x(h) <= 2h+1, "
          "x never 7 or 21 for h < 2^%d" % n1)
    res = {}
    for n in range(n0, n1 + 1):
        cnt, stats, mdr, mg, mw = census(n, q)
        Rn = int(np.ceil(np.log2(n))) + 1
        # cumulative max |s-t'| for words with no run > R
        cum = []
        cur = -1
        for r in range(n + 1):
            cur = max(cur, int(mdr[r]))
            cum.append(cur)
        entry = {
            'n': n, 'unordered_pairs': int(stats[0]), 'viol_odd': int(stats[1]), 'viol_eq': int(stats[2]),
            'viol_gamma': int(stats[3]), 'viol_sign': int(stats[4]), 'viol_bal4890': int(stats[5]),
            'viol_balhere': int(stats[6]), 'max_absdiff_minus_maxd': int(stats[7]),
            'max_abs_tp_minus_s_minus_1_minus_dd': int(stats[11]),
            'max_group': int(stats[9]), 'pairs_gamma_eq_2w': int(stats[10]),
            'max_gamma_ratio': [int(mg), int(mw)],
            'max_abs_s_minus_tp_for_runfree_R': cum,
            'cnt': {}
        }
        nz = np.argwhere(cnt.sum(axis=2) > 0)
        for tp, s in nz:
            entry['cnt']['%d,%d' % (tp, s)] = [int(v) for v in cnt[tp, s, :]]
        res[str(n)] = entry
        viol_R = [(R, cum[R]) for R in range(1, n + 1) if cum[R] > R + 3]
        tight = max((cum[R] - R) for R in range(1, n + 1) if cum[R] >= 0)
        print("n=%2d unordered colliding pairs=%9d  viol(odd,eq,gamma,sign,bal4890,balhere)=(%d,%d,%d,%d,%d,%d)  "
              "max|gamma|/2^w=%d/2^%d=%.4f  #|gamma|=2^w: %d  max(|s-t'|-max(dX,dXi))=%d  max(|t'-s-1|-|dX-dXi|)=%d  "
              "max group=%d" % (n, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5], stats[6], mg, mw,
                                mg / 2.0 ** mw, stats[10], stats[7], stats[11], stats[9]))
        print("      run-free: max|s-t'| for R=1..12: %s ; violations of |s-t'|<=R+3: %s ; max over R of (max|s-t'| - R) = %d ; "
              "theorem R=%d: max|s-t'|=%d" % (cum[1:13], viol_R, tight, Rn, cum[Rn]))
        sys.stdout.flush()
    with open(out_json, 'w') as f:
        json.dump(res, f)
    print("written", out_json, "time %.1fs" % (time.time() - t0))


if __name__ == '__main__':
    main()
