#!/usr/bin/env python3
"""Skeptic check 4: parts (C), (D), (E) of THM-5141 (draft) on ALL colliding pairs of n-bit words.

For each n: x(h) for every h in [2^(n-1), 2^n); every unordered pair h != h' with x(h) = x(h'):
  * t' (common prefix), s (common suffix), w = n - t' - s, windows x, x';
  * oddness of X, Xi; the collision equation c X - P Xi = gamma with the stated gamma; |gamma| <= 2^w;
  * the carry bound |s - t'| <= rho + 2 with rho minimal such that 2^(w-rho) <= |X|, |Xi|;
  * the run bound |s - t'| <= R' + 3 where R' = longest run of equal bits in h or h' (so for every R >= R');
  * histogram of (t', s) for all pairs and for pairs inside W (R = ceil(log2 n) + 1);
  * (D): ordered pairs with parameters (t', s) <= 4 * #odd(w) + 4 * 2^max(t'-1,s) S_w   (exact S_w);
         E = ordered pairs in W with common prefix >= t = ceil(n/4) versus the explicit sum of (D);
  * (E): number of distinct values on W versus |W|^2/(2^(t-1)(|W| + E)); |W| >= 2^(n-2).
  * (n <= NKEY) per window pair: number of (P, c) <= 4 (1 + 2^max(t'-1,s)/M).
"""
import sys, math, time, os
import numpy as np
from numba import njit, prange, set_num_threads

set_num_threads(4)

@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

@njit(cache=True)
def revw(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r

@njit(cache=True)
def maxrun(h, n):
    best = 1; cur = 1
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
def bitlen(v):
    l = 0
    while v:
        l += 1
        v >>= 1
    return l

@njit(cache=True)
def scan_pairs(n, order, vals, q, R, keys, want_keys):
    """order: indices (h - 2^(n-1)) sorted by value; vals: sorted values."""
    base = 1 << (n - 1)
    N = order.size
    histA = np.zeros((n + 1, n + 1), dtype=np.int64)   # unordered pairs, all words
    histW = np.zeros((n + 1, n + 1), dtype=np.int64)   # unordered pairs, both in W
    viol = np.zeros(8, dtype=np.int64)
    stats = np.zeros(8, dtype=np.int64)   # 0: pairs, 1: max |s-t'| - (rho+2), 2: max |s-t'| - (R'+3), 3: max |gamma|/2^w *1000, 4: max group
    stats[1] = -10**9; stats[2] = -10**9
    nk = 0
    i = 0
    while i < N:
        j = i + 1
        while j < N and vals[j] == vals[i]:
            j += 1
        m = j - i
        if m > stats[4]:
            stats[4] = m
        if m >= 2:
            for a in range(i, j):
                h = base + order[a]
                for b in range(a + 1, j):
                    hp = base + order[b]
                    d = h ^ hp
                    tp = n - bitlen(d)
                    s = 0
                    while ((d >> s) & 1) == 0:
                        s += 1
                    w = n - tp - s
                    x = (h >> s) & ((1 << w) - 1)
                    xp = (hp >> s) & ((1 << w) - 1)
                    Z = h & ((1 << s) - 1)
                    Pp = h >> (w + s)
                    X = x - xp
                    Xi = revw(x, w) - revw(xp, w)
                    stats[0] += 1
                    histA[tp, s] += 1
                    if tp < 1 or w < 1:
                        viol[0] += 1
                    if (X & 1) == 0 or (Xi & 1) == 0:
                        viol[1] += 1
                    c = (1 << (s + 1)) - revw(Z, s)
                    gam = 2 * (q[x] - q[xp]) + (Pp & 1) * (((x >> (w - 1)) & 1) - ((xp >> (w - 1)) & 1)) \
                          + (c & 1) * ((x & 1) - (xp & 1))
                    if c * X - Pp * Xi != gam:
                        viol[2] += 1
                    if abs(gam) > (1 << w):
                        viol[3] += 1
                    v1000 = (abs(gam) * 1000) >> w
                    if v1000 > stats[3]:
                        stats[3] = v1000
                    aX = abs(X); aXi = abs(Xi)
                    mn = aX if aX < aXi else aXi
                    rho = w - (bitlen(mn) - 1)
                    dd = abs(s - tp)
                    if dd - (rho + 2) > stats[1]:
                        stats[1] = dd - (rho + 2)
                    if dd > rho + 2:
                        viol[4] += 1
                    r1 = maxrun(h, n); r2 = maxrun(hp, n)
                    Rp = r1 if r1 > r2 else r2
                    if dd - (Rp + 3) > stats[2]:
                        stats[2] = dd - (Rp + 3)
                    if dd > Rp + 3:
                        viol[5] += 1
                    # |X|,|Xi| >= 2^(w-1-R') claimed for words without runs longer than R'
                    if w - 1 - Rp >= 0:
                        if mn < (1 << (w - 1 - Rp)):
                            viol[6] += 1
                    if Rp <= R:
                        histW[tp, s] += 1
                    if want_keys:
                        # ordered window pair key (both orders)
                        keys[nk] = (((tp << 6) | s) << 52) | (x << 26) | xp
                        keys[nk + 1] = (((tp << 6) | s) << 52) | (xp << 26) | x
                        nk += 2
        i = j
    return histA, histW, viol, stats, nk

@njit(cache=True)
def check_keys(keys, nk, n):
    """keys sorted; for each distinct (t', s, x, x') count solutions and compare with the lattice-line bound."""
    worst_num = 0; worst_den = 1   # max of count / (4 (1 + 2^max(t'-1,s)/M))  as a fraction
    viol = 0
    worst_cnt = 0
    i = 0
    nd = 0
    while i < nk:
        j = i + 1
        while j < nk and keys[j] == keys[i]:
            j += 1
        cnt = j - i
        k = keys[i]
        xp = k & ((1 << 26) - 1)
        x = (k >> 26) & ((1 << 26) - 1)
        ts = k >> 52
        s = ts & 63
        tp = ts >> 6
        w = n - tp - s
        X = x - xp
        Xi = revw(x, w) - revw(xp, w)
        aX = abs(X); aXi = abs(Xi)
        g = gcd(aX, aXi)
        M = (aX if aX > aXi else aXi) // g
        e = tp - 1 if tp - 1 > s else s
        # bound = 4 (1 + 2^e / M) = 4 (M + 2^e) / M ; cnt <= bound  <=>  cnt * M <= 4 (M + 2^e)
        if cnt * M > 4 * (M + (1 << e)):
            viol += 1
        # track the max ratio cnt*M / (4 (M + 2^e))
        num = cnt * M; den = 4 * (M + (1 << e))
        if num * worst_den > worst_num * den:
            worst_num = num; worst_den = den
        if cnt > worst_cnt:
            worst_cnt = cnt
        nd += 1
        i = j
    return viol, worst_num, worst_den, worst_cnt, nd

def S_exact(w, cache={}):
    """returns (#odd pairs, S_w as float from exact counts)."""
    if w in cache:
        return cache[w]
    fn = "cnt_w%d.npy" % w
    if w >= 3 and os.path.exists(fn):
        cnt = np.load(fn)
        nz = np.nonzero(cnt)[0]
        res = (int(cnt.sum()), math.fsum(int(cnt[M]) / int(M) for M in nz))
    else:
        assert w <= 10, "need cnt file for w=%d" % w
        n = 1 << w
        rv = [int(format(x, "0%db" % w)[::-1], 2) for x in range(n)]
        nodd = 0; S = 0.0
        for x in range(n):
            for xp in range(n):
                X = x - xp; Xi = rv[x] - rv[xp]
                if X % 2 and Xi % 2:
                    nodd += 1
                    S += math.gcd(X, Xi) / max(abs(X), abs(Xi))
        res = (nodd, S)
    cache[w] = res
    return res

def main():
    n0 = int(sys.argv[1]); n1 = int(sys.argv[2])
    outname = sys.argv[3]
    nkey = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    out = []
    def P(*a):
        s = " ".join(str(t) for t in a)
        print(s, flush=True)
        out.append(s)
    for n in range(n0, n1 + 1):
        t0 = time.time()
        Ntot = 1 << n
        q = np.zeros(Ntot, dtype=np.int64)
        for j in range(1, n):
            idx = np.arange(1 << j, 1 << (j + 1), dtype=np.int64)
            hh = idx >> 1
            q[idx] = q[hh] + (idx & 1) * ((hh + 1) >> 1)
        q[1] = 0
        hs = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
        xv = 2 * hs + 1 - 2 * q[hs]
        order = np.argsort(xv, kind="stable").astype(np.int64)
        vals = xv[order]
        R = math.ceil(math.log2(n)) + 1
        t = -(-n // 4)
        want_keys = n <= nkey
        # number of unordered pairs: sum m(m-1)/2
        uniq, mult = np.unique(vals, return_counts=True)
        npairs = int(np.sum(mult.astype(np.int64) * (mult.astype(np.int64) - 1) // 2))
        keys = np.zeros(2 * npairs if want_keys else 1, dtype=np.int64)
        histA, histW, viol, stats, nk = scan_pairs(n, order, vals, q, R, keys, want_keys)
        P("n=%d: words 2^(n-1)=%d, distinct values %d, unordered colliding pairs %d (energy incl. diagonal / words = %.4f), "
          "largest fibre %d" % (n, 1 << (n - 1), uniq.size, npairs, (2 * npairs + (1 << (n - 1))) / (1 << (n - 1)), int(stats[4])))
        P("   violations: [t'<1 or w<1]=%d [X or Xi even]=%d [collision equation]=%d [|gamma|>2^w]=%d "
          "[|s-t'|>rho+2]=%d [|s-t'|>R'+3]=%d [min(|X|,|Xi|)<2^(w-1-R')]=%d"
          % tuple(int(v) for v in viol[:7]))
        P("   max(|s-t'| - (rho+2)) = %d ; max(|s-t'| - (R'+3)) = %d ; max |gamma|/2^w = %.3f"
          % (int(stats[1]), int(stats[2]), int(stats[3]) / 1000.0))
        # W
        # count W directly
        Wsize = count_W(n, R)
        P("   R=%d t=%d: |W|=%d  2^(n-2)=%d  |W|>=2^(n-2): %s ; words with a run >= R+1: %d <= (n-R)2^(n-1-R)=%d: %s"
          % (R, t, Wsize, 1 << (n - 2), Wsize >= (1 << (n - 2)), (1 << (n - 1)) - Wsize, (n - R) * (1 << (n - 1 - R)),
             (1 << (n - 1)) - Wsize <= (n - R) * (1 << (n - 1 - R))))
        # (D) per (t', s) for ALL words: ordered count <= 4 #odd + 4 2^max(t'-1,s) S_w
        worst = 0.0; wargs = None
        badcells = 0
        for tp in range(1, n):
            for s in range(0, n - tp):
                w = n - tp - s
                cntord = 2 * int(histA[tp, s])
                if cntord == 0:
                    continue
                if w > 20:
                    continue
                nodd, S = S_exact(w)
                bound = 4 * nodd + 4 * 2 ** max(tp - 1, s) * S
                r = cntord / bound
                if r > worst:
                    worst = r; wargs = (tp, s, w, cntord, bound)
                if cntord > bound:
                    badcells += 1
        P("   (D) per-cell: ordered pairs(t',s) <= 4*#odd(w) + 4*2^max(t'-1,s) S_w for all cells with w<=20: violations %d; "
          "max ratio %.5f at (t',s,w,count,bound)=%s" % (badcells, worst, str(wargs)))
        # E for W and prefix >= t
        E = 2 * int(histW[t:, :].sum())
        Eall = 2 * int(histA[t:, :].sum())
        # check the balance claim on W cells: all nonzero cells of histW have |s-t'| <= R+3
        offbal = 0
        maxdd = 0
        for tp in range(1, n):
            for s in range(0, n - tp):
                if histW[tp, s] > 0:
                    maxdd = max(maxdd, abs(s - tp))
                    if abs(s - tp) > R + 3:
                        offbal += 1
        boundD = 0.0
        ncell = 0
        for tp in range(t, n):
            for s in range(max(0, tp - R - 3), min(n - tp - 1, tp + R + 3) + 1):
                w = n - tp - s
                if w < 1:
                    continue
                ncell += 1
                nodd, S = S_exact(w)
                boundD += 4 * nodd + 4 * 2 ** max(tp - 1, s) * S
        P("   W-pairs: max |s-t'| = %d (R+3 = %d), cells off balance: %d" % (maxdd, R + 3, offbal))
        P("   (D): E = ordered pairs in W with common prefix >= t: %d (all words: %d);  E/2^n = %.6f ; explicit (D) sum over %d cells "
          "(<= n(2R+7) = %d) = %.4e ; E <= sum: %s ; sum/2^n = %.2f ; n^(9/2) log n = %.3e"
          % (E, Eall, E / 2 ** n, ncell, n * (2 * R + 7), boundD, E <= boundD, boundD / 2 ** n, n ** 4.5 * math.log(n)))
        # (E): distinct values on W and the Cauchy-Schwarz chain
        dW, sum_dp, lower_chain = distinct_on_W(n, R, t, xv)
        lb = Wsize ** 2 / (2 ** (t - 1) * (Wsize + E))
        P("   (E): distinct values on W = %d ; sum_p d_p = %d ; |W|^2/(|W|+E) = %.1f <= sum_p d_p: %s ; "
          "lower bound |W|^2/(2^(t-1)(|W|+E)) = %.1f <= distinct: %s ; distinct/2^(3n/4) = %.3f ; max value %d <= 2^(n+1)-1: %s"
          % (dW, sum_dp, Wsize ** 2 / (Wsize + E), Wsize ** 2 / (Wsize + E) <= sum_dp, lb, lb <= dW,
             dW / 2 ** (0.75 * n), int(vals[-1]), int(vals[-1]) <= 2 ** (n + 1) - 1))
        # prefix profile
        prof = []
        tot_by_tp = histA.sum(axis=1)
        for tp in range(1, min(n, 14)):
            prof.append("%d:%.4f" % (tp, 2 * int(tot_by_tp[tp:].sum()) / (1 << (n - 1))))
        P("   ordered colliding pairs per word with common prefix >= t': " + " ".join(prof))
        if want_keys:
            keys = keys[:nk]
            keys.sort()
            v, wn, wd, wc, nd = check_keys(keys, nk, n)
            P("   lattice line (C): distinct ordered (t',s,x,x') = %d ; violations of #(P,c) <= 4(1 + 2^max(t'-1,s)/M): %d ; "
              "max ratio %.4f ; max #(P,c) for one window pair = %d" % (nd, v, wn / wd, wc))
        P("   time %.1fs" % (time.time() - t0))
        del q, xv, order, vals, keys
    with open(outname, "w") as f:
        f.write("\n".join(out) + "\n")

@njit(parallel=True, cache=True)
def _count_W(n, R):
    base = 1 << (n - 1)
    tot = 0
    for h in prange(base, 2 * base):
        if maxrun(h, n) <= R:
            tot += 1
    return tot

def count_W(n, R):
    return int(_count_W(n, R))

@njit(cache=True)
def _inW(n, R, flags):
    base = 1 << (n - 1)
    for i in range(base):
        flags[i] = maxrun(base + i, n) <= R

def distinct_on_W(n, R, t, xv):
    base = 1 << (n - 1)
    flags = np.zeros(base, dtype=np.bool_)
    _inW(n, R, flags)
    vW = xv[flags]
    hW = np.arange(base, 2 * base, dtype=np.int64)[flags]
    dW = np.unique(vW).size
    pref = hW >> (n - t)
    key = pref * (1 << (n + 2)) + vW
    sum_dp = np.unique(key).size
    return dW, sum_dp, None

if __name__ == "__main__":
    main()
