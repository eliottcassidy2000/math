#!/usr/bin/env python3
"""Independent check (verifier C) of THM-5141 (C) on ALL colliding pairs of n-bit words, n = 1 .. NMAX.

usage: python c2_collisions.py NMAX [NPY]
  NMAX : largest word length (default 20)
  NPY  : pure-Python string-based cross-check of every per-pair quantity for n <= NPY (default 12)

For every ORDERED pair h != h' of n-bit words (top bit 1) with x(h) = x(h') (x(h) = 2h + 1 - 2q(h)):
  t' = length of the longest common prefix, s = length of the longest common suffix, w = n - t' - s,
  h = P | x | Z, h' = P | x' | Z, c = 2^(s+1) - rev_s(Z), X = x - x', Xi = rev_w x - rev_w x',
  gamma = 2(q(x) - q(x')) + (P mod 2)(x_(w-1) - x'_(w-1)) + (c mod 2)(x_0 - x'_0).
Checks (violation counters, all must be 0):
  v0  common prefix / suffix decomposition consistent          v1  t' >= 1, s >= 0, w >= 1
  v2  P in [2^(t'-1), 2^t')                                      v3  c in (2^s, 2^(s+1)]
  v4  X and Xi odd                                               v5  c X - P Xi = gamma
  v6  |gamma| <= 2^w                                             v7  (information) X Xi > 0
  v8  |Xi| P <= |X| c + |gamma| < 2^(w+s+2)                      v9  |Xi| P >= 2^(w-rho+t'-1)
  v10 |X| c <= |Xi| P + |gamma| < 2^(w+t')                       v11 |X| c > 2^(w-rho+s)
  v12 t' <= s + rho + 2                                          v13 s <= t' + rho - 1
  v14 (information, THM-4890 (C) translated) t' - s - (dXi - dX) in {0, 1, 2}
  v15 run lemma: |X|, |Xi| >= 2^(w-1-m), strict if w >= m + 2   (m = longest run of equal bits in h or h')
  v16 |s - t'| <= m + 3  (this is the statement for every R >= m, i.e. for every R for which the pair is run-free)
where rho = w - floor(log2 min(|X|, |Xi|)) is the least integer with 2^(w-rho) <= |X|, |Xi|.
Lattice-line bound: the ordered pairs are grouped by (t', s, x, x'); the size of a group is exactly the number of
(P, c) of statement (C) for that odd pair, and is compared with 4 (1 + 2^max(t'-1,s) gcd/max(|X|,|Xi|)); the groups
refined by (P mod 2, c mod 2) are compared with 1 + g min(2^(t'-1)/|X|, 2^s/|Xi|).
"""
import sys
import time
import numpy as np
from numba import njit

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
NPY = int(sys.argv[2]) if len(sys.argv) > 2 else 12
NV = 24
NE = 24


@njit(cache=False)
def build_q(nmax):
    L = 1 << nmax
    q = np.zeros(L, dtype=np.int64)
    for h in range(1, L):
        if h & 1:
            q[h] = q[h >> 1] + (((h >> 1) + 1) >> 1)
        else:
            q[h] = q[h >> 1]
    return q


@njit(cache=False)
def bitlen(a):
    b = 0
    while a > 0:
        a >>= 1
        b += 1
    return b


@njit(cache=False)
def revw(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


@njit(cache=False)
def maxrun(h, n):
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


@njit(cache=False)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=False)
def collect(n, q):
    lo = 1 << (n - 1)
    m = lo
    vals = np.empty(m, dtype=np.int64)
    for i in range(m):
        h = lo + i
        vals[i] = 2 * h + 1 - 2 * q[h]
    order = np.argsort(vals)
    total = 0
    i = 0
    maxgroup = 1
    ndist = 0
    while i < m:
        j = i + 1
        while j < m and vals[order[j]] == vals[order[i]]:
            j += 1
        L = j - i
        total += L * (L - 1)
        if L > maxgroup:
            maxgroup = L
        ndist += 1
        i = j
    A = np.empty(total, dtype=np.int64)
    B = np.empty(total, dtype=np.int64)
    k = 0
    i = 0
    while i < m:
        j = i + 1
        while j < m and vals[order[j]] == vals[order[i]]:
            j += 1
        for a in range(i, j):
            for b in range(i, j):
                if a != b:
                    A[k] = lo + order[a]
                    B[k] = lo + order[b]
                    k += 1
        i = j
    return A, B, ndist, maxgroup


@njit(cache=False)
def analyse(n, A, B, q, viol, ext, exw, hist, hist90, fext, dump):
    K = A.shape[0]
    key = np.empty(K, dtype=np.int64)
    key2 = np.empty(K, dtype=np.int64)
    dodump = dump.shape[0] == K
    for k in range(K):
        h = A[k]
        hp = B[k]
        d = h ^ hp
        hi = bitlen(d) - 1
        lo = 0
        while ((d >> lo) & 1) == 0:
            lo += 1
        tp = n - 1 - hi
        s = lo
        w = hi - lo + 1
        P = h >> (hi + 1)
        Z = h & ((1 << s) - 1)
        x = (h >> s) & ((1 << w) - 1)
        xp = (hp >> s) & ((1 << w) - 1)
        if (hp >> (hi + 1)) != P or (hp & ((1 << s) - 1)) != Z:
            viol[0] += 1
        if ((P << (w + s)) | (x << s) | Z) != h or ((P << (w + s)) | (xp << s) | Z) != hp:
            viol[0] += 1
        if tp < 1 or s < 0 or w < 1 or tp + s + w != n:
            viol[1] += 1
        if not ((1 << (tp - 1)) <= P and P < (1 << tp)):
            viol[2] += 1
        c = (1 << (s + 1)) - revw(Z, s)
        if not ((1 << s) < c and c <= (1 << (s + 1))):
            viol[3] += 1
        X = x - xp
        Xi = revw(x, w) - revw(xp, w)
        if (X & 1) == 0 or (Xi & 1) == 0:
            viol[4] += 1
        dt = ((x >> (w - 1)) & 1) - ((xp >> (w - 1)) & 1)
        db = (x & 1) - (xp & 1)
        gam = 2 * (q[x] - q[xp]) + (P & 1) * dt + (c & 1) * db
        if c * X - P * Xi != gam:
            viol[5] += 1
        ag = abs(gam)
        if ag > (1 << w):
            viol[6] += 1
        if X * Xi <= 0:
            viol[7] += 1
        aX = abs(X)
        aXi = abs(Xi)
        mn = min(aX, aXi)
        rho = w - (bitlen(mn) - 1)
        # carries, line by line
        if not (aXi * P <= aX * c + ag and aX * c + ag < (1 << (w + s + 2))):
            viol[8] += 1
        if not ((aXi * P) << rho >= (1 << (w + tp - 1))):
            viol[9] += 1
        if not (aX * c <= aXi * P + ag and aXi * P + ag < (1 << (w + tp))):
            viol[10] += 1
        if not ((aX * c) << rho > (1 << (w + s))):
            viol[11] += 1
        if not (tp <= s + rho + 2):
            viol[12] += 1
        if not (s <= tp + rho - 1):
            viol[13] += 1
        dX = w - bitlen(aX)
        dXi = w - bitlen(aXi)
        val = tp - s - (dXi - dX)
        if val < 0 or val > 2:
            viol[14] += 1
        else:
            hist90[val] += 1
        m1 = maxrun(h, n)
        m2 = maxrun(hp, n)
        m = max(m1, m2)
        if w - 1 - m >= 0:
            lowb = 1 << (w - 1 - m)
            if aX < lowb or aXi < lowb:
                viol[15] += 1
            if w >= m + 2 and (aX <= lowb or aXi <= lowb):
                viol[15] += 1
        ad = abs(s - tp)
        if ad > m + 3:
            viol[16] += 1
        # extremes
        sc = ag << (30 - w)                     # |gamma| / 2^w scaled by 2^30
        if sc > ext[0]:
            ext[0] = sc
            exw[0, 0] = h
            exw[0, 1] = hp
        if tp - s - rho > ext[1]:
            ext[1] = tp - s - rho
            exw[1, 0] = h
            exw[1, 1] = hp
        if s - tp - rho > ext[2]:
            ext[2] = s - tp - rho
            exw[2, 0] = h
            exw[2, 1] = hp
        if ad - m > ext[3]:
            ext[3] = ad - m
            exw[3, 0] = h
            exw[3, 1] = hp
        if s - tp - m > ext[4]:
            ext[4] = s - tp - m
            exw[4, 0] = h
            exw[4, 1] = hp
        if tp - s - m > ext[5]:
            ext[5] = tp - s - m
            exw[5, 0] = h
            exw[5, 1] = hp
        if ad > ext[6]:
            ext[6] = ad
            exw[6, 0] = h
            exw[6, 1] = hp
        if rho > ext[7]:
            ext[7] = rho
            exw[7, 0] = h
            exw[7, 1] = hp
        if w > ext[8]:
            ext[8] = w
        # float tightness of the two strict upper inequalities
        r1 = (aX * c + ag) / float(1 << (w + s + 2))
        if r1 > fext[0]:
            fext[0] = r1
        r3 = (aXi * P + ag) / float(1 << (w + tp))
        if r3 > fext[1]:
            fext[1] = r3
        r2 = float(1 << (w + tp - 1)) / float((aXi * P) << rho)       # (ii): <= 1
        if r2 > fext[2]:
            fext[2] = r2
        r4 = float(1 << (w + s)) / float((aX * c) << rho)            # (iv): < 1
        if r4 > fext[3]:
            fext[3] = r4
        if m < hist.shape[0] and ad < hist.shape[1]:
            hist[m, ad] += 1
        key[k] = (((tp << 5) | s) << 44) | (x << 22) | xp
        key2[k] = (key[k] << 2) | ((P & 1) << 1) | (c & 1)
        if dodump:
            dump[k, 0] = tp
            dump[k, 1] = s
            dump[k, 2] = w
            dump[k, 3] = P
            dump[k, 4] = c
            dump[k, 5] = x
            dump[k, 6] = xp
            dump[k, 7] = X
            dump[k, 8] = Xi
            dump[k, 9] = gam
            dump[k, 10] = rho
            dump[k, 11] = m
    return key, key2


@njit(cache=False)
def lattice(n, key, key2, lv, lf, lex):
    """lv[0] groups, lv[1] violations of the stated bound, lv[2] parity-class groups, lv[3] violations of
    1 + g min(2^(t'-1)/|X|, 2^s/|Xi|), lv[4] violations of the sharper 1 + floor((2^(t'-1)-1) g/(2|X|)) and
    1 + floor((2^s - 1) g / (2|Xi|)), lv[5] largest group, lv[6] largest parity-class group, lv[7] groups with
    solutions in all 4 parity classes, lv[8] sum of group sizes; lf[0] max count/bound, lf[1] max class count/bound"""
    ks = np.sort(key)
    K = ks.shape[0]
    i = 0
    while i < K:
        j = i + 1
        while j < K and ks[j] == ks[i]:
            j += 1
        cnt = j - i
        kk = ks[i]
        xp = kk & ((1 << 22) - 1)
        x = (kk >> 22) & ((1 << 22) - 1)
        ts = kk >> 44
        s = ts & 31
        tp = ts >> 5
        w = n - tp - s
        X = x - xp
        Xi = revw(x, w) - revw(xp, w)
        aX = abs(X)
        aXi = abs(Xi)
        g = gcd(aX, aXi)
        M = max(aX, aXi)
        e = max(tp - 1, s)
        lv[0] += 1
        lv[8] += cnt
        # cnt <= 4 (1 + 2^e g / M)
        if cnt * M > 4 * (M + (1 << e) * g):
            lv[1] += 1
        ratio = cnt * M / float(4 * (M + (1 << e) * g))
        if ratio > lf[0]:
            lf[0] = ratio
            lex[0, 0] = tp
            lex[0, 1] = s
            lex[0, 2] = x
            lex[0, 3] = xp
            lex[0, 4] = cnt
        if cnt > lv[5]:
            lv[5] = cnt
            lex[1, 0] = tp
            lex[1, 1] = s
            lex[1, 2] = x
            lex[1, 3] = xp
            lex[1, 4] = cnt
        i = j
    ks2 = np.sort(key2)
    i = 0
    ncls_in_group = 0
    prevkey = -1
    while i < K:
        j = i + 1
        while j < K and ks2[j] == ks2[i]:
            j += 1
        cnt = j - i
        kk = ks2[i] >> 2
        if kk == prevkey:
            ncls_in_group += 1
        else:
            ncls_in_group = 1
            prevkey = kk
        if ncls_in_group == 4:
            lv[7] += 1
        xp = kk & ((1 << 22) - 1)
        x = (kk >> 22) & ((1 << 22) - 1)
        ts = kk >> 44
        s = ts & 31
        tp = ts >> 5
        w = n - tp - s
        X = x - xp
        Xi = revw(x, w) - revw(xp, w)
        aX = abs(X)
        aXi = abs(Xi)
        g = gcd(aX, aXi)
        lv[2] += 1
        # cnt <= 1 + g min(2^(t'-1)/|X|, 2^s/|Xi|)
        if (cnt - 1) * aX > g * (1 << (tp - 1)) or (cnt - 1) * aXi > g * (1 << s):
            lv[3] += 1
        b1 = 1 + ((1 << (tp - 1)) - 1) * g // (2 * aX)
        b2 = 1 + ((1 << s) - 1) * g // (2 * aXi)
        if cnt > b1 or cnt > b2:
            lv[4] += 1
        bb = 1.0 + g * min((1 << (tp - 1)) / float(aX), (1 << s) / float(aXi))
        if cnt / bb > lf[1]:
            lf[1] = cnt / bb
        if cnt > lv[6]:
            lv[6] = cnt
            lex[2, 0] = tp
            lex[2, 1] = s
            lex[2, 2] = x
            lex[2, 3] = xp
            lex[2, 4] = cnt
        i = j


# ------------------------------------------------------------------ pure Python reference (strings)
def py_q(h):
    t = 0
    while h:
        if h & 1:
            t += ((h >> 1) + 1) // 2
        h >>= 1
    return t


def py_pairs(n):
    from collections import defaultdict
    byval = defaultdict(list)
    for h in range(1 << (n - 1), 1 << n):
        byval[2 * h + 1 - 2 * py_q(h)].append(h)
    res = {}
    for v, hs in byval.items():
        for h in hs:
            for hp in hs:
                if h == hp:
                    continue
                a = format(h, "0%db" % n)
                b = format(hp, "0%db" % n)
                tp = 0
                while a[tp] == b[tp]:
                    tp += 1
                s = 0
                while a[n - 1 - s] == b[n - 1 - s]:
                    s += 1
                w = n - tp - s
                Pstr, xs, Zs = a[:tp], a[tp:tp + w], a[tp + w:]
                xps = b[tp:tp + w]
                assert b[:tp] == Pstr and b[tp + w:] == Zs
                P = int(Pstr, 2)
                x = int(xs, 2)
                xp = int(xps, 2)
                revZ = int(Zs[::-1], 2) if s else 0
                c = 2 ** (s + 1) - revZ
                X = x - xp
                Xi = int(xs[::-1], 2) - int(xps[::-1], 2)
                gam = 2 * (py_q(x) - py_q(xp)) + (P % 2) * (int(xs[0]) - int(xps[0])) + (c % 2) * (int(xs[-1]) - int(xps[-1]))
                mn = min(abs(X), abs(Xi))
                rho = w - (mn.bit_length() - 1)
                import itertools
                m = max(max(len(list(grp)) for _, grp in itertools.groupby(a)),
                        max(len(list(grp)) for _, grp in itertools.groupby(b)))
                res[(h, hp)] = (tp, s, w, P, c, x, xp, X, Xi, gam, rho, m)
    return res


def bits(h, n):
    return format(int(h), "0%db" % n)


def main():
    t0 = time.time()
    q = build_q(max(NMAX, 2))
    print("THM-5141 (C): all colliding pairs, n = 1 ..", NMAX, flush=True)
    grand = 0
    cum_unordered = 0
    allviol = np.zeros(NV, dtype=np.int64)
    gext = np.full(NE, -10 ** 9, dtype=np.int64)
    gfext = np.zeros(8, dtype=np.float64)
    ghist = np.zeros((40, 40), dtype=np.int64)
    ghist90 = np.zeros(3, dtype=np.int64)
    glv = np.zeros(12, dtype=np.int64)
    for n in range(1, NMAX + 1):
        A, B, ndist, maxgroup = collect(n, q)
        K = A.shape[0]
        viol = np.zeros(NV, dtype=np.int64)
        ext = np.full(NE, -10 ** 9, dtype=np.int64)
        exw = np.zeros((NE, 2), dtype=np.int64)
        hist = np.zeros((40, 40), dtype=np.int64)
        hist90 = np.zeros(3, dtype=np.int64)
        fext = np.zeros(8, dtype=np.float64)
        dump = np.zeros((K if n <= NPY else 0, 12), dtype=np.int64)
        key, key2 = analyse(n, A, B, q, viol, ext, exw, hist, hist90, fext, dump)
        lv = np.zeros(12, dtype=np.int64)
        lf = np.zeros(4, dtype=np.float64)
        lex = np.zeros((4, 5), dtype=np.int64)
        if K:
            lattice(n, key, key2, lv, lf, lex)
        pychk = ""
        if n <= NPY:
            ref = py_pairs(n)
            ok = (len(ref) == K)
            for k in range(K):
                tup = tuple(int(v) for v in dump[k])
                if ref.get((int(A[k]), int(B[k]))) != tup:
                    ok = False
            pychk = " | pure-Python string reference agrees on every pair and every quantity: %s" % ok
        grand += K
        cum_unordered += K // 2
        allviol += viol
        ghist += hist
        ghist90 += hist90
        for i in range(NE):
            gext[i] = max(gext[i], ext[i])
        for i in range(8):
            gfext[i] = max(gfext[i], fext[i])
        glv[:5] += lv[:5]
        glv[7] += lv[7]
        glv[8] += lv[8]
        glv[5] = max(glv[5], lv[5])
        glv[6] = max(glv[6], lv[6])
        print("n = %2d: distinct values %8d of %8d words (largest fibre %d); colliding pairs: %9d unordered = %9d ordered"
              % (n, ndist, 1 << (n - 1), maxgroup, K // 2, K) + pychk, flush=True)
        if K == 0:
            continue
        print("        violations v0..v16:", " ".join(str(int(v)) for v in viol[:17]),
              " (v7 and v14 are information only)", flush=True)
        print("        max |gamma|/2^w = %.6f; max (t' - s - rho) = %d (lemma: <= 2); max (s - t' - rho) = %d (lemma: <= -1);"
              % (ext[0] / float(1 << 30), ext[1], ext[2]),
              "max rho = %d; max w = %d" % (ext[7], ext[8]), flush=True)
        print("        run-free: max (|s - t'| - m) = %d (statement: <= 3); max (s - t' - m) = %d; max (t' - s - m) = %d; max |s - t'| = %d"
              % (ext[3], ext[4], ext[5], ext[6]), flush=True)
        print("        tightness: max (|X| c + |gamma|)/2^(w+s+2) = %.6f (< 1); max (|Xi| P + |gamma|)/2^(w+t') = %.6f (< 1);"
              % (fext[0], fext[1]),
              "max 2^(w-rho+t'-1)/(|Xi| P) = %.6f (<= 1); max 2^(w-rho+s)/(|X| c) = %.6f (< 1)" % (fext[2], fext[3]), flush=True)
        print("        THM-4890 (C) translated, t' - s - (dXi - dX) = 0, 1, 2 :", int(hist90[0]), int(hist90[1]), int(hist90[2]), flush=True)
        print("        lattice line: groups (t', s, x, x') with a solution: %d (sum of sizes %d); violations of 4(1 + 2^max(t'-1,s) g/max): %d;"
              % (lv[0], lv[8], lv[1]),
              "max count/bound = %.4f; largest group %d;" % (lf[0], lv[5]),
              "parity classes with a solution: %d; violations of 1 + g min(...): %d; of the halved-step bound: %d;"
              % (lv[2], lv[3], lv[4]),
              "max class count/bound = %.4f; largest class %d; groups with all 4 parity classes: %d" % (lf[1], lv[6], lv[7]), flush=True)
        w_ex = n - lex[0, 0] - lex[0, 1]
        print("        extremal group for count/bound: t' = %d, s = %d, x = %s, x' = %s, count %d"
              % (lex[0, 0], lex[0, 1], bits(lex[0, 2], w_ex), bits(lex[0, 3], w_ex), lex[0, 4]), flush=True)
        for (i, name) in ((1, "t' - s - rho"), (2, "s - t' - rho"), (3, "|s - t'| - m"), (0, "|gamma|/2^w")):
            print("        example for max %-13s: h = %s, h' = %s" % (name, bits(exw[i, 0], n), bits(exw[i, 1], n)), flush=True)
    print()
    print("TOTAL n <= %d: ordered colliding pairs %d (unordered %d)" % (NMAX, grand, cum_unordered))
    print("TOTAL violations v0..v16:", " ".join(str(int(v)) for v in allviol[:17]))
    print("TOTAL extremes: max |gamma|/2^w = %.6f; max (t' - s - rho) = %d; max (s - t' - rho) = %d; max (|s-t'| - m) = %d;"
          % (gext[0] / float(1 << 30), gext[1], gext[2], gext[3]),
          "max (s - t' - m) = %d; max (t' - s - m) = %d; max |s - t'| = %d; max rho = %d" % (gext[4], gext[5], gext[6], gext[7]))
    print("TOTAL tightness: %.6f %.6f %.6f %.6f" % (gfext[0], gfext[1], gfext[2], gfext[3]))
    print("TOTAL THM-4890 translated histogram (0, 1, 2):", [int(v) for v in ghist90])
    print("TOTAL lattice: groups %d, violations %d; classes %d, violations %d, halved-step violations %d; largest group %d,"
          % (glv[0], glv[1], glv[2], glv[3], glv[4], glv[5]), "largest class %d; groups with all four parity classes %d" % (glv[6], glv[7]))
    print("histogram of ordered pairs by m = longest run in h or h' (rows) and |s - t'| (columns 0..):")
    mmax = max(i for i in range(40) if ghist[i].sum() > 0)
    amax = max(j for j in range(40) if ghist[:, j].sum() > 0)
    print("   m \\ |s-t'| " + " ".join("%9d" % j for j in range(amax + 1)))
    for i in range(1, mmax + 1):
        print("   %2d         " % i + " ".join("%9d" % int(ghist[i, j]) for j in range(amax + 1)))
    print("for every R: run-free pairs (m <= R) and their max |s - t'| (statement: <= R + 3):")
    for R in range(1, mmax + 1):
        sub = ghist[1:R + 1].sum(axis=0)
        tot = int(sub.sum())
        mx = max((j for j in range(40) if sub[j] > 0), default=-1)
        print("   R = %2d: %10d ordered pairs, max |s - t'| = %2d, R + 3 = %2d  %s" % (R, tot, mx, R + 3, "OK" if mx <= R + 3 else "VIOLATION"))
    print("time %.1f s" % (time.time() - t0))


main()
