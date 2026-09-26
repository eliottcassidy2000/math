#!/usr/bin/env python3
"""
procgen_rank_20260926_q2.py -- Q2 of the rank lane: adaptive ("current center") ranks.

Candidates (a = 1, c = 1.1 log(3/2) unless stated):
  A_P   R = log m + c Lam_P(m),  Lam_P(m) = max v2(m - x) over expanding periodic x of period <= P
  A_inf the same with every period (Lam_inf); Proposition Q2.3: l*(m) <= Lam_inf(m) < 2 l*(m), where
        l*(m) = last time the multiplier 3^(a_l)/2^l of the orbit of m exceeds 1 (an excursion time)
  B_PQ  R = log m + c max_{0<=i<=Q} [i + Lam_P(T^i m)]   (preperiodic centers of preperiod <= Q)
  D_th  R = max_{0 <= i <= th log2 m} log T^i m          (windowed future envelope of log)
  C_th  R = max_{0 <= i <= th log2 m} R_A20(T^i m)       (windowed future envelope of A_20)
Test sets: hostile families; every step of the orbits of 3..5000; the trajectories of the path records
(OEIS A006884) and delay records (OEIS A006877) up to 10^12 (CITED lists, used only as inputs); 600 random
starts in [10^11, 10^12]; cycle controls for 3n-1 and 5n+1.
"""
import math
import random
from fractions import Fraction

import numpy as np

from procgen_rank_20260926_lib import (KAPPA, LN2, LN3, check, say, v2, T, periodic_point, parity_word,
                                       common_prefix, is_expanding_word)

# OEIS A006884 (path records) and A006877 (delay records), terms <= 10^12 (CITED; inputs only)
PATH_RECORDS = [1, 2, 3, 7, 15, 27, 255, 447, 639, 703, 1819, 4255, 4591, 9663, 20895, 26623, 31911, 60975,
                77671, 113383, 138367, 159487, 270271, 665215, 704511, 1042431, 1212415, 1441407, 1875711,
                1988859, 2643183, 2684647, 3041127, 3873535, 4637979, 5656191, 6416623, 6631675, 19638399,
                38595583, 80049391, 120080895, 210964383, 319804831, 1410123943, 8528817511, 12327829503,
                23035537407, 45871962271, 51739336447, 59152641055, 59436135663, 70141259775, 77566362559,
                110243094271, 204430613247, 231913730799, 272025660543, 446559217279, 567839862631,
                871673828443]
DELAY_RECORDS = [1, 2, 3, 6, 7, 9, 18, 25, 27, 54, 73, 97, 129, 171, 231, 313, 327, 649, 703, 871, 1161, 2223,
                 2463, 2919, 3711, 6171, 10971, 13255, 17647, 23529, 26623, 34239, 35655, 52527, 77031, 106239,
                 142587, 156159, 216367, 230631, 410011, 511935, 626331, 837799, 1117065, 1501353, 1723519,
                 2298025, 3064033, 3542887, 3732423, 5649499, 6649279, 8400511, 11200681, 14934241, 15733191,
                 31466382, 36791535, 63728127, 127456254, 169941673, 226588897, 268549803, 537099606, 670617279,
                 1341234558, 1412987847, 1674652263, 2610744987, 4578853915, 4890328815, 9780657630,
                 12212032815, 12235060455, 13371194527, 17828259369, 31694683323, 63389366646, 75128138247,
                 133561134663, 158294678119, 166763117679, 202485402111, 404970804222, 426635908975,
                 568847878633, 674190078379, 881715740415, 989345275647]

C_DEF = 1.1 * KAPPA


# ----------------------------------------------------------------------------------------------
# orbits, parity sequences, counters
# ----------------------------------------------------------------------------------------------
def orbit_with_tail(n, q=3, sgn=1, steps=None, tail=256):
    """values and parities of the orbit of n under x/2, (q x + sgn)/2; until 1 (+ tail) or `steps` steps"""
    vals = [n]
    x = n
    if steps is None:
        while x != 1:
            x = (q * x + sgn) >> 1 if x & 1 else x >> 1
            vals.append(x)
        for _ in range(tail):
            x = (q * x + sgn) >> 1 if x & 1 else x >> 1
            vals.append(x)
    else:
        for _ in range(steps):
            x = (q * x + sgn) >> 1 if x & 1 else x >> 1
            vals.append(x)
    return vals


def amin_table(Pmax, q=3):
    """least a with q^a > 2^p, for p = 0..Pmax"""
    out = []
    for p in range(Pmax + 1):
        a = 0
        while q ** a <= 2 ** p:
            a += 1
        out.append(a)
    return out


def lam_array(Q, P, q=3, amin=None):
    """Lam_P[t] = max over p <= P with Q[t:t+p] expanding of the length of the p-periodic run from t
    (= v2(T^t m - x) for the periodic point x of that block: parity vectors are a 2-adic isometry)"""
    Q = np.asarray(Q, dtype=np.int8)
    L = len(Q)
    if amin is None:
        amin = amin_table(P, q)
    cs = np.concatenate([[0], np.cumsum(Q, dtype=np.int64)])
    lam = np.zeros(L, dtype=np.int64)
    for p in range(1, min(P, L - 1) + 1):
        ones = cs[p:] - cs[:-p]                         # t = 0..L-p
        eq = Q[p:] == Q[:-p]                            # i = 0..L-p-1
        pos = np.arange(L - p)
        nf = np.where(eq, L - p, pos)
        nxt = np.minimum.accumulate(nf[::-1])[::-1]
        run = np.full(L - p + 1, p, dtype=np.int64)
        run[:L - p] = p + (nxt - pos)
        val = np.where(ones >= amin[p], run, 0)
        np.maximum(lam[:L - p + 1], val, out=lam[:L - p + 1])
    return lam


def lstar_array(Q, q=3):
    """l*(t) = max{l >= 1 : the block Q[t:t+l] is expanding} (0 if none), by a right-to-left scan"""
    Q = np.asarray(Q, dtype=np.int64)
    L = len(Q)
    g = np.where(Q == 1, math.log(q), 0.0) - LN2     # increments of log multiplier
    S = np.concatenate([[0.0], np.cumsum(g)])        # S[l] = log M at time l (from 0)
    out = np.zeros(L, dtype=np.int64)
    # l*(t) = max{l : S[t+l] - S[t] > 0}; scan positions of S in decreasing order of index keeping the max
    # over suffixes: for each t, the last index u > t with S[u] > S[t]
    # computed by sorting: process t from right to left with a monotone structure
    suffix_max = np.maximum.accumulate(S[::-1])[::-1]  # max_{u >= i} S[u]
    for t in range(L):
        if suffix_max[t + 1] <= S[t] + 1e-12:
            continue
        lo, hi = t + 1, L                              # find the last u with S[u] > S[t]: binary search on
        # suffix_max (non-increasing): last u with suffix_max[u] > S[t]
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if suffix_max[mid] > S[t] + 1e-12:
                lo = mid
            else:
                hi = mid - 1
        out[t] = lo - t
    return out


# ----------------------------------------------------------------------------------------------
# candidate increments along an orbit
# ----------------------------------------------------------------------------------------------
def increments(vals, R, tmax):
    return R[1:tmax + 1] - R[:tmax]


def orbit_stats(n, P_list=(1, 3, 20), c=C_DEF, thetas=(2, 4, 8, 16), Qwin=10):
    vals = orbit_with_tail(n, tail=64)
    Q = [v & 1 for v in vals]
    L = len(vals)
    # evaluate steps m -> T m for the points m >= 3 before the orbit first reaches 1
    tend = next(i for i, v in enumerate(vals) if v == 1)
    tmax = max([i for i in range(tend) if vals[i] >= 3], default=-1) + 1
    logs = np.log(np.array([float(v) for v in vals]))
    out = {}
    if tmax <= 0:
        return out, 0
    lam = {}
    for P in P_list:
        lam[P] = lam_array(Q, P)
        R = logs + c * lam[P]
        d = R[1:tmax + 1] - R[:tmax]
        dl = lam[P][1:tmax + 1] - lam[P][:tmax]
        out["A%d" % P] = (d, dl)
    li = lam_array(Q, L - 1)
    ls = lstar_array(Q)
    Rinf = logs + c * li
    out["Ainf"] = (Rinf[1:tmax + 1] - Rinf[:tmax], li[1:tmax + 1] - li[:tmax])
    out["lstar"] = (li[:tmax], ls[:tmax])
    out["big"] = (np.array([v >= 11 for v in vals[:tmax]]), None)
    # B_{3,Qwin}
    l3 = lam[3] if 3 in lam else lam_array(Q, 3)
    B = np.full(L, -10 ** 9, dtype=np.int64)
    for i in range(Qwin + 1):
        B[:L - i] = np.maximum(B[:L - i], i + l3[i:])
    RB = logs + c * B
    out["B3_%d" % Qwin] = (RB[1:tmax + 1] - RB[:tmax], None)
    # windowed envelopes
    RA20 = logs + c * (lam[20] if 20 in lam else lam_array(Q, 20))
    lg2 = logs / LN2
    for th in thetas:
        w = np.floor(th * lg2).astype(np.int64)
        RD = np.empty(tmax + 1)
        RC = np.empty(tmax + 1)
        for t in range(tmax + 1):
            e = min(L, t + w[t] + 1)
            RD[t] = logs[t:e].max()
            RC[t] = RA20[t:e].max()
        out["D%d" % th] = (RD[1:] - RD[:-1], None)
        out["C%d" % th] = (RC[1:] - RC[:-1], None)
    return out, tmax


def peak_time_ratio(vals, tmax):
    """max over t < tmax of (time from t to the global peak of the future orbit of vals[t]) / log2 vals[t]"""
    best = 0.0
    arr = np.array([float(v) for v in vals])
    # position of the last maximum of the suffix, computed right to left
    L = len(arr)
    arg = np.empty(L, dtype=np.int64)
    cur = L - 1
    for i in range(L - 1, -1, -1):
        if arr[i] >= arr[cur]:
            cur = i
        arg[i] = cur
    for t in range(tmax):
        r = (arg[t] - t) / math.log2(arr[t])
        best = max(best, r)
    return best


# ----------------------------------------------------------------------------------------------
# parts
# ----------------------------------------------------------------------------------------------
def part_isometry():
    say("== Q2.0 parity vectors are a 2-adic isometry: v2(m - x_w) = common prefix of Q(m) and w^inf ==")
    rng = random.Random(3)
    cnt = 0
    for _ in range(400):
        p = rng.randrange(1, 9)
        w = [rng.randrange(2) for _ in range(p)]
        if sum(w) == 0:
            continue
        x = periodic_point(w)
        for _ in range(5):
            m = rng.randrange(1, 10 ** 12)
            L = 80
            cp = common_prefix(parity_word(m, L), (w * (L // p + 2))[:L])
            val = v2(Fraction(m) - x)
            if cp < L and cp != val:
                check(False, "isometry fails: m=%d w=%s" % (m, w))
            cnt += 1
    check(cnt > 1500, "v2(m - x_w) equals the common-prefix length of the parity words in %d random cases" % cnt)


def lam_of_int(m, P, extra=64):
    vals = orbit_with_tail(m, steps=int(math.log2(m)) + 3 * P + extra)
    return lam_array([v & 1 for v in vals], P), vals


def part_Q21(c=C_DEF):
    say("== Q2.1 current-center ranks (period <= P) fail at the even preimages y_V = 2^(V+1) - 2 ==")
    for P in (1, 3, 11, 20):
        rows = []
        for V in (30, 60, 120, 240):
            y = 2 ** (V + 1) - 2
            lam, vals = lam_of_int(y, P)
            base = math.log(y) + c * lam[0]
            worst_min = min(math.log(vals[j]) + c * lam[j] - base for j in range(1, 11))
            rows.append((V, int(lam[0]), int(lam[1]), worst_min))
        lam_y = rows[0][1]
        ok = all(r[1] == (P if P >= 3 else 0) for r in rows) and all(r[2] >= r[0] for r in rows) \
            and all(r[3] >= c * (r[0] - 10 + 1 - P) - LN2 - 1e-9 for r in rows) and rows[-1][3] > 90
        check(ok, "P = %d: Lam_P(y_V) = %d = Lam_P(-2) for V = 30..240 while Lam_P(T y_V) >= V; "
              "min over 1 <= j <= 10 of R(T^j y_V) - R(y_V) = %s (>= c(V - 9 - P) - log 2): lookahead 10 fails" %
              (P, lam_y, ", ".join("%.1f" % r[3] for r in rows)))
    # the Kuratowski reset n_H = (2^(H+3) - 13)/9: Lam_P(n_H) bounded (n_H -> -13/9), Lam_P(T^3 n_H) >= H
    rows = []
    for j in (5, 10, 20, 40):
        H = 6 * j + 5
        n = (2 ** (H + 3) - 13) // 9
        lam, vals = lam_of_int(n, 20)
        rows.append((H, int(lam[0]), int(lam[3]), math.log(vals[3]) + c * lam[3] - math.log(n) - c * lam[0]))
    check(len(set(r[1] for r in rows)) == 1 and all(r[2] >= r[0] for r in rows) and rows[-1][3] > 90,
          "reset family n_H (H = 35..245): Lam_20(n_H) = %d constant, Lam_20(T^3 n_H) >= H, and "
          "R(T^3 n_H) - R(n_H) = %s" % (rows[0][1], ", ".join("%.1f" % r[3] for r in rows)))


def part_Q22(c=C_DEF):
    say("== Q2.2 bounded preperiod Q (period <= 3) fails on y = 2^(Q+1) (2^V - 1) ==")
    for Qw in (5, 10):
        incs = []
        for V in (40, 80, 160):
            y = 2 ** (Qw + 1) * (2 ** V - 1)
            vals = orbit_with_tail(y, steps=V + 3 * Qw + 80)
            lam = lam_array([v & 1 for v in vals], 3)
            Bv = [max(i + lam[t + i] for i in range(Qw + 1)) for t in range(2)]
            incs.append(math.log(vals[1]) + c * Bv[1] - math.log(vals[0]) - c * Bv[0])
        check(incs[0] > 0 and incs[1] > incs[0] and incs[2] > incs[1] and abs((incs[2] - incs[1]) - c * 80) < 1.0,
              "Q = %d: at y = 2^%d (2^V - 1) the first halving step increases R_B by %s at V = 40, 80, 160 "
              "(slope c: the preimage -2^%d of -1 lies just outside the window)" %
              (Qw, Qw + 1, ", ".join("%.2f" % v for v in incs), Qw + 1))


def lower_christoffel(a, p):
    return [(a * j) // p - (a * (j - 1)) // p for j in range(1, p + 1)]


def dlog2_mod_3pow(C, m):
    """x with 2^x = C (mod 3^m), 3 not dividing C (2 is a primitive root mod every 3^m); lifting mod phi(3^j)"""
    assert C % 3 != 0
    x = next(t for t in range(2) if pow(2, t, 3) == C % 3)
    for j in range(1, m):
        order = 2 * 3 ** (j - 1)
        x = next(x + t * order for t in range(3) if pow(2, x + t * order, 3 ** (j + 1)) == C % 3 ** (j + 1))
    return x


def hover_integer(w, r):
    """the positive integer k with parity word w^r followed by 0^s (T^(rp) k = 2^s), smallest s >= 1"""
    word = w * r
    A, B, C = 1, 1, 0
    for b in word:
        if b:
            A, C = 3 * A, 3 * C + B
        B *= 2
    rp = len(word)
    m = sum(word)
    x = dlog2_mod_3pow(C, m)                       # 2^(rp+s) = C (mod 3^m)
    per = 2 * 3 ** (m - 1)
    sft = (x - rp) % per
    if sft == 0:
        sft = per
    k = (2 ** (rp + sft) - C) // 3 ** m
    assert (2 ** (rp + sft) - C) % 3 ** m == 0 and k > 0
    return k, sft, word


def part_Q23():
    say("== Q2.3 the unbounded-period counter is an excursion time: l* <= Lam_inf < 2 l* ==")
    ok = True
    ncheck = 0
    rng = random.Random(5)
    starts = list(range(3, 400)) + [rng.randrange(10 ** 9, 10 ** 12) for _ in range(40)] + PATH_RECORDS[-8:]
    maxratio = 0.0
    odd_ok = True
    even_ok = True
    for n in starts:
        vals = orbit_with_tail(n)
        Q = [v & 1 for v in vals]
        tend = next(i for i, v in enumerate(vals) if v == 1)
        li = lam_array(Q, min(len(Q) - 1, 1200))
        ls = lstar_array(Q)
        g = np.where(np.array(Q) == 1, LN3, 0.0) - LN2
        S = np.concatenate([[0.0], np.cumsum(g)])
        for t in range(tend):
            if ls[t] > 0:
                ok &= ls[t] <= li[t] < 2 * ls[t]
                maxratio = max(maxratio, li[t] / ls[t])
                ncheck += 1
            if Q[t] == 1 and ls[t + 1] > 0:
                odd_ok &= ls[t + 1] <= ls[t] - 1
            if Q[t] == 0:
                above2 = [l for l in range(1, len(S) - t - 1) if S[t + 1 + l] - S[t + 1] > LN2 + 1e-12]
                even_ok &= ls[t] == (1 + max(above2) if above2 else 0)
    check(ok, "l*(m) <= Lam_inf(m) < 2 l*(m) at all %d orbit points with l* > 0 (max ratio %.3f)" % (ncheck, maxratio))
    check(odd_ok and even_ok, "on the same orbits: l*(T m) <= l*(m) - 1 at odd m, and l*(m) = 1 + max{l : M_l(m/2) > 2} "
          "(0 if none) at even m")
    agree = 0
    tot = 0
    for n in range(3, 3000):
        vals = orbit_with_tail(n, tail=4)
        Q = [v & 1 for v in vals]
        ls = int(lstar_array(Q)[0])
        last_above = max([i for i, v in enumerate(vals) if v > n], default=0)
        tot += 1
        agree += (ls == last_above)
        if ls > last_above:
            check(False, "l* exceeds the last passage above the start for n = %d" % n)
    check(agree / tot > 0.9, "l*(n) <= last passage of the orbit above n for every 3 <= n < 3000, with equality "
          "for %d of %d (the multiplier bound 3^a/2^l > 1 implies T^l n > n)" % (agree, tot))
    say("   hovering family (Proposition Q2.3(c)): k with parity word w^r 0^s, w = lower Christoffel word of an")
    say("   upper best approximation a/p of log_3 2, r <= log 2 / log(3^a/2^p); violation at 2k -> k >= c r p - log 2")
    rows = []
    for (a, pp, r) in ((2, 3, 1), (2, 3, 3), (2, 3, 5), (7, 11, 1), (7, 11, 2), (12, 19, 1)):
        w = lower_christoffel(a, pp)
        lam = 3 ** a / 2 ** pp
        assert lam > 1 and r <= math.log(2) / math.log(lam)
        k, sft, word = hover_integer(w, r)
        x = k
        for b in word:
            if (x & 1) != b:
                check(False, "hover integer does not follow its word")
            x = (3 * x + 1) >> 1 if x & 1 else x >> 1
        full = word + [0] * sft + [1, 0] * 8
        check_pow = (x == 2 ** sft)
        ls_k = int(lstar_array(full)[0])
        ls_2k = int(lstar_array([0] + full)[0])
        viol_lb = C_DEF * ls_k - LN2
        rows.append((a, pp, r, sft, k.bit_length(), ls_k, ls_2k, viol_lb))
        if not (check_pow and ls_k >= r * pp and ls_2k == 0):
            check(False, "hovering family fails for %d/%d r=%d" % (a, pp, r))
        say("     a/p = %d/%d, r = %d: s = %d, k has %d bits, T^(rp) k = 2^s, l*(k) = %d >= rp = %d, l*(2k) = 0, "
            "so Lam_inf(2k) = 0 and R_inf(k) - R_inf(2k) >= %.2f" % (a, pp, r, sft, k.bit_length(), ls_k, r * pp, viol_lb))
    check(len(rows) == 6, "the hovering integers exist and behave as proved for 2/3 (r = 1, 3, 5), 7/11 (r = 1, 2), 12/19 "
          "(r = 1); the proved bound c p floor(log 2/log lambda) - log 2 is %.1f, %.1f, %.1f, %.1f for 2/3, 7/11, 12/19, "
          "53/84, unbounded as a/p -> log_3 2" %
          tuple(C_DEF * pp * math.floor(math.log(2) / math.log(3 ** a / 2 ** pp)) - LN2
                for a, pp in ((2, 3), (7, 11), (12, 19), (53, 84))))


def collect(starts, label, P_list=(1, 3, 20), thetas=(2, 4, 8, 16)):
    agg = {}
    npts = 0
    reset_only = True
    for n in starts:
        out, tmax = orbit_stats(n, P_list=P_list, thetas=thetas)
        npts += tmax
        big = out["big"][0] if "big" in out else None
        for key, (d, dl) in out.items():
            if key in ("lstar", "big"):
                continue
            viol = d > 1e-9
            nv = int(viol.sum())
            mx = float(d.max()) if len(d) else -1e9
            a = agg.setdefault(key, [0, -1e9])
            a[0] += nv
            a[1] = max(a[1], mx)
            if dl is not None and key.startswith("A"):
                vb = viol & big
                if vb.any() and np.any(dl[vb] < 0):
                    reset_only = False
    return agg, npts, reset_only


def part_stats():
    say("== Q2.4 violation statistics (every step m -> T m with m >= 3 before reaching 1) ==")
    rng = random.Random(2026)
    sets = [("orbits of 3..5000", list(range(3, 5001))),
            ("path+delay records <= 10^12 (A006884, A006877)", sorted(set(PATH_RECORDS + DELAY_RECORDS) - {1, 2})),
            ("600 random n in [1e11, 1e12]", [rng.randrange(10 ** 11, 10 ** 12) for _ in range(600)])]
    keys = ["A1", "A3", "A20", "Ainf", "B3_10", "D2", "D4", "D8", "D16", "C2", "C4", "C8", "C16"]
    say("   %-48s %8s  " % ("set", "steps") + "  ".join("%-15s" % k for k in keys))
    results = {}
    for label, starts in sets:
        agg, npts, reset_only = collect(starts, label)
        results[label] = (agg, npts, reset_only)
        say("   %-48s %8d  " % (label, npts) + "  ".join("%5d/%-9.2f" % (agg[k][0], agg[k][1]) for k in keys))
    for label, (agg, npts, reset_only) in results.items():
        check(reset_only, "%s: every violation of A_1, A_3, A_20, A_inf (c = 1.1 log(3/2)) at a point m >= 11 happens "
              "at a reset (Delta Lam >= 0), never inside a shadow" % label)
        check(all(agg[k][0] > 0 for k in ("A1", "A3", "A20", "Ainf", "B3_10")),
              "%s: the static and bounded-window candidates A_1, A_3, A_20, A_inf, B_3,10 all have violations "
              "(max increments %s)" % (label, ", ".join("%.2f" % agg[k][1] for k in ("A1", "A3", "A20", "Ainf", "B3_10"))))
    return results


def D_violations(vals, tmax, th):
    logs = np.log(np.array([float(v) for v in vals]))
    L = len(vals)
    w = np.floor(th * logs / LN2).astype(np.int64)
    RD = np.array([logs[t:min(L, t + w[t] + 1)].max() for t in range(tmax + 1)])
    return int((RD[1:] - RD[:-1] > 1e-12).sum())


def part_theta():
    say("== Q2.5 windowed future envelopes: the horizon needed on the tested sets ==")
    rng = random.Random(99)
    sets = [("orbits of 3..5000", list(range(3, 5001))),
            ("records <= 10^12", sorted(set(PATH_RECORDS + DELAY_RECORDS) - {1, 2})),
            ("random [1e11,1e12]", [rng.randrange(10 ** 11, 10 ** 12) for _ in range(300)])]
    for label, starts in sets:
        th = 0.0
        data = []
        for n in starts:
            vals = orbit_with_tail(n, tail=8)
            tend = next(i for i, v in enumerate(vals) if v == 1)
            tmax = max([i for i in range(tend) if vals[i] >= 3], default=-1) + 1
            th = max(th, peak_time_ratio(vals, tmax))
            data.append((vals, tmax))
        thr = math.ceil(th * 100) / 100 + 0.01
        nv = sum(D_violations(v, t, thr) for v, t in data)
        nv_half = sum(D_violations(v, t, th / 2) for v, t in data)
        check(nv == 0 and nv_half > 0,
              "%s: max over orbit points m of (time to the future peak)/log2 m = %.3f; the windowed envelope D_theta "
              "is non-increasing on every step for theta = %.2f, and has %d violations for theta = %.3f" %
              (label, th, thr, nv_half, th / 2))


def part_controls(c=C_DEF):
    say("== Q2.6 cycle controls (3n-1 and 5n+1) ==")
    for q, sgn, start, name in ((5, 1, 13, "5n+1 cycle through 13"), (5, 1, 17, "5n+1 cycle through 17")):
        x = start
        per = 0
        while True:
            x = (q * x + sgn) >> 1 if x & 1 else x >> 1
            per += 1
            if x == start:
                break
        vals = orbit_with_tail(start, q=q, sgn=sgn, steps=per * 40)
        Q = [v & 1 for v in vals]
        lam = lam_array(Q, 20, q=q)
        R = np.log(np.array(vals[:per + 1], dtype=float)) + c * lam[:per + 1]
        d = R[1:] - R[:-1]
        check(lam[0] == lam[per] and abs(d.sum()) < 1e-9 and d.max() > 0,
              "%s (period %d, contracting: 5^a < 2^p): A_20 (with the 5n+1 expanding words) is finite and periodic "
              "on the cycle, its increments sum to 0, max %.3f > 0: flagged" % (name, per, d.max()))
    for start, name in ((5, "3n-1 cycle 5,7,10"), (17, "3n-1 cycle through 17")):
        lens = []
        for reps in (20, 40, 80):
            x = start
            per = 0
            while True:
                x = (3 * x - 1) >> 1 if x & 1 else x >> 1
                per += 1
                if x == start:
                    break
            vals = orbit_with_tail(start, q=3, sgn=-1, steps=per * reps)
            lam = lam_array([v & 1 for v in vals], 20)
            lens.append((len(vals), int(lam[0])))
        check(all(l0 == L for L, l0 in lens),
              "%s: the cycle is an expanding periodic point that is a positive integer, so Lam_20 = infinity there "
              "(the run equals the whole computed sequence: %s): the candidate is not finite on this map" %
              (name, ", ".join("%d/%d" % (l0, L) for L, l0 in lens)))


def part_worst(c=C_DEF):
    say("== Q2.7 where the largest violations of A_20 and A_inf occur (records <= 10^12) ==")
    best = {"A20": (-1e9, None), "Ainf": (-1e9, None)}
    for n in sorted(set(PATH_RECORDS + DELAY_RECORDS) - {1, 2}):
        vals = orbit_with_tail(n, tail=64)
        Q = [v & 1 for v in vals]
        L = len(vals)
        tend = next(i for i, v in enumerate(vals) if v == 1)
        logs = np.log(np.array([float(v) for v in vals]))
        l20 = lam_array(Q, 20)
        li = lam_array(Q, L - 1)
        ls = lstar_array(Q)
        for key, lam in (("A20", l20), ("Ainf", li)):
            R = logs + c * lam
            d = R[1:tend] - R[:tend - 1]
            t = int(np.argmax(d))
            if d[t] > best[key][0]:
                best[key] = (float(d[t]), (n, t, vals[t], Q[t], int(lam[t]), int(lam[t + 1]), int(ls[t]), int(ls[t + 1])))
    for key, (dv, info) in best.items():
        n, t, m, par, l0, l1, s0, s1 = info
        say("   %s: largest increment %.2f at step %d of the orbit of %d: m = %d (%s step), Lam %d -> %d, l* %d -> %d"
            % (key, dv, t, n, m, "odd" if par else "even", l0, l1, s0, s1))
    n, t, m, par, l0, l1, s0, s1 = best["Ainf"][1]
    check(par == 0 and l1 - l0 > 100 and s1 - s0 > 50,
          "the largest A_inf violation is an even step m -> m/2 at which the excursion time l* jumps by %d "
          "(the multiplier of m/2 stays in (1, 2] for a long stretch after its last visit above 2)" % (s1 - s0))
    n, t, m, par, l0, l1, s0, s1 = best["A20"][1]
    check(l1 - l0 >= 20, "the largest A_20 violation is a reset: Lam_20 jumps from %d to %d" % (l0, l1))


def run():
    part_isometry()
    part_Q21()
    part_Q22()
    part_Q23()
    part_stats()
    part_theta()
    part_controls()
    part_worst()


if __name__ == "__main__":
    run()
