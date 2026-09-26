#!/usr/bin/env python3
"""collatz_nodescent_order_20260926_audit.py -- independent adversarial audit of THM-4495
(the Spitzer/ballot identity k W_k = sum_n B_n W_(k-n) for the no-descent count W_k = |Bad_k|
and its exact order 2^(hk) k^(-3/2)).  Written blind to the audited script's code paths:
every object is re-implemented from the definitions in the note, and every comparison of two
partial sums S = o log_2 3 - j is done EXACTLY (sign of 3^o 2^j' - 3^o' 2^j), never in floats.

Audited: 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md and
         05-knowledge/results/collatz_nodescent_order_20260926_spitzer_ballot.md
         (script 04-computation/experiments/collatz_nodescent_order_20260926.py, output .out).

Sections (numbered as in the audit brief):
   1  no-ties setup, definitions, reversal formula S_i(w^R) = S_m(w) - S_(m-i)(w)
   2  Theorem A Step 1: minimum decomposition, bijection checked word by word (n <= 12)
   3  Step 2: reversal P_m = L_m, set equality of reversed min-ending and first-passage words
   4  Step 3: ladder blocks, E_(n,r) = [t^n] L(t)^r, direct counts vs polynomial powers (n <= 14)
   5  Step 4: rotation averaging, claims (a)-(d) word by word (n <= 12); B_n = n sum_r E_(n,r)/r exact;
      [t^n] log W = B_n/n from the DP series by an exact series logarithm (n <= 40)
   6  independent W_k (own ballot DP, brute force k <= 18, THM-4479 residue simulation k <= 16, and a
      first-passage DP through W = 1/(1-L)), B_n, recurrence k <= KMAX exact, THM-4479's table k <= 11
   7  Theorem B: (a) geometric tail, (b) Stirling bounds for ALL 1 <= j <= n-1 (n <= KMAX) + Robbins
      derivation, (c) h(j0/n) <= h, (d) b_n <= 2.05 n^(-1/2), (e) convolution lemma and induction step,
      (f) sigma and the constant 545, (g) the lower constant 0.26, (h) oscillation window and Burnside
   8  Corollary C3: carry bound |beta_j| <= (3/2)^j - 1, the four cases, brute-force dip counts
      D_+-(2^T, 1) on both sheets vs sum_(t<T) W_t, per-block classification
   9  C1/C2: the sandwich lines quoted from THM-4479 / THM-4485 / THM-4487 (and THM-4479's code)
  10  provenance: sha256 of the audited script and output vs the theorem header (raw and LF bytes)
Usage: python3 collatz_nodescent_order_20260926_audit.py [KMAX=3000] [TFULL=18] [TCOUNT=20]
"""
import hashlib
import math
import os
import sys
import time
from fractions import Fraction
from itertools import product
from math import comb, gcd

import mpmath

mpmath.mp.dps = 50
MP = mpmath.mpf
RHO = mpmath.log(2) / mpmath.log(3)                     # log_3 2
LOG2_3 = mpmath.log(3) / mpmath.log(2)                  # log_2 3
H = -(RHO * mpmath.log(RHO, 2) + (1 - RHO) * mpmath.log(1 - RHO, 2))   # h(rho), bits
RHOf, Hf = float(RHO), float(H)

FAILS = []
T0 = time.time()


def check(cond, msg):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAILS.append(msg)
    print("  [%s] %s" % (tag, msg))
    return cond


def hbits(p):
    p = float(p)
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def hbits_mp(p):
    p = MP(p)
    return -(p * mpmath.log(p, 2) + (1 - p) * mpmath.log(1 - p, 2))


# ---------------------------------------------------------------- exact walk arithmetic
def sign3_2(o, j):
    """sign of o*log_2(3) - j, exactly (0 only for o = j = 0)."""
    if o == 0 and j == 0:
        return 0
    if o >= 0 and j <= 0:
        return 1
    if o <= 0 and j >= 0:
        return -1
    if o > 0 and j > 0:
        a, b = 3 ** o, 2 ** j
        return 1 if a > b else (-1 if a < b else 0)
    return -sign3_2(-o, -j)


def cmpS(o1, j1, o2, j2):
    """sign of S(o1, j1) - S(o2, j2), S(o, j) = o log_2 3 - j."""
    return sign3_2(o1 - o2, j1 - j2)


def prefix_ones(w):
    o = [0]
    for x in w:
        o.append(o[-1] + x)
    return o


def is_positive(w):
    o = prefix_ones(w)
    return all(sign3_2(o[j], j) > 0 for j in range(1, len(w) + 1))


def argmin_index(w):
    """index m in [1, n] of the unique minimum of S_j over 1 <= j <= n."""
    o = prefix_ones(w)
    m = 1
    for j in range(2, len(w) + 1):
        if cmpS(o[j], j, o[m], m) < 0:
            m = j
    return m


def is_min_ending(w):
    """S_j > S_n for all 1 <= j < n."""
    o = prefix_ones(w)
    n = len(w)
    return all(cmpS(o[j], j, o[n], n) > 0 for j in range(1, n))


def is_first_passage(w):
    """S_j < 0 for 1 <= j < n and S_n > 0."""
    o = prefix_ones(w)
    n = len(w)
    return n >= 1 and all(sign3_2(o[j], j) < 0 for j in range(1, n)) and sign3_2(o[n], n) > 0


def ends_at_max(w):
    """S_n > S_j for all 0 <= j < n."""
    o = prefix_ones(w)
    n = len(w)
    return n >= 1 and all(cmpS(o[n], n, o[j], j) > 0 for j in range(0, n))


def ladder_epochs(w):
    """strict ascending ladder epochs: j in [1, n] with S_j > S_l for all 0 <= l < j."""
    o = prefix_ones(w)
    eps = []
    mo, mj = 0, 0
    for j in range(1, len(w) + 1):
        if cmpS(o[j], j, mo, mj) > 0:
            eps.append(j)
            mo, mj = o[j], j
    return eps


def all_words(n):
    return product((0, 1), repeat=n)


# ================================================================ 1. setup
def section1():
    print("== 1. setup: no ties, definitions, reversal formula ==")
    # 3^a = 2^b only for a = b = 0 (unique factorisation); exact search as a sanity check of the comparator
    ties = [(a, b) for a in range(0, 200) for b in range(0, 320) if 3 ** a == 2 ** b]
    check(ties == [(0, 0)], "3^a = 2^b only at a = b = 0 (a < 200, b < 320): no two partial sums of any word coincide, none vanishes")
    # no two partial sums S_0..S_n of any word of length <= 12 coincide (exact comparator)
    bad = 0
    for n in range(1, 13):
        for w in all_words(n):
            o = prefix_ones(w)
            vals = [(o[j], j) for j in range(n + 1)]
            for a in range(n + 1):
                for b in range(a + 1, n + 1):
                    if cmpS(*vals[a], *vals[b]) == 0:
                        bad += 1
    check(bad == 0, "exact comparator: 0 coinciding pairs among S_0..S_n over all words of length <= 12")
    # reversal formula S_i(w^R) = S_m(w) - S_(m-i)(w): as (o, j)-pairs, o_i(w^R) = o_m(w) - o_(m-i)(w)
    bad = 0
    for m in range(1, 11):
        for w in all_words(m):
            o = prefix_ones(w)
            oR = prefix_ones(w[::-1])
            for i in range(m + 1):
                if oR[i] != o[m] - o[m - i]:
                    bad += 1
    check(bad == 0, "reversal formula S_i(w^R) = S_m(w) - S_(m-i)(w) holds for all words of length <= 10 and all 0 <= i <= m")
    print("  definitions used below: positive = S_j > 0 for 1 <= j <= n; first-passage = S_j < 0 for 1 <= j < m and S_m > 0;")
    print("  ends-at-max = S_n > S_j for 0 <= j < n; min-ending = S_j > S_n for 1 <= j < n (index j = 0 excluded, as positivity forces).")


# ================================================================ 2. Step 1
def section2(NMAX=12):
    print("== 2. Theorem A Step 1: minimum decomposition W_k = sum_m P_m W_(k-m) ==")
    POS = {0: [()]}
    PM = {}
    for n in range(1, NMAX + 1):
        POS[n] = [w for w in all_words(n) if is_positive(w)]
        PM[n] = [w for w in POS[n] if is_min_ending(w)]
    ok_fwd = ok_back = ok_bij = ok_rec = True
    for k in range(1, NMAX + 1):
        seen = set()
        for w in POS[k]:
            m = argmin_index(w)
            a, b = w[:m], w[m:]
            if not (1 <= m <= k and is_positive(a) and is_min_ending(a) and (len(b) == 0 or is_positive(b))):
                ok_fwd = False
            seen.add((a, b))
        pairs = set()
        for m in range(1, k + 1):
            for a in PM[m]:
                for b in POS[k - m]:
                    ab = a + b
                    if not (is_positive(ab) and argmin_index(ab) == m):
                        ok_back = False
                    pairs.add((a, b))
        if pairs != seen or len(pairs) != len(POS[k]):
            ok_bij = False
        if sum(len(PM[m]) * len(POS[k - m]) for m in range(1, k + 1)) != len(POS[k]):
            ok_rec = False
    check(ok_fwd, "forward: every positive w of length k <= %d splits at its minimum index m in [1,k] into a positive min-ending a and a positive (possibly empty) b" % NMAX)
    check(ok_back, "converse: for every positive min-ending a and positive b, ab is positive with minimum index exactly |a|")
    check(ok_bij, "the two maps are mutually inverse (the set of pairs (a,b) from splitting equals PM_m x POS_(k-m) over m), so the decomposition is unique")
    check(ok_rec, "W_k = sum_(m=1)^k P_m W_(k-m) with W_0 = 1 (empty b allowed, m = k term) for k <= %d" % NMAX)
    print("  W_k direct, k = 1..%d:" % NMAX, [len(POS[k]) for k in range(1, NMAX + 1)])
    print("  P_m direct (positive min-ending), m = 1..%d:" % NMAX, [len(PM[m]) for m in range(1, NMAX + 1)])
    return POS, PM


# ================================================================ 3. Step 2
def section3(PM, NMAX=12):
    print("== 3. Theorem A Step 2: reversal, P_m = L_m ==")
    ok = True
    L = []
    for m in range(1, NMAX + 1):
        fp = set(w for w in all_words(m) if is_first_passage(w))
        rev = set(a[::-1] for a in PM[m])
        L.append(len(fp))
        if fp != rev:
            ok = False
    check(ok, "{a^R : a positive min-ending of length m} == {first-passage words of length m} for m <= %d (both directions)" % NMAX)
    print("  L_m direct (first-passage), m = 1..%d:" % NMAX, L)
    return L


# ================================================================ 4. Step 3
def section4(NMAX=14):
    print("== 4. Theorem A Step 3: ladder blocks, E_(n,r) = [t^n] L(t)^r ==")
    Lc = [0] * (NMAX + 1)
    E = {}
    ok_blocks = True
    for n in range(1, NMAX + 1):
        for w in all_words(n):
            if is_first_passage(w):
                Lc[n] += 1
            eps = ladder_epochs(w)
            if ends_at_max(w):
                if not eps or eps[-1] != n:
                    ok_blocks = False
                prev = 0
                for e in eps:
                    if not is_first_passage(w[prev:e]):
                        ok_blocks = False
                    prev = e
                E[(n, len(eps))] = E.get((n, len(eps)), 0) + 1
            else:
                if eps and eps[-1] == n:
                    ok_blocks = False
    check(ok_blocks, "a word ends at its maximum iff n is a ladder epoch, and then the blocks between consecutive ladder epochs are all first-passage (n <= %d)" % NMAX)
    # polynomial powers of L(t)
    ok_pow = True
    Lr = [0] + Lc[1:]
    Lpow = {1: Lr[:]}
    for r in range(2, NMAX + 1):
        prev = Lpow[r - 1]
        cur = [0] * (NMAX + 1)
        for i in range(1, NMAX + 1):
            if prev[i]:
                for m in range(1, NMAX + 1 - i):
                    cur[i + m] += prev[i] * Lr[m]
        Lpow[r] = cur
    for n in range(1, NMAX + 1):
        for r in range(1, n + 1):
            if E.get((n, r), 0) != Lpow[r][n]:
                ok_pow = False
    check(ok_pow, "E_(n,r) (direct count of words ending at max with r ladder epochs) == [t^n] L(t)^r for all n <= %d, r <= n" % NMAX)
    # converse: every concatenation of first-passage words ends at its max with r = number of blocks ladder epochs (n <= 10)
    FP = {m: [w for w in all_words(m) if is_first_passage(w)] for m in range(1, 11)}
    ok_conv = True

    def compositions(n):
        if n == 0:
            yield ()
            return
        for first in range(1, n + 1):
            for rest in compositions(n - first):
                yield (first,) + rest

    for n in range(1, 11):
        for comp in compositions(n):
            for tup in product(*[FP[m] for m in comp]):
                y = sum(tup, ())
                if not (ends_at_max(y) and len(ladder_epochs(y)) == len(comp)):
                    ok_conv = False
    check(ok_conv, "converse: every concatenation of r first-passage words ends at its maximum with exactly r ladder epochs (n <= 10, all tuples)")
    logW = [Fraction(0)] * (NMAX + 1)   # [t^n] log 1/(1-L) = sum_r E_(n,r)/r
    for n in range(1, NMAX + 1):
        logW[n] = sum(Fraction(E.get((n, r), 0), r) for r in range(1, n + 1))
    print("  sum_r E_(n,r)/r for n = 1..%d:" % NMAX, [str(x) for x in logW[1:]])
    return E, logW


# ================================================================ 5. Step 4
def section5(E, logW, NMAX=12):
    print("== 5. Theorem A Step 4: rotation averaging (the crux) ==")
    ok_a = ok_b = ok_c = ok_d = ok_win = True
    n_words = 0
    for n in range(1, NMAX + 1):
        Bn = 0
        total = Fraction(0)
        for x in all_words(n):
            o = prefix_ones(x)
            on = o[n]
            if sign3_2(on, n) <= 0:
                continue
            Bn += 1
            n_words += 1

            def ext(j):
                q, rr = divmod(j, n)
                return (o[rr] + q * on, j)

            def is_record(j, window=n):
                oj, jj = ext(j)
                return all(cmpS(oj, jj, *ext(l)) > 0 for l in range(j - window, j))

            R = [i for i in range(n) if is_record(i)]
            # (a) records exist and are periodic; the window [j-n, j) test agrees with a wider window [j-3n, j)
            if len(R) < 1:
                ok_a = False
            for i in range(-n, 2 * n):
                if is_record(i) != is_record(i + n):
                    ok_a = False
                if is_record(i) != is_record(i, 3 * n):
                    ok_win = False
            s = Fraction(0)
            for i in range(n):
                rot = x[i:] + x[:i]
                em = ends_at_max(rot)
                # (b)
                if em != (i in R):
                    ok_b = False
                if em:
                    eps = ladder_epochs(rot)
                    # (c)
                    if eps != [l for l in range(1, n + 1) if is_record(i + l)] or len(eps) != len(R):
                        ok_c = False
                    s += Fraction(1, len(eps))
            # (d) per word
            if s != 1:
                ok_d = False
            if ends_at_max(x):
                total += Fraction(1, len(ladder_epochs(x)))
        # (d) global: B_n = n sum_(y ends at max) 1/r(y) = n sum_r E_(n,r)/r
        if Bn != n * total or total != logW[n]:
            ok_d = False
    check(ok_a, "(a) for every x with S_n(x) > 0, n <= %d: records of the periodic extension exist and j is a record iff j + n is" % NMAX)
    check(ok_win, "    record test over the window [j-n, j) agrees with the window [j-3n, j) (the reduction to one period is sound)")
    check(ok_b, "(b) rot_i(x) ends at its maximum iff i is a record residue (checked for every i, every x)")
    check(ok_c, "(c) for record i, the ladder epochs of rot_i(x) are exactly the l in [1,n] with i + l a record; their number is r(x)")
    check(ok_d, "(d) sum_i [rot_i x ends at max]/r(rot_i x) = 1 for every x; hence B_n = n sum_r E_(n,r)/r exactly (Fraction) for n <= %d" % NMAX)
    print("  words with positive total examined: %d" % n_words)


# ================================================================ 6. independent W_k, B_n, recurrence
def W_ballot_dp(K):
    """own DP: states = number of ones o after j letters, kept only while 3^o > 2^j."""
    pow3 = [3 ** o for o in range(K + 2)]
    pow2 = [2 ** j for j in range(K + 2)]
    W = [1]
    cur = {0: 1}
    for j in range(1, K + 1):
        nxt = {}
        for o, c in cur.items():
            for s in (0, 1):
                o2 = o + s
                if pow3[o2] > pow2[j]:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
        W.append(sum(cur.values()))
    return W


def L_firstpassage_dp(K):
    """own DP for first-passage counts L_m: states o with 3^o < 2^j (walk strictly negative) until the first j with 3^o > 2^j."""
    pow3 = [3 ** o for o in range(K + 2)]
    pow2 = [2 ** j for j in range(K + 2)]
    L = [0] * (K + 1)
    cur = {0: 1}
    for j in range(1, K + 1):
        nxt = {}
        for o, c in cur.items():
            for s in (0, 1):
                o2 = o + s
                if pow3[o2] > pow2[j]:
                    L[j] += c
                else:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
    return L


def W_bruteforce(k):
    return sum(1 for w in all_words(k) if is_positive(w))


def W_residue_simulation(k):
    """THM-4479's definition literally: residues r mod 2^k (representatives 0..2^k-1), a_j = odd terms among
    r, T r, ..., T^(j-1) r, M_j = 3^(a_j)/2^j > 1 for j = 1..k (strict)."""
    pow3 = [3 ** o for o in range(k + 2)]
    cnt = 0
    for r in range(2 ** k):
        x = r
        a = 0
        ok = True
        for j in range(1, k + 1):
            if x & 1:
                a += 1
                x = (3 * x + 1) // 2
            else:
                x //= 2
            if not (pow3[a] > (1 << j)):
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def B_tails(K):
    B = [0] * (K + 1)
    J0 = [0] * (K + 1)
    for n in range(1, K + 1):
        j0 = 0
        while 3 ** j0 <= 2 ** n:
            j0 += 1
        J0[n] = j0
        c = comb(n, j0)
        s = 0
        for j in range(j0, n + 1):
            s += c
            c = c * (n - j) // (j + 1)
        B[n] = s
    return B, J0


def section6(KMAX):
    print("== 6. independent W_k, B_n and the recurrence ==")
    W = W_ballot_dp(KMAX)
    B, J0 = B_tails(KMAX)
    wb = [W_bruteforce(k) for k in range(1, 19)]
    check(wb == W[1:19], "own ballot DP == brute-force enumeration of positive words for k <= 18: %s" % wb)
    ws = [W_residue_simulation(k) for k in range(1, 17)]
    check(ws == W[1:17], "own ballot DP == THM-4479's residue simulation |{r mod 2^k : M_j(r) > 1, j = 1..k}| for k <= 16: %s" % ws)
    L = L_firstpassage_dp(300)
    WL = [1]
    for k in range(1, 301):
        WL.append(sum(L[m] * WL[k - m] for m in range(1, k + 1)))
    check(WL == W[:301], "W = 1/(1 - L) with L from an independent first-passage DP reproduces W_k for k <= 300 (Steps 1-2 at scale)")
    table = {2: 1, 3: 2, 4: 3, 5: 4, 6: 8, 7: 13, 8: 19, 9: 38, 10: 64, 11: 128}
    check(all(W[k] == v for k, v in table.items()) and W[1] == 1, "W_k == THM-4479's |Bad_k| table for k = 2..11 (1,2,3,4,8,13,19,38,64,128); W_1 = 1 = |Bad_1| = |{1 mod 2}|")
    print("  THM-4479 defines Bad_k = {r mod 2^k : M_j(r) > 1 for j = 1..k}, M_j = 3^(a_j)/2^j (strict; equality never occurs);")
    print("  its code tests `mul ** a > 2 ** j` (procgen_cubedist_20260925_lib.py) and `3 ** a2 > 2 ** j` (orchestrator ballot DP)")
    print("  over representatives 0..2^k-1: identical to W_k's definition 3^(o_j) > 2^j. Definitions AGREE.")
    ok250 = all(k * W[k] == sum(B[n] * W[k - n] for n in range(1, k + 1)) for k in range(1, 251))
    check(ok250, "k W_k == sum_(n=1)^k B_n W_(k-n) EXACTLY (own W by DP, own B by binomial tails) for all k <= 250")
    okK = all(k * W[k] == sum(B[n] * W[k - n] for n in range(1, k + 1)) for k in range(251, KMAX + 1))
    check(okK, "... and for all 251 <= k <= %d" % KMAX)
    print("  W_k, k = 1..20:", W[1:21])
    print("  B_n, n = 1..20:", B[1:21])
    print("  j0(n) = least j with 3^j > 2^n, n = 1..20:", J0[1:21])
    # exact series logarithm of W(t) from the DP: n c_n = n w_n - sum_(m<n) m c_m w_(n-m); claim c_n = B_n/n
    N = 40
    c = [Fraction(0)] * (N + 1)
    for n in range(1, N + 1):
        c[n] = Fraction(W[n]) - sum((m * c[m] * W[n - m] for m in range(1, n)), Fraction(0)) / n
    check(all(c[n] == Fraction(B[n], n) for n in range(1, N + 1)), "[t^n] log W(t) == B_n/n exactly for n <= %d (series log of the DP generating function)" % N)
    return W, B, J0


# ================================================================ 7. Theorem B
def necklaces_burnside(k, pow3):
    """own Burnside: # binary necklaces of length k with j ones, summed over j with 3^j > 2^k."""
    def phi(m):
        r, x, p = m, m, 2
        while p * p <= x:
            if x % p == 0:
                while x % p == 0:
                    x //= p
                r -= r // p
            p += 1
        if x > 1:
            r -= r // x
        return r
    p2 = 2 ** k
    tot = 0
    for j in range(0, k + 1):
        if pow3[j] > p2:
            g = gcd(k, j)
            s = 0
            for d in range(1, g + 1):
                if g % d == 0:
                    s += phi(d) * comb(k // d, j // d)
            assert s % k == 0
            tot += s // k
    return tot


def necklaces_direct(k):
    pow3 = [3 ** o for o in range(k + 2)]
    p2 = 2 ** k
    seen = set()
    cnt = 0
    for w in all_words(k):
        if w in seen:
            continue
        rots = [w[i:] + w[:i] for i in range(k)]
        for r in rots:
            seen.add(r)
        if pow3[sum(w)] > p2:
            cnt += 1
    return cnt


def section7(W, B, J0, KMAX):
    print("== 7. Theorem B ==")
    pow3 = [3 ** o for o in range(KMAX + 2)]
    q = (1 - RHO) / RHO
    tailconst = RHO / (2 * RHO - 1)
    print("  rho = %s, h = %s, (1-rho)/rho = %s, rho/(2rho-1) = %s" % (mpmath.nstr(RHO, 12), mpmath.nstr(H, 12), mpmath.nstr(q, 10), mpmath.nstr(tailconst, 10)))
    # (a) j0 = floor(n rho) + 1; C(n, j0) <= B_n <= (rho/(2rho-1)) C(n, j0); j0 <= n-1 for n >= 3
    ok_j0 = all(J0[n] == int(mpmath.floor(n * RHO)) + 1 for n in range(1, KMAX + 1))
    check(ok_j0, "(a) j0(n) = floor(n rho) + 1 for all n <= %d (least j with 3^j > 2^n)" % KMAX)
    check(all(J0[n] <= n - 1 for n in range(3, KMAX + 1)) and J0[1] == 1 and J0[2] == 2, "(a) j0 <= n - 1 for n >= 3 (p < 1); j0 = n for n = 1, 2")
    worst = MP(0)
    for n in range(1, KMAX + 1):
        r = MP(B[n]) / MP(comb(n, J0[n]))
        if r > worst:
            worst = r
        if comb(n, J0[n]) > B[n]:
            check(False, "(a) C(n,j0) <= B_n fails at n = %d" % n)
    check(worst <= tailconst, "(a) max_n B_n/C(n,j0) = %s <= rho/(2rho-1) = %s (ratio to the bound %s); C(n,j0) <= B_n trivially" % (mpmath.nstr(worst, 8), mpmath.nstr(tailconst, 8), mpmath.nstr(worst / tailconst, 6)))
    # analytic: (n-j)/(j+1) < (1-rho)/rho iff j > rho(n+1) - 1, implied by j >= j0 > n rho
    ok_ratio = all(MP(n - j) / (j + 1) < q for n in range(1, KMAX + 1) for j in range(J0[n], n))
    check(ok_ratio, "(a) C(n,j+1)/C(n,j) = (n-j)/(j+1) < (1-rho)/rho for all j >= j0, n <= %d (algebra: iff j > rho(n+1)-1, and j0 > n rho > rho(n+1)-1)" % KMAX)
    # (b) Stirling-type bounds for ALL 1 <= j <= n-1, n <= KMAX (exact binomials, float logs; margins are >> 1e-9)
    min_low, max_up = 1e9, -1e9
    argmin_low = argmax_up = None
    min_low_j0 = 1e9
    max_up_j0 = -1e9
    robbins_min = 1e9
    for n in range(2, KMAX + 1):
        c = n
        for j in range(1, n):
            lc = math.log2(c)
            p = j / n
            hp = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
            v = n * p * (1 - p)
            low = n * hp - 0.5 * math.log2(8 * v)
            up = n * hp - 0.5 * math.log2(2 * math.pi * v)
            dl, du = lc - low, lc - up
            if dl < min_low:
                min_low, argmin_low = dl, (n, j)
            if du > max_up:
                max_up, argmax_up = du, (n, j)
            if j == J0[n] and n >= 3:
                min_low_j0 = min(min_low_j0, dl)
                max_up_j0 = max(max_up_j0, du)
            if n >= 3:
                e = 1.0 / (12 * n + 1) - 1.0 / (12 * j) - 1.0 / (12 * (n - j))
                robbins_min = min(robbins_min, e)
            c = c * (n - j) // (j + 1)
    check(min_low >= -1e-9, "(b) lower Stirling bound C(n,pn) >= 2^(nh(p))/sqrt(8np(1-p)): min margin %.6f bits at (n,j) = %s over ALL 1 <= j <= n-1, n <= %d (equality at n = 2, j = 1)" % (min_low, argmin_low, KMAX))
    check(max_up <= 1e-9, "(b) upper Stirling bound C(n,pn) <= 2^(nh(p))/sqrt(2 pi np(1-p)): max margin %.6f bits at (n,j) = %s over ALL 1 <= j <= n-1, n <= %d" % (max_up, argmax_up, KMAX))
    print("      at the j = j0(n) that occur (n >= 3): min(log2 C - lower) = %.4f, max(log2 C - upper) = %.6f (note: 0.0376 and 0.0000)" % (min_low_j0, max_up_j0))
    check(robbins_min > math.log(math.sqrt(math.pi / 4)), "(b) Robbins derivation: 1/(12n+1) - 1/(12j) - 1/(12(n-j)) >= %.5f > ln sqrt(pi/4) = %.5f for n >= 3 (so the lower bound follows from n! = sqrt(2 pi n)(n/e)^n e^(r_n), 1/(12n+1) < r_n < 1/(12n));"
          " the upper bound follows from 1/(12j+1) > 1/(12n) for j <= n-1. Source: MacWilliams-Sloane, Theory of Error-Correcting Codes, Ch. 10 Lemma 7; Gallager 1968 Problem 5.8" % (robbins_min, math.log(math.sqrt(math.pi / 4))))
    # (c) h(j0/n) <= h(rho)
    worst_h = MP(-1)
    for n in range(3, KMAX + 1):
        d = hbits_mp(MP(J0[n]) / n) - H
        if d > worst_h:
            worst_h = d
    check(worst_h < 0, "(c) h(j0/n) < h(rho) for all 3 <= n <= %d: max difference %s (j0/n > rho > 1/2, h decreasing on [1/2,1])" % (KMAX, mpmath.nstr(worst_h, 4)))
    # (d) b_n = B_n 2^(-hn) <= 2.05 n^(-1/2) for all n
    b = [MP(0)] * (KMAX + 1)
    bs_max, bs_arg = MP(0), None
    bs_max30 = MP(0)
    for n in range(1, KMAX + 1):
        b[n] = MP(B[n]) * mpmath.power(2, -H * n)
        v = b[n] * mpmath.sqrt(n)
        if v > bs_max:
            bs_max, bs_arg = v, n
        if n >= 30 and v > bs_max30:
            bs_max30 = v
    pmax = RHO + MP(1) / 30
    const30 = tailconst / mpmath.sqrt(2 * mpmath.pi * pmax * (1 - pmax))
    print("      max_n b_n sqrt(n) = %s at n = %d; max over n >= 30 = %s; analytic constant for n >= 30: rho/(2rho-1)/sqrt(2 pi p(1-p)) with p = rho + 1/30 = %s, p(1-p) = %s, constant = %s"
          % (mpmath.nstr(bs_max, 6), bs_arg, mpmath.nstr(bs_max30, 6), mpmath.nstr(pmax, 6), mpmath.nstr(pmax * (1 - pmax), 6), mpmath.nstr(const30, 6)))
    check(bs_max <= MP("2.05") and const30 <= MP("2.05") and bs_max30 <= const30, "(d) b_n <= 2.05 n^(-1/2): numerically for n <= %d (max %s < 1.983 claimed? %s) and analytically for n >= 30 (constant %s <= 2.05)" % (KMAX, mpmath.nstr(bs_max, 5), bs_max <= MP("1.983"), mpmath.nstr(const30, 5)))
    # (e) convolution lemma and the induction step
    ok_alg = True
    for m in range(1, 41):
        lhs = 2 ** 1.5 * (1 + m * 2 ** (1.5 * (m - 1)))
        rhs = (m + 1) * 2 ** (1.5 * m)
        if lhs > rhs * (1 + 1e-12):
            ok_alg = False
    check(ok_alg, "(e) induction step 2^(3/2)(C s^m + C m (2^(3/2) s)^(m-1) s) <= C(m+1)(2^(3/2) s)^m reduces to 2^(3/2) <= 2^(3m/2), true for m >= 1 (checked m <= 40)")
    print("      convolution lemma: split at i = k/2; i <= k/2 gives b_(k-i) <= B (k/2)^(-3/2) and sums to <= B alpha (k/2)^(-3/2); i > k/2 symmetric with A beta: correct.")
    try:
        import numpy as np
        KN = min(1500, KMAX)
        sigma_full = sum(b[n] / n for n in range(1, KMAX + 1)) + 2 * MP("2.05") / mpmath.sqrt(KMAX)
        g = np.array([0.0] + [float(b[n] / n) for n in range(1, KN + 1)])
        C = 2.05
        s_up = float(sigma_full)
        conv = g.copy()
        ok_conv = True
        worst_ratio = 0.0
        exp_series = np.zeros(KN + 1)
        fact = 1.0
        for m in range(1, 41):
            if m > 1:
                conv = np.convolve(conv, g)[:KN + 1]
                fact *= m
            exp_series += conv / fact
            if m <= 6:
                for k in range(1, KN + 1):
                    bound = C * m * (2 ** 1.5 * s_up) ** (m - 1) * k ** -1.5
                    worst_ratio = max(worst_ratio, conv[k] / bound)
                    if conv[k] > bound * (1 + 1e-12):
                        ok_conv = False
        check(ok_conv, "(e) numerically (g^(*m))_k <= C m (2^(3/2) sigma)^(m-1) k^(-3/2) for m <= 6, k <= %d (max ratio to the bound %.4f)" % (KN, worst_ratio))
        wk = np.array([0.0] + [float(MP(W[k]) * mpmath.power(2, -H * k)) for k in range(1, KN + 1)])
        rel = max(abs(exp_series[k] - wk[k]) / wk[k] for k in range(1, KN + 1))
        check(rel < 1e-8, "(e) w_k = sum_(m<=40) (g^(*m))_k/m! reproduces W_k 2^(-hk) in floats for k <= %d (max relative error %.2e): the exponential form is consistent" % (KN, rel))
    except ImportError:
        print("      numpy not available: numerical convolution checks skipped")
    # (f) sigma and the constant 545
    sigma = sum(b[n] / n for n in range(1, KMAX + 1))
    tail = 2 * MP("2.05") / mpmath.sqrt(KMAX)
    const = MP("2.05") * mpmath.exp(2 * mpmath.sqrt(2) * (sigma + tail))
    tail199 = 2 * MP("1.99") / mpmath.sqrt(KMAX)
    const199 = MP("1.99") * mpmath.exp(2 * mpmath.sqrt(2) * (sigma + tail199))
    print("      sigma_%d = sum_(n<=%d) b_n/n = %s (note: 1.8989); tail bound sum_(n>%d) 2.05 n^(-3/2) <= 2*2.05/sqrt(%d) = %s (note: 0.0749; integral test);"
          % (KMAX, KMAX, mpmath.nstr(sigma, 8), KMAX, KMAX, mpmath.nstr(tail, 5)))
    print("      2.05 exp(2 sqrt2 (sigma + tail)) = %s (theorem: 545); the .out's 525.64 uses C = 1.99 (only numerically justified to n <= 3000): recomputed %s"
          % (mpmath.nstr(const, 6), mpmath.nstr(const199, 6)))
    check(const <= 545 and abs(sigma - MP("1.898923")) < MP("1e-5"), "(f) sigma = 1.898923 reproduced; 2.05 e^(2 sqrt2 (sigma + tail)) = %s <= 545" % mpmath.nstr(const, 6))
    # upper bound holds numerically
    ratio_max = MP(0)
    ratio_min = MP(1e9)
    argmin_w = None
    for k in range(1, KMAX + 1):
        v = MP(W[k]) * mpmath.power(2, -H * k) * mpmath.power(k, MP(3) / 2)
        if v > ratio_max:
            ratio_max = v
        if v < ratio_min:
            ratio_min, argmin_w = v, k
    check(ratio_max <= 545 and ratio_min >= MP("0.26"), "    W_k k^(3/2) 2^(-hk) in [%s (k = %d), %s] for 1 <= k <= %d: inside [0.26, 545]" % (mpmath.nstr(ratio_min, 5), argmin_w, mpmath.nstr(ratio_max, 6), KMAX))
    # (g) lower constant
    hp_max = mpmath.log((1 - (RHO + MP("0.1"))) / (RHO + MP("0.1")), 2)   # h'(p) = log2((1-p)/p), most negative at p = rho + 0.1
    lower_const_note = mpmath.power(2, MP("-1.4415")) / mpmath.sqrt(2)
    lower_const_true = mpmath.power(2, hp_max) / mpmath.sqrt(2)
    print("      max |h'| on [rho, rho+0.1] = |log2((1-p)/p)| at p = rho + 0.1 = %s (note and THM-4479: 1.4415, from rounding 0.7309/0.2691); 2^(-1.4415)/sqrt2 = %s (note: 0.2603); 2^(-%s)/sqrt2 = %s"
          % (mpmath.nstr(-hp_max, 7), mpmath.nstr(lower_const_note, 6), mpmath.nstr(-hp_max, 6), mpmath.nstr(lower_const_true, 6)))
    check(-hp_max > MP("1.4415") and lower_const_true >= MP("0.26"), "(g) |h'| <= 1.4415 is slightly FALSE at the right end (true max %s), but 2^(-%s)/sqrt2 = %s still >= 0.26: harmless rounding" % (mpmath.nstr(-hp_max, 6), mpmath.nstr(-hp_max, 5), mpmath.nstr(lower_const_true, 5)))
    ok_mv = True
    ok_mv_note = True
    for k in range(10, KMAX + 1):
        d = (H - hbits_mp(MP(J0[k]) / k)) * k
        if d > -hp_max:
            ok_mv = False
        if d > MP("1.4415"):
            ok_mv_note = False
    check(ok_mv, "(g) k (h - h(j0/k)) <= max|h'| for all 10 <= k <= %d (mean value theorem on [rho, rho + 1/k] subset [rho, rho+0.1]); with the note's 1.4415: %s" % (KMAX, ok_mv_note))
    check(all(0 < MP(0.25) - (MP(J0[k]) / k) * (1 - MP(J0[k]) / k) for k in range(3, KMAX + 1)), "(g) p(1-p) < 1/4 for p = j0/k (p > rho > 1/2), so sqrt(8kp(1-p)) < sqrt(2k)")
    rows = []
    cmin, bmin, nmin = MP(1e9), MP(1e9), MP(1e9)
    for k in range(1, KMAX + 1):
        f = mpmath.power(2, -H * k) * mpmath.sqrt(k)
        cv = MP(comb(k, J0[k])) * f
        bv = MP(B[k]) * f
        cmin, bmin = min(cmin, cv), min(bmin, bv)
        if k <= 400:
            Nk = necklaces_burnside(k, pow3)
            nv = MP(Nk) * f * k
            nmin = min(nmin, nv)
        if k <= 9:
            rows.append((k, J0[k], comb(k, J0[k]), B[k], Nk, mpmath.nstr(cv, 5), mpmath.nstr(bv, 5), mpmath.nstr(nv, 5)))
    print("      k < 10 (not all in the note's table): k, j0, C(k,j0), B_k, N_k, C k^.5 2^-hk, B_k k^.5 2^-hk, N_k k^1.5 2^-hk")
    for r in rows:
        print("        ", r)
    check(cmin >= MP("0.26") and bmin >= MP("0.26") and nmin >= MP("0.26"), "(g) C(k,j0) k^(1/2) 2^(-hk) >= %s, B_k k^(1/2) 2^(-hk) >= %s (k <= %d), N_k k^(3/2) 2^(-hk) >= %s (k <= 400): all >= 0.26, every k >= 1" % (mpmath.nstr(cmin, 5), mpmath.nstr(bmin, 5), KMAX, mpmath.nstr(nmin, 5)))
    ok_ck = all(MP(comb(k, J0[k])) >= mpmath.power(2, H * k + hp_max) / mpmath.sqrt(2 * k) for k in range(10, KMAX + 1))
    check(ok_ck, "(g) C(k,j0) >= 2^(hk - max|h'|)/sqrt(2k) for all 10 <= k <= %d" % KMAX)
    check(all(necklaces_burnside(k, pow3) == necklaces_direct(k) for k in range(1, 17)), "(g) own Burnside necklace count == direct enumeration for k <= 16; N_k >= C(k,j0)/k since a necklace has at most k rotations")
    check(all(necklaces_burnside(k, pow3) * k >= B[k] for k in range(1, 401)), "(g) N_k >= B_k/k for k <= 400 (so also N_k >= C(k,j0)/k)")
    # (h) oscillation window and Burnside error
    hprho = mpmath.log((1 - RHO) / RHO, 2)
    A = 1 / mpmath.sqrt(2 * mpmath.pi * RHO * (1 - RHO)) * tailconst
    print("      h'(rho) = log2((1-rho)/rho) = %s (note: -0.7731); window ends A = %s (theta -> 0) and A (1-rho)/rho = %s (theta = 1) (note: 1.992, 1.166); ratio rho/(1-rho) = %s"
          % (mpmath.nstr(hprho, 6), mpmath.nstr(A, 6), mpmath.nstr(A * q, 6), mpmath.nstr(1 / q, 6)))
    rmin, rmax = MP(1e9), MP(0)
    for k in range(max(1, KMAX - 1000), KMAX + 1):
        theta = mpmath.ceil(k * RHO) - k * RHO
        f = A * mpmath.power(q, theta)
        v = MP(B[k]) * mpmath.sqrt(k) * mpmath.power(2, -H * k)
        rmin, rmax = min(rmin, v / f), max(rmax, v / f)
    check(rmin > MP("0.95") and rmax < MP("1.01"), "(h) B_k k^(1/2) 2^(-hk) / [A ((1-rho)/rho)^theta] in [%s, %s] for %d <= k <= %d: the theta-window is the right asymptotic form (1 + o(1))" % (mpmath.nstr(rmin, 5), mpmath.nstr(rmax, 5), max(1, KMAX - 1000), KMAX))
    ok_bs = all(abs(Fraction(necklaces_burnside(k, pow3)) - Fraction(B[k], k)) <= Fraction(2 ** (k // 2 + 1)) for k in range(1, 401))
    check(ok_bs, "(h) |N_k - B_k/k| <= 2^(k/2) for k <= 400 (Burnside: the d >= 2 terms total at most (1/k) sum_(d|k) phi(d) 2^(k/2) = 2^(k/2)); the note's O(k 2^(k/2)) is weaker and correct")
    print("      note: the window [%s, %s] is asymptotic; the observed B_k k^0.5 2^-hk minimum 1.1355 (k <= 3000) lies below A q = %s by the finite-k correction (ratio bound (n-j)/(j+1) < q strict)." % (mpmath.nstr(A * q, 5), mpmath.nstr(A, 5), mpmath.nstr(A * q, 5)))
    return b


# ================================================================ 8. Corollary C3
def T_step(x, b):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


def section8(W, TFULL, TCOUNT):
    print("== 8. Corollary C3: carry bound, four cases, brute-force dip counts ==")
    # carry: T_b^j(n) = (3^(o_j) n + b c_j)/2^j with c_j = sum_(i<=j, w_i=1) 3^(o_j-o_i) 2^(i-1); c_j <= 3^j - 2^j (equality iff w = 1^j)
    ok_formula = True
    for bsh in (1, -1):
        for r in range(0, 2 ** 12):
            x = r
            o = 0
            c = 0
            for j in range(1, 13):
                if x % 2 == 1:
                    o += 1
                    c = 3 * c + 2 ** (j - 1)
                x = T_step(x, bsh)
                if (3 ** o * r + bsh * c) % 2 ** j != 0 or (3 ** o * r + bsh * c) // 2 ** j != x:
                    ok_formula = False
    check(ok_formula, "T_b^j(r) == (3^(o_j) r + b c_j)/2^j exactly for all r < 2^12, j <= 12, both sheets (c_j = 3 c_(j-1) + 2^(j-1) at odd steps): sign(beta_j) = b")
    ok_carry = True
    for n in range(1, 17):
        for w in all_words(n):
            c = 0
            for j, x in enumerate(w, 1):
                if x:
                    c = 3 * c + 2 ** (j - 1)
                if c > 3 ** j - 2 ** j or (c == 3 ** j - 2 ** j) != all(w[:j]):
                    ok_carry = False
                if (c == 0) != (sum(w[:j]) == 0):
                    ok_carry = False
    check(ok_carry, "|beta_j| = c_j/2^j <= (3/2)^j - 1 for all words of length <= 16, equality iff the first j letters are all 1; beta_j = 0 iff no odd step (induction: c_j <= 3(3^(j-1) - 2^(j-1)) + 2^(j-1))")
    # the four inequalities
    xs = [i / 10000 for i in range(0, 5001)]
    check(all(math.log2(1 - x) >= -2 * x - 1e-12 for x in xs[:-1]) and math.log2(0.5) == -1.0, "plus/upper: log2(1-x) >= -2x on [0, 1/2] (concave, chord through (0,0),(1/2,-1)); (3/4)^t <= 1/2 iff t >= 3; so S_j > -2(3/4)^t >= -0.844 > -1 for t >= 3")
    check(all(2 ** (t - 1) >= (Fraction(3, 2)) ** t for t in range(3, 60)) and 2 ** 1 < Fraction(9, 4), "minus/lower: 2^(t-1) >= (3/2)^t iff t >= 3 (fails at t = 2); then n >= 2^t gives n/2 >= (3/2)^t, i.e. (3/2)n - (3/2)^t >= n")
    print("      plus/upper chain: (1 - M_j) n <= beta_j <= (3/2)^j - 1 < (3/2)^t and n >= 2^t give 1 - 2^(S_j) < (3/4)^t, S_j > log2(1 - (3/4)^t) >= -2(3/4)^t;")
    print("      11w then has partial sums 0.585, 1.170, 1.170 + S_j(w) > 0.17: positive; w -> 11w injective: correct.")
    print("      minus/upper: beta_j <= 0 so T_-^j(n) <= M_j n; non-dipper forces M_j >= 1, hence M_j > 1 (no ties): w positive. correct.")
    print("      minus/lower: 11w with w positive has M_1 = 3/2, M_2 = 9/4, M_j = (9/4) 2^(S_(j-2)(w)) > 9/4 (j > 2), so M_j >= 3/2 for all j >= 1 and T_-^j(n) >= (3/2) n - ((3/2)^t - 1) > n for n >= 2^t, t >= 3. correct.")
    # the smallest |S_j| with S_j < 0, j <= t, against 2 (3/4)^t (when the plus-sheet argument already forces positivity)
    dmin = MP(10)
    first_t = None
    for t in range(1, 3001):
        o = int(mpmath.floor(t * RHO))
        d = t - o * LOG2_3      # the negative partial sum closest to 0 at step t
        dmin = min(dmin, d)
        if first_t is None and t >= 3 and 2 * MP(0.75) ** t < dmin:
            first_t = t
        if t >= 12 and not (2 * MP(0.75) ** t < dmin):
            first_t = -1
    print("      remark: min over j <= t of the negative S_j closest to 0 exceeds 2(3/4)^t from t = %d on (holds for all 12 <= t <= 3000), so on the plus sheet the note's own argument" % first_t)
    print("      makes every non-dipper in [2^t, 2^(t+1)) a positive word for 12 <= t <= 3000; the exact equality D_+ = sum W_t is then a finite check below t = 12.")
    # brute force D_b(2^T, 1) with full per-block classification (T < TFULL), definition T_b^j(n) >= n for 0 <= j <= floor(log2 n)
    pos_tab = {}
    for j in range(0, TFULL + 3):
        for o in range(0, j + 1):
            pos_tab[(o, j)] = 3 ** o > 2 ** j
    print("      per dyadic block [2^t, 2^(t+1)), t < %d: t, W_t | plus: non-dippers, of which non-positive word | minus: non-dippers, positive words that dip, 11w-classes (w positive) that are non-dippers, W_(t-2)" % TFULL)
    ok_plus = ok_minus = ok_11w = True
    cum = {1: [0] * (TFULL + 1), -1: [0] * (TFULL + 1)}
    for t in range(0, TFULL):
        stats = {}
        for bsh in (1, -1):
            nd = 0
            nd_nonpos = 0
            pos_dip = 0
            w11_nd = 0
            w11_tot = 0
            min_S_nd = 10.0
            for n in range(2 ** t, 2 ** (t + 1)):
                x = n
                o = 0
                positive = True
                nondip = True
                w1 = w2 = 0
                tailpos = True
                minS = 0.0
                for j in range(1, t + 1):
                    if x % 2 == 1:
                        o += 1
                        if j == 1:
                            w1 = 1
                        if j == 2:
                            w2 = 1
                    x = T_step(x, bsh)
                    if x < n:
                        nondip = False
                    if not pos_tab[(o, j)]:
                        positive = False
                    if j >= 3 and w1 and w2 and not pos_tab[(o - 2, j - 2)]:
                        tailpos = False
                    s = o * float(LOG2_3) - j
                    if s < minS:
                        minS = s
                is11w = (t >= 2 and w1 == 1 and w2 == 1 and tailpos)
                if nondip:
                    nd += 1
                    if not positive:
                        nd_nonpos += 1
                    if minS < min_S_nd:
                        min_S_nd = minS
                if positive and not nondip:
                    pos_dip += 1
                if is11w:
                    w11_tot += 1
                    if nondip:
                        w11_nd += 1
            stats[bsh] = (nd, nd_nonpos, pos_dip, w11_nd, w11_tot, min_S_nd)
            cum[bsh][t] = nd
        Wt = W[t]
        p, m = stats[1], stats[-1]
        print("        t=%2d W_t=%7d | plus: %7d nd, %d non-positive | minus: %7d nd, %d positive dippers, 11w nd %d of %d (W_(t-2) = %s); min S_j over plus non-dippers %.4f vs -2(3/4)^t = %.4f"
              % (t, Wt, p[0], p[1], m[0], m[2], m[3], m[4], W[t - 2] if t >= 2 else "-", p[5], -2 * 0.75 ** t))
        if t >= 1 and p[0] != Wt:
            ok_plus = False
        if t >= 1 and m[0] != Wt:
            ok_minus = False
        if t >= 3 and (m[3] != m[4] or m[4] != W[t - 2]):
            ok_11w = False
        if t >= 3 and p[5] <= -2 * 0.75 ** t:
            ok_plus = False
    check(ok_plus, "plus sheet: #non-dippers in [2^t, 2^(t+1)) == W_t for 1 <= t < %d (no non-positive non-dipper occurs), and all their S_j > -2(3/4)^t" % TFULL)
    check(ok_minus, "minus sheet: #non-dippers in [2^t, 2^(t+1)) == W_t for 1 <= t < %d (every positive word is a non-dipper, none dips)" % TFULL)
    check(ok_11w, "minus sheet: every class 11w with w positive (length t-2) is a non-dipper for 3 <= t < %d, and there are exactly W_(t-2) of them" % TFULL)
    # cumulative counts vs sum_(t<T) W_t, including the n = 1 term (t = 0, W_0 = 1)
    for bsh in (1, -1):
        for T in (8, 12, 16, min(20, TFULL)):
            if T > TFULL:
                continue
            D_all = sum(cum[bsh][t] for t in range(0, T))      # n = 1 .. 2^T - 1 (n = 2^T is a dipper: T(2^T) = 2^(T-1))
            D_ge2 = D_all - cum[bsh][0]
            S1 = sum(W[t] for t in range(1, T))
            S0 = S1 + W[0]
            print("      sheet %+d, T = %d: D_b(2^T - 1, 1) over n >= 1 = %d, over n >= 2 = %d; sum_(1<=t<T) W_t = %d, sum_(0<=t<T) W_t = %d" % (bsh, T, D_all, D_ge2, S1, S0))
            check(D_ge2 == S1 and D_all == S0, "      D_b(2^%d, 1) (n >= 2, THM-4487's script convention) == sum_(1<=t<%d) W_t = %d; with n = 1 it is %d = sum_(0<=t<%d) W_t (the theorem's 'sum_(t<T) W_t' omits W_0 = 1)" % (T, T, S1, S0, T))
    # count-only brute force to TCOUNT (odd n only; early break), both sheets
    if TCOUNT > TFULL:
        for bsh in (1, -1):
            cnt = 1   # n = 1
            for n in range(3, 2 ** TCOUNT, 2):
                x = n
                k = n.bit_length() - 1
                ok = True
                for _ in range(k):
                    x = T_step(x, bsh)
                    if x < n:
                        ok = False
                        break
                if ok:
                    cnt += 1
            S1 = sum(W[t] for t in range(1, TCOUNT))
            check(cnt - 1 == S1, "      sheet %+d: brute-force D_b(2^%d, 1) over n >= 2 = %d == sum_(1<=t<%d) W_t = %d (THM-4487's .out: %s)" % (bsh, TCOUNT, cnt - 1, TCOUNT, S1, {12: 281, 16: 2903, 20: 31730, 24: 367698}.get(TCOUNT, "n/a")))
    print("      sum_(1<=t<T) W_t for T = 8, 12, 16, 20, 24:", [sum(W[t] for t in range(1, T)) for T in (8, 12, 16, 20, 24)], "(THM-4487's gamma = 1 column: 281, 2903, 31730, 367698 on both sheets, n in [2, 2^T])")


# ================================================================ 9. quotes
def grep_lines(path, needles, maxlines=6):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().split("\n")
    except OSError as e:
        print("      (cannot open %s: %s)" % (path, e))
        return
    shown = 0
    for i, line in enumerate(lines, 1):
        if any(nd in line for nd in needles):
            print("      %s:%d: %s" % (os.path.basename(path), i, line.strip()))
            shown += 1
            if shown >= maxlines:
                break
    if shown == 0:
        print("      (no line of %s contains %s)" % (os.path.basename(path), needles))


def section9(root):
    print("== 9. Corollaries C1/C2: the sandwich lines as stated upstream ==")
    grep_lines(os.path.join(root, "01-canon", "theorems", "THM-4479-strategy-cube-distance-to-provability-sharp.md"),
               ["delta_k <= |Bad_k| <= 2^(hk)", "delta_k >= N_k >= 2^(hk)/(3k^2)", "Bad_k = {r mod 2^k : M_j(r) > 1"], 6)
    grep_lines(os.path.join(root, "01-canon", "theorems", "THM-4485-periodic-edit-price-feedback-sets.md"),
               ["(4) Chain: delta_k (sign flips) >= FVS^odd >= FVS >= nu (cycle packing)", ">= N (expanding necklaces)."], 3)
    grep_lines(os.path.join(root, "01-canon", "theorems", "THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md"),
               ["D_b(X, gamma) = #{n <= X : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n)}", "the script omits n = 1"], 3)
    grep_lines(os.path.join(root, "04-computation", "experiments", "collatz_dipspectrum_20260926.py"),
               ["for n in range(2, Xmax + 1):"], 2)
    grep_lines(os.path.join(root, "05-knowledge", "results", "collatz_dipspectrum_20260926.out"),
               ["X=2^24 D:", "X=2^12 D:", "X=2^16 D:", "X=2^20 D:"], 8)
    grep_lines(os.path.join(root, "04-computation", "experiments", "procgen_cubedist_20260925_lib.py"),
               ["mul ** a > 2 ** j", "if not mul ** a > 2 ** j"], 3)
    grep_lines(os.path.join(root, "04-computation", "experiments", "procgen_cubedist_20260925_orchestrator_check.py"),
               ["if 3 ** a2 > 2 ** j:", "BAD_TABLE = {"], 3)
    print("      read: THM-4479 states delta_k <= |Bad_k| (Theorem 1(d)) and delta_k >= N_k (Theorem 2), both for k >= 2; THM-4485 states the chain")
    print("      delta_k >= FVS^odd >= FVS >= nu >= N (item (4)); THM-4487 defines D_b with n <= X and 0 <= i <= floor(log_2 n), its script counts n in [2, 2^T].")


# ================================================================ 10. provenance
def sha256_of(path, lf=False):
    with open(path, "rb") as f:
        data = f.read()
    if lf:
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest(), data.count(b"\r\n")


def section10(root):
    print("== 10. provenance ==")
    exp_script = "9e3d674295ae19206c7ff5a3c9e053b372e9dc25273c4e8015f787299a1b148a"
    exp_out = "9f70677d1971708c6952ea21ada16abfb67a896cc243ba4e83104e06fa45a7e1"
    sp = os.path.join(root, "04-computation", "experiments", "collatz_nodescent_order_20260926.py")
    op = os.path.join(root, "04-computation", "experiments", "collatz_nodescent_order_20260926.out")
    for path, exp, name in ((sp, exp_script, "script"), (op, exp_out, "output")):
        try:
            raw, crlf = sha256_of(path)
            lfh, _ = sha256_of(path, lf=True)
            print("      %s: raw sha256 %s, LF-normalised %s, CRLF count %d; header claims %s" % (name, raw, lfh, crlf, exp))
            check(lfh == exp, "      %s sha256 (LF bytes) matches the theorem header" % name)
        except OSError as e:
            check(False, "      %s not readable: %s" % (name, e))
    me = os.path.abspath(__file__)
    raw, crlf = sha256_of(me)
    print("      this audit script: raw sha256 %s (CRLF count %d)" % (raw, crlf))


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    TFULL = int(sys.argv[2]) if len(sys.argv) > 2 else 18
    TCOUNT = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    print("collatz_nodescent_order_20260926_audit.py: independent audit of THM-4495 (KMAX = %d, TFULL = %d, TCOUNT = %d)" % (KMAX, TFULL, TCOUNT))
    print("rho = log_3 2 = %s, h(rho) = %s, log_2 3 = %s (mpmath, 50 digits)" % (mpmath.nstr(RHO, 15), mpmath.nstr(H, 15), mpmath.nstr(LOG2_3, 15)))
    section1()
    POS, PM = section2(12)
    section3(PM, 12)
    E, logW = section4(14)
    section5(E, logW, 12)
    W, B, J0 = section6(KMAX)
    section7(W, B, J0, KMAX)
    section8(W, TFULL, TCOUNT)
    section9(root)
    section10(root)
    print("== summary ==")
    if FAILS:
        print("  %d FAILED check(s):" % len(FAILS))
        for f in FAILS:
            print("    - " + f)
    else:
        print("  all checks passed")
    print("  elapsed %.1f s" % (time.time() - T0))


if __name__ == "__main__":
    main()
