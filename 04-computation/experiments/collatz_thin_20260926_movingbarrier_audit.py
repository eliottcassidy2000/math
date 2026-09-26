#!/usr/bin/env python3
"""collatz_thin_20260926_movingbarrier_audit.py -- independent adversarial audit of THM-4499
(thin divergence is o(X^(h*)): N(X) <= K X^(h*) (log_2 X)^a for every a > lambda*/h* - 3/2).

Audited: 01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md and
         05-knowledge/results/collatz_thin_20260926_little_o_thin_divergence.md
         (script 04-computation/experiments/collatz_thin_20260926_movingbarrier.py, output .out).
Prerequisites read: THM-4495's note (Theorem A, tilt identity, convolution lemma, W_k <= 545 2^(hk) k^(-3/2))
and THM-4476's note, sections 1.1, 1.3-1.5, 1.7 (recursion (R), counting lemma, Corollary 4).

Written blind to the audited script's code paths.  Every comparison of a partial sum
S_i = o log_2 3 - i with a rational barrier -y (y = p/q) is done EXACTLY: for o >= 1,
o log_2 3 - i > -p/q  <=>  floor(q o log_2 3) >= q i - p, with floor(m log_2 3) = bitlength(3^m) - 1
(3^m is never a power of 2 for m >= 1); for o = 0 the value -i is an integer and is compared exactly
with -p/q.  Comparisons of two partial sums use the sign of 3^o 2^(i2) - 3^(o2) 2^i.  Never a float.

Sections (numbered as in the audit brief):
   1  Lemma M, decomposition (D): bijection re-derived word by word (k <= 14), reversal algebra, uniqueness,
      boundary cases, positive words counted once; (D) exact by brute force for k <= 14, integer and
      half-integer y
   2  weighted Spitzer identity for NEGATIVE words: every step of the negated proof checked word by word
      (n <= 12), rotation step needs only weight-invariance under rotation; identity checked as an exact
      bivariate-polynomial identity (n <= 48, both signs) and with Fractions for four weight pairs (n <= 30)
   3  the tilt identity Q(w) = 2^(-hm) 2^(lambda* S_m(w)) (symbolic + 60-digit numerics)
   4  section 1.2: both inequalities, all signs, checked exactly/numerically for m <= 60
   5  sections 1.3-1.4: weights, spacing, local binomial bound to n = 5000, convolution lemma and the
      induction, p_n n^(3/2) bounded for s in {0.05, 0.1, 0.5}, n <= 2000 (weighted DP), recurrence
   6  section 1.5 assembly: convolution constant 14.78, m = k term, absorption, D_s at s = 0.1, the bound
      on the script's table
   7  Lemma 1.4c: carry bound (exact Fractions), the barrier for EVERY i <= k, strictness, direct enumeration
      of F_b(2^18, theta) on both sheets against 2 M_k(y), the numerical constants
   8  section 3 bootstrap: (R''), the exponent algebra, s, X_2, the induction bracket, well-foundedness,
      transfer to Corollary 4 (text of 1.7 quoted)
   9  consistency with THM-4487(2): concavity h(rho) <= h* + lambda* theta, theta in (0, 0.1]
  10  provenance: sha256 vs header, the float-margin barrier test (continued fraction of log_2 3, min gap),
      the (B) table and the (A) P_1..P_7 lists reproduced independently, W_1..W_20 of THM-4495
Usage: python3 collatz_thin_20260926_movingbarrier_audit.py   (no arguments; deterministic output)
"""
import hashlib
import math
import operator
import os
import re
import sys
import time
from fractions import Fraction
from functools import cmp_to_key
from itertools import product
from math import comb

import mpmath

mpmath.mp.dps = 60
MP = mpmath.mpf
ALPHA = mpmath.log(3) / mpmath.log(2)          # log_2 3
RHO = 1 / ALPHA                                 # rho* = log_3 2
H = -(RHO * mpmath.log(RHO, 2) + (1 - RHO) * mpmath.log(1 - RHO, 2))
LAM = mpmath.log(RHO / (1 - RHO), 2) / ALPHA    # lambda*
ASTAR = LAM / H - MP(3) / 2
LN2 = mpmath.log(2)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))

FAILS = []
NCHECK = [0]
T0 = time.time()


def check(cond, msg):
    NCHECK[0] += 1
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAILS.append(msg)
    print("  [%s] %s" % (tag, msg))
    return cond


def note(msg):
    print("  " + msg)


def elapsed(label):
    sys.stderr.write("   ... %s: %.1f s\n" % (label, time.time() - T0))


# ---------------------------------------------------------------- exact arithmetic on the walk
_POW3 = {0: 1}


def pow3(m):
    v = _POW3.get(m)
    if v is None:
        v = 3 ** m
        _POW3[m] = v
    return v


def floor_alpha(m):
    """floor(m log_2 3) exactly, m >= 0 (3^m is not a power of two for m >= 1)."""
    if m == 0:
        return 0
    return pow3(m).bit_length() - 1


def above(o, i, p, q):
    """exact test  o log_2 3 - i > -p/q  (barrier y = p/q >= 0, q >= 1)."""
    if o == 0:
        return q * i < p
    return floor_alpha(q * o) >= q * i - p


def jmax(o, p, q):
    """largest i with o log_2 3 - i > -p/q (all smaller i also satisfy it); -1 if none for o = 0."""
    if o == 0:
        return (p - 1) // q if p > 0 else -1
    return (floor_alpha(q * o) + p) // q


def is_negative_state(o, i):
    """S_i = o log_2 3 - i < 0 for i >= 1 (exact; no ties)."""
    if i == 0:
        return False
    if o == 0:
        return True
    return i >= floor_alpha(o) + 1


def cmp_S(o, i, o2, i2):
    """sign of (o log_2 3 - i) - (o2 log_2 3 - i2), exact, any integers o, i (also negative)."""
    om = min(o, o2)
    im = min(i, i2)
    lhs = pow3(o - om) << (i2 - im)
    rhs = pow3(o2 - om) << (i - im)
    return (lhs > rhs) - (lhs < rhs)


def build_ranks(nmax=64):
    """an exact total order of the values o log_2 3 - i on [0, nmax]^2 (ties only at equal pairs)."""
    keys = [(o, i) for o in range(nmax + 1) for i in range(nmax + 1)]
    keys.sort(key=cmp_to_key(lambda a, b: cmp_S(a[0], a[1], b[0], b[1])))
    rank = {}
    r = 0
    prev = None
    for kk in keys:
        if prev is not None and cmp_S(prev[0], prev[1], kk[0], kk[1]) == 0:
            rank[kk] = rank[prev]
        else:
            rank[kk] = r
        r += 1
        prev = kk
    return rank


RANK = build_ranks(64)
ZERO = RANK[(0, 0)]


def psums(w):
    """partial sums of word w as (o_i, i), i = 0..len(w)."""
    out = [(0, 0)]
    o = 0
    for i, c in enumerate(w, 1):
        o += c
        out.append((o, i))
    return out


def ranks(w):
    return [RANK[p] for p in psums(w)]


# ---------------------------------------------------------------- DP engines (exact, list based)
def barrier_dp(k, p, q):
    """M_j(y), y = p/q, for all j <= k, plus the final state vector (index o) at step k.
    State o alive at step j iff j <= jmax(o); the alive set is {o >= lo_j} (jmax increasing in o)."""
    cur = [1]
    Ms = [1]
    jm = {}
    for j in range(1, k + 1):
        nxt = list(map(operator.add, cur + [0], [0] + cur))
        for o in range(j + 1):
            if o not in jm:
                jm[o] = jmax(o, p, q)
            if j <= jm[o]:
                break
            nxt[o] = 0
        cur = nxt
        Ms.append(sum(cur))
    return Ms, cur


def negative_dp(n):
    """counts c[j][o] of words of length j with ALL partial sums < 0 and o ones, j <= n (exact)."""
    cur = [1]
    rows = [cur]
    for j in range(1, n + 1):
        nxt = list(map(operator.add, cur + [0], [0] + cur))
        for o in range(j, -1, -1):
            if is_negative_state(o, j):
                break
            nxt[o] = 0
        cur = nxt
        rows.append(cur)
    return rows


def positive_dp(n):
    """counts c[j][o] of words of length j with ALL partial sums > 0 and o ones, j <= n (exact)."""
    cur = [1]
    rows = [cur]
    for j in range(1, n + 1):
        nxt = list(map(operator.add, cur + [0], [0] + cur))
        for o in range(j + 1):
            if above(o, j, 0, 1):
                break
            nxt[o] = 0
        cur = nxt
        rows.append(cur)
    return rows


def N_m_y(negrows, m, p, q):
    """N_m(y): negative words of length m with total > -y (exact)."""
    return sum(c for o, c in enumerate(negrows[m]) if c and above(o, m, p, q))


# ====================================================================== section 1
def section1():
    print("== 1. Lemma M: the decomposition (D) at the minimum ==")
    # sanity of the exact machinery
    ok = all(cmp_S(o, i, o2, i2) != 0 for o in range(0, 40) for i in range(0, 40)
             for (o2, i2) in [(o + 1, i + 1), (o + 1, i + 2), (o + 2, i + 3), (o + 5, i + 8), (o + 12, i + 19)])
    check(ok, "exact comparator: S(o,i) != S(o+d,i+e) for the convergent shifts (d,e) in {(1,1),(1,2),(2,3),(5,8),(12,19)}, o,i < 40")
    ok = True
    for k in range(1, 15):
        for w in product((0, 1), repeat=k):
            r = ranks(w)
            if len(set(r[1:])) != k or ZERO in r[1:]:
                ok = False
    check(ok, "no ties: all S_1..S_k of every word of length k <= 14 are distinct and nonzero (minima unique)")
    # reversal formula
    ok = True
    for m in range(1, 11):
        for a in product((0, 1), repeat=m):
            pa = psums(a)
            pr = psums(a[::-1])
            for i in range(m + 1):
                # S_i(a^R) = S_m(a) - S_(m-i)(a): as (o, i) pairs
                if pr[i] != (pa[m][0] - pa[m - i][0], pa[m][1] - pa[m - i][1]):
                    ok = False
    check(ok, "reversal algebra S_i(a^R) = S_m(a) - S_(m-i)(a) exactly (as pairs (o, i)) for all words of length <= 10")

    ys = [(0, 1), (1, 2), (1, 1), (3, 2), (2, 1), (5, 2), (3, 1), (7, 2), (4, 1), (13, 2), (8, 1)]
    note("barriers y tested: " + ", ".join(str(Fraction(p, q)) for p, q in ys))
    KB = 14
    # brute-force objects for k <= KB
    words = {k: list(product((0, 1), repeat=k)) for k in range(0, KB + 1)}
    pos = {k: [w for w in words[k] if all(rr > ZERO for rr in ranks(w)[1:])] for k in range(0, KB + 1)}
    W = [len(pos[k]) for k in range(KB + 1)]
    neg = {k: [w for w in words[k] if all(rr < ZERO for rr in ranks(w)[1:])] for k in range(0, KB + 1)}
    check(W[1:] == [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734],
          "brute-force W_k, k <= 14, equals THM-4495's list 1,1,2,3,4,8,13,19,38,64,128,226,367,734")
    # positive min-ending words (THM-4495's P_m = L_m) and first-passage words
    def min_ending(w):
        r = ranks(w)
        m = len(w)
        return all(r[j] > r[m] for j in range(1, m))
    PM = [len([w for w in pos[m] if min_ending(w)]) for m in range(KB + 1)]
    L = [len([w for w in words[m] if all(rr < ZERO for rr in ranks(w)[1:m]) and ranks(w)[m] > ZERO]) for m in range(KB + 1)]
    check(PM[1:] == L[1:], "positive min-ending words P_m == first-passage words L_m (m <= 14): " + str(PM[1:]))
    check(all(W[k] == sum(PM[m] * W[k - m] for m in range(1, k + 1)) for k in range(1, KB + 1)),
          "THM-4495 Step 1: W_k = sum_m L_m W_(k-m) for k <= 14 (the W_k term of (D))")

    allD = True
    allbij = True
    allonce = True
    for (p, q) in ys:
        Nm = [0] + [len([w for w in neg[m] if above(sum(w), m, p, q)]) for m in range(1, KB + 1)]
        negset = {m: set(neg[m]) for m in range(KB + 1)}
        negy = {m: set(w for w in neg[m] if above(sum(w), m, p, q)) for m in range(KB + 1)}
        posset = {m: set(pos[m]) for m in range(KB + 1)}
        for k in range(1, KB + 1):
            Mk = [w for w in words[k] if all(above(o, i, p, q) for (o, i) in psums(w)[1:])]
            rhs = W[k] + sum(Nm[m] * W[k - m] for m in range(1, k + 1))
            if len(Mk) != rhs:
                allD = False
            # word-by-word bijection
            pairs_pos = []
            pairs_neg = []
            for w in Mk:
                r = ranks(w)
                m = min(range(1, k + 1), key=lambda j: r[j])      # unique argmin
                a, b = w[:m], w[m:]
                if b not in posset[k - m]:
                    allbij = False
                if r[m] > ZERO:
                    # positive min-ending a, whole word positive
                    if not (a in posset[m] and min_ending(a) and w in posset[k]):
                        allbij = False
                    pairs_pos.append((a, b))
                else:
                    ar = a[::-1]
                    if ar not in negy[m]:
                        allbij = False
                    # a itself ends at its strict minimum, below S_0 = 0, with S_m > -y
                    if not (all(r[j] > r[m] for j in range(1, m)) and r[m] < ZERO and above(sum(a), m, p, q)):
                        allbij = False
                    pairs_neg.append((ar, b))
            # injectivity and surjectivity onto the product sets
            if len(set(pairs_pos)) != len(pairs_pos) or len(set(pairs_neg)) != len(pairs_neg):
                allbij = False
            if len(pairs_pos) != W[k]:
                allonce = False
            if set(pairs_neg) != set((c, b) for m in range(1, k + 1) for c in negy[m] for b in posset[k - m]):
                allbij = False
            # gluing back (the converse): every (c, b) glues to a word of M_k(y) with argmin at |c|
            for m in range(1, k + 1):
                for c in negy[m]:
                    for b in pos[k - m]:
                        w = c[::-1] + b
                        r = ranks(w)
                        if not (all(above(o, i, p, q) for (o, i) in psums(w)[1:]) and min(range(1, k + 1), key=lambda j: r[j]) == m):
                            allbij = False
    check(allD, "(D) M_k(y) = W_k + sum_(m=1)^k N_m(y) W_(k-m) holds EXACTLY by brute force for k <= 14 and all 11 barriers")
    check(allbij, "the two-case decomposition at the argmin is a bijection onto {pos min-ending} x {pos} + sum_m {neg, total > -y} x {pos}; gluing is its inverse (argmin at |a|)")
    check(allonce, "words with all prefix sums positive are exactly the case S_m(a) > 0 and are counted once, by the W_k term")
    # boundary cases: m = k (b empty, W_0 = 1) and y = 0
    M0 = [len([w for w in words[k] if all(above(o, i, 0, 1) for (o, i) in psums(w)[1:])]) for k in range(KB + 1)]
    check(M0 == W, "y = 0: M_k(0) = W_k and N_m(0) = 0 (no negative word has total > 0)")
    Nk_half = [len([w for w in neg[k] if above(sum(w), k, 1, 2)]) for k in range(KB + 1)]
    check(all(Nk_half[k] == len([w for w in words[k] if all(above(o, i, 1, 2) for (o, i) in psums(w)[1:]) and ranks(w)[k] == min(ranks(w)[1:]) and ranks(w)[k] < ZERO]) for k in range(1, KB + 1)),
          "m = k term (b empty, W_0 = 1): words of M_k(1/2) whose minimum is at k and negative are in bijection with N_k(1/2)")
    note("integer barriers: the all-zero prefix of length y has S_y = -y exactly; the strict test S_i > -y excludes it, and (D) transports")
    note("the strict condition verbatim (checked above at y = 1, 2, 3, 4, 8 with the exact integer comparison), so no tie-breaking is needed")
    elapsed("section 1")


# ====================================================================== section 2
def poly_mul(A, B):
    """product of polynomials in x1 (dict o -> coeff), with the t-degree tracked outside."""
    out = {}
    for a, ca in A.items():
        for b, cb in B.items():
            out[a + b] = out.get(a + b, 0) + ca * cb
    return out


def section2():
    print("== 2. The weighted Spitzer identity for negative words (the negated proof, step by step) ==")
    NB = 12
    words = {k: list(product((0, 1), repeat=k)) for k in range(0, NB + 1)}

    def is_neg(w):
        return all(rr < ZERO for rr in ranks(w)[1:])

    def max_ending_neg(w):
        r = ranks(w)
        m = len(w)
        return is_neg(w) and all(r[j] < r[m] for j in range(1, m))

    def fpb(w):          # first-passage below: S_i > 0 for i < m, S_m < 0
        r = ranks(w)
        m = len(w)
        return m >= 1 and all(r[i] > ZERO for i in range(1, m)) and r[m] < ZERO

    def ends_at_min(w):  # S_n < S_j for all 0 <= j < n
        r = ranks(w)
        n = len(w)
        return n >= 1 and all(r[n] < r[j] for j in range(0, n))

    def desc_ladder(w):  # j in [1, n] with S_j < S_l for all 0 <= l < j
        r = ranks(w)
        return [j for j in range(1, len(w) + 1) if all(r[j] < r[l] for l in range(0, j))]

    NEG = {k: [w for w in words[k] if is_neg(w)] for k in range(NB + 1)}
    NEGset = {k: set(NEG[k]) for k in range(NB + 1)}
    MEN = {k: [w for w in NEG[k] if max_ending_neg(w)] for k in range(NB + 1)}
    FPB = {k: [w for w in words[k] if fpb(w)] for k in range(NB + 1)}
    Pn = [len(NEG[k]) for k in range(NB + 1)]
    note("P_n^- (negative words, unit weights), n = 1..12: " + str(Pn[1:]))
    # Step 1^-: decomposition at the maximum
    ok = True
    for n in range(1, NB + 1):
        pairs = []
        for w in NEG[n]:
            r = ranks(w)
            m = max(range(1, n + 1), key=lambda j: r[j])
            a, b = w[:m], w[m:]
            if not (max_ending_neg(a) and (b in NEGset[n - m])):
                ok = False
            pairs.append((a, b))
        target = set((a, b) for m in range(1, n + 1) for a in MEN[m] for b in NEG[n - m])
        if set(pairs) != target or len(set(pairs)) != len(pairs):
            ok = False
        for m in range(1, n + 1):
            for a in MEN[m]:
                for b in NEG[n - m]:
                    w = a + b
                    r = ranks(w)
                    if not (is_neg(w) and max(range(1, n + 1), key=lambda j: r[j]) == m):
                        ok = False
    check(ok, "Step 1^-: a negative word splits at its unique argmax into a max-ending negative word and a negative word; bijection with inverse gluing (n <= 12)")
    check(all(Pn[n] == sum(len(MEN[m]) * Pn[n - m] for m in range(1, n + 1)) for n in range(1, NB + 1)),
          "P_n^- = sum_m PM_m^- P_(n-m)^-  (n <= 12), PM_m^- = max-ending negative words: " + str([len(MEN[m]) for m in range(1, NB + 1)]))
    # Step 2^-: reversal
    check(all(set(a[::-1] for a in MEN[m]) == set(FPB[m]) for m in range(1, NB + 1)),
          "Step 2^-: reversal maps max-ending negative words of length m ONTO first-passage-below words (S_i > 0 for i < m, S_m < 0), m <= 12")
    # letter multisets preserved by reversal / concatenation (weights)
    check(all(sum(a) == sum(a[::-1]) for m in range(1, NB + 1) for a in MEN[m]), "reversal preserves the letter multiset, hence the weight x(w) = x_0^(#0) x_1^(#1)")
    # Step 3^-: ladder blocks, as polynomials in x1 (coefficient of x1^o), for each n, r
    Lpoly = {m: {} for m in range(NB + 1)}
    for m in range(1, NB + 1):
        for w in FPB[m]:
            o = sum(w)
            Lpoly[m][o] = Lpoly[m].get(o, 0) + 1
    ok = True
    for n in range(1, NB + 1):
        E = {}
        for w in words[n]:
            if ends_at_min(w):
                lad = desc_ladder(w)
                if lad[-1] != n:
                    ok = False
                # blocks between consecutive epochs must be first-passage-below
                prev = 0
                for j in lad:
                    if not fpb(w[prev:j]):
                        ok = False
                    prev = j
                key = (len(lad), sum(w))
                E[key] = E.get(key, 0) + 1
            else:
                lad = desc_ladder(w)
                if lad and lad[-1] == n:
                    ok = False
        # [t^n] L(t)^r as polynomial in x1: sum over compositions
        for r in range(1, n + 1):
            # coefficient of t^n in L^r: convolution over lengths
            acc = {}
            def rec(rem, parts, polyacc):
                nonlocal acc
                if parts == 0:
                    if rem == 0:
                        for o, c in polyacc.items():
                            acc[o] = acc.get(o, 0) + c
                    return
                for m in range(1, rem - parts + 2):
                    if Lpoly[m]:
                        rec(rem - m, parts - 1, poly_mul(polyacc, Lpoly[m]))
            rec(n, r, {0: 1})
            for o in range(n + 1):
                if acc.get(o, 0) != E.get((r, o), 0):
                    ok = False
    check(ok, "Step 3^-: a word ends at its minimum iff it is a concatenation of first-passage-below words (cut at descending ladder epochs); E^-_(n,r) = [t^n] L^-(t)^r as polynomials in x_1 (n <= 12)")
    # Step 4^-: rotation averaging with strict-minimum records, per word, then weighted regrouping
    ok = True
    okwin = True
    Bpoly_from_E = {}
    examined = 0
    for n in range(1, NB + 1):
        Erow = {}
        for w in words[n]:
            if ends_at_min(w):
                key = (len(desc_ladder(w)), sum(w))
                Erow[key] = Erow.get(key, 0) + 1
        for x in words[n]:
            o_n = sum(x)
            if not is_negative_state(o_n, n):
                continue
            examined += 1
            ps = psums(x)

            def S(l):  # periodic extension as (o, i) pairs
                t, rmd = divmod(l, n)
                return (ps[rmd][0] + t * o_n, rmd + t * n)

            def is_record(j, win):
                oj, ij = S(j)
                return all(cmp_S(oj, ij, *S(l)) < 0 for l in range(j - win, j))
            rec_res = set()
            for j in range(0, n):
                r1 = is_record(j, n)
                r3 = is_record(j, 3 * n)
                if r1 != r3:
                    okwin = False
                if is_record(j + n, n) != r1:
                    ok = False
                if r1:
                    rec_res.add(j)
            r_x = len(rec_res)
            if r_x == 0:
                ok = False
            total = Fraction(0)
            for i in range(n):
                rot = x[i:] + x[:i]
                em = ends_at_min(rot)
                if em != ((i % n) in rec_res):
                    ok = False
                if em:
                    lad = desc_ladder(rot)
                    if len(lad) != r_x or set((i + l) % n for l in lad) != rec_res:
                        ok = False
                    total += Fraction(1, len(lad))
            if total != 1:
                ok = False
        # B_n^- as polynomial in x1 equals n sum_r E_(n,r)/r
        Bp = {}
        for (r, o), c in Erow.items():
            Bp[o] = Bp.get(o, Fraction(0)) + Fraction(n * c, r)
        Bdirect = {o: comb(n, o) for o in range(n + 1) if is_negative_state(o, n)}
        if {o: c for o, c in Bp.items() if c} != {o: Fraction(c) for o, c in Bdirect.items()}:
            ok = False
    check(okwin, "Step 4^-: record test over the window [j-n, j) agrees with [j-3n, j) (S_(l-n) = S_l - S_n > S_l when S_n < 0)")
    check(ok, "Step 4^-: for every x with S_n(x) < 0 (n <= 12): records exist, j record iff j+n record, rot_i(x) ends at its min iff i is a record residue, "
              "its descending ladder epochs number r(x), sum_i [..]/r(rot_i x) = 1; hence B_n^-(x_1) = n sum_r E^-_(n,r)(x_1)/r as polynomials")
    note("words with negative total examined: %d.  The per-word identity is unweighted; multiplying by x(x) and regrouping by y = rot_i(x)" % examined)
    note("uses only x(rot_i x) = x(x) (same letter multiset) and that x -> rot_i(x) is a bijection of the words with negative total: the")
    note("weighted Step 4 needs exactly rotation-invariance of the weight, nothing else.  Steps 1-3 are letter-multiset-preserving bijections.")
    # exact bivariate polynomial identity, n <= 48, both signs
    NP = 48
    for sign, name, rows in ((-1, "negative", negative_dp(NP)), (+1, "positive", positive_dp(NP))):
        ok = True
        for n in range(1, NP + 1):
            for o in range(n + 1):
                lhs = n * rows[n][o]
                rhs = 0
                for m in range(1, n + 1):
                    for o2 in range(0, min(o, m) + 1):
                        if (is_negative_state(o2, m) if sign < 0 else above(o2, m, 0, 1)):
                            c = rows[n - m][o - o2] if o - o2 < len(rows[n - m]) else 0
                            if c:
                                rhs += comb(m, o2) * c
                if lhs != rhs:
                    ok = False
        check(ok, "n P_n = sum_m B_m P_(n-m) holds as an EXACT identity of bivariate polynomials in (x_0, x_1) for %s words, n <= %d (hence for ALL weights)" % (name, NP))
    # Fractions for four weight pairs, n <= 30, both signs; and the audited P_1..P_7 lists
    NF = 30
    negrows = negative_dp(NF)
    posrows = positive_dp(NF)
    heads = {}
    ok = True
    for (x0, x1) in ((Fraction(1), Fraction(1)), (Fraction(1, 2), Fraction(2)), (Fraction(3), Fraction(1, 3)), (Fraction(5, 7), Fraction(11, 2))):
        for sign, rows in ((+1, posrows), (-1, negrows)):
            P = [sum(Fraction(c) * x0 ** (n - o) * x1 ** o for o, c in enumerate(rows[n])) for n in range(NF + 1)]
            B = [Fraction(0)] + [sum(comb(n, o) * x0 ** (n - o) * x1 ** o for o in range(n + 1)
                                     if (is_negative_state(o, n) if sign < 0 else above(o, n, 0, 1))) for n in range(1, NF + 1)]
            if not all(n * P[n] == sum(B[m] * P[n - m] for m in range(1, n + 1)) for n in range(1, NF + 1)):
                ok = False
            heads[(str(x0), str(x1), sign)] = [str(v) for v in P[1:8]]
    check(ok, "n P_n(x) = sum_m B_m(x) P_(n-m)(x) exactly (Fractions) for weights (1,1), (1/2,2), (3,1/3), (5/7,11/2), both signs, n <= 30")
    elapsed("section 2")
    return heads


# ====================================================================== section 3
def section3():
    print("== 3. The tilt identity Q(w) = 2^(-hm) 2^(lambda* S_m(w)) ==")
    note("h* = %s, lambda* = %s, lambda*/h* = %s, a* = lambda*/h* - 3/2 = %s" % (mpmath.nstr(H, 10), mpmath.nstr(LAM, 10), mpmath.nstr(LAM / H, 8), mpmath.nstr(ASTAR, 10)))
    check(abs(H - MP('0.9499555')) < 1e-7 and abs(LAM - MP('0.488077')) < 1e-6 and abs(ASTAR - MP('-0.986211')) < 1e-6 and abs(LAM / H - MP('0.5138')) < 1e-4,
          "the note's constants h* = 0.9499555, lambda* = 0.488077, lambda*/h* = 0.5138, a* = -0.986211 are correct to the digits printed")
    check(abs(RHO * ALPHA - 1) < MP(10) ** -55, "zero drift: rho* alpha = 1 (the letter 1 steps alpha - 1, the letter 0 steps -1; E_Q step = rho* alpha - 1 = 0)")
    check(abs(LAM * ALPHA - mpmath.log(RHO / (1 - RHO), 2)) < MP(10) ** -55, "coefficient of o: lambda* alpha = log_2(rho*/(1-rho*))")
    check(abs(-H - LAM - mpmath.log(1 - RHO, 2)) < MP(10) ** -55, "coefficient of m: -h - lambda* = log_2(1 - rho*)  [uses rho* = 1/alpha]")
    note("symbolic: log_2 Q(w) = o log_2 rho* + (m-o) log_2(1-rho*) = m log_2(1-rho*) + o log_2(rho*/(1-rho*)) = -hm - lambda* m + lambda* alpha o = -hm + lambda* S_m(w).")
    worst = MP(0)
    for m in range(1, 201, 7):
        for o in range(0, m + 1, 3):
            lhs = o * mpmath.log(RHO, 2) + (m - o) * mpmath.log(1 - RHO, 2)
            rhs = -H * m + LAM * (o * ALPHA - m)
            worst = max(worst, abs(lhs - rhs))
    check(worst < MP(10) ** -50, "numerically (60 digits): max |log_2 Q(w) - (-hm + lambda* S_m)| over a grid of (o, m), m <= 200, is %s" % mpmath.nstr(worst, 3))
    check(LAM > 0 and RHO > MP(1) / 2, "lambda* > 0 (rho* = 0.6309 > 1/2), so S_m > -y gives Q(w) >= 2^(-hm) 2^(-lambda* y)")
    elapsed("section 3")


# ====================================================================== section 4
def weighted_negative_dp_mp(n, x0, x1):
    """p_j = sum over negative words of length j of x0^(#0) x1^(#1), j <= n (mpmath), plus final row."""
    cur = [MP(1)]
    ps = [MP(1)]
    rows = [cur]
    for j in range(1, n + 1):
        nxt = [MP(0)] * (j + 1)
        for o, c in enumerate(cur):
            if c:
                nxt[o] += c * x0
                nxt[o + 1] += c * x1
        for o in range(j + 1):
            if not is_negative_state(o, j):
                nxt[o] = MP(0)
        cur = nxt
        rows.append(cur)
        ps.append(sum(cur))
    return ps, rows


def section4(negrows):
    print("== 4. Section 1.2: the tilt inequality and the exponential weight ==")
    MB = 60
    ys = [(1, 2), (1, 1), (2, 1), (4, 1), (8, 1)]
    ok1 = ok1b = ok2 = ok3 = True
    worst1 = MP(10)
    worst2 = MP(10)
    pQ, rowsQ = weighted_negative_dp_mp(MB, 1 - RHO, RHO)
    for (p, q) in ys:
        y = MP(p) / q
        for m in range(1, MB + 1):
            # exact N_m(y) and the Q-probability of {negative, S_m > -y}
            Nm = N_m_y(negrows, m, p, q)
            Qev = sum(rowsQ[m][o] for o in range(m + 1) if negrows[m][o] and above(o, m, p, q))
            # pointwise: every counted word has Q(w) >= 2^(-hm) 2^(-lambda* y)
            for o in range(m + 1):
                if negrows[m][o] and above(o, m, p, q):
                    Qw = RHO ** o * (1 - RHO) ** (m - o)
                    if not Qw >= 2 ** (-H * m) * 2 ** (-LAM * y) * (1 - MP(10) ** -40):
                        ok1b = False
            bound1 = 2 ** (H * m) * 2 ** (LAM * y) * Qev
            if not Nm <= bound1 * (1 + MP(10) ** -40):
                ok1 = False
            if Nm:
                worst1 = min(worst1, bound1 / Nm)
            for s in (MP('0.05'), MP('0.1'), MP('0.5')):
                # E_Q[e^(s S_m); negative] via the weights x1 = rho e^(s(alpha-1)), x0 = (1-rho) e^(-s)
                Es = sum(rowsQ[m][o] * mpmath.e ** (s * (o * ALPHA - m)) for o in range(m + 1) if negrows[m][o])
                Es_y = sum(rowsQ[m][o] * mpmath.e ** (s * (o * ALPHA - m)) for o in range(m + 1) if negrows[m][o] and above(o, m, p, q))
                bound2 = mpmath.e ** (s * y) * Es
                if not Qev <= mpmath.e ** (s * y) * Es_y * (1 + MP(10) ** -40):
                    ok2 = False          # on the event, e^(s(S_m + y)) >= 1
                if not Es_y <= Es * (1 + MP(10) ** -40):
                    ok3 = False          # dropping the constraint S_m > -y (nonnegative integrand)
                if Nm:
                    worst2 = min(worst2, 2 ** (H * m) * 2 ** (LAM * y) * bound2 / Nm)
    check(ok1b, "pointwise: every negative word with S_m > -y has Q(w) >= 2^(-hm) 2^(-lambda* y) (m <= 60, y in {1/2,1,2,4,8})")
    check(ok1, "first inequality N_m(y) <= 2^(hm) 2^(lambda* y) Q(negative, S_m > -y) (exact N_m, 60-digit Q); min ratio bound/N = %s" % mpmath.nstr(worst1, 6))
    check(ok2 and ok3, "second inequality Q(neg, S_m > -y) <= e^(sy) E_Q[e^(sS_m); neg, S_m > -y] <= e^(sy) E_Q[e^(sS_m); neg] for s in {0.05, 0.1, 0.5}; min ratio full bound/N = %s" % mpmath.nstr(worst2, 6))
    note("signs: lambda* > 0 and S_m > -y give the first; s > 0 and S_m + y > 0 on the event give e^(s(S_m+y)) >= 1 for the second; the integrand")
    note("e^(sS_m) 1_(negative) is nonnegative so the constraint S_m > -y can be dropped.  Both directions correct.")
    elapsed("section 4")


# ====================================================================== section 5
def section5():
    print("== 5. Sections 1.3-1.4: the weights, the local bound, the convolution ==")
    rho = float(RHO)
    alpha = float(ALPHA)
    hf = float(H)
    # weights reproduce Q(w) e^(s S_n)
    ok = True
    for s in (MP('0.05'), MP('0.1'), MP('0.5')):
        x1 = RHO * mpmath.e ** (s * (ALPHA - 1))
        x0 = (1 - RHO) * mpmath.e ** (-s)
        for n in range(1, 120, 11):
            for o in range(0, n + 1, 5):
                lhs = x0 ** (n - o) * x1 ** o
                rhs = RHO ** o * (1 - RHO) ** (n - o) * mpmath.e ** (s * (o * ALPHA - n))
                if abs(lhs - rhs) > abs(rhs) * MP(10) ** -50:
                    ok = False
    check(ok, "x_1 = rho* e^(s(alpha-1)), x_0 = (1-rho*) e^(-s) give x(w) = Q(w) e^(s S_n(w)) exactly (60 digits, grid)")
    # spacing of the admissible values o alpha - n, o alpha < n
    ok = True
    for n in range(1, 2001):
        vals = []
        o = 0
        while is_negative_state(o, n):
            vals.append(o)
            o += 1
        # the admissible o are exactly 0..o_top with o_top = max o with o alpha < n
        o_top = vals[-1]
        if any(is_negative_state(oo, n) for oo in range(o_top + 1, n + 1)):
            ok = False
        top = o_top * float(ALPHA) - n
        if not (-alpha < top < 0):
            ok = False
    check(ok, "admissible o (o alpha < n) form an initial segment 0..o_top, values o alpha - n negative, spaced by alpha, the largest in (-alpha, 0) (n <= 2000)")
    note("hence the exponential factors are at most 1, e^(-s alpha), e^(-2 s alpha), ... and g_n <= (max_o Bin) * 1/(1 - e^(-s alpha)).")
    # local binomial bound to n = 5000 (log-space, full row)
    import numpy as np
    lg = np.array([math.lgamma(i + 1) for i in range(5002)])
    lp, lq = math.log(rho), math.log(1 - rho)
    worst = 0.0
    worst_n = 0
    okmode = True
    for n in range(1, 5001):
        o = np.arange(n + 1)
        logpmf = lg[n] - lg[o] - lg[n - o] + o * lp + (n - o) * lq
        mx = float(np.max(logpmf))
        v = math.exp(mx) * math.sqrt(n)
        if v > worst:
            worst, worst_n = v, n
        mode = int(np.argmax(logpmf))
        if abs(mode - rho * n) > 1 or mode != math.floor((n + 1) * rho):
            okmode = False
    check(worst <= 3.2, "max_o Bin(n, rho*)(o) <= 3.2 n^(-1/2) for all 1 <= n <= 5000: max of sqrt(n) max_o Bin = %.4f at n = %d" % (worst, worst_n))
    check(okmode, "the mode floor((n+1) rho*) is within 1 of rho* n for all n <= 5000")
    # Stirling argument for n >= 11: Bin(mode) <= 1/sqrt(2 pi n q(1-q)), q = mode/n in [rho - 1/n, rho + 1/n]
    qlo, qhi = rho - 1 / 11, rho + 1 / 11
    qq = min(qlo * (1 - qlo), qhi * (1 - qhi))
    stir = 1 / math.sqrt(2 * math.pi * qq)
    check(stir <= 3.2 and 3.2 / math.sqrt(10) >= 1,
          "Stirling (C(n,k) <= 2^(n h(k/n))/sqrt(2 pi n q(1-q)) for 1 <= k <= n-1, Robbins) gives Bin(mode) <= %.3f n^(-1/2) for n >= 11 (q(1-q) >= %.4f); for n <= 10, 3.2/sqrt(n) >= 3.2/sqrt(10) = %.3f >= 1 >= Bin" % (stir, qq, 3.2 / math.sqrt(10)))
    # g_n and p_n for s in {0.05, 0.1, 0.5}, n <= 2000, and the recurrence
    results = {}
    for s in (0.05, 0.1, 0.5):
        Cs = 3.2 / (1 - math.exp(-s * alpha))
        NMAX = 2000
        g = [0.0] * (NMAX + 1)
        for n in range(1, NMAX + 1):
            o = np.arange(n + 1)
            logpmf = lg[n] - lg[o] - lg[n - o] + o * lp + (n - o) * lq
            tops = [oo for oo in range(n + 1) if is_negative_state(oo, n)]
            o_top = tops[-1]
            sel = o[: o_top + 1]
            g[n] = float(np.sum(np.exp(logpmf[: o_top + 1] + s * (sel * alpha - n))))
        gmax = max(g[n] * math.sqrt(n) for n in range(1, NMAX + 1))
        # weighted DP for p_n (floats), states o <= hi_j
        x1 = rho * math.exp(s * (alpha - 1))
        x0 = (1 - rho) * math.exp(-s)
        cur = np.array([1.0])
        p = [1.0]
        for j in range(1, NMAX + 1):
            nxt = np.concatenate([cur, [0.0]]) * x0 + np.concatenate([[0.0], cur]) * x1
            hi = 0
            while hi + 1 <= j and is_negative_state(hi + 1, j):
                hi += 1
            nxt[hi + 1:] = 0.0
            cur = nxt
            p.append(float(np.sum(cur)))
        # recurrence n p_n = sum_m g_m p_(n-m) in floats
        maxrel = 0.0
        for n in range(1, NMAX + 1):
            rhs = sum(g[m] * p[n - m] for m in range(1, n + 1))
            maxrel = max(maxrel, abs(n * p[n] - rhs) / (n * p[n]))
        pn32 = [p[n] * n ** 1.5 for n in range(1, NMAX + 1)]
        sigma_num = sum(g[n] / n for n in range(1, NMAX + 1))
        tail = 2 * Cs / math.sqrt(NMAX)
        Dp_num = Cs * math.exp(2 * math.sqrt(2) * (sigma_num + tail))
        results[s] = dict(Cs=Cs, gmax=gmax, pmax=max(pn32), pmax_late=max(pn32[499:]), pmin_late=min(pn32[499:]),
                          p_last=pn32[-1], sigma_num=sigma_num, tail=tail, Dp_num=Dp_num, maxrel=maxrel, p=p, g=g)
        check(gmax <= Cs, "s = %.2f: g_n sqrt(n) <= C_s = %.3f for n <= 2000 (max %.4f); g_n = E_Q[e^(sS_n); S_n < 0]" % (s, Cs, gmax))
        check(maxrel < 1e-9, "s = %.2f: the weighted identity n p_n = sum_m g_m p_(n-m) holds numerically for the real weights, n <= 2000 (max rel err %.1e)" % (s, maxrel))
        check(max(pn32) <= Dp_num, "s = %.2f: p_n n^(3/2) in [%.4f, %.4f] for 500 <= n <= 2000 (max over all n <= 2000: %.4f), bounded; D'_s(numerical sigma) = %.3e >= it" % (s, min(pn32[499:]), max(pn32[499:]), max(pn32), Dp_num))
        note("s = %.2f: sigma_s numerical = %.4f (+ tail <= %.4f), proved bound C_s zeta(3/2) = %.3f; p_n n^(3/2) at n = 2000: %.4f" % (s, sigma_num, tail, Cs * float(mpmath.zeta(1.5)), pn32[-1]))
    # convolution lemma: statement and induction step re-derived numerically on a_n = n^(-3/2)
    K = 5000
    a = [0.0] + [n ** -1.5 for n in range(1, K + 1)]
    z32 = float(mpmath.zeta(1.5))
    worstc = 0.0
    for k in (2, 3, 5, 10, 50, 100, 500, 1000, 2000, 5000):
        conv = sum(a[i] * a[k - i] for i in range(1, k))
        worstc = max(worstc, conv * k ** 1.5)
    check(worstc <= 2 ** 1.5 * 2 * z32, "convolution lemma on a = b = n^(-3/2): k^(3/2) (a*a)_k <= 2^(3/2) * 2 zeta(3/2) = %.3f (max observed %.3f)" % (2 ** 1.5 * 2 * z32, worstc))
    note("lemma proof: split at i = k/2; i >= k/2: a_i <= A (k/2)^(-3/2), sum b <= beta; i < k/2: b_(k-i) <= B (k/2)^(-3/2), sum a <= alpha; total 2^(3/2)(A beta + B alpha) k^(-3/2).")
    note("induction (g^(*m))_k <= C m (2^(3/2) sigma)^(m-1) k^(-3/2): step needs 2^(3/2) sigma^m <= (2^(3/2) sigma)^m, true for m >= 1;")
    note("exp(G) = sum_m G^(*m)/m! with G_n = g_n/n <= C_s n^(-3/2), G_0 = 0, G >= 0, sum G = sigma_s <= C_s zeta(3/2): p_n <= C_s e^(2 sqrt2 sigma_s) n^(-3/2).")
    # numerical check of the induction bound on the actual G for s = 0.1
    s = 0.1
    Cs = results[s]['Cs']
    g = results[s]['g']
    NM = 400
    G = np.array([0.0] + [g[n] / n for n in range(1, NM + 1)])
    sig = float(np.sum(G))
    conv = G.copy()
    okind = True
    for m in range(1, 8):
        bound = Cs * m * (2 ** 1.5 * sig) ** (m - 1)
        for k in range(1, NM + 1):
            if conv[k] > bound * k ** -1.5 * (1 + 1e-12):
                okind = False
        conv = np.convolve(conv, G)[: NM + 1]
    check(okind, "s = 0.1: the induction bound (G^(*m))_k <= C_s m (2^(3/2) sigma)^(m-1) k^(-3/2) holds numerically for m <= 7, k <= 400 (sigma = numerical sum to 400)")
    elapsed("section 5")
    return results


# ====================================================================== section 6
def section6(results, Mtable, Wall):
    print("== 6. Section 1.5: the assembly ==")
    z32 = mpmath.zeta(MP(3) / 2)
    c2 = 545
    const = 2 ** MP(1.5) * 2 * z32
    check(const <= MP('14.8'), "2^(3/2) * 2 * zeta(3/2) = %s <= 14.8" % mpmath.nstr(const, 6))
    worst = 0.0
    for k in list(range(2, 200)) + [500, 1000, 1600, 3000]:
        v = sum(m ** -1.5 * (k - m) ** -1.5 for m in range(1, k)) * k ** 1.5
        worst = max(worst, v)
    check(worst <= 14.8, "k^(3/2) sum_(m=1)^(k-1) m^(-3/2) (k-m)^(-3/2) <= 14.8 for k <= 200 and k in {500, 1000, 1600, 3000} (max %.3f)" % worst)
    note("m = k term: N_k(y) W_0 with W_0 = 1 gives D'_s 2^(lambda* y) e^(sy) 2^(hk) k^(-3/2), the '1' inside the bracket [1 + 14.8 c_2];")
    note("the W_k term c_2 2^(hk) k^(-3/2) is absorbed using 2^(lambda* y) e^(sy) >= 1 (y >= 0, lambda*, s > 0): D_s = c_2 + D'_s (1 + 14.8 c_2).")
    # W_k <= 545 2^(hk) k^(-3/2) for the W_k computed here
    ok = all(Wall[k] <= 545 * 2 ** (H * k) * MP(k) ** -1.5 for k in range(1, len(Wall)))
    check(ok, "THM-4495's W_k <= 545 2^(hk) k^(-3/2) holds for the independently computed W_k, k <= %d" % (len(Wall) - 1))
    for s in (0.05, 0.1, 0.5):
        S = MP(s)
        Cs = MP('3.2') / (1 - mpmath.e ** (-S * ALPHA))
        sig = Cs * z32
        Dp = Cs * mpmath.e ** (2 * mpmath.sqrt(2) * sig)
        Ds = c2 + Dp * (1 + MP('14.8') * c2)
        okb = True
        worstlog = -MP(10) ** 9
        for (k, y), M in Mtable.items():
            bound = Ds * 2 ** (H * k) * MP(k) ** -1.5 * 2 ** (LAM * y) * mpmath.e ** (S * y)
            if not M <= bound:
                okb = False
            worstlog = max(worstlog, mpmath.log10(MP(M) / bound))
        check(okb, "s = %.2f: C_s = %s, sigma bound = %s, log10 D'_s = %s, log10 D_s = %s; M_k(y) <= D_s 2^(hk) k^(-3/2) 2^(lambda* y) e^(sy) on the (B) grid (max log10 M/bound = %s)"
              % (s, mpmath.nstr(Cs, 6), mpmath.nstr(sig, 6), mpmath.nstr(mpmath.log10(Dp), 5), mpmath.nstr(mpmath.log10(Ds), 5), mpmath.nstr(worstlog, 4)))
    # the same chain with the NUMERICAL p-constant (sanity, not proof): P_s = max_n p_n n^(3/2)
    for s in (0.05, 0.1, 0.5):
        Ps = results[s]['pmax']
        Dreal = c2 + Ps * (1 + 14.8 * c2)
        okb = True
        worst = -float('inf')
        for (k, y), M in Mtable.items():
            lnbound = mpmath.log(MP(Dreal)) + H * k * LN2 - MP(1.5) * mpmath.log(k) + LAM * y * LN2 + MP(s) * y
            r = float(mpmath.log(MP(M)) - lnbound)
            worst = max(worst, r)
            if r > 0:
                okb = False
        check(okb, "s = %.2f (sanity): with the numerical P_s = max_(n<=2000) p_n n^(3/2) = %.3f in place of D'_s the chain still bounds the (B) grid (max ln M/bound = %.2f)" % (s, Ps, worst))
    elapsed("section 6")


# ====================================================================== section 7
def T_step(x, b):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


def section7(Mfun):
    print("== 7. Lemma 1.4c: the counting lemma with the ballot factor ==")
    # carry bound, exact
    ok = True
    for b in (1, -1, 3, -5, 7):
        for k in range(0, 13):
            for w in product((0, 1), repeat=k):
                beta = Fraction(0)
                for i, c in enumerate(w, 1):
                    beta = beta / 2 if c == 0 else (3 * beta + b) / 2
                    if abs(beta) > abs(b) * (Fraction(3, 2) ** i - 1):
                        ok = False
                    if abs(beta) >= abs(b) * Fraction(3, 2) ** i:
                        ok = False
    check(ok, "carry bound |beta_i(w)| <= |b| ((3/2)^i - 1) < |b| (3/2)^i for all words of length <= 12, b in {1,-1,3,-5,7} (exact)")
    note("for every i <= k = floor(log_2 X): |beta_i| < |b| (3/2)^i <= |b| (3/2)^k <= |b| X^(log_2(3/2)) = Y_0 X^(-theta)/2 <= n X^(-theta)/2 for n >= Y_0;")
    note("no-dip T_b^i(n) >= n X^(-theta) then gives 3^(o_i) n / 2^i >= n X^(-theta) - |beta_i| > n X^(-theta)/2, i.e. S_i > -theta log_2 X - 1 STRICTLY")
    note("(the strictness, hidden in the note's 'i.e.', comes from the '-1' of the carry bound; it matters only when theta log_2 X is an integer and o_i = 0).")
    # direct enumeration at X = 2^18 on both sheets
    X = 2 ** 18
    Lx = 18
    k = 18
    for b in (1, -1):
        for theta in (Fraction(3, 100), Fraction(1, 10)):
            thf = float(theta)
            Y0 = 2 * abs(b) * X ** (math.log2(1.5) + thf)
            y = theta * Lx + 1                       # exact rational barrier
            p, q = y.numerator, y.denominator
            F = []
            for n in range(1, X + 1):
                x = n
                good = True
                thr = n * X ** (-thf)
                for i in range(1, k + 1):
                    x = T_step(x, b)
                    if x < thr:
                        good = False
                        break
                if good:
                    F.append(n)
            big = [n for n in F if n >= Y0]
            okbar = True
            for n in big:
                x = n
                o = 0
                for i in range(1, k + 1):
                    if x % 2:
                        o += 1
                    x = T_step(x, b)
                    if not above(o, i, p, q):
                        okbar = False
            M = Mfun(k, p, q)
            check(okbar, "b = %+d, theta = %s, X = 2^18: every n >= Y_0 = %.0f in F_b(X, theta) (%d of #F = %d) has S_i > -(theta log_2 X + 1) for ALL 1 <= i <= k (exact barrier y = %s)"
                  % (b, theta, Y0, len(big), len(F), y))
            check(len(F) <= Y0 + 2 * M, "b = %+d, theta = %s: #F_b = %d <= Y_0 + 2 M_k(y) = %.0f + 2 * %d" % (b, theta, len(F), Y0, M))
    # two representatives per class, exponents, constants
    ok = True
    for L in (MP(200), MP(201.5), MP(1000), MP(12345.678)):
        k = int(mpmath.floor(L))
        if not (k >= L - 1 and L - 1 >= MP('0.995') * L and 2 ** (H * k) <= 2 ** (H * L)):
            ok = False
    check(ok, "k = floor(log_2 X) >= log_2 X - 1 >= 0.995 log_2 X for X >= 2^200; 2^(hk) <= X^(h*)")
    note("a residue class mod 2^k meets [1, X] in at most floor((X-1)/2^k) + 1 <= 2 integers since 2^k > X/2.")
    c = 2 * 2 ** LAM * mpmath.e ** MP('0.1') * MP('0.995') ** MP(-1.5)
    check(c <= 4, "2^(lambda*) e^(s) 2 (0.995)^(-3/2) = %s <= 4 for s <= 0.1 (monotone in s)" % mpmath.nstr(c, 6))
    note("Lemma M needs no restriction on theta; the hypothesis theta <= theta_1 of Lemma 1.4c is inherited from 1.4b and unused here (harmless).")
    elapsed("section 7")


# ====================================================================== section 8
def section8():
    print("== 8. Section 3: the bootstrap ==")
    # (i) recursion (R'') = (R) of THM-4476 with Lemma 1.4c: quote the source lines
    thin = open(os.path.join(ROOT, '05-knowledge', 'results', 'collatz_thin_20260925_thin_divergent_orbits.md'), encoding='utf-8').read()
    r_line = [ln for ln in thin.splitlines() if '(R)' in ln and 'N(X) <= k + k N(X^(1-theta))' in ln]
    check(len(r_line) == 1, "THM-4476's recursion (R) found in its note: " + (r_line[0].strip() if r_line else "MISSING"))
    note("(R'') replaces the counting-lemma term by Lemma 1.4c's bound; log_2(3/2) = %s <= 0.585 so X^(log_2(3/2)+theta) <= X^(0.585+theta)." % mpmath.nstr(mpmath.log(MP(3) / 2, 2), 6))
    # (ii)-(v) for several a
    okall = True
    for a in (MP('-0.98'), MP('-0.9'), MP('-0.5'), MP(0), MP('0.5'), MP(2)):
        eta = (a - ASTAR) / 8
        c1 = (1 + eta) / H
        s = min(MP('0.1'), eta * LN2 / c1)
        lhs = LAM * c1 - MP(3) / 2 + eta
        ok_ii = (abs(lhs - (ASTAR + eta * LAM / H + eta)) < MP(10) ** -50) and (lhs <= ASTAR + MP('1.52') * eta) and (abs((ASTAR + MP('1.52') * eta) - (a - MP('6.48') * eta)) < MP(10) ** -50) and (a - MP('6.48') * eta <= a - 2 * eta)
        ok_iii = (s * c1 / LN2 <= eta) and (0 < s <= MP('0.1'))
        # (iv): X^(0.585+theta_X) L^(-(a-2eta)) = o(X^(h*)): check at L = 2^11 .. and the sign of a - 2 eta > a*
        okiv = (a - 2 * eta > ASTAR) and (a - 2 * eta > -1)
        # (v) the bracket: for L with 2^|a| L^(-eta) <= 1/4 and L^(-2eta) <= 1/4, and K >= 8 D_s + 1 (D_s := 1e75 here)
        Ds = MP(10) ** 75
        K = 8 * Ds + 1
        L3 = max((4 * 2 ** abs(a)) ** (1 / eta), 4 ** (1 / (2 * eta)), MP(200))
        okv = True
        for L in (L3 * (1 + MP(10) ** -30), 2 * L3, 10 * L3, 1000 * L3):   # L_3 itself is the equality case up to rounding
            if not (2 ** abs(a) * L ** (-eta) <= MP(1) / 4 and L ** (-2 * eta) <= MP(1) / 4):
                okv = False
            bracket = L ** (-2 * eta) + 2 ** abs(a) * K * L ** (-eta) + 4 * Ds * L ** (-2 * eta)
            if not bracket <= (1 + 4 * Ds) / 4 + K / 4 <= K:
                okv = False
            # Y = X L^(-c1) in [X^(1/2), X): needs c1 log_2 L <= L/2
            if not c1 * mpmath.log(L, 2) <= L / 2:
                okv = False
            logY = L - c1 * mpmath.log(L, 2)
            if not (L / 2 <= logY <= L and c1 * mpmath.log(L, 2) > 0):    # strictness log_2 Y < L from c_1 log_2 L > 0 (L ~ 10^1155 rounds)
                okv = False
            # (log_2 Y)^a <= 2^|a| L^a, either sign
            if not logY ** a <= 2 ** abs(a) * L ** a:
                okv = False
            # Y^(h*) = X^(h*) L^(-1-eta)
            if abs((H * logY) - (H * L - (1 + eta) * mpmath.log(L, 2))) > MP(10) ** -40 * H * L:   # relative: L can be 10^1155
                okv = False
            # the uniform gap for well-foundedness: Y <= X/200
            if not L ** c1 >= 200:
                okv = False
        okall = okall and ok_ii and ok_iii and okiv and okv
        note("a = %s: eta = %s, c_1 = %s, s = %s, lambda* c_1 - 3/2 + eta = %s <= a - 2 eta = %s; L_3 ~ %s; (ii) %s (iii) %s (iv) %s (v) %s"
             % (mpmath.nstr(a, 4), mpmath.nstr(eta, 5), mpmath.nstr(c1, 6), mpmath.nstr(s, 5), mpmath.nstr(lhs, 6), mpmath.nstr(a - 2 * eta, 6), mpmath.nstr(L3, 4), ok_ii, ok_iii, okiv, okv))
    check(okall, "(ii) lambda* c_1 - 3/2 + eta = a* + eta(1 + lambda*/h*) <= a* + 1.52 eta = a - 6.48 eta <= a - 2 eta; (iii) e^(s theta_X L) = L^(s c_1/ln 2) <= L^eta; "
                 "(iv) a - 2 eta > a* > -1; (v) Y in [X^(1/2), X), log_2 Y in [L/2, L), (log_2 Y)^a <= 2^|a| L^a, Y^(h*) = X^(h*) L^(-1-eta), bracket <= K, and L^(c_1) >= 200")
    note("(iv) X^(0.585+theta_X) L^(-(a-2eta)) = X^(0.585) L^(c_1 - a + 2 eta) = o(X^(h*)) since 0.585 < h* = 0.95: the polylog is harmless for either sign of a;")
    note("     the note's 'O(X^(0.7))' is correct once theta_X <= 0.115, and X^(0.7) <= X^(h*) L^(a-2eta) for large X because h* - 0.7 = 0.25 > 0.")
    note("(v)  base: N(X) <= X (distinct positive integers <= X) and K >= sup_(2<=X<X_3) X^(1-h*) (log_2 X)^(-a), finite on a bounded range; K >= 8 D_s + 1 >= (1+4D_s)/3 suffices.")
    note("(vi) well-founded: Y = X L^(-c_1) <= X/200 for L >= 200 (c_1 > 1), a uniform gap; if the set of counterexamples were nonempty with infimum X_0 >= X_3,")
    note("     a counterexample X < 2 X_0 has Y <= X/200 < X_0, so the claim holds at Y and hence at X: contradiction.  'Y < X' alone would not suffice for a real variable.")
    # (vii) Corollary 4's recursion in 1.7
    c4 = [ln for ln in thin.splitlines() if ln.startswith('So `N_A(X) <= k N_A(X^(1-theta))')]
    c4b = [ln for ln in thin.splitlines() if 'N_A(X) <= k N_A(X^(1-theta))' in ln and not ln.startswith('So `')]
    check(len(c4) == 1 and len(c4b) == 1, "Corollary 4's recursion found in 1.7 of THM-4476's note (quoted once more in its addendum 1.6b): " + (c4[0].strip()[:150] if c4 else "MISSING"))
    note("(vii) it is (R) with #F_b replaced by #F_b(X,theta) + #F_(-b)(X,theta) and the (E)-term k by (k+1)(|b|/3 + 1), base N_A(X) <= 2X: NOT verbatim (R),")
    note("      but Lemma 1.4c applies to both F_b and F_(-b) with the same |b|, so the bootstrap goes through with 8 D_s in place of 4 D_s and the O(L) term absorbed in X_2.")
    elapsed("section 8")


# ====================================================================== section 9
def section9():
    print("== 9. Consistency with THM-4487(2) ==")
    hp = mpmath.log((1 - RHO) / RHO, 2)
    check(abs(hp + ALPHA * LAM) < MP(10) ** -55, "h'(rho*) = log_2((1-rho*)/rho*) = -alpha lambda*")
    ok = True
    okgt = True
    worst = MP(-1)
    for i in range(1, 1001):
        th = MP(i) / 10000
        rho = (1 - th) / ALPHA
        hr = -(rho * mpmath.log(rho, 2) + (1 - rho) * mpmath.log(1 - rho, 2))
        if not hr <= H + LAM * th:
            ok = False
        if not hr > H:
            okgt = False
        worst = max(worst, hr - (H + LAM * th))
    check(ok, "concavity: h(rho(theta)) <= h* + lambda* theta for theta in (0, 0.1] (max h(rho) - h* - lambda* theta = %s < 0)" % mpmath.nstr(worst, 4))
    check(okgt, "h(rho(theta)) > h* for theta in (0, 0.1] (rho < rho* < ... h decreasing on [1/2, 1]): THM-4487's fixed-theta exponent exceeds h*")
    note("THM-4487(2): c X^(h(rho)) L^(-3/2) <= #F_b <= C X^(h(rho)) L^2.  Lemma 1.4c: #F_b <= 2|b| X^(0.585+theta) + 4 D_s X^(h* + lambda* theta + s theta/ln 2) L^(-3/2).")
    note("upper >= lower needs h(rho) <= h* + lambda* theta + s theta/ln 2: true by concavity (strictly, for theta > 0). No contradiction; for fixed theta the")
    note("lemma is WEAKER than THM-4487's upper bound in the exponent and only wins in the polynomial factor as theta -> 0: 'the ballot factor is a polynomial, not an exponent' is accurate.")
    elapsed("section 9")


# ====================================================================== section 10
def sha256_file(path):
    data = open(path, 'rb').read()
    return hashlib.sha256(data).hexdigest(), hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest(), data


def section10(Mtable, heads, Wall):
    print("== 10. Provenance, the float-margin barrier test, the (B) table and (A) lists reproduced ==")
    thm = open(os.path.join(ROOT, '01-canon', 'theorems', 'THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md'), encoding='utf-8').read()
    hs = re.search(r'script_sha256:\s*([0-9a-f]{64})', thm).group(1)
    ho = re.search(r'output_sha256:\s*([0-9a-f]{64})', thm).group(1)
    sp = os.path.join(ROOT, '04-computation', 'experiments', 'collatz_thin_20260926_movingbarrier.py')
    op = os.path.join(ROOT, '04-computation', 'experiments', 'collatz_thin_20260926_movingbarrier.out')
    raw_s, lf_s, _ = sha256_file(sp)
    raw_o, lf_o, out_bytes = sha256_file(op)
    check(raw_s == hs and lf_s == hs, "script sha256 (raw = LF-normalised, no CR bytes) matches the header: " + hs[:16] + "...")
    check(raw_o == ho and lf_o == ho, "output sha256 (raw = LF-normalised) matches the header: " + ho[:16] + "...")
    # the float-margin test: min |o alpha - j + y| for o >= 1
    cf = []
    x = ALPHA
    for _ in range(14):
        qd = int(mpmath.floor(x))
        cf.append(qd)
        x = 1 / (x - qd)
    p0, q0, p1, q1 = 1, 0, cf[0], 1
    conv = [(p1, q1)]
    for c in cf[1:]:
        p0, q0, p1, q1 = p1, q1, c * p1 + p0, c * q1 + q0
        conv.append((p1, q1))
    note("continued fraction of log_2 3: " + str(cf) + "; convergents " + str(conv[:11]))
    best = (MP(1), 0)
    for o in range(1, 3201):
        v = o * ALPHA
        d = abs(v - mpmath.nint(v))
        if d < best[0]:
            best = (d, o)
    best_even = (MP(1), 0)
    for o in range(2, 3201, 2):
        v = o * ALPHA
        d = abs(v - mpmath.nint(v))
        if d < best_even[0]:
            best_even = (d, o)
    # float error of o*ALPHA_float for o <= 1600
    af = float(ALPHA)
    ferr = max(abs(MP(o * af) - o * ALPHA) for o in range(1, 1601))
    check(best[1] == 665 and best[0] > MP('6e-5') and ferr < MP('1e-11'),
          "min_(1<=o<=3200) ||o log_2 3|| = %s at o = %d (denominator of the convergent 1054/665; next 24727/15601); max float error of o*ALPHA, o <= 1600: %s"
          % (mpmath.nstr(best[0], 4), best[1], mpmath.nstr(ferr, 3)))
    note("for integer y the script's v = o alpha - j + y has |v| >= 6.3e-5 for o in [1, 1600]; for y = 1/2, |v| = |2o alpha - (2j-1)|/2 >= %s/2 (o' = 2o even, min at o' = 1330);"
         % mpmath.nstr(best_even[0], 4))
    note("the 1e-9 margin is therefore never approached (gap/margin > 6e4) and the float sign is exact: the barrier test is justified, o = 0 being compared exactly.")
    # reproduce the (B) table
    Hf, Lf = float(H), float(LAM)
    rows = {}
    for ln in out_bytes.decode().splitlines():
        m = re.match(r'\s*k=\s*(\d+):\s*(.*?)\s*\|', ln)
        if m:
            rows[int(m.group(1))] = [float(v) for v in m.group(2).split()]
    ys = [0, 1, 2, 4, 8, 16, 32]
    ok = True
    worst = 0.0
    for k, vals in rows.items():
        for y, v in zip(ys, vals):
            M = Mtable[(k, y)]
            mine = float(2 ** (mpmath.log(M, 2) - H * k + MP(1.5) * mpmath.log(k, 2) - LAM * y))
            worst = max(worst, abs(mine - v))
            if abs(mine - v) > 0.0015:
                ok = False
    check(ok and len(rows) == 6, "the (B) table of the .out (6 rows x 7 barriers) is reproduced by the exact DP with 60-digit normalisation (max |diff| = %.5f, 3 printed decimals)" % worst)
    for k in sorted(rows):
        note("k = %4d: exact M_k(y) bit lengths %s; normalised %s" % (k, [Mtable[(k, y)].bit_length() for y in ys],
             ["%.3f" % float(2 ** (mpmath.log(Mtable[(k, y)], 2) - H * k + MP(1.5) * mpmath.log(k, 2) - LAM * y)) for y in ys]))
    # (A) P_1..P_7 lists
    okA = True
    for ln in out_bytes.decode().splitlines():
        m = re.match(r"\s*weights \(x0, x1\) = \(([^,]+), ([^)]+)\), (positive|negative) words: (\w+); P_1..P_7 = \[(.*)\]", ln)
        if m:
            x0, x1, kind, flag, lst = m.groups()
            sign = +1 if kind == 'positive' else -1
            lst = [v.strip().strip("'") for v in lst.split(',')]
            key = (str(Fraction(x0)), str(Fraction(x1)), sign)
            if flag != 'True' or heads.get(key) != lst:
                okA = False
    check(okA, "the (A) lists P_1..P_7 of the .out (three weight pairs, both signs) agree with the independent weighted DP; all six flags True")
    check([Wall[k] for k in range(1, 21)] == [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328],
          "W_1..W_20 from the exact barrier DP at y = 0 equal THM-4495's list")
    note("status-block claims: 'PROVED (elementary; rests on THM-4476's (R) and THM-4495's apparatus)', 'every a > a*', 'N(X) = o(X^(h*))', 'same for injective")
    note("invariant sets' -- all supported by sections 1-9 (the last with the non-verbatim recursion of 1.7, see section 8).  Overclaim in wording only: the note's")
    note("'the truth is M_k(y) ~ c (y+1) ...' is an observation (the /(y+1) column moves from 10.8 to 15.0 at k = 1600), not a proved statement; not used in the proof.")
    elapsed("section 10")


# ====================================================================== main
def main():
    try:
        sys.stdout.reconfigure(newline='\n')
    except Exception:
        pass
    print("collatz_thin_20260926_movingbarrier_audit.py: independent audit of THM-4499 (moving-barrier ballot bound, weighted Spitzer identity, bootstrap)")
    print("alpha = log_2 3 = %s, rho* = %s, h* = %s, lambda* = %s, a* = %s (mpmath, 60 digits)"
          % (mpmath.nstr(ALPHA, 15), mpmath.nstr(RHO, 15), mpmath.nstr(H, 15), mpmath.nstr(LAM, 15), mpmath.nstr(ASTAR, 15)))
    section1()
    heads = section2()
    section3()
    negrows = negative_dp(60)
    section4(negrows)
    results = section5()
    # exact M_k(y) on the (B) grid, and W_k for all k <= 1600
    Mtable = {}
    Wall = None
    Mcache = {}

    def Mfun(k, p, q):
        key = (k, p, q)
        if key not in Mcache:
            Mcache[key] = barrier_dp(k, p, q)[0][k]
        return Mcache[key]
    for y in (0, 1, 2, 4, 8, 16, 32):
        Ms, _ = barrier_dp(1600, y, 1)
        if y == 0:
            Wall = Ms
        for k in (100, 200, 400, 800, 1200, 1600):
            Mtable[(k, y)] = Ms[k]
    elapsed("(B) grid")
    section6(results, Mtable, Wall)
    section7(Mfun)
    section8()
    section9()
    section10(Mtable, heads, Wall)
    print("== summary ==")
    print("  checks: %d, failures: %d" % (NCHECK[0], len(FAILS)))
    for f in FAILS:
        print("  FAILED: " + f)
    print("  verdict of the numerical/exact part: %s" % ("all checks passed" if not FAILS else "FAILURES PRESENT"))


if __name__ == '__main__':
    main()
