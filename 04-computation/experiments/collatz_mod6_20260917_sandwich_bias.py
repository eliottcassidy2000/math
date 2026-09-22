#!/usr/bin/env python3
"""
collatz_mod6_20260917_sandwich_bias.py -- signed sandwich discrepancy census.

Lane: sandwich_bias (session collatz-mod6-20260917, machine mac-mini).
Finalized 2026-09-21 from the recovered 2026-09-17 draft: the draft lacked
section 4(iv) (pooled-gap control), the three-way split of the predicted
drift, and used an incomplete seven-stratum mod-35 partition (the audits
found the missing strata 35|left and 35|right, off by the SP pair (35,37));
its peak RSS was 1.7 GB.  All of that is repaired here; every check is an
explicit `raise` so it survives python -O.

Object.  Centers W = 6k, endpoints 6k-1 (class 5 mod 6) and 6k+1 (class 1
mod 6).  Omega = number of prime factors with multiplicity, omega = number
of distinct prime factors, b = number of prime factors congruent to 5 mod 6
with multiplicity.  chi(n) = +1 if n = 1 mod 6, -1 if n = 5 mod 6
(the nontrivial character mod 3 restricted to (n,6)=1).

Inheritance (read, not re-derived):
  05-knowledge/results/arithmetic_braids_20260917_divisors.md
    SW1 ordered matrix N_ij(K) and the exact CRT sandwich identity,
    SW2 k->-k side symmetry, SW3 K=10^6 4x4 Omega matrix,
    SW4 class-parity law b(6k-1) odd / b(6k+1) even.
This script re-derives the K=10^6 matrix only as a cross-check and raises
if it disagrees with the inherited table.

Universe.  Omega, omega, b sieved exactly on [1, 6*10^7 + 301] (uint8);
the +301 covers the pooled-gap partner n+g with g <= 300.

Sections printed:
  0  sieve + independent trial-division audit (positive control)
  1  4x4 Omega matrices at K = 10^5, 10^6, 10^7; antisymmetric part;
     running sign of N_PS - N_SP, N_PC - N_CP, N_SC - N_CS for K <= 10^7
  2  class decomposition: primes / semiprimes / 3-almost-primes by class
     pattern, chi-sums by Omega, PROVED identity for S_1 - S_5 checked
     exactly, class-labelled exponent-shape determination for Omega<=3 and
     its Omega=4 failures, independence (product-of-marginals) prediction,
     its two-mechanism split and the three-way split in both conventions
  3  3-almost-prime class patterns (omega, b) on each side
  4  hostile controls: (i) shifted pairs; (ii) the nine mod-35 strata;
     (iii) prime race; (iv) pooled-gap control with block-level residual test
  5  summary verdicts

Run:  python3 04-computation/experiments/collatz_mod6_20260917_sandwich_bias.py
"""
import math
import resource
import sys
import time
from fractions import Fraction

import numpy as np

T0 = time.time()
GMAX = 300                       # largest pooled gap
NMAX = 6 * 10**7 + 1 + GMAX      # sieve bound: largest partner endpoint used
KMAX = 10**7                     # largest number of centers
K_SCALES = [10**5, 10**6, 10**7]
K_TRACK = [10**4, 2 * 10**4, 5 * 10**4, 10**5, 2 * 10**5, 5 * 10**5,
           10**6, 2 * 10**6, 5 * 10**6, 10**7]
NB_FINE = 1000                   # blocks of 10^4 centers over K=10^7 (first 100 = K=10^6 blocks)

INHERITED_K1E6 = [
    [37915, 78689, 59706, 30192],
    [78277, 157420, 112416, 52182],
    [59992, 112305, 72363, 28755],
    [30161, 52125, 28736, 8766],
]


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def elapsed():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_gb = rss / 2**30 if sys.platform == "darwin" else rss / 2**20   # bytes on macOS, KiB on Linux
    return "[t=%.1fs maxrss=%.2fGB]" % (time.time() - T0, rss_gb)


def hr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# --------------------------------------------------------------------------
# 0. Sieve
# --------------------------------------------------------------------------
hr("0. SIEVE Omega / omega / b on [1, %d]  (uint8, exact)" % NMAX)
is_p = np.ones(NMAX + 1, dtype=np.bool_)
is_p[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if is_p[i]:
        is_p[i * i::i] = False
PRIMES = np.nonzero(is_p)[0].astype(np.int64)
del is_p
OM = np.zeros(NMAX + 1, np.uint8)     # Omega
OMD = np.zeros(NMAX + 1, np.uint8)    # omega (distinct)
BB = np.zeros(NMAX + 1, np.uint8)     # b = # prime factors 5 mod 6, with mult.
for p in map(int, PRIMES):
    OMD[p::p] += 1
    pj = p
    if p % 6 == 5:
        while pj <= NMAX:
            OM[pj::pj] += 1
            BB[pj::pj] += 1
            pj *= p
    else:
        while pj <= NMAX:
            OM[pj::pj] += 1
            pj *= p
print("primes <= %d : %d   %s" % (NMAX, len(PRIMES), elapsed()))
check(int(OM[1]) == 0 and int(OM[2]) == 1 and int(OM[12]) == 3 and int(OM[25]) == 2, "Omega spot")
check(int(BB[25]) == 2 and int(BB[35]) == 1 and int(BB[49]) == 0 and int(BB[125]) == 3, "b spot")


def trial(n):
    """Independent trial division: (Omega, omega, b)."""
    O = w = b = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            w += 1
            while n % d == 0:
                n //= d
                O += 1
                if d % 6 == 5:
                    b += 1
        d += 1
    if n > 1:
        O += 1
        w += 1
        if n % 6 == 5:
            b += 1
    return O, w, b


def shape(n):
    """Class-labelled exponent shape: sorted tuple of (class mod 6, exponent) over prime factors."""
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d % 6, e))
        d += 1
    if n > 1:
        out.append((n % 6, 1))
    return tuple(sorted(out))


# positive control: all n <= 30000 plus 3000 pseudo-random large n
rng = np.random.RandomState(20260917)
audit = list(range(1, 30001)) + [int(x) for x in rng.randint(30001, NMAX + 1, size=3000)]
for n in audit:
    O, w, b = trial(n)
    check((int(OM[n]), int(OMD[n]), int(BB[n])) == (O, w, b), "trial division mismatch at n=%d" % n)
print("trial-division audit: %d integers agree with the sieve  %s" % (len(audit), elapsed()))

# class-parity law b(6k-1) odd, b(6k+1) even (inherited SW4) -- re-checked on the full universe
kk = np.arange(1, KMAX + 1, dtype=np.int32)
BASE = 6 * kk                       # int32 centers 6k (6*10^7 fits int32)
check(bool(np.all(BB[BASE - 1] % 2 == 1)), "b(6k-1) odd fails")
check(bool(np.all(BB[BASE + 1] % 2 == 0)), "b(6k+1) even fails")
print("class-parity law b(6k-1) odd, b(6k+1) even: holds for all k <= %d (inherited SW4)" % KMAX)

CLS = np.minimum(OM, 4).astype(np.uint8)   # Omega classes 1,2,3,>=4 (0 for n=1)
LAB = {1: "P", 2: "S", 3: "C", 4: "Q"}     # Q = Omega>=4


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def pair_classes(c, g, K):
    """Omega classes of (6k+c, 6k+c+g) for k=1..K (uint8 arrays)."""
    return CLS[BASE[:K] + c], CLS[BASE[:K] + (c + g)]


def matrix4(left, right, K):
    """4x4 count matrix of (class(left), class(right)) over first K centers."""
    idx = (left[:K].astype(np.int32) - 1) * 4 + (right[:K].astype(np.int32) - 1)
    m = np.bincount(idx, minlength=16).reshape(4, 4)
    return [[int(m[i][j]) for j in range(4)] for i in range(4)]


def marginals(M):
    R = [sum(M[i]) for i in range(4)]
    C = [sum(M[i][j] for i in range(4)) for j in range(4)]
    return R, C


def print_matrix(M, rowlab="6k-1", collab="6k+1"):
    print("   rows: Omega(%s) in 1,2,3,>=4 ; cols: Omega(%s) in 1,2,3,>=4" % (rowlab, collab))
    print("   %10s %10s %10s %10s | %10s" % ("1", "2", "3", ">=4", "row sum"))
    for i in range(4):
        print("   %10d %10d %10d %10d | %10d" % (M[i][0], M[i][1], M[i][2], M[i][3], sum(M[i])))
    R, cs = marginals(M)
    print("   %10d %10d %10d %10d | %10d  (col sums)" % (cs[0], cs[1], cs[2], cs[3], sum(cs)))


def antisym_report(M, K, tag):
    """Actual antisymmetric part, independence prediction, residual (exact Fractions)."""
    R, C = marginals(M)
    check(sum(R) == K and sum(C) == K, "marginals must sum to K")
    print("   %s  pairs (i,j): actual N_ij-N_ji | independence (R_i C_j - R_j C_i)/K | residual" % tag)
    out = {}
    for (i, j) in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
        act = M[i][j] - M[j][i]
        pred = Fraction(R[i] * C[j] - R[j] * C[i], K)
        res = act - pred
        out[(i, j)] = (act, pred, res)
        print("     N_%s%s - N_%s%s = %7d | pred = %10.2f | residual = %10.2f"
              % (LAB[i + 1], LAB[j + 1], LAB[j + 1], LAB[i + 1], act, float(pred), float(res)))
    return R, C, out


def running_sign(left, right, a, b, K, tag):
    """Running N_ab(K') - N_ba(K') for K' <= K; returns summary dict.
    'sign changes' = transitions between + and - in the nonzero-sign subsequence (zeros skipped)."""
    d = np.cumsum((left[:K] == a) & (right[:K] == b), dtype=np.int32)
    d -= np.cumsum((left[:K] == b) & (right[:K] == a), dtype=np.int32)
    npos = int((d > 0).sum())
    nzero = int((d == 0).sum())
    nneg = int((d < 0).sum())
    s = np.sign(d).astype(np.int8)
    nzmask = s != 0
    sv = s[nzmask]
    ch = np.flatnonzero(sv[1:] != sv[:-1])            # positions in sv after which the sign flips
    changes = []
    if len(ch):
        cs = np.cumsum(nzmask, dtype=np.int32)         # cs[k-1] = # nonzero entries among the first k
        # the (m+1)-th nonzero entry (0-based m) sits at the first k with cs[k-1] == m+1
        first_k = int(np.searchsorted(cs, int(ch[0]) + 2)) + 1
        last_k = int(np.searchsorted(cs, int(ch[-1]) + 2)) + 1
        del cs
        changes = (len(ch), first_k, last_k)
    del s, nzmask, sv, ch
    nonpos = d <= 0
    lastnonpos = K - int(np.argmax(nonpos[::-1])) if nonpos.any() else 0
    del nonpos
    nonneg = d >= 0
    lastnonneg = K - int(np.argmax(nonneg[::-1])) if nonneg.any() else 0
    del nonneg
    minv = int(d.min()); maxv = int(d.max())
    neg = d < 0
    first_neg = int(np.argmax(neg)) + 1 if neg.any() else 0
    del neg
    nch, firstc, lastc = changes if changes else (0, "none", "none")
    print("   %s: K'<=%d with diff>0: %d, =0: %d, <0: %d ; sign changes: %d (first at K'=%s, last at K'=%s)"
          % (tag, K, npos, nzero, nneg, nch, firstc, lastc))
    print("        last K' with diff<=0: %d ; last K' with diff>=0: %d ; running min %d ; running max %d ; final %d ; first K' with diff<0: %s"
          % (lastnonpos, lastnonneg, minv, maxv, int(d[-1]), first_neg if first_neg else "none"))
    track = {Kt: int(d[Kt - 1]) for Kt in K_TRACK if Kt <= K}
    return dict(npos=npos, nzero=nzero, nneg=nneg, nchanges=nch, lastnonpos=lastnonpos,
                lastnonneg=lastnonneg, final=int(d[-1]), track=track, minv=minv, maxv=maxv,
                first_neg=first_neg)


# --------------------------------------------------------------------------
# 1. 4x4 matrices and running signs for (6k-1, 6k+1)
# --------------------------------------------------------------------------
hr("1. ORDERED SANDWICH MATRICES N_ij(K), pair (6k-1, 6k+1)")
L0, R0 = pair_classes(-1, 2, KMAX)
MAT = {}
ANTI = {}
for K in K_SCALES:
    M = matrix4(L0, R0, K)
    MAT[K] = M
    print()
    print(" K = %d  (centers 6..%d)" % (K, 6 * K))
    print_matrix(M)
    if K == 10**6:
        check(M == INHERITED_K1E6, "K=10^6 matrix disagrees with inherited SW3 table")
        print("   cross-check: identical to the inherited K=10^6 table (divisors note, SW3)  [PASS]")
    ANTI[K] = antisym_report(M, K, "K=%d" % K)

print()
print(" Running sign of the three signed discrepancies, K' = 1..%d" % KMAX)
print("   definition: 'sign changes' counts transitions between + and - in the nonzero-sign")
print("   subsequence of the running difference (values equal to 0 are skipped).")
RS = {}
RS["PS"] = running_sign(L0, R0, 1, 2, KMAX, "N_PS-N_SP")
RS["PC"] = running_sign(L0, R0, 1, 3, KMAX, "N_PC-N_CP")
RS["SC"] = running_sign(L0, R0, 2, 3, KMAX, "N_SC-N_CS")
check(RS["PS"]["first_neg"] == 1453, "minimal witness of N_PS<N_SP must be K=1453")
# independent confirmation of the minimal witness by direct trial division of the 2*1453 endpoints
nps = nsp = 0
first_wit = None
for k in range(1, 1454):
    Ol = trial(6 * k - 1)[0]; Or = trial(6 * k + 1)[0]
    if Ol == 1 and Or == 2:
        nps += 1
    if Ol == 2 and Or == 1:
        nsp += 1
    if nps < nsp and first_wit is None:
        first_wit = (k, nps, nsp)
check(first_wit == (1453, 245, 246), "trial-division witness: %s" % (first_wit,))
print("   minimal witness re-verified by trial division: K=1453, N_PS=245, N_SP=246  [PASS]")
print()
print("   tracked values of the three differences:")
print("   %10s %10s %10s %10s" % ("K", "PS-SP", "PC-CP", "SC-CS"))
for Kt in K_TRACK:
    print("   %10d %10d %10d %10d" % (Kt, RS["PS"]["track"][Kt], RS["PC"]["track"][Kt], RS["SC"]["track"][Kt]))
neg_tracked = [Kt for Kt in K_TRACK if RS["PS"]["track"][Kt] < 0]
print("   tracked K with N_PS-N_SP < 0: %s" % (neg_tracked,))
print(elapsed())

# --------------------------------------------------------------------------
# 2. Class decomposition
# --------------------------------------------------------------------------
hr("2. CLASS DECOMPOSITION: primes, semiprimes, 3-almost-primes by class mod 6")


def coprime6_upto(x):
    """n<=x with gcd(n,6)=1 and n>=5, as int32 (two residue classes)."""
    n1 = np.arange(7, x + 1, 6, dtype=np.int32)   # 1 mod 6, n>=7
    n5 = np.arange(5, x + 1, 6, dtype=np.int32)   # 5 mod 6
    return n1, n5


def class_counts(x):
    """Counts of n<=x, gcd(n,6)=1, n>=5, by (Omega class 1..4, residue class 1 or 5)."""
    n1, n5 = coprime6_upto(x)
    c1 = np.bincount(CLS[n1], minlength=5)
    c5 = np.bincount(CLS[n5], minlength=5)
    return [int(c1[j]) for j in range(1, 5)], [int(c5[j]) for j in range(1, 5)]


def pattern_counts(x, target_omega):
    """(omega, b) pattern counts of n<=x coprime to 6 with Omega=target_omega."""
    n1, n5 = coprime6_upto(x)
    cnt = {}
    tot = 0
    for arr in (n1, n5):
        sel = arr[OM[arr] == target_omega]
        tot += len(sel)
        key = OMD[sel].astype(np.int32) * 8 + BB[sel].astype(np.int32)
        bc = np.bincount(key, minlength=64)
        for w in range(1, 4):
            for b in range(0, 4):
                cnt[(w, b)] = cnt.get((w, b), 0) + int(bc[w * 8 + b])
    check(sum(cnt.values()) == tot, "pattern partition, Omega=%d" % target_omega)
    return cnt


def semiprime_patterns(x):
    c = pattern_counts(x, 2)
    return {
        "(1,1) distinct": c[(2, 0)],
        "(5,5) distinct": c[(2, 2)],
        "(1,5) mixed": c[(2, 1)],
        "p^2, p=1 mod 6": c[(1, 0)],
        "p^2, p=5 mod 6": c[(1, 2)],
    }


CUBE_NAMES = {
    (1, 0): "p1^3            [right,b=0]",
    (1, 3): "p5^3            [left, b=3]",
    (2, 0): "p1^2 q1         [right,b=0]",
    (2, 1): "p1^2 q5         [left, b=1]",
    (2, 2): "p5^2 q1         [right,b=2]",
    (2, 3): "p5^2 q5         [left, b=3]",
    (3, 0): "p1 q1 r1        [right,b=0]",
    (3, 1): "p1 q1 r5        [left, b=1]",
    (3, 2): "p1 q5 r5        [right,b=2]",
    (3, 3): "p5 q5 r5        [left, b=3]",
}


def cube_patterns(x):
    c = pattern_counts(x, 3)
    out = {nm: c[key] for key, nm in CUBE_NAMES.items()}
    check(sum(out.values()) == sum(c.values()), "3-almost-prime ten-pattern partition")
    return out


P5 = PRIMES[PRIMES >= 5]
CHI = np.where(P5 % 6 == 1, 1, -1).astype(np.int64)
CUMCHI = np.concatenate([[0], np.cumsum(CHI)])          # CUMCHI[i] = sum chi over first i primes>=5


def pi_chi(y):
    return int(CUMCHI[np.searchsorted(P5, y, side="right")])


def pi_chi_identity(x):
    """PROVED identity: S_1(x)-S_5(x) = (1/2)[ sum_{5<=p<=x/5} chi(p) pi_chi(x/p) + pi'(sqrt x) ],
    pi_chi(y)=#{5<=p<=y, p=1 mod 6} - #{5<=p<=y, p=5 mod 6}, pi'(y)=#{5<=p<=y}.
    Returns (rhs, ordered double sum O, pi'(sqrt x), primes<=x/5, terms, partial sums)."""
    Pm = P5[P5 <= x // 5]
    y = x // Pm
    idx = np.searchsorted(P5, y, side="right")
    pichi_y = CUMCHI[idx]
    terms = np.where(Pm % 6 == 1, 1, -1).astype(np.int64) * pichi_y
    ordered = int(terms.sum())
    r = math.isqrt(x)
    piprime_sqrt = int((P5 <= r).sum())
    check((ordered + piprime_sqrt) % 2 == 0, "identity parity")
    rhs = (ordered + piprime_sqrt) // 2
    partial = np.cumsum(terms)
    return rhs, ordered, piprime_sqrt, Pm, terms, partial


# small-x brute-force check of the identity (all x >= 1; both sides vanish for x < 25)
def s1_minus_s5_brute(x):
    tot = 0
    for n in range(5, x + 1):
        if n % 6 in (1, 5) and trial(n)[0] == 2:
            tot += 1 if n % 6 == 1 else -1
    return tot


for xs in [1, 24, 25, 26, 35, 49, 100, 1000, 10007]:
    rhs_s = pi_chi_identity(xs)[0] if xs >= 25 else 0
    if xs < 25:
        # no prime p with 5<=p<=x/5 and pi'(sqrt x)=0: rhs is 0 by definition
        check(int((P5 <= math.isqrt(xs)).sum()) == 0, "pi'(sqrt x) must vanish for x<25")
    lhs_s = s1_minus_s5_brute(xs)
    check(lhs_s == rhs_s, "identity brute force at x=%d: %d vs %d" % (xs, lhs_s, rhs_s))
print(" identity S_1-S_5 = (O + pi'(sqrt x))/2 brute-forced at x = 1,24,25,26,35,49,100,1000,10007"
      " (x=10007: %d)  [PASS]" % s1_minus_s5_brute(10007))

CLASSDATA = {}
IDENT = {}
for K in K_SCALES:
    x = 6 * K + 1
    c1, c5 = class_counts(x)
    CLASSDATA[K] = (c1, c5)
    print()
    print(" x = 6K+1 = %d  (K=%d)" % (x, K))
    print("   Omega class :      1 (P)      2 (S)      3 (C)    >=4 (Q)")
    print("   n = 1 mod 6 : %10d %10d %10d %10d   (right endpoints, column sums)" % tuple(c1))
    print("   n = 5 mod 6 : %10d %10d %10d %10d   (left  endpoints, row sums)" % tuple(c5))
    print("   chi-sum     : %10d %10d %10d %10d   (= #class1 - #class5 ; sign pattern by Omega)"
          % tuple(c1[j] - c5[j] for j in range(4)))
    M = MAT[K]
    Rr, Cc = marginals(M)
    check(Rr == c5 and Cc == c1, "matrix marginals must equal class counts at K=%d" % K)
    print("   check: matrix row sums == class-5 counts, column sums == class-1 counts  [PASS]")
    sp = semiprime_patterns(x)
    print("   semiprime patterns (n<=x, gcd(n,6)=1):")
    for nm, c in sp.items():
        print("      %-18s %10d" % (nm, c))
    same = sp["(1,1) distinct"] + sp["(5,5) distinct"] + sp["p^2, p=1 mod 6"] + sp["p^2, p=5 mod 6"]
    mixed = sp["(1,5) mixed"]
    check(same == c1[1] and mixed == c5[1], "class-parity law for semiprimes (b even <-> same class)")
    print("      same-class (incl. squares) = %d = S_1 ; mixed = %d = S_5   [inherited SW4 law, PASS]" % (same, mixed))
    rhs, ordered, pisq, Pm, terms, partial = pi_chi_identity(x)
    lhs = c1[1] - c5[1]
    check(lhs == rhs, "PROVED identity S_1-S_5 = (ordered + pi'(sqrt x))/2 fails: %d vs %d" % (lhs, rhs))
    distinct_diff = sp["(1,1) distinct"] + sp["(5,5) distinct"] - sp["(1,5) mixed"]
    check(sp["p^2, p=1 mod 6"] + sp["p^2, p=5 mod 6"] == pisq, "#squares = pi'(sqrt x)")
    check(2 * distinct_diff == ordered - pisq, "distinct-pair difference = (O - pi'(sqrt x))/2")
    IDENT[K] = dict(O=ordered, pisq=pisq, distinct=distinct_diff)
    print("   PROVED identity S_1-S_5 = (1/2)[sum_p chi(p) pi_chi(x/p) + pi'(sqrt x)]:")
    print("      S_1-S_5 = %d ; ordered double sum O = %d ; pi'(sqrt x) = %d ; rhs = %d  [PASS]"
          % (lhs, ordered, pisq, rhs))
    print("      split by convention: Ford-Sneed form  O/2 = %.1f  + squares pi'/2 = %.1f ;"
          % (ordered / 2, pisq / 2))
    print("                           squares-in-full  pi' = %d  + distinct pairs (O-pi')/2 = %d"
          " (= same-class distinct - mixed = %d)  [PASS]" % (pisq, (ordered - pisq) // 2, distinct_diff))
    print("      partial ordered sums through p = 5,7,11,13,17,19,23,29,31,37 and p<=100, p<=1000:")
    for pcut in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 100, 1000]:
        i = int(np.searchsorted(Pm, pcut, side="right"))
        print("         p<=%-5d : %8d" % (pcut, int(partial[i - 1]) if i > 0 else 0))
    for yy in [x, x // 5, x // 7, x // 11, x // 13]:
        print("      pi_chi(%d) = %d" % (yy, pi_chi(yy)))

# class-labelled exponent shape determination (Lemma 2.1) checked by full factorization
print()
print(" Class-labelled exponent shape vs (Omega, omega, b), all n <= 300000 coprime to 6 with Omega <= 4:")
SHAPES = {}
WITNESS = {}
n1s, n5s = coprime6_upto(300000)
for arr in (n1s, n5s):
    sel = arr[(OM[arr] >= 1) & (OM[arr] <= 4)]
    for n in sel.tolist():
        key = (int(OM[n]), int(OMD[n]), int(BB[n]))
        sh = shape(n)
        d = SHAPES.setdefault(key, {})
        if sh not in d:
            d[sh] = n
for Om in (1, 2, 3, 4):
    keys = sorted(k for k in SHAPES if k[0] == Om)
    multi = [k for k in keys if len(SHAPES[k]) > 1]
    print("   Omega=%d: %d keys (Omega,omega,b); keys with more than one shape: %d" % (Om, len(keys), len(multi)))
    if Om <= 3:
        check(len(multi) == 0, "shape must be determined for Omega=%d" % Om)
    for k in multi:
        items = sorted(SHAPES[k].items(), key=lambda t: t[1])
        print("      key %s : " % (k,) + " | ".join("%s (n=%d)" % (sh, n) for sh, n in items))
check(sum(1 for k in SHAPES if 1 <= k[0] <= 3) == 17, "17 keys for 1<=Omega<=3")
check(len(SHAPES[(4, 3, 2)]) == 2 and 2275 in SHAPES[(4, 3, 2)].values() and 2695 in SHAPES[(4, 3, 2)].values(),
      "Omega=4 collision (4,3,2) witnesses 2275, 2695")
check(len(SHAPES[(4, 2, 0)]) == 2 and 4459 in SHAPES[(4, 2, 0)].values() and 8281 in SHAPES[(4, 2, 0)].values(),
      "Omega=4 collision (4,2,0) witnesses 4459, 8281")
check(len(SHAPES[(4, 2, 4)]) == 2 and 1375 in SHAPES[(4, 2, 4)].values() and 3025 in SHAPES[(4, 2, 4)].values(),
      "Omega=4 collision (4,2,4) witnesses 1375, 3025")
# bare class multiset is always {1^(Omega-b), 5^b}: tautological from the definition of b; spot-print only
print("   (the bare class multiset is {1^(Omega-b), 5^b} for every n coprime to 6: determined by (Omega,b) alone)")
print("   Omega<=3: every key has exactly one shape (17 keys)  [PASS] ; Omega=4 minimal collisions: "
      "(4,3,2) 2275=5^2*7*13 vs 2695=5*7^2*11 ; (4,2,0) 4459=7^3*13 vs 8281=7^2*13^2 ; "
      "(4,2,4) 1375=5^3*11 vs 3025=5^2*11^2  [PASS]")

print()
print(" Independence (product-of-marginals) prediction for N_PS - N_SP and the two-mechanism split")
print("   pred = (R_1 C_2 - R_2 C_1)/K = 2(pibar*dS + Sbar*dpi)/K,  pibar=(R_1+C_1)/2, dpi=(R_1-C_1)/2,")
print("   Sbar=(R_2+C_2)/2, dS=(C_2-R_2)/2.   R = class-5 (left) counts, C = class-1 (right) counts.")
print("   %8s %8s %8s %8s %8s %10s %10s %10s %10s %10s %8s"
      % ("K", "pi_5", "pi_1", "S_5", "S_1", "actual", "pred", "semi-part", "prime-part", "residual", "act/pred"))
MECH = {}
for K in K_SCALES:
    c1, c5 = CLASSDATA[K]
    R1, C1, R2, C2 = c5[0], c1[0], c5[1], c1[1]
    pibar = Fraction(R1 + C1, 2); dpi = Fraction(R1 - C1, 2)
    Sbar = Fraction(R2 + C2, 2); dS = Fraction(C2 - R2, 2)
    pred = Fraction(R1 * C2 - R2 * C1, K)
    semi = Fraction(2) * pibar * dS / K
    prim = Fraction(2) * Sbar * dpi / K
    check(pred == semi + prim, "two-mechanism split must be exact")
    act = MAT[K][0][1] - MAT[K][1][0]
    O = IDENT[K]["O"]; pisq = IDENT[K]["pisq"]
    check(4 * dS == O + pisq, "dS = (O + pi'(sqrt x))/4")
    sq_fs = pibar * Fraction(pisq, 2 * K)          # Ford-Sneed form: squares pi'/2
    cross_fs = pibar * Fraction(O, 2 * K)          # Ford-Sneed form: O/2
    sq_full = pibar * Fraction(pisq, K)            # squares in full: pi'
    distinct = pibar * Fraction(O - pisq, 2 * K)   # distinct pairs (O-pi')/2
    check(semi == sq_fs + cross_fs == sq_full + distinct, "three-way splits must be exact")
    MECH[K] = dict(act=act, pred=pred, semi=semi, prim=prim, res=act - pred, sq_fs=sq_fs, cross_fs=cross_fs,
                   sq_full=sq_full, distinct=distinct)
    print("   %8d %8d %8d %8d %8d %10d %10.2f %10.2f %10.2f %10.2f %8.3f"
          % (K, R1, C1, R2, C2, act, float(pred), float(semi), float(prim), float(act - pred),
             act / float(pred) if pred != 0 else float('nan')))
print("   three-way split of pred (exact; the semiprime part is split by Theorem 2.2 in two conventions):")
print("   %8s %10s | %10s %10s %10s | %10s %10s" % ("K", "prime", "FS-square", "FS-cross", "FS-sq %", "full-sq", "distinct"))
print("   %8s %10s | %10s %10s %10s | %10s %10s" % ("", "2Sbar dpi/K", "pibar pi'/2K", "pibar O/2K", "of pred", "pibar pi'/K", "pibar(O-pi')/2K"))
for K in K_SCALES:
    m = MECH[K]
    print("   %8d %10.2f | %10.2f %10.2f %9.1f%% | %10.2f %10.2f   (full-sq = %5.1f%% of pred, prime = %4.1f%%)"
          % (K, float(m["prim"]), float(m["sq_fs"]), float(m["cross_fs"]), 100 * float(m["sq_fs"] / m["pred"]),
             float(m["sq_full"]), float(m["distinct"]), 100 * float(m["sq_full"] / m["pred"]),
             100 * float(m["prim"] / m["pred"])))
print("   ratio form: N_PS/N_SP actual vs (pi_5 S_1)/(pi_1 S_5):")
for K in K_SCALES:
    c1, c5 = CLASSDATA[K]
    a = Fraction(MAT[K][0][1], MAT[K][1][0])
    p = Fraction(c5[0] * c1[1], c1[0] * c5[1])
    print("      K=%d : actual %.6f  predicted %.6f  (pi_5/pi_1 = %.6f, S_1/S_5 = %.6f)"
          % (K, float(a), float(p), c5[0] / c1[0], c1[1] / c5[1]))
print("   heuristic size of pred against two candidate scales (x = 6K+1):")
for K in K_SCALES:
    x = 6 * K + 1
    s1 = math.sqrt(x) / math.log(x)
    s2 = math.sqrt(x) * math.log(math.log(x)) / math.log(x) ** 2
    print("      K=%d : pred/(sqrt x/log x) = %.2f ; pred/(sqrt x loglog x/log^2 x) = %.2f"
          % (K, float(MECH[K]["pred"]) / s1, float(MECH[K]["pred"]) / s2))
print(elapsed())

# --------------------------------------------------------------------------
# 3. 3-almost-prime patterns and the (P,C),(S,C) pairs
# --------------------------------------------------------------------------
hr("3. 3-ALMOST-PRIME CLASS PATTERNS (omega, b) AND THE (P,C)/(C,P), (S,C)/(C,S) DISCREPANCIES")
for K in K_SCALES:
    x = 6 * K + 1
    cp = cube_patterns(x)
    c1, c5 = CLASSDATA[K]
    print()
    print(" x = %d (K=%d): 3-almost-primes coprime to 6 by pattern" % (x, K))
    left = right = 0
    for nm, c in cp.items():
        print("      %-32s %10d" % (nm, c))
        if "left" in nm:
            left += c
        else:
            right += c
    check(left == c5[2] and right == c1[2], "3-almost class-parity partition")
    print("      left total (b odd) = %d = C_5 ; right total (b even) = %d = C_1   [PASS]" % (left, right))
    print("   two-mechanism split for the other two pairs (independence prediction, exact):")
    for (i, j, nm) in [(0, 2, "N_PC-N_CP"), (1, 2, "N_SC-N_CS")]:
        Ri, Ci, Rj, Cj = c5[i], c1[i], c5[j], c1[j]
        pred = Fraction(Ri * Cj - Rj * Ci, K)
        ibar = Fraction(Ri + Ci, 2); di = Fraction(Ri - Ci, 2)
        jbar = Fraction(Rj + Cj, 2); dj = Fraction(Cj - Rj, 2)
        part_j = Fraction(2) * ibar * dj / K   # driven by class bias of Omega=j+1 numbers
        part_i = Fraction(2) * jbar * di / K   # driven by class bias of Omega=i+1 numbers
        check(pred == part_i + part_j, "split")
        act = MAT[K][i][j] - MAT[K][j][i]
        print("      %s: actual %7d  pred %9.2f  [Omega=%d-bias part %9.2f, Omega=%d-bias part %9.2f]  residual %9.2f"
              % (nm, act, float(pred), j + 1, float(part_j), i + 1, float(part_i), float(act - pred)))

# --------------------------------------------------------------------------
# 4. Hostile controls
# --------------------------------------------------------------------------
hr("4. HOSTILE CONTROLS")
print(" (i) shifted pairs.  H1 predicts: class-reversed pair (6k+1,6k+5) -> sign of N_PS-N_SP REVERSED;")
print("     same-class pairs (6k-1,6k+5) and (6k+1,6k+7) -> no class-driven term (pred=0 exactly up to")
print("     the marginal mismatch between the two shifted copies), residual only.")
PAIRS = {
    "(6k-1,6k+1) [5|1]": (-1, 2),
    "(6k+1,6k+5) [1|5]": (1, 4),
    "(6k-1,6k+5) [5|5]": (-1, 6),
    "(6k+1,6k+7) [1|1]": (1, 6),
}
CTRL = {}
for nm, (c, g) in PAIRS.items():
    La, Rb = pair_classes(c, g, KMAX)
    print()
    print(" pair %s" % nm)
    for K in K_SCALES:
        M = matrix4(La, Rb, K)
        R, C = marginals(M)
        rows = []
        for (i, j, lab) in [(0, 1, "PS-SP"), (0, 2, "PC-CP"), (1, 2, "SC-CS")]:
            act = M[i][j] - M[j][i]
            pred = Fraction(R[i] * C[j] - R[j] * C[i], K)
            rows.append((lab, act, pred, act - pred))
        CTRL[(nm, K)] = dict(M=M, rows=rows)
        print("   K=%d  PP=%d PS=%d SP=%d SS=%d | " % (K, M[0][0], M[0][1], M[1][0], M[1][1]) +
              " ; ".join("%s: act %6d pred %9.2f res %9.2f" % (lab, act, float(pred), float(res))
                         for (lab, act, pred, res) in rows))
    rs = running_sign(La, Rb, 1, 2, KMAX, "   running N_PS-N_SP")
    CTRL[(nm, "rs")] = rs
    del La, Rb
check(CTRL[("(6k-1,6k+1) [5|1]", 10**7)]["M"] == MAT[10**7], "pair (i) row 1 must equal section 1")

print()
print(" (ii) the nine mod-35 strata of the original pair (6k-1,6k+1).")
print("     stratum = (5-state, 7-state) with each state in {none, left, right}: which endpoint (if any)")
print("     is divisible by 5, resp. by 7 (never both endpoints, since gcd(6k-1,6k+1)=1).")
print("     Symmetric strata: 'none' and the two 'cross' strata (one prime each side).")
print("     H1 predicts the sign persists in symmetric strata with the independence magnitude computed")
print("     from the stratum's own marginals; a purely local-sieve explanation would make it vanish there.")
d5m = ((BASE - 1) % 5 == 0); d5p = ((BASE + 1) % 5 == 0)
d7m = ((BASE - 1) % 7 == 0); d7p = ((BASE + 1) % 7 == 0)
STRATA = [
    ("none  (no 5, no 7)", ~(d5m | d5p | d7m | d7p), 15),
    ("cross (5|left,7|right)", d5m & d7p, 1),
    ("cross (7|left,5|right)", d7m & d5p, 1),
    ("left-only (5|left, no 7)", d5m & ~d7m & ~d7p, 5),
    ("right-only (5|right, no 7)", d5p & ~d7m & ~d7p, 5),
    ("left-only (7|left, no 5)", d7m & ~d5m & ~d5p, 3),
    ("right-only (7|right, no 5)", d7p & ~d5m & ~d5p, 3),
    ("35|left  (5,7|left)", d5m & d7m, 1),
    ("35|right (5,7|right)", d5p & d7p, 1),
]
del d5m, d5p, d7m, d7p
cover = np.zeros(KMAX, dtype=np.uint8)
for nm, mask, sz in STRATA:
    cover += mask
check(bool(np.all(cover == 1)), "the nine strata must partition the centers")
del cover
sizes = {nm: int(mask[:35].sum()) for nm, mask, sz in STRATA}
check(all(sizes[nm] == sz for nm, mask, sz in STRATA) and sum(sizes.values()) == 35,
      "CRT stratum sizes on one block of 35: %s" % sizes)
print("   partition check: the nine strata cover every center exactly once; on one block of 35 centers the")
print("   sizes are none=15, cross=1+1, one-sided 5,5,3,3, 35|left=1, 35|right=1 (sum 35)  [PASS]")
for K in K_SCALES:
    print()
    print("   K = %d" % K)
    tot_all_act = 0
    tot_sym_act = 0
    tot_sym_pred = Fraction(0)
    tot_seven_act = 0
    none_left_primes = 0
    for nm, mask, sz in STRATA:
        m = mask[:K]
        n = int(m.sum())
        Ls = L0[:K][m]; Rs = R0[:K][m]
        M = matrix4(Ls, Rs, n)
        R, C = marginals(M)
        act = M[0][1] - M[1][0]
        pred = Fraction(R[0] * C[1] - R[1] * C[0], n) if n else Fraction(0)
        tot_all_act += act
        if nm.startswith("none") or nm.startswith("cross"):
            tot_sym_act += act
            tot_sym_pred += pred
        if not nm.startswith("35|"):
            tot_seven_act += act
        if nm.startswith("none"):
            none_left_primes = R[0]
            check(M[0][0] == MAT[K][0][0] - 1, "all PP pairs except (5,7) lie in the 'none' stratum")
        if nm == "cross (5|left,7|right)":
            check(M[0][0] == 1 and M[0][1] == 0 and M[1][0] == 0, "cross stratum: exactly the PP pair (5,7) at k=1")
        if nm == "cross (7|left,5|right)":
            check(M[0][0] == 0 and M[0][1] == 0 and M[1][0] == 0, "cross stratum (7|left,5|right) has no P endpoints")
        if nm.startswith("35|left"):
            check(M[0][0] == 0 and M[0][1] == 0 and M[1][0] == 1, "35|left stratum: exactly the SP pair (35,37) at k=6")
        if nm.startswith("35|right"):
            check(M[0][0] == 0 and M[0][1] == 0 and M[1][0] == 0, "35|right stratum: no P endpoints (35m = 6k+1 has m>1)")
        print("     %-28s n=%8d  PP=%6d PS=%6d SP=%6d SS=%6d  PS-SP act %6d pred %8.2f res %8.2f"
              % (nm, n, M[0][0], M[0][1], M[1][0], M[1][1], act, float(pred), float(act - pred)))
    total = MAT[K][0][1] - MAT[K][1][0]
    check(tot_all_act == total, "nine strata must sum to the total N_PS-N_SP")
    check(tot_seven_act == total + 1, "the seven strata without 35|left/right sum to total+1 (the SP pair (35,37))")
    print("     nine strata sum: PS-SP = %d = total  [PASS] ; the seven strata without 35|left/35|right sum to %d"
          " (the missing -1 is the SP pair (35,37) at k=6)" % (tot_all_act, tot_seven_act))
    print("     symmetric strata (none + cross) combined: act %d  pred %.2f  res %.2f" %
          (tot_sym_act, float(tot_sym_pred), float(tot_sym_act - tot_sym_pred)))
    print("     left-endpoint primes in the 'none' stratum: %d of %d (%.1f%%); PP pairs there: %d of %d"
          " (all but (5,7))" % (none_left_primes, sum(MAT[K][0]), 100.0 * none_left_primes / sum(MAT[K][0]),
                                MAT[K][0][0] - 1, MAT[K][0][0]))
del STRATA

print()
print(" (iii) prime-only Chebyshev census pi(x;6,5)-pi(x;6,1) and the running sign of the class-1 vs class-5")
print("      prime race up to %d, to place the sandwich sign next to the classical bias." % NMAX)
cum = CUMCHI[1:]
zeros_at = P5[np.nonzero(cum == 0)[0]]
print("      primes >= 5 up to %d: %d ; pi_chi at end = %d ; #prefixes with pi_chi<0: %d, =0: %d, >0: %d"
      % (NMAX, len(P5), int(cum[-1]), int((cum < 0).sum()), int((cum == 0).sum()), int((cum > 0).sum())))
print("      primes at which pi_chi = 0: %s" % (zeros_at.tolist(),))
pos = np.nonzero(cum > 0)[0]
print("      largest prime at which pi_chi > 0 (class 1 ahead): %s" % (int(P5[pos[-1]]) if len(pos) else "never"))
print("      pi_chi(6*10^7+1) = %d" % pi_chi(6 * 10**7 + 1))
print(elapsed())

# --------------------------------------------------------------------------
# 4(iv). Pooled-gap control with block-level residual test
# --------------------------------------------------------------------------
print()
print(" (iv) pooled-gap control.  For each orientation [class n | class n+g] the 50 gaps g in 2..300 with both")
print("      endpoints coprime to 6: [5|1] n=6k-1, g=2 mod 6 ; [1|5] n=6k+1, g=4 mod 6 ; [5|5] n=6k-1, g=0 mod 6 ;")
print("      [1|1] n=6k+1, g=0 mod 6.  D_g = N_ij - N_ji over k<=K ; pred_g = (R_i C_j - R_j C_i)/K from the pair's")
print("      own marginals.  Block analysis: K split into %d blocks (size K/%d); p_(g,b) = block product of"
      % (100, 100))
print("      marginals, D_(g,b) = block discrepancy, r_(g,b) = D_(g,b) - p_(g,b).  Reported per orientation:")
print("      pooled mean/sd(ddof=1)/#D_g>0 of D_g, mean pred_g ; corr and slope of block means m_b=mean_g D_(g,b)")
print("      against p_b=mean_g p_(g,b) ; mean off-diagonal inter-gap correlation of D_(g,b) and of r_(g,b) over")
print("      blocks ; residual test: RES = sum_b mean_g r_(g,b), SE = sd_b(mean_g r_(g,b))*sqrt(#blocks), z=RES/SE ;")
print("      per-gap residual sd = sd_g(D_g - sum_b p_(g,b)).  The K=10^6 blocks are the first 100 blocks of")
print("      10^4 centers; the K=10^7 blocks are 100 blocks of 10^5 centers (merged from 1000 blocks of 10^4).")
ORIENT = [
    ("[5|1]", -1, [g for g in range(2, GMAX + 1) if g % 6 == 2]),
    ("[1|5]", 1, [g for g in range(2, GMAX + 1) if g % 6 == 4]),
    ("[5|5]", -1, [g for g in range(2, GMAX + 1) if g % 6 == 0]),
    ("[1|1]", 1, [g for g in range(2, GMAX + 1) if g % 6 == 0]),
]
for nm, c, gl in ORIENT:
    check(len(gl) == 50, "50 gaps per orientation")
    for g in gl:
        check(all(x % 6 in (1, 5) for x in (6 * 7 + c, 6 * 7 + c + g)), "endpoints coprime to 6")
check(ORIENT[0][2][0] == 2 and ORIENT[1][2][0] == 4 and ORIENT[2][2][0] == 6, "first gaps 2,4,6")
BLK = np.repeat(np.arange(NB_FINE, dtype=np.int32), KMAX // NB_FINE)   # block id of each k (0..999)
check(len(BLK) == KMAX, "block ids")
PAIRIDX = [(0, 1, "PS-SP"), (0, 2, "PC-CP"), (1, 2, "SC-CS")]
POOL = {}   # (orient, K, lab) -> dict of statistics
CNT = {}    # orient -> array [50, 1000, 4, 4] of block counts
for nm, c, gl in ORIENT:
    arr = np.zeros((len(gl), NB_FINE, 4, 4), dtype=np.int64)
    for gi, g in enumerate(gl):
        L, R = pair_classes(c, g, KMAX)
        idx = (L.astype(np.int32) - 1) * 4 + (R.astype(np.int32) - 1)
        cnt = np.bincount(BLK * 16 + idx, minlength=NB_FINE * 16).reshape(NB_FINE, 4, 4)
        arr[gi] = cnt
        del L, R, idx, cnt
    CNT[nm] = arr
# consistency with the single-pair controls
for nm, c, gl in ORIENT:
    key = [k for k in PAIRS if k.endswith(nm)][0]
    g0 = PAIRS[key][1]
    gi = gl.index(g0)
    for K in (10**6, 10**7):
        M = CNT[nm][gi, :K // 10**4].sum(axis=0)
        check([[int(M[i][j]) for j in range(4)] for i in range(4)] == CTRL[(key, K)]["M"],
              "pooled block counts must reproduce the single-pair matrix %s at K=%d" % (key, K))
print("      consistency: the block counts at the four control gaps reproduce the section 4(i) matrices  [PASS]")


def mean_offdiag_corr(X):
    """X: [ngap, nblock]; mean off-diagonal Pearson correlation across gaps (over blocks)."""
    Cm = np.corrcoef(X)
    n = Cm.shape[0]
    return float((Cm.sum() - np.trace(Cm)) / (n * (n - 1)))


for K in (10**6, 10**7):
    nb = 100
    bsz = K // nb
    merge = bsz // 10**4                  # fine blocks per coarse block
    print()
    print("   K = %d  (%d blocks of %d centers)" % (K, nb, bsz))
    for (i, j, lab) in PAIRIDX:
        print("    %s:" % lab)
        print("      %-6s %9s %8s %6s %9s | %6s %6s | %7s %7s | %8s %7s %6s | %9s"
              % ("orient", "mean D_g", "sd", "#>0", "mean pred", "corr", "slope", "cor D", "cor r",
                 "RES", "SE", "z", "gap-res sd"))
        for nm, c, gl in ORIENT:
            A = CNT[nm][:, :K // 10**4].reshape(len(gl), nb, merge, 4, 4).sum(axis=2)   # [50, nb, 4, 4]
            Rb = A.sum(axis=3).astype(np.float64)      # row marginals per block [50, nb, 4]
            Cb = A.sum(axis=2).astype(np.float64)      # col marginals per block [50, nb, 4]
            Dgb = (A[:, :, i, j] - A[:, :, j, i]).astype(np.float64)
            pgb = (Rb[:, :, i] * Cb[:, :, j] - Rb[:, :, j] * Cb[:, :, i]) / bsz
            Dg = Dgb.sum(axis=1)
            Rt = Rb.sum(axis=1); Ct = Cb.sum(axis=1)
            predg = (Rt[:, i] * Ct[:, j] - Rt[:, j] * Ct[:, i]) / K
            mean_D = float(Dg.mean()); sd_D = float(Dg.std(ddof=1)); npos = int((Dg > 0).sum())
            mean_pred = float(predg.mean())
            mb = Dgb.mean(axis=0); pb = pgb.mean(axis=0)
            corr = float(np.corrcoef(mb, pb)[0, 1])
            slope = float(np.polyfit(pb, mb, 1)[0])
            corD = mean_offdiag_corr(Dgb)
            rgb = Dgb - pgb
            corr_r = mean_offdiag_corr(rgb)
            Rbm = rgb.mean(axis=0)                     # per-block residual mean over gaps
            RES = float(Rbm.sum())
            SE = float(Rbm.std(ddof=1) * math.sqrt(nb))
            z = RES / SE if SE > 0 else float("nan")
            gapres = Dg - pgb.sum(axis=1)
            gapres_sd = float(gapres.std(ddof=1))
            POOL[(nm, K, lab)] = dict(mean_D=mean_D, sd_D=sd_D, npos=npos, mean_pred=mean_pred, corr=corr,
                                      slope=slope, corD=corD, corr_r=corr_r, RES=RES, SE=SE, z=z,
                                      gapres_sd=gapres_sd, Dg=Dg, predg=predg, sumblock=pgb.sum(axis=1),
                                      gaps=gl)
            print("      %-6s %9.1f %8.1f %3d/50 %9.1f | %6.2f %6.2f | %7.2f %7.2f | %8.1f %7.1f %6.2f | %9.1f"
                  % (nm, mean_D, sd_D, npos, mean_pred, corr, slope, corD, corr_r, RES, SE, z, gapres_sd))
    # the sandwich gap itself inside its orientation
    ps = POOL[("[5|1]", K, "PS-SP")]
    D2 = ps["Dg"][0]; pred2 = ps["predg"][0]; sb2 = ps["sumblock"][0]
    rank = int((ps["Dg"][1:] > D2).sum())
    check(int(D2) == MAT[K][0][1] - MAT[K][1][0], "g=2 pooled datum equals section 1")
    print("    sandwich gap g=2 within [5|1] PS-SP: D_2 = %d ; global pred %.2f ; sum of block preds %.2f ;"
          " residual vs block preds %.1f (%.2f per-gap residual sd) ; %d of the other 49 gaps have larger D_g"
          % (int(D2), pred2, sb2, D2 - sb2, (D2 - sb2) / ps["gapres_sd"], rank))
    # falsifier bookkeeping for PS-SP: same-class pooled |mean| vs opposite-class pooled |mean|
    m51 = POOL[("[5|1]", K, "PS-SP")]["mean_D"]; m15 = POOL[("[1|5]", K, "PS-SP")]["mean_D"]
    m55 = POOL[("[5|5]", K, "PS-SP")]["mean_D"]; m11 = POOL[("[1|1]", K, "PS-SP")]["mean_D"]
    print("    falsifier check PS-SP: |same-class means| = %.1f, %.1f vs |opposite-class means| = %.1f, %.1f -> %s"
          % (abs(m55), abs(m11), abs(m51), abs(m15),
             "FIRES (a same-class |mean| exceeds an opposite-class |mean|)"
             if max(abs(m55), abs(m11)) > min(abs(m51), abs(m15)) else "does not fire"))
    for lab in ("PS-SP", "SC-CS"):
        for nm in ("[5|1]", "[1|5]"):
            q = POOL[(nm, K, lab)]
            print("    %s %s: pooled mean / mean pred = %.2f (%.0f%% off)"
                  % (lab, nm, q["mean_D"] / q["mean_pred"], 100 * abs(q["mean_D"] / q["mean_pred"] - 1)))
    # z-values in units of block SE, listed for the record
    big = [(lab, nm, POOL[(nm, K, lab)]["z"]) for (i, j, lab) in PAIRIDX for nm, c, gl in ORIENT
           if abs(POOL[(nm, K, lab)]["z"]) >= 2]
    print("    cells with |z| >= 2 in the residual test: %s"
          % (", ".join("%s %s z=%.2f" % t for t in big) if big else "none"))
del CNT, BLK
print(elapsed())

# --------------------------------------------------------------------------
# 5. Verdicts
# --------------------------------------------------------------------------
hr("5. VERDICTS (finite-exact facts; the labels are for the note)")
for K in K_SCALES:
    m = MECH[K]
    print(" K=%d: N_PS-N_SP = %d ; independence pred %.2f (semiprime-class part %.2f, prime-class part %.2f) ; residual %.2f"
          % (K, m["act"], float(m["pred"]), float(m["semi"]), float(m["prim"]), float(m["res"])))
rs = RS["PS"]
print(" running sign K'<=%d: N_PS>N_SP for %d values, = for %d, < for %d; last K' with N_PS<=N_SP: %d ;"
      " minimal witness K=%d" % (KMAX, rs["npos"], rs["nzero"], rs["nneg"], rs["lastnonpos"], rs["first_neg"]))
rev = CTRL[("(6k+1,6k+5) [1|5]", "rs")]
print(" class-reversed pair (6k+1,6k+5): final N_PS-N_SP = %d ; K'<=%d with diff<0: %d, >0: %d"
      % (rev["final"], KMAX, rev["nneg"], rev["npos"]))
for nm in ["(6k-1,6k+5) [5|5]", "(6k+1,6k+7) [1|1]"]:
    r = CTRL[(nm, "rs")]
    print(" same-class pair %s: final N_PS-N_SP = %d ; K'<=%d with diff>0: %d, =0: %d, <0: %d ;"
          " running min %d, running max %d ; extremes over the ten tracked K only: min %d, max %d"
          % (nm, r["final"], KMAX, r["npos"], r["nzero"], r["nneg"], r["minv"], r["maxv"],
             min(r["track"].values()), max(r["track"].values())))
for K in (10**6, 10**7):
    q = ["%s %+.1f (pred %+.1f)" % (nm, POOL[(nm, K, "PS-SP")]["mean_D"], POOL[(nm, K, "PS-SP")]["mean_pred"])
         for nm, c, gl in ORIENT]
    print(" pooled PS-SP at K=%d: %s" % (K, " ; ".join(q)))
print(" total time %s" % elapsed())
print("ALL CHECKS PASSED")
