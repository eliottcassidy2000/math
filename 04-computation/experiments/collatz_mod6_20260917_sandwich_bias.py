#!/usr/bin/env python3
"""
collatz_mod6_20260917_sandwich_bias.py -- signed sandwich discrepancy census.

Lane: sandwich_bias (session collatz-mod6-20260917).

Object.  Centers W = 6k, endpoints 6k-1 (class 5 mod 6) and 6k+1 (class 1
mod 6).  Omega = number of prime factors with multiplicity, omega = number
of distinct prime factors, b = number of prime factors congruent to 5 mod 6
with multiplicity.  chi(n) = +1 if n = 1 mod 6, -1 if n = 5 mod 6
(the nontrivial character mod 3 restricted to (n,6)=1).

Inheritance (read, not re-derived):
  05-knowledge/results/arithmetic_braids_20260917_divisors.md
    SW1 exact CRT sandwich identity, SW2 k->-k side symmetry and class-parity
    law b(6k-1) odd / b(6k+1) even, SW3 K=10^6 4x4 Omega matrix.
This script re-derives the K=10^6 matrix only as a cross-check and raises
if it disagrees with the inherited table.

Universe.  Omega, omega, b sieved exactly on [1, 6*10^7 + 7] (uint8).
Everything load-bearing is an exact integer or Fraction.  Explicit `raise`
is used everywhere so the checks stay active under python -O.

Sections printed:
  0  sieve + independent trial-division audit (positive control)
  1  4x4 Omega matrices at K = 10^5, 10^6, 10^7; antisymmetric part;
     running sign of N_PS - N_SP, N_PC - N_CP, N_SC - N_CS for K <= 10^7
  2  class decomposition: primes / semiprimes / 3-almost-primes by class
     pattern, chi-sums by Omega, PROVED identity for S_1 - S_5 checked
     exactly, independence (product-of-marginals) prediction and the
     two-mechanism split of the predicted discrepancy
  3  3-almost-prime class patterns (omega, b) on each side
  4  hostile controls: shifted pairs (6k+1,6k+5), (6k-1,6k+5), (6k+1,6k+7);
     mod-35 local strata (Chebyshev-free-in-the-local-sense control)
  5  summary verdicts

Run:  python3 04-computation/experiments/collatz_mod6_20260917_sandwich_bias.py
"""
import sys
import time
from fractions import Fraction

import numpy as np

T0 = time.time()
NMAX = 6 * 10**7 + 7          # sieve bound; 6*10^7+7 is the largest endpoint used
KMAX = 10**7                  # largest number of centers
K_SCALES = [10**5, 10**6, 10**7]
K_TRACK = [10**4, 2 * 10**4, 5 * 10**4, 10**5, 2 * 10**5, 5 * 10**5,
           10**6, 2 * 10**6, 5 * 10**6, 10**7]

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
    return "[t=%.1fs]" % (time.time() - T0)


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
PRIMES = np.nonzero(is_p)[0]
del is_p
OM = np.zeros(NMAX + 1, np.uint8)     # Omega
OMD = np.zeros(NMAX + 1, np.uint8)    # omega (distinct)
BB = np.zeros(NMAX + 1, np.uint8)     # b = # prime factors 5 mod 6, with mult.
for p in PRIMES.tolist():
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


# positive control: all n <= 30000 plus 3000 pseudo-random large n
rng = np.random.RandomState(20260917)
audit = list(range(1, 30001)) + [int(x) for x in rng.randint(30001, NMAX + 1, size=3000)]
for n in audit:
    O, w, b = trial(n)
    check((int(OM[n]), int(OMD[n]), int(BB[n])) == (O, w, b), "trial division mismatch at n=%d" % n)
print("trial-division audit: %d integers agree with the sieve  %s" % (len(audit), elapsed()))

# class-parity law b(6k-1) odd, b(6k+1) even (inherited SW4) -- re-checked on the full universe
kk = np.arange(1, KMAX + 1, dtype=np.int64)
check(bool(np.all(BB[6 * kk - 1] % 2 == 1)), "b(6k-1) odd fails")
check(bool(np.all(BB[6 * kk + 1] % 2 == 0)), "b(6k+1) even fails")
print("class-parity law b(6k-1) odd, b(6k+1) even: holds for all k <= %d (inherited SW4)" % KMAX)

CLS = np.minimum(OM, 4).astype(np.uint8)   # Omega classes 1,2,3,>=4 (0 for n=1)
LAB = {1: "P", 2: "S", 3: "C", 4: "Q"}     # Q = Omega>=4


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def matrix4(left, right, K):
    """4x4 count matrix of (class(left), class(right)) over first K centers."""
    idx = (left[:K].astype(np.int64) - 1) * 4 + (right[:K].astype(np.int64) - 1)
    m = np.bincount(idx, minlength=16).reshape(4, 4)
    return [[int(m[i][j]) for j in range(4)] for i in range(4)]


def print_matrix(M, rowlab="6k-1", collab="6k+1"):
    print("   rows: Omega(%s) in 1,2,3,>=4 ; cols: Omega(%s) in 1,2,3,>=4" % (rowlab, collab))
    print("   %10s %10s %10s %10s | %10s" % ("1", "2", "3", ">=4", "row sum"))
    for i in range(4):
        print("   %10d %10d %10d %10d | %10d" % (M[i][0], M[i][1], M[i][2], M[i][3], sum(M[i])))
    cs = [sum(M[i][j] for i in range(4)) for j in range(4)]
    print("   %10d %10d %10d %10d | %10d  (col sums)" % (cs[0], cs[1], cs[2], cs[3], sum(cs)))


def antisym_report(M, K, tag):
    """Actual antisymmetric part, independence prediction, residual (exact Fractions)."""
    R = [sum(M[i]) for i in range(4)]
    C = [sum(M[i][j] for i in range(4)) for j in range(4)]
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
    """Running N_ab(K') - N_ba(K') for K' <= K; returns summary dict."""
    ab = ((left[:K] == a) & (right[:K] == b)).astype(np.int32)
    ba = ((left[:K] == b) & (right[:K] == a)).astype(np.int32)
    d = np.cumsum(ab) - np.cumsum(ba)
    npos = int((d > 0).sum())
    nzero = int((d == 0).sum())
    nneg = int((d < 0).sum())
    s = np.sign(d)
    nz = np.nonzero(s)[0]
    changes = []
    if len(nz):
        sv = s[nz]
        ch = np.nonzero(sv[1:] != sv[:-1])[0]
        changes = [int(nz[c + 1]) + 1 for c in ch]   # K' (1-based) at which new sign first holds
    lastnonpos = np.nonzero(d <= 0)[0]
    lastnonpos = int(lastnonpos[-1]) + 1 if len(lastnonpos) else 0
    lastnonneg = np.nonzero(d >= 0)[0]
    lastnonneg = int(lastnonneg[-1]) + 1 if len(lastnonneg) else 0
    minv = int(d.min()); maxv = int(d.max())
    print("   %s: K'<=%d with diff>0: %d, =0: %d, <0: %d ; sign changes: %d (first at K'=%s, last at K'=%s)"
          % (tag, K, npos, nzero, nneg, len(changes),
             changes[0] if changes else "none", changes[-1] if changes else "none"))
    print("        last K' with diff<=0: %d ; last K' with diff>=0: %d ; min diff %d ; max diff %d ; final %d"
          % (lastnonpos, lastnonneg, minv, maxv, int(d[-1])))
    track = {Kt: int(d[Kt - 1]) for Kt in K_TRACK if Kt <= K}
    return dict(npos=npos, nzero=nzero, nneg=nneg, changes=changes, lastnonpos=lastnonpos,
                lastnonneg=lastnonneg, final=int(d[-1]), track=track)


# --------------------------------------------------------------------------
# 1. 4x4 matrices and running signs for (6k-1, 6k+1)
# --------------------------------------------------------------------------
hr("1. ORDERED SANDWICH MATRICES N_ij(K), pair (6k-1, 6k+1)")
L0 = CLS[6 * kk - 1]
R0 = CLS[6 * kk + 1]
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
RS = {}
RS["PS"] = running_sign(L0, R0, 1, 2, KMAX, "N_PS-N_SP")
RS["PC"] = running_sign(L0, R0, 1, 3, KMAX, "N_PC-N_CP")
RS["SC"] = running_sign(L0, R0, 2, 3, KMAX, "N_SC-N_CS")
print()
print("   tracked values of the three differences:")
print("   %10s %10s %10s %10s" % ("K", "PS-SP", "PC-CP", "SC-CS"))
for Kt in K_TRACK:
    print("   %10d %10d %10d %10d" % (Kt, RS["PS"]["track"][Kt], RS["PC"]["track"][Kt], RS["SC"]["track"][Kt]))
print(elapsed())

# --------------------------------------------------------------------------
# 2. Class decomposition
# --------------------------------------------------------------------------
hr("2. CLASS DECOMPOSITION: primes, semiprimes, 3-almost-primes by class mod 6")


def class_counts(x):
    """Counts of n<=x, gcd(n,6)=1, n>=5, by (Omega class 1..4, residue class 1 or 5)."""
    n1 = np.arange(7, x + 1, 6, dtype=np.int64)   # 1 mod 6, n>=7
    n5 = np.arange(5, x + 1, 6, dtype=np.int64)   # 5 mod 6
    c1 = np.bincount(CLS[n1], minlength=5)
    c5 = np.bincount(CLS[n5], minlength=5)
    return [int(c1[j]) for j in range(1, 5)], [int(c5[j]) for j in range(1, 5)]


def semiprime_patterns(x):
    """Semiprimes n<=x coprime to 6 by pattern: (1,1),(5,5),(1,5),p^2 p=1,p^2 p=5."""
    n = np.arange(5, x + 1, dtype=np.int64)
    n = n[(n % 6 == 1) | (n % 6 == 5)]
    sel = n[OM[n] == 2]
    w = OMD[sel]; b = BB[sel]
    out = {
        "(1,1) distinct": int(((w == 2) & (b == 0)).sum()),
        "(5,5) distinct": int(((w == 2) & (b == 2)).sum()),
        "(1,5) mixed": int(((w == 2) & (b == 1)).sum()),
        "p^2, p=1 mod 6": int(((w == 1) & (b == 0)).sum()),
        "p^2, p=5 mod 6": int(((w == 1) & (b == 2)).sum()),
    }
    check(sum(out.values()) == len(sel), "semiprime pattern partition")
    return out


def cube_patterns(x):
    """3-almost-primes n<=x coprime to 6 by (omega,b) pattern."""
    n = np.arange(5, x + 1, dtype=np.int64)
    n = n[(n % 6 == 1) | (n % 6 == 5)]
    sel = n[OM[n] == 3]
    w = OMD[sel]; b = BB[sel]
    names = {
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
    out = {}
    tot = 0
    for (wi, bi), nm in names.items():
        c = int(((w == wi) & (b == bi)).sum())
        out[nm] = c
        tot += c
    check(tot == len(sel), "3-almost-prime pattern partition")
    return out


def pi_chi_identity(x):
    """PROVED identity: S_1(x)-S_5(x) = (1/2)[ sum_{5<=p<=x/5} chi(p) pi_chi(x/p) + pi'(sqrt x) ],
    pi_chi(y)=#{5<=p<=y, p=1 mod 6} - #{5<=p<=y, p=5 mod 6}, pi'(y)=#{5<=p<=y}.
    Returns (lhs from sieve, rhs from prime list, partial sums by p)."""
    P = PRIMES[PRIMES >= 5]
    chi = np.where(P % 6 == 1, 1, -1).astype(np.int64)
    cum = np.concatenate([[0], np.cumsum(chi)])          # cum[i] = sum chi over first i primes>=5
    Pm = P[P <= x // 5]
    y = x // Pm
    idx = np.searchsorted(P, y, side="right")
    pichi_y = cum[idx]
    terms = np.where(Pm % 6 == 1, 1, -1).astype(np.int64) * pichi_y
    ordered = int(terms.sum())
    r = int(np.sqrt(x))
    while (r + 1) * (r + 1) <= x:
        r += 1
    while r * r > x:
        r -= 1
    piprime_sqrt = int((P <= r).sum())
    check((ordered + piprime_sqrt) % 2 == 0, "identity parity")
    rhs = (ordered + piprime_sqrt) // 2
    partial = np.cumsum(terms)
    return rhs, ordered, piprime_sqrt, Pm, terms, partial


CLASSDATA = {}
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
    Rr = [sum(M[i]) for i in range(4)]
    Cc = [sum(M[i][j] for i in range(4)) for j in range(4)]
    check(Rr == c5 and Cc == c1, "matrix marginals must equal class counts at K=%d" % K)
    print("   check: matrix row sums == class-5 counts, column sums == class-1 counts  [PASS]")
    sp = semiprime_patterns(x)
    print("   semiprime patterns (n<=x, gcd(n,6)=1):")
    for nm, c in sp.items():
        print("      %-18s %10d" % (nm, c))
    same = sp["(1,1) distinct"] + sp["(5,5) distinct"] + sp["p^2, p=1 mod 6"] + sp["p^2, p=5 mod 6"]
    mixed = sp["(1,5) mixed"]
    check(same == c1[1] and mixed == c5[1], "class-parity law for semiprimes (b even <-> same class)")
    print("      same-class (incl. squares) = %d = S_1 ; mixed = %d = S_5   [PROVED law, PASS]" % (same, mixed))
    rhs, ordered, pisq, Pm, terms, partial = pi_chi_identity(x)
    lhs = c1[1] - c5[1]
    check(lhs == rhs, "PROVED identity S_1-S_5 = (ordered + pi'(sqrt x))/2 fails: %d vs %d" % (lhs, rhs))
    print("   PROVED identity S_1-S_5 = (1/2)[sum_p chi(p) pi_chi(x/p) + pi'(sqrt x)]:")
    print("      S_1-S_5 = %d ; ordered double sum = %d ; pi'(sqrt x) = %d ; rhs = %d  [PASS]"
          % (lhs, ordered, pisq, rhs))
    print("      partial ordered sums through p = 5,7,11,13,17,19,23,29,31,37 and p<=100, p<=1000:")
    Pl = Pm.tolist()
    for pcut in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 100, 1000]:
        i = int(np.searchsorted(Pm, pcut, side="right"))
        print("         p<=%-5d : %8d" % (pcut, int(partial[i - 1]) if i > 0 else 0))
    # pi_chi at x itself and at x/5, x/7
    P = PRIMES[PRIMES >= 5]
    chi = np.where(P % 6 == 1, 1, -1).astype(np.int64)
    cum = np.concatenate([[0], np.cumsum(chi)])
    for yy in [x, x // 5, x // 7, x // 11, x // 13]:
        print("      pi_chi(%d) = %d" % (yy, int(cum[np.searchsorted(P, yy, side='right')])))

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
    MECH[K] = dict(act=act, pred=pred, semi=semi, prim=prim, res=act - pred)
    print("   %8d %8d %8d %8d %8d %10d %10.2f %10.2f %10.2f %10.2f %8.3f"
          % (K, R1, C1, R2, C2, act, float(pred), float(semi), float(prim), float(act - pred),
             act / float(pred) if pred != 0 else float('nan')))
print("   ratio form: N_PS/N_SP actual vs (pi_5 S_1)/(pi_1 S_5):")
for K in K_SCALES:
    c1, c5 = CLASSDATA[K]
    a = Fraction(MAT[K][0][1], MAT[K][1][0])
    p = Fraction(c5[0] * c1[1], c1[0] * c5[1])
    print("      K=%d : actual %.6f  predicted %.6f  (pi_5/pi_1 = %.6f, S_1/S_5 = %.6f)"
          % (K, float(a), float(p), c5[0] / c1[0], c1[1] / c5[1]))
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
    "(6k-1,6k+1) [5|1]": (6 * kk - 1, 6 * kk + 1),
    "(6k+1,6k+5) [1|5]": (6 * kk + 1, 6 * kk + 5),
    "(6k-1,6k+5) [5|5]": (6 * kk - 1, 6 * kk + 5),
    "(6k+1,6k+7) [1|1]": (6 * kk + 1, 6 * kk + 7),
}
CTRL = {}
for nm, (a_idx, b_idx) in PAIRS.items():
    La = CLS[a_idx]; Rb = CLS[b_idx]
    print()
    print(" pair %s" % nm)
    for K in K_SCALES:
        M = matrix4(La, Rb, K)
        R = [sum(M[i]) for i in range(4)]
        C = [sum(M[i][j] for i in range(4)) for j in range(4)]
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

print()
print(" (ii) local mod-35 strata of the original pair (6k-1,6k+1).")
print("     stratum = (h5-,h5+,h7-,h7+) divisibility pattern of (6k-1,6k+1) by 5 and 7.")
print("     Symmetric strata: 'none' (no 5,7 divisor either side) and 'cross' (one each side).")
print("     A local-sieve explanation would make the sign vanish in symmetric strata; H1 predicts it persists")
print("     with the same independence-predicted magnitude computed from the stratum's own marginals.")
d5m = ((6 * kk - 1) % 5 == 0); d5p = ((6 * kk + 1) % 5 == 0)
d7m = ((6 * kk - 1) % 7 == 0); d7p = ((6 * kk + 1) % 7 == 0)
STRATA = {
    "none  (h-=h+=0)": ~(d5m | d5p | d7m | d7p),
    "cross (5|left,7|right)": d5m & d7p,
    "cross (7|left,5|right)": d7m & d5p,
    "left-only (5|left, 7 none)": d5m & ~d7m & ~d7p,
    "right-only (5|right, 7 none)": d5p & ~d7m & ~d7p,
    "left-only (7|left, 5 none)": d7m & ~d5m & ~d5p,
    "right-only (7|right, 5 none)": d7p & ~d5m & ~d5p,
}
for K in K_SCALES:
    print()
    print("   K = %d" % K)
    tot_sym_act = 0
    tot_sym_pred = Fraction(0)
    for nm, mask in STRATA.items():
        m = mask[:K]
        n = int(m.sum())
        Ls = L0[:K][m]; Rs = R0[:K][m]
        M = matrix4(Ls, Rs, n)
        R = [sum(M[i]) for i in range(4)]
        C = [sum(M[i][j] for i in range(4)) for j in range(4)]
        act = M[0][1] - M[1][0]
        pred = Fraction(R[0] * C[1] - R[1] * C[0], n) if n else Fraction(0)
        if nm.startswith("none") or nm.startswith("cross"):
            tot_sym_act += act
            tot_sym_pred += pred
        print("     %-30s n=%8d  PP=%6d PS=%6d SP=%6d SS=%6d  PS-SP act %6d pred %8.2f res %8.2f"
              % (nm, n, M[0][0], M[0][1], M[1][0], M[1][1], act, float(pred), float(act - pred)))
    print("     symmetric strata combined: act %d  pred %.2f  res %.2f" %
          (tot_sym_act, float(tot_sym_pred), float(tot_sym_act - tot_sym_pred)))
    # CRT exactness check of the stratum sizes: on a full block of 35 consecutive k, sizes are fixed
    if K % 35 == 0:
        pass
blk = 35
sizes = {nm: int(mask[:blk].sum()) for nm, mask in STRATA.items()}
check(sizes["none  (h-=h+=0)"] == 3 * 5 and sizes["cross (5|left,7|right)"] == 1
      and sizes["cross (7|left,5|right)"] == 1, "CRT stratum sizes on one block of 35: %s" % sizes)
print("   CRT check: on one block of 35 centers the strata have sizes none=15, cross=1+1, left/right-only=%d,%d,%d,%d  [PASS]"
      % (sizes["left-only (5|left, 7 none)"], sizes["right-only (5|right, 7 none)"],
         sizes["left-only (7|left, 5 none)"], sizes["right-only (7|right, 5 none)"]))

print()
print(" (iii) prime-only Chebyshev census pi(x;6,5)-pi(x;6,1) and the running sign of the class-1 vs class-5")
print("      prime race up to 6*10^7+1, to place the sandwich sign next to the classical bias.")
P = PRIMES[PRIMES >= 5]
chi = np.where(P % 6 == 1, 1, -1).astype(np.int64)
cum = np.cumsum(chi)   # pi_chi after each prime
print("      primes >= 5 up to %d: %d ; pi_chi at end = %d ; #prefixes with pi_chi<0: %d, =0: %d, >0: %d"
      % (NMAX, len(P), int(cum[-1]), int((cum < 0).sum()), int((cum == 0).sum()), int((cum > 0).sum())))
pos = np.nonzero(cum > 0)[0]
print("      largest prime at which pi_chi > 0 (class 1 ahead): %s" % (int(P[pos[-1]]) if len(pos) else "never"))
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
print(" running sign K'<=%d: N_PS>N_SP for %d values, = for %d, < for %d; last K' with N_PS<=N_SP: %d"
      % (KMAX, rs["npos"], rs["nzero"], rs["nneg"], rs["lastnonpos"]))
rev = CTRL[("(6k+1,6k+5) [1|5]", "rs")]
print(" class-reversed pair (6k+1,6k+5): final N_PS-N_SP = %d ; K'<=%d with diff<0: %d, >0: %d"
      % (rev["final"], KMAX, rev["nneg"], rev["npos"]))
for nm in ["(6k-1,6k+5) [5|5]", "(6k+1,6k+7) [1|1]"]:
    r = CTRL[(nm, "rs")]
    print(" same-class pair %s: final N_PS-N_SP = %d ; K'<=%d with diff>0: %d, =0: %d, <0: %d ; min %s max %s"
          % (nm, r["final"], KMAX, r["npos"], r["nzero"], r["nneg"],
             min(r["track"].values()), max(r["track"].values())))
print(" total time %s" % elapsed())
print("ALL CHECKS PASSED")
