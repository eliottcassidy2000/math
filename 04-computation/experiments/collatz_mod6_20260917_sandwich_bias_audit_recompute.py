#!/usr/bin/env python3
"""
collatz_mod6_20260917_sandwich_bias_audit_recompute.py -- INDEPENDENT recomputation
for the adversarial audit (lane sandwich_bias, lens recompute).

Does NOT import the explorer's script.  Different algorithm: smallest-prime-factor
sieve stored compactly on n coprime to 6 (index n//3), then vectorised repeated
division to get (Omega, omega, b).  Everything load-bearing raises on failure.

Recomputed: 4x4 matrices at K=1e5,1e6,1e7; antisymmetric parts; running-sign
statistics of N_PS-N_SP (minimal witness, sign changes); class counts / chi-sums;
semiprime and 3-almost-prime patterns; the S_1-S_5 identity (brute force at small x
and at 6K+1); prime race; product-of-marginals prediction + three-way split; shifted
pairs; mod-35 strata; pooled-gap control with an HONEST block-based standard error
(the explorer's 50 gaps share the same left endpoints and are correlated);
sympy check of K=1453; the Omega=4 counterexample of Lemma 2.1.

Run: python3 04-computation/experiments/collatz_mod6_20260917_sandwich_bias_audit_recompute.py
"""
import sys
import time
from fractions import Fraction
from math import isqrt

import numpy as np

T0 = time.time()
GMAX = 300
NMAX = 6 * 10**7 + 1 + GMAX
KMAX = 10**7
SCALES = [10**5, 10**6, 10**7]


def check(c, msg):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def t():
    return "[t=%.1fs]" % (time.time() - T0)


def hr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ---------------------------------------------------------------- primes
hr("A. primes and compact smallest-prime-factor sieve on n coprime to 6")
isp = np.ones(NMAX + 1, dtype=np.bool_)
isp[:2] = False
for i in range(2, isqrt(NMAX) + 1):
    if isp[i]:
        isp[i * i::i] = False
PR = np.flatnonzero(isp).astype(np.int64)
del isp
print("pi(%d) = %d  %s" % (NMAX, len(PR), t()))
check(len(PR) == 3562133, "prime count differs from explorer's 3562133")

# compact index: n coprime to 6  ->  i = n//3 ; n(i) = 3i+1 (i even) or 3i+2 (i odd)
NC = NMAX // 3 + 1
spfc = np.zeros(NC, dtype=np.int32)
small = [int(p) for p in PR if p >= 5 and p * p <= NMAX]
for p in reversed(small):                      # descending: smallest prime wins
    for r in (1, 5):                           # n = p*m, m = 6j+r >= p
        j0 = (p - r + 5) // 6                  # smallest j with 6j+r >= p
        m0 = 6 * j0 + r
        start = (p * m0) // 3
        if p * m0 <= NMAX:
            spfc[start::2 * p] = p
print("spf sieve filled (%d primes)  %s" % (len(small), t()))

# vectorised factorisation of every n coprime to 6, n <= NMAX
idx = np.arange(NC, dtype=np.int64)
m = 3 * idx + 1 + (idx & 1)                    # n(i)
del idx
m = m.astype(np.int32)
OMc = np.zeros(NC, np.uint8)
OMDc = np.zeros(NC, np.uint8)
Bc = np.zeros(NC, np.uint8)
act = np.flatnonzero(m > 1).astype(np.int32)
prevp = np.zeros(len(act), np.int32)
it = 0
while act.size:
    mm = m[act]
    p = spfc[mm // 3]
    p = np.where(p == 0, mm, p)                # 0 in the spf table means prime
    OMc[act] += 1
    Bc[act] += (p % 6 == 5).astype(np.uint8)
    OMDc[act] += (p != prevp).astype(np.uint8)
    mm //= p
    m[act] = mm
    keep = mm > 1
    act = act[keep]
    prevp = p[keep]
    it += 1
del m, spfc, act, prevp
print("factorisation loop: %d iterations  %s" % (it, t()))


def om(n):   # Omega for n coprime to 6 (scalar)
    return int(OMc[n // 3])


# spot checks + independent trial division
def trial(n):
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


rng = np.random.RandomState(12345)
aud = [n for n in range(1, 40001) if n % 6 in (1, 5)]
aud += [int(x) for x in rng.randint(1, NMAX // 6, size=4000) * 6 + 1]
aud += [int(x) for x in rng.randint(1, NMAX // 6, size=4000) * 6 - 1]
for n in aud:
    check((int(OMc[n // 3]), int(OMDc[n // 3]), int(Bc[n // 3])) == trial(n), "trial mismatch n=%d" % n)
print("trial-division audit on %d integers coprime to 6: agree  %s" % (len(aud), t()))
CLSc = np.minimum(OMc, 4).astype(np.uint8)

# class parity law on the whole universe
check(bool(np.all(Bc[1:2 * KMAX:2] % 2 == 1)), "b(6k-1) odd")
check(bool(np.all(Bc[2:2 * KMAX + 1:2] % 2 == 0)), "b(6k+1) even")
print("class-parity law b(6k-1) odd / b(6k+1) even re-verified for k<=%d" % KMAX)


# ---------------------------------------------------------------- helpers
def endpoint(c, K):
    """CLS of 6k+c for k=1..K, c = +-1 mod 6 (c may be >6)."""
    a = 2 + c // 3
    return CLSc[a:a + 2 * K:2]


def mat(L, R):
    code = (L.astype(np.int16) - 1) * 4 + (R.astype(np.int16) - 1)
    cnt = np.bincount(code, minlength=16)
    return [[int(cnt[4 * i + j]) for j in range(4)] for i in range(4)]


def anti(M, K, i, j):
    R = [sum(M[a]) for a in range(4)]
    C = [sum(M[a][b] for a in range(4)) for b in range(4)]
    act = M[i][j] - M[j][i]
    pred = Fraction(R[i] * C[j] - R[j] * C[i], K)
    return act, pred


INHERITED_K1E6 = [[37915, 78689, 59706, 30192], [78277, 157420, 112416, 52182],
                  [59992, 112305, 72363, 28755], [30161, 52125, 28736, 8766]]
EXPL = {  # explorer's claimed matrices
    10**5: [[5330, 10091, 6537, 2614], [10050, 17794, 10719, 3851], [6561, 10751, 5600, 1664], [2583, 3844, 1666, 345]],
    10**6: INHERITED_K1E6,
    10**7: [[280557, 632905, 539140, 328636], [632766, 1380659, 1118019, 632028],
            [539061, 1118744, 837851, 416682], [328491, 631688, 416902, 165871]],
}
EXPL_ANTI = {10**5: (41, -24, -32), 10**6: (412, -286, 111), 10**7: (139, 79, -725)}
EXPL_PRED = {10**5: "36.58", 10**6: "113.23", 10**7: "229.95"}

hr("B. SB1: ordered matrices (6k-1,6k+1) and antisymmetric parts")
MAT = {}
for K in SCALES:
    L = endpoint(-1, K); R = endpoint(1, K)
    M = mat(L, R)
    MAT[K] = M
    check(M == EXPL[K], "matrix at K=%d differs from explorer" % K)
    a = [anti(M, K, 0, 1), anti(M, K, 0, 2), anti(M, K, 1, 2)]
    check(tuple(x[0] for x in a) == EXPL_ANTI[K], "antisym at K=%d" % K)
    check("%.2f" % float(a[0][1]) == EXPL_PRED[K], "pred at K=%d" % K)
    print("K=%d matrix matches explorer; PS-SP=%d (pred %.2f) PC-CP=%d (pred %.2f) SC-CS=%d (pred %.2f)"
          % (K, a[0][0], float(a[0][1]), a[1][0], float(a[1][1]), a[2][0], float(a[2][1])))

hr("C. SB2: running sign of N_PS-N_SP, K'<=1e7")
L = endpoint(-1, KMAX); R = endpoint(1, KMAX)
step = ((L == 1) & (R == 2)).astype(np.int8) - ((L == 2) & (R == 1)).astype(np.int8)
d = np.cumsum(step, dtype=np.int32)
del step
npos, nzero, nneg = int((d > 0).sum()), int((d == 0).sum()), int((d < 0).sum())
s = np.sign(d).astype(np.int8)
nzpos = np.flatnonzero(s)
nz = s[nzpos]
chg = np.flatnonzero(nz[1:] != nz[:-1])          # index into nz of the entry BEFORE the change
chgK = nzpos[chg + 1] + 1                        # K' at which the new sign first holds
firstneg = int(np.argmax(d < 0)) + 1
print("diff>0: %d ; =0: %d ; <0: %d ; sign changes (nonzero-subsequence transitions): %d" % (npos, nzero, nneg, len(chg)))
print("first change at K'=%d, last at K'=%d ; first K' with diff<0: %d (diff=%d)" % (chgK[0], chgK[-1], firstneg, int(d[firstneg - 1])))
lastnonpos = int(np.flatnonzero(d <= 0)[-1]) + 1
print("last K' with diff<=0: %d ; min %d ; max %d ; final %d ; positive fraction %.4f"
      % (lastnonpos, int(d.min()), int(d.max()), int(d[-1]), npos / KMAX))
check((npos, nzero, nneg, len(chg), int(chgK[0]), int(chgK[-1]), firstneg, lastnonpos, int(d.min()), int(d.max()))
      == (7716244, 14947, 2268809, 992, 1453, 8714332, 1453, 8714401, -436, 650), "running-sign stats differ from explorer")
print("running-sign statistics agree with explorer  [PASS]")
# arcsine-law tail for the positive fraction
from math import asin, pi, sqrt
print("arcsine P(frac>=%.4f) = %.4f (explorer says about 0.32)" % (npos / KMAX, 1 - 2 / pi * asin(sqrt(npos / KMAX))))
# independent sympy witness check at K=1453
try:
    import sympy
    nps = nsp = 0
    for k in range(1, 1454):
        a = sum(sympy.factorint(6 * k - 1).values()); b = sum(sympy.factorint(6 * k + 1).values())
        nps += (a == 1 and b == 2); nsp += (a == 2 and b == 1)
    print("sympy: K=1453 N_PS=%d N_SP=%d" % (nps, nsp))
    check((nps, nsp) == (245, 246), "sympy witness")
    # also confirm no earlier K has N_PS<N_SP via sympy running counts
    nps = nsp = 0; first = None
    for k in range(1, 1454):
        a = sum(sympy.factorint(6 * k - 1).values()); b = sum(sympy.factorint(6 * k + 1).values())
        nps += (a == 1 and b == 2); nsp += (a == 2 and b == 1)
        if nps < nsp and first is None:
            first = k
    check(first == 1453, "sympy minimal witness")
    print("sympy: minimal K with N_PS<N_SP is %d  [PASS]" % first)
except ImportError:
    print("sympy unavailable; witness checked by sieve only")
del d, s, nzpos, nz

hr("D. SB3/SB5: class counts, chi-sums, semiprime and 3-almost-prime patterns")
CLASS = {}
for K in SCALES:
    c5 = [int((endpoint(-1, K) == j).sum()) for j in range(1, 5)]   # n=5 mod 6, n<=6K-1
    c1 = [int((endpoint(1, K) == j).sum()) for j in range(1, 5)]    # n=1 mod 6, 7<=n<=6K+1
    CLASS[K] = (c1, c5)
    M = MAT[K]
    check([sum(M[i]) for i in range(4)] == c5 and [sum(M[i][j] for i in range(4)) for j in range(4)] == c1, "marginals")
    chi = [c1[j] - c5[j] for j in range(4)]
    print("K=%d  class1 %s  class5 %s  chi-sums %s" % (K, c1, c5, chi))
EXPL_CHI = {10**5: [-48, 66, -54, 36], 10**6: [-157, 244, -194, 107], 10**7: [-363, 524, -426, 265]}
for K in SCALES:
    c1, c5 = CLASS[K]
    check([c1[j] - c5[j] for j in range(4)] == EXPL_CHI[K], "chi-sums K=%d" % K)
print("chi-sums agree with explorer  [PASS]")
# patterns at the largest scale
K = 10**7
sl5 = slice(1, 2 * K, 2); sl1 = slice(2, 2 * K + 1, 2)
w5, b5, o5 = OMDc[sl5], Bc[sl5], OMc[sl5]
w1, b1, o1 = OMDc[sl1], Bc[sl1], OMc[sl1]
semi = {
    "(1,1)": int(((o1 == 2) & (w1 == 2) & (b1 == 0)).sum()),
    "(5,5)": int(((o1 == 2) & (w1 == 2) & (b1 == 2)).sum()),
    "mixed": int(((o5 == 2) & (w5 == 2) & (b5 == 1)).sum()),
    "p1^2": int(((o1 == 2) & (w1 == 1)).sum()),
    "p5^2": int(((o1 == 2) & (w1 == 1) & (b1 == 2)).sum()),
}
# every semiprime on the right must be same-class or square; on the left mixed
check(int(((o1 == 2) & (b1 == 1)).sum()) == 0 and int(((o5 == 2) & (b5 != 1)).sum()) == 0, "semiprime side law")
# p^2 with p = 1 mod 6 has b=0; recount to split squares
sq1 = int(((o1 == 2) & (w1 == 1) & (b1 == 0)).sum())
semi["p1^2"] = sq1
print("x=%d semiprimes: %s" % (6 * K + 1, semi))
check(semi == {"(1,1)": 1732047, "(5,5)": 2030969, "mixed": 3763472, "p1^2": 484, "p5^2": 496}, "semiprime patterns")
cube = {}
for (wi, bi) in [(1, 0), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]:
    side = o1 if bi % 2 == 0 else o5
    w = w1 if bi % 2 == 0 else w5
    b = b1 if bi % 2 == 0 else b5
    cube[(wi, bi)] = int(((side == 3) & (w == wi) & (b == bi)).sum())
print("x=%d 3-almost-prime (omega,b) census: %s" % (6 * K + 1, cube))
check(cube == {(1, 0): 36, (1, 3): 39, (2, 0): 84521, (2, 1): 85178, (2, 2): 138521, (2, 3): 139167,
               (3, 0): 522639, (3, 1): 1909337, (3, 2): 2166195, (3, 3): 778617}, "3-almost patterns")
check(sum(v for (w_, b_), v in cube.items() if b_ % 2 == 0) == CLASS[K][0][2], "right total = C_1")
check(sum(v for (w_, b_), v in cube.items() if b_ % 2 == 1) == CLASS[K][1][2], "left total = C_5")
print("semiprime / 3-almost pattern censuses agree with explorer  [PASS]")
# Lemma 2.1 wording: multiset of classes is determined by (Omega,b) ALONE for every Omega
# (b fives and Omega-b ones).  The Omega=4 'counterexample' has the SAME class multiset.
n_a, n_b = 7 * 7 * 5 * 11, 7 * 13 * 5 * 5
ta, tb = trial(n_a), trial(n_b)
print("Omega=4 example: %d=7^2*5*11 -> (Omega,omega,b)=%s ; %d=7*13*5^2 -> %s ; class multisets both {1,1,5,5}"
      % (n_a, ta, n_b, tb))
check(ta == tb == (4, 3, 2), "counterexample triple")

hr("E. SB4: identity S_1(x)-S_5(x) = (O(x) + pi'(sqrt x))/2")


def brute_identity(x):
    # brute force S_1-S_5 by trial factorisation; O and pi' by explicit double loop over primes
    s1 = s5 = 0
    for n in range(5, x + 1):
        if n % 6 in (1, 5) and trial(n)[0] == 2:
            if n % 6 == 1:
                s1 += 1
            else:
                s5 += 1
    ps = [int(p) for p in PR if 5 <= p <= x]
    O = 0
    for p in ps:
        for q in ps:
            if p * q > x:
                break
            O += (1 if p % 6 == 1 else -1) * (1 if q % 6 == 1 else -1)
    pisq = sum(1 for p in ps if p * p <= x)
    return s1 - s5, O, pisq


for x in [1, 24, 25, 26, 35, 49, 100, 1000, 10007]:
    lhs, O, pisq = brute_identity(x)
    check((O + pisq) % 2 == 0 and lhs == (O + pisq) // 2, "identity fails at x=%d" % x)
    print("x=%d: S_1-S_5=%d ; O=%d ; pi'(sqrt x)=%d  [PASS]" % (x, lhs, O, pisq))
print("identity holds for all x>=1 (both sides vanish for x<25); the 'x>=25' hypothesis is unnecessary")

P5 = PR[PR >= 5]
chiP = np.where(P5 % 6 == 1, 1, -1).astype(np.int32)
cum = np.concatenate([[0], np.cumsum(chiP)])


def pichi(y):
    return int(cum[np.searchsorted(P5, y, side="right")])


EXPL_ID = {10**5: (66, -3, 135), 10**6: (244, 127, 361), 10**7: (524, 68, 980)}
for K in SCALES:
    x = 6 * K + 1
    Pm = P5[P5 <= x // 5]
    y = x // Pm
    O = int((np.where(Pm % 6 == 1, 1, -1) * cum[np.searchsorted(P5, y, side="right")]).sum())
    pisq = int((P5 <= isqrt(x)).sum())
    lhs = CLASS[K][0][1] - CLASS[K][1][1]
    check(lhs == (O + pisq) // 2 and (O + pisq) % 2 == 0, "identity at x=%d" % x)
    check((lhs, O, pisq) == EXPL_ID[K], "identity numbers K=%d" % K)
    print("x=%d: S_1-S_5=%d = (%d + %d)/2  [PASS]" % (x, lhs, O, pisq))

hr("F. SB5: prime race pi_chi over primes >= 5 up to %d" % NMAX)
c = cum[1:]
print("pi_chi at end %d ; prefixes <0: %d, =0: %d, >0: %d ; zero prefixes at primes %s"
      % (int(c[-1]), int((c < 0).sum()), int((c == 0).sum()), int((c > 0).sum()), P5[c == 0].tolist()))
check(int((c > 0).sum()) == 0 and int((c == 0).sum()) == 9, "prime race")

hr("G. SB6: three-way split of the product-of-marginals prediction")
for K in SCALES:
    c1, c5 = CLASS[K]
    R1, C1, R2, C2 = c5[0], c1[0], c5[1], c1[1]
    pred = Fraction(R1 * C2 - R2 * C1, K)
    pibar = Fraction(R1 + C1, 2); dpi = Fraction(R1 - C1, 2); Sbar = Fraction(R2 + C2, 2); dS = Fraction(C2 - R2, 2)
    check(pred == 2 * (pibar * dS + Sbar * dpi) / K, "algebraic split")
    lhs, O, pisq = EXPL_ID[K]
    check(4 * dS == O + pisq, "dS = (O+pi')/4")
    prim = 2 * Sbar * dpi / K; sq = pibar * pisq / (2 * K); cr = pibar * O / (2 * K)
    check(prim + sq + cr == pred, "three-way sum")
    act = MAT[K][0][1] - MAT[K][1][0]
    print("K=%d: actual %d pred %.2f = prime %.2f (%.1f%%) + square %.2f (%.1f%%) + cross %.2f (%.1f%%) ; act/pred %.3f ; ratio N_PS/N_SP %.6f vs pred %.6f"
          % (K, act, float(pred), float(prim), 100 * float(prim / pred), float(sq), 100 * float(sq / pred),
             float(cr), 100 * float(cr / pred), act / float(pred),
             MAT[K][0][1] / MAT[K][1][0], (R1 * C2) / (C1 * R2)))
    # PC and SC predictions
    for (i, j, nm) in [(0, 2, "PC-CP"), (1, 2, "SC-CS")]:
        a, p = anti(MAT[K], K, i, j)
        print("      %s actual %d pred %.2f" % (nm, a, float(p)))

hr("H. SB8: shifted single pairs")
EXPL_SHIFT = {(1, 5): (-5, -43, 461), (-1, 5): (39, -170, 129), (1, 7): (109, -187, 108)}
for (ca, cb), exp in EXPL_SHIFT.items():
    vals = []
    for K in SCALES:
        M = mat(endpoint(ca, K), endpoint(cb, K))
        a, p = anti(M, K, 0, 1)
        vals.append(a)
        print("pair (6k%+d,6k%+d) K=%d: PS-SP=%d pred %.2f" % (ca, cb, K, a, float(p)))
    check(tuple(vals) == exp, "shifted pair (%d,%d)" % (ca, cb))
print("shifted-pair values agree with explorer  [PASS]")

hr("I. SB9: mod-35 strata of (6k-1,6k+1); CRT sizes; involution check")
kk = np.arange(1, KMAX + 1, dtype=np.int64)
d5m = (6 * kk - 1) % 5 == 0; d5p = (6 * kk + 1) % 5 == 0; d7m = (6 * kk - 1) % 7 == 0; d7p = (6 * kk + 1) % 7 == 0
none = ~(d5m | d5p | d7m | d7p)
crossA = d5m & d7p; crossB = d7m & d5p
blk = 35
print("block of 35: none=%d crossA=%d crossB=%d 5|left-only=%d 5|right-only=%d 7|left-only=%d 7|right-only=%d"
      % (int(none[:blk].sum()), int(crossA[:blk].sum()), int(crossB[:blk].sum()),
         int((d5m & ~d7m & ~d7p)[:blk].sum()), int((d5p & ~d7m & ~d7p)[:blk].sum()),
         int((d7m & ~d5m & ~d5p)[:blk].sum()), int((d7p & ~d5m & ~d5p)[:blk].sum())))
check(int(none[:blk].sum()) == 15 and int(crossA[:blk].sum()) == 1 and int(crossB[:blk].sum()) == 1, "CRT sizes")
EXPL_NONE = {10**5: (21, "5.60"), 10**6: (182, "139.37"), 10**7: (-334, "255.22")}
for K in SCALES:
    mk = none[:K]
    M = mat(endpoint(-1, K)[mk], endpoint(1, K)[mk])
    a, p = anti(M, int(mk.sum()), 0, 1)
    print("K=%d 'none' stratum n=%d PS-SP=%d pred %.2f" % (K, int(mk.sum()), a, float(p)))
    check(a == EXPL_NONE[K][0] and "%.2f" % float(p) == EXPL_NONE[K][1], "none stratum K=%d" % K)
    # cross strata: primes present?
    for nm, mk2 in [("crossA", crossA[:K]), ("crossB", crossB[:K])]:
        Lc = endpoint(-1, K)[mk2]; Rc = endpoint(1, K)[mk2]
        print("      %s: #P on left %d, #P on right %d (only k=1: (5,7))" % (nm, int((Lc == 1).sum()), int((Rc == 1).sum())))
del kk, d5m, d5p, d7m, d7p, none, crossA, crossB
# involution: k* = -k - (2c+g)*inv6 mod M sends (6k+c, 6k+c+g) to (-(6k+c+g), -(6k+c)) mod M
for Mmod in [5, 7, 35, 385, 1001]:
    inv6 = pow(6, -1, Mmod)
    for c, g in [(-1, 2), (1, 4), (-1, 6), (1, 6), (5, 8), (1, 300)]:
        for k in range(Mmod):
            ks = (-k - (2 * c + g) * inv6) % Mmod
            check((6 * ks + c - (-(6 * k + c + g))) % Mmod == 0 and (6 * ks + c + g - (-(6 * k + c))) % Mmod == 0,
                  "involution M=%d c=%d g=%d k=%d" % (Mmod, c, g, k))
            check((-ks - (2 * c + g) * inv6) % Mmod == k, "involution is an involution")
print("involution k -> -k-(2c+g)/6 mod M verified exhaustively for M in {5,7,35,385,1001}, six (c,g)  [PASS]")

hr("J. SB7: pooled-gap control with block-based standard errors")
NB = 100
EXPL_POOL = {  # explorer: (orientation) -> (mean PS-SP, mean pred) at K=1e6, 1e7
    10**6: {"[5|1]": (100.7, 117.5), "[1|5]": (-113.0, -109.1), "[5|5]": (9.7, 4.1), "[1|1]": (12.5, 4.4)},
    10**7: {"[5|1]": (195.1, 234.2), "[1|5]": (-66.4, -225.6), "[5|5]": (96.7, 4.3), "[1|1]": (-32.7, 4.4)},
}
for K in [10**6, 10**7]:
    print("\nK=%d  (%d blocks of %d centres)" % (K, NB, K // NB))
    res = {}   # (lab, ori) -> list of (g, D_total, pred, D_blocks)
    for c in (-1, 1):
        Lc = endpoint(c, K)
        Lm = {1: Lc == 1, 2: Lc == 2, 3: Lc == 3}
        for g in range(2, GMAX + 1, 2):
            e = (c + g) % 6
            if e == 3:
                continue
            ori = "[%d|%d]" % (c % 6, e)
            Rc = endpoint(c + g, K)
            Rm = {1: Rc == 1, 2: Rc == 2, 3: Rc == 3}
            M = mat(Lc, Rc)
            for (i, j, lab) in [(1, 2, "PS-SP"), (1, 3, "PC-CP"), (2, 3, "SC-CS")]:
                blocks = (Lm[i] & Rm[j]).reshape(NB, -1).sum(1) - (Lm[j] & Rm[i]).reshape(NB, -1).sum(1)
                a, p = anti(M, K, i - 1, j - 1)
                check(int(blocks.sum()) == a, "block sum")
                res.setdefault((lab, ori), []).append((g, a, p, blocks.astype(np.int64)))
    for lab in ["PS-SP", "PC-CP", "SC-CS"]:
        print("  %s:" % lab)
        for ori in ["[5|1]", "[1|5]", "[5|5]", "[1|1]"]:
            lst = res[(lab, ori)]
            n = len(lst)
            acts = np.array([a for (_, a, _, _) in lst], dtype=np.float64)
            preds = float(sum(p for (_, _, p, _) in lst) / n)
            mean = acts.mean(); sd = acts.std(ddof=1)
            Bm = np.stack([b for (_, _, _, b) in lst])          # (n gaps, NB blocks)
            mb = Bm.mean(0)                                     # per-block pooled mean over gaps
            se_block = sqrt(NB * mb.var(ddof=1))                # SE of the pooled mean if blocks ~independent
            # per-gap noise from blocks (should be ~ sd across gaps if gaps were exchangeable)
            pergap_sd_blocks = float(np.sqrt(NB * Bm.var(1, ddof=1)).mean())
            # average pairwise correlation between gaps, estimated from blocks
            Cm = np.corrcoef(Bm)
            rho = float((Cm.sum() - n) / (n * (n - 1)))
            nposb = int((mb > 0).sum())
            # paired sign test over blocks: two-sided binomial p-value for #positive blocks
            from math import comb
            kpos = min(nposb, NB - nposb)
            pval = 2 * sum(comb(NB, i) for i in range(kpos + 1)) / 2**NB
            print("    %s n=%d mean=%8.1f (naive sd %6.1f, naive SE %5.1f) pred %7.1f | block SE %6.1f  z_block=%5.2f | blocks>0: %d/%d (sign p=%.3f) | mean inter-gap corr %.3f | per-gap sd from blocks %.1f"
                  % (ori, n, mean, sd, sd / sqrt(n), preds, se_block, mean / se_block, nposb, NB, pval, rho, pergap_sd_blocks))
            if lab == "PS-SP":
                em, ep = EXPL_POOL[K][ori]
                check(abs(mean - em) < 0.06 and abs(preds - ep) < 0.06, "pooled mean %s K=%d" % (ori, K))
    print("  pooled PS-SP means and predictions agree with explorer  [PASS]  %s" % t())

hr("K. Sanity: per-gap sd vs sqrt(K) scaling, and orientation falsifier at 1e7")
print("Explorer's falsifier text: 'if same-class pooled means are as large as opposite-class ones, H1 is refuted'.")
print("At K=1e7 PS-SP: |[5|5]| = 96.7 > |[1|5]| = 66.4 ; at K=1e6 both same-class means are < both opposite-class means.")
print("total time %s" % t())
print("ALL AUDIT CHECKS PASSED")
