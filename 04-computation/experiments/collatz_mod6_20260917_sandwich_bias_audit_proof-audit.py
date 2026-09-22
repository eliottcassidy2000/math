#!/usr/bin/env python3
"""
collatz_mod6_20260917_sandwich_bias_audit_proof-audit.py

INDEPENDENT recomputation for the proof-audit of lane sandwich_bias.
Does not import the explorer's script.  Uses a smallest-prime-factor (spf)
sieve (uint16: every composite <= 6e7 has spf <= 7745; primes carry a 0
sentinel) and per-number factor chains (the explorer used prime-power
slicing), restricted to n = +-1 mod 6, then recomputes every load-bearing
number and probes the hostile cases found while reading the proofs:

  A. the mod-35 strata listed in the explorer's section 4(ii) do NOT
     partition the centers (35|6k-1 and 35|6k+1 strata are missing);
  B. Lemma 2.1 as worded (multiset of classes) is trivially true for ALL n
     -- the content is the class-labelled exponent shape;
  C. the three-way split of Theorem 2.2 is convention dependent: with the
     squares counted in full, distinct-pair semiprimes favour class 5.

Memory: spf uint16 120 MB, CLS uint8 60 MB, class-indexed omega/b 40 MB.
"""
import sys
import time
import math
import resource
from fractions import Fraction

import numpy as np

T0 = time.time()
GMAX = 300
NMAX = 6 * 10**7 + 1 + GMAX
KMAX = 10**7
K_SCALES = [10**5, 10**6, 10**7]


def el():
    return "[t=%.1fs]" % (time.time() - T0)


def check(c, msg):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def hr(t):
    print(); print("=" * 78); print(t)
    print("(peak RSS so far: %.2f GB)" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30))
    print("=" * 78)


# ---------------------------------------------------------------- spf sieve
hr("0. INDEPENDENT spf SIEVE on [1, %d]" % NMAX)
r = math.isqrt(NMAX)
check(r < 65535, "uint16 spf needs sqrt(NMAX) < 65535")
spf = np.zeros(NMAX + 1, np.uint16)
for p in range(2, r + 1):
    if spf[p] == 0:
        v = spf[p::p]
        v[v == 0] = p          # stamps spf[p]=p for the prime p itself too
isp = spf == 0
isp[:2] = False
for p in range(2, r + 1):      # small primes were stamped with themselves
    if spf[p] == p:
        isp[p] = True
PR = np.nonzero(isp)[0].astype(np.int64)
del isp
print("primes <= %d : %d  %s" % (NMAX, len(PR), el()))
check(len(PR) == 3562133, "prime count differs from explorer (3562133)")

# CLS[n] = min(Omega(n),4) on the full range (only for n = +-1 mod 6);
# omega and b are stored per residue class: index (n-5)//6 for n=5 mod 6, (n-7)//6 for n=1 mod 6
CLS = np.zeros(NMAX + 1, np.uint8)
NCLS = (NMAX - 5) // 6 + 1
OMD5 = np.zeros(NCLS, np.uint8); BB5 = np.zeros(NCLS, np.uint8)
OMD1 = np.zeros(NCLS, np.uint8); BB1 = np.zeros(NCLS, np.uint8)
CHUNK = 2 * 10**6
itmax = 0
for start, OMDc, BBc in ((5, OMD5, BB5), (7, OMD1, BB1)):
    nsel = (NMAX - start) // 6 + 1
    for c0 in range(0, nsel, CHUNK):
        sel = np.arange(start + 6 * c0, min(start + 6 * (c0 + CHUNK), NMAX + 1), 6, dtype=np.int32)
        work = sel.copy()
        OM_s = np.zeros(len(sel), np.uint8)
        OMD_s = np.zeros(len(sel), np.uint8)
        BB_s = np.zeros(len(sel), np.uint8)
        prev = np.zeros(len(sel), np.int32)
        active = np.nonzero(work > 1)[0].astype(np.int32)
        it = 0
        while active.size:
            w = work[active]
            p = spf[w].astype(np.int32)
            p = np.where(p == 0, w, p)          # 0 sentinel: w itself is prime
            OM_s[active] += 1
            BB_s[active] += (p % 6 == 5).astype(np.uint8)
            OMD_s[active] += (p != prev[active]).astype(np.uint8)
            prev[active] = p
            work[active] = w // p
            active = active[work[active] > 1]
            it += 1
        itmax = max(itmax, it)
        CLS[sel] = np.minimum(OM_s, 4)
        OMDc[c0:c0 + len(sel)] = OMD_s; BBc[c0:c0 + len(sel)] = BB_s
        del work, prev, active, OM_s, OMD_s, BB_s, sel
print("factor chains done (max %d iterations per chunk)  %s" % (itmax, el()))
del spf


def wb(n):
    """(omega, b) for n = +-1 mod 6 (n an int or an array of a SINGLE residue class), uint8 views."""
    arr = np.asarray(n, dtype=np.int32)
    cls = int(arr.flat[0]) % 6
    if arr.size:
        check(bool(np.all(arr % 6 == cls)), "wb: mixed residue classes")
    if cls == 5:
        idx = (arr - 5) // 6
        return OMD5[idx], BB5[idx]
    idx = (arr - 7) // 6
    return OMD1[idx], BB1[idx]


# independent spot audit against sympy on random n = +-1 mod 6
try:
    from sympy import factorint
    rng = np.random.RandomState(4242)
    cnt = 0
    for n in rng.randint(5, NMAX + 1, size=1500):
        n = int(n)
        if n % 6 not in (1, 5):
            continue
        f = factorint(n)
        O = sum(f.values()); w = len(f); b = sum(e for q, e in f.items() if q % 6 == 5)
        ww, bb = wb(n)
        check((int(CLS[n]), int(ww), int(bb)) == (min(O, 4), w, b), "sympy mismatch at n=%d" % n)
        cnt += 1
    print("sympy factorint audit: %d random n = +-1 mod 6 agree  %s" % (cnt, el()))
    HAVE_SYMPY = True
except ImportError:
    HAVE_SYMPY = False
    print("sympy not available; skipped random audit")

kk = np.arange(1, KMAX + 1, dtype=np.int32)
L0 = CLS[6 * kk - 1]
R0 = CLS[6 * kk + 1]
check(bool(np.all(BB5[:KMAX] % 2 == 1)) and bool(np.all(BB1[:KMAX] % 2 == 0)), "parity law SW4")
print("class-parity law SW4 re-verified for k <= %d" % KMAX)


def mat4(L, R, K):
    Lk = L[:K]; Rk = R[:K]
    return [[int(np.count_nonzero((Lk == i) & (Rk == j))) for j in range(1, 5)] for i in range(1, 5)]


NOTE_1E5 = [[5330, 10091, 6537, 2614], [10050, 17794, 10719, 3851],
            [6561, 10751, 5600, 1664], [2583, 3844, 1666, 345]]
NOTE_1E6 = [[37915, 78689, 59706, 30192], [78277, 157420, 112416, 52182],
            [59992, 112305, 72363, 28755], [30161, 52125, 28736, 8766]]
NOTE_1E7 = [[280557, 632905, 539140, 328636], [632766, 1380659, 1118019, 632028],
            [539061, 1118744, 837851, 416682], [328491, 631688, 416902, 165871]]

hr("1. MATRICES (SB1) AND RUNNING SIGNS (SB2)")
MAT = {}
for K, ref in zip(K_SCALES, [NOTE_1E5, NOTE_1E6, NOTE_1E7]):
    M = mat4(L0, R0, K)
    MAT[K] = M
    check(M == ref, "matrix K=%d differs from note" % K)
    print("K=%d matrix identical to the note  [PASS]  PS-SP=%d PC-CP=%d SC-CS=%d"
          % (K, M[0][1] - M[1][0], M[0][2] - M[2][0], M[1][2] - M[2][1]))
check([MAT[K][0][1] - MAT[K][1][0] for K in K_SCALES] == [41, 412, 139], "PS-SP triple")
check([MAT[K][0][2] - MAT[K][2][0] for K in K_SCALES] == [-24, -286, 79], "PC-CP triple")
check([MAT[K][1][2] - MAT[K][2][1] for K in K_SCALES] == [-32, 111, -725], "SC-CS triple")

K_TRACK = [10**4, 2 * 10**4, 5 * 10**4, 10**5, 2 * 10**5, 5 * 10**5, 10**6, 2 * 10**6, 5 * 10**6, 10**7]


def running(L, R, a, b, K):
    d = np.cumsum(((L[:K] == a) & (R[:K] == b)).astype(np.int8)
                  - ((L[:K] == b) & (R[:K] == a)).astype(np.int8), dtype=np.int32)
    s = np.sign(d).astype(np.int8)
    nz = s[s != 0]
    changes = int(np.count_nonzero(nz[1:] != nz[:-1]))
    del s, nz
    firstneg = int(np.argmax(d < 0)) + 1 if bool((d < 0).any()) else 0
    lastnonpos = K - int(np.argmax((d <= 0)[::-1])) if bool((d <= 0).any()) else 0
    out = dict(npos=int((d > 0).sum()), nzero=int((d == 0).sum()), nneg=int((d < 0).sum()),
               changes=changes, firstneg=firstneg, lastnonpos=lastnonpos,
               minv=int(d.min()), maxv=int(d.max()), final=int(d[-1]),
               track={Kt: int(d[Kt - 1]) for Kt in K_TRACK if Kt <= K})
    del d
    return out


rs = running(L0, R0, 1, 2, KMAX)
print("PS-SP running: >0 %d =0 %d <0 %d changes %d firstneg %d lastnonpos %d min %d max %d"
      % (rs["npos"], rs["nzero"], rs["nneg"], rs["changes"], rs["firstneg"], rs["lastnonpos"], rs["minv"], rs["maxv"]))
check((rs["npos"], rs["nzero"], rs["nneg"], rs["changes"], rs["firstneg"], rs["lastnonpos"], rs["minv"], rs["maxv"])
      == (7716244, 14947, 2268809, 992, 1453, 8714401, -436, 650), "PS-SP running summary vs note")
print("PS-SP at tracked K: " + ", ".join("%d:%d" % (Kt, v) for Kt, v in rs["track"].items()))
neg_tracked = [Kt for Kt, v in rs["track"].items() if v < 0]
print("tracked K with PS-SP<0: %s  (the note's 'sign predicted at all three scales' uses K=10^5,10^6,10^7 only)" % neg_tracked)
rs2 = running(L0, R0, 1, 3, KMAX)
rs3 = running(L0, R0, 2, 3, KMAX)
print("PC-CP running: changes %d firstneg %d ; SC-CS: <0 %d changes %d min %d max %d"
      % (rs2["changes"], rs2["firstneg"], rs3["nneg"], rs3["changes"], rs3["minv"], rs3["maxv"]))
check(rs2["changes"] == 639 and rs2["firstneg"] == 21, "PC-CP running vs note")
check(rs3["nneg"] == 9114559 and rs3["changes"] == 470 and rs3["minv"] == -1323 and rs3["maxv"] == 319, "SC-CS running vs note")
del rs2, rs3

# K=1453 witness, brute force with sympy (independent of every sieve)
if HAVE_SYMPY:
    nps = nsp = 0
    first_bad = None
    for k in range(1, 1454):
        a = sum(factorint(6 * k - 1).values()); b = sum(factorint(6 * k + 1).values())
        if a == 1 and b == 2:
            nps += 1
        elif a == 2 and b == 1:
            nsp += 1
        if nps < nsp and first_bad is None:
            first_bad = (k, nps, nsp)
    print("sympy brute force: first K with N_PS<N_SP is %s ; at K=1453: N_PS=%d N_SP=%d" % (first_bad, nps, nsp))
    check(first_bad == (1453, 245, 246), "K=1453 witness")

frac = rs["npos"] / KMAX
p_arcsine = 1 - (2 / math.pi) * math.asin(math.sqrt(frac))
print("fraction of K' with PS-SP>0 = %.4f ; arcsine-law P(fraction >= this) = %.3f (note: about 0.32)" % (frac, p_arcsine))
del rs
print(el())

# ---------------------------------------------------------------- class counts
hr("2. CLASS COUNTS, SEMIPRIME PATTERNS, THEOREM 2.2 (SB3-SB6)")
chi_p = np.where(PR % 6 == 1, 1, -1).astype(np.int8)
chi_p[PR < 5] = 0
cum = np.concatenate([[0], np.cumsum(chi_p, dtype=np.int32)])   # cum[i] = sum chi over first i primes
cnt_ge5 = np.concatenate([[0], np.cumsum((PR >= 5).astype(np.int32))])


def pi_chi(y):
    return int(cum[int(np.searchsorted(PR, y, side="right"))])


def pi_prime(y):
    return int(cnt_ge5[int(np.searchsorted(PR, y, side="right"))])


def ordered_O(x):
    """O(x) = sum_{5<=p<=x/5} chi(p) pi_chi(x/p), computed my way."""
    Pm = PR[(PR >= 5) & (PR <= x // 5)]
    idx = np.searchsorted(PR, x // Pm, side="right")
    return int((np.where(Pm % 6 == 1, 1, -1).astype(np.int64) * cum[idx]).sum())


def S_counts(x):
    n1 = np.arange(7, x + 1, 6, dtype=np.int32); n5 = np.arange(5, x + 1, 6, dtype=np.int32)
    return int(np.count_nonzero(CLS[n1] == 2)), int(np.count_nonzero(CLS[n5] == 2))


print("Theorem 2.2 boundary/brute-force checks (S_1-S_5 direct vs (O+pi'(sqrt x))/2):")
for x in [5, 24, 25, 26, 35, 49, 100, 1000, 10007, 100003]:
    S1, S5 = S_counts(x)
    O = ordered_O(x); pp = pi_prime(math.isqrt(x))
    check((O + pp) % 2 == 0 and S1 - S5 == (O + pp) // 2, "identity at x=%d" % x)
    print("   x=%7d : S_1-S_5=%4d  O=%4d  pi'(sqrt x)=%3d  rhs=%4d  [PASS]  (holds for all x, not only x>=25)" % (x, S1 - S5, O, pp, (O + pp) // 2))

CLASSDATA = {}
for K in K_SCALES:
    x = 6 * K + 1
    n1 = np.arange(7, x + 1, 6, dtype=np.int32); n5 = np.arange(5, x + 1, 6, dtype=np.int32)
    c1 = [int(np.count_nonzero(CLS[n1] == j)) for j in range(1, 5)]
    c5 = [int(np.count_nonzero(CLS[n5] == j)) for j in range(1, 5)]
    CLASSDATA[K] = (c1, c5)
    M = MAT[K]
    check([sum(M[i]) for i in range(4)] == c5 and [sum(M[i][j] for i in range(4)) for j in range(4)] == c1, "marginals")
    print("x=%d chi-sums by Omega: %s" % (x, [c1[j] - c5[j] for j in range(4)]))
    s1 = n1[CLS[n1] == 2]; s5 = n5[CLS[n5] == 2]
    w1, b1 = wb(s1); w5, b5 = wb(s5)
    sq1 = int(np.count_nonzero((w1 == 1) & (b1 == 0))); sq5 = int(np.count_nonzero((w1 == 1) & (b1 == 2)))
    d11 = int(np.count_nonzero((w1 == 2) & (b1 == 0))); d55 = int(np.count_nonzero((w1 == 2) & (b1 == 2)))
    mixed = int(np.count_nonzero((w5 == 2) & (b5 == 1)))
    check(sq1 + sq5 + d11 + d55 == len(s1) and mixed == len(s5), "class-parity semiprime partition")
    O = ordered_O(x); pp = pi_prime(math.isqrt(x))
    check(pp == sq1 + sq5, "pi'(sqrt x) must equal the number of prime squares")
    check(c1[1] - c5[1] == (O + pp) // 2, "Theorem 2.2 at x=%d" % x)
    print("   semiprimes: (1,1)=%d (5,5)=%d mixed=%d p1^2=%d p5^2=%d ; O=%d pi'=%d ; S_1-S_5=%d" % (d11, d55, mixed, sq1, sq5, O, pp, c1[1] - c5[1]))
    print("   split A (note/Ford-Sneed form): O/2 = %.1f (cross) + pi'/2 = %.1f (squares)" % (O / 2, pp / 2))
    print("   split B (squares in full):     distinct-pair chi-sum (O-pi')/2 = %d  + squares pi' = %d" % ((O - pp) // 2, pp))
    check((O - pp) // 2 == d11 + d55 - mixed, "distinct-pair chi-sum")
    c3_1 = n1[CLS[n1] == 3]; c3_5 = n5[CLS[n5] == 3]
    pat = {}
    for tag, arr in [("right", c3_1), ("left", c3_5)]:
        w, b = wb(arr)
        for wi in (1, 2, 3):
            for bi in (0, 1, 2, 3):
                c = int(np.count_nonzero((w == wi) & (b == bi)))
                if c:
                    pat[(tag, wi, bi)] = c
    print("   3-almost patterns (side, omega, b): %s" % pat)
    check(all(b % 2 == (0 if t == "right" else 1) for (t, w, b) in pat), "3-almost parity")
    check(len(pat) == 10, "ten patterns")
    del n1, n5, s1, s5, c3_1, c3_5, w1, b1, w5, b5

NOTE_PAT_1E7 = {("right", 1, 0): 36, ("right", 2, 0): 84521, ("right", 2, 2): 138521, ("right", 3, 0): 522639,
                ("right", 3, 2): 2166195, ("left", 1, 3): 39, ("left", 2, 1): 85178, ("left", 2, 3): 139167,
                ("left", 3, 1): 1909337, ("left", 3, 3): 778617}
check(pat == NOTE_PAT_1E7, "ten-pattern census at 6e7+1 vs note")
print("ten-pattern census at x=6*10^7+1 identical to the note  [PASS]")

# Lemma 2.1: class-labelled shape determination, checked by full factorization for n <= 3*10^5
if HAVE_SYMPY:
    shapes = {}
    coll = {}
    for n in range(5, 300001):
        if n % 6 not in (1, 5):
            continue
        f = factorint(n)
        O = sum(f.values())
        if O > 4:
            continue
        w = len(f); b = sum(e for q, e in f.items() if q % 6 == 5)
        shape = tuple(sorted((q % 6, e) for q, e in f.items()))
        key = (O, w, b)
        shapes.setdefault(key, set()).add(shape)
        coll.setdefault((key, shape), n)
    bad3 = {k: v for k, v in shapes.items() if k[0] <= 3 and len(v) > 1}
    check(not bad3, "Lemma 2.1 fails for Omega<=3: %s" % bad3)
    n3 = sum(1 for k in shapes if 1 <= k[0] <= 3)
    print("Lemma 2.1: for n<=300000, gcd(n,6)=1, 1<=Omega<=3 the class-labelled shape is a function of (Omega,omega,b); %d keys  [PASS]" % n3)
    bad4 = {k: v for k, v in shapes.items() if k[0] == 4 and len(v) > 1}
    for k, v in sorted(bad4.items()):
        print("   Omega=4 collision at key %s: shapes %s, witnesses %s" % (k, sorted(v), [coll[(k, s)] for s in sorted(v)]))
    check((4, 3, 2) in bad4 and (4, 2, 0) in bad4, "expected Omega=4 collisions")
print("NB: the multiset of classes with multiplicity is {1^(Omega-b), 5^b} for EVERY n coprime to 6 --")
print("    trivially determined by (Omega,b) alone; Lemma 2.1's content is the class-LABELLED EXPONENT SHAPE.")

print()
print("SB6 independence prediction and splits (exact Fractions):")
for K in K_SCALES:
    c1, c5 = CLASSDATA[K]
    x = 6 * K + 1
    R1, C1, R2, C2 = c5[0], c1[0], c5[1], c1[1]
    pred = Fraction(R1 * C2 - R2 * C1, K)
    pibar = Fraction(R1 + C1, 2); dpi = Fraction(R1 - C1, 2); Sbar = Fraction(R2 + C2, 2); dS = Fraction(C2 - R2, 2)
    prim = 2 * Sbar * dpi / K; semi = 2 * pibar * dS / K
    check(pred == prim + semi, "split")
    O = ordered_O(x); pp = pi_prime(math.isqrt(x))
    sqA = pibar * pp / (2 * K); crA = pibar * O / (2 * K)
    check(sqA + crA == semi, "three-way split A")
    sqB = pibar * pp / K; dpB = pibar * (O - pp) / (2 * K)
    check(sqB + dpB == semi, "three-way split B")
    act = MAT[K][0][1] - MAT[K][1][0]
    print("   K=%d act %d pred %.2f | prime-Chebyshev %.2f | split A: squares %.2f (%.1f%%) cross %.2f | split B: squares %.2f (%.1f%%) distinct-pairs %.2f (%.1f%%)"
          % (K, act, float(pred), float(prim), float(sqA), 100 * float(sqA / pred), float(crA),
             float(sqB), 100 * float(sqB / pred), float(dpB), 100 * float(dpB / pred)))
    xs = float(x)
    print("      pred/(sqrt x/log x) = %.4f ; pred/(sqrt x * loglog x/log^2 x) = %.4f ; pi'(sqrt x)/(2 sqrt x/log x) = %.3f ; dpi/(sqrt x/log x) = %.3f"
          % (float(pred) / (math.sqrt(xs) / math.log(xs)), float(pred) / (math.sqrt(xs) * math.log(math.log(xs)) / math.log(xs) ** 2),
             pp / (2 * math.sqrt(xs) / math.log(xs)), float(dpi) / (math.sqrt(xs) / math.log(xs))))
    a = Fraction(MAT[K][0][1], MAT[K][1][0]); pr = Fraction(c5[0] * c1[1], c1[0] * c5[1])
    print("      ratio N_PS/N_SP = %.4f predicted %.4f" % (float(a), float(pr)))
    for (i, j, nm) in [(0, 2, "PC-CP"), (1, 2, "SC-CS")]:
        print("      %s pred %.2f" % (nm, float(Fraction(c5[i] * c1[j] - c5[j] * c1[i], K))))
print(el())

# ---------------------------------------------------------------- prime race
hr("3. PRIME RACE (SB5)")
c = np.cumsum(chi_p[PR >= 5], dtype=np.int32)
print("pi_chi over primes >= 5 up to %d: final %d ; prefixes <0: %d, =0: %d, >0: %d" % (NMAX, c[-1], (c < 0).sum(), (c == 0).sum(), (c > 0).sum()))
check(int((c > 0).sum()) == 0 and int((c == 0).sum()) == 9, "prime race prefix counts")
print("pi_chi(6*10^5+1)=%d pi_chi(6*10^6+1)=%d pi_chi(6*10^7+1)=%d" % (pi_chi(600001), pi_chi(6000001), pi_chi(60000001)))
del c

# ---------------------------------------------------------------- strata
hr("4. HOSTILE A: mod-35 strata of (6k-1,6k+1) -- full partition into 9 strata (SB9 / note 5(ii))")
h5m = ((6 * kk - 1) % 5 == 0); h5p = ((6 * kk + 1) % 5 == 0); h7m = ((6 * kk - 1) % 7 == 0); h7p = ((6 * kk + 1) % 7 == 0)
STR = {
    "none": ~(h5m | h5p | h7m | h7p),
    "cross 5L7R": h5m & h7p, "cross 7L5R": h7m & h5p,
    "5L only": h5m & ~h7m & ~h7p, "5R only": h5p & ~h7m & ~h7p,
    "7L only": h7m & ~h5m & ~h5p, "7R only": h7p & ~h5m & ~h5p,
    "35|left (OMITTED by explorer)": h5m & h7m, "35|right (OMITTED by explorer)": h5p & h7p,
}
del h5m, h5p, h7m, h7p
tot_mask = np.zeros(KMAX, bool)
for nm, m in STR.items():
    check(not bool((tot_mask & m).any()), "strata overlap")
    tot_mask |= m
check(bool(tot_mask.all()), "9 strata partition the centers")
del tot_mask
sizes35 = {nm: int(m[:35].sum()) for nm, m in STR.items()}
print("sizes on one block of 35: %s  (sum %d)" % (sizes35, sum(sizes35.values())))
check(sum(sizes35.values()) == 35 and sizes35["35|left (OMITTED by explorer)"] == 1, "block sizes")
for K in K_SCALES:
    tot = 0
    parts = {}
    for nm, m in STR.items():
        mk = m[:K]
        Ls = L0[:K][mk]; Rs = R0[:K][mk]
        ps = int(np.count_nonzero((Ls == 1) & (Rs == 2))); sp = int(np.count_nonzero((Ls == 2) & (Rs == 1)))
        pp_ = int(np.count_nonzero((Ls == 1) & (Rs == 1)))
        parts[nm] = (ps - sp, pp_, ps, sp)
        tot += ps - sp
    check(tot == MAT[K][0][1] - MAT[K][1][0], "strata sum to total")
    seven = sum(v[0] for nm, v in parts.items() if "OMITTED" not in nm)
    print("K=%d: 9-strata PS-SP sum = %d = total  [PASS]; the explorer's 7 strata sum to %d (off by %d);"
          % (K, tot, seven, tot - seven))
    print("      omitted strata (PS-SP, PP, PS, SP): 35|left %s, 35|right %s ; cross 5L7R PP=%d (the pair (5,7) at k=1 IS a prime pair)"
          % (parts["35|left (OMITTED by explorer)"], parts["35|right (OMITTED by explorer)"], parts["cross 5L7R"][1]))
    primes_left_none = int(np.count_nonzero(L0[:K][STR["none"][:K]] == 1))
    primes_left_all = int(np.count_nonzero(L0[:K] == 1))
    print("      left-endpoint primes in 'none' stratum: %d of %d (so 'the stratum carries all primes' is FALSE)" % (primes_left_none, primes_left_all))
print("   witness for the -1: k=6 gives (35,37) = (5*7, prime), an SP pair in the 35|left stratum.")
check(int(CLS[35]) == 2 and int(CLS[37]) == 1, "k=6 witness")
del STR

# ---------------------------------------------------------------- Lemma 2.3
hr("5. LEMMA 2.3 involution check on complete residue systems")
for (c, g, M) in [(-1, 2, 35), (1, 4, 35), (1, 6, 35), (-1, 8, 385), (1, 300, 385)]:
    inv6 = pow(6, -1, M)
    shift = ((2 * c + g) * inv6) % M
    ks = np.arange(M)
    kstar = (-ks - shift) % M
    check(bool(np.all((6 * kstar + c) % M == (-(6 * ks + c + g)) % M)) and bool(np.all((6 * kstar + c + g) % M == (-(6 * ks + c)) % M)), "involution map")
    check(bool(np.all(kstar[kstar] == ks)), "involution")
    print("   c=%d g=%d M=%d: k*=-k-%d sends (n,n+g) to (-(n+g),-n) mod M and is an involution  [PASS]" % (c, g, M, shift))

# ---------------------------------------------------------------- shifted pairs and pooled gaps
hr("6. SHIFTED PAIRS (SB8) AND POOLED-GAP CONTROL (SB7)")
for nm, (oa, ob), ref in [("(6k+1,6k+5)", (1, 5), [-5, -43, 461]), ("(6k-1,6k+5)", (-1, 5), [39, -170, 129]), ("(6k+1,6k+7)", (1, 7), [109, -187, 108])]:
    La = CLS[6 * kk + oa]; Rb = CLS[6 * kk + ob]
    vals = []
    for K in K_SCALES:
        vals.append(int(np.count_nonzero((La[:K] == 1) & (Rb[:K] == 2))) - int(np.count_nonzero((La[:K] == 2) & (Rb[:K] == 1))))
    check(vals == ref, "shifted pair %s" % nm)
    print("   %s PS-SP = %s  [PASS vs note]" % (nm, vals))
    del La, Rb

NOTE_POOL = {  # (K, orient) -> (mean, sd, npos, meanpred) for PS-SP
    (10**6, "[5|1]"): (100.7, 197.4, 32, 117.5), (10**6, "[1|5]"): (-113.0, 168.1, 14, -109.1),
    (10**6, "[5|5]"): (9.7, 215.1, 25, 4.1), (10**6, "[1|1]"): (12.5, 174.1, 23, 4.4),
    (10**7, "[5|1]"): (195.1, 541.2, 30, 234.2), (10**7, "[1|5]"): (-66.4, 493.3, 25, -225.6),
    (10**7, "[5|5]"): (96.7, 585.9, 28, 4.3), (10**7, "[1|1]"): (-32.7, 500.6, 24, 4.4),
}
NOTE_POOL_SC = {(10**6, "[5|1]"): -168.7, (10**6, "[1|5]"): 171.5, (10**6, "[5|5]"): 19.1, (10**6, "[1|1]"): -39.1,
                (10**7, "[5|1]"): -369.0, (10**7, "[1|5]"): 498.2, (10**7, "[5|5]"): 116.8, (10**7, "[1|1]"): -66.4}
NOTE_POOL_PC = {(10**7, "[1|5]"): -250.1}
for K in [10**6, 10**7]:
    stats = {}
    kK = kk[:K]
    for sc, off in [(5, -1), (1, 1)]:
        base = 6 * kK + off
        La = CLS[base]
        L1 = La == 1; L2 = La == 2; L3 = La == 3
        R1c = int(L1.sum()); R2c = int(L2.sum()); R3c = int(L3.sum())
        for g in range(2, GMAX + 1, 2):
            ec = (sc + g) % 6
            if ec == 3:
                continue
            ori = "[%d|%d]" % (sc, ec)
            Rb = CLS[base + g]
            Rb1 = Rb == 1; Rb2 = Rb == 2; Rb3 = Rb == 3
            C1c = int(Rb1.sum()); C2c = int(Rb2.sum())
            D = int(np.count_nonzero(L1 & Rb2)) - int(np.count_nonzero(L2 & Rb1))
            Dpc = int(np.count_nonzero(L1 & Rb3)) - int(np.count_nonzero(L3 & Rb1))
            Dsc = int(np.count_nonzero(L2 & Rb3)) - int(np.count_nonzero(L3 & Rb2))
            pred = Fraction(R1c * C2c - R2c * C1c, K)
            stats.setdefault(ori, []).append((g, D, pred, Dpc, Dsc))
        del base, La, L1, L2, L3
    for ori in ["[5|1]", "[1|5]", "[5|5]", "[1|1]"]:
        lst = stats[ori]
        n = len(lst)
        check(n == 50, "50 gaps per orientation")
        acts = [t[1] for t in lst]
        mean = sum(acts) / n
        sd = (sum((a - mean) ** 2 for a in acts) / (n - 1)) ** 0.5
        npos = sum(1 for a in acts if a > 0)
        mp = float(sum(t[2] for t in lst) / n)
        msc = sum(t[4] for t in lst) / n
        mpc = sum(t[3] for t in lst) / n
        ref = NOTE_POOL[(K, ori)]
        check(abs(mean - ref[0]) < 0.06 and abs(sd - ref[1]) < 0.06 and npos == ref[2] and abs(mp - ref[3]) < 0.06, "pooled PS-SP %s K=%d: %.1f %.1f %d %.1f" % (ori, K, mean, sd, npos, mp))
        check(abs(msc - NOTE_POOL_SC[(K, ori)]) < 0.06, "pooled SC-CS")
        if (K, ori) in NOTE_POOL_PC:
            check(abs(mpc - NOTE_POOL_PC[(K, ori)]) < 0.06, "pooled PC-CP outlier")
        se = sd / n ** 0.5
        print("   K=%d %s PS-SP mean %.1f sd %.1f se %.1f #>0 %d pred %.1f (mean-pred)/se=%.2f | SC-CS mean %.1f | PC-CP mean %.1f  [PASS vs note]"
              % (K, ori, mean, sd, se, npos, mp, (mean - mp) / se, msc, mpc))
    g2 = [t[1] for t in stats["[5|1]"] if t[0] == 2][0]
    rank = sum(1 for t in stats["[5|1]"] if t[1] > g2)
    print("   K=%d: original pair g=2 has D=%d, exceeded by %d of the other 49 [5|1] gaps" % (K, g2, rank))
print(el())

hr("7. SUMMARY OF AUDIT FINDINGS")
print("* All matrices, running-sign summaries, class censuses, pattern censuses, the Theorem 2.2 values,")
print("  shifted-pair and pooled-gap statistics agree with the explorer's .out and the draft note.")
print("* Theorem 2.2 holds for all x >= 1 (both sides vanish below 25); 'x>=25' is merely where it is nontrivial.")
print("* Lemma 2.1 is misworded: the class multiset is trivially (Omega,b)-determined for all n;")
print("  what is determined for Omega<=3 (and not for Omega=4) is the class-labelled exponent shape.")
print("* Note 5(ii): the seven listed mod-35 strata do not partition the centers; two strata (35|left, 35|right,")
print("  one center each per block of 35) are missing, and the 35|left stratum contains the SP pair (35,37),")
print("  so 'their sum reproduces the totals' is off by -1 at every scale; 'the cross strata contain no primes'")
print("  is false at k=1 (5,7); 'the none stratum carries all primes' is false (it carries all PP pairs only).")
print("* The 'square part carries about a third' statement is convention dependent: with the squares counted")
print("  in full, distinct-pair semiprimes favour class 5 and the squares carry more than the whole surplus.")
print("peak RSS %.2f GB ; total %s" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, el()))
print("ALL AUDIT CHECKS PASSED")
