#!/usr/bin/env python3
"""
collatz_mod6_20260917_cell_ordering_scale_audit_proof-audit.py

INDEPENDENT recomputation for the proof-audit of lane cell_ordering_scale.
Does not import the explorer's script.  Different Omega method: smallest-prime-
factor sieve + repeated division (the explorer adds 1 on prime-power progressions).
Different census code: full-length cumulative arrays instead of chunked one-hot.
Selberg-Delange coefficients by hand from Mertens/Euler/prime-zeta values instead
of mp.taylor of the full product.

Everything load-bearing raises on failure.  Expected values are the explorer's
frozen numbers; a mismatch is printed as AUDIT-FAIL, not hidden.
"""
import math
import sys
import time
import itertools
from fractions import Fraction

import numpy as np

T0 = time.time()
KMAX = 10**7
NMAX = 6 * KMAX + 1
DEC = [10**3, 10**4, 10**5, 10**6, 10**7]
FAILS = []


def check(cond, msg):
    if cond:
        print("  ok   ", msg)
    else:
        print("  AUDIT-FAIL", msg)
        FAILS.append(msg)


def t():
    return "[%.1fs]" % (time.time() - T0)


# ---------------------------------------------------------------- spf sieve
print("== 0. smallest-prime-factor sieve to", NMAX, t())
spf = np.zeros(NMAX + 1, dtype=np.int32)
r = math.isqrt(NMAX)
for p in range(2, r + 1):
    if spf[p] == 0:
        v = spf[p * p::p]
        v[v == 0] = p
print("spf sieve done", t())


def _omega_block(n):
    """Omega for an int64 array n>=1 via repeated division by spf (0 => prime)."""
    n = n.astype(np.int64).copy()
    om = np.zeros(len(n), dtype=np.uint8)
    while True:
        idx = np.flatnonzero(n > 1)
        if len(idx) == 0:
            break
        sub = n[idx]
        q = spf[sub].astype(np.int64)
        q[q == 0] = sub[q == 0]
        n[idx] = sub // q
        om[idx] += 1
    return om


def omega_arr(n, blk=2 * 10**6):
    """Chunked wrapper so int64 temporaries stay ~16 MB each (RAM budget < 1 GB)."""
    return np.concatenate([_omega_block(n[s:s + blk]) for s in range(0, len(n), blk)])


def omega_trial(n):
    c = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            c += 1
        d += 1
    return c + (1 if n > 1 else 0)


# global pi_k counts at 10^6, 10^7, 6*10^7 (chunked over all n)
bounds = [10**6, 10**7, 6 * 10**7]
GLOB = {}
acc = np.zeros(64, dtype=np.int64)
lo = 1
for hi in [10**6, 10**7] + [10**7 * j for j in range(2, 7)]:
    om = omega_arr(np.arange(lo, hi + 1, dtype=np.int64))
    acc += np.bincount(om, minlength=64)
    if hi in bounds:
        GLOB[hi] = acc.copy()
    lo = hi + 1
    del om
print("global Omega histogram done", t())
check(GLOB[10**6][1] == 78498, "pi(10^6)=78498")
check(GLOB[10**7][1] == 664579, "pi(10^7)=664579")
check(GLOB[6 * 10**7][1] == 3562115, "pi(6e7)=3562115 (explorer prime count)")
check(GLOB[10**6][2] == 210035, "pi_2(10^6)=210035")
check(GLOB[10**6][3] == 250853 and GLOB[10**6][4] == 198062, "pi_3(10^6)=250853, pi_4(10^6)=198062")
check(GLOB[10**7][2] == 1904324 and GLOB[10**7][3] == 2444359 and GLOB[10**7][4] == 2050696, "pi_k(10^7) k=2,3,4")
check(GLOB[6 * 10**7][2] == 10655932 and GLOB[6 * 10**7][3] == 14333707 and GLOB[6 * 10**7][4] == 12525750, "pi_k(6e7) k=2,3,4")
# GLOB[1] at 6e7 counts n=1 as Omega=0; fine.

BLK = 2 * 10**6
OmL = np.concatenate([omega_arr(6 * np.arange(s, min(s + BLK - 1, KMAX) + 1, dtype=np.int64) - 1) for s in range(1, KMAX + 1, BLK)])
OmR = np.concatenate([omega_arr(6 * np.arange(s, min(s + BLK - 1, KMAX) + 1, dtype=np.int64) + 1) for s in range(1, KMAX + 1, BLK)])
check(len(OmL) == KMAX and len(OmR) == KMAX, "endpoint array lengths")
del spf
print("endpoint Omega done", t())
# audits of my own Omega
for k in range(1, 5001):
    if OmL[k - 1] != omega_trial(6 * k - 1) or OmR[k - 1] != omega_trial(6 * k + 1):
        raise RuntimeError("trial division mismatch at k=%d" % k)
print("  ok    trial-division audit k<=5000")
try:
    from sympy import factorint
    rng0 = np.random.default_rng(917)
    for k in rng0.integers(1, KMAX + 1, size=200).tolist():
        if OmL[k - 1] != sum(factorint(6 * k - 1).values()) or OmR[k - 1] != sum(factorint(6 * k + 1).values()):
            raise RuntimeError("sympy mismatch at k=%d" % k)
    print("  ok    sympy.factorint audit on 200 random k")
except ImportError:
    print("  (sympy unavailable)")

Lc = np.minimum(OmL, 4).astype(np.int32)
Rc = np.minimum(OmR, 4).astype(np.int32)
cell = ((Lc - 1) * 4 + (Rc - 1)).astype(np.int32)
del OmL, OmR


def mat(K, c=cell):
    return np.bincount(c[:K], minlength=16).reshape(4, 4)


EXP = {
    10**3: [[142, 182, 60, 13], [167, 192, 71, 7], [65, 71, 7, 2], [10, 9, 2, 0]],
    10**4: [[810, 1338, 687, 206], [1323, 1942, 915, 231], [685, 935, 367, 62], [196, 233, 61, 9]],
    10**5: [[5330, 10091, 6537, 2614], [10050, 17794, 10719, 3851], [6561, 10751, 5600, 1664], [2583, 3844, 1666, 345]],
    10**6: [[37915, 78689, 59706, 30192], [78277, 157420, 112416, 52182], [59992, 112305, 72363, 28755], [30161, 52125, 28736, 8766]],
    10**7: [[280557, 632905, 539140, 328636], [632766, 1380659, 1118019, 632028], [539061, 1118744, 837851, 416682], [328491, 631688, 416902, 165871]],
}
print("== C6. 4x4 matrices at the decades")
MATS = {}
for K in DEC:
    MATS[K] = mat(K)
    check(MATS[K].tolist() == EXP[K], "4x4 matrix at K=%d equals explorer/inherited table" % K)
M7 = MATS[10**7]
check(M7.sum(axis=1).tolist() == [1781238, 3763472, 2912338, 1542952], "row marginals 10^7")
check(M7.sum(axis=0).tolist() == [1780875, 3763996, 2911912, 1543217], "col marginals 10^7")
check(int(M7.sum(axis=1)[0] + M7.sum(axis=0)[0]) + 2 == GLOB[6 * 10**7][1], "primes 5 mod 6 + primes 1 mod 6 + {2,3} = pi(6e7)")

# ---------------------------------------------------------------- rankings
NAMES9 = {(i, j): "PSC"[i - 1] + "PSC"[j - 1] for i in (1, 2, 3) for j in (1, 2, 3)}


def rank_str(items):
    items = sorted(items, key=lambda x: -x[1])
    s = items[0][0]
    for a, b in zip(items, items[1:]):
        s += (" = " if a[1] == b[1] else " > ") + b[0]
    return s


def six(M):
    return [("PP", int(M[0, 0])), ("PS+SP", int(M[0, 1] + M[1, 0])), ("SS", int(M[1, 1])),
            ("PC+CP", int(M[0, 2] + M[2, 0])), ("SC+CS", int(M[1, 2] + M[2, 1])), ("CC", int(M[2, 2]))]


print("== C9/C13. rankings at decades")
EXP9 = {10**3: "SS > PS > SP > PP > SC = CS > CP > PC > CC",
        10**4: "SS > PS > SP > CS > SC > PP > PC > CP > CC",
        10**5: "SS > CS > SC > PS > SP > CP > PC > CC > PP",
        10**6: "SS > SC > CS > PS > SP > CC > CP > PC > PP",
        10**7: "SS > CS > SC > CC > PS > SP > PC > CP > PP"}
EXP6 = {10**3: "PS+SP > SS > PP = SC+CS > PC+CP > CC",
        10**4: "PS+SP > SS > SC+CS > PC+CP > PP > CC",
        10**5: "SC+CS > PS+SP > SS > PC+CP > CC > PP",
        10**6: "SC+CS > SS > PS+SP > PC+CP > CC > PP",
        10**7: "SC+CS > SS > PS+SP > PC+CP > CC > PP"}
for K in DEC:
    M = MATS[K]
    r9 = rank_str([(NAMES9[(i, j)], int(M[i - 1, j - 1])) for (i, j) in NAMES9])
    r6 = rank_str(six(M))
    print("  K=%-8d nine: %s" % (K, r9))
    print("  K=%-8d six : %s" % (K, r6))
    check(r9 == EXP9[K], "nine-cell ranking K=%d" % K)
    check(r6 == EXP6[K], "six-class ranking K=%d" % K)
    # position of SS and CC in the six-class ranking (draft note claims 'SS first, CC fourth of six')
    order6 = [n for n, _ in sorted(six(M), key=lambda x: -x[1])]
    print("     six-class position of SS = %d, of CC = %d" % (order6.index("SS") + 1, order6.index("CC") + 1))
for K, exp in [(100, "PS > PP > SP > SS"), (270, "PS > SP > PP > SS"), (500, "PS > SP > SS > PP"), (1000, "SS > PS > SP > PP")]:
    M = mat(K)
    r = rank_str([("PP", int(M[0, 0])), ("PS", int(M[0, 1])), ("SP", int(M[1, 0])), ("SS", int(M[1, 1]))])
    check(r == exp, "2x2 ranking K=%d: %s" % (K, r))
M272 = mat(272)
check(int(M272[2, 2]) == 1 and int(M272[1, 1]) == 43, "N33(272)=1, N22(272)=43")
M270 = mat(270)
print("  N33(270)=%d N22(270)=%d ; exp(e^2)/6 = %.2f" % (M270[2, 2], M270[1, 1], math.exp(math.exp(2)) / 6))

# ---------------------------------------------------------------- cumulative census
print("== C7/C8/C9. running census over all K", t())
CUM = {}
for (i, j) in NAMES9:
    CUM[(i, j)] = np.cumsum((cell == (i - 1) * 4 + (j - 1)).astype(np.int32), dtype=np.int32)
PP, PS, SP, SS = CUM[(1, 1)], CUM[(1, 2)], CUM[(2, 1)], CUM[(2, 2)]
PC, CP, SC, CS, CC = CUM[(1, 3)], CUM[(3, 1)], CUM[(2, 3)], CUM[(3, 2)], CUM[(3, 3)]
check(int(PP[-1]) == 280557 and int(SS[-1]) == 1380659 and int(CC[-1]) == 837851, "cumulative totals")

v1 = ~(SS > np.maximum(PS, SP))
v2 = ~(np.minimum(PS, SP) > PP)
n1, last1 = int(v1.sum()), int(np.flatnonzero(v1)[-1]) + 1
n2, last2 = int(v2.sum()), int(np.flatnonzero(v2)[-1]) + 1
print("  SS>max(PS,SP) violations: %d, last K=%d ; all K<=761 violate: %s" % (n1, last1, bool(v1[:761].all())))
print("  min(PS,SP)>PP violations: %d, last K=%d ; non-violating K<160: %s" % (n2, last2, (np.flatnonzero(~v2[:160]) + 1).tolist()))
check(n1 == 761 and last1 == 761, "C7: 761 violations of SS>max(PS,SP), last at 761")
check(n2 == 151 and last2 == 160, "C7: 151 violations of min(PS,SP)>PP, last at 160")

EXPM = {"PS-SP": (139, 14947, 992, 650, 5512997), "PC-CP": (79, 12330, 639, 570, 7331193), "SC-CS": (-725, 4654, 470, 1323, 9539838)}
for name, (A, B) in [("PS-SP", (PS, SP)), ("PC-CP", (PC, CP)), ("SC-CS", (SC, CS))]:
    D = A - B
    zeros = int((D == 0).sum())
    sg = np.sign(D)
    nz = sg[sg != 0]
    flips = int((nz[1:] != nz[:-1]).sum())
    am = int(np.argmax(np.abs(D)))
    fin = int(D[-1])
    tot = int(A[-1] + B[-1])
    print("  %s final=%+d sqrt(tot)=%.0f zeros=%d flips=%d max|D|=%d at K=%d (first K with D=0: %d)"
          % (name, fin, math.sqrt(tot), zeros, flips, abs(int(D[am])), am + 1, int(np.flatnonzero(D == 0)[0]) + 1))
    e = EXPM[name]
    check((fin, zeros, flips, abs(int(D[am])), am + 1) == e, "C8 mirror stats %s" % name)

# crossovers: last K with A<=B
CR = {("SS", "PS"): (SS, PS, 761), ("SS", "SP"): (SS, SP, 672), ("SS", "PS+SP"): (SS, PS + SP, 920807),
      ("SC+CS", "SS"): (SC + CS, SS, 16815), ("SC+CS", "PS+SP"): (SC + CS, PS + SP, 66959),
      ("CC", "PS"): (CC, PS, 1591506), ("CC", "SP"): (CC, SP, 1559904), ("PC+CP", "PP"): (PC + CP, PP, 1445),
      ("CC", "PP"): (CC, PP, 83683)}
for (a, b), (A, B, exp) in CR.items():
    le = A <= B
    last = int(np.flatnonzero(le)[-1]) + 1
    check(last == exp, "C9 crossover %s > %s for good after K=%d (expected %d)" % (a, b, last, exp))
for (a, b), (A, B) in {("CC", "SS"): (CC, SS), ("CC", "PS+SP"): (CC, PS + SP), ("CC", "PC+CP"): (CC, PC + CP)}.items():
    check(bool((A <= B).all()), "C9 %s never exceeds %s below 10^7" % (a, b))

# ranking-change census (weak order incl. ties) via pairwise sign vectors, full arrays
def change_census(groups):
    G = [CUM[g[0]] if len(g) == 1 else sum(CUM[c] for c in g) for g in groups]   # int32 views/sums
    ng = len(G)
    changed = np.zeros(KMAX, dtype=bool)
    for a in range(ng):
        for b in range(a + 1, ng):
            s = np.sign(G[a] - G[b]).astype(np.int8)
            changed[1:] |= (s[1:] != s[:-1])
            del s
    idx = np.flatnonzero(changed)
    return int(len(idx)), (int(idx[-1]) + 1 if len(idx) else 0)


nine_g = [[(i, j)] for (i, j) in NAMES9]
six_g = [[(1, 1)], [(1, 2), (2, 1)], [(2, 2)], [(1, 3), (3, 1)], [(2, 3), (3, 2)], [(3, 3)]]
two_g = [[(1, 1)], [(1, 2)], [(2, 1)], [(2, 2)]]
twom_g = [[(1, 1)], [(1, 2), (2, 1)], [(2, 2)]]
for name, g, exp in [("nine", nine_g, (8918, 8714402)), ("six", six_g, (259, 920808)), ("2x2", two_g, (3885, 8714402)), ("2x2m", twom_g, (33, 920808))]:
    res = change_census(g)
    check(res == exp, "ranking-change census %s: changes=%d last=%d (expected %s)" % (name, res[0], res[1], exp))
print(t())

# ---------------------------------------------------------------- independence ratios, shuffle
print("== C12. independence ratios")
Mf = M7.astype(np.float64)
Rm = Mf * KMAX / np.outer(Mf.sum(1), Mf.sum(0))
print("  R at 10^7:\n" + "\n".join("   " + " ".join("%.4f" % v for v in row) for row in Rm))
check(abs(Rm[0, 0] - 0.8844) < 6e-5 and abs(Rm[1, 1] - 0.9746) < 6e-5 and abs(Rm[2, 2] - 0.9880) < 6e-5 and abs(Rm[3, 3] - 0.6966) < 6e-5, "R diagonal at 10^7")
check(abs(Rm[0, 1] - 0.9440) < 6e-5 and abs(Rm[0, 2] - 1.0394) < 6e-5 and abs(Rm[1, 2] - 1.0202) < 6e-5 and abs(Rm[0, 3] - 1.1955) < 6e-5, "R off-diagonal at 10^7")
rpp = []
for K in DEC:
    Mk = MATS[K].astype(np.float64)
    rpp.append(Mk[0, 0] * K / (Mk.sum(1)[0] * Mk.sum(0)[0]))
print("  R_PP decades:", ["%.4f" % v for v in rpp])
check(all(abs(a - b) < 6e-5 for a, b in zip(rpp, [0.9315, 0.8837, 0.8845, 0.8898, 0.8844])), "R_PP across decades")
# twin prime constant with an own sieve to 10^6
sv = np.ones(10**6 + 1, dtype=bool)
sv[:2] = False
for p in range(2, 1001):
    if sv[p]:
        sv[p * p::p] = False
pr = np.flatnonzero(sv)
C2 = float(np.prod(1 - 1.0 / (pr[1:].astype(np.float64) - 1) ** 2))
print("  C_2 = %.7f, 4C_2/3 = %.5f ; R_PP(10^7)/(4C_2/3) = %.4f" % (C2, 4 * C2 / 3, Rm[0, 0] / (4 * C2 / 3)))
check(abs(C2 - 0.6601618) < 3e-6, "C_2")
# HL derivation: twins (6k-1,6k+1), k<=K ~ 2C2*(6K)/log^2(6K); product (3K/log 6K)^2/K = 9K/log^2 6K; ratio 12C2/9=4C2/3
check(abs(Fraction(12, 9) - Fraction(4, 3)) == 0, "4C_2/3 algebra: 2C_2*6 / 9 = 4C_2/3")
del CUM, PP, PS, SP, SS, PC, CP, SC, CS, CC, v1, v2
rng = np.random.default_rng(20260917)
Rsh = Rc[rng.permutation(KMAX)]
csh = ((Lc - 1) * 4 + (Rsh - 1)).astype(np.int32)
Msh = mat(KMAX, csh).astype(np.float64)
Rmsh = Msh * KMAX / np.outer(Msh.sum(1), Msh.sum(0))
zmax = float(np.max(np.abs(Rmsh - 1) * np.sqrt(Msh)))
print("  shuffle control max |R-1|sqrt(N) = %.2f ; shuffled R_PP = %.4f" % (zmax, Rmsh[0, 0]))
check(abs(zmax - 1.33) < 0.01, "shuffle control z-max 1.33")

# ---------------------------------------------------------------- 3x3 conditionals
print("== C13. 3x3 conditionals")
for K in DEC:
    blk = MATS[K][:3, :3]
    tot = int(blk.sum())
    fr = blk / tot
    print("  K=%d block=%d (%.4f)  " % (K, tot, tot / K) + " ".join("%.5f" % v for v in fr.reshape(9)))
    check(sum(Fraction(int(v), tot) for v in blk.reshape(9)) == 1, "fractions sum to 1 at K=%d" % K)
blk7 = MATS[10**7][:3, :3]
check(int(blk7.sum()) == 7079702, "block total 7079702 at 10^7")
fr7 = blk7 / blk7.sum()
check(abs(fr7[0, 0] - 0.03963) < 6e-6 and abs(fr7[1, 1] - 0.19502) < 6e-6 and abs(fr7[2, 2] - 0.11835) < 6e-6 and abs(fr7[1, 2] - 0.15792) < 6e-6, "3x3 fractions at 10^7")
blk6 = MATS[10**6][:3, :3]
check(int(blk6.sum()) == 769083, "block total 769083 at 10^6")
fr6 = blk6 / blk6.sum()
check(abs(fr6[0, 0] - 0.04930) < 6e-6 and abs(fr6[1, 1] - 0.20469) < 6e-6 and abs(fr6[2, 2] - 0.09409) < 6e-6, "3x3 fractions at 10^6")
check([round(int(MATS[K][:3, :3].sum()) / K, 3) for K in DEC] == [0.957, 0.900, 0.834, 0.769, 0.708], "block fractions")

# ---------------------------------------------------------------- Selberg-Delange by hand
print("== C4. Selberg-Delange coefficients by hand")
from mpmath import mp, mpf, euler, mertens, primezeta, pi as mpi, rgamma, taylor

mp.dps = 30
h1 = mertens - euler                       # sum_p [log(1-1/p)+1/p]
P2, P3, P4 = primezeta(2), primezeta(3), primezeta(4)
h2, h3, h4 = P2 / 2, P3 / 3, P4 / 4        # log H(1,z) = h1 z + h2 z^2 + h3 z^3 + ...  (z^m coefficient = P(m)/m for m>=2)
# check h1 independently: -sum_{m>=2} P(m)/m
h1_alt = -mp.fsum(primezeta(m) / m for m in range(2, 80))
check(abs(h1 - h1_alt) < mpf(10) ** -20, "sum_p[log(1-1/p)+1/p] = M - gamma (Mertens' theorem, numerically to 1e-20)")
a = taylor(rgamma, 0, 5)                   # 1/Gamma(z) coefficients
check(abs(a[1] - 1) < mpf(10) ** -25 and abs(a[2] - euler) < mpf(10) ** -25, "1/Gamma(z) = z + gamma z^2 + ...")
check(abs(a[3] - (euler**2 / 2 - mpi**2 / 12)) < mpf(10) ** -25, "[z^3] 1/Gamma = gamma^2/2 - pi^2/12")
# exp(h1 z + h2 z^2 + h3 z^3 + h4 z^4) coefficients e0..e4
e0 = mpf(1)
e1 = h1
e2 = h1**2 / 2 + h2
e3 = h1**3 / 6 + h1 * h2 + h3
e4 = h1**4 / 24 + h1**2 * h2 / 2 + h2**2 / 2 + h1 * h3 + h4
E = [e0, e1, e2, e3, e4]
f = [sum(a[i] * E[m - i] for i in range(0, m + 1)) for m in range(0, 5)]
print("  f =", [mp.nstr(v, 12) for v in f])
check(abs(f[0]) < mpf(10) ** -25 and abs(f[1] - 1) < mpf(10) ** -25, "f_0=0, f_1=1")
check(abs(f[2] - mertens) < mpf(10) ** -25, "f_2 = M exactly (algebra: gamma + (M - gamma))")
check(abs(f[3] - mpf("-0.562152927240038")) < mpf(10) ** -12 and abs(f[4] - mpf("0.305977961616984")) < mpf(10) ** -12, "f_3, f_4 numerical")
# F_6 = F * (1 - 5z/6 + z^2/6)
g = [f[m] - mpf(5) / 6 * (f[m - 1] if m >= 1 else 0) + (f[m - 2] if m >= 2 else 0) / 6 for m in range(0, 5)]
print("  g =", [mp.nstr(v, 12) for v in g])
check(abs(g[2] - (mertens - mpf(5) / 6)) < mpf(10) ** -25, "g_2 = M - 5/6")
check(abs(g[2] - mpf("-0.571836120485691")) < mpf(10) ** -12 and abs(g[3] - mpf("-0.613400604613073")) < mpf(10) ** -12, "g_2, g_3 numerical")
print("  M = %s  M-5/6 = %s" % (mp.nstr(mertens, 12), mp.nstr(mertens - mpf(5) / 6, 12)))


def sdpoly(c, k, L):
    return float(sum(c[m] * mpf(L) ** (k - m) / mp.factorial(k - m) for m in range(0, k + 1)))


# C5 Landau/SD ratios
print("== C5. Landau and SD ratios")
EXP5 = {(10**6, 2): (1.105, 1.005), (10**7, 2): (1.104, 1.009), (6 * 10**7, 2): (1.102, 1.011),
        (10**6, 3): (1.005, 0.970), (10**7, 3): (1.020, 0.978), (6 * 10**7, 3): (1.028, 0.983),
        (10**6, 4): (0.907, 0.996), (10**7, 4): (0.923, 0.991), (6 * 10**7, 4): (0.934, 0.990)}
for x in bounds:
    L = math.log(math.log(x))
    for k in (2, 3, 4):
        act = int(GLOB[x][k])
        la = x / math.log(x) * L ** (k - 1) / math.factorial(k - 1)
        sd = x / math.log(x) * sdpoly(f, k, L)
        ok = abs(act / la - EXP5[(x, k)][0]) < 6e-4 and abs(act / sd - EXP5[(x, k)][1]) < 6e-4
        check(ok, "global pi_%d(%d): actual/Landau=%.4f actual/SD=%.4f" % (k, x, act / la, act / sd))
EXPm = {2: ([0.586, 0.675, 0.727, 0.758, 0.779], [0.796, 0.886, 0.933, 0.957, 0.971]),
        3: ([0.180, 0.261, 0.325, 0.377, 0.418], [0.859, 0.844, 0.868, 0.894, 0.915])}
for i in (2, 3):
    la_r, sd_r = [], []
    for K in DEC:
        x = 6 * K
        L = math.log(math.log(x))
        Ni = int(MATS[K].sum(1)[i - 1])
        la_r.append(Ni / (3 * K / math.log(x) * L ** (i - 1) / math.factorial(i - 1)))
        sd_r.append(Ni / (3 * K / math.log(x) * sdpoly(g, i, L)))
    print("  i=%d N_i./Landau %s  N_i./SD-F6 %s" % (i, ["%.3f" % v for v in la_r], ["%.3f" % v for v in sd_r]))
    check(all(abs(a - b) < 6e-4 for a, b in zip(la_r, EXPm[i][0])) and all(abs(a - b) < 6e-4 for a, b in zip(sd_r, EXPm[i][1])), "class-marginal ratios i=%d" % i)

# ---------------------------------------------------------------- C3 Landau weights, C11 models
print("== C3/C11. product-model predictions")
L2 = 2.0
w = lambda i, j, L: L ** (i + j - 2) / (math.factorial(i - 1) * math.factorial(j - 1))
check(abs(w(3, 3, L2) - w(2, 3, L2)) < 1e-12 and abs(w(2, 3, L2) - w(2, 2, L2)) < 1e-12, "Landau thresholds (3,3)=(2,3)=(2,2) at L=2")
L = 2.5
order = sorted([("CC", w(3, 3, L)), ("CS", w(2, 3, L)), ("SS", w(2, 2, L)), ("CP", w(1, 3, L)), ("PS", w(1, 2, L)), ("PP", w(1, 1, L))], key=lambda x: -x[1])
check([n for n, _ in order] == ["CC", "CS", "SS", "CP", "PS", "PP"], "Landau order for L>2 (L=2.5): " + " > ".join(n for n, _ in order))
FIT = [10**3, 2 * 10**3, 5 * 10**3, 10**4, 2 * 10**4, 5 * 10**4, 10**5, 2 * 10**5, 5 * 10**5, 10**6, 2 * 10**6, 5 * 10**6, 10**7]
rows = []
for K in FIT:
    M = mat(K)
    x = 6 * K
    Lv = math.log(math.log(x))
    r = M.sum(1)
    c = M.sum(0)
    N33, N22 = int(M[2, 2]), int(M[1, 1])
    rm, rp = r[2] / r[1], c[2] / c[1]
    c33 = N33 * K / (r[2] * c[2])
    c22 = N22 * K / (r[1] * c[1])
    rows.append((K, Lv, N33, N22, N33 / N22, rm, rp, c33 / c22))
for K, Lv, N33, N22, rat, rm, rp, cr in rows:
    if K in DEC:
        print("  K=%-8d L=%.4f N33/N22=%.4f r-=%.4f r+=%.4f c33/c22=%.4f" % (K, Lv, rat, rm, rp, cr))
r7 = rows[-1]
check(abs(r7[4] - 0.6069) < 6e-5 and abs(r7[5] - 0.774) < 6e-4 and abs(r7[7] - 1.0137) < 6e-5, "10^7 ratios N33/N22, r-, c33/c22")
check(abs(rows[-4][7] - 0.987) < 6e-4 and rows[-4][0] == 10**6, "c33/c22 at 10^6 = 0.987")


def rth(Lv):
    return sdpoly(g, 3, Lv) / sdpoly(g, 2, Lv)


def solve(fn, lo, hi):
    from scipy.optimize import brentq
    return brentq(fn, lo, hi, xtol=1e-10)


try:
    import scipy  # noqa
    have_scipy = True
except ImportError:
    have_scipy = False


def bis(fn, lo, hi):
    if have_scipy:
        return solve(fn, lo, hi)
    for _ in range(200):
        mid = (lo + hi) / 2
        if (fn(mid) > 0) == (fn(lo) > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


c_last, c_prev = r7[7], rows[-4][7]
LB = bis(lambda Lv: c_last * rth(Lv) ** 2 - 1, 2.0, 8.0)
sel = [r for r in rows if r[0] >= 10**5]
Ls = np.array([r[1] for r in sel])
ac, bc = np.polyfit(Ls, np.array([math.sqrt(r[5] * r[6]) for r in sel]), 1)
LC = bis(lambda Lv: c_last * (ac * Lv + bc) ** 2 - 1, 1.0, 12.0)
ad, bd = np.polyfit(Ls, np.array([math.sqrt(r[4]) for r in sel]), 1)
LD = (1 - bd) / ad
L6, L7 = rows[-4][1], r7[1]
slope = (c_last - c_prev) / (1 / L7 - 1 / L6)
cofL = lambda Lv: c_last + slope * (1 / Lv - 1 / L7)
LE = bis(lambda Lv: cofL(Lv) * rth(Lv) ** 2 - 1, 2.0, 8.0)
Kof = lambda Lv: math.exp(math.exp(Lv)) / 6
for tag, Lv, expL, expK in [("B", LB, 3.159, 2.8e9), ("C", LC, 3.220, 1.2e10), ("D", LD, 3.187, 5.4e9), ("E", LE, 3.127, 1.3e9)]:
    Kst = Kof(Lv)
    check(abs(Lv - expL) < 6e-4 and abs(math.log10(Kst) - math.log10(expK)) < 0.02, "Model %s: L*=%.4f K*=%.3e (log10 %.2f)" % (tag, Lv, Kst, math.log10(Kst)))
print("  fit C: %.3fL%+.3f ; fit D: %.3fL%+.3f ; c_inf(E)=%.4f ; r_theory(L7)=%.4f vs r-=%.4f" % (ac, bc, ad, bd, cofL(1e9), rth(L7), r7[5]))
check(abs(rth(L7) - 0.8209) < 6e-4, "SD marginal ratio overshoot 0.821 at 10^7")
wS = lambda i, Lv: sdpoly(g, i, Lv)
Lcc_ss = bis(lambda Lv: wS(3, Lv) - wS(2, Lv), 2.0, 9.0)
Lcc_sccs = bis(lambda Lv: wS(3, Lv) - 2 * wS(2, Lv), 2.0, 9.0)
Lss_pssp = bis(lambda Lv: wS(2, Lv) - 2 * wS(1, Lv), 2.0, 9.0)
for tag, Lv, expL in [("CC>SS", Lcc_ss, 3.1699), ("CC>SC+CS", Lcc_sccs, 4.9285), ("SS>PS+SP", Lss_pssp, 2.5718)]:
    check(abs(Lv - expL) < 6e-4, "SD-F6 c=1 crossing %s at L=%.4f, x~10^%.1f" % (tag, Lv, math.exp(Lv) / math.log(10)))
# six-class model predictions vs actual at decades
def merged(wv):
    return [("PP", wv[1] * wv[1]), ("PS+SP", 2 * wv[1] * wv[2]), ("SS", wv[2] * wv[2]), ("PC+CP", 2 * wv[1] * wv[3]), ("SC+CS", 2 * wv[2] * wv[3]), ("CC", wv[3] * wv[3])]


sd_match, mg_match, la_match = [], [], []
for K in DEC:
    Lv = math.log(math.log(6 * K))
    act = rank_str(six(MATS[K]))
    rs = rank_str(merged({i: wS(i, Lv) for i in (1, 2, 3)}))
    rl = rank_str(merged({i: Lv ** (i - 1) / math.factorial(i - 1) for i in (1, 2, 3)}))
    rr, cc = MATS[K].sum(1), MATS[K].sum(0)
    rm = rank_str([("PP", rr[0] * cc[0]), ("PS+SP", rr[0] * cc[1] + rr[1] * cc[0]), ("SS", rr[1] * cc[1]), ("PC+CP", rr[0] * cc[2] + rr[2] * cc[0]), ("SC+CS", rr[1] * cc[2] + rr[2] * cc[1]), ("CC", rr[2] * cc[2])])
    sd_match.append(rs == act)
    mg_match.append(rm == act)
    la_match.append(rl == act)
    print("  K=%d Landau: %s" % (K, rl))
check(sd_match == [False, True, False, True, True], "SD-F6 six-class match at 10^4, 10^6, 10^7 only")
check(mg_match == [False, True, True, False, True], "measured-marginal c=1 match at 10^4, 10^5, 10^7 only")
check(la_match == [False] * 5, "Landau order matches at no decade")

# ---------------------------------------------------------------- C1 tournaments
print("== C1. 4-vertex tournaments")
pairs = list(itertools.combinations(range(4), 2))
scores = set()
for orient in itertools.product((0, 1), repeat=6):
    s = [0] * 4
    for (u, v), o in zip(pairs, orient):
        s[u if o else v] += 1
    scores.add(tuple(sorted(s)))
scores = sorted(scores)
print("  score sequences:", scores)
check(scores == [(0, 1, 2, 3), (0, 2, 2, 2), (1, 1, 1, 3), (1, 1, 2, 2)], "Landau n=4 score sequences")
check(not any(len(set(s)) == 3 and s[1] == s[2] for s in scores), "no (min, mid, mid, max) profile")
# parity proof: scores sum to 6; (a,m,m,b) with a<m<b integers in [0,3] -> a+2m+b=6 has no solution
sols = [(a, m, b) for a in range(4) for m in range(4) for b in range(4) if a < m < b and a + 2 * m + b == 6]
check(sols == [], "arithmetic proof: a+2m+b=6 with 0<=a<m<b<=3 has no solution (no enumeration needed)")
# Landau's theorem check on the four sequences
def landau_ok(s):
    s = sorted(s)
    return all(sum(s[:k]) >= k * (k - 1) // 2 for k in range(1, 5)) and sum(s) == 6
allseq = [s for s in itertools.combinations_with_replacement(range(4), 4) if landau_ok(s)]
check(sorted(allseq) == scores, "Landau (1953) criterion reproduces exactly the enumerated set")

print()
print("AUDIT-FAIL count:", len(FAILS))
for m in FAILS:
    print("  ", m)
print(t(), "done")
