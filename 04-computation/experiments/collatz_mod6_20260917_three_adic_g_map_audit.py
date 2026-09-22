#!/usr/bin/env python3
"""
collatz_mod6_20260917_three_adic_g_map_audit.py  (adversarial audit of lane three_adic_g_map)

Independent recomputation of the key numbers of
  04-computation/experiments/collatz_mod6_20260917_three_adic_g_map.py
with separately written code (no reduction tricks shared with the lane script), plus targeted
hunts for the issues found while reading the note:
  A1  k-table, transfer law J=1..7 (J=1 special), sharpness witnesses
  A2  the 2-block recoding: the note's iff "k_(n+1) in {0,2} iff X_n in {1,2,4,5}" is REFUTED;
      the true resolving statement is "k_(n+1) in {0,2} iff X_n in {1,4,7} iff X_(n-1) in {1,2,4,5}"
  A3  genuine Markov property at level J=3 (conditional law of X_2 given (X_0,X_1)), mu-invariance,
      Perron data, admissible k-word counts
  A4  d_J, g_J for J<=12 by unreduced integer iteration; Terras F(J) for J<=12
  A5  sigma vs sigma_res on m<=10^6 (pure Python), max sigma
  A6  exact moments for u in {1/2, 2, 3, 5, 7} at J<=4, u*, sharp factor
  A7  threshold enumeration i<=12: max m*(w), candidates, and the cycle boundary m = m*(w)
      (finds every positive G-cycle of length <=12)
  A8  census |b|<=49: cycle counts, content d = gcd(cycle, b), primitive (d=1) vs non-universal (d<|b|),
      fixed points, F3 counts
  A9  hostiles: 3^j+1, peak bound on m<=10^6, full sweep m<=10^7 (independent numpy), (1,3,2^(n-1),0)
      family: number of classes mod 3^(n+3) (two, not one), smallest members
  A10 reversed G-orbit is an E-path (m<=3000)
All checks are explicit raises (active under python -O).
"""
import sys, time
from fractions import Fraction as Fr
from math import gcd, log, sqrt
import numpy as np

T0 = time.time()


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def banner(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def kof(m):
    """minimal k>=0 with 2^k m in {4,7} mod 9 (independent implementation)."""
    for k in range(0, 7):
        if (pow(2, k, 9) * m) % 9 in (4, 7):
            return k
    raise RuntimeError("no k")


def G(m):
    k = kof(m)
    return ((2 ** k) * m - 1) // 3, k


def Gb(m, b):
    for j in range(0, 7):
        x = (pow(2, j, 9) * m) % 9
        if x == (b + 3) % 9 or x == (b + 6) % 9:
            return ((2 ** j) * m - b) // 3, j
    raise RuntimeError("no j")


def Tb(n, b):
    y = 3 * n + b
    e = 0
    while y % 2 == 0:
        y //= 2
        e += 1
    return y, e


def word(m, J):
    w = []
    x = m
    for _ in range(J):
        x, k = G(x)
        w.append(k)
    return w, x


# ============================================================================
banner("A1  k-table, transfer law, sharpness")
# ============================================================================
KT = {r: kof(r) for r in (1, 2, 4, 5, 7, 8)}
print("k-table:", KT)
check(KT == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "k-table")
for m in range(1, 2000):
    if m % 3 == 0:
        continue
    g, k = G(m)
    check(3 * g + 1 == 2 ** k * m and g % 3 != 0, "E-path identity / unit image at m=%d" % m)
# J=1: lifts of 1 and 2 mod 3 to mod 9
for a in (1, 2):
    ks = [G(a + 3 * t)[1] for t in range(3)]
    ims = sorted(G(a + 3 * t)[0] % 3 for t in range(3))
    print("  J=1 a=%d: k along lifts %s, images mod 3 %s" % (a, ks, ims))
    check(ims == [1, 1, 2], "J=1 images (1,1,2)")
check([G(1 + 3 * t)[1] for t in range(3)] == [2, 0, 0] and [G(2 + 3 * t)[1] for t in range(3)] == [1, 3, 1], "J=1 k pattern")
for J in range(2, 8):
    MJ = 3 ** J
    for a in range(MJ):
        if a % 3 == 0:
            continue
        ks = set()
        ims = []
        for t in range(3):
            g, k = G(a + t * MJ)
            ks.add(k)
            ims.append(g % MJ)
        check(len(ks) == 1, "same k J=%d a=%d" % (J, a))
        check(len(set(ims)) == 3 and len(set(x % (MJ // 3) for x in ims)) == 1, "bijective lifts J=%d a=%d" % (J, a))
        # G mod 3^J is a function of m mod 3^(J+1): compare with a+7*3^(J+1)
        for t in range(3):
            check(G(a + t * MJ)[0] % MJ == G(a + t * MJ + 7 * 3 * MJ)[0] % MJ, "function of m mod 3^(J+1)")
print("  transfer law verified J=2..7 (all unit classes)")
# sharpness witnesses quoted in the note
for J, (a, a2, i1, i2) in {1: (1, 7, 1, 2), 2: (1, 10, 1, 4), 3: (1, 28, 1, 10), 4: (1, 82, 1, 28)}.items():
    MJ = 3 ** J
    check(a2 % MJ == a % MJ and G(a)[0] % MJ == i1 and G(a2)[0] % MJ == i2, "sharpness witness J=%d" % J)
print("  sharpness witnesses (1,7->1,2), (1,10->1,4), (1,28->1,10), (1,82->1,28) confirmed")

# ============================================================================
banner("A2  the 2-block recoding: which iff resolves the pairs {4,7}, {2,8}?")
# ============================================================================
# residues X_(n-1) with k=0: {4,7}; with k=1: {2,8}.  Test both statements over all m mod 3^4.
bad_note = 0
bad_true = 0
wit = None
for m in range(1, 81):
    if m % 3 == 0:
        continue
    X0 = m % 9
    X1 = G(m)[0] % 9
    k2 = KT[X1]
    # note's statement: k_2 in {0,2} <=> X_1 in {1,2,4,5}
    if (k2 in (0, 2)) != (X1 in (1, 2, 4, 5)):
        bad_note += 1
        if wit is None:
            wit = (m, X0, X1, k2)
    # true statement: k_2 in {0,2} <=> X_1 in {1,4,7} <=> X_0 in {1,2,4,5}
    if (k2 in (0, 2)) != (X1 in (1, 4, 7)) or (X1 in (1, 4, 7)) != (X0 in (1, 2, 4, 5)):
        bad_true += 1
print("  note's clause 'k_(n+1) in {0,2} iff X_n in {1,2,4,5}': violations among units mod 81 = %d, first witness (m, X0, X1, k2) = %s" % (bad_note, wit))
print("  corrected clause 'k_(n+1) in {0,2} iff X_n in {1,4,7} iff X_(n-1) in {1,2,4,5}': violations = %d" % bad_true)
check(bad_note > 0, "the note's clause is refuted (e.g. X_n = 7 has k=0 but 7 not in {1,2,4,5})")
check(bad_true == 0, "corrected clause holds")
# the pair resolution itself: (k_n, k_(n+1)) determines X_(n-1)
for J in range(1, 7):
    M1 = 3 ** (J + 1)
    seen = {}
    for m in range(M1):
        if m % 3 == 0:
            continue
        w, _ = word(m, J + 1)
        seen.setdefault(tuple(w), set()).add(m)
    check(all(len(v) == 1 for v in seen.values()), "(J+1)-letter k-word pins m mod 3^(J+1), J=%d" % J)
    nJ = len(set(w[:J] for w in seen))
    check(nJ == 4 * 3 ** (J - 1), "J-letter k-words 4*3^(J-1)")
print("  (J+1)-letter word pins m mod 3^(J+1) and #J-letter words = 4*3^(J-1) for J<=6  [checked]")

# ============================================================================
banner("A3  Markov property at J=3, mu-invariance, Perron data")
# ============================================================================
# conditional law of X_2 = G^2(m) mod 27 given (X_0, X_1) under uniform m mod 3^5
J = 3
MJ = 27
hist = {}
for m in range(3 ** 5):
    if m % 3 == 0:
        continue
    x0 = m % MJ
    g1 = G(m)[0]
    x1 = g1 % MJ
    x2 = G(g1)[0] % MJ
    hist.setdefault((x0, x1), {}).setdefault(x2, 0)
    hist[(x0, x1)][x2] += 1
for (x0, x1), d in hist.items():
    tot = sum(d.values())
    check(len(d) == 3 and all(3 * c == tot for c in d.values()), "X_2 uniform on 3 classes given history")
    check(all(x2 % 9 == G(x1)[0] % 9 for x2 in d), "X_2 lies over G(X_1) mod 9")
print("  J=3: given (X_0,X_1) the residue X_2 is uniform on the three lifts of G(X_1) mod 9 (%d histories)  [checked]" % len(hist))
# mu-invariance by direct pushforward on classes mod 3^(J+1), J<=6
for J in range(1, 7):
    MJ, M1 = 3 ** J, 3 ** (J + 1)
    push = {}
    for a in range(M1):
        if a % 3 == 0:
            continue
        w = Fr(4, 3) if a % 3 == 1 else Fr(2, 3)
        b = G(a)[0] % MJ
        push[b] = push.get(b, 0) + w / (2 * 3 ** J)
    for b in range(MJ):
        if b % 3 == 0:
            continue
        w = Fr(4, 3) if b % 3 == 1 else Fr(2, 3)
        check(push[b] == w / (2 * 3 ** (J - 1)), "mu invariance J=%d b=%d" % (J, b))
print("  mu = w(m mod 3) Haar is G-invariant on cylinders of level J<=6  [checked]")
# preimage multiplicity: 4 classes over 1 mod 3, 2 over 2 mod 3
pre = {1: 0, 2: 0}
for a in (1, 2, 4, 5, 7, 8):
    pre[G(a)[0] % 3] += 1
check(pre == {1: 4, 2: 2}, "4-to-1 over 1 mod 3, 2-to-1 over 2 mod 3")
# Perron data
STATES = (1, 2, 4, 5, 7, 8)
A = [[1 if b % 3 == G(a)[0] % 3 else 0 for b in STATES] for a in STATES]
pi = [Fr(2, 9), Fr(1, 9), Fr(2, 9), Fr(1, 9), Fr(2, 9), Fr(1, 9)]
check(all(sum(r) == 3 for r in A), "row sums 3")
check(all(sum(pi[i] * A[i][j] for i in range(6)) == 3 * pi[j] for j in range(6)), "left Perron vector")
Ek = sum(pi[i] * KT[STATES[i]] for i in range(6))
check(Ek == 1, "E_pi[k]=1")
print("  A row sums 3, pi=(2,1,2,1,2,1)/9 left eigenvector, E_pi[k]=%s, uniform E[k]=%s, drift log(2/3)=%.6f" % (Ek, Fr(sum(KT.values()), 6), log(2 / 3)))
# primitivity at J<=5 by matrix powers
for J in range(1, 6):
    MJ = 3 ** J
    U = [a for a in range(MJ) if a % 3]
    idx = {a: i for i, a in enumerate(U)}
    P = np.zeros((len(U), len(U)), dtype=np.int64)
    for a in U:
        for t in range(3):
            P[idx[a], idx[G(a + t * MJ)[0] % MJ]] += 1
    Q = np.linalg.matrix_power(P, J)
    check(np.all(Q > 0), "P_J^J positive at J=%d" % J)
    if J >= 2:
        Q1 = np.linalg.matrix_power(P, J - 1)
        check(not np.all(Q1 > 0), "P_J^(J-1) has zeros at J=%d" % J)
print("  P_J^J > 0 and P_J^(J-1) not > 0 for J<=5: primitive with index exactly J  [checked]")

# ============================================================================
banner("A4  d_J, g_J (J<=12) by unreduced iteration; Terras F(J), J<=12")
# ============================================================================
t1 = time.time()
JMAX = 12
mod = 3 ** (JMAX + 1)
r = np.arange(mod, dtype=np.int64)
r = r[r % 3 != 0]
NU = len(r)
x = r.copy()
K = np.zeros(NU, dtype=np.int64)
done = np.zeros(NU, dtype=bool)
KT_ARR = np.zeros(9, dtype=np.int64)
for rr, kk in KT.items():
    KT_ARR[rr] = kk
dJ = {}
gJ = {}
expected_d = {1: "2/3", 2: "8/9", 3: "25/27", 4: "26/27", 5: "236/243", 6: "239/243", 7: "241/243", 8: "2173/2187",
              9: "19609/19683", 10: "58868/59049", 11: "176809/177147", 12: "530885/531441"}
expected_g = {1: "2/3", 2: "5/6", 3: "43/54", 4: "73/81", 5: "214/243", 6: "686/729", 7: "2126/2187", 8: "12647/13122",
              9: "19339/19683", 10: "6413/6561", 11: "58396/59049", 12: "1057277/1062882"}
for i in range(1, JMAX + 1):
    k = KT_ARR[x % 9]
    x = (x * (2 ** k).astype(np.int64) - 1) // 3   # no modular reduction: true integer values of G^i(a), a < 3^13
    K += k
    now = (2.0 ** K) < (3.0 ** i)      # float comparison is exact here (K <= 2i+1 <= 25)
    done |= now
    dJ[i] = Fr(int(done.sum()), NU)
    gJ[i] = Fr(int(now.sum()), NU)
    check(str(dJ[i]) == expected_d[i], "d_%d = %s (note says %s)" % (i, dJ[i], expected_d[i]))
    check(str(gJ[i]) == expected_g[i], "g_%d = %s (note says %s)" % (i, gJ[i], expected_g[i]))
    check(1 - dJ[i] <= Fr(7, 9) ** (i - 1), "Markov bound")
# exactness of the float comparison: check 2^K vs 3^i on the boundary values
for i in range(1, 30):
    for KK in range(0, 2 * i + 3):
        check((2.0 ** KK < 3.0 ** i) == (2 ** KK < 3 ** i), "float boundary")
print("  d_J and g_J for J<=12 agree with the note's table; max int64 value used %d  (%.1fs)" % (int(x.max()), time.time() - t1))
del r, x, K, done
ratios = ["%.3f" % float((1 - dJ[i + 1]) / (1 - dJ[i])) for i in range(1, 12)]
print("  ratios (1-d_(J+1))/(1-d_J): %s; geometric mean J=8..12: %.4f" % (ratios, float((1 - dJ[12]) / (1 - dJ[8])) ** 0.25))
# Terras densities mod 2^J, J<=12, per-residue simulation
terras_expected = {1: "1/2", 2: "3/4", 3: "3/4", 4: "13/16", 5: "7/8", 6: "7/8", 7: "115/128", 8: "237/256", 9: "237/256",
                   10: "15/16", 11: "15/16", 12: "1935/2048"}
JC = 12
cnt = [0] * (JC + 1)
for n0 in range(2 ** JC):
    n = n0
    a = 0
    stopped = None
    for i in range(1, JC + 1):
        if n % 2:
            n = (3 * n + 1) // 2
            a += 1
        else:
            n //= 2
        if 3 ** a < 2 ** i:
            stopped = i
            break
    if stopped is not None:
        for i in range(stopped, JC + 1):
            cnt[i] += 1
for i in range(1, JC + 1):
    F = Fr(cnt[i], 2 ** JC)
    check(str(F) == terras_expected[i], "Terras F(%d)=%s" % (i, F))
print("  Terras F(J), J<=12, agree with the note (F(12)=%s=%.4f)" % (Fr(cnt[12], 2 ** JC), cnt[12] / 2 ** JC))

# ============================================================================
banner("A5  sigma vs sigma_res on m<=10^6 (pure Python)")
# ============================================================================
t1 = time.time()
KMAX = {i: max(k for k in range(0, 3 * i + 2) if 2 ** k < 3 ** i) for i in range(1, 200)}
mism = []
max_sig = (0, 0)
for m in range(2, 10 ** 6 + 1):
    if m % 3 == 0:
        continue
    x = m
    Kc = 0
    sig = 0
    sres = 0
    i = 0
    while not (sig and sres):
        i += 1
        r9 = x % 9
        k = KT[r9]
        x = (x * (1 << k) - 1) // 3
        Kc += k
        if not sres and Kc <= KMAX[i]:
            sres = i
        if not sig and x < m:
            sig = i
    if sig != sres:
        mism.append((m, sig, sres))
    if sig > max_sig[0]:
        max_sig = (sig, m)
print("  m<=10^6: mismatches {sigma != sigma_res} = %d %s; max sigma = %d at m=%d  (%.1fs)" % (len(mism), mism[:10], max_sig[0], max_sig[1], time.time() - t1))
check(len(mism) == 0 and max_sig == (31, 128669), "no mismatch to 10^6; max sigma 31 at 128669")

# ============================================================================
banner("A6  exact moments for every tilt (incl. u=1/2), u*, sharp factor")
# ============================================================================
for u in (Fr(1, 2), Fr(2), Fr(3), Fr(5), Fr(7)):
    for J in (1, 2, 3, 4):
        M1 = 3 ** (J + 1)
        tot = Fr(0)
        n = 0
        for a in range(M1):
            if a % 3 == 0:
                continue
            w, _ = word(a, J)
            tot += u ** sum(w)
            n += 1
        lhs = tot / n
        rhs = (u + 1) * (u * u + 2) / 6 * ((u * u + u + 1) / 3) ** (J - 1)
        check(lhs == rhs, "moment identity u=%s J=%d: %s vs %s" % (u, J, lhs, rhs))
    print("  u=%s: E[u^K_J] = (u+1)(u^2+2)/6 * ((u^2+u+1)/3)^(J-1) verified J<=4" % u)
c = log(3) / log(2)
# Lambda'(theta) = u(2u+1)/(u^2+u+1); solve = c by bisection on u
lo, hi = 1.0, 10.0
for _ in range(200):
    mid = (lo + hi) / 2
    if mid * (2 * mid + 1) / (mid * mid + mid + 1) < c:
        lo = mid
    else:
        hi = mid
u_star = (lo + hi) / 2
rho = (u_star ** 2 + u_star + 1) / 3
fac = rho / u_star ** c
print("  c=log_2 3=%.6f; u*=%.6f; rho(u*)=%.6f; exp(-I(c))=%.6f; 7/9=%.6f; Lambda'(log 2)=%.4f" % (c, u_star, rho, fac, 7 / 9, 10 / 7))
check(abs(u_star - 2.782079) < 2e-6 and abs(rho - 3.840680) < 2e-6 and abs(fac - 0.758751) < 2e-6, "u*, rho, factor")
check(fac < 7 / 9, "sharp factor below 7/9")

# ============================================================================
banner("A7  threshold enumeration i<=12; cycle boundary m = m*(w)")
# ============================================================================
t1 = time.time()
note_table = {1: ("1", 0), 2: ("11/7", 0), 3: ("53/5", 1), 4: ("239/47", 1), 5: ("973/13", 3), 6: ("827/59", 2), 7: ("17269/1909", 2),
              8: ("63071/1631", 4), 9: ("51769/2617", 4), 10: ("940375/6487", 8), 11: ("3820549/84997", 8), 12: ("15459343/517135", 8)}
cands_all = []
cycles_found = []
for i in range(1, 13):
    M1 = 3 ** (i + 1)
    maxr = Fr(0)
    n_grow = 0
    n_cand = 0
    n_mis = 0
    for a in range(1, M1):
        if a % 3 == 0:
            continue
        x = a
        Kt = 0
        B = 0
        for s in range(i):
            x, k = G(x)
            B = B * 2 ** k + 3 ** s
            Kt += k
        check((2 ** Kt) * a - (3 ** i) * x == B, "carry B_i")
        if 2 ** Kt <= 3 ** i:
            continue
        n_grow += 1
        mstar = Fr(B, 2 ** Kt - 3 ** i)
        maxr = max(maxr, mstar)
        # cycle boundary: m = m*(w) integer in the class a
        if mstar.denominator == 1 and mstar % M1 == a and mstar >= 1:
            cycles_found.append((i, int(mstar)))
            xx = int(mstar)
            for _ in range(i):
                xx = G(xx)[0]
            check(xx == int(mstar), "boundary point is a cycle")
        mm = a
        while mm < mstar:
            if mm >= 2:
                n_cand += 1
                cands_all.append((i, mm))
                y = mm
                early = False
                for l in range(1, i):
                    y = G(y)[0]
                    if y < mm:
                        early = True
                        break
                if not early:
                    n_mis += 1
            mm += M1
    check(str(maxr) == note_table[i][0] and n_cand == note_table[i][1] and n_mis == 0, "level %d: max m* %s cand %d mis %d" % (i, maxr, n_cand, n_mis))
    print("  i=%2d: growth classes %5d, max m* = %s = %.4f, candidates %d, true mismatches %d" % (i, n_grow, maxr, float(maxr), n_cand, n_mis))
print("  candidates (i, m): %s" % cands_all)
print("  integer boundary points m = m*(w) in their own class (= positive G-cycles of length i): %s" % cycles_found)
check(all(m == 1 for _, m in cycles_found) and len(cycles_found) == 12, "the only positive G-cycle of length <=12 is {1} (word 2^i at every level)")
# the note's stopping-boundary conjecture 'm*(w) < least positive member of the class' is refuted:
w3, _ = word(2, 3)
check(w3 == [1, 2, 2] and Fr(37, 5) > 2, "witness: word (1,2,2) at level 3 has class 2 mod 81, m* = 37/5 > 2")
print("  REFUTED (stopping-boundary conjecture): word (1,2,2), class 2 mod 81, m* = 37/5 = 7.4 > least member 2 (m=2 descends at step 1).  (%.1fs)" % (time.time() - t1))

# ============================================================================
banner("A8  census |b|<=49: counts, content d = gcd(cycle,b), primitive vs non-universal, fixed points")
# ============================================================================
t1 = time.time()


def cycles_of(step, starts, ESC=10 ** 18):
    label = {}
    cyc = []
    for s in starts:
        if s in label:
            continue
        path = []
        seen = {}
        x = s
        while x not in label and x not in seen:
            if abs(x) > ESC:
                raise RuntimeError("escape")
            seen[x] = len(path)
            path.append(x)
            x = step(x)
        if x in seen:
            cyc.append(path[seen[x]:])
            lab = len(cyc) - 1
        else:
            lab = label[x]
        for y in path:
            label[y] = lab
    return cyc


BS = [b for b in range(-49, 50) if b % 2 and b % 3]
F3 = {1: 4, 5: 9, 7: 5, 11: 7, 13: 13}
resG = {}
resT = {}
for b in BS:
    startsG = sorted((m for m in range(-10 ** 5, 10 ** 5 + 1) if m and m % 3), key=abs)
    resG[b] = cycles_of(lambda m, b=b: Gb(m, b)[0], startsG)
    startsT = sorted((n for n in range(-2 * 10 ** 5, 2 * 10 ** 5 + 1) if n % 2), key=abs)
    resT[b] = cycles_of(lambda n, b=b: Tb(n, b)[0], startsT)
print("  census done (%.1fs)" % (time.time() - t1))
note_counts = {1: (3, 4), 5: (5, 9), 7: (4, 5), 11: (5, 7), 13: (7, 13), 17: (6, 8), 19: (6, 6), 23: (8, 15), 25: (7, 12), 29: (4, 9),
               31: (4, 6), 35: (7, 12), 37: (4, 7), 41: (6, 5), 43: (5, 5), 47: (6, 11), 49: (8, 7)}
print("  b | #G_b | #T_b | G_b non-universal (d<|b|) | G_b primitive (d=1) | G_b contents d of non-universal cycles | T_b non-universal | T_b primitive")
for b in BS:
    cG, cT = resG[b], resT[b]
    if b > 0:
        check((len(cG), len(cT)) == note_counts[b], "counts b=%d: %d %d" % (b, len(cG), len(cT)))
    if abs(b) in F3:
        check(len(cT) == F3[abs(b)], "F3 count")
    # negation symmetry
    if b > 0:
        check(set(frozenset(-x for x in c) for c in cG) == set(frozenset(c) for c in resG[-b]), "G negation")
        check(set(frozenset(-x for x in c) for c in cT) == set(frozenset(c) for c in resT[-b]), "T negation")
    # gate identity on every G_b cycle and content
    dG = []
    for cyc in cG:
        L = len(cyc)
        js = [Gb(cyc[i], b)[1] for i in range(L)]
        Jt = sum(js)
        Bp = 0
        Ji = 0
        for i in range(1, L + 1):
            Ji += js[i - 1]
            Bp += 3 ** (i - 1) * 2 ** (Jt - Ji)
        check(cyc[0] * (2 ** Jt - 3 ** L) == b * Bp, "G gate b=%d" % b)
        d = abs(b)
        for x in cyc:
            d = gcd(d, x)
        check(all(gcd(abs(b), x) == d for x in cyc), "content constant along the cycle")
        dG.append(d)
    dT = []
    for cyc in cT:
        d = abs(b)
        for x in cyc:
            d = gcd(d, x)
        dT.append(d)
    # universal cycles: content |b| exactly three for G (b{1}, b{-1}, b{-4,-11}) and four for T
    check(sum(1 for d in dG if d == abs(b)) == 3, "three universal G_b cycles")
    check(sum(1 for d in dT if d == abs(b)) == 4, "four universal T_b cycles")
    # common cycles = the two fixed points
    common = set(frozenset(c) for c in cG) & set(frozenset(c) for c in cT)
    check(common == {frozenset([b]), frozenset([-b])}, "common cycles are {b},{-b}")
    if b > 0:
        nonu = [(min(c, key=abs), d) for c, d in zip(cG, dG) if d < b]
        print("  %2d | %d | %2d | %d | %d | %s | %d | %d" % (b, len(cG), len(cT), sum(1 for d in dG if d < b), sum(1 for d in dG if d == 1),
                                                          sorted(nonu, key=lambda t: abs(t[0])), sum(1 for d in dT if d < b), sum(1 for d in dT if d == 1)))
# the composite-b discrepancy: 25, 35, 49 have non-universal cycles of content > 1
for b, exp_prim, exp_nonu in ((25, 2, 4), (35, 1, 4), (49, 4, 5)):
    dG = []
    for cyc in resG[b]:
        d = b
        for x in cyc:
            d = gcd(d, x)
        dG.append(d)
    check(sum(1 for d in dG if d == 1) == exp_prim and sum(1 for d in dG if d < b) == exp_nonu, "b=%d primitive %d vs non-universal %d" % (b, exp_prim, exp_nonu))
print("  b=25: non-universal 4 but primitive (d=1) 2 (cycles -20, 80 are 5 x the G_5 cycles -4, 16); b=35: 4 vs 1; b=49: 5 vs 4  [checked]")
# fixed points of G_b for all |b|<=49: exactly {b, -b}
for b in BS:
    fps = set()
    for j in range(0, 4):
        num = b
        den = 2 ** j - 3
        if den != 0 and num % den == 0:
            m = num // den
            if m % 3 and Gb(m, b) == (m, j):
                fps.add(m)
    check(fps == {b, -b}, "fixed points of G_b are {b,-b}, b=%d: %s" % (b, fps))
print("  fixed points of G_b: exactly {b, -b} for every odd |b|<=49, 3 not | b  [checked]")
# scaling lemma
for b in BS:
    for m in range(-300, 301):
        if m % 3 == 0:
            continue
        g1, j1 = Gb(m, 1)
        gb, jb = Gb(b * m, b)
        check(gb == b * g1 and jb == j1, "scaling lemma")
print("  scaling lemma G_b(bm) = b G_1(m), same j, verified |b|<=49, |m|<=300")
# sample b=13 positive cycles quoted in the note
for cyc, w in (([337, 445, 589, 781, 256], [2, 2, 2, 0, 2]), ([209, 553, 733, 973, 320], [3, 2, 2, 0, 1]), ([236, 625, 829, 272, 721], [3, 2, 0, 3, 0])):
    ww = [Gb(cyc[i], 13)[1] for i in range(5)]
    check(ww == w and all(Gb(cyc[i], 13)[0] == cyc[(i + 1) % 5] for i in range(5)), "b=13 cycle words")
print("  b=13 positive cycle words (2,2,2,0,2), (3,2,2,0,1), (3,2,0,3,0) confirmed; 2^8-3^5 = %d" % (2 ** 8 - 3 ** 5))

# ============================================================================
banner("A9  hostiles: 3^j+1, peak bound, full sweep m<=10^7, the (1,3,2^(n-1),0) family")
# ============================================================================
for j in range(1, 16):
    m = 3 ** j + 1
    x = m
    for i in range(j):
        if i <= j - 1:
            check(x == 4 ** i * 3 ** (j - i) + 1, "3^j+1 orbit value")
        x = G(x)[0]
    check(x == 4 ** (j - 1), "G^j(3^j+1) = 4^(j-1), j=%d" % j)
check(3 * 4 ** 14 + 1 == 805306369 and abs(805306369 / 14348908 - 56.123) < 1e-3, "j=15 numbers")
print("  3^j+1: G^i = 4^i 3^(j-i)+1 (i<j), G^j = 4^(j-1) for j<=15; j=15 ratio %.3f" % (805306369 / 14348908))
# peak bound along orbits, m<=10^6 (pure Python), and full sweep m<=10^7 (independent numpy, single pass in 2 halves)
t1 = time.time()
viol = 0
for m in range(1, 10 ** 6 + 1, 1):
    if m % 3 == 0:
        continue
    x = m
    Kc = 0
    n = 0
    f3 = 1 if m % 9 == 5 else 0
    while x != 1:
        x, k = G(x)
        Kc += k
        n += 1
        if Kc > 2 * n + f3:
            viol += 1
            break
check(viol == 0, "K_n <= 2n + [k_1=3] on m<=10^6")
print("  K_n <= 2n + [k_1=3] along every orbit m<=10^6 (pure Python)  (%.1fs)" % (time.time() - t1))
t1 = time.time()
best_steps = (0, 0)
best_ratio = (0.0, 0, 0)
total = 0
for lo, hi in ((1, 5 * 10 ** 6), (5 * 10 ** 6 + 1, 10 ** 7)):
    m0 = np.arange(lo, hi + 1, dtype=np.int64)
    m0 = m0[m0 % 3 != 0]
    total += len(m0)
    cur = m0.copy()
    pk = m0.copy()
    st = np.zeros(len(m0), dtype=np.int64)
    alive = np.ones(len(m0), dtype=bool)
    alive[m0 == 1] = False
    n = 0
    while alive.any():
        n += 1
        ii = np.flatnonzero(alive)
        v = cur[ii]
        k = KT_ARR[v % 9]
        v = (v * (2 ** k).astype(np.int64) - 1) // 3
        cur[ii] = v
        pk[ii] = np.maximum(pk[ii], v)
        st[ii] = n
        alive[ii[v == 1]] = False
        check(n < 500, "termination")
    i1 = int(st.argmax())
    if int(st[i1]) > best_steps[0]:
        best_steps = (int(st[i1]), int(m0[i1]))
    rat = pk / m0
    i2 = int(rat.argmax())
    if float(rat[i2]) > best_ratio[0]:
        best_ratio = (float(rat[i2]), int(m0[i2]), int(pk[i2]))
    # ties for the max step count?
    ties = m0[st == st[i1]]
    print("  range [%d,%d]: max steps %d at m=%s (ties: %s), max peak/m %.4f at m=%d" % (lo, hi, int(st[i1]), int(m0[i1]), ties.tolist(), float(rat[i2]), int(m0[i2])))
    del m0, cur, pk, st, alive, rat
check(total == 6666667 and best_steps == (93, 8751065) and best_ratio[1] == 4847486 and best_ratio[2] == 644874241, "sweep extremes")
print("  all %d non-multiples of 3 in [1,10^7] reach 1; max steps %d at %d; max peak/m %.3f at %d (peak %d)  (%.1fs)"
      % (total, best_steps[0], best_steps[1], best_ratio[0], best_ratio[1], best_ratio[2], time.time() - t1))
# the (8,5,1^(n-1),{4,7}) family = word 1,3,2^(n-1),0: number of classes mod 3^(n+3), smallest member
note_smallest = [8, 143, 62, 1277, 548, 11483, 4922, 103337, 44288, 930023, 398582]
for n in range(1, 8):
    M = 3 ** (n + 3)
    target = [1, 3] + [2] * (n - 1) + [0]
    classes = [a for a in range(M) if a % 3 and word(a, n + 2)[0] == target]
    check(len(classes) == 2, "word 1,3,2^(n-1),0 is realized by exactly TWO classes mod 3^(n+3), n=%d: %s" % (n, classes))
    check(all(a % 27 == 8 for a in classes), "both classes are 8 mod 27")
    check(len(set(a % (3 ** (n + 2)) for a in classes)) == 1, "one class mod 3^(n+2) (the prefix 1,3,2^(n-1))")
    check(min(classes) == note_smallest[n - 1], "smallest member n=%d" % n)
    print("  n=%d: word %s realized by classes %s mod 3^%d (two lifts of one class mod 3^%d); smallest member %d" % (n, target, classes, n + 3, n + 2, min(classes)))
# residue count check: residues 8,5 then (n-1) residues 1 then {4,7}
w11, _ = word(398582, 13)
check(w11 == [1, 3] + [2] * 10 + [0], "n=11 word")
res = []
x = 398582
for _ in range(13):
    res.append(x % 9)
    x = G(x)[0]
check(res == [8, 5] + [1] * 10 + [4] or res == [8, 5] + [1] * 10 + [7], "n=11 residues 8,5,1^10,{4,7}: %s" % res)
print("  n=11: residues X_0..X_12 = %s (ten residues 1 = n-1, not n)" % res)
# the max-steps tie: 9454814 also needs 93 steps
for m_t in (8751065, 9454814):
    x = m_t
    n = 0
    while x != 1:
        x = G(x)[0]
        n += 1
    check(n == 93 and m_t % 3 != 0, "93 steps at m=%d" % m_t)
print("  steps-to-1 = 93 is attained at BOTH m=8751065 (first) and m=9454814 (tie) in [1,10^7]")
# the worst m: 3-adic reason
m = 4847486
orb = [m]
x = m
w = []
while x != 1:
    x, k = G(x)
    orb.append(x)
    w.append(k)


def v3(x):
    e = 0
    while x % 3 == 0:
        x //= 3
        e += 1
    return e


check(len(w) == 76 and max(orb) == 644874241 and orb.index(max(orb)) == 17 and sum(w[:17]) == 34, "worst orbit data")
check(v3(orb[1] - 1) == 6 and v3(orb[8] - 1) == 10 and w[:18] == [3] + [2] * 5 + [0, 3] + [2] * 9 + [0], "3-adic reason")
print("  m=4847486: 76 letters, peak 644874241 at step 17, K_17=34, v_3(G(m)-1)=6, v_3(G^8(m)-1)=10, 2^34/3^17=%.2f" % (2 ** 34 / 3 ** 17))

# ============================================================================
banner("A10  reversed G-orbit is an E-path")
# ============================================================================
for m in range(1, 3001):
    if m % 3 == 0:
        continue
    x = m
    while x != 1:
        y, k = G(x)
        # E-path from y to x: y -> 3y+1 = 2^k x -> halvings -> x
        z = 3 * y + 1
        check(z == 2 ** k * x, "E-step")
        x = y
print("  for every m<=3000 (3 not | m) the reversed G-orbit is an E-path 1 -> m  [checked]")

print("\ntotal time %.1fs" % (time.time() - T0))
print("ALL AUDIT CHECKS PASSED")
