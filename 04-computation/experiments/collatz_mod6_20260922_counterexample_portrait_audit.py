#!/usr/bin/env python3
"""Adversarial audit of lane counterexample_portrait (session collatz-mod6, 2026-09-22).

Independent recomputation of every key number of
collatz_mod6_20260922_counterexample_portrait.{py,out,md}: different code paths
(brute-force residue counts instead of the DP, Burnside instead of orbit
enumeration, a sieve written from scratch, mpmath continued fraction by the
standard algorithm), plus the checks that the audit found necessary:
the Hardy-Wright 171 spacing inequality verified numerically on every
convergent used, the exact L_A boundary (L_A passes, L_A+1 fails), and the
count of |2^K-3^L| <= 100 clocks in the boxes that the note compares with the
pillai lane.  Every check raises on failure.  Run time about 60 s, < 300 MB.
"""
import math
import sys
import time
from fractions import Fraction
from itertools import combinations

import mpmath

mpmath.mp.dps = 400
ALPHA = mpmath.log(3) / mpmath.log(2)
LN2 = mpmath.log(2)
T0 = time.time()


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def hdr(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


def v2(x):
    return (x & -x).bit_length() - 1


def Tb(n, b):
    m = 3 * n + b
    return m >> v2(m)


# ---------------------------------------------------------------- A1
hdr("A1. conjugacy T_+(-n) = -T_-(n), odd n in [-20001, 20001]")
bad = sum(1 for n in range(-20001, 20002, 2) if Tb(-n, +1) != -Tb(n, -1))
print(f"violations: {bad}")
check(bad == 0, "A1")

# ---------------------------------------------------------------- A2
hdr("A2. gate data of the four known cycles, recomputed from the words")
CYCLES = [(+1, [1]), (-1, [1]), (-1, [5, 7]), (-1, [17, 25, 37, 55, 41, 61, 91])]


def word_of_cycle(cyc, b):
    w = []
    for i, x in enumerate(cyc):
        m = 3 * x + b
        k = v2(m)
        check(m >> k == cyc[(i + 1) % len(cyc)], "cycle closes")
        w.append(k)
    return w


def gate(w):
    L = len(w)
    K = sum(w)
    Ks = [0]
    for k in w:
        Ks.append(Ks[-1] + k)
    B = sum(3 ** (L - 1 - i) * 2 ** Ks[i] for i in range(L))
    Delta = 2 ** K - 3 ** L
    q = abs(Delta) // math.gcd(B, abs(Delta))
    return K, L, B, Delta, q


for b, cyc in CYCLES:
    w = word_of_cycle(cyc, b)
    K, L, B, Delta, q = gate(w)
    n0 = Fraction(b * B, Delta)
    print(f"b={b:+d} cycle {cyc}: word {w} K/L={K}/{L} B={B} Delta={Delta} q={q} n0={n0}")
    check(q == 1 and n0 == cyc[0], "A2 gate")
    check((Delta > 0) == (b > 0), "A2 sign law: positive cycle has sign(Delta) = sign(b)")

# ---------------------------------------------------------------- A3
hdr("A3. continued fraction of log_2 3 (standard algorithm), convergents, 11/7 placement")
x = ALPHA
pq = []
a_list = []
p0, q0, p1, q1 = 0, 1, 1, 0
for i in range(60):
    a = int(mpmath.floor(x))
    a_list.append(a)
    p0, p1 = p1, a * p1 + p0
    q0, q1 = q1, a * q1 + q0
    pq.append((p1, q1))
    x = 1 / (x - a)
print("a_0..a_29:", a_list[:30])
check(a_list[:30] == [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8], "A3 pq")
for n in range(1, 60):
    check(abs(pq[n - 1][0] * pq[n][1] - pq[n][0] * pq[n - 1][1]) == 1, "A3 determinant")
signs = []
for n in range(25):
    p, q = pq[n]
    if q <= 80000:
        D = 2 ** p - 3 ** q
        s = 1 if D > 0 else -1
        tag = str(D) if abs(D) < 10 ** 6 else f"{'+' if D > 0 else '-'}({len(str(abs(D)))} digits)"
    else:
        s = 1 if p - q * ALPHA > 0 else -1
        tag = "by log"
    signs.append(s)
    check(s == (1 if n % 2 else -1), "A3 alternation")
    if n <= 11:
        print(f"  n={n:2d} {p}/{q} Delta {tag}")
print("alternation even->Delta<0, odd->Delta>0 holds for n<=24: True")
print(f"q_21={pq[21][1]} q_23={pq[23][1]} q_25={pq[25][1]} q_27={pq[27][1]}")
check(pq[23] == (217976794617, 137528045312) and pq[21][1] == 6586818670, "A3 q_21 q_23")
print(f"q_23 = {pq[23][1]} = {pq[23][1] / 1e11:.5f}e11  (Hercher's 'next bound' K >= 1.375e11 is this clock)")
# 11/7 as mediant
p2, q2 = pq[2]
p3, q3 = pq[3]
check((p2 + p3, q2 + q3) == (11, 7) and a_list[4] == 2, "A3 11/7 mediant, a_4=2")
print(f"11/7 = ({p2}+{p3})/({q2}+{q3}), a_4 = {a_list[4]} -> the only intermediate between 3/2 and 19/12")
for K, L in [(1, 1), (2, 1), (3, 2), (11, 7)]:
    D = abs(2 ** K - 3 ** L)
    print(f"  clock {K}/{L}: |Delta|={D} tierA cap {3 ** L / (4 * L):.3f} {'in' if D <= 3 ** L / (4 * L) else 'OUT'}"
          f" tierB cap {3 ** L / (2 * L):.3f} {'in' if D <= 3 ** L / (2 * L) else 'OUT'}")
check(139 > 2187 / 28 and 139 <= 2187 / 14, "A3 11/7 tier B")

# ---------------------------------------------------------------- A4
hdr("A4. Eliahou product identity and N_min(L) (tier A thresholds)")
for b, cyc in CYCLES:
    w = word_of_cycle(cyc, b)
    K, L, B, Delta, q = gate(w)
    prod = Fraction(1)
    for n in cyc:
        prod *= 1 + Fraction(b, 3 * n)
    check(prod == Fraction(2 ** K, 3 ** L), "A4 identity")
    N = min(cyc)
    rel = abs(Delta) / 3 ** L
    bound = (1 + 1 / (3 * N)) ** L - 1 if b > 0 else 1 - (1 - 1 / (3 * N)) ** L
    print(f"b={b:+d} {cyc}: |Delta|/3^L={rel:.4f} bound={bound:.4f}")
    check(rel <= bound + 1e-12, "A4 bound")


def nmin(L, sign):
    N = 1
    while True:
        val = (1 + 1 / (3 * N)) ** L - 1 if sign > 0 else 1 - (1 - 1 / (3 * N)) ** L
        if val <= 1 / (4 * L):
            return N
        N += 1


table = {}
for L in (2, 5, 7, 12, 41, 53):
    table[L] = (nmin(L, +1), nmin(L, -1))
    print(f"  L={L}: N_min plus={table[L][0]} minus={table[L][1]} 4L^2/3={4 * L * L / 3:.2f}")
check(table == {2: (6, 6), 5: (34, 33), 7: (67, 65), 12: (194, 191), 41: (2248, 2235), 53: (3754, 3737)}, "A4 N_min")

# ---------------------------------------------------------------- A5
hdr("A5. minus-sheet census to 10^7 (independent implementation)")
X = 10 ** 7
cyc_of = {}
for i, (b, cyc) in enumerate(CYCLES[1:], start=1):
    for c in cyc:
        cyc_of[c] = i
tag = bytearray(X // 2 + 1)
cnt = [0, 0, 0, 0]
best = (0, 0)
t = time.time()
for n in range(1, X + 1, 2):
    if n in cyc_of:
        r = cyc_of[n]
        steps = 0
    else:
        m = n
        steps = 0
        while True:
            y = 3 * m - 1
            m = y >> v2(y)
            steps += 1
            if m < n:
                r = tag[m >> 1]
                break
            if m in cyc_of:
                r = cyc_of[m]
                break
    check(r != 0, "A5 unclassified")
    tag[n >> 1] = r
    cnt[r] += 1
    if steps > best[0]:
        best = (steps, n)
tot = sum(cnt)
print(f"odd starts {tot}; basins {cnt[1:]} fractions {[round(c / tot, 6) for c in cnt[1:]]}; census {time.time() - t:.1f}s (timing)")
print(f"longest first-drop-or-cycle-hit segment: {best[0]} steps at n={best[1]}")
check(cnt[1:] == [1636054, 1623149, 1740797] and best == (172, 5960769), "A5 census")
# sample of full orbits (no memoisation) as a second control
import random
random.seed(1)
for _ in range(2000):
    n = random.randrange(1, X) | 1
    m = n
    for _ in range(10000):
        if m in cyc_of:
            break
        y = 3 * m - 1
        m = y >> v2(y)
    check(m in cyc_of, "A5 sample orbit reaches a cycle")
print("2000 random odd starts < 10^7 iterated to a cycle without memoisation: all reach a known cycle")

# ---------------------------------------------------------------- A6
hdr("A6. L_A, spacing thresholds, HW171 numeric check, admissible convergents")


def tierA_ok(L, N):
    # exp(L/(3N)) - 1 <= 1/(4L)  (sufficient on both sheets)
    return mpmath.expm1(mpmath.mpf(L) / (3 * N)) <= mpmath.mpf(1) / (4 * L)


for label, N, sheet, LA_claim, thr_claim, first_claim in [
        ("plus", 1 << 68, +1, 14878203146, 6.137428e20, 23),
        ("minus", 10 ** 7 + 1, -1, 2738, 1.946585e7, 10)]:
    check(tierA_ok(LA_claim, N) and not tierA_ok(LA_claim + 1, N), f"A6 L_A boundary {label}")
    if sheet > 0:
        factor = mpmath.exp(-mpmath.mpf(LA_claim) / (3 * N))
    else:
        y0 = -mpmath.log(mpmath.mpf(7) / 8)
        factor = (1 - mpmath.exp(-y0)) / y0
    thr = 3 * N * LN2 * factor
    print(f"{label}: L_A={LA_claim} passes, L_A+1={LA_claim + 1} fails; factor={float(factor):.6f}; threshold={float(thr):.6e}")
    check(abs(float(thr) / thr_claim - 1) < 1e-6, f"A6 threshold {label}")
    adm = []
    for n in range(0, 40):
        p, q = pq[n]
        qn = pq[n + 1][1]
        if (n % 2 == 1) != (sheet > 0):
            continue
        # HW 171 numeric check: |alpha - p/q| > 1/(q(q+q_next))
        gap = abs(ALPHA - mpmath.mpf(p) / q)
        check(gap > 1 / (mpmath.mpf(q) * (q + qn)) and gap < 1 / (mpmath.mpf(q) * qn), f"A6 HW171 n={n}")
        if q * (q + qn) >= thr:
            adm.append((n, p, q))
    print(f"  HW171 spacing 1/(q_n(q_n+q_(n+1))) < |alpha-p_n/q_n| < 1/(q_n q_(n+1)) verified for n<=39 on this side")
    print(f"  first admissible convergent: n={adm[0][0]} {adm[0][1]}/{adm[0][2]}; survivors {[a[0] for a in adm[:3]]}")
    check(adm[0][0] == first_claim, f"A6 first admissible {label}")
    check(all(q > LA_claim for (_, _, q) in adm), f"A6 no admissible with q<=L_A {label}")
    # the convergents with q_n <= L_A on the correct side all fail the threshold:
    small = [(n, pq[n][1]) for n in range(40) if ((n % 2 == 1) == (sheet > 0)) and pq[n][1] <= LA_claim]
    print(f"  convergents on this side with q_n <= L_A: {small}; all fail the threshold: "
          f"{all(q * (q + pq[n + 1][1]) < thr for n, q in small)}")
    check(all(q * (q + pq[n + 1][1]) < thr for n, q in small), "A6 all small fail")
print("NOTE (audit): tier A, hence 'clock is a convergent', is only available for L <= L_A; a cycle with")
print("L > L_A may sit on a non-convergent clock.  The bound L > L_A is what the chain proves; the")
print("convergent survivors describe only the case in which the clock happens to be a convergent.")
print(f"n=21 spacing {pq[21][1] * (pq[21][1] + pq[22][1]):.3e} vs plus threshold 6.137e20: excluded")

# ---------------------------------------------------------------- A7
hdr("A7. full gate census L <= 10, L <= K <= 2L (independent enumeration)")


def comps(K, L):
    # compositions of K into L positive parts via cut positions
    for cuts in combinations(range(1, K), L - 1):
        prev = 0
        w = []
        for c in cuts:
            w.append(c - prev)
            prev = c
        w.append(K - prev)
        yield tuple(w)


nwords = 0
found = {+1: set(), -1: set()}
for L in range(1, 11):
    for K in range(L, 2 * L + 1):
        for w in comps(K, L):
            nwords += 1
            Kt, Lt, B, Delta, q = gate(w)
            if q == 1:
                for b in (+1, -1):
                    n0 = b * B // Delta
                    orbit = [n0]
                    for _ in range(L - 1):
                        orbit.append(Tb(orbit[-1], b))
                    check(Tb(orbit[-1], b) == n0, "A7 closes")
                    found[b].add(tuple(sorted(set(orbit))))
print(f"words: {nwords} (sum C(2L,L) = {sum(math.comb(2 * L, L) for L in range(1, 11))})")
check(nwords == 250952, "A7 count")
for b in (+1, -1):
    print(f"  b={b:+d}: cycles (as sets) {sorted(found[b])}")
check(found[-1] == {(1,), (5, 7), (17, 25, 37, 41, 55, 61, 91), (-1,)}, "A7 minus")
check(found[+1] == {(1,), (-1,), (-7, -5), (-91, -61, -55, -41, -37, -25, -17)}, "A7 plus")

# ---------------------------------------------------------------- A8
hdr("A8. 1-cycle and 2-cycle word counts and outcomes")
one = {+1: [], -1: []}
for L in range(1, 3001):
    Kc = int(mpmath.floor(L * ALPHA))
    for K in range(Kc - 1, Kc + 3):
        k = K - (L - 1)
        if k < 1:
            continue
        w = (1,) * (L - 1) + (k,)
        Kt, Lt, B, Delta, q = gate(w)
        if q == 1:
            for b in (+1, -1):
                n0 = b * B // Delta
                if n0 > 0:
                    orbit = [n0]
                    for _ in range(L - 1):
                        orbit.append(Tb(orbit[-1], b))
                    if len(set(orbit)) == L:
                        one[b].append((L, k))
print(f"1-cycles L<=3000: plus {one[+1]} minus {one[-1]}")
check(one[+1] == [(1, 2)] and one[-1] == [(1, 1), (2, 2)], "A8 one-cycles")
n2 = 0
two = {+1: set(), -1: set()}
for L in range(2, 25):
    Kc = int(mpmath.floor(L * ALPHA))
    for K in range(Kc - 1, Kc + 3):
        for a in range(0, L - 1):
            bb = L - 2 - a
            for k1 in range(2, K - a - bb - 1):
                k2 = K - a - bb - k1
                if k2 < 2:
                    continue
                n2 += 1
                w = (1,) * a + (k1,) + (1,) * bb + (k2,)
                Kt, Lt, B, Delta, q = gate(w)
                if q == 1:
                    for b in (+1, -1):
                        n0 = b * B // Delta
                        if n0 > 0:
                            orbit = [n0]
                            for _ in range(L - 1):
                                orbit.append(Tb(orbit[-1], b))
                            two[b].add(tuple(sorted(set(orbit))))
print(f"2-cycle words L<=24 in the range K in floor(L a)+{{-1..2}}: {n2}; plus {sorted(two[+1])} minus {sorted(two[-1])}")
check(n2 == 9699 and two[+1] == {(1,)} and two[-1] == {(5, 7), (17, 25, 37, 41, 55, 61, 91)}, "A8 two-cycles")
print("range completeness (audit): a positive cycle NOT already known has all members > 10^7 (minus, A5)")
print("or >= 2^68 (plus, Barina), so |2^K/3^L - 1| < 1e-6 at L <= 24 and K = floor(L a) or +1 is forced.")

# ---------------------------------------------------------------- A9
hdr("A9. |Delta| <= 100 clocks in three boxes; necklace counts by Burnside")
restricted, allK1, allK0 = [], [], []
for L in range(1, 401):
    Kc = int(mpmath.floor(L * ALPHA))
    for K in range(0, 2 * L + 8):
        D = 2 ** K - 3 ** L
        if abs(D) <= 100:
            allK0.append((K, L, D))
            if K >= 1:
                allK1.append((K, L, D))
            if K in (Kc, Kc + 1):
                restricted.append((K, L, D))
print(f"K in {{floor,floor+1}}: {len(restricted)} {restricted}")
print(f"K >= 1 (pillai box): {len(allK1)};  K >= 0: {len(allK0)}")
check(len(restricted) == 9 and len(allK1) == 26 and len(allK0) == 30, "A9 clocks")
check([c for c in restricted if abs(c[2]) == 1] == [(1, 1, -1), (2, 1, 1), (3, 2, -1)], "A9 unit gaps")


def necklaces_burnside(K, L):
    # rotation classes of compositions of K into L parts = (1/L) sum_{d | gcd-compatible} phi(d) * #comps fixed
    tot = 0
    for d in range(1, L + 1):
        if L % d == 0 and K % d == 0:
            # rotation by L/d positions has period d; fixed words are d-fold repeats of a comp of K/d into L/d
            phi = sum(1 for r in range(1, d + 1) if math.gcd(r, d) == 1)
            tot += phi * math.comb(K // d - 1, L // d - 1)
    check(tot % L == 0, "Burnside integrality")
    return tot // L


for K, L in [(1, 1), (2, 1), (3, 2), (11, 7), (19, 12), (65, 41)]:
    print(f"  clock {K}/{L}: compositions {math.comb(K - 1, L - 1)}, necklaces {necklaces_burnside(K, L)}")
check(necklaces_burnside(11, 7) == 30 and necklaces_burnside(19, 12) == 2652
      and necklaces_burnside(65, 41) == 6113392816333320 and math.comb(64, 40) == 250649105469666120, "A9 necklaces")
q1 = sum(1 for w in comps(11, 7) if gate(w)[4] == 1)
print(f"  q=1 words on 11/7: {q1}")
check(q1 == 7, "A9 q=1 words")

# ---------------------------------------------------------------- A10
hdr("A10. D per period, budget sums, q_(7m)")
print("(the strict budget inequality on an infinite minus orbit is glued_xor (B10), which uses Garcia-Tal 1999 (CITED there); not re-proved here)")
for b, cyc in CYCLES:
    w = word_of_cycle(cyc, b)
    K, L, B, Delta, q = gate(w)
    D = float(K - L * ALPHA)
    if b < 0:
        Ks = [0]
        for k in w:
            Ks.append(Ks[-1] + k)
        period_sum = sum(Fraction(2 ** Ks[i], 3 ** i) for i in range(L))
        ratio = Fraction(2 ** K, 3 ** L)
        total = period_sum / (1 - ratio)
        print(f"b=-1 {cyc}: D/period={D:+.6f} period sum={period_sum} ratio={ratio} total={total} 3n0={3 * cyc[0]}")
        check(total == 3 * cyc[0], "A10 budget")
    else:
        print(f"b=+1 {cyc}: D/period={D:+.6f} ratio={Fraction(2 ** K, 3 ** L)} (sum diverges)")
r = Fraction(2048, 2187)
print("q_(7m) m=1..10:", " ".join(f"{float(r ** m):.4f}" for m in range(1, 11)))
check(abs(float(r ** 10) - 0.5186) < 6e-5, "A10 q_70")

# ---------------------------------------------------------------- A11
hdr("A11. Terras densities by brute force over residues (no DP)")
# odd-step convention: coefficient stopping <= J iff some j<=J has K_j > j*alpha; determined by n mod 2^13 for J<=8
M = 1 << 14
dJ = {}
for J in range(1, 9):
    cnt = 0
    for n in range(1, M, 2):
        m, Ksum = n, 0
        for j in range(1, J + 1):
            y = 3 * m + 1
            k = v2(y)
            m = y >> k
            Ksum += k
            if 2 ** Ksum > 3 ** j:
                cnt += 1
                break
    dJ[J] = Fraction(cnt, M // 2)
print("d_J (odd-step convention), J=1..8:", [str(dJ[J]) for J in range(1, 9)])
check([dJ[J] for J in range(1, 9)] == [Fraction(1, 2), Fraction(5, 8), Fraction(3, 4), Fraction(51, 64),
                                       Fraction(109, 128), Fraction(7, 8), Fraction(911, 1024), Fraction(3729, 4096)], "A11 d_J")
# Terras convention: U(n)=(3n+1)/2 or n/2, stop at first j with 3^o < 2^j; determined by n mod 2^J
FJ = {}
for J in (1, 2, 4, 8, 12):
    cnt = 0
    for n in range(1, (1 << J) + 1):
        m, o = n, 0
        for j in range(1, J + 1):
            if m & 1:
                m = (3 * m + 1) // 2
                o += 1
            else:
                m //= 2
            if 3 ** o < (1 << j):
                cnt += 1
                break
    FJ[J] = Fraction(cnt, 1 << J)
print("F_J (Terras convention) J=1,2,4,8,12:", [f"{FJ[J]} = {float(FJ[J]):.6f}" for J in (1, 2, 4, 8, 12)])
check(FJ[1] == Fraction(1, 2) and FJ[2] == Fraction(3, 4) and FJ[4] == Fraction(13, 16)
      and FJ[8] == Fraction(237, 256) and abs(float(FJ[12]) - 0.944824) < 1e-6, "A11 F_J")

# ---------------------------------------------------------------- A12
hdr("A12. squarefree growth prefixes L=5, q<=100000, sieve written from scratch")
Q = 100000
L = 5
for b in (+1, -1):
    ok = bytearray([1]) * (Q + 1)
    for j in range(L + 1):
        c = 3 ** j * 2 ** (L + 1 - j)
        maxv = c * Q + 1
        p = 2
        while p * p <= maxv:
            # p prime test by trial division (p small)
            if all(p % r for r in range(2, int(p ** 0.5) + 1)):
                pp = p * p
                if math.gcd(c, pp) == 1:
                    inv = pow(c, -1, pp)
                    r0 = (b * inv) % pp  # c q - b = 0 mod pp  <=>  q = b c^{-1}
                    for qq in range(r0 if r0 else pp, Q + 1, pp):
                        ok[qq] = 0
            p += 1
    cnt = sum(ok[1:])
    print(f"b={b:+d}: all 6 nodes squarefree for {cnt} of {Q} starts ({cnt / Q:.4f})")
    check(cnt == (49167 if b > 0 else 49173), "A12 squarefree")

# ---------------------------------------------------------------- A13
hdr("A13. Lyapunov obstruction witnesses")
for b in (+1, -1):
    n0 = 2 ** 13 * 6 - b
    m = n0
    res = set()
    for _ in range(12):
        res.add(m % 6)
        y = 3 * m + b
        check(v2(y) == 1, "A13 k=1")
        m = y >> 1
    print(f"b={b:+d}: n0={n0}, T^12={m}, residues mod 6 {sorted(res)}")
    check(m == 3 ** 12 * 12 - b and m > n0, "A13")

# ---------------------------------------------------------------- A14
hdr("A14. greedy map G with the g_negatives definition (2^k m in {4,7} mod 9 / {2,5} mod 9)")


def Gp(m):
    x = m
    while (x % 9) not in (4, 7):
        x *= 2
    return (x - 1) // 3


def Gm(m):
    x = m
    while (x % 9) not in (2, 5):
        x *= 2
    return (x + 1) // 3


bad = sum(1 for m in range(1, 3001) if m % 3 and Gp(-m) != -Gm(m))
print(f"violations of G(-m) = -G_-(m), m<=3000: {bad}")
check(bad == 0, "A14")

# ---------------------------------------------------------------- A16
hdr("A16. entropy of 8/pi^2 and the Paley tournament on F_7")
pv = 8 / math.pi ** 2
H = -pv * math.log2(pv) - (1 - pv) * math.log2(1 - pv)
print(f"8/pi^2 = {pv:.5f}, H = {H:.5f} bits")
check(abs(H - 0.70028) < 1e-5 and abs(pv - 0.81057) < 1e-5, "A16 entropy")
QR = {1, 2, 4}
arcs = {(i, j) for i in range(7) for j in range(7) if (j - i) % 7 in QR}
cyc3, trans3 = [], 0
for tri in combinations(range(7), 3):
    a, b_, c = tri
    outs = [sum(1 for y in tri if (x, y) in arcs) for x in tri]
    if sorted(outs) == [1, 1, 1]:
        cyc3.append(frozenset(tri))
    else:
        trans3 += 1
dev = {frozenset({(s + d) % 7 for d in base}) for base in ((0, 1, 3), (0, 1, 5)) for s in range(7)}
per_arc = {a: sum(1 for t in cyc3 if a[0] in t and a[1] in t) for a in arcs}
print(f"arcs {len(arcs)}, cyclic triples {len(cyc3)} ((7^3-7)/24 = {(343 - 7) // 24}), transitive {trans3}, "
      f"cyclic = dev{{0,1,3}} u dev{{0,1,5}}: {set(cyc3) == dev}, arcs per cyclic triple count {set(per_arc.values())}")
check(len(arcs) == 21 and len(cyc3) == 14 and trans3 == 21 and set(cyc3) == dev and set(per_arc.values()) == {2}, "A16 Paley")
# octonion index rule on the lines {s,s+1,s+3}: s->s+1->s+3->s
rule = all(((s, (s + 1) % 7) in arcs and ((s + 1) % 7, (s + 3) % 7) in arcs and ((s + 3) % 7, s) in arcs) for s in range(7))
print(f"on every line {{s,s+1,s+3}} the orientation is s->s+1->s+3->s: {rule}")
check(rule, "A16 rule")

print(f"\nall audit checks passed; total {time.time() - T0:.1f}s (timing line, not a result)")
