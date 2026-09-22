#!/usr/bin/env python3
"""Adversarial audit of lane pillai_convergents_cycle_gates (2026-09-22).

Independent recomputation of every key number of
05-knowledge/results/collatz_mod6_20260921_pillai_convergents_cycle_gates.md
with different code paths: mpmath (not PARI, not decimal) for the continued
fraction; an explicit mediant enumeration for the intermediate fractions; a
Horner recursion for the carry B(w) and itertools cut positions for the
compositions (the lane uses a recursive DFS with the closed sum); a residue
DP mod |Delta| as a second, enumeration-free count of the q=1 words on the
small-gap clocks; sympy isprime/products for the factorizations; exact
Fractions for the Eliahou bounds; an empirical test of the Legendre and
Fatou-Grace placement statements on five irrationals; and a wider box for
the two Pillai censuses.  All checks raise explicitly; stdout is identical
under python3 -O.
"""
import itertools
import math
import sys
import time
from fractions import Fraction
from math import comb, gcd

import mpmath
from sympy import isprime

T0 = time.time()


def out(*a):
    print(*a)
    sys.stdout.flush()


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


out("=" * 78)
out("A1. Continued fraction of log_2 3 by mpmath (400 dps), convergents, intermediates")
out("=" * 78)
mpmath.mp.dps = 400
alpha = mpmath.log(3) / mpmath.log(2)
NT = 45
cf = []
x = alpha
for _ in range(NT):
    a = int(mpmath.floor(x))
    cf.append(a)
    x = 1 / (x - a)
LANE_CF = [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5,
           7, 1, 1, 4, 8, 1, 11, 1, 20, 2, 1, 10, 1, 4, 1, 1, 1, 1, 1, 37]
check(cf == LANE_CF, "continued fraction differs from the lane: %s" % cf)
out("45 partial quotients agree with the lane:", cf == LANE_CF)
out("largest among the 45: %d at index %d; second largest: %d at index %d; a_20=%d, a_44=%d"
    % (max(cf), cf.index(max(cf)), sorted(set(cf))[-2], cf.index(sorted(set(cf))[-2]), cf[20], cf[44]))


def convergents(cf):
    p0, q0, p1, q1 = 1, 0, cf[0], 1
    res = [(p1, q1)]
    for a in cf[1:]:
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        res.append((p1, q1))
    return res


conv = convergents(cf)
LANE_CONV = [(1, 1), (2, 1), (3, 2), (8, 5), (19, 12), (65, 41), (84, 53), (485, 306), (1054, 665),
             (24727, 15601), (50508, 31867), (125743, 79335), (176251, 111202), (301994, 190537),
             (16785921, 10590737)]
check(conv[:15] == LANE_CONV, "convergents differ")
for n in range(1, 15):
    (p0, q0), (p1, q1) = conv[n - 1], conv[n]
    check(abs(p0 * q1 - p1 * q0) == 1, "determinant at n=%d" % n)
out("convergents n=0..14 agree with the lane and satisfy p_(n-1) q_n - p_n q_(n-1) = +-1: True")
# exact integer check of |Delta| L < 3^L on convergents with L <= 200000 (indices 0..13)
nchk = 0
for (K, L) in conv:
    if L > 200000:
        break
    D = 2 ** K - 3 ** L
    check(L == 1 or abs(D) * L < 3 ** L, "convergent property at %d/%d" % (K, L))
    nchk += 1
out("convergents with L <= 200000 checked exactly for |Delta|*L < 3^L: %d" % nchk)


def intermediates(cf, conv, nmax):
    """(K, L, n, j, a_{n+1}) for the mediant chain between p_{n-1}/q_{n-1} and p_{n+1}/q_{n+1}."""
    res = {}
    for n in range(0, nmax):
        pm, qm = (1, 0) if n == 0 else conv[n - 1]
        pn, qn = conv[n]
        p, q = pm, qm
        for j in range(1, cf[n + 1]):
            p, q = p + pn, q + qn          # repeated mediant
            res[(p, q)] = (n, j, cf[n + 1])
    return res


inter = intermediates(cf, conv, 40)
conv_idx = {pq: n for n, pq in enumerate(conv)}


def classify(K, L):
    g = gcd(K, L)
    r = (K // g, L // g)
    tag = "" if g == 1 else " (x%d)" % g
    if r in conv_idx:
        return "convergent n=%d%s" % (conv_idx[r], tag)
    if r in inter:
        n, j, a = inter[r]
        return "intermediate n=%d j=%d/%d %s%s" % (n, j, a - 1, "extreme" if j in (1, a - 1) else "INTERIOR", tag)
    return "NEITHER"


for (K, L, expect) in ((5, 3, "intermediate n=2 j=1/1 extreme"), (11, 7, "intermediate n=3 j=1/1 extreme"),
                       (27, 17, "intermediate n=4 j=1/2 extreme"), (46, 29, "intermediate n=4 j=2/2 extreme"),
                       (149, 94, "intermediate n=6 j=1/4 extreme"), (233, 147, "intermediate n=6 j=2/4 INTERIOR"),
                       (7, 4, "NEITHER"), (19, 12, "convergent n=4"), (6, 4, "convergent n=2 (x2)")):
    check(classify(K, L) == expect, "classify(%d,%d)=%s" % (K, L, classify(K, L)))
out("classification of 5/3, 11/7, 27/17, 46/29, 149/94, 233/147, 7/4, 19/12, 6/4 agrees with the lane")
for (K, L, D) in ((5, 3, 5), (8, 5, 13), (11, 7, -139), (19, 12, -7153), (27, 17, 5077565),
                  (46, 29, 1738366812781), (65, 41, 420491770248316829), (84, 53, -40432553845953101497907)):
    check(2 ** K - 3 ** L == D, "Delta at %d/%d" % (K, L))
out("Delta at 5/3, 8/5, 11/7, 19/12, 27/17, 46/29, 65/41, 84/53 agree with the lane's table")
out("|Delta|/3^L at 1/1, 3/2, 5/3, 11/7, 19/12: %s" % ", ".join(
    "%.4e" % (abs(2 ** K - 3 ** L) / 3 ** L) for (K, L) in ((1, 1), (3, 2), (5, 3), (11, 7), (19, 12))))

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A2. Factorizations quoted in the note: products and primality (sympy isprime)")
out("=" * 78)
FACTS = {
    (19, 12): [23, 311],
    (27, 17): [5, 71, 14303],
    (46, 29): [39409, 44110909],
    (65, 41): [19, 29, 17021, 44835377399],
    (84, 53): [11, 467, 11979433, 657030219467],
    (149, 94): [7, 7, 7, 30809, 1765553, 1009273943, 353617911534038093791],
    (485, 306): [929, 84958721, 1437465479, 46777127526357837196396057,
                 19231970699168568692206159641463898527274405039282219231295668859629511743697206424938341838460889],
}
for (K, L), fs in FACTS.items():
    prod = 1
    for f in fs:
        prod *= f
        check(isprime(f), "non-prime factor %d at %d/%d" % (f, K, L))
    check(prod == abs(2 ** K - 3 ** L), "product mismatch at %d/%d" % (K, L))
    out("  |2^%d-3^%d| = %s : product and primality OK (%d-digit largest factor)" % (K, L, "*".join(map(str, fs)), len(str(max(fs)))))
for (K, L) in ((5, 3), (8, 5), (11, 7)):
    check(isprime(abs(2 ** K - 3 ** L)), "prime gap at %d/%d" % (K, L))
out("  |Delta| = 5, 13, 139 prime at 5/3, 8/5, 11/7: OK")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A3. Tier constants, the elementary log inequality, and the L <= 400 hit list")
out("=" * 78)
ln2 = math.log(2)
out("3 ln 2 = %.6f > 2: %s;  2/(3 ln 2) = %.6f < 1: %s;  2/ln 2 = %.6f < 2.9: %s"
    % (3 * ln2, 3 * ln2 > 2, 2 / (3 * ln2), 2 / (3 * ln2) < 1, 2 / ln2, 2 / ln2 < 2.9))
check(3 * ln2 > 2 and 2 / (3 * ln2) < 1 and 2 / ln2 < 2.9, "tier constants")
mpmath.mp.dps = 50
worst = mpmath.mpf(0)
for i in range(1, 100000):
    e = mpmath.mpf(i) / 400000          # |e| < 0.25
    for s in (1, -1):
        v = abs(mpmath.log(1 + s * e, 2))
        bnd = e / ((1 - e) * mpmath.log(2))
        check(v <= bnd, "log inequality at %s" % (s * e))
        worst = max(worst, v / bnd)
out("|log_2(1+e)| <= |e|/((1-|e|) ln 2) on 199998 grid points |e| < 1/4 (mpmath 50 dps): max ratio %s" % mpmath.nstr(worst, 6))
out("(the lemma is also proved by ln(1+e) <= e for e >= 0 and -ln(1-|e|) = int_0^|e| dt/(1-t) <= |e|/(1-|e|))")

mpmath.mp.dps = 60
alpha60 = mpmath.log(3) / mpmath.log(2)
hits = []
for L in range(1, 401):
    P = 3 ** L
    Kc = int(mpmath.floor(L * alpha60))
    for K in range(max(1, Kc - 3), Kc + 5):
        D = 2 ** K - P
        if abs(D) * L < P:
            hits.append((K, L, D, 4 * L * abs(D) <= P, 2 * L * abs(D) <= P, classify(K, L)))
out("tier-C hits (|Delta| L < 3^L), L <= 400: %d" % len(hits))
for (K, L, D, tA, tB, cls) in hits:
    out("  %3d/%-3d A=%d B=%d %s" % (K, L, tA, tB, cls))
check(len(hits) == 19, "hit count")
check(sum(1 for h in hits if h[3]) == 3 and all(h[5].startswith("convergent") for h in hits if h[3]), "tier A")
check(sum(1 for h in hits if h[4]) == 11 and all(h[5].startswith("convergent") or "extreme" in h[5] for h in hits if h[4]), "tier B")
check(all(h[5].startswith("convergent") or "extreme" in h[5] for h in hits), "tier C empirical")
check([(h[0], h[1]) for h in hits if h[3]] == [(3, 2), (19, 12), (84, 53)], "tier A list")
out("tier A: 3 hits (3/2, 19/12, 84/53), all convergents, all Delta < 0: %s" % all(2 ** K - 3 ** L < 0 for (K, L) in ((3, 2), (19, 12), (84, 53))))
out("tier B: 11 hits, all convergents or extreme intermediates; tier C: all 19 convergent-or-extreme (no INTERIOR, no NEITHER)")
out("hits with L > 60 (not listed in the lane's .out): %s" % ", ".join("%d/%d" % (K, L) for (K, L, D, tA, tB, cls) in hits if L > 60))
# multiple K per L can only happen for small L: for L >= 6, 2.885/L < 1/2
per_L = {}
for h in hits:
    per_L.setdefault(h[1], []).append(h[0])
out("L values with two or more hits: %s (bound 2/(L ln 2) < 1/2 for L >= 6)" % {L: Ks for L, Ks in per_L.items() if len(Ks) > 1})

# Legendre and Fatou-Grace, empirical, on five irrationals with q <= 2000
out()
out("empirical placement test, all reduced p/q with q <= 2000 and |x - p/q| < 1/q^2 (resp. < 1/(2q^2)):")
tests = [("log_2 3", alpha60), ("sqrt 2", mpmath.sqrt(2)), ("pi", mpmath.pi), ("e", mpmath.e), ("2^(1/3)", mpmath.cbrt(2))]
for name, xval in tests:
    cfx = []
    y = xval
    while True:
        a = int(mpmath.floor(y))
        cfx.append(a)
        if convergents(cfx)[-1][1] > 4000:
            break
        y = 1 / (y - a)
    cvx = convergents(cfx)
    cvx_set = {pq: n for n, pq in enumerate(cvx)}
    intx = intermediates(cfx, cvx, len(cfx) - 1)
    n_fg = n_leg = n_bad_fg = n_bad_leg = n_interior_between = 0
    for q in range(1, 2001):
        for p in (int(mpmath.floor(q * xval)), int(mpmath.floor(q * xval)) + 1):
            if gcd(p, q) != 1:
                continue
            err = abs(xval - mpmath.mpf(p) / q)
            if err < mpmath.mpf(1) / q ** 2:
                n_fg += 1
                ok = (p, q) in cvx_set or ((p, q) in intx and intx[(p, q)][1] in (1, intx[(p, q)][2] - 1))
                if not ok:
                    n_bad_fg += 1
            if err < mpmath.mpf(1) / (2 * q ** 2):
                n_leg += 1
                if (p, q) not in cvx_set:
                    n_bad_leg += 1
            if mpmath.mpf(1) / q ** 2 <= err < mpmath.mpf(3) / q ** 2 and (p, q) in intx and intx[(p, q)][1] not in (1, intx[(p, q)][2] - 1):
                n_interior_between += 1
    check(n_bad_fg == 0 and n_bad_leg == 0, "placement failure for %s" % name)
    out("  %-8s  |x-p/q|<1/q^2: %3d fractions, all convergent-or-extreme: %s;  <1/(2q^2): %3d, all convergents: %s;  interior intermediates with 1/q^2 <= err < 3/q^2: %d"
        % (name, n_fg, n_bad_fg == 0, n_leg, n_bad_leg == 0, n_interior_between))

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A4. Eliahou product bound on the eight known cycles; thresholds; L* at N = 2^68")
out("=" * 78)


def horner_B(w):
    B, Ki = 0, 0
    for k in w:
        B = 3 * B + 2 ** Ki
        Ki += k
    return B, Ki


def word_of_cycle(nodes, b):
    w = []
    for i, n in enumerate(nodes):
        m = 3 * n + b
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        check(m == nodes[(i + 1) % len(nodes)], "not a cycle")
        w.append(k)
    return w


KNOWN = [((1,), 1), ((-1,), 1), ((-5, -7), 1), ((-17, -25, -37, -55, -41, -61, -91), 1),
         ((1,), -1), ((5, 7), -1), ((17, 25, 37, 55, 41, 61, 91), -1), ((-1,), -1)]
for nodes, b in KNOWN:
    w = word_of_cycle(nodes, b)
    B, K = horner_B(w)
    L = len(w)
    D = 2 ** K - 3 ** L
    q = abs(D) // gcd(B, abs(D))
    check(q == 1 and Fraction(b * B, D) == nodes[0], "gate on %s" % (nodes,))
    N = min(abs(v) for v in nodes)
    rel = Fraction(abs(D), 3 ** L)
    bnd = (1 + Fraction(1, 3 * N)) ** L - 1 if nodes[0] * b > 0 else 1 - (1 - Fraction(1, 3 * N)) ** L
    check(rel <= bnd, "Eliahou bound on %s" % (nodes,))
    # sign purity: 3n+b has the sign of n for every odd n (|3n| >= 3 > 1)
    out("  %-32s b=%2d w=%s K/L=%d/%d Delta=%d B=%d q=%d N=%d rel=%s=%.4f bound=%s=%.4f"
        % (nodes, b, w, K, L, D, B, q, N, rel, float(rel), bnd, float(bnd)))
check(1 - (1 - Fraction(1, 15)) ** 2 == Fraction(29, 225), "bound 29/225")
out("  seven-cycle word rotations: B(rot w) = (3B + Delta)/2^k_1 (inherited signed_cycles sec.3), n_0 = -B/Delta:")
w = [1, 1, 1, 2, 1, 1, 4]
B, K = horner_B(w)
D = 2 ** 11 - 3 ** 7
rots = []
for i in range(7):
    rots.append((tuple(w), B, -B // D))
    check((-B) % D == 0, "integrality")
    B = (3 * B + D) // 2 ** w[0]
    w = w[1:] + w[:1]
out("  " + "; ".join("%s B=%d n_0=%d" % r for r in rots))
check([r[1] for r in rots] == [2363, 3475, 5143, 7645, 5699, 8479, 12649], "rotation carries")
check(sorted(r[2] for r in rots) == [17, 25, 37, 41, 55, 61, 91], "rotation nodes")
out("  (the lane's .out lists the same seven (B, n_0) pairs in lexicographic word order)")

mpmath.mp.dps = 40
out("Legendre threshold N_min(L) = 1/(3((1+1/(4L))^(1/L)-1)) versus 4L^2/3 (mpmath):")
for L in (5, 12, 41, 53, 306, 665, 15601, 31867, 79335):
    Nmin = 1 / (3 * ((1 + mpmath.mpf(1) / (4 * L)) ** (mpmath.mpf(1) / L) - 1))
    out("  L=%-6d N_min = %s   4L^2/3 = %s   ratio = %s" % (L, mpmath.nstr(Nmin, 5), mpmath.nstr(mpmath.mpf(4 * L * L) / 3, 5), mpmath.nstr(Nmin / (mpmath.mpf(4 * L * L) / 3), 6)))
N68 = mpmath.mpf(2) ** 68
Lstar = 14878195707


def tierA_ok(L):
    return mpmath.expm1(mpmath.mpf(L) / (3 * N68)) <= mpmath.mpf(1) / (4 * L)


check(tierA_ok(Lstar), "L* at 2^68")
lo, hi = Lstar, 2 * Lstar
while hi - lo > 1:
    mid = (lo + hi) // 2
    if tierA_ok(mid):
        lo = mid
    else:
        hi = mid
check(tierA_ok(lo) and not tierA_ok(lo + 1), "exact threshold")
out("exp(L/(3*2^68)) - 1 <= 1/(4L) holds at the lane's L = %d (a sufficient bound found with a 1e-6 float margin);"
    % Lstar)
out("  the exact largest L is %d (mpmath 40 dps bisection), i.e. %d above the lane's; isqrt(3*2^68/4) = %d"
    % (lo, lo - Lstar, math.isqrt(3 * 2 ** 68 // 4)))
check(math.isqrt(3 * 2 ** 68 // 4) == 14878203147, "isqrt")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A5. Cycle-gate census, independent enumeration (cut positions + Horner carry)")
out("=" * 78)


def census(K, L):
    D = abs(2 ** K - 3 ** L)
    n_words = n_q1 = 0
    min_q = None
    hist = {}
    q1_words = []
    pw2 = [2 ** i for i in range(K + 1)]
    for cuts in itertools.combinations(range(1, K), L - 1):
        B = 0
        prev = 0
        for c in cuts:
            B = 3 * B + pw2[prev]
            prev = c
        B = 3 * B + pw2[prev]
        q = D // gcd(B, D)
        n_words += 1
        if q == 1:
            n_q1 += 1
            if len(q1_words) < 40:
                parts = [c - p for p, c in zip((0,) + cuts, cuts + (K,))]
                q1_words.append((tuple(parts), B))
        if min_q is None or q < min_q:
            min_q = q
        if q <= 20:
            hist[q] = hist.get(q, 0) + 1
    check(n_words == comb(K - 1, L - 1), "word count")
    return n_words, n_q1, min_q, hist, q1_words


def census_dp(K, L):
    """q=1 count and min q from a residue DP mod D (no enumeration of words)."""
    D = abs(2 ** K - 3 ** L)
    # state after placing i parts: (K_i, B_i mod D) -> count, where B_i = Horner carry over the first i parts
    layer = {(0, 0): 1}
    for i in range(L):
        nxt = {}
        for (Ki, Bm), c in layer.items():
            Bn = (3 * Bm + pow(2, Ki, D)) % D
            parts_left = L - 1 - i
            if parts_left == 0:
                key = (K, Bn)
                nxt[key] = nxt.get(key, 0) + c
            else:
                for k in range(1, K - Ki - parts_left + 1):
                    key = (Ki + k, Bn)
                    nxt[key] = nxt.get(key, 0) + c
        layer = nxt
    total = sum(layer.values())
    n_q1 = sum(c for (Ki, Bm), c in layer.items() if Bm == 0)
    min_q = min(D // gcd(Bm, D) for (Ki, Bm) in layer)
    return total, n_q1, min_q


# expected rows (words, #q=1, min_q) from the lane's .out, keyed by (K, L)
EXPECT = {
    (1, 1): (1, 1, 1), (2, 1): (1, 1, 1), (2, 2): (1, 1, 1), (3, 2): (2, 2, 1), (4, 2): (3, 1, 1),
    (3, 3): (1, 1, 1), (5, 3): (6, 0, 5), (6, 3): (10, 1, 1), (4, 4): (1, 1, 1), (6, 4): (10, 2, 1),
    (8, 4): (35, 1, 1), (5, 5): (1, 1, 1), (8, 5): (35, 0, 13), (10, 5): (126, 1, 1), (6, 6): (1, 1, 1),
    (9, 6): (56, 2, 1), (10, 6): (126, 0, 5), (12, 6): (462, 1, 1), (7, 7): (1, 1, 1), (11, 7): (210, 7, 1),
    (14, 7): (1716, 1, 1), (8, 8): (1, 1, 1), (12, 8): (330, 2, 1), (16, 8): (6435, 1, 1), (9, 9): (1, 1, 1),
    (15, 9): (3003, 0, 5), (18, 9): (24310, 1, 1), (10, 10): (1, 1, 1), (15, 10): (2002, 2, 1),
    (16, 10): (5005, 0, 13), (20, 10): (92378, 1, 1), (11, 11): (1, 1, 1), (22, 11): (352716, 1, 1),
    (12, 12): (1, 1, 1), (18, 12): (12376, 2, 1), (19, 12): (31824, 0, 23), (20, 12): (75582, 0, 5),
    (24, 12): (1352078, 1, 1), (13, 13): (1, 1, 1), (26, 13): (5200300, 1, 1),
}
clocks = [(K, L) for L in range(1, 14) for K in range(L, 2 * L + 1) if classify(K, L) != "NEITHER"]
check(len(clocks) == 40 and set(clocks) == set(EXPECT), "clock list")
out("convergent/intermediate clocks with L <= 13: %d (same set as the lane)" % len(clocks))
out("clock   words    #q=1  min_q  hist(q<=20)             DP(mod |Delta|, run when |Delta| <= 30000) total/#q=1/min_q")
q1_pos, q1_neg = [], []
n_dp = 0
for (K, L) in clocks:
    D = 2 ** K - 3 ** L
    n_words, n_q1, min_q, hist, q1w = census(K, L)
    check((n_words, n_q1, min_q) == EXPECT[(K, L)], "row %d/%d: %s vs %s" % (K, L, (n_words, n_q1, min_q), EXPECT[(K, L)]))
    dp = ""
    if abs(D) <= 30000:
        t2, n2, m2 = census_dp(K, L)
        check((t2, n2, m2) == (n_words, n_q1, min_q), "DP disagrees at %d/%d" % (K, L))
        dp = "%d/%d/%d agree" % (t2, n2, m2)
        n_dp += 1
    out("  %2d/%-2d %-8d %-5d %-6d %-24s %s" % (K, L, n_words, n_q1, min_q, " ".join("%d:%d" % (q, hist[q]) for q in sorted(hist)), dp))
    for (wd, B) in q1w:
        (q1_pos if D > 0 else q1_neg).append((K, L, wd, B))
check(len(q1_pos) == 13 and all(set(wd) == {2} for (K, L, wd, B) in q1_pos), "3n+1 side")
check(len(q1_neg) == 32, "3n-1 side count")
out("residue DP agrees with the enumeration on all %d clocks with |Delta| <= 30000" % n_dp)
out("q=1 words: 3n+1 side %d (all (2,...,2)), 3n-1 side %d" % (len(q1_pos), len(q1_neg)))
# closed forms for the repeated fixed points and the alternating two-cycle words (PROVED by geometric sums)
for d in range(1, 14):
    check(horner_B([2] * d)[0] == 4 ** d - 3 ** d == 2 ** (2 * d) - 3 ** d, "all-2 carry")
    check(horner_B([1] * d)[0] == 3 ** d - 2 ** d, "all-1 carry")
for d in range(1, 7):
    check(horner_B([1, 2] * d)[0] == 5 * (9 ** d - 8 ** d), "(1,2)^d carry")
    check(horner_B([2, 1] * d)[0] == 7 * (9 ** d - 8 ** d), "(2,1)^d carry")
out("closed forms verified: B((2)^d) = 4^d - 3^d = Delta (n_0 = 1); B((1)^d) = 3^d - 2^d = -Delta (n_0 = 1);")
out("  B((1,2)^d) = 5 (9^d - 8^d) and B((2,1)^d) = 7 (9^d - 8^d) with Delta = 8^d - 9^d (n_0 = 5, 7), d <= 6")
grouped = {}
for (K, L, wd, B) in q1_neg:
    grouped[(K, L)] = grouped.get((K, L), 0) + 1
check(grouped[(11, 7)] == 7, "seven on 11/7")
out("3n-1 q=1 words per clock: %s" % {k: grouped[k] for k in sorted(grouped)})
# the six q=5 words on (10,6), (15,9), (20,12) are the d-fold repeats of the six words on 5/3
for d in (2, 3, 4):
    n_words, n_q1, min_q, hist, q1w = census(5 * d, 3 * d)
    check(hist.get(5, 0) == 6 and n_q1 == 0, "(5d,3d) at d=%d" % d)
    D = abs(2 ** (5 * d) - 3 ** (3 * d))
    reps = 0
    for wd in ((1, 1, 3), (1, 2, 2), (1, 3, 1), (2, 1, 2), (2, 2, 1), (3, 1, 1)):
        B = horner_B(list(wd) * d)[0]
        reps += (D // gcd(B, D) == 5)
    check(reps == 6, "repeats")
out("on (5d,3d), d=2,3,4: the six q=5 words are exactly the d-fold repeats of the six words of 5/3 (checked)")
# 19/12: full q spectrum
n_words, n_q1, min_q, hist, q1w = census(19, 12)
D = 7153
spec = {}
for cuts in itertools.combinations(range(1, 19), 11):
    B, prev = 0, 0
    for c in cuts:
        B = 3 * B + 2 ** prev
        prev = c
    B = 3 * B + 2 ** prev
    q = D // gcd(B, D)
    spec[q] = spec.get(q, 0) + 1
out("19/12 (|Delta| = 7153 = 23*311): q spectrum over all 31824 words: %s" % {q: spec[q] for q in sorted(spec)})
check(min(spec) == 23 and 1 not in spec, "19/12 spectrum")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A6. Full-clock census L <= 11 (all K in [L, 2L])")
out("=" * 78)
LANE_Q1 = [2, 4, 2, 4, 2, 4, 9, 4, 2, 4, 2]
LANE_NC = [2, 1, 2, 2, 4, 2, 3, 2, 5, 2, 6]
tot = 0
q1_per_L = []
nc_per_L = []
for L in range(1, 12):
    q1 = nc = 0
    mins = []
    words11 = []
    for K in range(L, 2 * L + 1):
        n_words, n_q1, min_q, hist, q1w = census(K, L)
        tot += n_words
        q1 += n_q1
        nc += (min_q == abs(2 ** K - 3 ** L))
        mins.append((min_q, K))
        if L == 11:
            words11 += [wd for (wd, B) in q1w]
    mins.sort()
    q1_per_L.append(q1)
    nc_per_L.append(nc)
    out("  L=%2d q=1 words=%d no-cancel clocks=%d/%d three smallest min_q: %s" % (L, q1, nc, L + 1, mins[:3]))
    if L == 11:
        check(sorted(words11) == [tuple([1] * 11), tuple([2] * 11)], "L=11 q=1 words")
        out("  L=11 q=1 words are exactly (1^11) and (2^11): no signed cycle of least period 11 at b=+-1")
check(q1_per_L == LANE_Q1 and nc_per_L == LANE_NC, "per-L counts")
check(tot == 956384 and sum(q1_per_L[:10]) == 37, "totals")
out("total words %d; q=1 per L %s (sum over L<=10: %d); no-cancellation clocks per L %s (L=1 degenerate: |Delta|=1)"
    % (tot, q1_per_L, sum(q1_per_L[:10]), nc_per_L))

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("A7. Pillai censuses in wider boxes")
out("=" * 78)
sols = [(a, b, 2 ** a - 10 ** b) for a in range(0, 301) for b in range(0, 91) if abs(2 ** a - 10 ** b) <= 100]
out("|2^a - 10^b| <= 100, 0 <= a <= 300, 0 <= b <= 90: %d solutions; with a, b >= 1: %d; largest b: %d; largest a: %d"
    % (len(sols), sum(1 for s in sols if s[0] >= 1 and s[1] >= 1), max(s[1] for s in sols), max(s[0] for s in sols)))
check(len(sols) == 23 and (10, 3, 24) in sols, "2-10 census")
kl = [(K, L, 2 ** K - 3 ** L) for K in range(1, 701) for L in range(1, 401) if abs(2 ** K - 3 ** L) <= 100]
out("|2^K - 3^L| <= 100, 1 <= K <= 700, 1 <= L <= 400: %d solutions (lane box K<=200, L<=130 had 26); largest L: %d; NEITHER: %d; tier-C: %d"
    % (len(kl), max(s[1] for s in kl), sum(1 for s in kl if classify(s[0], s[1]) == "NEITHER"),
       sum(1 for s in kl if abs(s[2]) * s[1] < 3 ** s[1])))
check(len(kl) == 26 and (7, 4, 47) in kl, "2-3 census")
# completeness inside L <= 400: |Delta| <= 100 forces 2^K <= 3^400 + 100, i.e. K <= 634 < 700
check(2 ** 634 > 3 ** 400 + 100 and 2 ** 633 < 3 ** 400, "K bound at L = 400")
out("for L <= 400, |Delta| <= 100 forces K <= 634 < 700, so this list is complete for ALL K and all L <= 400 (FINITE-EXACT);")
check(all(100 * 4 * L <= 3 ** L for L in range(8, 401)), "100 <= 3^L/(4L) for L >= 8")
out("  beyond that, 100 <= 3^L/(4L) for every L >= 8, so any further solution would be in tier A and sit on a convergent")
out("  multiple (PROVED via Legendre); no elementary bound on L follows, which is why the full list stays CITED context.")

print("elapsed %.1f s" % (time.time() - T0), file=sys.stderr)
out("AUDIT DONE")
