#!/usr/bin/env python3
"""Lane pillai_convergents_cycle_gates (session collatz-mod6, 2026-09-21).

Continued fraction of log_2 3, the power gap Delta = 2^K - 3^L along its
convergents and intermediate fractions, the Legendre / Fatou-Grace placement
of small-gap clocks, the Eliahou cycle-to-convergent bridge, an exact
cycle-gate census (q = |Delta| / gcd(B, |Delta|)) on every clock (K, L) with
L <= 13 sitting on a convergent or intermediate fraction, a full-clock
minimum-q census for L <= 11, and the complete solution set of
|2^a - 10^b| <= 100.

Inherited (cited by path in the note, never re-derived here):
  - cycle gate n_0 = bB/(2^K-3^L), q = |Delta|/gcd(B,|Delta|), q | b criterion:
    05-knowledge/results/arithmetic_braids2_20260917_signed_cycles.md
  - unit-gap classification, clock bound L <= K <= 2L, period-10 certificate:
    05-knowledge/results/catalan_elliptic_20260921_catalan.md
  - c(b) census and necklace theorem:
    05-knowledge/results/collatz_mod6_20260917_wild_typing.md

All checks raise explicitly; the script must produce identical stdout under
python3 and python3 -O (timing goes to stderr only).
"""
import math
import subprocess
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb, gcd

T0 = time.time()


def out(*a):
    print(*a)
    sys.stdout.flush()


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def gp(cmd, timeout=300):
    r = subprocess.run(["gp", "-q"], input=cmd + "\n", capture_output=True,
                       text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError("gp failed: " + r.stderr)
    return r.stdout.strip()


# ---------------------------------------------------------------------------
out("=" * 78)
out("S1. Continued fraction of alpha = log_2 3 (PARI/GP at 320 digits, cross-checked")
out("    by an independent Python decimal computation at 320 digits)")
out("=" * 78)
NTERMS = 45
gp_out = gp("default(realprecision,320); v=contfrac(log(3)/log(2)); "
            "print(vector(%d,i,v[i]))" % NTERMS)
cf_gp = [int(t) for t in gp_out.strip("[]").split(",")]
check(len(cf_gp) == NTERMS, "gp returned %d terms" % len(cf_gp))

getcontext().prec = 330
alpha_dec = Decimal(3).ln() / Decimal(2).ln()
x = alpha_dec
cf_py = []
for _ in range(NTERMS):
    a = int(x.to_integral_value(rounding="ROUND_FLOOR"))
    cf_py.append(a)
    frac = x - a
    x = 1 / frac
check(cf_py == cf_gp, "gp and decimal continued fractions disagree")
out("partial quotients a_0..a_%d:" % (NTERMS - 1))
out(" ", cf_gp)
out("index: value  (index 0 is the integer part)")
for i in (9, 14, 20, 44):
    out("  a_%d = %d" % (i, cf_gp[i]))
out("NOTE (numerology flag): a_9=23 and a_14=55 are the two largest early")
out("  partial quotients (0-indexed); in 1-indexed counting they are the 10th")
out("  and 15th terms.  23 = 3n+1 reset residue (guards lane) and 55 = F_10 are")
out("  coincidences of small integers; no map from partial quotients to Collatz")
out("  residues or Fibonacci numbers is claimed.")


def convergents(cf):
    p0, q0, p1, q1 = 1, 0, cf[0], 1
    res = [(p1, q1)]
    for a in cf[1:]:
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        res.append((p1, q1))
    return res


conv = convergents(cf_gp)
out("convergents p_n/q_n, n=0..14:")
for n in range(15):
    out("  n=%2d  %d/%d" % (n, conv[n][0], conv[n][1]))

# exact integer sanity check of the convergent property |Delta| * L < 3^L
# (for a convergent with L>=2, |L alpha - K| < 1/L gives |2^K/3^L - 1| < 1/L)
out("exact integer check |2^K-3^L|*L < 3^L on convergents with L <= 200000:")
n_checked = 0
for (K, L) in conv:
    if L > 200000:
        break
    D = 2 ** K - 3 ** L
    check(L == 1 or abs(D) * L < 3 ** L, "convergent property fails at %d/%d" % (K, L))
    n_checked += 1
out("  convergents checked exactly:", n_checked, "(indices 0..%d)" % (n_checked - 1))

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S2. Intermediate fractions (semiconvergents) and Delta = 2^K - 3^L")
out("=" * 78)


def intermediates(cf, conv, max_index):
    """All (K, L, n, j) with K/L = (p_{n-1} + j p_n)/(q_{n-1} + j q_n),
    1 <= j <= a_{n+1}-1, n >= 1 (family index n, semiconvergent between
    p_{n-1}/q_{n-1} and p_{n+1}/q_{n+1}).  For n=0 the family is
    (1 + j*p_0)/(0 + j*q_0) = (1+j)/j, 1 <= j <= a_1 - 1 (empty here since a_1=1)."""
    res = []
    for n in range(0, max_index):
        a_next = cf[n + 1]
        if n == 0:
            pm, qm = 1, 0
        else:
            pm, qm = conv[n - 1]
        pn, qn = conv[n]
        for j in range(1, a_next):
            res.append((pm + j * pn, qm + j * qn, n, j, a_next))
    return res


inter = intermediates(cf_gp, conv, 12)
rows = []
for n in range(12):
    K, L = conv[n]
    rows.append((L, K, "convergent n=%d" % n))
for (K, L, n, j, a_next) in inter:
    if L <= 10 ** 7:
        extreme = (j == 1 or j == a_next - 1)
        rows.append((L, K, "intermediate n=%d j=%d/%d%s" % (n, j, a_next - 1,
                     " (extreme)" if extreme else "")))
rows.sort()
out("clock  K/L      kind                          Delta=2^K-3^L (digits if long)   |Delta|/3^L")
delta_table = {}
for (L, K, kind) in rows:
    D = 2 ** K - 3 ** L
    delta_table[(K, L)] = D
    rel = abs(D) / 3 ** L if L < 1000 else float(Fraction(abs(D), 3 ** L))
    ds = str(D) if abs(D) < 10 ** 30 else ("%s(%d digits)" % ("-" if D < 0 else "+", len(str(abs(D)))))
    out("  %5d/%-5d %-32s %-32s %.3e" % (K, L, kind, ds, rel))
check(delta_table[(11, 7)] == -139, "11/7 gap")
check(delta_table[(19, 12)] == -7153, "19/12 gap")
check(delta_table[(5, 3)] == 5 and delta_table[(8, 5)] == 13, "5/3, 8/5 gaps")
out("CHECK: 11/7 is the j=1 (and only) intermediate of family n=3 between 8/5 and 19/12: Delta=-139")
out("CHECK: the task's sequence 1,-1,5,13,-7153 = Delta at 2/1, 3/2, 5/3, 8/5, 19/12;")
out("       5/3 (Delta=5) and 11/7 (Delta=-139) are intermediates, not convergents;")
out("       1/1 (Delta=-1), 2/1, 3/2, 8/5, 19/12 are convergents.")

out("factorizations of |Delta| along convergents and intermediates with L <= 100, plus 485/306 (gp factor):")
fact_rows = []
for (L, K, kind) in rows:
    if L <= 100 or (K, L) == (485, 306):
        D = delta_table[(K, L)]
        f = gp("print(factor(%d))" % abs(D))
        fact_rows.append((K, L, kind, D, f))
        out("  %d/%d  %s  |Delta| = %s" % (K, L, kind.split()[0], f))
out("  (each row printed as gp's factor matrix [p, e; ...])")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S3. Legendre / Fatou-Grace placement of small power gaps (exact integers)")
out("=" * 78)
out("For each L <= 400 take every K with |2^K - 3^L| < 3^L / L (at most one K for")
out("L >= 4) and classify the reduced fraction of K/L: convergent, intermediate")
out("(with j and whether j is extreme), or neither.  Three tiers of the hypothesis:")
out("  tier A: |Delta| <= 3^L/(4L)  -> PROVED to force a convergent (Legendre)")
out("  tier B: |Delta| <= 3^L/(2L)  -> PROVED to force |alpha-K/L| < 1/L^2, hence")
out("          convergent or extreme intermediate (Fatou-Grace, CITED)")
out("  tier C: |Delta| <  3^L/L     -> the task's hypothesis; only |alpha-K/L| <")
out("          2.9/L^2 follows, which is NOT enough for Fatou-Grace.")

conv_set = {}
for n, (p, q) in enumerate(conv):
    conv_set[(p, q)] = n
inter_all = intermediates(cf_gp, conv, 30)
inter_set = {}
for (K, L, n, j, a_next) in inter_all:
    inter_set[(K, L)] = (n, j, a_next)


def classify(K, L):
    g = gcd(K, L)
    Kr, Lr = K // g, L // g
    if (Kr, Lr) in conv_set:
        return "convergent n=%d" % conv_set[(Kr, Lr)] + ("" if g == 1 else " (x%d)" % g)
    if (Kr, Lr) in inter_set:
        n, j, a_next = inter_set[(Kr, Lr)]
        ext = (j == 1 or j == a_next - 1)
        return "intermediate n=%d j=%d/%d%s" % (n, j, a_next - 1, " extreme" if ext else " INTERIOR") + \
               ("" if g == 1 else " (x%d)" % g)
    return "NEITHER"


LMAX_TIER = 400
hits = []
tierC_not_conv_or_ext = []
tierB_fail = []
tierA_fail = []
for L in range(1, LMAX_TIER + 1):
    P = 3 ** L
    Kc = round(L * float(alpha_dec))
    for K in range(max(1, Kc - 2), Kc + 3):
        D = 2 ** K - P
        if abs(D) * L < P:
            cls = classify(K, L)
            tA = 4 * L * abs(D) <= P
            tB = 2 * L * abs(D) <= P
            hits.append((L, K, D, tA, tB, cls))
            if tA and not cls.startswith("convergent"):
                tierA_fail.append((K, L, cls))
            if tB and not (cls.startswith("convergent") or "extreme" in cls):
                tierB_fail.append((K, L, cls))
            if not (cls.startswith("convergent") or "extreme" in cls):
                tierC_not_conv_or_ext.append((K, L, D, cls))
check(not tierA_fail, "tier A witness: %s" % tierA_fail)
check(not tierB_fail, "tier B witness: %s" % tierB_fail)
out("hits with L <= 60 (K/L, Delta, tierA, tierB, class):")
for (L, K, D, tA, tB, cls) in hits:
    if L <= 60:
        out("  %3d/%-3d Delta=%-22d A=%d B=%d  %s" % (K, L, D, tA, tB, cls))
out("total tier-C hits for L <= %d: %d" % (LMAX_TIER, len(hits)))
out("tier-A hits: %d, all convergents (checked)" % sum(1 for h in hits if h[3]))
out("tier-B hits: %d, all convergents or extreme intermediates (checked)" % sum(1 for h in hits if h[4]))
out("tier-C hits that are NOT convergents and NOT extreme intermediates:")
if tierC_not_conv_or_ext:
    for (K, L, D, cls) in tierC_not_conv_or_ext:
        out("   %d/%d Delta=%d  %s" % (K, L, D, cls))
else:
    out("   none for L <= %d" % LMAX_TIER)
n_neither = sum(1 for h in hits if h[5] == "NEITHER")
out("tier-C hits classified NEITHER (neither convergent nor any intermediate): %d" % n_neither)

# the elementary inequality |log2(1+e)| <= |e|/((1-|e|) ln 2) checked numerically
worst = 0.0
for i in range(1, 2500):
    e = i / 10000.0
    for s in (1, -1):
        v = abs(math.log2(1 + s * e))
        bound = e / ((1 - e) * math.log(2))
        check(v <= bound, "log inequality at e=%g" % (s * e))
        worst = max(worst, v / bound)
out("elementary inequality |log_2(1+e)| <= |e|/((1-|e|) ln 2) checked on a grid |e| < 0.25; max ratio %.4f" % worst)

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S4. Eliahou bridge: cycle min-element N forces the clock onto a convergent")
out("=" * 78)
out("For a signed odd cycle at b=+-1 with all |n_i| >= N:  Product(3n_i+b)=2^K Product n_i")
out("gives 2^K/3^L = Product(1 + b/(3 n_i)), so |Delta|/3^L <= (1+1/(3N))^L - 1 (sign +)")
out("or <= 1 - (1-1/(3N))^L (sign -).  Tier A follows once L^2 < roughly 3N/4.")


def cycle_word(nodes, b):
    """actual halving word of a signed odd cycle; raises if not a cycle"""
    w = []
    L = len(nodes)
    for i in range(L):
        m = 3 * nodes[i] + b
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        check(m == nodes[(i + 1) % L], "not a cycle: %s at %d" % (nodes, i))
        w.append(k)
    return w


def BKL(w):
    L = len(w)
    B = 0
    Ki = 0
    for i in range(L):
        B += 3 ** (L - 1 - i) * 2 ** Ki
        Ki += w[i]
    return B, Ki, L


known = [((1,), 1), ((-1,), 1), ((-5, -7), 1), ((-17, -25, -37, -55, -41, -61, -91), 1),
         ((1,), -1), ((5, 7), -1), ((17, 25, 37, 55, 41, 61, 91), -1), ((-1,), -1)]
out("known cycles: nodes, b, word, (K,L), Delta, B, q, N=min|n|, |Delta|/3^L, Eliahou bound, tierA, tierB, class")
for nodes, b in known:
    w = cycle_word(list(nodes), b)
    B, K, L = BKL(w)
    D = 2 ** K - 3 ** L
    q = abs(D) // gcd(B, abs(D))
    check(q == 1, "known cycle has q != 1")
    check(Fraction(b * B, D) == nodes[0], "fixed point mismatch")
    N = min(abs(v) for v in nodes)
    rel = Fraction(abs(D), 3 ** L)
    if nodes[0] * b > 0:
        bnd = (1 + Fraction(1, 3 * N)) ** L - 1
    else:
        bnd = 1 - (1 - Fraction(1, 3 * N)) ** L
    check(rel <= bnd, "Eliahou bound fails on %s" % (nodes,))
    tA = 4 * L * abs(D) <= 3 ** L
    tB = 2 * L * abs(D) <= 3 ** L
    out("  %s b=%d w=%s (K,L)=(%d,%d) Delta=%d B=%d q=%d N=%d rel=%.4f bound=%.4f A=%d B=%d %s"
        % (nodes, b, w, K, L, D, B, q, N, float(rel), float(bnd), tA, tB, classify(K, L)))

out("Legendre threshold: smallest N such that (1+1/(3N))^L - 1 <= 1/(4L), for L = convergent denominators:")
for n in range(1, 12):
    L = conv[n][1]
    # solve (1+1/(3N))^L <= 1 + 1/(4L): N >= 1/(3((1+1/(4L))^(1/L) - 1))
    Nmin = 1.0 / (3.0 * ((1.0 + 1.0 / (4.0 * L)) ** (1.0 / L) - 1.0))
    out("  L=%-9d N_min ~ %.4e  (about 4L^2/3 = %.4e)" % (L, Nmin, 4 * L * L / 3.0))
N68 = 2 ** 68
# largest L with exp(L/(3N)) - 1 <= 1/(4L) (exp(Lx) >= (1+x)^L, so this is a
# rigorous sufficient condition up to float rounding; a 1e-6 relative margin is kept)
Lstar = int(math.isqrt(3 * N68 // 4))
while math.expm1(Lstar / (3.0 * N68)) * (1 + 1e-6) > 1.0 / (4.0 * Lstar):
    Lstar -= 1
check(math.expm1((Lstar + 1) / (3.0 * N68)) * (1 + 1e-6) > 1.0 / (4.0 * (Lstar + 1)), "Lstar maximal")
out("With N >= 2^68 (Barina 2020 verification, CITED) tier A holds for every L <= %d (about sqrt(3*2^68/4) = %d);" % (Lstar, int(math.isqrt(3 * N68 // 4))))
out("  every positive 3n+1 cycle of odd length L <= %d has K/L (reduced) a convergent of log_2 3 (Eliahou 1993 mechanism)." % Lstar)

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S5. Cycle-gate census on convergent / intermediate clocks with L <= 13")
out("=" * 78)
out("For every (K,L), L<=13, whose reduced K/L is a convergent or intermediate")
out("fraction, all C(K-1,L-1) compositions w of K into L positive parts are")
out("enumerated; B(w) = sum 3^(L-1-i) 2^(K_i), q = |Delta|/gcd(B,|Delta|).")
out("q=1 means an exact integer odd cycle of 3n+1 (Delta>0) or 3n-1 (Delta<0).")


def census(K, L):
    """returns (count, n_q1, min_q, hist of q for q<=20, list of q=1 words, B stats)"""
    D = abs(2 ** K - 3 ** L)
    pow3 = [3 ** e for e in range(L)]
    total = 0
    n_q1 = 0
    min_q = None
    hist = {}
    q1_words = []
    # iterative DFS: w = (k_1..k_L), sum K; partial B and K_i
    w = [1] * L
    # enumerate compositions via recursion with explicit stack
    def rec(i, Ki, B, remaining):
        nonlocal total, n_q1, min_q
        # place k_i for position i (0-based); B already includes terms 0..i
        if i == L - 1:
            # last part is forced: remaining
            w[i] = remaining
            q = D // gcd(B, D)
            total += 1
            if q == 1:
                n_q1 += 1
                if len(q1_words) < 40:
                    q1_words.append((tuple(w), B))
            if min_q is None or q < min_q:
                min_q = q
            if q <= 20:
                hist[q] = hist.get(q, 0) + 1
            return
        parts_left = L - 1 - i
        for k in range(1, remaining - parts_left + 1):
            w[i] = k
            rec(i + 1, Ki + k, B + pow3[L - 2 - i] * 2 ** (Ki + k), remaining - k)
    rec(0, 0, pow3[L - 1], K)
    check(total == comb(K - 1, L - 1), "composition count")
    return total, n_q1, min_q, hist, q1_words


clocks = []
for L in range(1, 14):
    for K in range(L, 2 * L + 1):
        cls = classify(K, L)
        if cls != "NEITHER":
            clocks.append((L, K, cls))
out("clock  class                     Delta      #words   #q=1  min_q  q-histogram (q<=20)")
grand_q1 = {1: [], -1: []}
for (L, K, cls) in clocks:
    D = 2 ** K - 3 ** L
    total, n_q1, min_q, hist, q1w = census(K, L)
    hs = " ".join("%d:%d" % (q, hist[q]) for q in sorted(hist))
    out("  %2d/%-2d %-26s %-10d %-8d %-5d %-6d %s" % (K, L, cls, D, total, n_q1, min_q, hs))
    if n_q1:
        sign = 1 if D > 0 else -1
        for (w, B) in q1w:
            grand_q1[sign].append((K, L, w, B))
out("q=1 words (exact integer cycles), by sign of Delta:")
for sign in (1, -1):
    out("  3n%s1 (Delta %s 0):" % ("+" if sign > 0 else "-", ">" if sign > 0 else "<"))
    for (K, L, w, B) in grand_q1[sign]:
        D = 2 ** K - 3 ** L
        n0 = Fraction(sign * B, D)   # b = sign so that b/Delta > 0 ... see note
        out("    (K,L)=(%d,%d) w=%s B=%d n_0=bB/Delta with b=%d: %s" % (K, L, w, B, sign, n0))
n_pos = len(grand_q1[1])
n_neg = len(grand_q1[-1])
out("count of q=1 words: 3n+1 side %d, 3n-1 side %d" % (n_pos, n_neg))
# expectations from the inherited classification: 3n+1: (2,1) once, (4,2) once, (6,3) once, ... (2d,d) all-2 words (repeats of {1});
# 3n-1: (d,d) all-1 words (repeats of {1}); (3,2) two words, (6,4) [rotations/repeats of (1,2)]... ; (11,7) seven rotations
exp_pos = sorted((2 * d, d) for d in range(1, 14) if classify(2 * d, d) != "NEITHER")
got_pos = sorted((K, L) for (K, L, w, B) in grand_q1[1])
check(got_pos == exp_pos, "3n+1 q=1 clocks: got %s expected %s" % (got_pos, exp_pos))
for (K, L, w, B) in grand_q1[1]:
    check(all(k == 2 for k in w), "3n+1 q=1 word not all-2: %s" % (w,))
out("CHECK: every 3n+1 q=1 word with L<=13 on these clocks is (2,2,...,2) = repeated fixed point {1}: PASS")
neg_clocks = {}
for (K, L, w, B) in grand_q1[-1]:
    neg_clocks.setdefault((K, L), []).append(w)
out("3n-1 q=1 words grouped by clock (count):", {k: len(v) for k, v in sorted(neg_clocks.items())})
# (d,d): 1 word; (3d,2d): repeats/rotations of (1,2): number of words of length 2d over {(1,2),(2,1)} pattern = 2 for each d? check by reconstruction
for (K, L), ws in neg_clocks.items():
    for w in ws:
        # reconstruct the cycle nodes from the word and verify with the actual map
        B, K2, L2 = BKL(list(w))
        D = 2 ** K2 - 3 ** L2
        n0 = Fraction(-B, D)
        check(n0.denominator == 1, "not integral")
        n = int(n0)
        nodes = [n]
        for k in w:
            m = 3 * nodes[-1] - 1
            check(m % 2 ** k == 0 and (m // 2 ** k) % 2 == 1, "valuation mismatch")
            nodes.append(m // 2 ** k)
        check(nodes[-1] == nodes[0], "cycle does not close")
        least = set(nodes[:-1])
        check(least <= {1} or least <= {5, 7} or least <= {17, 25, 37, 55, 41, 61, 91},
              "unexpected 3n-1 cycle %s" % (nodes,))
out("CHECK: every 3n-1 q=1 word with L<=13 on these clocks reconstructs to {1}, {5,7} or the 17-cycle: PASS")
check(len(neg_clocks.get((11, 7), [])) == 7, "seven-cycle should give exactly 7 words on 11/7")
out("CHECK: clock 11/7 carries exactly 7 q=1 words = the 7 rotations of (1,1,1,2,1,1,4): PASS")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S6. Full-clock minimum-q census, all (K,L) with L<=11, L<=K<=2L")
out("=" * 78)
out("For each L: min over K and over all words of q, the clock attaining it, and")
out("whether that clock is a convergent/intermediate; plus the number of q=1 words")
out("per L (all clocks) and the number of clocks with min_q = |Delta| (no cancellation).")
FULL_L = 11
tot_words = 0
for L in range(1, FULL_L + 1):
    best = None
    q1_total = 0
    no_cancel = 0
    per_clock = []
    for K in range(L, 2 * L + 1):
        D = 2 ** K - 3 ** L
        total, n_q1, min_q, hist, q1w = census(K, L)
        tot_words += total
        q1_total += n_q1
        if min_q == abs(D):
            no_cancel += 1
        per_clock.append((K, min_q, n_q1))
        if best is None or min_q < best[0]:
            best = (min_q, K, n_q1)
    # second-best distinct clock by min_q for context
    per_clock_sorted = sorted(per_clock, key=lambda t: t[1])
    out("  L=%2d  words=%-8d q=1 words=%-3d  min_q=%d at K=%d (%s)  clocks with no cancellation: %d/%d  three smallest min_q: %s"
        % (L, sum(comb(K - 1, L - 1) for K in range(L, 2 * L + 1)), q1_total, best[0], best[1],
           classify(best[1], L), no_cancel, L + 1,
           ", ".join("K=%d:q=%d" % (K, mq) for (K, mq, nq) in per_clock_sorted[:3])))
out("total words enumerated in S6: %d = sum_{L<=%d} C(2L,L)" % (tot_words, FULL_L))
check(tot_words == sum(comb(2 * L, L) for L in range(1, FULL_L + 1)), "word total")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S7. Pillai-type census |2^a - 10^b| <= 100 (complete, with proof of completeness)")
out("=" * 78)
sols = []
for a in range(0, 201):
    for b in range(0, 62):
        d = 2 ** a - 10 ** b
        if abs(d) <= 100:
            sols.append((a, b, d))
out("all (a,b) with 0<=a<=200, 0<=b<=61 and |2^a-10^b|<=100:")
for (a, b, d) in sols:
    out("  2^%d - 10^%d = %d" % (a, b, d))
out("count:", len(sols))
# completeness proof check: for a >= b >= 1, 2^a-10^b = 2^b (2^(a-b) - 5^b) so |.| >= 2^b unless 2^(a-b)=5^b (impossible)
# hence b <= 6; for a < b, 10^b - 2^a > 10^b - 2^b > 100 for b >= 3.
maxb = max(b for (a, b, d) in sols)
check(maxb <= 6, "b bound")
check((10, 3, 24) in sols, "2^10-10^3=24")
out("2^10 = 10^3 + 24 = 10^3 + 2^3*3; the residual 24 is 2^3*3 = 23+1 = 4! = 5^2-1.")
out("largest b among solutions: %d (proof: 2^a-10^b = 2^b(2^(a-b)-5^b) forces |.| >= 2^b > 100 for b >= 7)" % maxb)
# also the honest Collatz-relevant family: |2^K - 3^L| <= 100
out("for contrast, all (K,L) with 1<=K<=200, 1<=L<=130 and |2^K - 3^L| <= 100 (coprime bases; NOT finite by an elementary 2-adic trick):")
kl = []
for K in range(1, 201):
    for L in range(1, 131):
        d = 2 ** K - 3 ** L
        if abs(d) <= 100:
            kl.append((K, L, d))
for (K, L, d) in kl:
    out("  2^%d - 3^%d = %d   [%s]" % (K, L, d, classify(K, L)))
out("count:", len(kl))
off_tree = [(K, L, d) for (K, L, d) in kl if classify(K, L) == "NEITHER"]
out("clocks in this list that are NEITHER convergent nor intermediate (reduced): %d" % len(off_tree))
for (K, L, d) in off_tree:
    out("  2^%d - 3^%d = %d : |Delta|*L = %d >= 3^L = %d, so outside tier C (the absolute bound 100 is not small relative to 3^L here)"
        % (K, L, d, abs(d) * L, 3 ** L))
    check(abs(d) * L >= 3 ** L, "off-tree clock inside tier C: %d/%d" % (K, L))
in_tierC = [(K, L, d) for (K, L, d) in kl if abs(d) * L < 3 ** L]
check(all(classify(K, L) != "NEITHER" for (K, L, d) in in_tierC), "tier-C clock off the convergent tree")
out("CHECK: every |2^K-3^L| <= 100 clock in the box with |Delta|*L < 3^L (%d of %d) is a convergent or intermediate fraction of log_2 3: PASS"
    % (len(in_tierC), len(kl)))
out("  (the absolute-size condition |Delta|<=100 and the relative-size condition |Delta|<3^L/L differ at small L;")
out("   only the relative one is a Legendre-type hypothesis.)")
out("  (Pillai 1931 / de Weger: |2^K-3^L| grows; the complete list of |2^K-3^L|<=100 is CITED context,")
out("   the box above is FINITE-EXACT only.)")

# ---------------------------------------------------------------------------
out()
out("=" * 78)
out("S8. Same-clock carry obstruction on 19/12 and 5/3 (how close q comes to 1)")
out("=" * 78)
for (K, L) in ((5, 3), (19, 12), (11, 7), (8, 5)):
    D = 2 ** K - 3 ** L
    total, n_q1, min_q, hist, q1w = census(K, L)
    out("  (K,L)=(%d,%d) Delta=%d words=%d min_q=%d #q=1=%d #q<=20=%d" % (K, L, D, total, min_q, n_q1, sum(hist.values())))
    if (K, L) == (5, 3):
        # list all six words with B and q
        for w in ((1, 1, 3), (1, 2, 2), (1, 3, 1), (2, 1, 2), (2, 2, 1), (3, 1, 1)):
            B, _, _ = BKL(list(w))
            q = abs(D) // gcd(B, abs(D))
            out("     w=%s B=%d q=%d  -> cycle at b=5 iff q|5: %s" % (w, B, q, q in (1, 5)))

print("elapsed %.1f s" % (time.time() - T0), file=sys.stderr)
out("DONE")
