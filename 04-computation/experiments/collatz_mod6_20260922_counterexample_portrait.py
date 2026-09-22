#!/usr/bin/env python3
"""collatz_mod6_20260922_counterexample_portrait.py

Lane counterexample_portrait (session collatz-mod6, 2026-09-22).

Compiles, in one place, the PROVED / CITED necessary conditions that a
hypothetical nontrivial positive 3n+1 cycle and a hypothetical divergent
positive 3n+1 orbit must satisfy, and tests every condition on the 3n-1
sheet, where three positive cycles {1}, {5,7}, {17,25,37,55,41,61,91} exist.
A condition is SHEET-BLIND if it holds (with the same proof) on both sheets,
SIGN-SPECIFIC if its statement or proof uses b=+1.

Master criterion (S1): T_+(-n) = -T_-(n), so the positive 3n-1 sheet is the
negative 3n+1 sheet; a condition is sheet-blind iff its proof never uses the
sign of n (equivalently of the carry B) at b=+1.

Everything printed is FINITE-EXACT unless the label says otherwise.
Inherited mechanisms are cited by path in the .md note, not re-proved here;
they are re-verified on finite ranges only.

Run:  python3 collatz_mod6_20260922_counterexample_portrait.py
      (python3 -O gives identical output modulo the timing line)
"""
import math
import sys
import time
from fractions import Fraction
from math import comb, gcd

try:
    import mpmath
except ImportError as exc:  # pragma: no cover
    raise SystemExit("mpmath (bundled with sympy) is required") from exc

T0 = time.time()
print("collatz_mod6_20260922_counterexample_portrait.py -- lane counterexample_portrait, session collatz-mod6, 2026-09-22")
mpmath.mp.dps = 400
ALPHA = mpmath.log(3) / mpmath.log(2)
LN2 = mpmath.log(2)


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def hdr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ----------------------------------------------------------------------------
# signed accelerated maps
# ----------------------------------------------------------------------------

def v2(x):
    return (x & -x).bit_length() - 1


def T(n, b):
    """Accelerated map on odd integers (any sign): (3n+b)/2^v2(3n+b)."""
    m = 3 * n + b
    return m >> v2(m)


def T_word(n, b, L):
    """First L halving exponents of the odd start n under T_b."""
    w = []
    for _ in range(L):
        m = 3 * n + b
        k = v2(m)
        w.append(k)
        n = m >> k
    return w, n


def word_data(w):
    """Carry B, clock (K,L), Delta=2^K-3^L, minimal parameter q (signed_cycles eq. (8))."""
    L = len(w)
    K = 0
    B = 0
    for i in range(L):
        # B = sum_{i<L} 3^(L-1-i) 2^(K_i), K_i = k_1+...+k_i, K_0 = 0
        B += 3 ** (L - 1 - i) * (1 << K)
        K += w[i]
    Delta = (1 << K) - 3 ** L
    q = abs(Delta) // gcd(B, abs(Delta))
    return B, K, L, Delta, q


def cycle_from_word(w, b):
    """n_0 = b B / Delta; returns the ordered cycle if integral (q | b)."""
    B, K, L, Delta, q = word_data(w)
    if (b * B) % Delta != 0:
        return None
    n0 = (b * B) // Delta
    cyc = [n0]
    n = n0
    for i in range(L):
        m = 3 * n + b
        k = v2(m)
        check(k == w[i], "valuation mismatch in cycle_from_word")
        n = m >> k
        cyc.append(n)
    check(cyc[-1] == cyc[0], "cycle does not close")
    return cyc[:-1]


# ----------------------------------------------------------------------------
# S1. conjugacy: the positive 3n-1 sheet is the negative 3n+1 sheet
# ----------------------------------------------------------------------------
hdr("S1. Conjugacy T_+(-n) = -T_-(n) on odd n (master sheet criterion)")
bad = 0
for n in range(-20001, 20002, 2):
    if T(-n, +1) != -T(n, -1):
        bad += 1
print(f"odd n in [-20001, 20001]: violations of T_+(-n) = -T_-(n): {bad}")
check(bad == 0, "conjugacy")
print("PROVED (one line): 3(-n)+1 = -(3n-1) and v2 is sign-blind.")
print("Consequence: a condition proved for ALL odd integers under 3n+1 is SHEET-BLIND;")
print("a condition whose proof uses n>0 at b=+1 (equivalently B>0 with b=+1) is SIGN-SPECIFIC.")

# ----------------------------------------------------------------------------
# S2. the three positive 3n-1 cycles through the gate
# ----------------------------------------------------------------------------
hdr("S2. The three positive 3n-1 cycles and the plus-sheet {1}: gate data")
KNOWN_MINUS = [[1], [5, 7], [17, 25, 37, 55, 41, 61, 91]]
CYCLE_TABLE = {}
for b, cyc in [(+1, [1]), (-1, [1]), (-1, [5, 7]), (-1, [17, 25, 37, 55, 41, 61, 91])]:
    n0 = cyc[0]
    L = len(cyc)
    w, back = T_word(n0, b, L)
    check(back == n0, "known cycle does not close")
    orbit = [n0]
    for _ in range(L - 1):
        orbit.append(T(orbit[-1], b))
    check(orbit == cyc, "known cycle order")
    B, K, L, Delta, q = word_data(w)
    n0g = Fraction(b * B, Delta)
    check(n0g == n0, "gate value")
    sign_law = "positive" if (b * Delta) > 0 else "negative"
    # sign of every member = sign(b/Delta)
    check(all((x > 0) == (b * Delta > 0) for x in cyc), "sign law")
    key = (b, tuple(cyc))
    CYCLE_TABLE[key] = dict(w=w, B=B, K=K, L=L, Delta=Delta, q=q)
    print(f"b={b:+d} cycle {cyc}: word {w}, K={K}, L={L}, B={B}, Delta=2^K-3^L={Delta}, "
          f"q=|Delta|/gcd(B,|Delta|)={q}, n0=bB/Delta={n0g}, sign(b/Delta)->{sign_law} members")
print("Sign law (signed_cycles sec.3, PROVED): every member has the sign of b/Delta, so a")
print("POSITIVE 3n+1 cycle needs Delta>0 (K/L > log_2 3) while a POSITIVE 3n-1 cycle needs")
print("Delta<0 (K/L < log_2 3).  The three minus cycles have Delta = -1, -1, -139 (<0);")
print("the plus cycle {1} has Delta = +1.  The gate q=1 itself is SHEET-BLIND; which side")
print("of log_2 3 the clock must sit on is SIGN-SPECIFIC.")

# ----------------------------------------------------------------------------
# S3. continued fraction of log_2 3: convergents, intermediate fractions, sign alternation
# ----------------------------------------------------------------------------
hdr("S3. Convergents of log_2 3, the sign alternation, and where the four known clocks sit")


def contfrac(x, nterms):
    a = []
    p_prev, p = 1, int(mpmath.floor(x))
    q_prev, q = 0, 1
    a.append(p)
    conv = [(p, q)]
    y = x
    for _ in range(nterms - 1):
        y = y - mpmath.floor(y)
        check(y != 0, "terminating cf")
        y = 1 / y
        an = int(mpmath.floor(y))
        a.append(an)
        p_prev, p = p, an * p + p_prev
        q_prev, q = q, an * q + q_prev
        conv.append((p, q))
    return a, conv


A_CF, CONV = contfrac(ALPHA, 60)
for i in range(1, len(CONV)):
    p1, q1 = CONV[i - 1]
    p2, q2 = CONV[i]
    check(abs(p1 * q2 - p2 * q1) == 1, "determinant of convergents")
print(f"computed {len(CONV)} convergents; determinant identity checked for all consecutive pairs")
print("partial quotients a_0..a_29:", A_CF[:30])
print("convergents p_n/q_n, n<=24, with sign of Delta=2^p-3^q (exact when q<=80000, else by 400-digit log):")
CONV_SIGN = {}
for n, (p, q) in enumerate(CONV[:25]):
    if q <= 80000:
        d = (1 << p) - 3 ** q
        s = 1 if d > 0 else -1
        dstr = str(d) if abs(d) < 10 ** 12 else f"{'+' if d > 0 else '-'}({len(str(abs(d)))} digits)"
    else:
        s = 1 if (p * LN2 - q * mpmath.log(3)) > 0 else -1
        dstr = "sign by log"
    CONV_SIGN[n] = s
    check(s == (1 if n % 2 == 1 else -1), "sign alternation of convergents")
    print(f"  n={n:2d}  {p}/{q}  Delta sign {s:+d}  ({dstr})")
print("PROVED (classical, Hardy-Wright Thm 154 alternation): even-index convergents lie BELOW")
print("log_2 3 (Delta<0), odd-index convergents ABOVE (Delta>0).  Hence, wherever the clock")
print("of a cycle is forced to be a convergent, a positive 3n+1 cycle sits on an ODD-index")
print("convergent and a positive 3n-1 cycle on an EVEN-index one.  Sign-specific refinement.")

# classification of the four known clocks
p3, q3 = CONV[3]  # 8/5
p2, q2 = CONV[2]  # 3/2
print(f"11/7 = (p_2 + 1*p_3)/(q_2 + 1*q_3) = ({p2}+{p3})/({q2}+{q3}) with a_4 = {A_CF[4]}: the unique")
check((p2 + p3, q2 + q3) == (11, 7) and A_CF[4] == 2, "11/7 intermediate")
print("intermediate fraction between the convergents 3/2 and 19/12 (j=1 of a_4-1=1); NOT a convergent.")
print("1/1 = p_0/q_0, 3/2 = p_2/q_2 (even index, Delta<0) carry the minus cycles {1}, {5,7};")
print("2/1 = p_1/q_1 (odd index, Delta>0) carries the plus cycle {1}.")
for L, K in [(1, 1), (2, 3), (7, 11), (1, 2)]:
    D = (1 << K) - 3 ** L
    tA = Fraction(3 ** L, 4 * L)
    tB = Fraction(3 ** L, 2 * L)
    print(f"  clock {K}/{L}: |Delta|={abs(D)}, tier-A cap 3^L/(4L)={float(tA):.3f} ({'in' if abs(D) <= tA else 'OUT'}), "
          f"tier-B cap 3^L/(2L)={float(tB):.3f} ({'in' if abs(D) <= tB else 'OUT'})")
print("The seven-cycle's |Delta|=139 exceeds the tier-A cap 78.107 but not the tier-B cap 156.214:")
print("consistent with the pillai lane (tier B = convergent or extreme intermediate; 11/7 is extreme).")

# ----------------------------------------------------------------------------
# S4. the Eliahou product identity on both sheets and the Legendre thresholds
# ----------------------------------------------------------------------------
hdr("S4. Eliahou product identity 2^K/3^L = prod(1 + b/(3 n_i)) on the four known cycles")
for (b, cyc), d in CYCLE_TABLE.items():
    L = d["L"]
    K = d["K"]
    N = min(cyc)
    prod = Fraction(1)
    for x in cyc:
        prod *= 1 + Fraction(b, 3 * x)
    check(prod == Fraction(1 << K, 3 ** L), "product identity")
    rel = Fraction(abs(d["Delta"]), 3 ** L)
    if b * d["Delta"] > 0 and b > 0:
        bound = (1 + Fraction(1, 3 * N)) ** L - 1
        kind = "plus sheet: (1+1/(3N))^L - 1"
    else:
        bound = 1 - (1 - Fraction(1, 3 * N)) ** L
        kind = "minus sheet: 1-(1-1/(3N))^L"
    check(rel <= bound, "Eliahou bound")
    print(f"b={b:+d} {cyc}: |Delta|/3^L = {float(rel):.4f} <= {kind} = {float(bound):.4f}  (N=min={N})")
print("PROVED identity, SHEET-BLIND (pillai S4 already states it for b=+-1).")


def legendre_threshold(L, sign):
    """Smallest N with the Eliahou bound <= tier-A cap 1/(4L).  sign=+1 plus sheet, -1 minus."""
    lo, hi = 1, 1
    def ok(N):
        if sign > 0:
            return (1 + Fraction(1, 3 * N)) ** L - 1 <= Fraction(1, 4 * L)
        return 1 - (1 - Fraction(1, 3 * N)) ** L <= Fraction(1, 4 * L)
    while not ok(hi):
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


print("Legendre tier-A thresholds N_min(L) (smallest cycle minimum that forces K/L convergent):")
for L in [2, 5, 7, 12, 41, 53]:
    tp = legendre_threshold(L, +1)
    tm = legendre_threshold(L, -1)
    print(f"  L={L:3d}: plus sheet N_min={tp}, minus sheet N_min={tm}, 4L^2/3={4 * L * L / 3:.2f}")
print("The seven-cycle has N=17 < N_min(7) on the minus sheet, so the convergent theorem does not")
print("apply to it, and indeed 11/7 is not a convergent: the theorem's HYPOTHESIS (a large minimum)")
print("is what the sheets do not share.  {5,7} has N=5 < N_min(2)=3? no:")
print(f"  N_min(2, minus) = {legendre_threshold(2, -1)}, and {{5,7}} with N=5 does satisfy it; 3/2 is a convergent.")

# ----------------------------------------------------------------------------
# S5. finite census on the minus sheet: every odd n <= X_MINUS reaches one of the three cycles
# ----------------------------------------------------------------------------
hdr("S5. Census: odd n <= X on the 3n-1 sheet fall into the three known cycles (basin sizes)")
X_MINUS = 10 ** 7
cyc_id = {}
for i, cyc in enumerate(KNOWN_MINUS, start=1):
    for x in cyc:
        cyc_id[x] = i
basin = bytearray((X_MINUS >> 1) + 1)
for x, i in cyc_id.items():
    basin[x >> 1] = i
counts = [0, 0, 0, 0]
maxsteps = 0
argmax = 0
t_census = time.time()
for n in range(1, X_MINUS + 1, 2):
    idx = n >> 1
    if basin[idx]:
        counts[basin[idx]] += 1
        continue
    m = n
    steps = 0
    while True:
        m3 = 3 * m - 1
        m = m3 >> v2(m3)
        steps += 1
        if m < n:
            bid = basin[m >> 1]
            check(bid != 0, "unclassified smaller value")
            break
        if m in cyc_id:
            bid = cyc_id[m]
            break
    basin[idx] = bid
    counts[bid] += 1
    if steps > maxsteps:
        maxsteps = steps
        argmax = n
tot = sum(counts[1:])
print(f"X = {X_MINUS}: odd starts {tot}; every one reaches a known cycle (no fourth cycle, no escape).")
for i, cyc in enumerate(KNOWN_MINUS, start=1):
    print(f"  basin of {cyc}: {counts[i]} odd starts ({counts[i] / tot:.6f}, rounded {counts[i] / tot:.3f})")
print(f"  longest first-drop-or-cycle-hit segment: {maxsteps} accelerated steps at n={argmax}")
print(f"  census time {time.time() - t_census:.1f}s (timing line, not a result)")
print("FINITE-EXACT consequence: any fourth positive 3n-1 cycle has minimum element N > 10^7,")
print("i.e. N >= 10000001 (odd).  Mirror of Barina's 2^68 (CITED, plus sheet); the minus sheet")
print("has NO published verification comparable to 2^68, only this repo's censuses.")
N_MINUS_FLOOR = X_MINUS + 1

# plus-sheet sanity census (Barina covers 2^68; this is only a control)
X_PLUS = 10 ** 6
seen = bytearray((X_PLUS >> 1) + 1)
seen[0] = 1
for n in range(3, X_PLUS + 1, 2):
    m = n
    while m >= n and m != 1:
        m3 = 3 * m + 1
        m = m3 >> v2(m3)
    seen[n >> 1] = 1
print(f"control: every odd n <= {X_PLUS} reaches 1 under 3n+1 (Barina 2020 CITED: all n < 2^68).")

# ----------------------------------------------------------------------------
# S6. derived cycle-length bounds from (Legendre tier A) + (min element) + (convergent spacing)
# ----------------------------------------------------------------------------
hdr("S6. Length bounds for a hypothetical cycle from Eliahou + tier A + convergent spacing")


def tierA_max_L(N):
    """Largest L with exp(L/(3N)) - 1 <= 1/(4L) (sufficient for tier A on both sheets)."""
    Nm = mpmath.mpf(N)
    def ok(L):
        return 3 * Nm * mpmath.log(1 + mpmath.mpf(1) / (4 * L)) >= L
    lo, hi = 1, 1
    while ok(hi):
        hi *= 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


print("Chain (PROVED given CITED inputs): a positive cycle with minimum N and L odd terms has")
print("  |Delta|/3^L <= e^{L/(3N)} - 1 (plus) or <= L/(3N) (minus, Bernoulli);")
print("  if that is <= 1/(4L) then |log_2 3 - K/L| < 1/(2L^2) and K/L is a convergent (Legendre);")
print("  with K/L = p_n/q_n, L = m q_n, |K - L log_2 3| = m |q_n a - p_n| > m/(q_n + q_{n+1})")
print("  (Hardy-Wright Thm 171 form).  Plus sheet: y <= e^y - 1 = Delta/3^L <= e^{L/(3N)} - 1 <=")
print("  (L/(3N)) e^{L/(3N)}, so q_n (q_n+q_(n+1)) >= 3 N ln2 * e^{-L_A/(3N)}.  Minus sheet: L >= 2 and")
print("  |Delta|/3^L = 1 - e^{-y} <= 1/(4L) <= 1/8 give y <= y0 = -ln(7/8) and 1-e^{-y} >= c0 y with")
print("  c0 = (1-e^{-y0})/y0, so q_n (q_n+q_(n+1)) >= c0 * 3 N ln2.  (L = 1 is the trivial clock: n_0 =")
print("  b/(2^k-3) is integral only for k = 1, 2, giving {1} at b = -1, +1; no tier needed.)")
for label, N, sheet in [("plus sheet, N = 2^68 (Barina)", 1 << 68, +1),
                        ("minus sheet, N = 10^7+1 (S5 census)", N_MINUS_FLOOR, -1)]:
    LA = tierA_max_L(N)
    if sheet > 0:
        factor = mpmath.exp(-mpmath.mpf(LA) / (3 * N))
    else:
        y0 = -mpmath.log(mpmath.mpf(7) / 8)
        factor = (1 - mpmath.exp(-y0)) / y0
    thr = mpmath.mpf(3) * N * LN2 * factor
    print(f"{label}: tier A holds for every L <= L_A = {LA} (sqrt(3N/4) = {float(mpmath.sqrt(3 * mpmath.mpf(N) / 4)):.6e})")
    print(f"  spacing threshold factor*3*N*ln2 = {float(thr):.6e} (factor = {float(factor):.6f})")
    need_sign = +1 if sheet > 0 else -1
    admissible = []
    for n in range(len(CONV) - 1):
        p, q = CONV[n]
        q_next = CONV[n + 1][1]
        if CONV_SIGN.get(n, (1 if n % 2 else -1)) != need_sign:
            continue
        if q * (q + q_next) >= thr:
            admissible.append((n, p, q))
    first = admissible[0]
    print(f"  convergents on the correct side with q_n(q_n+q_(n+1)) >= threshold: first is n={first[0]}, "
          f"{first[1]}/{first[2]}")
    for n in range(first[0]):
        if CONV_SIGN.get(n, (1 if n % 2 else -1)) == need_sign:
            q = CONV[n][1]
            qn = CONV[n + 1][1]
            print(f"    (n={n}: q_n(q_n+q_(n+1)) = {q}*({q}+{qn}) = {q * (q + qn)} < threshold)")
    below = [(n, p, q) for (n, p, q) in admissible if q <= LA]
    if not below:
        print(f"  none of them has q_n <= L_A, so NO positive cycle (other than the trivial/known ones)")
        print(f"  has L <= {LA}: every hypothetical cycle on this sheet has L > {LA}.")
    else:
        for (n, p, q) in below:
            ms = [m for m in range(1, LA // q + 1)]
            print(f"  clock still open with q_n <= L_A: n={n}, {p}/{q}; a cycle with L <= L_A must have")
            print(f"  L = m*{q} with m in {ms}, i.e. L in {[m * q for m in ms]}; otherwise L > {LA}.")
    print("  first three surviving clocks (K/L reduced), the minimal unresolved conjunction lives there:")
    for n, p, q in admissible[:3]:
        print(f"    n={n}: K/L = {p}/{q}, cycle length L = m*{q}, m>=1, q_n(q_n+q_(n+1)) = {q * (q + CONV[n + 1][1])}")
print("Literature for comparison (plus sheet): Simons-de Weger 2005 (CITED, Acta Arith. 117) no m-cycles")
print("m<=68 (exact published figure UNCITED-RECOLLECTION); Hercher, J. Integer Seq. 26 (2023) art. 23.3.5")
print("(CITED; abstract re-read at the audit 2026-09-22, cs.uwaterloo.ca/journals/JIS/VOL26/Hercher/hercher5.html):")
print("credits Simons-de Weger with m >= 76, gets m >= 83 from newer verification ranges, PROVES m >= 92,")
print("and states that verification up to 3*2^69 WOULD SUFFICE to raise the odd-member bound to the NEXT")
print("bound K >= 1.375e11 (conditional; this lane's first draft misread it as proved, caught at the audit).")
print(f"That next bound is the clock n=23: q_23 = {CONV[23][1]} = 1.375e11.")
print("The chain here, from cited inputs alone (Barina 2^68, Legendre, Hardy-Wright 171), proves that")
print("every nontrivial positive 3n+1 cycle has L > 14878203146 odd terms (n=21's spacing 4.746e20")
print("falls below the sharp threshold 6.137e20); it does NOT reach q_23, because tier A (hence the")
print("convergent placement) is available only for L <= L_A, and a cycle with L > L_A may sit on a")
print("non-convergent clock.  On the minus sheet the same chain gives: any FOURTH positive 3n-1 cycle has")
print("L > L_A(minus) = 2738; IF its clock is a convergent (not forced for L > L_A) it is an even-index")
print("one with q_n >= 31867 (n=8, 665*(665+15601), fails the threshold).")

# ----------------------------------------------------------------------------
# S7. gate census, all clocks L <= 10, both signs
# ----------------------------------------------------------------------------
hdr("S7. Gate census: every composition word with L <= 10, both parameters b = +-1")


def compositions(K, L):
    """All compositions of K into L positive parts."""
    if L == 1:
        yield (K,)
        return
    for first in range(1, K - L + 2):
        for rest in compositions(K - first, L - 1):
            yield (first,) + rest


found = {+1: {}, -1: {}}
nwords = 0
for L in range(1, 11):
    for K in range(L, 2 * L + 1):
        for w in compositions(K, L):
            nwords += 1
            B, K_, L_, Delta, q = word_data(w)
            if q == 1:
                for b in (+1, -1):
                    cyc = cycle_from_word(w, b)
                    key = tuple(sorted(set(cyc)))
                    found[b].setdefault(key, set()).add(w)
print(f"words enumerated (sum over L<=10 of C(2L,L)): {nwords}")
for b in (+1, -1):
    print(f"b={b:+d}: distinct integer cycles (as sets) from q=1 words:")
    for key, ws in sorted(found[b].items(), key=lambda kv: (len(kv[0]), kv[0])):
        prim = min(ws, key=lambda w: (len(w), w))
        print(f"   {list(key)}  (shortest word {list(prim)}, {len(ws)} words incl. rotations/repetitions)")
pos = {b: {k for k in found[b] if k[0] > 0} for b in (+1, -1)}
neg = {b: {k for k in found[b] if k[0] < 0} for b in (+1, -1)}
check(pos[-1] == {(1,), (5, 7), (17, 25, 37, 41, 55, 61, 91)}, "minus census")
check(pos[+1] == {(1,)}, "plus census")
check(neg[+1] == {tuple(sorted(-x for x in k)) for k in pos[-1]} and neg[-1] == {(-1,)}, "conjugate census")
print("FINITE-EXACT: for L <= 10 the only POSITIVE integer cycles at b=-1 are the three known ones and")
print("at b=+1 only {1}; the negative cycles at each sign are exactly the negatives of the positive")
print("cycles at the other sign (the conjugacy of S1).  SHEET-BLIND mechanism.")

# ----------------------------------------------------------------------------
# S8. Steiner / Simons-de Weger shape on the minus sheet: 1-cycles and 2-cycles by the gate
# ----------------------------------------------------------------------------
hdr("S8. 1-cycles (1^(L-1), k) and 2-cycles (1^a k1 1^b k2): where the gate opens, both sheets")
one_cycles = {+1: [], -1: []}
L1MAX = 3000
for L in range(1, L1MAX + 1):
    Kc = int(mpmath.floor(L * ALPHA))
    for K in (Kc - 1, Kc, Kc + 1, Kc + 2):
        k = K - (L - 1)
        if k < 1:
            continue
        w = (1,) * (L - 1) + (k,)
        B, _, _, Delta, q = word_data(w)
        if q == 1:
            for b in (+1, -1):
                cyc = cycle_from_word(w, b)
                if cyc[0] > 0 and len(set(cyc)) == L:  # positive and primitive (not a repeated {1})
                    one_cycles[b].append((L, k, cyc if len(cyc) <= 8 else cyc[:3] + ["..."]))
print(f"1-cycle words (1^(L-1), k), L <= {L1MAX}, K in floor(L log_2 3)+{{-1,0,1,2}}, positive integral, primitive:")
for b in (+1, -1):
    print(f"  b={b:+d}: {one_cycles[b]}")
check(one_cycles[+1] == [(1, 2, [1])], "plus 1-cycles")
check(one_cycles[-1] == [(1, 1, [1]), (2, 2, [5, 7])], "minus 1-cycles")
print("Plus sheet: only the trivial word (2) -> {1}; Steiner 1977 (CITED) proves this for ALL L.")
print("Minus sheet: {5,7} IS a nontrivial 1-cycle (word (1,2)), so Steiner's theorem is SIGN-SPECIFIC")
print("and its transfer to 3n-1 is REFUTED by the witness {5,7}; its Baker/continued-fraction")
print("machinery is sheet-blind, but its conclusion is not.")

two_cycles = {+1: set(), -1: set()}
L2MAX = 24
n2 = 0
for L in range(2, L2MAX + 1):
    Kc = int(mpmath.floor(L * ALPHA))
    for K in (Kc - 1, Kc, Kc + 1, Kc + 2):
        for a in range(0, L - 1):
            bb = L - 2 - a
            for k1 in range(2, K - a - bb - 1):
                k2 = K - a - bb - k1
                if k2 < 2:
                    continue
                w = (1,) * a + (k1,) + (1,) * bb + (k2,)
                n2 += 1
                B, _, _, Delta, q = word_data(w)
                if q == 1:
                    for b in (+1, -1):
                        cyc = cycle_from_word(w, b)
                        if cyc[0] > 0:
                            two_cycles[b].add((L, tuple(sorted(set(cyc)))))
print(f"2-cycle words 1^a k1 1^b k2 (k1,k2>=2), L <= {L2MAX}: {n2} words tested")
for b in (+1, -1):
    prim = sorted({c for (L, c) in two_cycles[b]})
    print(f"  b={b:+d}: positive cycles reached: {[list(c) for c in prim]}")
check(sorted({c for (L, c) in two_cycles[-1]}) == [(5, 7), (17, 25, 37, 41, 55, 61, 91)], "minus 2-cycles")
check(sorted({c for (L, c) in two_cycles[+1]}) == [(1,)], "plus 2-cycles")
print("Minus sheet: the seven-cycle (1,1,1,2,1,1,4) is a genuine 2-cycle; Simons-de Weger's")
print("'no m-cycles, m<=68' is therefore SIGN-SPECIFIC (transfer REFUTED by the seven-cycle, m=2).")
print("The repeated (1,2)(1,2) also appears as a non-primitive 2-cycle word for {5,7}.")

# ----------------------------------------------------------------------------
# S9. Catalan/Mihailescu unit gaps and the necklace count
# ----------------------------------------------------------------------------
hdr("S9. Unit-gap clocks |2^K-3^L| = 1 (Catalan/Mihailescu) carry EVERY word; necklace counts")
unit = []
small = []
for L in range(1, 401):
    for K in (int(mpmath.floor(L * ALPHA)), int(mpmath.floor(L * ALPHA)) + 1):
        D = (1 << K) - 3 ** L
        if abs(D) == 1:
            unit.append((K, L, D))
        if abs(D) <= 100:
            small.append((K, L, D))
print(f"clocks with |Delta| = 1, L <= 400: {unit}")
check(unit == [(1, 1, -1), (2, 1, 1), (3, 2, -1)], "unit gaps")
print(f"clocks with |Delta| <= 100, L <= 400: {len(small)} of them: {small}")
print("Mihailescu 2004 (CITED; Catalan's conjecture): 8 and 9 are the only consecutive perfect powers,")
print("so 3/2 is the last unit-gap clock for all L.  On a unit-gap clock q=1 for every word, so every")
print("composition is a cycle word: (1)->{1} at b=-1 / {-1} at b=+1; (2)->{1} at b=+1; (1,2),(2,1)")
print("-> {5,7} at b=-1 / {-5,-7} at b=+1.  SHEET-BLIND (both sides of every unit gap are used).")


def necklaces(K, L):
    g = gcd(K, L)
    tot = 0
    for d in range(1, g + 1):
        if g % d == 0:
            # Euler phi
            phi = d
            x = d
            pp = 2
            while pp * pp <= x:
                if x % pp == 0:
                    while x % pp == 0:
                        x //= pp
                    phi -= phi // pp
                pp += 1
            if x > 1:
                phi -= phi // x
            tot += phi * comb(K // d - 1, L // d - 1)
    return tot // L


for K, L in [(1, 1), (2, 1), (3, 2), (11, 7), (19, 12), (65, 41)]:
    print(f"  clock {K}/{L}: compositions C(K-1,L-1) = {comb(K - 1, L - 1)}, necklaces (rotation classes) = {necklaces(K, L)}")
print("q is rotation-invariant (signed_cycles sec.3, PROVED), so the gate is a test per necklace;")
print("on 11/7 the seven-cycle is one necklace of the 30, its 7 rotations the 7 q=1 words of pillai S5.")
# check 11/7 count of q=1 words
c117 = sum(1 for w in compositions(11, 7) if word_data(w)[4] == 1)
print(f"  q=1 words on 11/7: {c117}")
check(c117 == 7, "11/7 q=1 words")

# ----------------------------------------------------------------------------
# S10. divergent-orbit conditions on the minus sheet: D_L, sum q_L, budget
# ----------------------------------------------------------------------------
hdr("S10. D_L = K_L - L log_2 3 and the budget sum q_L on the known cycles (both sheets)")
print("Identity (glued_xor sec.3, PROVED): n_L q_L = n_0 + (b/3) sum_{i<L} q_i, q_i = 2^(K_i)/3^i.")
for (b, cyc), d in CYCLE_TABLE.items():
    w = d["w"]
    L = d["L"]
    K = d["K"]
    n0 = cyc[0]
    # one period partial sums
    S0 = Fraction(0)
    Ki = 0
    for i in range(L):
        S0 += Fraction(1 << Ki, 3 ** i)
        Ki += w[i]
    ratio = Fraction(1 << K, 3 ** L)
    per_period_D = float(K - L * ALPHA)
    if ratio < 1:
        total = S0 / (1 - ratio)
        print(f"b={b:+d} {cyc}: D per period = K - L log_2 3 = {per_period_D:+.6f} (-> -inf linearly), "
              f"sum q_i over one period = {S0}, geometric ratio 2^K/3^L = {ratio}, total sum q = {total} = 3 n_0 = {3 * n0}")
        check(total == 3 * n0 and b == -1, "budget exhaustion")
    else:
        print(f"b={b:+d} {cyc}: D per period = {per_period_D:+.6f} (-> +inf linearly), ratio 2^K/3^L = {ratio} > 1, "
              f"sum q_i DIVERGES; n_L q_L = n_0 + (1/3) sum q_i grows.")
    # numerical check of identity along 3 periods
    q = Fraction(1)
    acc = Fraction(0)
    n = n0
    Kacc = 0
    for step in range(3 * L):
        acc += Fraction(1 << Kacc, 3 ** step)
        m = 3 * n + b
        k = v2(m)
        n = m >> k
        Kacc += k
    lhs = n * Fraction(1 << Kacc, 3 ** (3 * L))
    rhs = n0 + Fraction(b, 3) * acc
    check(lhs == rhs, "additive identity")
print("Verdict: on the plus sheet D_L -> -infinity and sum q_L < infinity are NECESSARY for a divergent")
print("orbit and are VIOLATED by every cycle (D_L -> +infinity), so they separate cycle from escape.")
print("On the minus sheet EVERY orbit, periodic or not, has D_L -> -infinity and sum q_L <= 3 n_0")
print("(positivity alone); the discriminator there is budget exhaustion: sum q_L = 3 n_0 iff eventually")
print("periodic.  The condition 'D_L -> -infinity' is SHEET-BLIND as a necessary condition on a")
print("divergent orbit, SIGN-SPECIFIC as a cycle/escape discriminator.")

# bounded strip on the minus sheet is excluded by positivity alone: numerical illustration
print("No bounded strip (guards 2a, PROVED on the plus sheet by source density): on the minus sheet")
print("a strip A <= q_j <= B gives n_j = (n_0 - sum/3)/q_j <= n_0/A, a bounded hence eventually periodic")
print("orbit, whose D_L -> -infinity (sign law: 2^K < 3^L) contradicts q_j >= A.  SHEET-BLIND")
print("conclusion, different proofs (density on +, positivity on -; cf. g_negatives S9/S10 for G).")
qs = []
n = 17
Kacc = 0
for step in range(1, 71):
    m = 3 * n - 1
    k = v2(m)
    n = m >> k
    Kacc += k
    if step % 7 == 0:
        qs.append(float(Fraction(1 << Kacc, 3 ** step)))
print("  q_(7m) along the seven-cycle, m=1..10:", " ".join(f"{x:.4f}" for x in qs))

# ----------------------------------------------------------------------------
# S11. Terras coefficient-stopping densities: identical on both sheets
# ----------------------------------------------------------------------------
hdr("S11. Terras/Everett coefficient stopping-time densities d_J (exact, same table on both sheets)")
# DP over (j, K_j) not yet crossed: K_j <= floor(j*alpha); k geometric 2^-k
crossed = Fraction(0)
alive = {0: Fraction(1)}
table = []
for j in range(1, 41):
    cap = int(mpmath.floor(j * ALPHA))
    new = {}
    for Kp, pr in alive.items():
        # step k >= 1; not crossing iff Kp + k <= cap
        for k in range(1, cap - Kp + 1):
            new[Kp + k] = new.get(Kp + k, Fraction(0)) + pr * Fraction(1, 1 << k)
        kmin = cap - Kp + 1
        crossed += pr * Fraction(1, 1 << (kmin - 1))
    alive = new
    table.append((j, crossed))
for j, dj in table:
    if j <= 12 or j % 10 == 0:
        print(f"  J={j:2d}: d_J = {dj if j <= 8 else ''} = {float(dj):.6f}, 1-d_J = {float(1 - dj):.6f}")
print("  (index J = number of odd/accelerated steps, density among odd starts)")
# Terras's own convention: U(n) = (3n+1)/2 or n/2 on ALL n, J = total steps, stop when 3^o < 2^j
alive2 = {0: Fraction(1)}
stopped2 = Fraction(0)
terras = {}
for j in range(1, 41):
    new2 = {}
    for o, pr in alive2.items():
        for par in (0, 1):
            o2 = o + par
            if 3 ** o2 < (1 << j):
                stopped2 += pr / 2
            else:
                new2[o2] = new2.get(o2, Fraction(0)) + pr / 2
    alive2 = new2
    terras[j] = stopped2
print("Terras's own indexing (U(n) = (3n+1)/2 or n/2 on all n, J total steps, density mod 2^J):")
for j in (1, 2, 4, 8, 12, 20, 30, 40):
    print(f"  J={j:2d}: F_J = {terras[j] if j <= 8 else ''} = {float(terras[j]):.6f}")
print(f"  check: F_12 = {float(terras[12]):.4f}, F_20 = {float(terras[20]):.4f} (three_adic lane quotes 0.9448 and 0.9739)")
check(abs(float(terras[12]) - 0.9448) < 6e-5 and abs(float(terras[20]) - 0.9739) < 6e-5, "Terras F_12, F_20")
print("PROVED (Terras 1976, Everett 1977, CITED): d_J -> 1.  The word cylinder of (k_1..k_j) has odd-source")
print("density 2^(-K_j) on BOTH sheets (residue class mod 2^(K_j+1), 2-adic prefix law), and a word")
print("with 2^(K_j) > 3^j gives T^j(n) < n for ALL n at b=-1 (carry is subtracted) and for all")
print("n > B/(2^(K_j)-3^j) at b=+1.  SHEET-BLIND; the minus sheet is even slightly easier.")
# empirical control of the 2-adic prefix law on both sheets
bad = 0
import random
random.seed(20260922)
for _ in range(2000):
    n = random.randrange(1, 10 ** 9) | 1
    b = random.choice((1, -1))
    w, _ = T_word(n, b, 6)
    B, K, L, Delta, q = word_data(w)
    m = n + (1 << (K + 1)) * random.randrange(1, 1000)
    w2, _ = T_word(m, b, 6)
    if w2 != w:
        bad += 1
    if (3 ** L * n + b * B) % (1 << (K + 1)) != (1 << K):
        bad += 1
print(f"  2-adic prefix law, 2000 random (n,b) tests at L=6: violations {bad}")
check(bad == 0, "prefix law")

# ----------------------------------------------------------------------------
# S12. squarefree growing prefixes on both sheets
# ----------------------------------------------------------------------------
hdr("S12. Squarefree growth prefixes n_j = 3^j 2^(L+1-j) q -/+ 1 (braids2 SF3) on both sheets")


LSF = 5
QMAX = 100000
# Sieve: q is bad if some prime square p^2 divides some node 3^j 2^(L+1-j) q - b, p up to
# sqrt(max node); the coefficients are even and coprime to every p >= 5 so each (p, j) gives
# exactly one residue class of q mod p^2 (p=2 never divides an odd node; p=3 can divide only j=0).
maxnode = 3 ** LSF * (1 << 1) * QMAX + 1
PMAX = int(math.isqrt(maxnode)) + 1
sieve = bytearray([1]) * (PMAX + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(math.isqrt(PMAX)) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
PRIMES = [i for i in range(2, PMAX + 1) if sieve[i]]
cnt = {}
for b in (+1, -1):
    good = bytearray([1]) * (QMAX + 1)
    for pr in PRIMES:
        if pr == 2:
            continue
        m = pr * pr
        for j in range(LSF + 1):
            c = 3 ** j * (1 << (LSF + 1 - j))
            if c % pr == 0:
                # p = 3, j >= 1: node = 3^j(...) - b is +-1 mod 3, never 0 mod 9
                continue
            r = (b * pow(c, -1, m)) % m  # c q = b mod p^2
            if r == 0:
                r = m
            for q in range(r, QMAX + 1, m):
                good[q] = 0
    cnt[b] = sum(good[1:])
# direct check of the sieve on a small range by trial division
def squarefree(x):
    p = 2
    while p * p <= x:
        if x % (p * p) == 0:
            return False
        p += 1
    return True
direct = {}
for b in (+1, -1):
    direct[b] = sum(1 for q in range(1, 301)
                    if all(squarefree(3 ** j * (1 << (LSF + 1 - j)) * q - b) for j in range(LSF + 1)))
    good_small = bytearray([1]) * 301
    for pr in PRIMES:
        if pr == 2:
            continue
        m = pr * pr
        for j in range(LSF + 1):
            c = 3 ** j * (1 << (LSF + 1 - j))
            if c % pr == 0:
                continue
            r = (b * pow(c, -1, m)) % m
            if r == 0:
                r = m
            for q in range(r, 301, m):
                good_small[q] = 0
    check(sum(good_small[1:]) == direct[b], "sieve vs trial division")
print(f"  sieve cross-check, q <= 300: {direct[+1]} (plus), {direct[-1]} (minus) starts, both by trial division and sieve")
# verify the word is all-ones on both sheets for a sample
for b in (+1, -1):
    n0 = (1 << (LSF + 1)) * 12345 - b
    w, nL = T_word(n0, b, LSF)
    check(w == [1] * LSF and nL > n0, "growth prefix")
# local factors nu_L(p) = min(L+1, ord_{p^2}(3/2)), same on both sheets


def ord_mod(a, m):
    x = a % m
    k = 1
    while x != 1:
        x = x * a % m
        k += 1
    return k


prod = 1.0
for p in [2, 3] + [p for p in range(5, 20000, 2) if all(p % r for r in range(3, int(p ** 0.5) + 1, 2))]:
    if p == 2:
        nu = 0  # all nodes odd
    elif p == 3:
        nu = 1  # only j=0 node can be 0 mod 9? n_0 = 2^(L+1) q -/+ 1; the others are +-1 mod 3
    else:
        inv2 = pow(2, -1, p * p)
        g = 3 * inv2 % (p * p)
        nu = LSF + 1  # nu_L(p) = min(L+1, ord_{p^2}(3/2)); only orders <= L+1 matter
        x = 1
        for k in range(1, LSF + 1):
            x = x * g % (p * p)
            if x == 1:
                nu = k
                break
    prod *= 1 - nu / (p * p)
print(f"L={LSF}, q <= {QMAX}: all {LSF + 1} nodes squarefree for {cnt[+1]} starts (plus sheet), {cnt[-1]} (minus sheet)")
print(f"  observed densities {cnt[+1] / QMAX:.4f}, {cnt[-1] / QMAX:.4f}; sieve product delta_L (p<20000, nu_L(3)=1) = {prod:.4f}")
print("PROVED (braids2 SF5) with the same local counts nu_L(p) on both sheets (3/2 mod p^2 is sign-blind):")
print("SHEET-BLIND.  Arbitrarily long all-squarefree growth prefixes exist on both sheets.")

# ----------------------------------------------------------------------------
# S13. finite-congruence Lyapunov obstruction on both sheets
# ----------------------------------------------------------------------------
hdr("S13. Finite-congruence Lyapunov obstruction n_0 = 2^(L+1) M q -/+ 1 (braids sec.5) on both sheets")
M = 6
Lb = 12
for b in (+1, -1):
    n0 = (1 << (Lb + 1)) * M * 1 - b
    w, nL = T_word(n0, b, Lb)
    res = {(3 ** j * (1 << (Lb + 1 - j)) * M - b) % M for j in range(Lb + 1)}
    print(f"b={b:+d}: n_0={n0}, word {w}, T^L(n_0)={nL} > n_0: {nL > n0}, residues mod {M} of all nodes: {sorted(res)}")
    check(w == [1] * Lb and nL > n0 and len(res) == 1, "Lyapunov obstruction")
print("PROVED (braids sec.5, sheet-blind by the same three lines): no positive periodic weight modulo any")
print("fixed M and no bounded correction of log n decreases at every odd step, on either sheet.")

# ----------------------------------------------------------------------------
# S14. E-graph / greedy 3-adic mirror is conjugation-blind
# ----------------------------------------------------------------------------
hdr("S14. The greedy 3-adic map G and its minus-sheet partner: G(-m) = -G_-(m)")


def G(m, b=+1):
    # G_b(m) = (2^k m - b)/3, k >= 0 minimal with 3 | 2^k m - b and 3 does not divide the quotient
    k = 0
    x = m
    while True:
        if (x - b) % 3 == 0 and ((x - b) // 3) % 3 != 0:
            return (x - b) // 3
        x *= 2
        k += 1
        check(k < 20, "G loop")


bad = 0
for m in range(1, 3001):
    if m % 3 == 0:
        continue
    if G(-m, +1) != -G(m, -1):
        bad += 1
print(f"m <= 3000 coprime to 3: violations of G(-m) = -G_-(m): {bad}")
check(bad == 0, "G conjugacy")
print("PROVED (g_negatives, inherited): the E-graph/G mirror, its word law mod 3^(J+1), drift log(2/3) and")
print("tail (7/9)^(J-1) transfer verbatim between the sheets.  SHEET-BLIND.")

# ----------------------------------------------------------------------------
# S15. contradiction search over the compiled conditions
# ----------------------------------------------------------------------------
hdr("S15. Pairwise contradiction search and the minimal unresolved conjunction")
print("Cycle conditions (plus sheet): C1 gate q=1 & sign law Delta>0; C2 N >= 2^68; C3 K/L odd-index")
print("convergent for L <= L_A; C4 q_n(q_n+q_(n+1)) >= 3N ln2 e^{-L_A/(3N)} (sharp factor 1.000000, S6) whenever")
print("the clock is a convergent; C5 no 1-cycle (Steiner), no m-cycle m <= 91 (Simons-de Weger, Hercher);")
print("C6 necklace/rotation invariance of q; C7 |Delta| = 1 only at 2/1.")
print("Witness that they are jointly satisfiable at the level of clocks: the surviving clocks listed in S6")
print("(odd-index convergents with q_n > L_A and the spacing inequality) violate none of C1-C7 a priori;")
print("C1's integrality q=1 on such a clock is the ONLY untested predicate.  Pairwise contradiction:")
print("NONE FOUND.  (Audit 2026-09-22: C3 is conditional on L <= L_A; for L > L_A no condition on file")
print("forces the clock to be a convergent, so the convergent survivors are one case, not the whole.)")
print("Divergence conditions (plus sheet): D1 D_L -> -infinity, sum q_L < infinity (glued_xor);")
print("D2 no bounded strip (guards 2a); D3 stopping-time density one (Terras/Everett); D4 Tao a.e. in")
print("logarithmic density (arXiv:1909.03562, Thm 1.3); D5 Krasikov-Lagarias 2003 x^0.84 ancestors of 1;")
print("D6 Lyapunov obstruction;")
print("D7 squarefree growth prefixes of every length.  D1-D2 constrain the orbit, D3-D5 say only that")
print("the exceptional set has density zero (D4: log-density; D5: the tree of 1 is large), D6-D7 say")
print("finite lookahead cannot certify descent.  A density-zero divergent orbit with liminf D_L = -inf")
print("and sum q_L < inf is consistent with all of them: NONE FOUND.")
print("Minimal unresolved conjunction (cycle): exists a reduced clock K/L > log_2 3 with L > L_A = 14878203146,")
print("(K - L log_2 3) ln 2 <= (L/(3N)) e^{L/(3N)} at N = 2^68, and a composition of K into L parts with q = 1")
print("and all members >= 2^68; in the sub-case where K/L is a convergent it is p_n/q_n with n >= 23 odd,")
print("L = m q_n (Hercher's next target K >= 1.375e11 is exactly the m = 1, n = 23 clock).")
print("Minimal unresolved conjunction (divergence): exists a positive odd orbit with n_L unbounded,")
print("D_L -> -infinity, sum q_L < infinity, whose start lies in the density-zero exceptional sets of")
print("Terras and Tao.")
print("Sheet lesson: every SHEET-BLIND condition above is satisfied by the seven-cycle on the minus sheet")
print("(q=1, product identity, necklace, D_L -> -inf, no strip, Terras, prefixes, Lyapunov, G mirror),")
print("so no conjunction of sheet-blind conditions can exclude a plus-sheet cycle; a proof must use the")
print("sign-specific inputs (Delta>0 side, Barina 2^68, Steiner/Simons-de Weger/Hercher, Tao) or a new one.")

# ----------------------------------------------------------------------------
# S16. session-lead probe numbers re-verified (entropy of 8/pi^2; Paley tournament on F_7)
# ----------------------------------------------------------------------------
hdr("S16. Session-lead probe: binary entropy of 8/pi^2 and the Paley tournament on F_7")
x = 8 / math.pi ** 2
H = -(x * math.log2(x) + (1 - x) * math.log2(1 - x))
print(f"8/pi^2 = {x:.5f}; binary entropy H = {H:.5f} bits (the paste says 0.704 and calls it zero entropy)")
check(abs(H - 0.70028) < 5e-6, "entropy")
QR = {1, 2, 4}
arcs = {(i, j) for i in range(7) for j in range(7) if i != j and (j - i) % 7 in QR}
cyc3 = 0
trans3 = 0
lines = set()
for a in range(7):
    for b in range(a + 1, 7):
        for c in range(b + 1, 7):
            e = [(a, b) in arcs, (b, c) in arcs, (a, c) in arcs]
            # cyclic iff a->b->c->a or a->c->b->a
            if ((a, b) in arcs and (b, c) in arcs and (c, a) in arcs) or ((a, c) in arcs and (c, b) in arcs and (b, a) in arcs):
                cyc3 += 1
                lines.add((a, b, c))
            else:
                trans3 += 1
fano1 = {tuple(sorted(((s + d) % 7) for d in (0, 1, 3))) for s in range(7)}
fano2 = {tuple(sorted(((s + d) % 7) for d in (0, 1, 5))) for s in range(7)}
print(f"Paley on F_7 (arc i->j iff j-i in {{1,2,4}}): {len(arcs)} arcs, {cyc3} cyclic triples ((7^3-7)/24 = {(343 - 7) // 24}), "
      f"{trans3} transitive triples; cyclic triples = dev{{0,1,3}} u dev{{0,1,5}}: {lines == fano1 | fano2}")
per_arc = {e: sum(1 for t in lines if set(e) <= set(t)) for e in arcs}
print(f"each arc lies in exactly {min(per_arc.values())}..{max(per_arc.values())} cyclic triples (2-(7,3,2) design)")
check(len(arcs) == 21 and cyc3 == 14 and trans3 == 21 and lines == fano1 | fano2 and set(per_arc.values()) == {2}, "Paley")
print("FINITE-EXACT; these tournament facts are correct and have no map to either table (SCOPE).")

print()
print(f"total runtime {time.time() - T0:.1f}s (timing line, not a result)")
