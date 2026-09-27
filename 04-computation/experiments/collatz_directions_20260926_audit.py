#!/usr/bin/env python3
"""collatz_directions_20260926_audit.py -- independent adversarial audit of the directions note
05-knowledge/results/collatz_directions_20260926.md and its script collatz_directions_20260926.py.

Written blind to the note's script (own DP, own scans, own identities). Every check prints PASS/FAIL and a
counter is kept; the final line reports the totals. Sections:

 (A1) W_k by an own no-descent DP (strict 3^o > 2^j as in THM-4495; the note's >= is equivalent since 3^o = 2^j
      is impossible for j >= 1), cross-checked against THM-4495's published W_1..W_20 and against its Spitzer
      identity k W_k = sum B_n W_(k-n) for k <= 300; the box-dimension digits of the note at k = 60, 200, 1000,
      2000 (plain and with the +1.5 log_2 k / k ballot correction); h* = h(log_3 2).
 (A2) Proposition 2: every no-descent T-word u of length k <= 12 has fixed point x_u = S_u/(2^k - 3^p) that is
      negative, 2-adically integral, has parity word u^infinity (checked by 2-adic simulation on rationals) and
      all heights >= 0; E_inf cap Z_<0 to 10^6 by the HEIGHT criterion (independent of the note's actual-descent
      criterion); first coefficient-descent time versus first actual-descent time for every negative odd
      |n| <= 10^6 (the note's "actual and coefficient descent coincide" claim); the rotations of the three
      negative cycles; why -7, -25 are excluded; |U(m)| < |m| iff v >= 2.
 (A3) Proposition 3: Terras equality sigma = sigma_inf checked for 2 <= n <= 10^6 (n = 1 is the known exception,
      THM-4512); the note's "exactly" needs it.
 (A4) Proposition 4: heights of (1,1,2)^N in both codings; the stricter witness (1,1,1,1,2,3)^N whose periods
      contain excursions returning strictly below their own start while staying above the source.
 (B)  Section 2: the lsb-map identities; every n <= 10^5 reaches a power of two; maximal odd-step count and its
      argmax with the total stopping time cross-check; n = 27. Section 5: best lower approximations of log_2 3 by
      exact integer comparison for p <= 700, the continued fraction of log_2 3 (Decimal, 200 digits) and which
      listed fractions are convergents; exact gaps; 8/5 as an upper approximation; 3^p - 2^A = +-1 enumeration
      (Catalan check); 139 | S for -17. Section 3: the identities 2^(d_l) m_l/3^l = n C_l and
      S_(l-1)/3^l = n (C_l - 1) exactly to l = 60 for 27, -1, -5 and 200 random n of both signs; the exact real
      value of the Bernstein series on the three cycles (= -n) and the 2-adic convergence v_2(S_(l-1) + 3^l n) = d_l.
Usage: python3 collatz_directions_20260926_audit.py > ../../05-knowledge/results/collatz_directions_20260926_audit.out
"""
import hashlib
import math
import os
import random
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction

T0 = time.time()
PASSES = 0
FAILS = 0


def check(name, ok, detail=""):
    global PASSES, FAILS
    if ok:
        PASSES += 1
        print("  PASS  %s%s" % (name, ("  [" + detail + "]") if detail else ""))
    else:
        FAILS += 1
        print("  FAIL  %s%s" % (name, ("  [" + detail + "]") if detail else ""))


def elapsed():
    return "%.1fs" % (time.time() - T0)


LOG23 = math.log(3) / math.log(2)
LOG32 = math.log(2) / math.log(3)


def binent(p):
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


HSTAR = binent(LOG32)

# ----------------------------------------------------------------------------------------------------------------
# powers, exact comparisons
KMAX = 2000
POW2 = [1 << j for j in range(KMAX + 4000)]
POW3 = [3 ** o for o in range(KMAX + 4000)]


def no_descent_ok(o, j):
    """THM-4495's condition 3^o > 2^j (strict)."""
    return POW3[o] > POW2[j]


# ================================================================================================================
print("=== (A1) W_k by an own DP, THM-4495 cross-checks, box-dimension digits ===")
# equality 3^o = 2^j impossible for j >= 1 (parity); so >= and > agree for j >= 1
eq_found = any(POW3[o] == POW2[j] for j in range(1, 3001) for o in range(0, 3001) if POW3[o] <= POW2[j] and POW3[o] >= POW2[j])
check("3^o = 2^j never holds for 1 <= j <= 3000 (so the note's >= equals THM-4495's strict >)", not eq_found)

# own DP: state = number of odd letters o after j letters; keep only o with 3^o > 2^j
omin = [0] * (KMAX + 1)
o = 0
for j in range(1, KMAX + 1):
    while POW3[o] <= POW2[j]:
        o += 1
    omin[j] = o  # minimal o with 3^o > 2^j
cnt = {0: 1}
W = [1]
for j in range(1, KMAX + 1):
    nxt = {}
    for oo, c in cnt.items():
        # letter 0: o stays; letter 1: o + 1
        if oo >= omin[j]:
            nxt[oo] = nxt.get(oo, 0) + c
        if oo + 1 >= omin[j]:
            nxt[oo + 1] = nxt.get(oo + 1, 0) + c
    cnt = nxt
    W.append(sum(cnt.values()))
print("  W_1..W_20 (own DP): %s" % W[1:21])
THM4495_W = [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328]
check("W_1..W_20 equal THM-4495's published list", W[1:21] == THM4495_W)


def binom_tail(n):
    return sum(math.comb(n, j) for j in range(0, n + 1) if POW3[j] > POW2[n])


B = [0] + [binom_tail(n) for n in range(1, 301)]
spitzer_ok = all(k * W[k] == sum(B[n] * W[k - n] for n in range(1, k + 1)) for k in range(1, 301))
check("Spitzer identity k W_k = sum_n B_n W_(k-n) holds for the own W_k, k <= 300 (independent of the DP)", spitzer_ok)
print("  h* = h(log_3 2) = %.10f  (note: 0.9499555 / 0.94996)" % HSTAR)
check("h* rounds to 0.9499555 and 0.94996", abs(HSTAR - 0.9499555) < 5e-8 and abs(HSTAR - 0.94996) < 5e-6)
NOTE_PLAIN = {60: 0.8496, 200: 0.9084, 1000: 0.9384, 2000: 0.9435}
NOTE_CORR = {60: 0.9973, 200: 0.9657, 1000: 0.9534, 2000: 0.9517}
NOTE_DIGITS = {60: 16, 200: 55, 1000: 283, 2000: 569}
for k in (10, 20, 41, 60, 100, 200, 500, 1000, 1500, 2000):
    lw = math.log2(W[k]) if W[k] < 2 ** 1000 else (W[k].bit_length() - 1 + math.log2(W[k] / (1 << (W[k].bit_length() - 1))))
    plain = lw / k
    corr = plain + 1.5 * math.log2(k) / k
    c_k = W[k] * k ** 1.5 / 2 ** (HSTAR * k) if k <= 1000 else float(Fraction(W[k]) / Fraction(2 ** int(HSTAR * k))) * k ** 1.5 / 2 ** (HSTAR * k - int(HSTAR * k))
    print("  k=%5d  digits=%4d  log2(W_k)/k=%.6f  +1.5log2(k)/k=%.6f  W_k k^1.5 2^(-h*k)=%.4f" % (k, len(str(W[k])), plain, corr, c_k))
    if k in NOTE_PLAIN:
        check("note's digits at k=%d: plain %.4f, corrected %.4f, %d digits" % (k, NOTE_PLAIN[k], NOTE_CORR[k], NOTE_DIGITS[k]),
              round(plain, 4) == NOTE_PLAIN[k] and round(corr, 4) == NOTE_CORR[k] and len(str(W[k])) == NOTE_DIGITS[k],
              "own: %.4f, %.4f, %d digits" % (plain, corr, len(str(W[k]))))
# W_k = Theta(2^(h k) k^-3/2): the corrected sequence converges like log_2(c_k)/k with c_k in [9.66, 11.05] (THM-4495)
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (A2) Proposition 2: fixed points of no-descent words; E_inf cap Z_<0 by the height criterion ===")


def t_step_frac(x):
    """Collatz T on a rational with odd denominator; parity = parity of the numerator."""
    if x.numerator & 1:
        return (3 * x + 1) / 2
    return x / 2


def words_no_descent(k):
    """All T-words of length k (tuples of 0/1) with 3^(o_j) > 2^j for all 1 <= j <= k."""
    out = []

    def rec(prefix, o):
        j = len(prefix)
        if j == k:
            out.append(tuple(prefix))
            return
        for letter in (0, 1):
            oo = o + letter
            if POW3[oo] > POW2[j + 1]:
                prefix.append(letter)
                rec(prefix, oo)
                prefix.pop()

    rec([], 0)
    return out


def affine_S(word):
    """T^k(x) = (3^p x + S)/2^k for the T-word `word`; returns (p, S)."""
    p = sum(word)
    S = 0
    o_so_far = 0
    for i, letter in enumerate(word):  # position i (0-based): step i+1
        if letter:
            o_so_far += 1
            S += POW3[p - o_so_far] * POW2[i]
    return p, S


allgood = True
fixed_negative = True
integral = True
word_ok = True
heights_ok = True
counts_ok = True
example_1110 = None
for k in range(1, 13):
    ws = words_no_descent(k)
    counts_ok &= (len(ws) == W[k])
    for u in ws:
        p, S = affine_S(u)
        x = Fraction(S, POW2[k] - POW3[p])
        fixed_negative &= (x < 0)
        integral &= (x.denominator % 2 == 1)
        # simulate the parity word of x for 3k letters and check periodicity u^infinity and T^k(x) = x
        y = x
        letters = []
        for _ in range(3 * k):
            letters.append(y.numerator & 1)
            y = t_step_frac(y)
        word_ok &= (tuple(letters) == u * 3)
        yk = x
        for _ in range(k):
            yk = t_step_frac(yk)
        word_ok &= (yk == x)
        # heights of u^infinity: h_(mk+i) = m h_k + h_i >= 0
        oo = 0
        for j, letter in enumerate(u * 3, start=1):
            oo += letter
            heights_ok &= (POW3[oo] > POW2[j])
        if u == (1, 1, 1, 0):
            example_1110 = x
check("number of no-descent words of length k equals own W_k for k <= 12", counts_ok)
check("fixed points x_u of all no-descent words (k <= 12) are negative", fixed_negative)
check("fixed points x_u have odd denominators (lie in Z_2)", integral)
check("parity word of x_u is u^infinity and T^k(x_u) = x_u (2-adic simulation)", word_ok)
check("u^infinity has all prefix heights >= 0 (extension argument)", heights_ok)
check("x_(1110) = -19/11 (the (1,1,2) word; S12 note agrees)", example_1110 == Fraction(-19, 11), str(example_1110))


def scan(n, want_actual=True, maxsteps=200000):
    """Follow the T-orbit of integer n. Returns (jc, ja, status) where jc = first coefficient-descent step
    (3^(o_j) < 2^j), ja = first actual-descent step (|T^j n| < |n|), status in {'cycle+', 'cycle-', 'both',
    'open'} describing how the scan ended: 'cycle+' = closed a cycle with positive period height before any
    coefficient descent (so n in E_inf), 'cycle-' = closed a cycle with negative period height."""
    x = n
    j = 0
    o = 0
    jc = None
    ja = None
    seen = {x: (0, 0)}
    absn = abs(n)
    while j < maxsteps:
        if x & 1:
            x = (3 * x + 1) >> 1
            o += 1
        else:
            x >>= 1
        j += 1
        if jc is None and POW3[o] < POW2[j]:
            jc = j
        if ja is None and abs(x) < absn:
            ja = j
        if jc is not None and (ja is not None or not want_actual):
            return jc, ja, "both"
        if x in seen:
            j1, o1 = seen[x]
            P = j - j1
            dO = o - o1
            if jc is None:
                return jc, ja, ("cycle+" if POW3[dO] > POW2[P] else "cycle-")
            return jc, ja, "cycle"
        seen[x] = (j, o)
    return jc, ja, "open"


# E_inf cap Z_<0 for |n| <= 10^6 by the height criterion
members = []
coincide_fail = []
status_count = {}
LIM = 10 ** 6
for n in range(-1, -LIM - 1, -2):
    jc, ja, st = scan(n)
    status_count[st] = status_count.get(st, 0) + 1
    if st == "cycle+":
        members.append(n)
    if jc != ja:
        coincide_fail.append((n, jc, ja, st))
print("  scan statuses (negative odd |n| <= 10^6): %s" % status_count)
check("E_inf cap Z_<0 (height criterion, |n| <= 10^6) = {-1, -5, -17}", members == [-1, -5, -17], str(members))
check("first coefficient-descent step == first actual-descent step for every negative odd |n| <= 10^6",
      not coincide_fail, "failures: %s" % coincide_fail[:5])
# even negative n: word starts with 0, h_1 = -1 < 0: never in E_inf
check("even negative n are never in E_inf (first letter 0 gives h_1 = -1)", all(scan(n)[0] == 1 for n in range(-2, -2000, -2)))
# why -7 and -25 are excluded
for n in (-7, -25, -3, -9):
    jc, ja, st = scan(n)
    x = n
    hs = []
    o = 0
    for j in range(1, (jc or 12) + 1):
        if x & 1:
            x = (3 * x + 1) >> 1
            o += 1
        else:
            x >>= 1
        hs.append(round(o * LOG23 - j, 3))
    print("  n=%d: first coefficient descent at j=%s, first actual descent at j=%s; heights %s" % (n, jc, ja, hs))
check("-7 is excluded (descent at j=2: word 1,0 -> -10 -> -5)", scan(-7)[0] == 2 and scan(-7)[1] == 2)
check("-25 is excluded (descent at j=10, landing on -17)", scan(-25)[0] == 10 and scan(-25)[1] == 10)
# rotations of the three cycles: which elements are in E_inf
for start in (-1, -5, -17):
    cyc = [start]
    x = start
    while True:
        x = (3 * x + 1) >> 1 if x & 1 else x >> 1
        if x == start:
            break
        cyc.append(x)
    inE = [c for c in cyc if scan(c)[2] == "cycle+"]
    print("  cycle of %d: %s; elements in E_inf: %s" % (start, cyc, inE))
    check("only the height-minimal rotation of the %d-cycle lies in E_inf" % start, inE == [start])
# |U(m)| < |m| iff v >= 2 for negative odd m
one_step_ok = True
for m in range(-1, -200001, -2):
    y = 3 * m + 1
    v = (y & -y).bit_length() - 1
    U = y >> v
    one_step_ok &= ((abs(U) < abs(m)) == (v >= 2))
check("|U(m)| < |m| iff v >= 2 for all negative odd |m| <= 2*10^5", one_step_ok)
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (A3) Proposition 3: Terras equality sigma = sigma_inf on positive integers ===")
terras_fail = []
for n in range(2, LIM + 1):
    jc, ja, st = scan(n)
    if jc != ja:
        terras_fail.append((n, jc, ja, st))
check("sigma(n) = sigma_inf(n) for all 2 <= n <= 10^6 (T-coded; n = 1 is the known exception)", not terras_fail, str(terras_fail[:5]))
jc1, ja1, st1 = scan(1, maxsteps=50)
print("  n = 1: first coefficient descent j=%s, first actual descent j=%s (never strictly below 1): sigma_inf(1) = 2, sigma(1) = infinity" % (jc1, ja1))
print("  NOTE: Collatz => Z^+ cap E_inf = empty; the converse needs sigma = sigma_inf for ALL n (open; THM-4512 reports it to 10^7).")
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (A4) Proposition 4: the witness words ===")


def heights_T(valword, N):
    """T-coded prefix heights of (valword)^N."""
    hs = []
    o = 0
    j = 0
    for _ in range(N):
        for v in valword:
            o += 1
            for _ in range(v):
                j += 1
                hs.append(o * LOG23 - j)
    return hs


def heights_syr(valword, N):
    hs = []
    o = 0
    d = 0
    for _ in range(N):
        for v in valword:
            o += 1
            d += v
            hs.append(o * LOG23 - d)
    return hs


h1 = heights_syr((1, 1, 2), 1)
hT = heights_T((1, 1, 2), 1)
print("  (1,1,2): Syracuse-coded heights %s; T-coded heights %s; period height %.3f" % ([round(h, 3) for h in h1], [round(h, 3) for h in hT], 3 * LOG23 - 4))
check("(1,1,2) prefix heights are 0.585, 1.170, 0.755 (Syracuse coding) and period height +0.755", [round(h, 3) for h in h1] == [0.585, 1.17, 0.755] and round(3 * LOG23 - 4, 3) == 0.755)
check("(1,1,2)^N is no-descent for N <= 200 in both codings", min(heights_T((1, 1, 2), 200)) >= 0 and min(heights_syr((1, 1, 2), 200)) >= 0)
# does any excursion of (1,1,2) return to or below its own start? Syracuse coding: bases 0 -> ends .755 (no), .585 -> .755 (no)
print("  (1,1,2): in Syracuse coding no rise-then-fall returns to or below its base (0 -> .755, .585 -> .755); in T coding the base 1.170 -> 1.755 -> 0.755 does")
# the stricter witness
w2 = (1, 1, 1, 1, 2, 3)
hs2 = heights_syr(w2, 1)
hT2 = heights_T(w2, 1)
per2 = 6 * LOG23 - 9
print("  (1,1,1,1,2,3): Syracuse-coded heights %s; T-coded %s; period height %.3f" % ([round(h, 3) for h in hs2], [round(h, 3) for h in hT2], per2))
check("(1,1,1,1,2,3)^N is no-descent (N <= 200) with period height +0.510", min(heights_T(w2, 200)) >= 0 and round(per2, 3) == 0.51)
check("(1,1,1,1,2,3): the excursion based at height 0.585 rises to 2.925 and ends the period at 0.510 < 0.585 while staying >= 0",
      hs2[0] > hs2[-1] >= 0 and max(hs2) > hs2[0])
p2, S2 = affine_S(tuple(int(c) for c in "111110100"))
x2 = Fraction(S2, POW2[9] - POW3[p2])
print("  fixed point of (1,1,1,1,2,3): x = %s (negative rational point of E_inf)" % x2)
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (B1) Section 2: the lowest-set-bit map ===")


def lsb(r):
    return r & -r


def rho(r):
    return 3 * r + lsb(r)


def syr(m):
    y = 3 * m + 1
    v = (y & -y).bit_length() - 1
    return y >> v, v


ident_ok = True
for n in range(1, 2001, 2):
    r = n
    m = n
    d = 0
    S = 0
    for l in range(201):
        ident_ok &= (r == (m << d)) and (r == POW3[l] * n + S)
        # lsb(rho(r)) = lsb(r) 2^v
        m2, v = syr(m)
        r2 = rho(r)
        ident_ok &= (lsb(r2) == lsb(r) << v) and (r2 == (m2 << (d + v)))
        S = 3 * S + (1 << d)
        r, m, d = r2, m2, d + v
check("r_l = 2^(d_l) m_l = 3^l n + S_(l-1) and lsb(rho(r)) = lsb(r) 2^v for odd n <= 2000, l <= 200", ident_ok)
# algebra: 3(2^d m) + 2^d = 2^d (3m+1) = 2^(d+v) U(m): exact identity, also checked above via r2 == m2 << (d+v)
maxsteps = 0
argmax = None
all_reach = True
for n in range(1, 100001):
    r = n
    steps = 0
    while r & (r - 1):
        r = rho(r)
        steps += 1
        if steps > 5000:
            all_reach = False
            break
    if steps > maxsteps:
        maxsteps, argmax = steps, n
check("every n <= 10^5 reaches a power of two under rho", all_reach)
check("maximal number of rho-steps for n <= 10^5 is 129", maxsteps == 129, "max %d at n = %s" % (maxsteps, argmax))


def odd_steps_and_total(n):
    odd = 0
    tot = 0
    while n != 1:
        if n & 1:
            n = 3 * n + 1
            odd += 1
        else:
            n >>= 1
        tot += 1
    return odd, tot


mo = max((odd_steps_and_total(n)[0], n) for n in range(1, 100001))
print("  maximal Syracuse (odd) step count for n <= 10^5: %d at n = %d; its (3n+1, n/2) total stopping time = %d" % (mo[0], mo[1], odd_steps_and_total(mo[1])[1]))
check("rho-step maximum equals the odd-step maximum (129) and the record holder is 77031 with total stopping time 350", mo == (129, 77031) and odd_steps_and_total(77031)[1] == 350)
r = 27
steps = 0
while r & (r - 1):
    r = rho(r)
    steps += 1
check("27 reaches 2^70 after 41 rho-steps", steps == 41 and r == 1 << 70, "reaches 2^%d after %d" % (r.bit_length() - 1, steps))
check("27 has 41 odd steps and 70 halvings", odd_steps_and_total(27) == (41, 111))
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (B2) Section 5: lower approximations of log_2 3, continued fraction, gaps, Catalan ===")
records = []
best = Fraction(0)
A = 0
for p in range(1, 701):
    while POW2[A + 1] < POW3[p]:
        A += 1
    # now 2^A < 3^p < 2^(A+1)
    assert POW2[A] < POW3[p] < POW2[A + 1]
    if Fraction(A, p) > best:
        best = Fraction(A, p)
        records.append((A, p, POW3[p] - POW2[A]))
print("  best lower approximations (A/p records, p <= 700): %s" % [(a, p) for a, p, g in records])


def sci(g):
    """Scientific notation for an integer of any size (floats overflow beyond 1.8e308)."""
    return format(Decimal(g), ".4e")


for a, p, g in records:
    print("    %4d/%3d = %.9f   gap 3^p - 2^A = %s" % (a, p, a / p, g if g < 10 ** 8 else sci(g)))
NOTE_LIST = [(1, 1), (3, 2), (11, 7), (19, 12), (84, 53), (569, 359)]
check("the note's list 1/1, 3/2, 11/7, 19/12, 84/53, 569/359 is exactly the records for p < 400", [(a, p) for a, p, g in records if p < 400] == NOTE_LIST)
gaps = {(a, p): g for a, p, g in records}
check("gaps 1, 1, 139, 7153 exact", [gaps[(1, 1)], gaps[(3, 2)], gaps[(11, 7)], gaps[(19, 12)]] == [1, 1, 139, 7153])
check("gap at 84/53 = 4.043e22 and at 569/359 = 2.061e168", "%.3e" % gaps[(84, 53)] == "4.043e+22" and "%.3e" % gaps[(569, 359)] == "2.061e+168",
      "%.3e, %.3e" % (gaps[(84, 53)], gaps[(569, 359)]))
check("84/53 and 569/359 lie below log_2 3 (2^84 < 3^53, 2^569 < 3^359)", POW2[84] < POW3[53] and POW2[569] < POW3[359])
check("8/5 is an UPPER approximation (2^8 = 256 > 243 = 3^5), correctly absent", POW2[8] > POW3[5])
# continued fraction of log_2 3 with 200-digit Decimal
getcontext().prec = 220
alpha = Decimal(3).ln() / Decimal(2).ln()
cf = []
x = alpha
for _ in range(16):
    a = int(x)
    cf.append(a)
    x = 1 / (x - a)
print("  continued fraction of log_2 3: %s" % cf)
check("CF of log_2 3 begins [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55]", cf[:15] == [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55])
conv = []
p0, q0, p1, q1 = 1, 0, cf[0], 1
conv.append((p1, q1))
for a in cf[1:12]:
    p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
    conv.append((p1, q1))
print("  convergents: %s" % conv)
lower_conv = [(p, q) for p, q in conv if Fraction(p, q) < Fraction(alpha.__floor__() if False else 0) or POW2[p] < POW3[q]] if False else [(p, q) for p, q in conv if p < 5000 and POW2[p] < POW3[q]]
print("  lower convergents (2^p < 3^q): %s; the listed 11/7 and 569/359 are intermediate fractions, not convergents" % lower_conv)
check("1/1, 3/2, 19/12, 84/53 are convergents; 11/7 and 569/359 are semiconvergents",
      all(t in conv for t in [(1, 1), (3, 2), (19, 12), (84, 53)]) and (11, 7) not in conv and (569, 359) not in conv)
# Catalan: 3^p - 2^A = +-1
sols_plus = []
sols_minus = []
A = 0
for p in range(1, 3001):
    while POW2[A + 1] < POW3[p]:
        A += 1
    if POW3[p] - POW2[A] == 1:
        sols_plus.append((p, A))
    if POW2[A + 1] - POW3[p] == 1:
        sols_minus.append((p, A + 1))
print("  3^p - 2^A = +1 for p <= 3000: %s;  2^A - 3^p = +1: %s" % (sols_plus, sols_minus))
check("3^p - 2^A = 1 only at (p,A) = (1,1), (2,3) [3-2, 9-8]: the note's 'gap 1 twice' is complete for the lower side", sols_plus == [(1, 1), (2, 3)])
check("2^A - 3^p = 1 only at (p,A) = (1,2) [4-3], an UPPER approximation irrelevant to negative cycles", sols_minus == [(1, 2)])
# 139 | S for -17
u17 = (1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0)
p17, S17 = affine_S(u17)
check("-17: T-word 11110111000, p = 7, S = 2363 = 17 * 139, x = S/(2^11 - 3^7) = -17", (p17, S17) == (7, 2363) and Fraction(S17, POW2[11] - POW3[7]) == -17, "p=%d S=%d" % (p17, S17))
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== (B3) Section 3: two-place identities, cycle series ===")


def two_place(n, L):
    """Check 2^(d_l) m_l/3^l = n C_l and S_(l-1)/3^l = n (C_l - 1) exactly for l <= L; return the real partial sums."""
    m = n
    C = Fraction(1)
    S = 0
    d = 0
    sums = []
    ok = True
    for l in range(1, L + 1):
        C *= 1 + Fraction(1, 3 * m)
        S = 3 * S + (1 << d)
        m, v = syr(m)
        d += v
        ok &= (Fraction(m << d, POW3[l]) == n * C) and (Fraction(S, POW3[l]) == n * (C - 1))
        # 2-adic: v_2(S + 3^l n) = d_l
        t = S + POW3[l] * n
        ok &= (t != 0 and (t & -t).bit_length() - 1 == d)
        sums.append(Fraction(S, POW3[l]))
    return ok, sums


for n in (27, -1, -5, -17):
    ok, sums = two_place(n, 60)
    check("identities 2^(d_l) m_l/3^l = n C_l, S_(l-1)/3^l = n(C_l - 1), v_2(S_(l-1) + 3^l n) = d_l exact for n=%d, l <= 60" % n, ok)
    print("    n=%d real partial sums l=1..5: %s ... l=22,23,24: %s ... l=60: %.6f" % (n, ["%.4f" % float(s) for s in sums[:5]], ["%.4f" % float(s) for s in sums[21:24]], float(sums[59])))
ok27, s27 = two_place(27, 24)
check("27: partial sums 0.3333, 0.5556, 0.8519, 1.0494, 1.1811 and 2.1370 at l = 24 (note: 0.33, 0.56, 0.85, ..., 2.14)",
      ["%.4f" % float(s) for s in s27[:5]] == ["0.3333", "0.5556", "0.8519", "1.0494", "1.1811"] and "%.4f" % float(s27[23]) == "2.1370")
random.seed(20260926)
rand_ok = True
for _ in range(200):
    n = random.choice([1, -1]) * (2 * random.randrange(1, 10 ** 6) + 1)
    rand_ok &= two_place(n, 40)[0]
check("the identities hold exactly for 200 random odd n of both signs, l <= 40", rand_ok)


def cycle_series_exact(valword):
    """Exact real value of sum_t 2^(A_t)/3^(t+1) for the periodic valuation word: (one period)/(1 - 2^A/3^P)."""
    P = len(valword)
    A = sum(valword)
    Acum = 0
    s = Fraction(0)
    for t, v in enumerate(valword):
        s += Fraction(POW2[Acum], POW3[t + 1])
        Acum += v
    return s / (1 - Fraction(POW2[A], POW3[P]))


for n, w in ((-1, (1,)), (-5, (1, 2)), (-17, (1, 1, 1, 2, 1, 1, 4))):
    val = cycle_series_exact(w)
    check("real Bernstein series of the cycle point %d (word %s) equals -n = %d exactly" % (n, w, -n), val == -n, str(val))
print("  (the 2-adic value is -n by v_2(S_(l-1) + 3^l n) = d_l -> infinity, checked above; so real and 2-adic values coincide on cycles)")
# eta_l >= 1/3: C_inf/C_l - 1 >= 1/(3 m_l) since the first omitted factor alone gives it (algebraic; trivial check on truncations)
print("  eta_l = m_l (C_inf/C_l - 1) >= m_l * (1/(3 m_l)) = 1/3: the first factor of the tail product; algebraic, no computation needed")
print("  [%s]" % elapsed())

# ================================================================================================================
print("=== hashes of the audited files (raw bytes) ===")
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", ".."))
for rel in ("04-computation/experiments/collatz_directions_20260926.py", "04-computation/experiments/collatz_directions_20260926.out", "05-knowledge/results/collatz_directions_20260926.md"):
    path = os.path.join(root, rel)
    try:
        with open(path, "rb") as f:
            print("  %s  %s" % (hashlib.sha256(f.read()).hexdigest(), rel))
    except OSError as e:
        print("  (cannot read %s: %s)" % (rel, e))

print("=== SUMMARY: %d PASS, %d FAIL  [%s] ===" % (PASSES, FAILS, elapsed()))
