#!/usr/bin/env python3
"""collatz_shadow_error_20260927_audit.py -- independent adversarial audit of
05-knowledge/results/collatz_shadow_error_flp_20260927.md, HYP-9164 and
04-computation/experiments/collatz_shadow_error_20260927.py (auditor agent, 2026-09-27).

Everything is exact (fractions / integers); no float is used except for display.
The code is written from the definitions, not copied from the script under audit.

Sections
 (A1) Proposition 1: recursion, eta_l = m_l (C_K/C_l - 1) (truncated form, exact), x_l = xi 3^l/2^(d_l),
      and the general solution of the affine recursion (the 'bounded-ratio uniqueness' wording).
 (A2) Proposition 2: eta >= 1 on exact tails; largest partial eta_0 over no-descent words of length 8..16.
 (A3) Proposition 3: budget (note's form and the sharper 2^(d_(l+k)-d_l) <= 3^k eta_l), dip bounds.
 (A4) Proposition 4: periodic values, R = R_2 for periodic and eventually periodic words, cycle cancellation.
 (A5) Proposition 5: minus sheet recursion / cancellation / orbits reaching 1; the minus-sheet
      counterexamples to the 'equivalently no positive integer orbit ... has sup eta < inf' phrasing.
 (A6) Proposition 6: carries c_l = 3{eta_l} - 2^v {eta_(l+1)}, exact, range; N-coordinate formula; Mahler ceil map.
 (A7) Proposition 7 / bounded shadow: valuation bound, f maps [1,2) to itself, independent enumeration of the
      periodic f-itineraries (letters {1,2} forced), maxima of (1^a 2), trivial bound 5/3 and its exceptions.
 (B)  All quoted numbers of the note re-derived (orbits 27, 703, 871, 6171; -5 families; xi for 27; carries).
"""
import math
import itertools
from fractions import Fraction as Fr

# ----------------------------------------------------------------------------- basic maps
def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v

def Um(m):
    m = 3 * m - 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v

def orbit(n, step, K):
    ms = [n]; w = []
    for _ in range(K):
        m, v = step(ms[-1]); ms.append(m); w.append(v)
    return ms, w

def steps_to_one(n, step):
    m = n; K = 0
    while m != 1:
        m, _ = step(m); K += 1
    return K

def S_of(w):
    """S_w = sum_(t<p) 3^(p-1-t) 2^(d_t), d_0 = 0"""
    p = len(w); S = 0; d = 0
    for t in range(p):
        S += 3 ** (p - 1 - t) * 2 ** d; d += w[t]
    return S

def R_periodic(w):
    """real value of sum_(k>=0) 2^(d_k)/3^(k+1) for the word w^inf; None if 2^A >= 3^p (divergent)"""
    p = len(w); A = sum(w)
    if 2 ** A >= 3 ** p:
        return None
    return Fr(S_of(w), 3 ** p - 2 ** A)

def R_head(w):
    s = Fr(0); d = 0
    for t, v in enumerate(w):
        s += Fr(2 ** d, 3 ** (t + 1)); d += v
    return s

def R_evper(head, w):
    """real value of the eventually periodic word head w^inf (None if divergent)"""
    Rw = R_periodic(w)
    if Rw is None:
        return None
    return R_head(head) + Fr(2 ** sum(head), 3 ** len(head)) * Rw

def R2_mod(head, w, N):
    """2-adic value of head w^inf modulo 2^N (sum of 2^(d_t) 3^(-t-1) until 2^(d_t) = 0 mod 2^N)"""
    mod = 2 ** N; inv3 = pow(3, -1, mod)
    s = 0; d = 0; t = 0
    while d < N:
        v = head[t] if t < len(head) else w[(t - len(head)) % len(w)]
        s = (s + 2 ** d * pow(inv3, t + 1, mod)) % mod
        d += v; t += 1
    return s

def x_plus(w):
    """2-adic point whose plus-sheet halving word is w^inf: 2^A x = 3^p x + S"""
    return Fr(S_of(w), 2 ** sum(w) - 3 ** len(w))

def x_minus(w):
    """2-adic point whose minus-sheet halving word is w^inf: 2^A x = 3^p x - S"""
    return Fr(S_of(w), 3 ** len(w) - 2 ** sum(w))

def eta_trunc(w):
    """truncated errors eta_l^(K) = sum_(k<K-l) 2^(d_(l+k)-d_l)/3^(k+1), l = 0..K-1, and the d_l"""
    K = len(w); d = [0]
    for v in w:
        d.append(d[-1] + v)
    etas = [sum(Fr(2 ** (d[l + k] - d[l]), 3 ** (k + 1)) for k in range(K - l)) for l in range(K)]
    return etas, d

def eventually_periodic_word(n, step, maxlen=400):
    """(head, period) of the halving word of the orbit of n under step (assumes the orbit enters a cycle)"""
    seen = {}; ms = [n]; w = []
    while ms[-1] not in seen and len(w) < maxlen:
        seen[ms[-1]] = len(w)
        m, v = step(ms[-1]); ms.append(m); w.append(v)
    if ms[-1] not in seen:
        return None
    s = seen[ms[-1]]
    return w[:s], w[s:]

def eta_exact_along(n, step, L):
    """exact eta_l, l < L, along an eventually periodic orbit (None entries if the tail series diverges)"""
    hp = eventually_periodic_word(n, step)
    if hp is None:
        return None
    head, per = hp
    full = head + per * (L // len(per) + 2)
    out = []
    for l in range(L):
        if l < len(head):
            out.append(R_evper(full[l:len(head)], per))
        else:
            rot = (l - len(head)) % len(per)
            out.append(R_periodic(per[rot:] + per[:rot]))
    return out

def sha_line(title):
    print("\n== %s ==" % title)

ok_all = True
def report(flag, text):
    global ok_all
    ok_all &= bool(flag)
    print(" [%s] %s" % ("ok" if flag else "FAIL", text))

# ============================================================================= (A1) Proposition 1
sha_line("(A1) Proposition 1: recursion, closed form, shadow, and the uniqueness wording")
for n in (27, 703, 871, 6171):
    K = steps_to_one(n, U)
    ms, w = orbit(n, U, K)
    etas, d = eta_trunc(w)
    rec = all(2 ** w[l] * etas[l + 1] == 3 * etas[l] - 1 for l in range(K - 1))
    # truncated closed form: eta_l^(K) = m_l (C_K/C_l - 1), C_l = prod_(j<l)(1 + 1/(3 m_j)); exact
    C = [Fr(1)]
    for j in range(K):
        C.append(C[-1] * (1 + Fr(1, 3 * ms[j])))
    closed = all(etas[l] == ms[l] * (C[K] / C[l] - 1) for l in range(K))
    xi = n * C[K]
    shadow = all(xi * 3 ** l / 2 ** d[l] == ms[l] + etas[l] for l in range(K))
    xi2 = n + etas[0]
    report(rec and closed and shadow and xi == xi2,
           "n=%d (K=%d): recursion exact; eta_l = m_l(C_K/C_l - 1) exact; xi 3^l/2^(d_l) = m_l + eta_l exact; xi = n C_K = n + eta_0 = %.6f" % (n, K, float(xi)))
# general solution of 2^v y' = 3y - 1: y_l = eta_l + lambda 3^l/2^(d_l); its ratio to m_l stays bounded (-> lambda/(n C_inf))
n = 27; K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w)
lam = Fr(7, 3)
ys = [etas[l] + lam * Fr(3 ** l, 2 ** d[l]) for l in range(K)]
rec2 = all(2 ** w[l] * ys[l + 1] == 3 * ys[l] - 1 for l in range(K - 1))
ratios = [float(ys[l] / ms[l]) for l in range(K)]
report(rec2, "general solution y_l = eta_l + lambda 3^l/2^(d_l) (lambda = 7/3) also satisfies 2^v y' = 3y - 1; y_l/m_l ranges in [%.3f, %.3f] (bounded, not -> 0): uniqueness needs y_l = o(m_l), 'bounded ratio' is not enough" % (min(ratios), max(ratios)))

# ============================================================================= (A2) Proposition 2
sha_line("(A2) Proposition 2: eta >= 1 with equality iff all later valuations are 1; no-descent maxima")
neg_cycles = {(-1,): [(1,)], (-5,): [(1, 2), (2, 1)], (-17,): [(1, 1, 1, 2, 1, 1, 4)]}
allw = [(1,), (1, 2), (2, 1), (1, 1, 1, 2, 1, 1, 4), (1, 1, 2, 1, 1, 4, 1), (1, 1, 2), (1, 2, 1, 2, 2, 1), (1, 1, 1, 1, 2), (2, 1, 1, 1, 1)]
report(all(R_periodic(w) >= 1 for w in allw) and R_periodic((1,)) == 1 and all(R_periodic(w) > 1 for w in allw if w != (1,)),
       "eta >= 1 on the exact periodic tails listed, = 1 exactly for (1)^inf, > 1 otherwise")
# lower bound sum 2^k/3^(k+1) = 1 : exact partial sums
report(all(sum(Fr(2 ** k, 3 ** (k + 1)) for k in range(N)) == 1 - Fr(2 ** N, 3 ** N) for N in (5, 12, 30)),
       "all-ones partial sums are 1 - (2/3)^N exactly (so truncated tails may fall below 1: eta >= 1 is a statement about the exact tail)")

def nodescent_words(L):
    """all words (v_1..v_L) with 2^(d_j) < 3^j for j = 1..L (exact integer test)"""
    out = []
    def rec(prefix, d):
        j = len(prefix)
        if j == L:
            out.append(tuple(prefix)); return
        v = 1
        while 2 ** (d + v) < 3 ** (j + 1):
            rec(prefix + [v], d + v); v += 1
    rec([], 0)
    return out

print(" largest partial eta_0 = sum_(t<L) 2^(d_t)/3^(t+1) over no-descent words of length L (2^(d_j) < 3^j for all j <= L):")
for L in (8, 10, 12, 14, 16):
    words = nodescent_words(L)
    best = max(words, key=lambda w: R_head(w))
    val = R_head(best)
    print("   L=%2d: %6d words, max = %.4f = %.3f L (word %s)" % (L, len(words), float(val), float(val) / L, best))
    if L == 12:
        report(abs(float(val) - 2.9518) < 5e-4, "L = 12 maximum %.4f matches the note's 2.95" % float(val))

# ============================================================================= (A3) Proposition 3
sha_line("(A3) Proposition 3: budget and dip bounds")
# exact periodic tails: at every point of the three negative cycles and of the 8 itineraries, 2^v <= 3 eta - 1 and
# the sharper 2^(d_(l+k)-d_l) <= 3^k eta_l (k >= 0), the note's 2^(d_(l+k)-d_l) <= 3^(k+1) eta_l
def check_budget_periodic(w):
    p = len(w); res = True
    for r in range(p):
        wr = w[r:] + w[:r]; eta = R_periodic(wr)
        d = 0
        for k in range(3 * p + 1):
            if 2 ** d > 3 ** k * eta:
                res = False
            d += wr[k % p]
        if 2 ** wr[0] > 3 * eta - 1:
            res = False
    return res
report(all(check_budget_periodic(w) for w in allw), "on every rotation of the exact periodic tails: 2^v <= 3 eta - 1 and 2^(d_(l+k)-d_l) <= 3^k eta_l for all k (sharper than the note's 3^(k+1) eta_l)")
# derivation used: 2^(d_(l+k)-d_l) eta_(l+k) = 3^k eta_l - S_k with S_k = sum_(j<k) 3^(k-1-j) 2^(d_(l+j)-d_l) >= 0 -- check exactly on truncated data
n = 27; K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w)
ident = True
for l in range(K):
    for k in range(K - l):
        Sk = sum(3 ** (k - 1 - j) * 2 ** (d[l + j] - d[l]) for j in range(k))
        ident &= (2 ** (d[l + k] - d[l]) * etas[l + k] == 3 ** k * etas[l] - Sk)
        ident &= (2 ** (d[l + k] - d[l]) * ms[l + k] == 3 ** k * ms[l] + Sk)
report(ident, "identities 2^(d_(l+k)-d_l) eta_(l+k) = 3^k eta_l - S_k and 2^(d_(l+k)-d_l) m_(l+k) = 3^k m_l + S_k hold exactly on the orbit of 27 (so x_(l+k) = 3^k x_l/2^(d_(l+k)-d_l))")
for n in (27, 703, 871, 6171):
    K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w)
    note_dip = all(ms[l + k] >= Fr(ms[l], 3 * etas[l]) for l in range(K) for k in range(K - l))
    sharp_dip_fail = sum(1 for l in range(K) for k in range(K - l) if ms[l + k] < ms[l] / etas[l])
    exact_budget_trunc = all(2 ** w[l] <= 3 * etas[l] - 1 for l in range(K - 1))
    exact_budget_30 = all(2 ** w[l] <= 3 * etas[l] - 1 for l in range(K - 1) if K - l >= 30)
    always_budget = all(2 ** w[l] <= 9 * etas[l] - 3 for l in range(K - 1))
    print("   n=%d: note's dip bound m_(l+k) >= m_l/(3 eta_l) on truncated data: %s; the sharper m_(l+k) >= m_l/eta_l (valid for exact tails, needs eta_(l+k) >= 1) fails at %d truncated pairs (expected near the truncation); exact-tail budget 2^v <= 3 eta_l - 1 on all truncated values: %s (expected to fail near the end), on tails of >= 30 terms: %s; the truncation-proof 2^v <= 9 eta_l - 3: %s" % (n, note_dip, sharp_dip_fail, exact_budget_trunc, exact_budget_30, always_budget))

# ============================================================================= (A4) Proposition 4
sha_line("(A4) Proposition 4: periodic tails, eta = |x_w|, real = 2-adic, cancellation on the negative cycles")
expected = {(1,): Fr(1), (1, 2): Fr(5), (2, 1): Fr(7), (1, 1, 1, 2, 1, 1, 4): Fr(17), (1, 1, 2, 1, 1, 4, 1): Fr(25), (1, 1, 2): Fr(19, 11), (1, 2, 1, 2, 2, 1): Fr(1213, 217)}
for w, val in expected.items():
    R = R_periodic(w); x = x_plus(w)
    r2 = R2_mod([], list(w), 200)
    two_adic_agree = (R.numerator * pow(R.denominator, -1, 2 ** 200) - r2) % 2 ** 200 == 0
    report(R == val and R + x == 0 and two_adic_agree, "w=%s: R(w^inf) = %s (note: %s), x_w = %s, R + x_w = 0, R == R_2 mod 2^200: %s" % (w, R, val, x, two_adic_agree))
for n in (-1, -5, -17):
    L = 24
    ms, w = orbit(n, U, L)
    etas = eta_exact_along(n, U, L)
    report(all(ms[l] + etas[l] == 0 for l in range(L)), "negative cycle through %d: m_l + eta_l = 0 at all %d points checked (values %s)" % (n, L, sorted(set(ms[:L]))))

# ============================================================================= (A5) Proposition 5
sha_line("(A5) Proposition 5: minus sheet")
for n in (1, 5, 17):
    L = 30
    ms, w = orbit(n, Um, L)
    etas = eta_exact_along(n, Um, L)
    rec = all(2 ** w[l] * etas[l + 1] == 3 * etas[l] - 1 for l in range(L - 1))
    xs = [ms[l] - etas[l] for l in range(L)]
    shadow = all(2 ** w[l] * xs[l + 1] == 3 * xs[l] for l in range(L - 1))
    report(rec and shadow and all(x == 0 for x in xs), "3x-1 cycle through %d (points %s): recursion 2^v eta' = 3 eta - 1 exact, x = m - eta obeys 2^v x' = 3x, and m - eta = 0 everywhere" % (n, sorted(set(ms[:L]))))
for n in (3, 11, 29):
    hp = eventually_periodic_word(n, Um); head, per = hp
    R = R_evper(head, per)
    report(R == n and per == [1], "3x-1 orbit of %d: word head %s then (1)^inf; R(d) = %s = n" % (n, head, R))
# by hand for n = 3: 1/3 + 8/9 + 16/27 + ... = 1/3 + (8/9)/(1 - 2/3) = 3
report(Fr(1, 3) + Fr(8, 9) / (1 - Fr(2, 3)) == 3, "hand check n=3: 1/3 + (8/9)/(1-2/3) = 3")
# general claim: eventually periodic words (3^p > 2^A) have R = R_2 as rationals (also for non-integer values)
gen_ok = True; cnt = 0
for head in ([], [3], [2, 5], [1, 7], [4, 1, 2], [6, 6]):
    for per in ([1], [1, 2], [1, 1, 2], [1, 2, 1, 2, 2, 1], [1, 1, 1, 2, 1, 1, 4], [2, 1, 1, 1]):
        R = R_evper(head, per)
        if R is None:
            continue
        r2 = R2_mod(head, per, 300)
        gen_ok &= (R.numerator * pow(R.denominator, -1, 2 ** 300) - r2) % 2 ** 300 == 0
        cnt += 1
report(gen_ok, "R(d) = R_2(d) (mod 2^300, as rationals with odd denominator) for %d eventually periodic words with convergent real part, integer-valued or not" % cnt)
# minus-sheet foundry Prop 6(4): R(d) <= n for every positive minus-sheet orbit, = n iff the orbit is eventually periodic: check R = n on all odd n <= 2001 (all reach a cycle)
allR = True
for n in range(1, 2002, 2):
    hp = eventually_periodic_word(n, Um, 5000)
    if hp is None:
        allR = False; continue
    allR &= (R_evper(*hp) == n)
report(allR, "R(d) = n for every odd n <= 2001 on the minus sheet (all these orbits enter a 3x-1 cycle): equality is the eventually periodic case (foundry Prop. 6(4), not cited by the note)")
# the HYP's 'equivalently no positive integer orbit on either sheet has sup eta < infinity': minus-sheet counterexamples
print(" minus-sheet positive integers with BOUNDED shadow error (eventually periodic orbits; eta_l exact):")
for n in (1, 3, 5, 7, 11, 17, 25, 29, 37, 41, 55, 61, 91, 3073, 65537):
    L = 60
    etas = eta_exact_along(n, Um, L)
    if etas is None:
        print("   n=%5d: orbit did not enter a cycle within the step cap (skipped)" % n); continue
    print("   n=%5d: sup_(l<%d) eta_l = %-10s  (all eta_l < 2: %s)" % (n, L, max(etas), all(e < 2 for e in etas)))
report(all(e == 1 for e in eta_exact_along(1, Um, 40)),
       "n = 1 on the minus sheet: eta_l = 1 for all l -- a positive integer orbit with all eta_l < 2 (the word (1)^inf IS an f-itinerary and its minus-sheet 2-adic point is the positive integer 1)")
# plus sheet: every orbit reaching 1 has a divergent tail series (tail (2)^inf, 4 > 3)
report(R_periodic((2,)) is None, "plus sheet: the 1-cycle word (2)^inf has 2^A = 4 > 3 = 3^p, so eta = infinity for every orbit reaching 1 (the HYP's parenthesis is right on the plus sheet only)")

# ============================================================================= (A6) Proposition 6
sha_line("(A6) FLP carry structure, N-coordinates, Mahler's map")
def frac(q):
    return q - (q.numerator // q.denominator)
for n in (27, 703, 871, 6171):
    K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w)
    Ms = [ms[l] + etas[l].numerator // etas[l].denominator for l in range(K)]
    cs = [2 ** w[l] * Ms[l + 1] - 3 * Ms[l] for l in range(K - 1)]
    formula = all(cs[l] == 3 * frac(etas[l]) - 2 ** w[l] * frac(etas[l + 1]) for l in range(K - 1))
    rng = all(-(2 ** w[l]) + 1 <= cs[l] <= 2 for l in range(K - 1))
    vmax = max(w)
    report(formula and rng, "n=%d: c_l = 3{eta_l} - 2^v{eta_(l+1)} exactly, integer, in {-2^v+1..2}; values seen %s; max valuation %d (so the widest range is {%d..2})" % (n, sorted(set(cs)), vmax, -(2 ** vmax) + 1))
# N-coordinates: N = (m+1)/2, N' = (3N - 1 + 2^(v-1))/2^v
ncoord = True
for m in range(1, 200001, 2):
    m2, v = U(m); N = (m + 1) // 2; N2 = (m2 + 1) // 2
    num = 3 * N - 1 + 2 ** (v - 1)
    ncoord &= (num % 2 ** v == 0 and num // 2 ** v == N2)
    if v == 1:
        ncoord &= (N % 2 == 0 and N2 * 2 == 3 * N)
    if v == 2:
        ncoord &= (N2 * 4 == 3 * N + 1)
report(ncoord, "N = (m+1)/2 -> (3N - 1 + 2^(v-1))/2^v for all odd m <= 200001; v = 1 iff N even and then N' = 3N/2; v = 2 gives (3N+1)/4")
# Mahler: xi (3/2)^n = g_n + f_n with f_n < 1/2 forces g_(n+1) = ceil(3 g_n / 2); check the algebra on a grid of (g, f)
mahler = True
for g in range(0, 40):
    for k in range(0, 60):
        f = Fr(k, 120)  # f in [0, 1/2)
        y = Fr(3, 2) * (g + f); g2 = y.numerator // y.denominator; f2 = y - g2
        if f2 < Fr(1, 2):
            mahler &= (g2 == -((-3 * g) // 2))  # ceil(3g/2)
            mahler &= ((g % 2 == 0 and f < Fr(1, 3)) or (g % 2 == 1 and f >= Fr(1, 3)))
report(mahler, "Mahler: if {xi(3/2)^n} < 1/2 at two consecutive times then g_(n+1) = ceil(3 g_n/2), with f < 1/3 at even g and f in [1/3, 1/2) at odd g")
# FLP carry count in the same normalisation: q g' = p g + c, c = p f - q f' in (-q, p): p + q - 1 integers
p, q = 3, 2
report(len(range(-q + 1, p)) == p + q - 1 == 4, "FLP carry c = p f_n - q f_(n+1) lies in (-q, p): %d integer values for p/q = 3/2 (the note says 'ranges over p values')" % len(range(-q + 1, p)))
# FLP identity xi (3/2)^l = 2^(e_l) (m_l + eta_l), e_l = d_l - l >= 0, exact on 27's truncated shadow
n = 27; K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w); xi = n + etas[0]
report(all(xi * Fr(3, 2) ** l == 2 ** (d[l] - l) * (ms[l] + etas[l]) and d[l] - l >= 0 for l in range(K)),
       "xi (3/2)^l = 2^(e_l)(m_l + eta_l) with e_l = d_l - l >= 0 exactly (so {xi (3/2)^l} = {2^(e_l) eta_l})")

# ============================================================================= (A7) bounded shadow
sha_line("(A7) bounded-shadow regime, the interval map f, the eight itineraries, the 5/3 bound")
for B, vb in ((2, 2), (5, 3), (17, 5)):
    report(max(v for v in range(1, 40) if 2 ** v <= 3 * B - 1) == vb and 2 ** (vb + 1) > 3 * B - 1, "eta < B = %d forces 2^v <= 3 eta - 1 < %d, i.e. v <= %d" % (B, 3 * B - 1, vb))
# branch determination: eta in [1,2), eta' = (3 eta - 1)/2^v in [1,2) forces v = 1 iff eta < 5/3, v = 2 iff eta >= 5/3, never v >= 3
def f(e):
    return (3 * e - 1) / 2 if e < Fr(5, 3) else (3 * e - 1) / 4
branch = True
for k in range(0, 3000):
    e = 1 + Fr(k, 3000)
    admissible = [v for v in range(1, 6) if 1 <= (3 * e - 1) / 2 ** v < 2]
    branch &= (admissible == [1 if e < Fr(5, 3) else 2]) and 1 <= f(e) < 2
report(branch, "for eta in [1,2) exactly one v keeps eta' in [1,2): v = 1 iff eta < 5/3, v = 2 iff eta >= 5/3; f maps [1,2) into itself (second branch lands in [1,5/4))")

def necklaces(p, letters=(1, 2)):
    """primitive necklace representatives (lexicographically least rotation) of length p"""
    out = []
    for w in __import__('itertools').product(letters, repeat=p):
        rots = [w[i:] + w[:i] for i in range(p)]
        if w == min(rots) and all(r != w for r in rots[1:]):
            out.append(w)
    return out

def itinerary_words(pmax, letters=(1, 2)):
    found = []
    for p in range(1, pmax + 1):
        for w in necklaces(p, letters):
            R = R_periodic(w)
            if R is None:
                continue
            vals = [R_periodic(w[r:] + w[:r]) for r in range(p)]   # eta at each position, each from its own geometric series
            if not all(1 <= e < 2 for e in vals):
                continue
            # recursion consistency and itinerary of f
            e = R; it = []
            for r in range(p):
                it.append(1 if e < Fr(5, 3) else 2)
                e = f(e)
            found.append((w, R, max(vals), tuple(it) == w and e == R, vals))
    return found
found10 = itinerary_words(10)
report(len(found10) == 8 and all(fl for _, _, _, fl, _ in found10), "independent enumeration: %d primitive periodic words with p <= 10 over {1,2} with all eta in [1,2); itinerary of f = word in every case" % len(found10))
quoted = {(1,): (Fr(1), 1.0), (1, 1, 1, 1, 2): (Fr(211, 179), 1.905), (1, 1, 1, 1, 1, 2): (Fr(665, 601), 1.809), (1, 1, 1, 1, 1, 1, 2): (Fr(2059, 1931), 1.755),
          (1, 1, 1, 1, 1, 1, 1, 2): (Fr(6305, 6049), 1.723), (1, 1, 1, 1, 1, 1, 1, 1, 2): (Fr(19171, 18659), 1.703), (1, 1, 1, 1, 1, 1, 1, 1, 1, 2): (Fr(58025, 57001), 1.691),
          (1, 1, 1, 1, 1, 2, 1, 1, 1, 2): (Fr(62185, 54953), 1.9994)}
for w, R, mx, fl, vals in found10:
    qR, qmax = quoted.get(w, (None, None))
    xw = x_plus(w); xm = x_minus(w)
    report(qR == R and abs(float(mx) - qmax) < 6e-4, "w=%s eta_0 = %s max = %.4f (quoted %s, %.4f); plus-sheet point %s (integer: %s), minus-sheet point %s (integer: %s)" % (w, R, float(mx), qR, qmax, xw, xw.denominator == 1, xm, xm.denominator == 1))
# letters >= 3 are excluded in the [1,2) regime (2^v <= 3 eta - 1 < 5); confirm that allowing letter 3 adds nothing
found10b = itinerary_words(10, letters=(1, 2, 3))
report(len(found10b) == 8, "allowing the letter 3 in the enumeration adds nothing (%d words): the restriction to {1,2} is justified by 2^v <= 3 eta - 1 < 5" % len(found10b))
found14 = itinerary_words(14)
print("   with p <= 14 there are %d such words: %s" % (len(found14), [w for w, *_ in found14 if len(w) > 10]))
# maxima of (1^a 2) decrease to 5/3
prev = None; mono = True
print("   maxima of the family (1^a 2): ", end="")
for a in range(4, 16):
    w = (1,) * a + (2,)
    mx = max(R_periodic(w[r:] + w[:r]) for r in range(len(w)))
    print("a=%d: %.5f  " % (a, float(mx)), end="")
    mono &= (prev is None or mx < prev) and mx > Fr(5, 3); prev = mx
print()
report(mono, "maxima of (1^a 2) strictly decrease and stay above 5/3 for a = 4..15; the limit is 1/3 + 4/3 = 5/3 (value before the 2 with a long run of ones behind it)")
# the trivial bound: every word other than 1^inf has sup eta >= 5/3 (a position before a letter v >= 2 has eta >= (2^v + 1)/3 >= 5/3)
report(all(min(max(R_periodic(w[r:] + w[:r]) for r in range(len(w))) for w in necklaces(p) if R_periodic(w) is not None) >= Fr(5, 3) for p in range(2, 9)),
       "every periodic word with p in 2..8 (letters {1,2}, convergent) has max eta >= 5/3; the only word with sup eta < 5/3 is (1)^inf (point -1 on the plus sheet, +1 on the minus sheet)")

# ============================================================================= (B) quoted numbers
sha_line("(B) quoted numbers: truncated maxima, eta_0, carries, -5 families, xi for 27")
quoted_max = {27: (332, 32, 3077), 703: (17293, 34, 83501), 871: (9909, 14, 63665), 6171: (61981, 32, 325133)}
quoted_eta0 = {27: 5.3689, 703: 147.5929, 871: 137.1101, 6171: 1179.3637}
for n in (27, 703, 871, 6171):
    K = steps_to_one(n, U); ms, w = orbit(n, U, K); etas, d = eta_trunc(w)
    imax = max(range(K), key=lambda l: etas[l])
    qm, ql, qmm = quoted_max[n]
    report(round(float(etas[imax])) == qm and imax == ql and ms[imax] == qmm and abs(float(etas[0]) - quoted_eta0[n]) < 1e-4,
           "n=%d: K=%d, max truncated eta = %.3f at l=%d, m_l=%d (note: %d at m=%d); eta_0^(K) = %.4f; min of the orbit after l=%d is %d" % (n, K, float(etas[imax]), imax, ms[imax], qm, qmm, float(etas[0]), imax, min(ms[imax:K])))
    if n == 27:
        report(abs(float(27 + etas[0]) - 32.368923) < 1e-6, "xi(27, truncated) = 27 + eta_0 = %.6f (note: 32.368923)" % float(27 + etas[0]))
fam = True
for k in range(1, 9):
    n = 4 * 8 ** k - 5
    ms, w = orbit(n, U, 2 * k)
    fam &= (w == [1, 2] * k) and (R_head(w) == 5 * (1 - Fr(8, 9) ** k))
report(fam, "-5 families n_k = 4 8^k - 5, k = 1..8: word (1,2)^k and partial eta_0 = 5(1 - (8/9)^k) exactly")

# ============================================================================= (A8) Proposition 7 (added in commit 76db2bffb) and the itineraries exploration
sha_line("(A8) Proposition 7: runs in the B = 2 regime; replication of the 490-itinerary exploration")
g = lambda e: (3 * e - 1) / 2
# images of [1, 5/4) under g: endpoints map to endpoints (g increasing)
ends = [Fr(5, 4)]
for _ in range(3):
    ends.append(g(ends[-1]))
report(ends == [Fr(5, 4), Fr(11, 8), Fr(25, 16), Fr(59, 32)] and Fr(25, 16) < Fr(5, 3) < Fr(59, 32) and g(Fr(5, 3)) == 2 and (3 * Fr(2) - 1) / 4 == Fr(5, 4),
       "after v = 2 the error is in [1, 5/4); g-images [1,11/8), [1,25/16) lie below 5/3, the third image [1,59/32) = [1,1.84375) reaches 5/3: at least three v = 1 after every v = 2")
def runs_ok(w, cyclic=True):
    """every 2 is followed by at least three 1's (cyclically for periodic words)"""
    p = len(w); ww = w + w[:3] if cyclic else w
    for i, v in enumerate(w):
        if v == 2:
            nxt = ww[i + 1:i + 4]
            if len(nxt) == 3 and any(x != 1 for x in nxt):
                return False
    return True
report(all(runs_ok(w) for w, *_ in found14) and all(w.count(2) * 4 <= len(w) for w, *_ in found14),
       "all %d periodic f-itineraries with p <= 14 obey the three-ones rule and have frequency of 2 at most 1/4" % len(found14))
expo = math.log2(3) - Fr(5, 4)
report(abs(expo - 0.33496) < 1e-5 and Fr(81, 32) == Fr(3 ** 4, 2 ** 5), "growth exponent log_2 3 - 5/4 = %.5f (note: 0.335); a block (1,1,1,2) has 3^4/2^5 = 81/32" % expo)
h = lambda r: -r * math.log2(r) - (1 - r) * math.log2(1 - r)
print("   with frequency of 2 at most 1/4, the parity word has odd-fraction rho = l/d_l >= 4/5 = 0.8; 1 - h(0.8) = %.4f (THM-4476's counting-lemma exponent; the note's '(1 - h(rho)) with rho >= 0.8')" % (1 - h(0.8)))
# the exploration: rationals in [1,2) with denominator <= 40
def phi(b):
    return sum(1 for a in range(1, b + 1) if math.gcd(a, b) == 1)
samples = sorted({Fr(a, b) for b in range(1, 41) for a in range(b, 2 * b)})
report(len(samples) == 490 == sum(phi(b) for b in range(1, 41)), "rational starting errors in [1,2) with denominators <= 40: %d = sum_(b<=40) phi(b)" % len(samples))
def itin(e0, L):
    e = e0; w = []
    for _ in range(L):
        if e < Fr(5, 3):
            e = (3 * e - 1) / 2; w.append(1)
        else:
            e = (3 * e - 1) / 4; w.append(2)
    return w, e
def prefix_class(w):
    """the class mod 2^(d_K+1) of integers whose first K valuations are w: from 2^(d_K) m_K = 3^K n + S_K with m_K odd"""
    K = len(w); dK = sum(w); S = S_of(w); mod = 2 ** (dK + 1)
    return ((2 ** dK - S) * pow(3, -K, mod)) % mod, mod
# sanity: the prefix classes of the orbit of 27 contain 27, and every odd n <= 999 lies in its own prefix class
K27 = steps_to_one(27, U); _, w27 = orbit(27, U, K27)
pc_ok = all(27 % prefix_class(w27[:K])[1] == prefix_class(w27[:K])[0] for K in range(1, K27 + 1))
pc_ok &= all((lambda r: n % r[1] == r[0])(prefix_class(orbit(n, U, 12)[1])) for n in range(1, 1000, 2))
report(pc_ok, "prefix-class formula rho_K = (2^(d_K) - S_K) 3^(-K) mod 2^(d_K+1) verified on the orbit of 27 (all K) and on all odd n < 1000 (K = 12)")
L = 80
freqs = []; worst = (Fr(0), None); stable = 0; three_ones = True; growth = 0; tail_bound_ok = True
for e0 in samples:
    w, eL = itin(e0, L)
    freqs.append(w.count(2) / L)
    three_ones &= runs_ok(w, cyclic=False)
    d = [0]
    for v in w:
        d.append(d[-1] + v)
    growth += (2 ** d[L] < 3 ** L)
    partial = R_head(w)
    dev = e0 - partial                     # exact: equals (2^(d_L)/3^L) eta_L by the recursion
    tail_bound_ok &= (dev == Fr(2 ** d[L], 3 ** L) * eL) and 0 < dev < 2 * Fr(2 ** d[L], 3 ** L)
    if dev > worst[0]:
        worst = (dev, e0)
    rhos = [prefix_class(w[:K])[0] for K in range(L // 2, L + 1)]
    stable += (len(set(rhos)) == 1)
report(min(freqs) == 0 and "%.3f" % (sum(freqs) / len(freqs)) == "0.152" and "%.3f" % max(freqs) == "0.188" and max(freqs) == Fr(15, 80) and three_ones and growth == 490,
       "490 itineraries of length 80: frequency of 2 min %.3f mean %.4f max %.4f = 15/80 (note: 0, 0.152, 0.188); three-ones rule holds in all; all 490 are growth words" % (min(freqs), sum(freqs) / len(freqs), max(freqs)))
report(abs(float(worst[0]) - 3.761e-10) < 1e-12 and worst[1] == Fr(53, 38) and tail_bound_ok,
       "eta_0 - (partial sum) = (2^(d_L)/3^L) eta_L exactly; max %.3e at eta_0 = %s (note: 3.761e-10 at 53/38); it is < 2 2^(d_L)/3^L <= 2 2^(-0.335 L) = %.1e, NOT < (2/3)^L = %.1e as the itineraries .out prints" % (float(worst[0]), worst[1], 2 * 2 ** (-float(expo) * L), (2 / 3) ** L))
report(stable == 0, "prefix residues rho_K constant over K in [40, 80] (a positive integer point below 2^41 would show up): %d of 490 (note: 0)" % stable)
report(prefix_class([1] * 40)[0] == 2 ** 41 - 1, "eta_0 = 1 (word 1^inf): rho_40 = 2^41 - 1, i.e. the point is -1: an integer, but negative, so it does not stabilise (the test detects positive integers only)")

# ============================================================================= wrap-up
sha_line("summary")
print(" all automated checks passed: %s" % ok_all)
print(" findings that are not automated (see the audit report): the 'bounded-ratio' uniqueness wording (A1),")
print(" the minus-sheet counterexamples to 'no positive integer orbit on either sheet has sup eta < inf' (A5),")
print(" the FLP carry count (A6), the 2-adic point -1 / +1 of the itinerary (1) (A7).")
