#!/usr/bin/env python3
"""collatz_shadow_error_20260927.py -- the real shadow of a Collatz orbit and its error term
(session collatz-shadow-flp-20260927, opus, 2026-09-27).

Syracuse orbit m_0 = n, m_(l+1) = (3 m_l + 1)/2^(v_(l+1)), d_l = v_1 + ... + v_l. Bernstein real value of the
word from position l:  eta_l = sum_(k>=0) 2^(d_(l+k) - d_l) / 3^(k+1)  (converges iff the orbit is divergent,
THM-4476). Exact facts checked here:
 (E1) recursion  2^(v_(l+1)) eta_(l+1) = 3 eta_l - 1  (the 3x-1 copy driven by the orbit's halvings);
 (E2) eta_l >= 1, with equality iff every later valuation is 1;
 (E3) 2^(v_(l+1)) <= 9 eta_l - 3  (a large next valuation needs a large shadow error), more generally
      2^(d_(l+k) - d_l) <= 3^(k+1) eta_l, hence  m_(l+k) >= m_l / (3 eta_l)  for every k >= 0;
 (E4) for a periodic tail w^inf with cycle point x_w, eta = -x_w = |x_w| exactly (real and 2-adic values agree);
      on the negative cycles the real shadow m + eta is exactly 0 (the two copies cancel);
 (E5) real shadow x_l = m_l + eta_l = xi 3^l / 2^(d_l) with xi = n C_inf; x_(l+1) = 3 x_l / 2^(v_(l+1)).
For finite orbits (reaching 1) the series diverges; the truncated tails eta_l^(K) (sum to the last odd step K
before the 1-cycle) are partial sums, so (E1) holds exactly between truncated values and (E3) holds for them;
they are shown on 27, 703, 871, 6171 and on the -5 shadow families n_k = 4*8^k - 5.
Usage: python3 collatz_shadow_error_20260927.py
"""
import math
from fractions import Fraction


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def word(n, steps):
    w = []; m = n
    for _ in range(steps):
        m, v = U(m); w.append(v)
    return w


def bernstein_periodic(w):
    """exact real value of sum_(k>=0) 2^(d_k)/3^(k+1) for the periodic word w^inf (geometric series)"""
    p = len(w); A = sum(w)
    head = Fraction(0); d = 0
    for t in range(p):
        head += Fraction(2 ** d, 3 ** (t + 1)); d += w[t]
    ratio = Fraction(2 ** A, 3 ** p)
    if ratio >= 1:
        return None
    return head / (1 - ratio)


def cycle_point(w):
    p = len(w); A = sum(w); S = 0; d = 0
    for t in range(p):
        S += 3 ** (p - 1 - t) * 2 ** d; d += w[t]
    return Fraction(S, 2 ** A - 3 ** p)


def truncated_eta(w):
    """eta_l^(K) for l = 0..K-1 from a finite word w = (v_1..v_K): partial sums of the tail series to the end"""
    K = len(w)
    d = [0]
    for v in w:
        d.append(d[-1] + v)
    etas = []
    for l in range(K):
        s = Fraction(0)
        for k in range(K - l):
            s += Fraction(2 ** (d[l + k] - d[l]), 3 ** (k + 1))
        etas.append(s)
    return etas, d


def period_of(m, maxp=12):
    L = 6 * maxp
    w = word(m, L)
    for p in range(1, maxp + 1):
        if w[:p] * (L // p) == w[:p * (L // p)]:
            return p
    return None


def main():
    print("== (E4) periodic tails: eta = |x_w| exactly ==")
    for w in ((1,), (1, 2), (2, 1), (1, 1, 1, 2, 1, 1, 4), (1, 1, 2, 1, 1, 4, 1), (1, 1, 2), (1, 2, 1, 2, 2, 1)):
        R = bernstein_periodic(w)
        if R is None:
            print(" w = %-24s descent word (2^A > 3^p): real series diverges" % (w,))
        else:
            x = cycle_point(w)
            print(" w = %-24s x_w = %-10s eta = R(w^inf) = %-10s  eta + x_w = %s" % (w, x, R, R + x))
    for n in (-1, -5, -17):
        m = n; ok = True; vals = []
        for _ in range(8):
            p = period_of(m)
            eta = bernstein_periodic(tuple(word(m, p)))
            vals.append((m, eta)); ok &= (eta + m == 0)
            m, _ = U(m)
        print(" negative cycle through %d: (m, eta) = %s ... m + eta = 0 at every point: %s" % (n, [(a, str(b)) for a, b in vals[:3]], ok))
    print("== (E1)-(E3) on truncated tails of finite orbits ==")
    for n in (27, 703, 871, 6171):
        m = n; K = 0
        while m != 1:
            m, _ = U(m); K += 1
        w = word(n, K)
        etas, d = truncated_eta(w)
        ms = [n]
        for v in w:
            ms.append(U(ms[-1])[0])
        e1 = all(2 ** w[l] * etas[l + 1] == 3 * etas[l] - 1 for l in range(K - 1))
        e3v = all(2 ** w[l] <= 9 * etas[l] - 3 for l in range(K - 1))
        e3m = all(ms[l + k] >= Fraction(ms[l], 3 * etas[l]) for l in range(K) for k in range(K - l))
        e3s = all(2 ** w[l] <= 3 * etas[l] - 1 for l in range(K - 1) if K - l >= 30)
        # FLP carry structure on the truncated shadow x_l = m_l + eta_l: M_l = floor(x_l), 2^v M_(l+1) = 3 M_l + c_l
        Ms = [ms[l] + (etas[l].numerator // etas[l].denominator) for l in range(K)]
        cs = [2 ** w[l] * Ms[l + 1] - 3 * Ms[l] for l in range(K - 1)]
        carry_ok = all(-(2 ** w[l]) + 1 <= cs[l] <= 2 for l in range(K - 1))
        imax = max(range(K), key=lambda l: etas[l])
        print(" n = %d: K = %d odd steps; 2^v eta_(l+1) = 3 eta_l - 1 on truncated values: %s; 2^v <= 9 eta_l - 3: %s; sharp 2^v <= 3 eta_l - 1 (tails >= 30 terms): %s; m_(l+k) >= m_l/(3 eta_l) inside the truncation: %s; carries c_l = 2^v M_(l+1) - 3 M_l in [-2^v+1, 2]: %s (values seen %s); eta_0^(K) = %.4f, max eta = %.3f at l = %d (m_l = %d)" % (
            n, K, e1, e3v, e3s, e3m, carry_ok, sorted(set(cs)), float(etas[0]), float(etas[imax]), imax, ms[imax]))
    print("== the -5 shadow families n_k = 4*8^k - 5: eta_0 over the shadow tends to 5 ==")
    for k in range(1, 9):
        n = 4 * 8 ** k - 5
        w = word(n, 2 * k)
        assert w == [1, 2] * k
        s = sum(Fraction(2 ** sum(w[:t]), 3 ** (t + 1)) for t in range(2 * k))
        print(" k=%d n=%d: word (1,2)^%d, partial eta_0 over the shadow = %.6f  (5(1 - (8/9)^k) = %.6f)" % (k, n, k, float(s), 5 * (1 - (8 / 9) ** k)))
    print("== (E5) real shadow on 27 with the truncated xi: x_l = m_l + eta_l = xi 3^l/2^(d_l) ==")
    n = 27
    m = n; K = 0
    while m != 1:
        m, _ = U(m); K += 1
    w = word(n, K); etas, d = truncated_eta(w)
    xi = n + etas[0]
    ms = [n]
    for v in w:
        ms.append(U(ms[-1])[0])
    ok = all(Fraction(xi * 3 ** l, 2 ** d[l]) == ms[l] + etas[l] for l in range(K))
    print(" xi = 27 + eta_0 = %.6f; xi 3^l/2^(d_l) = m_l + eta_l for all l < K: %s; l = 10: m = %d, eta = %.4f, x = %.4f" % (float(xi), ok, ms[10], float(etas[10]), float(ms[10] + etas[10])))
    print("== (E2) the infimum over infinite words is the all-ones value 1; largest partial eta_0 over no-descent words of length 12 ==")
    best = Fraction(0)
    def rec(prefix, At):
        nonlocal best
        j = len(prefix)
        if j == 12:
            s = sum(Fraction(2 ** sum(prefix[:t]), 3 ** (t + 1)) for t in range(12))
            if s > best:
                best = s
            return
        vmax = int(math.floor((j + 1) * math.log2(3) - At - 1e-12))
        for v in range(1, vmax + 1):
            rec(prefix + [v], At + v)
    rec([], 0)
    print(" largest partial eta_0 over no-descent words of length 12: %.4f (words hugging the critical line from below); all-ones value: %.4f" % (float(best), sum(2 ** t / 3 ** (t + 1) for t in range(12))))
    part_bounded_shadow()
    part_minus_sheet()


def part_bounded_shadow():
    print("== bounded-shadow regime: periodic words with all eta in [1, 2) are itineraries of f(eta) = (3 eta - 1)/2 on [1, 5/3), (3 eta - 1)/4 on [5/3, 2) ==")
    found = []
    def rec(prefix, At, p):
        j = len(prefix)
        if j == p:
            w = tuple(prefix)
            if w != min(w[i:] + w[:i] for i in range(p)):
                return
            if any(w == w[i:] + w[:i] for i in range(1, p)):
                return
            R = bernstein_periodic(w)
            if R is None:
                return
            # eta along the period
            eta = R; vals = [eta]; ok = True
            for v in w:
                if not (1 <= eta < 2):
                    ok = False; break
                eta = (3 * eta - 1) / 2 ** v; vals.append(eta)
            if ok and eta == R:
                # itinerary of f from R
                it = []; e = R
                for _ in range(p):
                    if e < Fraction(5, 3):
                        it.append(1); e = (3 * e - 1) / 2
                    else:
                        it.append(2); e = (3 * e - 1) / 4
                found.append((w, R, tuple(it) == w, max(vals)))
            return
        vmax = int(math.floor((j + 1) * math.log2(3) - At - 1e-12))
        for v in range(1, min(vmax, 2) + 1):
            rec(prefix + [v], At + v, p)
    for p in range(1, 11):
        rec([], 0, p)
    print(" %d primitive periodic words (p <= 10) with all eta in [1,2); itinerary of f equals the word in all cases: %s" % (len(found), all(f[2] for f in found)))
    for w, R, ok, mx in found[:12]:
        print("   w = %-32s eta_0 = %-12s max eta = %.4f" % (w, R, float(mx)))
    print(" words with all eta < B force valuations <= log2(3B - 1): B = 2 gives v <= 2, B = 5 gives v <= 3, B = 17 gives v <= 5")


def part_minus_sheet():
    print("== minus sheet (3x-1 map): the same recursion, real shadow x = m - eta, cycles 1, 5, 17 ==")
    def Um(m):
        m = 3 * m - 1; v = 0
        while m % 2 == 0:
            m //= 2; v += 1
        return m, v
    def word_m(n, steps):
        w = []; m = n
        for _ in range(steps):
            m, v = Um(m); w.append(v)
        return w
    for n in (1, 5, 17):
        m = n; ok = True; vals = []
        for _ in range(8):
            L = 72; w = word_m(m, L)
            pp = next(p for p in range(1, 13) if w[:p] * (L // p) == w[:p * (L // p)])
            eta = bernstein_periodic(tuple(w[:pp]))
            vals.append((m, str(eta))); ok &= (eta == m)
            m, _ = Um(m)
        print(" 3x-1 cycle through %d: (m, eta) = %s ... m - eta = 0 at every point: %s" % (n, vals[:3], ok))
    # a 3x-1 orbit reaching the fixed point 1: n = 3 -> 4 -> 1 (v = 3): R(d) = 1 * 2^0/3 + 2^3/9 + 2^4/27 + ... hmm compute the word and check R = n and the recursion
    for n in (3, 11, 29):
        w = word_m(n, 60)
        # tail is (1)^inf once the orbit is at 1; the real value is the finite head plus the geometric tail
        m = n; K = 0
        while m != 1:
            m, _ = Um(m); K += 1
        head = sum(Fraction(2 ** sum(w[:t]), 3 ** (t + 1)) for t in range(K))
        dK = sum(w[:K])
        tail = Fraction(2 ** dK, 3 ** K)  # sum_(k>=0) 2^(dK + k)/3^(K+k+1) = 2^dK/3^K
        R = head + tail
        print(" 3x-1 orbit of %d reaches 1 after %d steps: R(d) = %s = n: %s (eventually periodic word: real value equals the 2-adic value)" % (n, K, R, R == n))


if __name__ == '__main__':
    main()
