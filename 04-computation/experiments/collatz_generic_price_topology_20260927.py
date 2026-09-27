#!/usr/bin/env python3
"""collatz_generic_price_topology_20260927.py -- D14 (pricing the generic words) and the topological reframe
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, sixth note).

 (1) The empirical price of generic growth: for odd n < 10^6, peak(n) = max of the Syracuse orbit, net growth
     h(n) = log_2(peak/n); the cheapness ratio log_2 n / h(n) (bits of source per bit of net growth), its record lows,
     and the record excursions' exponent log(peak)/log(n) (the Lagarias-Weiss heuristic says peak <= n^(2+o(1))).
     Against the cycle-shadow price 1/c_(-1) = 1.71 bits of source per bit of growth.
 (2) Why non-integer points of E_inf carry no size price: for a 2-adic point x, the least positive integer in the
     depth-K approach class x + 2^K Z is rho_K(x) = x mod 2^K, so the price pi_K(x) = log_2 rho_K(x) = K - z_K(x) where
     z_K(x) is the run of zero bits of x just below position K. For the negative integer cycle points z_K = 0 (price K);
     for the trivial cycle the price is 0; for the Sturmian point of E_inf and for fixed points of random no-descent
     words the deficit K - pi_K is a geometric-like fluctuation, and an integer m approaches x to depth K > log_2 m + 1
     iff x's bits between log_2 m + 1 and K - 1 are zero (x looks like m up to bit K): the deep approaches of small
     integers are to the points that look like them, i.e. the integer's own word -- D14 is circular below the size.
 (3) Fixed points of the conjugacy map Q (parity vector as a 2-adic number): N_k = #{x mod 2^k : Q(x) = x mod 2^k} for
     k <= 22, the solutions that lift through every level (the true 2-adic fixed points seen to depth 22), the
     inclusion {0} u {-2^j} subset Fix(Q) (the doubling chain of the fixed point -1 of T), and the 2-cycle {1, -1/3}.
 (4) Necklace splitting of a cycle word: the doubled -17 word (14 beads, types 1, 2, 4) splits by <= 3 cuts into two
     shares with equal type counts (Alon's bound k(q-1) = 3): an explicit split, and the remark that the shares are
     not orbit segments.
Usage: python3 collatz_generic_price_topology_20260927.py
"""
import math, itertools

LOG23 = math.log2(3)


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def part1():
    print("== (1) the empirical price of generic growth (odd n < 10^6) ==")
    N = 10 ** 6
    best_ratio = []; rec_peak = 0; records = []
    for n in range(3, N, 2):
        m = n; peak = n
        while m != 1:
            m, _ = U(m)
            if m > peak:
                peak = m
        if peak > rec_peak:
            rec_peak = peak
            records.append((n, peak))
        if peak > n:
            h = math.log2(peak / n)
            if h >= 4:
                ratio = math.log2(n) / h
                best_ratio.append((ratio, n, round(h, 2)))
    best_ratio.sort()
    print(" cheapest generic growth (>= 4 bits of net growth), ratio log_2 n / log_2(peak/n):", [(round(r, 3), n, h) for r, n, h in best_ratio[:8]])
    print(" record excursions (n, peak, log(peak)/log(n), bits of source per bit of net growth):")
    for n, peak in records[-14:]:
        print("   n = %8d  peak = %16d  exponent %.3f  price %.3f" % (n, peak, math.log(peak) / math.log(n), math.log2(n) / math.log2(peak / n)))
    print(" reading: the cheapest generic growth costs about 0.5-1 bit of source per bit of growth (record excursions approach 1 as")
    print(" peak ~ n^2), against 1/0.585 = 1.71 bits per bit for the -1 shadow: generic words are cheaper than any cycle shadow")


def residue_of_word(word_d, K):
    # 2-adic point with cumulative valuations d_0 = 0, d_1, ..., modulo 2^K: x = -sum_(k) 2^(d_k)/3^(k+1) (first note, Prop. 9)
    mod = 2 ** K
    s = 0
    for k, d in enumerate(word_d):
        if d >= K + 2:
            break
        s = (s + pow(2, d, mod) * pow(3, -(k + 1), mod)) % mod
    return (-s) % mod


def part2():
    print("== (2) approach classes and the price pi_K(x) = log_2(least positive member of x + 2^K Z) = K - z_K(x) ==")
    # Sturmian point: d_k = ceil(k log_2 3)
    dS = [0] + [math.ceil(k * LOG23) for k in range(1, 400)]
    rows = []
    for K in (10, 20, 30, 40, 60, 80, 120, 160):
        rho = residue_of_word(dS, K)
        rows.append((K, round(math.log2(rho), 2), K - math.floor(math.log2(rho)) - 1))
    print(" Sturmian point x_S of E_inf: (K, pi_K, zero run z_K):", rows)
    for x in (-1, -5, -17):
        print(" negative integer cycle point %d: rho_K = 2^K - %d, pi_K = K - O(2^-K); trivial cycle 1: rho_K = 1, pi_K = 0" % (x, -x) if x == -1 else " negative integer cycle point %d: rho_K = 2^K - %d" % (x, -x))
    # random no-descent words: fixed points x_w of w^inf, deficits K - pi_K
    import random
    random.seed(11)
    deficits = {}
    for trial in range(2000):
        # random no-descent Syracuse word of length 30 with sum d < 30 log_2 3 (coefficient no-descent at the end)
        while True:
            w = [1 if random.random() < 0.7 else 2 for _ in range(30)]
            d = [0]
            for v in w:
                d.append(d[-1] + v)
            if all(3 ** j > 2 ** d[j] for j in range(1, 31)):
                break
        K = 40
        # x = fixed point of w^inf: word d extends periodically
        dd = [0]; t = 0
        while dd[-1] < K + 40:
            dd.append(dd[-1] + w[t % 30]); t += 1
        rho = residue_of_word(dd, K)
        z = K - math.floor(math.log2(rho)) - 1 if rho else K
        deficits[z] = deficits.get(z, 0) + 1
    tot = sum(deficits.values())
    print(" fixed points of 2000 random no-descent words, K = 40: P(deficit K - pi_K >= z) for z = 0..6:", [round(sum(c for zz, c in deficits.items() if zz >= z) / tot, 3) for z in range(7)], "(uniform bits would give 2^-z)")
    # the tautology: an integer m approaches x to depth K > log_2 m + 1 iff x = m mod 2^K
    print(" an integer m lies in the depth-K class of x with K > log_2 m + 1 iff x = m mod 2^K, i.e. the bits of x from log_2 m + 1 to K - 1 vanish:")
    print(" the only 2-adic points a small integer approaches deeply are the points that look like it, so the price of a generic word")
    print(" is the integer's own word; the size price exists exactly for the negative integer cycle points, whose classes have no small members")


def Qk(x, k):
    # parity vector of the T-orbit of x (T(x) = x/2 or (3x+1)/2) as an integer mod 2^k; x taken as an integer representative
    q = 0
    for i in range(k):
        b = x & 1
        q |= b << i
        x = (3 * x + 1) // 2 if b else x // 2
    return q


def part3():
    print("== (3) fixed points of the conjugacy map Q ==")
    kmax = 22
    sols = {}
    for k in range(1, kmax + 1):
        mod = 1 << k
        S = [x for x in range(mod) if Qk(x, k) % mod == x]
        sols[k] = S
    print(" N_k = #{x mod 2^k : Q(x) = x mod 2^k} for k = 1..%d:" % kmax, [len(sols[k]) for k in range(1, kmax + 1)])
    # true fixed points seen to depth kmax: solutions mod 2^kmax whose projections are solutions at every level (automatic) -- list them
    S = sols[kmax]; mod = 1 << kmax
    named = []
    for x in S:
        if x == 0:
            named.append("0")
        elif (mod - x) & (mod - x - 1) == 0:
            named.append("-2^%d" % ((mod - x).bit_length() - 1))
        else:
            named.append(str(x))
    print(" solutions mod 2^%d: %s" % (kmax, named))
    # inclusion {0} u {-2^j}: Q(-2^j) = -2^j (word 0^j 1^inf)
    ok = all(Qk((mod - (1 << j)) % mod, kmax) % mod == (mod - (1 << j)) % mod for j in range(kmax))
    print(" Q(-2^j) = -2^j for j < %d (the doubling chain of the T-fixed point -1): %s" % (kmax, ok))
    # the 2-cycle {1, -1/3}
    inv3 = pow(3, -1, mod); m13 = (-inv3) % mod
    print(" Q(1) = -1/3 and Q(-1/3) = 1 mod 2^%d: %s" % (kmax, Qk(1, kmax) % mod == m13 and Qk(m13, kmax) % mod == 1))
    # which solutions mod 2^k lift to mod 2^(k+1)? count spurious ones
    spur = []
    for k in range(1, kmax):
        lifted = {x % (1 << k) for x in sols[k + 1]}
        spur.append(len([x for x in sols[k] if x not in lifted]))
    print(" solutions mod 2^k that do not lift to 2^(k+1), k = 1..%d: %s" % (kmax - 1, spur))


def part4():
    print("== (4) necklace splitting of the doubled -17 cycle word ==")
    w = (1, 1, 1, 2, 1, 1, 4) * 2
    types = sorted(set(w)); n = len(w)
    target = {t: w.count(t) // 2 for t in types}
    found = None
    for cuts in itertools.combinations(range(1, n), 3):
        pieces = []; prev = 0
        for c in list(cuts) + [n]:
            pieces.append(w[prev:c]); prev = c
        for assign in itertools.product((0, 1), repeat=len(pieces)):
            share = [p for p, a in zip(pieces, assign) if a == 0]
            cnt = {t: sum(p.count(t) for p in share) for t in types}
            if cnt == target:
                found = (cuts, assign, share); break
        if found:
            break
    cuts, assign, share = found
    print(" word %s: 14 beads, 3 types; a split with cuts at %s and pieces assigned %s gives thief A the pieces %s" % (w, cuts, assign, share))
    print(" both thieves get %s; the shares are unions of intervals, not orbit segments, so they carry no dynamics (Alon's bound k(q-1) = 3 cuts attained here with 3)" % target)


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
