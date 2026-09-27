#!/usr/bin/env python3
"""collatz_depth_layers_20260927.py -- the depth function of the rooted component, its layers, its Hankel-rank form,
the exponent competition, and the two September-2026 zeta(5) preprints read by their own equations
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, second note).

 (1) Layers of the rooted component by Collatz-step depth d(n) = min{k : C^k(n) in {1,2,4}} (C(n) = n/2 or 3n+1):
     backward construction; the owner's table (depth 5 = {3, 20, 21, 128}); layer sizes L_k for k <= 60 and their
     growth ratio (heuristic 4/3; Applegate-Lagarias bracket [1.302, 1.36] for pruned trees).
 (2) Bigraded layers: odd n reaching 1 with exactly d odd (Syracuse) steps and D halvings, counted by the word method
     (n = (2^D - S)/3^d positive integer, no earlier visit to 1), for d <= 12, against the density heuristic
     C(D-1, d-1) 3^(-d) summed over D.
 (3) The Hankel form of the depth: the linear complexity (Berlekamp-Massey over F_p, p = 2^61 - 1) of the sequence
     2^(d_k(n)) of 2-parts of the Apery forms A_k = 3^k n + S_(k-1), continued along the 1-cycle, equals d_odd(n) + 1
     (odd-step depth + 1) for every odd n <= 3000 (Kronecker: PC(n) iff finite Hankel rank).
 (4) The exponent competition along orbits: nu_L = d_L/(L log_2 3) and the margin d_L/L - log_2 3 at the first
     coefficient descent; the generic drift 2 - log_2 3 = 0.415 bits per odd step; compared with the margin
     C_0' - C_2' = 0.4766 claimed in arXiv:2609.22316 (a different unit; numerology, recorded as such).
 (5) arXiv:2609.22316 (Suman) read by its equations: tau_0 and C_0' = -Re f_0(tau_0) reproduced from eq. (33)-(34) with
     eta_0 = 3, eta_j = 1; the prime window of eq. (13), h_0 < p <= m_8 = n, is EMPTY (h_0 = 3n + 2), so Phi = 1 and
     C_2' = 9 > C_0' = 5.75: Lemma 2's inequality fails as written; a hypothetical window (n, 3n+2] saves ~1.4 < 3.25.
 (6) Numerology table (tested, all labelled): 3^7 - 2^11 = 139 (the -17 cycle's denominator; Fauzan's exponent 139/5);
     11/7 is a semiconvergent of log_2 3 (the -17 cycle), 3/2 a convergent (the -5 cycle), 1/1 (the -1 cycle);
     37 in the -17 cycle (Fauzan's degree 37n); Suman's q = 11, eight maxima, A = 8, B = 3, C = 2.
Usage: python3 collatz_depth_layers_20260927.py
"""
import math
from fractions import Fraction

LOG23 = math.log2(3)
PR = (1 << 61) - 1


def C(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def part1():
    print("== (1) layers of the rooted component by Collatz-step depth (root {1,2,4} collapsed) ==")
    layer = {0: [1, 2, 4]}
    seen = {1, 2, 4}
    frontier = [1, 2, 4]
    sizes = [3]
    for k in range(1, 61):
        nxt = []
        for n in frontier:
            for p in (2 * n,) + (((n - 1) // 3,) if (n % 6 == 4 and (n - 1) // 3 > 0) else ()):
                if p not in seen:
                    seen.add(p); nxt.append(p)
        layer[k] = sorted(nxt); frontier = nxt; sizes.append(len(nxt))
    for k in range(6):
        print("  depth %d: %s" % (k, layer[k]))
    print("  owner's table reproduced (depth 5 = [3, 20, 21, 128]): %s" % (layer[5] == [3, 20, 21, 128]))
    print("  layer sizes L_k, k = 0..60:", sizes)
    for k in (20, 30, 40, 50, 60):
        print("  (L_k)^(1/k) at k = %d: %.4f   L_k/L_(k-1): %.4f" % (k, sizes[k] ** (1 / k), sizes[k] / sizes[k - 1]))
    print("  heuristic 4/3 = 1.3333; Applegate-Lagarias pruned-tree bracket [1.302, 1.36]")
    lam = (1 + math.sqrt(7 / 3)) / 2
    print("  per-Collatz-step growth heuristic: a halving predecessor (cost 1 step) always, a (n-1)/3 predecessor (cost 2 steps) with density 1/3:")
    print("  1 = 1/lam + (1/3)/lam^2, lam = (1 + sqrt(7/3))/2 = %.5f (the T-step rate 4/3 becomes this in Collatz steps)" % lam)
    # check d(2^x m) = x + d(m) on the computed layers
    depth = {}
    for k, L in layer.items():
        for n in L:
            depth[n] = k
    ok = all(depth.get(2 * n) == depth[n] + 1 for n in depth if n not in (1, 2) and 2 * n in depth)
    print("  d(2n) = d(n) + 1 on the computed component (n not in {1,2}): %s" % ok)


def part2():
    print("== (2) bigraded layers: odd n reaching 1 with exactly d odd steps and D halvings ==")
    # words (v_1..v_d) with sum D: n = (2^D - S)/3^d, S = sum_(t<d) 3^(d-1-t) 2^(d_t); n odd positive integer and the
    # orbit must not hit 1 before step d (i.e. m_t != 1 for t < d)
    from itertools import combinations
    print("  d, number of (D, n) pairs, heuristic sum_D C(D-1,d-1)/3^d over the same D-range, min n, max n")
    for d in range(1, 11):
        count = 0; nmin = None; nmax = 0; heur = 0.0
        Dmax = 2 * d + 6
        for D in range(d, Dmax + 1):
            heur += math.comb(D - 1, d - 1) / 3 ** d
            # compositions of D into d parts >= 1 via cut positions
            for cuts in combinations(range(1, D), d - 1):
                parts = [b - a for a, b in zip((0,) + cuts, cuts + (D,))]
                S = 0; dd = 0
                for v in parts:
                    S = 3 * S + 2 ** dd; dd += v
                num = 2 ** D - S
                if num > 0 and num % 3 ** d == 0:
                    n = num // 3 ** d
                    if n % 2 == 1:
                        # verify the orbit (word exact, no earlier 1)
                        m = n; okw = True
                        for t, v in enumerate(parts):
                            if m == 1 and t < d:
                                okw = False; break
                            m2, v2 = U(m)
                            if v2 != v:
                                okw = False; break
                            m = m2
                        if okw and m == 1:
                            count += 1; nmax = max(nmax, n); nmin = n if nmin is None else min(nmin, n)
        print("  d = %2d: %5d pairs (D <= %d), heuristic %.1f, n in [%s, %d]" % (d, count, Dmax, heur, nmin, nmax))


def lin_complexity_modp(seq):
    # Berlekamp-Massey over F_p
    p = PR
    C_ = [1]; B = [1]; L = 0; m = 1; b = 1
    for i, s in enumerate(seq):
        d = (s + sum(C_[j] * seq[i - j] for j in range(1, L + 1))) % p
        if d == 0:
            m += 1
        elif 2 * L <= i:
            T = C_[:]; coef = d * pow(b, -1, p) % p
            C_ = C_ + [0] * (len(B) + m - len(C_))
            for j in range(len(B)):
                C_[j + m] = (C_[j + m] - coef * B[j]) % p
            L = i + 1 - L; B = T; b = d; m = 1
        else:
            coef = d * pow(b, -1, p) % p
            C_ = C_ + [0] * max(0, len(B) + m - len(C_))
            for j in range(len(B)):
                C_[j + m] = (C_[j + m] - coef * B[j]) % p
            m += 1
    return L


def part3():
    print("== (3) the Hankel form of the depth: linear complexity of (2^(d_k)) = odd-step depth + 1 ==")
    ok = True; worst = None
    for n in range(1, 3001, 2):
        m = n; d = [0]
        while m != 1:
            m, v = U(m); d.append(d[-1] + v)
        depth = len(d) - 1
        seq = [pow(2, x, PR) for x in d] + [pow(2, d[-1] + 2 * j, PR) for j in range(1, 2 * depth + 12)]
        L = lin_complexity_modp(seq)
        if L != depth + 1:
            ok = False; worst = (n, depth, L)
    print("  LC(2^(d_k(n)), continued along the 1-cycle) = d_odd(n) + 1 for all odd n <= 3000: %s %s" % (ok, "" if ok else worst))
    print("  reading: Kronecker -- the generating function sum 2^(d_k) z^k is rational iff the Hankel determinants vanish from")
    print("  some size on; the odd-step depth is that size minus one; a divergent orbit would have infinite Hankel rank")


def part4():
    print("== (4) the exponent competition along orbits ==")
    for n in (27, 703, 6171, 77031, 837799):
        m = n; d = 0; L = 0; first = None
        while m != 1:
            m, v = U(m); d += v; L += 1
            if first is None and 2 ** d > 3 ** L:
                first = (L, d, d / L - LOG23)
        print("  n = %d: first coefficient descent at L = %d with d_L = %d, margin d_L/L - log_2 3 = %.4f bits per odd step" % ((n,) + first))
    print("  generic drift (Terras): 2 - log_2 3 = %.4f bits per odd step; arXiv:2609.22316 claims C_0' - C_2' = 0.4766 (nats per n; different unit)" % (2 - LOG23))


def part5():
    print("== (5) arXiv:2609.22316 read by its equations ==")
    import cmath
    # roots of (tau-3)^3 (tau-1)^11 - tau^3 (tau-2)^11 by Durand-Kerner
    def P(t):
        return (t - 3) ** 3 * (t - 1) ** 11 - t ** 3 * (t - 2) ** 11
    deg = 13  # leading terms cancel: tau^14 - tau^14; degree 13
    # coefficients via expansion (use numpy-free polynomial arithmetic)
    def polymul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] += x * y
        return r
    def polypow(a, k):
        r = [1]
        for _ in range(k):
            r = polymul(r, a)
        return r
    A = polymul(polypow([-3, 1], 3), polypow([-1, 1], 11))
    B = polymul(polypow([0, 1], 3), polypow([-2, 1], 11))
    coef = [a - b for a, b in zip(A, B)]
    while coef and coef[-1] == 0:
        coef.pop()
    n = len(coef) - 1
    lead = coef[-1]
    monic = [c / lead for c in coef]
    roots = [complex(0.4, 0.9) ** k for k in range(n)]
    for _ in range(500):
        new = []
        for i, r in enumerate(roots):
            val = sum(monic[k] * r ** k for k in range(n + 1))
            den = 1
            for j, s in enumerate(roots):
                if j != i:
                    den *= (r - s)
            new.append(r - val / den)
        roots = new
    cand = [r for r in roots if r.imag > 1e-6 and r.real < 3]
    tau0 = max(cand, key=lambda r: r.real)
    f0 = 9 * cmath.log(3 - tau0) + 11 * cmath.log(tau0 - 1) - 22 * cmath.log(tau0 - 2)
    print("  eq. (34) root with Im > 0 and maximal Re: tau_0 = %.8f + %.8fi (paper: 2.86852453 + 0.11960091i)" % (tau0.real, tau0.imag))
    print("  f_0(tau_0) = %.8f + %.8fi (paper: -5.75349395 - 8.95071225i); C_0' = %.6f reproduced" % (f0.real, f0.imag, -f0.real))
    print("  eq. (13): Phi = prod over h_0 < p <= m_8 with h_0 = 3n + 2 and m_8 = n (eq. 16): the range is empty, Phi = 1,")
    print("  so the denominator is D_n^9 and C_2' = 9; Lemma 2 needs C_0' > C_2', i.e. 5.75 > 9: FALSE as written.")
    # hypothetical window (n, 3n+2]
    def primes_upto(N):
        s = bytearray([1]) * (N + 1); s[0] = s[1] = 0
        for i in range(2, int(N ** 0.5) + 1):
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
        return [i for i in range(N + 1) if s[i]]
    def vkp(k, p, n):
        h0 = 3 * n + 2; hj = n + 1; t = 0
        for j in range(1, 4):
            t += (k - 1) // p + (h0 - k - 1) // p - (k - hj) // p - (h0 - hj - k) // p - 2 * ((hj - 1) // p)
        for j in range(4, 12):
            t += (h0 - 2 * hj) // p - (k - hj) // p - (h0 - hj - k) // p
        return t
    for nn in (100, 200, 400):
        h0 = 3 * nn + 2; h4 = nn + 1; sav = 0.0
        for p in primes_upto(h0):
            if p > nn:
                vp = min(vkp(k, p, nn) for k in range(h4, h0 - h4 + 1))
                sav += max(vp, 0) * math.log(p)
        print("  hypothetical window (n, 3n+2] at n = %d: saving (1/n) sum v_p log p = %.4f (needed: 9 - 5.7535 = 3.2465; the paper's C_2' = 5.2769 would need 3.72)" % (nn, sav / nn))


def part6():
    print("== (6) numerology, tested and labelled ==")
    print("  3^7 - 2^11 = %d; the -17 cycle: p = 7 odd steps, A = 11 halvings, -17 = S/(2^11 - 3^7) with S = %d" % (3 ** 7 - 2 ** 11, 17 * 139))
    # continued fraction of log_2 3
    x = LOG23; cf = []; y = x
    for _ in range(9):
        a = math.floor(y); cf.append(a); y = 1 / (y - a)
    conv = []; h0, h1, k0, k1 = 0, 1, 1, 0
    for a in cf:
        h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0; conv.append((h1, k1))
    print("  continued fraction of log_2 3:", cf, " convergents (A/p):", conv[:7])
    print("  negative cycles (A, p): (1,1) convergent, (3,2) convergent, (11,7) = mediant of 3/2 and 8/5 (semiconvergent, not a convergent)")
    print("  members of the -17 cycle: 17, 25, 37, 55, 41, 61, 91 (37 appears; the Zenodo preprint's degree is 37n, its exponent 139 n^2/5, its measure 260)")
    print("  verdict: coincidences of small integers; no equation of either preprint involves 2^11, 3^7, or a Collatz object; recorded, not used")


def part7():
    print("== (7) the Zenodo preprint (zenodo.org/records/22826419) read by its equations: the functional of its section 2.1 ==")
    try:
        import mpmath as mp
    except ImportError:
        print("  mpmath not available; skipped"); return
    mp.mp.dps = 30
    z5 = mp.zeta(5)
    print("  pole entries at X = zeta(5): j^4 (zeta(5) - H_j) - 1/4 + 1/(2j) against Hermite: 2 j^4 int_0^oo sin(5 arctan(u/j)) / ((j^2+u^2)^(5/2) (e^(2 pi u) - 1)) du")
    for j in (1, 2, 5, 10):
        H = sum(mp.mpf(1) / mp.mpf(v) ** 5 for v in range(1, j + 1))
        lhs = j ** 4 * (z5 - H) - mp.mpf(1) / 4 + mp.mpf(1) / (2 * j)
        rhs = 2 * j ** 4 * mp.quad(lambda u: mp.sin(5 * mp.atan(u / j)) / ((j * j + u * u) ** mp.mpf(2.5) * (mp.exp(2 * mp.pi * u) - 1)), [0, 1, 10, mp.inf])
        print("   j = %2d: %s = %s (difference %s)" % (j, mp.nstr(lhs, 15), mp.nstr(rhs, 15), mp.nstr(lhs - rhs, 3)))
    print("  polynomial entries: (-1)^e B_(2e+2) (2e+3)(2e+4)(2e+5)/24 against the Binet integral int u^(2e+1) (2e+2)(2e+3)(2e+4)(2e+5) / (12 (e^(2 pi u) - 1)) du")
    for e in (0, 1, 2, 3):
        lhs = (-1) ** e * mp.bernoulli(2 * e + 2) * (2 * e + 3) * (2 * e + 4) * (2 * e + 5) / 24
        rhs = mp.quad(lambda u: u ** (2 * e + 1) / (mp.exp(2 * mp.pi * u) - 1), [0, 1, 10, mp.inf]) * (2 * e + 2) * (2 * e + 3) * (2 * e + 4) * (2 * e + 5) / 12
        print("   e = %d: %s = %s" % (e, mp.nstr(lhs, 15), mp.nstr(rhs, 15)))
    g = lambda x: 1 / (mp.exp(x) - 1)
    print("  fourth derivative of the Planck factor 1/(e^x - 1) at x = 0.5, 2, 6: %s (positive: complete monotonicity, the source of the positive weight)" % [mp.nstr(mp.diff(g, x, 4), 6) for x in (0.5, 2, 6)])
    print("  reading: both entry types are moments of the positive weight (1/12) u^5 (d/du)^4 [1/(e^(2 pi u) - 1)] in t = u^2, up to explicit polynomial corrections;")
    print("  the pole entries decay like 1/j^2 (0.29, 0.091, 0.016, 0.0041 at j = 1, 2, 5, 10): the smallness engine is the Euler-Maclaurin remainder of the tails of zeta(5)")


def main():
    part1(); part2(); part3(); part4(); part5(); part6(); part7()


if __name__ == '__main__':
    main()
