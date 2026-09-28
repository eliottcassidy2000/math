#!/usr/bin/env python3
"""collatz_two_carries_typology_20260927.py -- D24, D25, D26 and the typology of arithmetic iterations
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, tenth note).

 (1) D24, pointwise memorylessness. For m = 2^(K+1) t - 1 (odd, t >= 1) the Syracuse orbit follows v = 1 for K steps,
     U^K(m) = 2 * 3^K t - 1, and the next valuation is 1 + v_2(3^(K+1) t - 1): the valuation after the -1 shadow is a
     function of the free cofactor t, and regeneration of a depth-J valuation needs t = 3^(-(K+1)) mod 2^(J-1), i.e.
     m >= 2^(K+J) up to a unit: the price sheet with the cofactor explicit. The two-carries identities: Collatz
     d_j = v_2(3^j n + S_(j-1)) with v_2(3^j n) = 0 (all 2-adic content is created by the carry), aliquot
     v_2(s(2^a m)) = min(a, v_2(sigma(m))) unless equal (the carry -n caps the content of sigma at the driver).
 (2) D25, average-case persistence of the growth event on both sides: Collatz P(v_(j+1) = 1 | v_j = 1) along orbits
     (Terras: 1/2); aliquot P(s_(k+1)(n) > s_k(n) for all k < K | s(n) > n) on a sample of abundant n <= 10^6 (Erdos 1976:
     the growth ratio persists for almost all abundant n over any fixed number of steps).
 (3) D26, in-degrees: aliquot preimage counts of n <= 632 (exact) and the touched density; Syracuse: infinite in-degree
     off the multiples of 3 (density 1/3 of leaves), descent tree in-degree c_D from the parallel note (cited).
 (4) Typology: drift and memory for Collatz T_3, T_5, aliquot (cited from the ninth note), Juggler (n -> floor(n^(3/2))
     odd, floor(sqrt n) even) and reverse-and-add (Lychrel), with the provable invariant-driven iterations (Ducci,
     Kaprekar, look-and-say) as controls.
Usage: python3 collatz_two_carries_typology_20260927.py
"""
import math, random, sys

sys.setrecursionlimit(10000)


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def v2(x):
    c = 0
    while x % 2 == 0 and x:
        x //= 2; c += 1
    return c


def part1():
    print("== (1) D24: the valuation after a -1 shadow is a function of the free cofactor ==")
    ok = True
    for K in range(1, 12):
        for t in range(1, 60):
            m = 2 ** (K + 1) * t - 1
            x = m; word = []
            for _ in range(K):
                x, v = U(x); word.append(v)
            ok &= (word == [1] * K) and (x == 2 * 3 ** K * t - 1)
            _, vnext = U(x)
            ok &= (vnext == 1 + v2(3 ** (K + 1) * t - 1))
    print(" for m = 2^(K+1) t - 1 (K <= 11, t < 60): word (1)^K, U^K(m) = 2 * 3^K t - 1, next valuation = 1 + v_2(3^(K+1) t - 1): %s" % ok)
    # regeneration: next valuation >= J iff t = 3^(-(K+1)) mod 2^(J-1)
    ok2 = True
    for K in range(1, 8):
        for J in range(2, 8):
            inv = pow(3, -(K + 1), 2 ** (J - 1))
            for t in range(1, 2 ** (J - 1) * 4):
                m = 2 ** (K + 1) * t - 1
                x = m
                for _ in range(K):
                    x, _v = U(x)
                _, vnext = U(x)
                ok2 &= ((vnext >= J) == (t % 2 ** (J - 1) == inv))
    print(" next valuation >= J iff t = 3^(-(K+1)) mod 2^(J-1) (K < 8, J < 8): %s -> the least such m is 2^(K+1) rho_J - 1 with rho_J = 3^(-(K+1)) mod 2^(J-1): regeneration is priced by the cofactor" % ok2)
    # free regeneration (t = 1, the Mersenne number 2^(K+1) - 1): next valuation = 2 (K even) or 3 + v_2(K+1) (K odd), by LTE
    ok5 = True
    for K in range(1, 200):
        m = 2 ** (K + 1) - 1; x = m
        for _ in range(K):
            x, _v = U(x)
        _, vnext = U(x)
        ok5 &= (vnext == (2 if K % 2 == 0 else 3 + v2(K + 1)))
    print(" free regeneration after the run of the Mersenne number 2^(K+1) - 1 (cofactor t = 1): next valuation = 2 for K even, 3 + v_2(K+1) for K odd (K < 200): %s -- logarithmic in the run length" % ok5)
    # two carries
    ok3 = True
    for n in range(1, 4001, 2):
        m = n; S = 0; d = 0
        for j in range(1, 40):
            S = 3 * S + 2 ** d
            m, v = U(m); d += v
            A = 3 ** j * n + S
            ok3 &= (v2(A) == d) and (v2(3 ** j * n) == 0)
            if m == 1:
                break
    print(" Collatz: d_j = v_2(3^j n + S_(j-1)) with v_2(3^j n) = 0 (odd n < 4000, all steps to 1): %s -- the 2-adic content is created by the carry" % ok3)
    # aliquot cap
    def sigma(n):
        s = 1; p = 2; nn = n
        while p * p <= nn:
            if nn % p == 0:
                e = 0
                while nn % p == 0:
                    nn //= p; e += 1
                s *= (p ** (e + 1) - 1) // (p - 1)
            p += 1
        if nn > 1:
            s *= nn + 1
        return s
    ok4 = True; equal_cases = 0
    for n in range(2, 60001, 2):
        a = v2(n); m = n >> a
        sm = sigma(m); vs = v2(sm)
        s = sigma(n) - n
        if vs != a:
            ok4 &= (v2(s) == min(a, vs))
        else:
            equal_cases += 1
            ok4 &= (v2(s) >= a + 1)
    print(" aliquot: v_2(s(2^a m)) = min(a, v_2 sigma(m)) when they differ, >= a+1 when equal (even n <= 60000, %d equal cases): %s -- the carry -n caps the content" % (equal_cases, ok4))


def part2():
    print("== (2) D25: persistence of the growth event, both sides ==")
    # Collatz: P(v_(j+1) = 1 | v_j = 1) along orbits of odd n < 2*10^5, and runs
    both = 0; first = 0
    for n in range(3, 2 * 10 ** 5, 2):
        m = n; prev = None
        while m != 1:
            m, v = U(m)
            if prev == 1:
                first += 1; both += (v == 1)
            prev = v
    print(" Collatz: P(v_(j+1) = 1 | v_j = 1) = %.4f (Terras: 0.5; the growth event does not persist)" % (both / first))
    # aliquot: sample abundant n <= 10^6, follow K steps with sympy
    try:
        from sympy import factorint
    except ImportError:
        print(" sympy not available"); return
    def sig(n):
        s = 1
        for p, e in factorint(n).items():
            p = int(p); e = int(e)
            s *= (p ** (e + 1) - 1) // (p - 1)
        return int(s)
    random.seed(4)
    K = 6; runs = [0] * (K + 1); tot = 0
    while tot < 3000:
        n = random.randrange(2, 10 ** 6)
        s1 = sig(n) - n
        if s1 <= n:
            continue
        tot += 1
        prev = n; cur = s1; k = 1
        while k < K:
            nxt = sig(cur) - cur
            if nxt > cur:
                k += 1; cur = nxt
            else:
                break
        for j in range(1, k + 1):
            runs[j] += 1
    print(" aliquot: among %d random abundant n <= 10^6, P(s_2 > s_1 > n) = %.3f, P(s_3 > s_2 > s_1 > n) = %.3f, ..., P(six increases) = %.3f" % (tot, runs[2] / tot, runs[3] / tot, runs[6] / tot))
    print(" (Erdos 1976: for almost all n with s(n) > n the ratio s_(k+1)/s_k stays close to s(n)/n for any fixed number of steps: the growth event persists)")


def part3():
    print("== (3) D26: in-degrees ==")
    N = 4 * 10 ** 5
    sig = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            sig[m] += d
    exact = math.isqrt(N)
    indeg = {}
    for m in range(2, N + 1):
        s = sig[m] - m
        if 1 < s <= exact:
            indeg[s] = indeg.get(s, 0) + 1
    dist = {}
    for n in range(2, exact + 1):
        k = indeg.get(n, 0); dist[k] = dist.get(k, 0) + 1
    print(" aliquot: in-degree distribution of n in [2, %d] (exact; preimages m <= n^2 scanned): %s; mean in-degree %.3f" % (exact, dict(sorted(dist.items())), sum(indeg.get(n, 0) for n in range(2, exact + 1)) / (exact - 1)))
    print(" Syracuse: in-degree infinite for m not divisible by 3 (predecessors (2^v m - 1)/3), zero for multiples of 3 (density 1/3); descent tree in-degree c_D in [1.67, 1.70] (parallel note)")


def juggler_stats():
    drift = 0.0; steps = 0; same = 0; pairs = 0; maxlen = 0; maxpeak = 0
    for n in range(2, 1501):
        x = n; prev = None; L = 0
        while x > 1 and L < 400:
            nxt = math.isqrt(x ** 3) if x % 2 == 1 else math.isqrt(x)
            if nxt > 1 and x > 1:
                drift += math.log2(math.log2(nxt) / math.log2(x)); steps += 1
            par = x % 2
            if prev is not None:
                pairs += 1; same += (par == prev)
            prev = par; x = nxt; L += 1
            maxpeak = max(maxpeak, x)
        maxlen = max(maxlen, L)
    return drift / steps, same / pairs, maxlen, len(str(maxpeak))


def lychrel_stats():
    N = 10 ** 5
    drift = 0.0; steps = 0; nonterm = []
    for n in range(10, N):
        x = n; k = 0; done = False
        while k < 200:
            s = str(x)
            if s == s[::-1] and k > 0:
                done = True; break
            y = x + int(s[::-1])
            drift += math.log10(y) - math.log10(x); steps += 1
            x = y; k += 1
        if not done:
            nonterm.append(n)
    return drift / steps, nonterm[:12], len(nonterm)


def part4():
    print("== (4) typology: invariant, drift, memory ==")
    jd, jsame, jlen, jpk = juggler_stats()
    print(" Juggler n -> floor(n^(3/2)) (odd) / floor(sqrt n) (even), n <= 1500: a walk on log log n; mean log_2 of the ratio of log-sizes per step %+.3f (fair coin on (3/2, 1/2) would give %+.3f), P(parity persists) = %.3f (memoryless would be 0.5); longest run %d steps, largest peak %d digits" % (jd, (math.log2(1.5) + math.log2(0.5)) / 2, jsame, jlen, jpk))
    ld, lyc, cnt = lychrel_stats()
    print(" reverse-and-add, n < 10^5, 200 steps: mean digit growth per step %+.3f (doubling would be log10 2 = 0.301); starts not palindromic within 200 steps: %d, first %s (196 is the smallest Lychrel candidate; the delays below 10^5 are at most 55 steps, so these are the candidates)" % (ld, cnt, lyc))
    print(" table (drift sign, memory of the driving quantity, exact invariant?, conjectured fate, status):")
    rows = [
        ("Collatz 3n+1", "negative (-0.415 bits/odd step)", "memoryless (Terras)", "none", "all terminate", "OPEN; Terras density 1, Tao almost-bounded"),
        ("3n-1", "negative", "memoryless", "none", "three cycles, all bounded", "OPEN (bounded conj.)"),
        ("5n+1", "positive (+0.16)", "memoryless", "none", "most diverge", "OPEN; no divergence proved"),
        ("aliquot s(n)", "negative in class 1, positive in classes >= 2", "sticky (drivers)", "none", "many diverge (Guy-Selfridge); Lehmer five", "OPEN; Catalan-Dickson vs Guy-Selfridge"),
        ("Juggler", "negative (-0.21 bits/step)", "memoryless (parity persists ~0.5)", "none", "all terminate", "OPEN"),
        ("reverse-and-add", "positive (+0.3 digits/step)", "carries", "none", "Lychrel numbers never palindromic", "OPEN (196)"),
        ("Ducci (length 2^k)", "-", "-", "linear over GF(2), nilpotent", "terminates", "PROVED"),
        ("Kaprekar (fixed digit count)", "-", "-", "finite state space", "terminates", "PROVED"),
        ("look-and-say", "positive (Conway's constant 1.3036)", "-", "linear over 92 atoms", "diverges", "PROVED"),
    ]
    for r in rows:
        print("   ", r)
    print(" reading: provability tracks the presence of an exact invariant (finite state, linear algebra); among the invariant-free maps the conjectured fate tracks the sign of the drift and the memory of the driving quantity")


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
