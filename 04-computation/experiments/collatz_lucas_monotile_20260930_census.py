#!/usr/bin/env python3
"""collatz_lucas_monotile_20260930_census.py -- the torsion-point reading of the Collatz cycle equation, made exact
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, eleventh note; companion of collatz_lucas_monotile_discrepancy_20260930.py).

 (5) Torsion-point census. For a shape (A, p) (A halvings, p odd steps) the cycle map is x -> (3^p x + S_w)/2^A with fixed
     point S_w/(2^A - 3^p); the class of S_w in Z/(2^A - 3^p) is a torsion point of the group of fixed points of the linear
     part (the solenoid analogue of the monotile paper's O_j = Fix(Q^j) on the torus); an integer cycle of that shape exists
     iff some word of the shape has S_w = 0 in that group.  Counted exactly for every shape with A <= 22 (direct enumeration
     of the C(A-1, p-1) words) and by a residue DP for the near-critical shapes with 23 <= A <= 27 and |2^A - 3^p| <= 10^7.
     Expected hits: Gersonides (1,1), (2,1), (3,2), the sporadic (11,7), and their multiples (repeated cycles).
 (6) The uniform-residue heuristic (words / (p |clock|)), its no-descent refinement, h(log_3 2) against the repo value
     0.9499555272, and the heuristic tail along the near-critical shapes.
 (7) Prime occurrence (the analogue of Proposition 7.4 of the monotile paper): a prime l >= 5 divides 2^A - 3^p iff (A, p)
     lies in the kernel lattice of (A, p) -> 2^A 3^(-p) in F_l^x, of index |<2, 3>|; tabulated, with l = 139 and (11, 7).
 (8) The cycle census over Z for the shortcut map T (|x| <= 2 * 10^5) and the conjectural Artin-Mazur zeta function.
 (9) Cayley-Hamilton for Q^4 and Q^5 (the owner's identities), the associated Mersenne numbers |O_j|, the doubling chain
     |O_(2^k)| = 5 F_(2^(k-1))^2, the resultant reading of the three clocks, the spectra of the skew tower and of Sylvester,
     and Zsigmondy along the cycle rays.
Usage: python3 collatz_lucas_monotile_20260930_census.py
"""
import math, time
from fractions import Fraction
from itertools import combinations
from math import comb
import numpy as np

LOG23 = math.log2(3)


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def zero_words_direct(A, p, D, want_residues=False):
    """words (v_1..v_p), v_i >= 1, sum A  <->  cut sets {d_1 < ... < d_(p-1)} in {1..A-1};
    S_w = sum_{i=0}^{p-1} 3^(p-1-i) 2^(d_i) with d_0 = 0.  Returns (#words with S_w = 0 mod |D|, #distinct residues)."""
    D = abs(D)
    base = pow(3, p - 1, D)
    tab = [[pow(3, p - 1 - i, D) * pow(2, d, D) % D for d in range(A)] for i in range(p)]
    zeros = 0; res = set()
    for cuts in combinations(range(1, A), p - 1):
        s = base
        for i, d in enumerate(cuts, start=1):
            s += tab[i][d]
        s %= D
        if s == 0:
            zeros += 1
        if want_residues:
            res.add(s)
    return zeros, (len(res) if want_residues else None)


def zero_words_dp(A, p, D):
    """the same count by a DP over residues mod |D| (3 is invertible mod D, so each odd step permutes the residues)."""
    D = abs(D)
    idx = np.arange(D, dtype=np.int64)
    dp = [None] * (A + 1); dp[0] = np.zeros(D, dtype=np.int64); dp[0][0] = 1
    for j in range(p):
        new = [None] * (A + 1); acc = None
        for d in range(A + 1):
            if acc is not None:
                new[d] = acc.copy()
            if dp[d] is not None and j < p:
                perm = (3 * idx + pow(2, d, D)) % D
                vec = np.zeros(D, dtype=np.int64); vec[perm] = dp[d]
                acc = vec if acc is None else acc + vec
        dp = new
    return int(dp[A][0]) if dp[A] is not None else 0


def nodescent_words(A, p):
    """words of shape (A, p) with d_j <= floor(j log_2 3) for all j < p (never below the start, at the word level)."""
    dp = {0: 1}
    for j in range(1, p + 1):
        cap = A if j == p else int(math.floor(j * LOG23))
        new = {}
        for d, c in dp.items():
            for v in range(1, cap - d + 1):
                new[d + v] = new.get(d + v, 0) + c
        dp = new
    return dp.get(A, 0)


def part5_6():
    print("== (5) torsion-point census: words of shape (A, p) with S_w = 0 in Z/(2^A - 3^p) ==")
    t0 = time.time()
    # cross-check of the two counters on small shapes
    bad = [(A, p) for A in range(1, 13) for p in range(1, A + 1) if 2 ** A != 3 ** p and zero_words_direct(A, p, 2 ** A - 3 ** p)[0] != zero_words_dp(A, p, 2 ** A - 3 ** p)]
    print(" DP vs direct enumeration agree on all shapes with A <= 12: %s" % (bad == []))
    hits = []; expect = Fraction(0); nd_sum = 0.0; shapes = 0; coverage = []
    for A in range(1, 23):
        for p in range(1, A + 1):
            D = 2 ** A - 3 ** p
            shapes += 1
            z, cov = zero_words_direct(A, p, D, want_residues=(abs(D) <= 2000))
            words = comb(A - 1, p - 1); nd = nodescent_words(A, p)
            if abs(D) <= 2000 and abs(D) > 1:
                coverage.append((A, p, abs(D), cov, words))
            gers = (A, p) in ((1, 1), (2, 1), (3, 2))
            if not gers:
                expect += Fraction(words, p * abs(D)); nd_sum += nd / abs(D)
            if z:
                hits.append((A, p, D, z, words, nd))
    for A, p in ((23, 14), (23, 15), (24, 15), (25, 16), (27, 17)):
        D = 2 ** A - 3 ** p; shapes += 1
        z = zero_words_dp(A, p, D); words = comb(A - 1, p - 1); nd = nodescent_words(A, p)
        expect += Fraction(words, p * abs(D)); nd_sum += nd / abs(D)
        print("  DP shape (%d,%d): clock %d, zero-words %d of %d (no-descent %d)" % (A, p, D, z, words, nd))
        if z:
            hits.append((A, p, D, z, words, nd))
    print(" shapes examined: %d (every shape with A <= 22; five near-critical shapes with 23 <= A <= 27, |clock| <= 10^7); time %.0fs" % (shapes, time.time() - t0))
    print(" hits (A, p, clock, zero-residue words, all words, no-descent words):")
    prim = set()
    for A, p, D, z, w, nd in hits:
        g = math.gcd(A, p); a, b = A // g, p // g
        tag = "primitive" if g == 1 else "repeat of (%d,%d) x%d" % (a, b, g)
        prim.add((a, b))
        print("   (%2d,%2d) clock %9d zero-words %3d of %6d (no-descent %5d)  %s" % (A, p, D, z, w, nd, tag))
    print(" primitive hit shapes: %s -> the cycles -1, +1, -5, -17 and nothing else in range" % sorted(prim))
    print(" coverage of the torsion group by the words of the shape, |clock| <= 2000 (A, p, |clock|, residues hit, words):")
    print("  ", [c for c in coverage if c[2] > 1][:24])
    print("== (6) the uniform-residue heuristic ==")
    print(" sum over examined non-Gersonides shapes of words/(p |clock|) = %.4f; actual primitive hits beyond Gersonides: 1 (the -17 cycle; its own term is 210/973 = %.3f)" % (float(expect), 210 / 973))
    print(" the no-descent refinement, sum of no-descent words/|clock| = %.4f" % nd_sum)
    rho = math.log(2) / math.log(3)
    h = -(rho * math.log2(rho) + (1 - rho) * math.log2(1 - rho))
    print(" h(log_3 2) = %.10f (repo: 0.9499555272), 1 - h = %.10f = the codimension of the no-descent set (sibling dimension ladder, Thm 1a) = the decay exponent of the heuristic per halving" % (h, 1 - h))
    tail = 0.0; top = []
    for A in range(23, 61):
        p0 = round(A / LOG23)
        for q in (p0 - 1, p0, p0 + 1):
            if 1 <= q <= A:
                D = abs(2 ** A - 3 ** q); t = nodescent_words(A, q) / D; tail += t; top.append((t, A, q))
    top.sort(reverse=True)
    print(" heuristic tail over the three nearest-critical shapes for 23 <= A <= 60: %.5f; largest terms %s (the convergents 65/41 and 84/53 of log_2 3 dominate)" % (tail, [(round(t, 5), A, q) for t, A, q in top[:4]]))


def part7():
    print("== (7) prime occurrence: l | 2^A - 3^p iff (A, p) lies in the kernel lattice of (A,p) -> 2^A 3^(-p) in F_l^x ==")
    def order(a, l):
        k = 1; x = a % l
        while x != 1:
            x = x * a % l; k += 1
        return k
    rows = []; full = 0; cnt = 0
    for l in primes_upto(2000):
        if l < 5:
            continue
        o2, o3 = order(2, l), order(3, l)
        sub = set(); x = 1
        for i in range(o2):
            y = x
            for j in range(o3):
                sub.add(y); y = y * 3 % l
            x = x * 2 % l
        size = len(sub); index = (l - 1) // size
        cnt += 1; full += (index == 1)
        if l <= 60 or l == 139:
            rows.append((l, o2, o3, size, index))
    print(" (l, ord 2, ord 3, |<2,3>|, index) for l <= 60 and l = 139:", rows)
    print(" fraction of primes 5 <= l < 2000 with <2,3> = F_l^x: %d/%d = %.4f (two-generator Artin density; the GRH-conditional constant is CITED from memory as Matthews 1976, unverified here)" % (full, cnt, full / cnt))
    l = 139
    lat = sorted((A, p) for A in range(0, 30) for p in range(0, 30) if pow(2, A, l) == pow(3, p, l))
    print(" l = 139: kernel lattice points with A, p < 30: %s" % lat[:16])
    print(" (11,7) in the lattice: %s; 2^11 mod 139 = %d, 3^7 mod 139 = %d; index of the lattice = |<2,3>| = %d = 139 - 1: both 2 and 3 generate F_139^x? ord2 = %d, ord3 = %d" % ((11, 7) in lat, pow(2, 11, 139), pow(3, 7, 139), (139 - 1) // ((139 - 1) // len({pow(2, a, l) * pow(3, b, l) % l for a in range(138) for b in range(138)})), order(2, 139), order(3, 139)))


def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def part8():
    print("== (8) cycle census over Z for the shortcut map T (|x| <= 2 * 10^5) and the conjectural zeta function ==")
    N = 2 * 10 ** 5
    seen = {}; cycles = []
    for x0 in range(-N, N + 1):
        x = x0; path = []
        while abs(x) <= 4 * N and x not in seen:
            seen[x] = x0; path.append(x); x = T(x)
        if x in path:
            i = path.index(x); cyc = path[i:]
            m = min(cyc, key=lambda t: (abs(t), t))
            cycles.append((m, len(cyc), sum(1 for t in cyc if t % 2)))
    cycles = sorted(set(cycles), key=lambda c: (abs(c[0]), c[0]))
    print(" cycles (element of least modulus, T-period = halvings A, odd steps p):", cycles)
    per = sorted(c[1] for c in cycles)
    print(" Artin-Mazur zeta of T on Z, conjecturally exact: zeta(z) = 1 / %s" % " ".join("(1 - z^%d)" % k for k in per))
    print(" T-periods %s: 1, 1, 2, 3 are F_1..F_4 and 11 = L_5 (golden numerology, explained by the shared prefix [1;1,1] of log_2 3 and phi and the Lucas mediant 11/7)" % per)


def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def part9():
    print("== (9) Cayley-Hamilton, associated Mersenne numbers, the resultant reading, spectra, Zsigmondy ==")
    Q = np.array([[0, 1], [1, 1]], dtype=object)
    I = np.array([[1, 0], [0, 1]], dtype=object)
    def mpow(M, k):
        R = I.copy()
        for _ in range(k):
            R = R.dot(M)
        return R
    for n in (4, 5):
        Qn = mpow(Q, n); ch = mpow(Q, 2 * n) - lucas(n) * Qn + (-1) ** n * I
        det = int(Qn[0][0] * Qn[1][1] - Qn[0][1] * Qn[1][0])
        print(" Q^%d = %s, trace %d = L_%d, det %d; Q^%d - L_%d Q^%d + (-1)^%d I = 0: %s  (the owner: phi^%d %s 1 = %d phi^%d)" % (
            n, Qn.tolist(), int(Qn[0][0] + Qn[1][1]), n, det, 2 * n, n, n, n, bool((ch == 0).all()), 2 * n, "+" if n % 2 == 0 else "-", lucas(n), n))
    am = [abs(lucas(j) - 1 - (-1) ** j) for j in range(0, 17)]
    print(" associated Mersenne numbers |O_j| = |L_j - 1 - (-1)^j| = |det(Q^j - I)|, j = 0..16 (OEIS A001350):", am)
    print(" doubling chain j = 2^k: |O_j| = 5 F_(j/2)^2:", [(j, abs(lucas(j) - 2), 5 * fib(j // 2) ** 2) for j in (4, 8, 16, 32)], "; traces L_(2^k) = 3, 7, 47, 2207 (the historian chain), L_(2n) = L_n^2 - 2 for even n")
    from sympy import symbols, resultant, factorint
    t = symbols("t")
    r1 = resultant(t ** 7 - 1, t - 2, t); r2 = resultant(t ** 5 - 1, t ** 2 - t - 1, t); r3 = resultant(t ** 3 - 1, 4 * t - 3, t)
    print(" resultants: Res(t^7 - 1, t - 2) = %s = 2^7 - 1; Res(t^5 - 1, t^2 - t - 1) = %s = |O_5| = L_5; Res(t^3 - 1, 4t - 3) = %s = 2^6 - 3^3 (constant-valuation ray v = 2, p = 3)" % (r1, r2, r3))
    M = np.array([[0, 1], [-1, 0]], dtype=float)
    for k in range(4):
        n = len(M); S = M + np.eye(n)
        ev = np.linalg.eigvals(S)
        print(" skew tower order %2d: eigenvalues of S = M + I: %s (1 +- i sqrt(n-1): non-real, no real eigenvector attains sqrt(n))" % (n, sorted(set(np.round(ev, 6).tolist()), key=lambda z: (z.real, z.imag))[:2]))
        D = np.zeros((2 * n, 2 * n))
        D[:n, :n] = M; D[:n, n:] = M + np.eye(n); D[n:, :n] = M - np.eye(n); D[n:, n:] = -M
        M = D
    H = np.array([[1.0]])
    for m in range(1, 5):
        H = np.block([[H, H], [H, -H]])
        ev = np.linalg.eigvals(H)
        print(" Sylvester H_%d: eigenvalues %s (+- sqrt(n), real; multiplicities %d / %d): Reis-Song attain sqrt(n) at even m via H_4 y = 2 y" % (len(H), sorted(set(np.round(ev.real, 6).tolist())), int(np.sum(ev.real > 0)), int(np.sum(ev.real < 0))))
    for (a, b) in ((1, 1), (2, 1), (3, 2), (11, 7)):
        seen = set(); out = []
        kmax = 12 if a <= 3 else 5
        for k in range(1, kmax + 1):
            nval = abs(2 ** (k * a) - 3 ** (k * b))
            fs = {int(q): int(e) for q, e in factorint(nval).items()} if nval > 1 else {}
            new = [q for q in fs if q not in seen]
            seen |= set(fs)
            out.append((k, "primitive %s" % new if new else "none"))
        print(" ray (%d,%d), primitive prime divisors of 2^(ka) - 3^(kb) (Zsigmondy: present for every k >= 2 here):" % (a, b), out)


if __name__ == "__main__":
    part5_6(); part7(); part8(); part9()
