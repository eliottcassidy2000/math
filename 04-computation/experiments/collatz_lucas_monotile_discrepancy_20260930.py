#!/usr/bin/env python3
"""collatz_lucas_monotile_discrepancy_20260930.py -- the owner's Lucas identities, the Fibonacci-monodromy monotile's
arithmetic, the discrepancy paper's Hadamard rigidity, and the Collatz cycles (session collatz-posets-zeta5-20260927,
opus, 2026-09-30, eleventh note).

 (1) phi^(2n) + (-1)^n = L_n phi^n (the owner's phi^8 + 1 = 7 phi^4 and phi^10 = 1 + 11 phi^5 are n = 4, 5); the norm
     N(phi^j - 1) = (-1)^j - L_j + 1, so |Z[phi]/(phi^j - 1)| = |L_j - 1 - (-1)^j| (the monotile paper's Theorem E(ii)):
     5, 11 at j = 4, 5 and 121 = 11^2 at j = 10 -- the Brocard roots 5, 11 and the square 5! + 1; the field criterion
     (Z[phi]/(phi^p - 1) a field iff L_p prime) tabulated.
 (2) The shared prefix of two continued fractions: log_2 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...], phi = [1; 1, 1, 1, ...];
     shared convergents 1/1, 2/1, 3/2 and the common convergent 8/5; the negative Collatz cycles' (A, p) = (1,1), (3,2), (11,7)
     are F_2/F_1, F_4/F_3 and L_5/L_4 = mediant(3/2, 8/5) = (F_4 + F_6)/(F_3 + F_5); Gersonides: |2^A - 3^p| = 1 exactly for the
     three shortest cycles (-1: (1,1); +1: (2,1); -5: (3,2)) and 139 (prime) for -17: the 'field case' Z/139.
 (3) The tower's skew-Hadamard matrices S = M + I (THM-447 doubling from the Paley heptagon's and from H_2): S S^T = nI, so every
     signing has ||Sx||_oo >= sqrt(n) (Parseval); exact discrepancy min_x ||Sx||_oo for orders 2, 4, 8, 16 (and Sylvester H_n),
     heuristic for 32; the Sylvester attainment disc(H_(2^m)) = sqrt(n) at even m via H_4 y = 2y reproduced.
 (4) The cycle ring: for each rotation class of the negative cycle words, S_w mod (2^A - 3^p) = 0 (integrality), and for the
     Lucas quadratics x^2 - L_n x - (-1)^n the roots phi^n, psi^n; the 'unit^j - 1' torsion L_j - 1 - (-1)^j against the
     Collatz 'torsion' 3^p - 2^A along the convergents of log_2 3.
Usage: python3 collatz_lucas_monotile_discrepancy_20260930.py
"""
import math, itertools
from fractions import Fraction

PHI = (1 + 5 ** 0.5) / 2


def lucas(n):
    a, b = 2, 1
    if n == 0:
        return 2
    for _ in range(n - 1):
        a, b = b, a + b
    return b


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def part1():
    print("== (1) Lucas identities, norms in Z[phi], the monotile's torsion orders ==")
    for n in (4, 5):
        lhs = PHI ** (2 * n) + (-1) ** n; rhs = lucas(n) * PHI ** n
        print(" phi^%d + (-1)^%d = %.10f, L_%d phi^%d = %.10f -> phi^(2n) + (-1)^n = L_n phi^n with L_%d = %d" % (2 * n, n, lhs, n, n, rhs, n, lucas(n)))
    # exact: phi^n = (L_n + F_n sqrt5)/2; phi^(2n) + (-1)^n - L_n phi^n = 0 in Q(sqrt5)
    def phi_pow(n):  # (a, b) with phi^n = (a + b sqrt5)/2
        return Fraction(lucas(n), 2), Fraction(fib(n), 2)
    ok = True
    for n in range(1, 30):
        a2, b2 = phi_pow(2 * n); a1, b1 = phi_pow(n)
        ok &= (a2 + (-1) ** n - lucas(n) * a1 == 0) and (b2 - lucas(n) * b1 == 0)
    print(" exact in Q(sqrt5) for n < 30: %s" % ok)
    print(" |Z[phi]/(phi^j - 1)| = |N(phi^j - 1)| = |L_j - 1 - (-1)^j| (Theorem E(ii) of the monotile paper):")
    rows = []
    for j in range(2, 21):
        order = abs(lucas(j) - 1 - (-1) ** j)
        field = (j in (3, 4)) or (j >= 5 and is_prime(j) and is_prime(lucas(j)))
        rows.append((j, lucas(j), order, "field" if field else ""))
    print("  (j, L_j, order, field?):", rows)
    print(" j = 4 -> 5 (F_5), j = 5 -> 11 (F_11), j = 10 -> 121 = 11^2 = 5! + 1: the Brocard roots of the eighth note; L_4 = 7 and L_5 = 11 are the owner's coefficients")


def contfrac(x, k):
    cf = []
    for _ in range(k):
        a = math.floor(x); cf.append(a); x = 1 / (x - a)
    return cf


def convergents(cf):
    h0, h1, k0, k1 = 0, 1, 1, 0
    out = []
    for a in cf:
        h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0; out.append((h1, k1))
    return out


def part2():
    print("== (2) the shared prefix of the continued fractions of log_2 3 and phi, and the negative cycles ==")
    cf3 = contfrac(math.log2(3), 10); cfphi = [1] * 10
    print(" log_2 3 = %s..., phi = %s..." % (cf3, cfphi))
    c3 = convergents(cf3); cphi = convergents(cfphi)
    print(" convergents of log_2 3:", c3[:7]); print(" convergents of phi:    ", cphi[:7])
    shared = [c for c in c3 if c in cphi]
    print(" common convergents:", shared, "(1/1, 2/1, 3/2 from the shared prefix [1; 1, 1]; 8/5 = F_6/F_5 is a convergent of both although log_2 3 skips 5/3)")
    cycles = {"-1": (1, 1, (1,)), "+1": (2, 1, (2,)), "-5": (3, 2, (1, 2)), "-17": (11, 7, (1, 1, 1, 2, 1, 1, 4))}
    for name, (A, p, w) in cycles.items():
        S = 0; d = 0
        for v in w:
            S = 3 * S + 2 ** d; d += v
        den = 2 ** A - 3 ** p
        print(" cycle %s: (A, p) = (%d, %d), A/p = %s, 2^A - 3^p = %d (%s), x = S/(2^A - 3^p) = %d/%d = %s" % (name, A, p, Fraction(A, p), den, "unit" if abs(den) == 1 else ("prime" if is_prime(abs(den)) else "composite"), S, den, Fraction(S, den)))
    print(" (1,1) = F_2/F_1, (3,2) = F_4/F_3, (11,7) = L_5/L_4 = (F_4 + F_6)/(F_3 + F_5) = mediant(3/2, 8/5): %s" % (Fraction(11, 7) == Fraction(fib(4) + fib(6), fib(3) + fib(5)) and lucas(5) == 11 and lucas(4) == 7))
    print(" Gersonides (1343): 3^p - 2^A = +-1 only for (p, A) in {(1,1), (1,2), (2,3)} -- checked for p, A <= 60: %s" % ([(p, A) for p in range(1, 61) for A in range(1, 61) if abs(3 ** p - 2 ** A) == 1] == [(1, 1), (1, 2), (2, 3)]))
    print(" reading: the golden numbers in the Collatz cycle data are the shared initial segment of two continued fractions; the -17 cycle's (p, A) = (L_4, L_5) = (7, 11) is the Lucas mediant one step past the shared segment, with 3^7 - 2^11 = 139 prime")


def skew_double(M):
    n = len(M)
    D = [[0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            D[i][j] = M[i][j]
            D[i][j + n] = M[i][j] + (1 if i == j else 0)
            D[i + n][j] = M[i][j] - (1 if i == j else 0)
            D[i + n][j + n] = -M[i][j]
    return D


def sylvester(m):
    H = [[1]]
    for _ in range(m):
        n = len(H)
        H = [[H[i % n][j % n] * (-1 if (i >= n and j >= n) else 1) for j in range(2 * n)] for i in range(2 * n)]
    return H


def disc_exact(A):
    n = len(A); best = None
    for bits in range(1 << (n - 1)):   # x_0 = +1 by symmetry
        x = [1] + [1 if (bits >> k) & 1 else -1 for k in range(n - 1)]
        val = max(abs(sum(A[i][j] * x[j] for j in range(n))) for i in range(n))
        if best is None or val < best:
            best = val
    return best


def disc_heuristic(A, trials=4000, seed=1):
    import random
    random.seed(seed); n = len(A); best = None
    for _ in range(trials):
        x = [random.choice((1, -1)) for _ in range(n)]
        # local search: flip a coordinate if it reduces the max
        improved = True
        while improved:
            improved = False
            cur = max(abs(sum(A[i][j] * x[j] for j in range(n))) for i in range(n))
            for k in range(n):
                x[k] = -x[k]
                val = max(abs(sum(A[i][j] * x[j] for j in range(n))) for i in range(n))
                if val < cur:
                    cur = val; improved = True
                else:
                    x[k] = -x[k]
        if best is None or cur < best:
            best = cur
    return best


def part3():
    print("== (3) discrepancy of the tower's skew-Hadamard matrices against Sylvester (Reis-Song, arXiv:2609.34471) ==")
    # tower from H_2: M_2 = [[0,1],[-1,0]] dominance of the 2-vertex tournament; S = M + I
    M = [[0, 1], [-1, 0]]
    for k in range(1, 5):
        n = len(M)
        S = [[M[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
        gram_ok = all(sum(S[i][t] * S[j][t] for t in range(n)) == (n if i == j else 0) for i in range(n) for j in range(n))
        d = disc_exact(S) if n <= 16 else disc_heuristic(S)
        print(" skew tower order %2d: S S^T = nI: %s; disc(S) = %d (%s), sqrt(n) = %.3f, ratio %.3f" % (n, gram_ok, d, "exact" if n <= 16 else "local search", math.sqrt(n), d / math.sqrt(n)))
        M = skew_double(M)
    for m in range(1, 6):
        H = sylvester(m); n = len(H)
        d = disc_exact(H) if n <= 16 else disc_heuristic(H)
        print(" Sylvester H_%d (m = %d): disc = %d (%s), sqrt(n) = %.3f; paper: = sqrt(n) at even m via H_4 y = 2y, y = (1,1,1,-1)" % (n, m, d, "exact" if n <= 16 else "local search", math.sqrt(n)))
    print(" reading: the L^2 identity S S^T = nI forces every signing to have a coordinate of size sqrt(n) -- a rigidity valid for all 2^n signings at once,")
    print(" which is the kind of universal statement the Collatz thread lacks for orbits: Terras's bijection says every word occurs, so no residue is forced to descend")


def part4():
    print("== (4) the cycle ring and the Lucas quadratics ==")
    for name, w in (("-5", (1, 2)), ("-17", (1, 1, 1, 2, 1, 1, 4))):
        p = len(w); A = sum(w); den = 2 ** A - 3 ** p
        pts = []
        for r in range(p):
            wr = w[r:] + w[:r]; S = 0; d = 0
            for v in wr:
                S = 3 * S + 2 ** d; d += v
            pts.append((S % abs(den) == 0, S // den))
        print(" cycle %s: all %d rotations have S_w = 0 mod |2^A - 3^p| = %d: %s; points %s" % (name, p, abs(den), all(t for t, _ in pts), [x for _, x in pts]))
    for n in (3, 4, 5, 6):
        print(" Lucas quadratic x^2 - %d x %s 1 = 0 (x^2 - L_n x + (-1)^n) has roots phi^%d, psi^%d: %s" % (lucas(n), "-" if n % 2 else "+", n, n, abs((PHI ** n) ** 2 - lucas(n) * PHI ** n + (-1) ** n) < 1e-9))
    print(" torsion of the monotile tower: |L_j - 1 - (-1)^j| for j = 3, 6, 9, 12: %s (orders L_(3n) or L_(3n) - 2); Collatz cycle 'torsion' 3^p - 2^A along the convergents of log_2 3: %s" % ([abs(lucas(j) - 1 - (-1) ** j) for j in (3, 6, 9, 12)], [(A, p, 3 ** p - 2 ** A) for A, p in ((1, 1), (2, 1), (3, 2), (8, 5), (19, 12), (65, 41))]))


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
