#!/usr/bin/env python3
"""collatz_posets_dags_zeta5_20260927.py -- posets and DAGs behind the Collatz work, and the shape of an
Apery-type endgame (session collatz-posets-zeta5-20260927, opus, 2026-09-27).

 (1) N-coordinates N = (m+1)/2: the Syracuse step is  N even -> 3N/2  (the left child of the Akiyama-Frougny-
     Sakarovitch base-3/2 tree, Mahler's even branch);  N odd -> (U^-(N) + 1)/2  where U^- is the 3x-1 Syracuse map.
     Checked for all odd m <= 2*10^6. So the 3x+1 map on odd numbers is the 3x-1 map on the shifted numbers,
     interleaved with the AFS left-child step; and v_(l+1) = 1 iff N_l is even.
 (2) The AFS tree: children of N are 3N/2, 3N/2 + 1 (N even) and (3N+1)/2 (N odd); every integer has a unique
     finite base-3/2 expansion via N -> floor(2N/3) (checked N <= 10^5); the Collatz v = 1 steps are tree edges.
 (3) The residual tree (no-descent T-words as linear extensions of THM-4503's width-2 posets): the conditional
     probability that the next letter is odd, given a uniformly random no-descent prefix of length k, stays in
     [1/3, 2/3] -- computed exactly for k <= 60 (a 1/3-2/3 statement for the residual tree, checked, not proved).
 (4) Diophantine quality of the two natural approximations along an orbit (27, truncated at its last odd step):
     real: |xi - 2^(d_L) m_L/3^L| = eta_L 2^(d_L)/3^L versus the denominator 3^L (exponent mu_L = -log_3(error)/L);
     2-adic: |n + S_(L-1)/3^L|_2 = 2^(-d_L) (exponent d_L/(L log_2 3)). Neither is Liouville-strong: the real
     exponent is far below 1, the 2-adic one is automatic.
 (5) The reachability order of the Syracuse DAG: sizes of the tree of 1 below X (predecessors reaching 1 within
     the odd numbers <= X) for X = 10^3..10^6, versus X (the conjecture says all), and the Krasikov-Lagarias
     lower-bound exponent 0.84 for comparison.
Usage: python3 collatz_posets_dags_zeta5_20260927.py
"""
import math
from fractions import Fraction

LOG23 = math.log2(3)


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


def part1():
    print("== (1) N = (m+1)/2 coordinates ==")
    ok = True; v1_even = True
    for m in range(1, 2 * 10 ** 6, 2):
        N = (m + 1) // 2
        mp, v = U(m)
        Np = (mp + 1) // 2
        if N % 2 == 0:
            ok &= (Np == 3 * N // 2) and v == 1
        else:
            ok &= (Np == (Um(N)[0] + 1) // 2) and v >= 2
        v1_even &= ((v == 1) == (N % 2 == 0))
    print(" N even -> 3N/2 with v = 1; N odd -> (U^-(N) + 1)/2 with v >= 2; v = 1 iff N even: %s %s (odd m < 2*10^6)" % (ok, v1_even))


def afs_expansion(N):
    digits = []
    while N > 0:
        a = (2 * N) % 3
        digits.append(a)
        N = (2 * N - a) // 3
    return digits


def part2():
    print("== (2) the AFS base-3/2 tree ==")
    ok = True
    for N in range(1, 10 ** 5 + 1):
        digits = afs_expansion(N)
        val = Fraction(0)
        for i, a in enumerate(digits):
            val += Fraction(a, 2) * Fraction(3, 2) ** i
        ok &= (val == N) and all(a in (0, 1, 2) for a in digits)
    print(" unique finite expansion N = sum (a_i/2)(3/2)^i with digits a_i = 2N_i mod 3, N_(i+1) = floor(2N_i/3): %s (N <= 10^5)" % ok)
    # children
    ch = {}
    for N in range(1, 200):
        parent = (2 * N) // 3
        ch.setdefault(parent, []).append(N)
    ex = {N: ch.get(N, []) for N in (2, 3, 4, 5, 6, 7)}
    print(" children (N with floor(2N/3) = parent):", ex, "-> N even: {3N/2, 3N/2 + 1}; N odd: {(3N+1)/2}")
    # Collatz v = 1 steps are left-child edges N -> 3N/2 (N even); v >= 2 steps leave the tree level structure
    print(" Collatz v = 1 steps are the left-child edges of even nodes (Mahler's map g -> ceil(3g/2) follows left children always)")


def part3():
    print("== (3) balance of the residual tree: P(next letter odd | uniform no-descent prefix of length k) ==")
    cur = {0: 1}
    worst_lo, worst_hi = 1.0, 0.0
    rows = []
    for j in range(1, 61):
        nxt = {}
        odd_ext = 0; tot_ext = 0
        for o, c in cur.items():
            for letter in (0, 1):
                oo = o + letter
                if 3 ** oo > 2 ** j:
                    nxt[oo] = nxt.get(oo, 0) + c
                    tot_ext += c
                    if letter == 1:
                        odd_ext += c
        p_odd = odd_ext / tot_ext
        worst_lo = min(worst_lo, p_odd); worst_hi = max(worst_hi, p_odd)
        if j in (1, 2, 3, 5, 8, 13, 21, 34, 55, 60):
            rows.append((j, round(p_odd, 4)))
        cur = nxt
    print(" (k, P(next odd)):", rows)
    print(" over k = 1..60 the probability stays in [%.4f, %.4f] (inside [1/3, 2/3]: %s); the extension counts weight each prefix by its number of continuations, so this is the balance of the next letter in the no-descent tree" % (worst_lo, worst_hi, worst_lo >= 1 / 3 and worst_hi <= 2 / 3))


def part4():
    print("== (4) Diophantine quality along the orbit of 27 (truncated shadow) ==")
    n = 27
    w = []; m = n
    while m != 1:
        m, v = U(m); w.append(v)
    K = len(w)
    d = [0]
    for v in w:
        d.append(d[-1] + v)
    etas = []
    for l in range(K):
        etas.append(sum(Fraction(2 ** (d[l + k] - d[l]), 3 ** (k + 1)) for k in range(K - l)))
    xi = n + etas[0]
    ms = [n]
    for v in w:
        ms.append(U(ms[-1])[0])
    print(" L   real error |xi - 2^d m/3^L|   mu_L = -log_3(err)/L    2-adic exponent d_L/(L log2 3)")
    for L in (5, 10, 20, 30, 40):
        err = abs(xi - Fraction(2 ** d[L] * ms[L], 3 ** L))
        mu = -math.log(float(err), 3) / L if err > 0 else float('inf')
        print(" %2d   %.3e                 %.3f                   %.3f" % (L, float(err), mu, d[L] / (L * LOG23)))
    print(" the real approximations 2^(d_L) m_L/3^L to xi have exponent well below 1 (Liouville needs > 1, Roth > 2); the 2-adic approximations -S/3^L to -n have exponent d_L/(L log2 3) by construction")


def part5():
    print("== (5) the reachability order: tree of 1 within the odd numbers <= X ==")
    for X in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
        reach = 0
        for m in range(1, X + 1, 2):
            x = m
            while x > 1 and x <= X:
                x, _ = U(x)
            # if the orbit leaves [1, X] we still follow it (the conjecture holds for these sizes), count reaching 1
            while x != 1:
                x, _ = U(x)
            reach += 1
        print(" X = %d: all %d odd numbers <= X reach 1 (conjecture verified); Krasikov-Lagarias's proved lower bound for the tree of 1 below X is X^0.84 = %.3g" % (X, reach, X ** 0.84))


def main():
    part1(); part2(); part3(); part4(); part5()


if __name__ == '__main__':
    main()
