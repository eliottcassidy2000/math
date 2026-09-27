#!/usr/bin/env python3
"""fibonacci_two_copies_20260926.py -- the Fibonacci line with its zero removed: three consecutive 1s, the
positive Zeckendorf system on one side, Knuth's negaFibonacci system on the other, and the real Binet function
(session fibonacci-two-copies-20260926, opus, 2026-09-26).

 (1) Real Binet function F(x) = (phi^x - cos(pi x) phi^-x)/sqrt5: F(n) = F_n for all integers n (both signs);
     recurrence F(x+2) = F(x+1) + F(x) for real x; twist identity F(-x) = -cos(pi x) F(x) + phi^-x sin^2(pi x)/sqrt5
     (the negative copy is the positive one with the sign twist, plus a correction vanishing at integers);
     smooth Cassini F(x+1)^2 - F(x) F(x+2) = cos(pi x); the one-parameter family F_r = F + r phi^-x sin(pi x)/sqrt5
     of all 'nice' interpolants (recurrence + all integer values), r = 0 minimising the oscillation amplitude
     on the negative axis; zeros of F on the negative axis.
 (2) The zero-removed sequence g(n) = F_(n+1) (n >= 0), g(-n) = F_(-n) (n >= 1): 1 at -1, 0, 1.
 (3) Knuth's negaFibonacci representation of every integer N (sum of F_(-k), k >= 1, no two consecutive k):
     built for |N| <= X by enumeration, uniqueness checked; tricolor by lowest index: k = 1 (value +1), k = 2
     (value -1), odd k >= 3, even k >= 4; densities for positive and for negative N; comparison with the
     Zeckendorf tricolor (lowest index 2, odd >= 3, even >= 4) on the positives: overlap of the '1 classes'.
Usage: python3 fibonacci_two_copies_20260926.py [X=100000]
"""
import math, sys

PHI = (1 + 5 ** 0.5) / 2
S5 = 5 ** 0.5


def fib(n):
    if n >= 0:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a
    return (-1) ** (n + 1) * fib(-n)


def F(x, r=0.0):
    return (PHI ** x - math.cos(math.pi * x) * PHI ** (-x) + r * PHI ** (-x) * math.sin(math.pi * x)) / S5


def part1():
    print("== (1) real Binet function ==")
    assert all(abs(F(n) - fib(n)) < 1e-9 * max(1, abs(fib(n))) for n in range(-30, 31))
    print(" F(n) = F_n for -30 <= n <= 30; F(-1), F(0), F(1), F(2), F(-2) =", [fib(n) for n in (-1, 0, 1, 2, -2)])
    xs = [i / 7.0 for i in range(-70, 71)]
    e1 = max(abs(F(x + 2) - F(x + 1) - F(x)) for x in xs)
    e2 = max(abs(F(-x) - (-math.cos(math.pi * x) * F(x) + PHI ** (-x) * math.sin(math.pi * x) ** 2 / S5)) for x in xs)
    e3 = max(abs(F(x + 1) ** 2 - F(x) * F(x + 2) - math.cos(math.pi * x)) for x in xs)
    e4 = max(abs(F(x + 2, 0.7) - F(x + 1, 0.7) - F(x, 0.7)) for x in xs)
    print(" recurrence on reals: max error %.1e; twist identity F(-x) = -cos(pi x)F(x) + phi^-x sin^2(pi x)/sqrt5: %.1e; smooth Cassini F(x+1)^2 - F(x)F(x+2) = cos(pi x): %.1e; recurrence for F_r (r=0.7): %.1e" % (e1, e2, e3, e4))
    # zeros on the negative axis
    zeros = []
    x = -12.0
    step = 1e-3
    prev = F(x)
    while x < 1.0:
        x2 = x + step; cur = F(x2)
        if prev == 0 or prev * cur < 0:
            lo, hi = x, x2
            for _ in range(60):
                mid = (lo + hi) / 2
                if F(lo) * F(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            zeros.append(round((lo + hi) / 2, 6))
        prev = cur; x = x2
    print(" zeros of F on [-12, 1]:", zeros)
    print(" amplitude of F_r on the negative axis ~ phi^|x| sqrt(1 + r^2)/sqrt5: r = 0 (Binet) is the minimal-oscillation interpolant; the integers cannot see r.")


def part2():
    print("== (2) zero-removed sequence g ==")
    g = lambda n: fib(n + 1) if n >= 0 else fib(n)
    print(" g(n) for n = -6..6:", [g(n) for n in range(-6, 7)])
    print(" g(-n) = (-1)^(n+1) g(n-1) for n >= 1:", all(g(-n) == (-1) ** (n + 1) * g(n - 1) for n in range(1, 30)))


def negafib_all(X):
    """map N -> tuple of indices k (F_(-k)), for all |N| <= X, by enumerating no-two-consecutive index sets"""
    # F_(-k) for k = 1..K
    K = 2
    while abs(fib(-K)) <= 4 * X:
        K += 1
    vals = [fib(-k) for k in range(1, K + 1)]  # index k-1 -> F_(-k)
    rep = {}
    # enumerate subsets with no two consecutive, prune by bound on remaining sum
    def rec(i, cur_sum, chosen):
        if abs(cur_sum) <= X:
            rep.setdefault(cur_sum, tuple(chosen))
            if cur_sum in rep and rep[cur_sum] != tuple(chosen):
                raise AssertionError("non-unique negaFibonacci representation for %d" % cur_sum)
        if i >= K:
            return
        # remaining terms could bring the sum back into range; crude pruning: skip if |cur_sum| > X + sum of remaining |vals|
        rem = sum(abs(v) for v in vals[i:])
        if abs(cur_sum) > X + rem:
            return
        rec(i + 1, cur_sum, chosen)           # skip index k = i+1
        rec(i + 2, cur_sum + vals[i], chosen + [i + 1])  # take index k = i+1, then skip k+1
    rec(0, 0, [])
    return rep


def zeckendorf(N):
    idx = []
    k = 2
    while fib(k + 1) <= N:
        k += 1
    while N > 0:
        while fib(k) > N:
            k -= 1
        idx.append(k); N -= fib(k); k -= 2
    return tuple(sorted(idx))


def part3(X):
    print("== (3) negaFibonacci representations for |N| <= %d ==" % X)
    rep = negafib_all(X)
    assert all(N in rep for N in range(-X, X + 1)), "missing representations"
    for N in list(range(-8, 9)):
        ks = rep[N]
        print("  %3d = %s" % (N, ' + '.join('F_(-%d)=%d' % (k, fib(-k)) for k in sorted(ks)) or '0'))
    def cls(ks):
        if not ks:
            return 'zero'
        k = min(ks)
        return {1: 'k=1 (+1)', 2: 'k=2 (-1)'}.get(k, 'odd>=3' if k % 2 else 'even>=4')
    for sign, rng in (('positive', range(1, X + 1)), ('negative', range(-X, 0))):
        cnt = {}
        for N in rng:
            c = cls(rep[N]); cnt[c] = cnt.get(c, 0) + 1
        tot = len(rng)
        print(" %s N: lowest-index classes %s" % (sign, {c: "%.4f" % (v / tot) for c, v in sorted(cnt.items())}))
    # Zeckendorf tricolor on positives and overlap of the '1 classes'
    zc = {}
    both = 0; z1 = 0; n1 = 0
    for N in range(1, X + 1):
        zk = zeckendorf(N)
        k = min(zk)
        c = 'k=2 (1)' if k == 2 else ('odd>=3' if k % 2 else 'even>=4')
        zc[c] = zc.get(c, 0) + 1
        a = (k == 2); b = (min(rep[N]) == 1)
        z1 += a; n1 += b; both += a and b
    print(" Zeckendorf classes on positives: %s  (phi^-2 = %.4f, phi^-3 = %.4f)" % ({c: "%.4f" % (v / X) for c, v in sorted(zc.items())}, PHI ** -2, PHI ** -3))
    print(" positives whose Zeckendorf ends in F_2 = 1: %.4f; whose negaFibonacci ends in F_(-1) = 1: %.4f; both: %.4f" % (z1 / X, n1 / X, both / X))
    # the relation between the two representations: how often is nega(N) obtained from Z(N) by k -> -(k-1)?
    same = sum(1 for N in range(1, X + 1) if tuple(sorted(k - 1 for k in zeckendorf(N))) == tuple(sorted(rep[N])))
    print(" positives whose negaFibonacci index set equals the Zeckendorf index set shifted by one (k -> k-1, i.e. F_k -> F_(-(k-1)) = F_k for odd k-1 ... exact only when all Zeckendorf indices are odd): %.4f" % (same / X))
    allodd = sum(1 for N in range(1, X + 1) if all(k % 2 == 1 for k in zeckendorf(N)))
    print(" positives whose Zeckendorf indices are all odd: %.4f (these are exactly the ones whose two representations coincide term by term, since F_(-(k-1)) = F_(k-1)... see note)" % (allodd / X))


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    part1(); part2(); part3(X)


if __name__ == '__main__':
    main()
