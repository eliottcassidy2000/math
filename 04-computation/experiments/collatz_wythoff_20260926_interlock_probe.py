#!/usr/bin/env python3
"""collatz_wythoff_20260926_interlock_probe.py -- the owner's three-colouring decoded (Wythoff classes) and
tested against the Collatz step (session collatz-crossings-20260926, opus, 2026-09-26).

The owner's diagonal stripes (red, black, blue, red, black, red, black, blue, ...) are the classes
   red   = AA = { floor(floor(k phi) phi) }  (Zeckendorf representation ends in F_2 = 1),
   black = B  = { floor(k phi^2) }           (lowest Zeckendorf index odd, >= 3),
   blue  = AB = { floor(floor(k phi^2) phi) } (lowest Zeckendorf index even, >= 4),
which partition the positive integers (N = A + B, A = AA + AB). The owner's list matches this for n <= 34 and
puts 35 in blue where AB has 37 (35 is in AA). Targets as in the earlier probe; here also the colour of the
odd part U(n) = (3n+1)/2^v against the colour of n, and the colour pair (n, 2n) (halving preserves?).
Usage: python3 collatz_wythoff_20260926_interlock_probe.py [N=300000]
"""
import math, sys
from collections import Counter


def fib_upto(N):
    F = [1, 2]
    while F[-1] <= N:
        F.append(F[-1] + F[-2])
    return F


def lowest_index(n, F):
    i = len(F) - 1; low = None
    while n > 0:
        while F[i] > n:
            i -= 1
        low = i + 2; n -= F[i]; i -= 2
    return low


def colour(n, F):
    k = lowest_index(n, F)
    return 0 if k == 2 else (2 if k % 2 == 0 else 1)      # 0 red, 1 black, 2 blue


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def mi(pairs):
    N = len(pairs); cx, cy, cxy = Counter(), Counter(), Counter()
    for x, y in pairs:
        cx[x] += 1; cy[y] += 1; cxy[(x, y)] += 1
    I = sum((c / N) * math.log2(c * N / (cx[x] * cy[y])) for (x, y), c in cxy.items())
    Hy = -sum((c / N) * math.log2(c / N) for c in cy.values())
    return I, Hy


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    F = fib_upto(4 * N)
    phi = (1 + 5 ** 0.5) / 2
    # sanity: the three classes are AA, B, AB
    A = set(math.floor(k * phi) for k in range(1, N)); B = set(math.floor(k * phi * phi) for k in range(1, N))
    AA = set(math.floor(a * phi) for a in A); AB = set(math.floor(b * phi) for b in B)
    ok = all((colour(n, F) == 0) == (n in AA) and (colour(n, F) == 1) == (n in B) and (colour(n, F) == 2) == (n in AB) for n in range(1, 20000))
    print("colour classes equal (AA, B, AB) for n < 20000: %s; densities: red %.4f black %.4f blue %.4f (phi^-2 = %.4f, phi^-1 = %.4f... expected red 1/phi^3=%.4f? black 1/phi^2=%.4f blue 1/phi^3=%.4f)" % (
        ok, sum(1 for n in range(1, N) if colour(n, F) == 0) / N, sum(1 for n in range(1, N) if colour(n, F) == 1) / N, sum(1 for n in range(1, N) if colour(n, F) == 2) / N,
        phi ** -2, phi ** -1, phi ** -3, phi ** -2, phi ** -3))
    stop = {1: 0}
    def stopping(n):
        path = []; m = n
        while m not in stop:
            path.append(m); m = T(m)
        s = stop[m]
        for x in reversed(path):
            s += 1; stop[x] = s
        return stop[n]
    rows = []
    for n in range(2, N + 1):
        c = colour(n, F)
        m = n; desc = 0
        for _ in range(20):
            m = T(m)
            if m < n: desc = 1; break
        st = stopping(n) % 3
        if n % 2 == 1:
            x = 3 * n + 1; v = v2(x); U = x >> v
            rows.append((n, c, n % 2, v % 3, desc, st, colour(U, F) if U <= 4 * N else -1, colour(2 * n, F) if 2 * n <= 4 * N else -1))
        else:
            rows.append((n, c, n % 2, -1, desc, st, -1, colour(2 * n, F) if 2 * n <= 4 * N else -1))
    print("N = %d: mutual information I(colour(n); target) [target entropy]" % N)
    tg = {"parity": 2, "v2(3n+1) mod 3 (odd)": 3, "descent<=20": 4, "stop mod 3": 5, "colour(U(n)) (odd)": 6, "colour(2n)": 7}
    line = ""
    for name, ti in tg.items():
        pairs = [(r[1], r[ti]) for r in rows if r[ti] >= 0]
        I, Hy = mi(pairs)
        line += "   %s: %.5f [%.3f]" % (name, I, Hy)
    print(line)
    M = [[0] * 3 for _ in range(3)]
    for r in rows:
        if r[7] >= 0: M[r[1]][r[7]] += 1
    print("colour(n) -> colour(2n) (doubling), row-normalised:")
    for i in range(3):
        s = sum(M[i]) or 1
        print("   %s: " % ["red", "black", "blue"][i] + "  ".join("%.4f" % (M[i][j] / s) for j in range(3)))
    M = [[0] * 3 for _ in range(3)]
    for r in rows:
        if r[6] >= 0: M[r[1]][r[6]] += 1
    print("colour(n) -> colour(U(n)) for odd n (Syracuse step), row-normalised:")
    for i in range(3):
        s = sum(M[i]) or 1
        print("   %s: " % ["red", "black", "blue"][i] + "  ".join("%.4f" % (M[i][j] / s) for j in range(3)))


if __name__ == '__main__':
    main()
