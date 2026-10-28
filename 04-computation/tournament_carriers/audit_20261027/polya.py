"""Referee: Polya/Burnside count of oriented graphs (no 2-cycles) with n vertices and m arcs.
Independent check of census totals (j=1: m=n-1; j=2: m=2n-3; j=3: m=3n-6)."""
import itertools
import sys
from fractions import Fraction


def count_oriented(n, m):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tot = Fraction(0)
    nf = 0
    for p in itertools.permutations(range(n)):
        nf += 1
        seen = set()
        poly = [1]  # polynomial in x
        for (a, b) in pairs:
            if (a, b) in seen:
                continue
            # orbit of the ordered pair (a,b) under <p>
            orb = []
            x, y = a, b
            reverses = False
            while True:
                key = (min(x, y), max(x, y))
                if key in seen:
                    break
                seen.add(key)
                orb.append(key)
                x, y = p[x], p[y]
            # check whether some power maps (a,b) to (b,a)
            x, y = a, b
            for _ in range(len(orb) * 2):
                x, y = p[x], p[y]
                if (x, y) == (b, a):
                    reverses = True
                    break
                if (x, y) == (a, b):
                    break
            L = len(orb)
            factor = [1] + [0] * (L - 1) + ([0] if reverses else [2])
            if reverses:
                factor = [1]
            new = [0] * (len(poly) + len(factor) - 1)
            for i, c in enumerate(poly):
                for j, d in enumerate(factor):
                    new[i + j] += c * d
            poly = new
        tot += poly[m] if m < len(poly) else 0
    return tot / nf


if __name__ == "__main__":
    for n in range(3, 8):
        print("j=1 n=%d m=%d: %s" % (n, n - 1, count_oriented(n, n - 1)))
    for n in range(4, 7):
        print("j=2 n=%d m=%d: %s" % (n, 2 * n - 3, count_oriented(n, 2 * n - 3)))
    for n in range(5, 7):
        print("j=3 n=%d m=%d: %s" % (n, 3 * n - 6, count_oriented(n, 3 * n - 6)))
