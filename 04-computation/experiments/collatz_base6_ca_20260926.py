#!/usr/bin/env python3
"""collatz_base6_ca_20260926.py -- the Collatz map is a radius-1 cellular automaton on base-6 digit strings
(Cloney-Goles-Vichniac 1987), verified by brute force (session collatz-oscillation-20260926, opus, 2026-09-26).

For x with base-6 digits x_0 (least significant), x_1, ..., the digit i of T(x) = x/2 (x even) or (3x+1)/2 (x odd)
is claimed to be a function of (x_(i-1), x_i, x_(i+1), parity of x_0)... more precisely of (x_(i-1), x_i, x_(i+1))
together with the global parity bit x_0 mod 2 (which decides which branch runs). We test: for random x, collect
the map (x_0 mod 2, x_(i-1), x_i, x_(i+1)) -> T(x)_i and report any inconsistency; then the same with radius 2
to confirm that radius 1 already suffices, and print the local rule tables for the odd branch. The reason:
multiplying a base-6 digit by 3 gives 3 x_i in {0,3,6,9,12,15}, and an incoming carry c in {0,1,2} never crosses a
multiple of 6, so the outgoing carry floor((3 x_i + c)/6) depends on x_i only; halving reads one digit above.
Usage: python3 collatz_base6_ca_20260926.py [N=200000]
"""
import random, sys


def digits6(x, width):
    d = []
    for _ in range(width):
        d.append(x % 6); x //= 6
    return d


def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    random.seed(1)
    W = 14
    for radius in (0, 1, 2):
        table = {}; bad = 0
        for _ in range(N):
            x = random.randrange(1, 6 ** (W - 3))
            dx = digits6(x, W); dt = digits6(T(x), W); par = x % 2
            for i in range(W - 3):
                key = (par, tuple(dx[max(0, i - radius):i + radius + 1]), i - max(0, i - radius))
                if key in table and table[key] != dt[i]:
                    bad += 1
                table[key] = dt[i]
        print("radius %d: %d inconsistencies over %d samples x %d digits (0 means digit i of T(x) is a function of the neighbourhood and the parity bit)" % (radius, bad, N, W - 3))
    # print the odd-branch rule for interior digits: digit i of (3x+1)/2 from (x_(i-1), x_i, x_(i+1)), i >= 1
    table = {}
    for _ in range(N):
        x = random.randrange(1, 6 ** (W - 3)) | 1
        dx = digits6(x, W); dt = digits6(T(x), W)
        for i in range(1, W - 3):
            table[(dx[i - 1], dx[i], dx[i + 1])] = dt[i]
    print("odd branch, interior digit i of (3x+1)/2 as a function of (x_(i-1), x_i, x_(i+1)): %d of 216 neighbourhoods seen" % len(table))
    # show dependence structure: does it depend on x_(i-1)? on x_(i+1)?
    dep_left = any(table.get((a, b, c)) != table.get((a2, b, c)) for a in range(6) for a2 in range(6) for b in range(6) for c in range(6) if (a, b, c) in table and (a2, b, c) in table)
    dep_right = any(table.get((a, b, c)) != table.get((a, b, c2)) for a in range(6) for b in range(6) for c in range(6) for c2 in range(6) if (a, b, c) in table and (a, b, c2) in table)
    print("   depends on the lower digit x_(i-1): %s; on the upper digit x_(i+1): %s" % (dep_left, dep_right))
    print("   rule for x_i = 0..5 with x_(i-1) = 0 and x_(i+1) = 0..5 (rows x_i, columns x_(i+1)):")
    for b in range(6):
        print("     x_i=%d: " % b + " ".join(str(table.get((0, b, c), '?')) for c in range(6)))


if __name__ == '__main__':
    main()
