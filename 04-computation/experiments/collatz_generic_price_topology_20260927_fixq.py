#!/usr/bin/env python3
"""collatz_generic_price_topology_20260927_fixq.py -- the fixed points of the conjugacy map Q, level by level
(session collatz-posets-zeta5-20260927, opus, 2026-09-27; companion of collatz_generic_price_topology_20260927.py part 3).

Q(x) = sum_i (parity of T^i x) 2^i. Solutions of Q(x) = x mod 2^k are lifted level by level (each solution mod 2^k has two
candidate lifts mod 2^(k+1)); N_k counts them. A solution mod 2^k is persistent to depth D if some lift chain reaches
depth D. The inclusion {0} u {-2^j} subset Fix(Q) is a theorem (Q(2y) = 2Q(y); Q(-1) = -1); the computation shows
which other solutions persist. The odd fixed points satisfy the self-referential tree of equations Q(a z + b) = z
starting at (a, b) = (3, 2): a node (a, b) with b odd is dead; with b even its children are (a, b/2) and
(3a, (3a + 3b + 1)/2); the -1 chain is (3^i, 3^i - 1).
Usage: python3 collatz_generic_price_topology_20260927_fixq.py [DEPTH=40]
"""
import sys


def Qk(x, k):
    q = 0
    for i in range(k):
        b = x & 1
        q |= b << i
        x = (3 * x + 1) // 2 if b else x // 2
    return q


def main():
    D = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    sols = {1: [x for x in range(2) if Qk(x, 1) % 2 == x]}
    for k in range(2, D + 1):
        mod = 1 << k
        nxt = []
        for x in sols[k - 1]:
            for y in (x, x + (1 << (k - 1))):
                if Qk(y, k) % mod == y:
                    nxt.append(y)
        sols[k] = nxt
    print(" N_k for k = 1..%d:" % D, [len(sols[k]) for k in range(1, D + 1)])
    # persistence: solutions mod 2^k that have a lift chain to depth D
    persistent = {}
    live = set(sols[D])
    for k in range(D, 0, -1):
        persistent[k] = sorted({x % (1 << k) for x in live})
        live = {x % (1 << k) for x in live}
    for k in (10, 16, 22, 28, 34, D):
        if k > D:
            continue
        mod = 1 << k
        names = []
        for x in persistent[k]:
            if x == 0:
                names.append("0")
            elif (mod - x) & (mod - x - 1) == 0:
                names.append("-2^%d" % ((mod - x).bit_length() - 1))
            else:
                names.append(str(x))
        print(" solutions mod 2^%d persistent to depth %d: %d of %d: %s" % (k, D, len(persistent[k]), len(sols[k]), names))
    # the odd-fixed-point tree of equations, live nodes per level
    level = [(3, 2)]; counts = []
    for depth in range(1, 31):
        nxt = []
        for a, b in level:
            if b % 2 == 0:
                nxt.append((a, b // 2))
                nxt.append((3 * a, (3 * a + 3 * b + 1) // 2))
        level = [(a, b) for a, b in nxt if b % 2 == 0]
        counts.append(len(level))
    print(" live nodes of the odd-fixed-point equation tree Q(az + b) = z at depths 1..30:", counts)
    print(" the -1 chain (3^i, 3^i - 1) is live at every depth: %s" % all(((3 ** i) % 2 == 1 and (3 ** i - 1) % 2 == 0) for i in range(1, 40)))


if __name__ == '__main__':
    main()
