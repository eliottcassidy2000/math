#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  Add the raw lines of x2_box.c and print the exact values
B(n, k) = E (-1)^(N_k) (uniform labelled tournament on N = n + k vertices, k disjoint marked pairs), the deletion-pair
side of (4.3) (k = 1: B(m+1, 1) = E_R[b_R^2] at order m), D_N and D_(N-1), and the chain B(n,k)^2 <= B(n,k+1).
Usage: python3 x2_merge.py file [file ...]"""
import sys, re
from fractions import Fraction
from math import factorial

lab = {}; cls = {}
for path in sys.argv[1:]:
    for line in open(path):
        g = re.match(r'N (\d+) labelled tournaments (\d+) pointwise mismatches route1/route2 (\d+) ;(.*)', line)
        if g:
            N = int(g.group(1)); tot = int(g.group(2))
            assert tot == 2 ** (N * (N - 1) // 2)
            lab[N] = (int(g.group(3)), {int(k): (Fraction(int(a), tot), Fraction(int(b), tot))
                                        for k, a, b in re.findall(r' k (\d+) n \d+ sum_route1 (-?\d+) sum_route2 (-?\d+)', g.group(4))})
            continue
        g = re.match(r'N (\d+) classes (\d+) W (\d+) sumD (-?\d+) sumDel (-?\d+) route1_compared (\d+) route1_mismatch (\d+) ;(.*)', line)
        if g:
            N = int(g.group(1)); t = cls.setdefault(N, dict(classes=0, W=0, sD=0, sDel=0, cmp=0, mis=0, sums={}, places={}))
            t['classes'] += int(g.group(2)); t['W'] += int(g.group(3)); t['sD'] += int(g.group(4)); t['sDel'] += int(g.group(5))
            t['cmp'] += int(g.group(6)); t['mis'] += int(g.group(7))
            for k, pl, s in re.findall(r' k (\d+) n \d+ placements (\d+) sum (-?\d+)', g.group(8)):
                t['sums'][int(k)] = t['sums'].get(int(k), 0) + int(s); t['places'][int(k)] = int(pl)

B = {}
print("labelled brute force (marked pairs {0,1},{2,3},...; route 1 = explicit odd-cycle count, route 2 = product of mod-4 signs over transversal deletions):")
for N in sorted(lab):
    mis, d = lab[N]
    print("  N = %d: pointwise mismatches between the routes: %d ; " % (N, mis) + " ; ".join("B(%d,%d) = %s%s" % (N - k, k, a, "" if a == b else " (ROUTE 2: %s)" % b) for k, (a, b) in sorted(d.items())))
    for k, (a, b) in d.items(): B[(N - k, k)] = a
print("class representatives with weights N!/|Aut| (all placements of the marked pairs, route 2):")
for N in sorted(cls):
    t = cls[N]; W = t['W']; full = 2 ** (N * (N - 1) // 2)
    line = "  N = %d: %d classes, weight %s 2^C(N,2) ; route 1 compared on %d placements, mismatches %d ; D_%d = %s ; D_%d (from deletions) = %s" % (
        N, t['classes'], "=" if W == full else "NOT", t['cmp'], t['mis'], N, Fraction(t['sD'], W), N - 1, Fraction(t['sDel'], N * W))
    print(line)
    for k in sorted(t['sums']):
        nm = factorial(N) // (2 ** k * factorial(k) * factorial(N - 2 * k))
        assert nm == t['places'][k], (N, k, nm, t['places'][k])
        v = Fraction(t['sums'][k], nm * W)
        if (N - k, k) in B and B[(N - k, k)] != v: print("    DISAGREES with the labelled value: B(%d,%d)" % (N - k, k))
        B[(N - k, k)] = v
        print("    B(%d,%d) = %s = %.9f" % (N - k, k, v, float(v)))
print("chain B(n,k)^2 <= B(n,k+1) (0 <= k <= n-1) and |D_n|^(2^k) <= B(n,k):")
ok = True
for n in sorted({n for (n, k) in B}):
    ks = sorted(k for (nn, k) in B if nn == n)
    for k in ks:
        if (n, k + 1) in B:
            c = B[(n, k)] ** 2 <= B[(n, k + 1)]; ok &= c
            print("    n = %d, k = %d: B(n,k)^2 = %.3e <= B(n,k+1) = %.3e : %s%s" % (n, k, float(B[(n, k)] ** 2), float(B[(n, k + 1)]), c,
                  "  (equality)" if B[(n, k)] ** 2 == B[(n, k + 1)] else ""))
        if k >= 1 and (n, 0) in B:
            c = abs(B[(n, 0)]) ** (2 ** k) <= B[(n, k)]; ok &= c
print("ALL CHAIN INEQUALITIES HOLD" if ok else "A CHAIN INEQUALITY FAILS")
print("deletion-pair values B(m+1, 1) = E_R[b_R^2] at order m:", {n - 1: str(B[(n, 1)]) for (n, k) in sorted(B) if k == 1})
