#!/usr/bin/env python3
"""collatz_procgen_20260922_loops_escapes.py -- escapes from the hostile 3-adic points 1, 1/2, 43/32, 59/64.

Part A (exact, proof check).  For a positive dyadic rational q (a 3-adic unit) the reverse moves are
x -> (2^k x - 1)/3, legal iff the result is a 3-adic unit (and nonzero).  A path of s moves with K halvings
ending at x_s has multiplier 2^K/3^s = (x_s + B/3^s)/q with B/3^s >= 1/3 (s >= 1) and >= 1/3 + 1/9 (s >= 2).
The script enumerates the finite NON-INTEGER part of the reverse tree of q (moves with 2^k q < 2^b stay
non-integral; large k give integers), records every node with its exact (s,K), and checks:
  (1) every non-integer node (s >= 1) has multiplier > 1;
  (2) every first integer node x reached at depth 1 has (x + 1/3)/q > 1, at depth >= 2 the bound
      (1 + 4/9)/q > 1 holds (q < 13/9); from an integer node all later nodes are integers >= 1, so every
      later prefix has multiplier >= 13/(9q) > 1.
  Hence q is in Bad_inf (no certificate at any level).  Also prints the cheapest routes from q to 1 and 1/2.
Part B (exact integer checks of the transfer lemmas on random members of the threads):
  1-escape  m = 1 + 3^k u  --(loop of length k-1, then move 0)-->  2^K u < m;
  1/2-transfer m = (1 + 3^k w)/2 --(k_1+1, k_2..k_s, 0)--> 2^K w, with m < 2^K w < 2m;
  43/32 route m = (43 + 3^k W)/32 --(3,2,0) then 1-escape--> 2^{K'} W;  59/64 route to the 1/2 thread.
Part C (tables): transfer prices 2c(k-1)/3, the two-move exit conditions, and the route inequalities.
Loops for Part B are read from a family file ('s=.. k:..' lines, e.g. family1000r.txt).
Part D (with the level-12 dump of loops_seesaw): after the 1/2-transfer the landing point y = 2^K w is an
arbitrary unit class; m descends through one level-12 certificate of y iff its factor f < 3/(2c(k-1)).
The table counts, for each k, the unit classes mod 3^13 (out of 2*3^12) where this fails.
Usage: python3 ..._loops_escapes.py FAMILYFILE [SEESAW_DUMP]
"""
import sys, re, random
from fractions import Fraction as F

def unit3(x):          # x a Fraction: is it a nonzero 3-adic unit (v_3(x) = 0)?
    return x != 0 and x.numerator % 3 != 0 and x.denominator % 3 != 0

def rev(x, k):
    return (x * 2 ** k - 1) / 3

def nonint_tree(q, kmax=40):
    """BFS over reverse moves from q, expanding only non-integer nodes. Returns list of
    (x, s, K, path) for all nodes (non-integer nodes and the first integer nodes)."""
    out = []; frontier = [(q, 0, 0, [])]
    while frontier:
        nf = []
        for x, s, K, path in frontier:
            for k in range(0, kmax + 1):
                y = rev(x, k)
                if not unit3(y) or y < 0:
                    continue
                node = (y, s + 1, K + k, path + [k])
                out.append(node)
                if y.denominator != 1:
                    nf.append(node)
        frontier = nf
    return out

def partA(q, name):
    """Proof check that q is in Bad_inf.  Every node (x,s) has multiplier (x + B_s/3^s)/q with
    B_s/3^s >= beta_s = (1 - 3^-s)/2.  Non-integer nodes: exact check.  Integer nodes x >= 2 (s >= 1):
    multiplier >= (2 + 1/3)/q > 1 as q < 7/3.  Integer node x = 1 at depth s: multiplier >= (1 + beta_s)/q;
    let s* = least s with 1 + beta_s > q; the value 1 at depth s needs a node x = 4/2^k (x in {4,2,1,1/2,..})
    at depth s-1, so it suffices that no node of depth <= s*-2 has that form (checked on the explicit nodes of
    depth <= 1, which are all listed; we require s* <= 3)."""
    nodes = nonint_tree(q)
    nonint = [n for n in nodes if n[0].denominator != 1]
    ints = [n for n in nodes if n[0].denominator == 1]
    ok = all(F(2 ** K, 3 ** s) > 1 for x, s, K, path in nonint)
    c7 = q < F(7, 3)
    sstar = 1
    while not (1 + (1 - F(1, 3 ** sstar)) / 2 > q):
        sstar += 1
        if sstar > 60: break
    def pow2le4(x):
        return x <= 4 and x.numerator in (1, 2, 4) and (x.denominator & (x.denominator - 1)) == 0
    small = [n for n in nodes if pow2le4(n[0])]
    mind = min([n[1] for n in small], default=99)
    # nodes in {1,2,4} can only be integers; integer nodes created at depth d have value >= min over them;
    # a node in {1,2,4} at depth <= sstar-2 would have to be one of the listed nodes (depth <= 1 fully listed
    # for integer children of q and of non-integer nodes) -- we require sstar <= 3 so only depth <= 1 matters.
    c_small = sstar <= 3 and mind > sstar - 2
    print(f"{name}: non-integer reverse-tree nodes: {sorted(set(str(n[0]) for n in nonint))}")
    for x, s, K, path in sorted(nonint, key=lambda n: (n[1], n[2])):
        print(f"   node {str(x):>6} at depth {s}, moves {path}, multiplier 2^{K}/3^{s} = {float(F(2**K,3**s)):.4f}")
    d1 = [n for n in ints if n[1] == 1]
    print(f"   non-integer nodes all > 1: {ok}; q < 7/3: {c7}; s* = {sstar} (value 1 must not occur before depth s*),"
          f" first node of the form 4/2^k at depth {mind if mind < 99 else '>1'}; smallest depth-1 integer node {min(n[0] for n in d1)}")
    res = ok and c7 and c_small
    print(f"   => {name} in Bad_inf: {res}")
    return res

def K0(s):
    return (3 ** (s + 1)).bit_length() - 1

def c_of(s):
    return F(2 ** K0(s), 3 ** s)

def run_rev(m, ks):
    x = m
    for k in ks:
        t = x * 2 ** k - 1
        assert t % 3 == 0
        x = t // 3
        assert x >= 1 and x % 3 != 0, "illegal"
    return x

def partB(loops):
    random.seed(20260922)
    n1 = n2 = n3 = n4 = 0
    for k in range(2, 60):
        s = k - 1
        if s not in loops and s != 1:
            continue
        ks = [2] if s == 1 else loops[s]
        K = sum(ks)
        for _ in range(40):
            u = random.randrange(1, 10 ** 30)
            if u % 3 == 0: continue
            m = 1 + 3 ** k * u
            y = run_rev(m, ks + [0]); assert y == 2 ** K * u and y < m; n1 += 1
            w = u if u % 2 == 1 else u + 3        # w odd, 3 !| w  -> m integer
            if w % 3 == 0: continue
            m2 = (1 + 3 ** k * w) // 2
            y2 = run_rev(m2, [ks[0] + 1] + ks[1:] + [0]); assert y2 == 2 ** K * w
            if k >= 3:
                assert m2 < y2 < 2 * m2
            n2 += 1
        # 43/32 route: 32 m = 43 + 3^k W, then moves 3,2,0 lead to 1 + 3^(k-3) W, then 1-escape (loop k-4)
        if k >= 6 and (k - 4) in loops:
            ks4 = loops[k - 4]; K4 = sum(ks4)
            for _ in range(20):
                W = random.randrange(1, 10 ** 30)
                if W % 3 == 0 or (43 + 3 ** k * W) % 32: continue
                m3 = (43 + 3 ** k * W) // 32
                y3 = run_rev(m3, [3, 2, 0]); assert y3 == 1 + 3 ** (k - 3) * W
                y4 = run_rev(m3, [3, 2, 0] + ks4 + [0]); assert y4 == 2 ** K4 * W
                assert (y4 < m3) == (2 ** (K4 + 5) < 3 ** k); n3 += 1
            for _ in range(20):
                W = random.randrange(1, 10 ** 30)
                if W % 3 == 0 or (59 + 3 ** k * W) % 64: continue
                m5 = (59 + 3 ** k * W) // 64
                y5 = run_rev(m5, [3, 2, 0])
                assert 2 * y5 == 1 + 3 ** (k - 3) * W        # lands on the 1/2 thread, precision k-3
                y6 = run_rev(m5, [3, 2, 0, ks4[0] + 1] + ks4[1:] + [0]); assert y6 == 2 ** K4 * W
                n4 += 1
    print(f"Part B: exact integer checks passed: 1-escape {n1}, 1/2-transfer {n2}, 43/32 route {n3}, 59/64 route {n4}")

def partC():
    print("Part C: prices along the threads (c = c(s) = 2^K0(s)/3^s, s = k-1 for 1 and 1/2, s = k-4 for 43/32, 59/64)")
    print("  k   c(k-1)   1-escape c/3   1/2 price 2c/3   2-move exit y=2,8 mod 9 (c<9/4)   43/32: 32c(k-4)/81  59/64 price 64c(k-4)/81")
    for k in range(3, 41):
        c = c_of(k - 1); c4 = c_of(k - 4) if k >= 6 else None
        e43 = f"{float(32*c4/81):.4f}{' <1' if 32*c4 < 81 else ' >1'}" if c4 else "   -"
        e59 = f"{float(64*c4/81):.4f}" if c4 else "  -"
        print(f"  {k:2d}  {float(c):.4f}   {float(c/3):.4f}        {float(2*c/3):.4f}          {'yes' if 4*c < 9 else 'no '}"
              f"                               {e43}            {e59}")
    # asymptotic frequencies (Weyl): c(s) = (3/2) 2^eps, eps uniform on (0,1)
    import math
    print("  frequencies over k (eps(k) equidistributed): 1/2 two-move exit on y=2,8 (mod 9) needs c<9/4: "
          f"{math.log2(1.5):.4f};  43/32 direct route needs c<81/32: {math.log2(27/16):.4f}")

def main():
    loops = {}
    for line in open(sys.argv[1]):
        m = re.search(r's=\s*(\d+).*?k:(\S+)', line)
        if m:
            loops[int(m.group(1))] = [int(t) for t in m.group(2).split(',')]
    allok = True
    for q, name in [(F(1, 2), '1/2'), (F(43, 32), '43/32'), (F(59, 64), '59/64'), (F(145, 128), '145/128'),
                    (F(209, 256), '209/256'), (F(371, 256), '371/256'), (F(499, 512), '499/512')]:
        try:
            allok &= partA(q, name)
        except ValueError:
            print(name, 'failed')
    partB(loops)
    partC()
    if len(sys.argv) > 2:
        rows = [l.split() for l in open(sys.argv[2])]
        print("Part D: 1/2-transfer followed by one level-12 certificate of the landing class y = 2^K w")
        print("  k   c(k-1)   threshold 3/(2c)   failing classes (f >= threshold or exceptional) of 2*3^12   fraction")
        for k in range(3, 41):
            c = c_of(k - 1); th = 3 / (2 * c)
            bad = [r for r in rows if r[4] == 'EXC' or F(2 ** int(r[2]), 3 ** int(r[1])) >= th]
            print(f"  {k:2d}  {float(c):.4f}   {float(th):.4f}             {len(bad):4d}                                   {len(bad)/(2*3**12):.2e}")

if __name__ == '__main__':
    main()
