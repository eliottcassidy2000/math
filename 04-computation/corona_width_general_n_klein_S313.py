#!/usr/bin/env python3
"""
corona_width_general_n_klein_S313.py — klein-2026-07-15-S313 (cont.3)

THE LANDAU CORONA AT GENERAL n (THM-868 §3 open item: "is the width floor(n/2)+1?").

Setup (THM-868 conventions): score deviations d_v = 2s_v - (n-1), all ≡ n-1 (mod 2),
sum 0, |d_v| <= n-1.  Shells: x = sum d^2 (step 8).  Per shell:
  A(x) = # sorted deviation multisets in the lattice box (parity + sum 0 + |d|<=n-1)
  B(x) = # of those that are LANDAU (s_(1)<=...<=s_(n), partial sums >= C(k,2), total C(n,2))
Landau filter: trivial zone (B=A), CORONA (0 < B < A), dead zone (B=0 < A, above ceiling).
Corona width W(n) = #{x : 0 < B(x) < A(x)}.  THM-868: W(8) = 5 with bites
10/12, 10/14, 11/15, 6/12, 1/10 at x = 136..168.  Question: W(n) = floor(n/2)+1?
"""
from math import comb
import sys

def corona(n):
    par = (n - 1) % 2
    vals = [v for v in range(-(n - 1), n) if abs(v) % 2 == par or (par == 0 and v % 2 == 0)]
    vals = [v for v in range(-(n - 1), n) if (v - (n - 1)) % 2 == 0]   # v ≡ n-1 mod 2
    vals.sort()
    A, B = {}, {}
    ms = []
    def rec(idx, left, tot, sq):
        if left == 0:
            if tot != 0: return
            x = sq
            A[x] = A.get(x, 0) + 1
            s = sorted((d + n - 1) // 2 for d in ms)
            ps = 0; ok = True
            for k, sv in enumerate(s, 1):
                ps += sv
                if ps < comb(k, 2): ok = False; break
            if ok and ps == comb(n, 2):
                B[x] = B.get(x, 0) + 1
            return
        if idx == len(vals): return
        v = vals[idx]
        # prune: remaining values range
        lo = tot + left * vals[idx]           # if all remaining = smallest available (vals sorted asc, choose from idx..)
        hi = tot + left * vals[-1]
        if lo > 0 or hi < 0: return
        for c in range(left + 1):
            if c: ms.append(v)
            rec(idx + 1, left - c, tot + c * v, sq + c * v * v)
        for c in range(left, 0, -1):
            if ms and ms[-1] == v: ms.pop()
    rec(0, n, 0, 0)
    shells = sorted(A)
    ceiling = (n ** 3 - n) // 3
    cor = [x for x in shells if 0 < B.get(x, 0) < A[x]]
    trivial_max = max([x for x in shells if B.get(x, 0) == A[x]], default=None)
    return A, B, cor, trivial_max, ceiling

print("n | shells | ceiling | corona levels (x: B/A) | width | floor(n/2)+1")
widths = {}
for n in range(4, 17):
    A, B, cor, triv, ceil_ = corona(n)
    bites = ", ".join(f"{x}:{B.get(x,0)}/{A[x]}" for x in cor)
    W = len(cor)
    widths[n] = W
    pred = n // 2 + 1
    mark = "MATCH" if W == pred else "DIFFER"
    print(f"n={n:2d} | {len(A):4d} | {ceil_:5d} | [{bites}] | W={W} vs {pred} {mark}")
    # cross-checks
    assert max(x for x in A if B.get(x, 0) > 0) == ceil_, "ceiling mismatch"
    if n == 8:
        expect = {136: (10, 12), 144: (10, 14), 152: (11, 15), 160: (6, 12), 168: (1, 10)}
        got = {x: (B.get(x, 0), A[x]) for x in cor}
        print("   THM-868 n=8 cross-check:", "PASS" if got == expect else f"FAIL {got}")
    # contiguity of the corona
    if cor:
        contig = all(cor[i + 1] - cor[i] == 8 for i in range(len(cor) - 1)) and cor[-1] == ceil_
        print(f"   corona contiguous ending at ceiling: {contig}; span x = {cor[0]}..{cor[-1]}")

print()
print("widths:", widths)
even_w = {n: w for n, w in widths.items() if n % 2 == 0}
odd_w = {n: w for n, w in widths.items() if n % 2 == 1}
print("even n:", even_w, " -> floor(n/2)+1?", all(w == n // 2 + 1 for n, w in even_w.items()))
print("odd  n:", odd_w, " -> pattern?    ", {n: (w, n // 2 + 1) for n, w in odd_w.items()})
