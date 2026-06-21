#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_residue_law_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

CRUCIAL DISCOVERY: e(11/10)=e(4/3), e(10/9)=e(3/2), e(12/11)=e(5/4).
HYPOTHESIS:  e_j = 7 r_j - pq  depends ONLY on (p mod 7, q mod 7)  (up to the cyclic
relabeling), i.e. e is a function of the residue pair (P7, Q7) = (p mod7, q mod7).

If TRUE, then sum_j|e_j| depends ONLY on (p mod 7, q mod 7), so it is bounded by a
FINITE table (at most 7*7 = 49 residue pairs, fewer after coprimality/window), and
the sharp constant 12 is just the max of sum|e|/p... no wait, sum|e| is residue-only
but 12p grows. Let's check: is sum|e| itself residue-only, OR sum|e|/something?

Actually if e is residue-only then sum|e| is BOUNDED (<= some constant ~44), independent
of p! Then D = sum|e|/(7pq) -> 0 like 1/(pq). That's even stronger. Let's test carefully.
"""
from math import gcd
from collections import defaultdict

P = 7

def evec(p, q):
    r = [0] * P
    for k in range(q):
        base = 14 * ((p * k) % q)
        for t in range(p):
            j = ((base + 2 * t + 1) // (2 * q)) % P
            r[j] += 1
    return tuple(7 * x - p * q for x in r)

def main():
    print("THREAD C: does e depend only on (p mod 7, q mod 7)?")
    print("=" * 72)
    # group ratios by (p%7, q%7) and see if e is constant within each group
    groups = defaultdict(list)
    for q in range(1, 120):
        for p in range(q + 1, (43 * q) // 20 + 1):
            if 20 * p > 43 * q or gcd(p, q) != 1:
                continue
            groups[(p % 7, q % 7)].append((p, q, evec(p, q)))
    residue_only = True
    sumabs_only = True
    for key, lst in sorted(groups.items()):
        evals = set(e for _, _, e in lst)
        svals = set(sum(abs(x) for x in e) for _, _, e in lst)
        if len(evals) > 1:
            residue_only = False
        if len(svals) > 1:
            sumabs_only = False
    print(f"e vector constant within (p%7,q%7) group: {'YES' if residue_only else 'NO'}")
    print(f"sum|e| constant within (p%7,q%7) group:  {'YES' if sumabs_only else 'NO'}")

    if not residue_only:
        # maybe e depends on (p%7,q%7) but with a SHIFT (cyclic). Check sum|e| only.
        print("\n(e itself varies -- likely cyclic relabel. Reporting sum|e| per residue pair.)")
    print("\nTABLE: (p%7,q%7) -> sum|e| values seen, sample ratios, max p in group")
    table = {}
    for key, lst in sorted(groups.items()):
        svals = sorted(set(sum(abs(x) for x in e) for _, _, e in lst))
        sample = lst[0]
        table[key] = svals
        print(f"  (p%7={key[0]}, q%7={key[1]}): sum|e| in {svals}  "
              f"e.g. {sample[0]}/{sample[1]} e={sample[2]}")

    # The bound sum|e| <= 12p: since sum|e| only depends on residues (if sumabs_only),
    # the SMALLEST p in each residue class is the binding case for sum|e|/p.
    print("\n" + "=" * 72)
    print("If sum|e| is residue-only (bounded), the sharp face sum|e|<=12p is tightest")
    print("at the SMALLEST p with that residue pair. Find max sum|e| / p_min per class:")
    print("=" * 72)
    worst = (0, 1, None)
    for key, lst in sorted(groups.items()):
        S = sum(abs(x) for x in lst[0][2])
        # smallest p in this group within window
        pmin = min(p for p, q, e in lst)
        # find the actual binding ratio (smallest p, which maximizes S/p)
        binding = min(lst, key=lambda t: t[0])
        Sb = sum(abs(x) for x in binding[2])
        from fractions import Fraction as Fr
        ratio = Fr(Sb, binding[0])
        if ratio > Fr(worst[0], worst[1]):
            worst = (Sb, binding[0], binding)
        print(f"  (p%7={key[0]},q%7={key[1]}): binding {binding[0]}/{binding[1]} "
              f"sum|e|={Sb} S/p={float(ratio):.4f}")
    from fractions import Fraction as Fr
    print(f"\nGLOBAL max sum|e|/p over all residue classes = {Fr(worst[0],worst[1])} "
          f"= {worst[0]/worst[1]:.5f}  at {worst[2][0]}/{worst[2][1]}  (target sup = 12)")

if __name__ == "__main__":
    main()
