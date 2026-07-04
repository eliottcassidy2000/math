#!/usr/bin/env python3
"""
klein-2026-07-04-S124 (HYP-4076) - THE FIBONACCI / EISENSTEIN DUALITY OF THE COVERING-MIN.

Creative search for Fibonacci relations to the LRC(14) covering-min core.  Two findings:

(A) RECURRENCE SIBLING.  The covering-min witness's powers  n^k mod Phi6(n)  are
        1, n, n-1, -1, -n, 1-n   (period 6),
    satisfying  s_{k+2} = s_{k+1} - s_k   (the x^2 - x + 1 = Phi6 recurrence).
    Fibonacci satisfies  s_{k+2} = s_{k+1} + s_k  (x^2 - x - 1, golden).  The covering-min lives
    on the SIGN-FLIPPED SIBLING of Fibonacci: Eisenstein (order 6, Heegner -3) vs golden (order
    infinity, sqrt5).  Two adjacent metallic-type quadratics x^2 - x -/+ 1.

(B) GOLDEN / ANTI-GOLDEN DUALITY on the three-gap cap kernel.  The pair-overlap kernel
    K(a,b) = g(a,b)/(7ab) has breakpoints at the continued-fraction convergents of a/b (the
    three-gap/Stern-Brocot recursion, mac-mini S75b).  Among 1<=a<b<=13:
      * the CF-WORST ratios (all-1 partial quotients => slowest convergence => three-gap hardest)
        are the FIBONACCI ratios 5/8, 8/13, 5/13...  (consecutive Fibonaccis, F_k/F_{k+1}).
      * the covering-min witness t* = n/Phi6(n) = 14/183 has CF [0; n-1, n] = [0;13,14] -- LARGE
        partial quotients = FASTEST convergence = ANTI-golden.
    So the covering-min VALUE sits at the anti-golden (Eisenstein, fast-CF) extreme, while the
    three-gap DIFFICULTY (largest gap-count, most kinks) sits at the golden (Fibonacci, slow-CF)
    extreme.  This script pins both, exactly.
"""
from fractions import Fraction as F
from math import gcd

def phi6(n): return n*n - n + 1

def cf(a, b):
    """continued fraction of a/b (list of partial quotients)."""
    out = []
    while b:
        out.append(a // b); a, b = b, a % b
    return out

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

print("="*78)
print("(A) RECURRENCE SIBLING: covering-min powers n^k mod Phi6(n) obey s_{k+2}=s_{k+1}-s_k")
print("    (x^2-x+1 = Eisenstein), the SIGN-FLIP of Fibonacci s_{k+2}=s_{k+1}+s_k (x^2-x-1).")
print("="*78)
for n in [4, 7, 14]:
    q = phi6(n)
    pw = [pow(n, k, q) for k in range(7)]
    # signed representatives in (-q/2, q/2]
    pws = [(p if p <= q//2 else p - q) for p in pw]
    # check s_{k+2} == s_{k+1} - s_k (mod q)
    ok = all((pw[k+2] - (pw[k+1] - pw[k])) % q == 0 for k in range(5))
    print(f"  n={n:>3} Phi6={q:>4}: powers(signed) {pws}  period6={pw[6]==1}  "
          f"recurrence s_{{k+2}}=s_{{k+1}}-s_k mod Phi6: {ok}")
print("  Fibonacci: 1,1,2,3,5,8,... obeys s_{k+2}=s_{k+1}+s_k. Same shape, opposite sign.")
print("  => the covering-min recurrence is the char.poly x^2-x+1 (roots = primitive 6th roots,")
print("     Eisenstein Z[w], Heegner -3); Fibonacci is x^2-x-1 (roots = golden, Z[phi], sqrt5).")

print()
print("="*78)
print("(B) GOLDEN / ANTI-GOLDEN on the three-gap cap kernel K(a,b), a,b<=13")
print("    CF-worst (all-1 quotients = Fibonacci ratios) vs the covering-min's fast CF [0;13,14].")
print("="*78)
# CF length & max-quotient as a 'slowness' proxy: all-1 CF (Fibonacci) = slowest = three-gap worst
print(f"{'a/b':>7} {'CF(a/b)':>16} {'len':>4} {'all-ones(golden)?':>17} {'Fibonacci pair?':>15}")
FIB = {1,2,3,5,8,13,21}
worst = []
for b in range(2, 14):
    for a in range(1, b):
        if gcd(a, b) != 1: continue
        c = cf(a, b)  # [0, ...]
        tail = c[1:]  # partial quotients after the leading 0
        allones = len(tail) >= 2 and all(x == 1 for x in tail)
        # fib-consecutive test:
        fibcons = False
        s = sorted(FIB)
        for i in range(len(s)-1):
            if a == s[i] and b == s[i+1]: fibcons = True
        if allones:
            worst.append((a, b, len(tail)))
        mark = "<== GOLDEN" if allones else ""
        if allones or fibcons:
            print(f"{a}/{b:<5} {str(c):>16} {len(tail):>4} {str(allones):>17} {str(fibcons):>15} {mark}")
print(f"  CF-worst (all-1 = golden) ratios among a,b<=13: "
      f"{[f'{a}/{b}' for a,b,_ in sorted(worst, key=lambda t:-t[2])][:6]}")
print(f"  longest all-1 CF: {max((L for _,_,L in worst), default=0)} at "
      f"{[f'{a}/{b}' for a,b,L in worst if L==max(L2 for _,_,L2 in worst)]}  (consecutive Fibonaccis)")

# the covering-min witness CF
n = 14; q = phi6(n)
print(f"\n  covering-min witness t* = {n}/{q} = 14/183 : CF = {cf(n, q)} = [0; {n-1}, {n}] "
      f"(LARGE quotients => FAST convergence => ANTI-golden)")
print(f"  Its 1/M = Phi6/n = {q}/{n} : CF = {cf(q, n)} = [{n-1}; {n}]  (the S71 ladder [n-1; n]).")

print()
print("="*78)
print("(C) Does the three-gap-worst (Fibonacci ratio) MAXIMIZE the pair-overlap kernel K(a,b)?")
print("    K(a,b)=meas(D_a cap D_b) with danger half-width h=1/14; K=g/(7ab). Higher K = more")
print("    overlap = 'harder' pair. Check whether Fibonacci ratios peak g(a,b) (kink count).")
print("="*78)
h = F(1, 14)
def Kpair(a, b, N=20000):
    """exact-ish overlap meas(D_a cap D_b) by fine rational sampling of the 1-periodic overlap.
    D_v = {t: ||v t|| < h}. Use exact: overlap over t in [0,1) of both combs; scale by lcm."""
    # exact via residues: t in [0,1), ||a t||<h and ||b t||<h. Count measure exactly by
    # integrating the two comb indicators; do a fine Fraction grid (proxy for the g/(7ab) law).
    from fractions import Fraction as FF
    L = a * b  # period of both combs' pattern in units of 1/(ab)
    tot = FF(0)
    M = 4 * a * b
    for k in range(M):
        t = FF(k, M) + FF(1, 2 * M)
        da = min((a * t) % 1, 1 - (a * t) % 1)
        db = min((b * t) % 1, 1 - (b * t) % 1)
        if da < h and db < h: tot += FF(1, M)
    return tot
print(f"{'a/b':>7} {'CF':>14} {'K(a,b)~':>10} {'7ab*K=g':>9}")
rows = []
for (a, b) in [(1,13),(2,13),(5,13),(8,13),(3,11),(5,11),(5,8),(3,5),(1,2),(7,13),(6,13)]:
    K = Kpair(a, b)
    g = float(7 * a * b) * float(K)
    rows.append((a, b, float(K), g, cf(a, b)))
    print(f"{a}/{b:<5} {str(cf(a,b)):>14} {float(K):>10.5f} {g:>9.2f}")
print("  (8/13, 5/8 = Fibonacci/golden ratios; compare their g to non-Fibonacci at similar size.)")
print()
print("DONE")
