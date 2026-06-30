#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CV(H)^2 over labeled tournaments: the EXACT even-run closed form  (klein-2026-06-29-S5).

The owner's identity: CV(H)^2 = (1/n!) sum_{pi': no descending consecutive adjacency} 2^{j(pi')} - 1,
j = #ascending consecutive adjacencies. mac-mini (HYP-3560) proved CV(H)^2 -> 0 via Poisson(1).

NEW (this script): the EXACT structure behind it. Inclusion-exclusion on descending successions
(binomial inversion: forcing a directed-linear-forest L of consecutive-adjacencies gives (n-|L|)!
permutations) collapses the permutation sum to an EDGE-SUBSET sum on the integer path 1-2-...-n.
Each maximal run of chosen consecutive edges contributes a per-run factor (1 + (-1)^m) = 2 if its
length m is EVEN, 0 if ODD (ascending vs descending orientation cancel for odd runs). Hence

    A_n(2) := sum_{pi'} c(pi') 2^{j(pi')}
            = sum_{ edge-subsets S of [1..n-1] with ALL maximal runs EVEN } 2^{#runs(S)} (n - |S|)!,

a PARITY phenomenon: two Hamiltonian paths' shared consecutive-integer arcs contribute only when the
overlap runs are EVEN-length (odd overlap cancels by orientation). Leading term (one length-2 run):
CV(H)^2 ~ 2/n exactly. Verified vs the permutation form (n<=9) and direct tournament enumeration (n<=5).
"""
from __future__ import annotations
import itertools
from fractions import Fraction as F
from math import comb, factorial

# ---- (1) permutation form: A_n(2) = sum_{no desc consec adj} 2^{#asc consec adj} ----
def A_perm(n):
    tot = 0
    for sigma in itertools.permutations(range(1, n+1)):
        desc = asc = 0
        ok = True
        for a, b in zip(sigma, sigma[1:]):
            if b == a - 1:
                ok = False; break
            if b == a + 1:
                asc += 1
        if ok:
            tot += 2**asc
    return tot

# ---- (2) even-run DP: A_n(2) = sum_t W(n-1,t)(n-t)!, W from the parity-cancellation ----
def W_polys(m):
    """polynomials (lists, index=t) for states after m edges; 2^#runs baked in, runs must end even."""
    # state g (gap), o (run-odd len), e (run-even len)
    g, o, e = [1], [0], [0]
    def add(p, q):
        r = [0]*max(len(p), len(q))
        for i, c in enumerate(p): r[i] += c
        for i, c in enumerate(q): r[i] += c
        return r
    def shift(p):  # multiply by x (advance t by 1)
        return [0] + p
    for _ in range(m):
        ng = add(g, e)                       # edge0 from g or e -> gap (x^0)
        no = add(shift([2*c for c in g]),    # edge1 from g -> run-odd, NEW run weight 2
                 shift(e))                   # edge1 from e -> run-odd
        ne = shift(o)                        # edge1 from o -> run-even ; (edge0 from o invalid)
        g, o, e = ng, no, ne
    return add(g, e)                         # accept gap or run-even

def A_evenrun(n):
    W = W_polys(n-1)
    return sum(W[t]*factorial(n-t) for t in range(len(W)) if t <= n)

# ---- (3) direct over labeled tournaments (cross-check, small n) ----
def H_count(n, bits, pairs, idx):
    cnt = 0
    for perm in itertools.permutations(range(n)):
        good = True
        for a, b in zip(perm, perm[1:]):
            i, j = (a, b) if a < b else (b, a)
            bit = (bits >> idx[(i, j)]) & 1     # 1 means i beats j
            beats = (i if bit else j)           # who wins edge (i,j)
            if beats != a:                      # need a -> b
                good = False; break
        if good: cnt += 1
    return cnt

def CV2_direct(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    d = len(pairs); S = E = E2 = 0
    for bits in range(2**d):
        h = H_count(n, bits, pairs, idx)
        E += h; E2 += h*h; S += 1
    EH = F(E, S); EH2 = F(E2, S)
    return EH2/(EH*EH) - 1

if __name__ == "__main__":
    print("="*78)
    print(" CV(H)^2 over labeled tournaments -- exact even-run closed form (klein-S5)")
    print("="*78)
    print(" cross-checks (perm form vs even-run DP vs direct tournament enumeration):")
    print(f"   {'n':>2} {'A_perm':>10} {'A_evenrun':>10} {'CV^2 (frac)':>14} {'CV^2':>9} {'n*CV^2':>8}")
    for n in range(2, 10):
        ap = A_perm(n)
        ae = A_evenrun(n)
        assert ap == ae, f"perm vs even-run mismatch n={n}: {ap} {ae}"
        cv2 = F(ae, factorial(n)) - 1
        flag = ""
        if n <= 5:
            cd = CV2_direct(n)
            flag = "  direct=OK" if cd == cv2 else f"  direct MISMATCH {cd}"
        print(f"   {n:>2} {ap:>10} {ae:>10} {str(cv2):>14} {float(cv2):>9.5f} {float(n*cv2):>8.4f}{flag}")

    print("\n A_n(2) sequence (for OEIS):",
          [A_evenrun(n) for n in range(1, 12)])
    print("\n n*CV(H)^2 -> 2  (the dominant fluctuation is a single length-2 overlap run):")
    for n in [10, 20, 40, 80, 160, 320]:
        cv2 = F(A_evenrun(n), factorial(n)) - 1
        print(f"   n={n:>4}:  CV^2 = {float(cv2):.6f}   n*CV^2 = {float(n*cv2):.5f}")
    print("\n => CV(H)^2 ~ 2/n exactly: H concentrates, std/mean ~ sqrt(2/n). The even-run parity")
    print("    (odd overlap runs cancel) is the S_n-side mirror of the LRC 2-adic descent (THM-580).")
