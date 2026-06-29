#!/usr/bin/env python3
"""
signed-cycle-index-metagraph-spectrum.py   (klein-2026-06-29-S2)

Closes HYP-3540 with an explicit closed form, and turns THM-584 into a tool.

THM-584: the iso-class metagraph adjacency A is the S_n-quotient of the arc-hypercube
Q_d (d=C(n,2)); eigenvalues are d-2k with multiplicity mult(k) = dim of S_n-invariants
at level k. complement R = antipodal = flips all bits, acting (-1)^k on level k.

mac-mini (HYP-3543) + klein: the S_n action on Q_d is SIGNED (a vertex swap (i j)
reverses the pair {i,j}), so mult(k) is the per-level evaluation of the CYCLE INDEX of
this signed S_n action (the vertex-induced subgroup of B_{C(n,2)}).

CLOSED FORM (klein-S2):  for sigma in S_n, its action on the d=C(n,2) unordered pairs
decomposes into signed cycles (length ell, orientation sign s_c = product of per-step
order-reversals around the cycle). A sigma-invariant subset S of pairs is a union of whole
cycles; on the level-k (symmetric-subset / Boolean-Fourier) module, including a cycle
contributes sign s_c at degree ell. Hence the PER-LEVEL SIGNED CYCLE INDEX is

    P_n(x) = sum_k mult(k) x^k = (1/n!) * sum_{sigma in S_n}  prod_{cycles c} (1 + s_c x^{ell_c}).

(NB: this is the SYMMETRIC-subset trace, not the exterior-power det(I+xM_sigma); the two
differ by reordering signs and only the former matches the verified metagraph spectrum.)

Evaluations (verified):
    P_n(1)  = A000568(n)   (total iso classes; = the all-odd-cycle Burnside, CLAUDE.md)
    P_n(-1) = SC(n)        (self-converse count; the ANTIPODAL Lefschetz/Euler number)
    dim R-even = (P_n(1)+P_n(-1))/2 = V_merged ;  dim R-odd = (P_n(1)-P_n(-1))/2 = #NS.

PAYOFF: P_n needs only n! permutations, NOT 2^{C(n,2)} tournaments -- it delivers the full
metagraph spectrum (eigenvalues d-2k with multiplicities) at n where the metagraph itself
is far too large to enumerate (n=7: 2^21; n=8: 2^28; n=9: 2^36).
"""
import itertools
from math import comb, factorial
from fractions import Fraction

A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:9733056}

# ---- integer polynomial helpers (coeff lists, index = power of x) ----
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] += ai*bj
    return r

def signed_pair_cycles(perm, n):
    """Decompose sigma's action on unordered pairs into (length, sign) signed cycles.
       sign = product of per-step orientation flips around the cycle (+1 preserve, -1 reverse)."""
    pair_id = {}
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pair_id[(i, j)] = len(pairs)
            pairs.append((i, j))
    # step map: pair (i<j) -> (canonical target, step_sign)
    def step(p):
        i, j = p
        a, b = perm[i], perm[j]
        if a < b:
            return (a, b), +1
        else:
            return (b, a), -1
    seen = [False]*len(pairs)
    cycles = []
    for start in range(len(pairs)):
        if seen[start]:
            continue
        p = pairs[start]
        s = 1
        length = 0
        while True:
            idx = pair_id[p]
            if seen[idx]:
                break
            seen[idx] = True
            p, st = step(p)
            s *= st
            length += 1
        cycles.append((length, s))
    return cycles

def poly_for_sigma(perm, n):
    """Level generating polynomial sum_k chi_k(sigma) x^k for sigma's action on the
       level-k subspaces span{chi_S:|S|=k} of the SYMMETRIC-SUBSET (Boolean-Fourier)
       module -- NOT the exterior power. A subset S invariant under sigma is a union of
       whole pair-cycles; including a cycle (length ell, orientation sign s_c) contributes
       sign s_c at degree ell, so the per-cycle factor is (1 + s_c x^ell)."""
    poly = [1]
    for (ell, s) in signed_pair_cycles(perm, n):
        fac = [0]*(ell+1)
        fac[0] = 1
        fac[ell] = s                # 1 + s_c * x^ell   (s_c = +-1)
        poly = pmul(poly, fac)
    return poly

def signed_cycle_index(n):
    """P_n(x) = (1/n!) sum_sigma det(I+xM_sigma); returns integer coeff list (mult(k))."""
    d = comb(n, 2)
    total = [0]*(d+1)
    for perm in itertools.permutations(range(n)):
        p = poly_for_sigma(perm, n)
        for i, c in enumerate(p):
            total[i] += c
    nf = factorial(n)
    mult = []
    for c in total:
        assert c % nf == 0, f"non-integer multiplicity at n={n}: {c}/{nf}"
        mult.append(c // nf)
    return mult

def report(n):
    d = comb(n, 2)
    mult = signed_cycle_index(n)
    total = sum(mult)
    P1 = sum(mult)                       # P_n(1)
    Pm1 = sum(((-1)**k)*m for k, m in enumerate(mult))   # P_n(-1)
    Veven = (P1+Pm1)//2
    Vodd = (P1-Pm1)//2
    even_levels = [mult[k] for k in range(0, d+1, 2)]
    odd_levels  = [mult[k] for k in range(1, d+1, 2)]
    print(f"\n n={n}  d=C(n,2)={d}")
    print(f"   level multiplicities mult(k), k=0..{d}: {mult}")
    print(f"   P_n(1)  = {P1:>9}   A000568({n}) = {A000568[n]:>9}   {'OK' if P1==A000568[n] else 'MISMATCH'}")
    print(f"   P_n(-1) = {Pm1:>9}   = SC(n) (self-converse count) = ANTIPODAL Euler number")
    print(f"   R-even dim (even levels) = {Veven} = V_merged ;  R-odd dim (odd levels) = {Vodd} = #NS pairs")
    print(f"   even-level mults: {even_levels}")
    print(f"   odd-level  mults: {odd_levels}")
    # metagraph spectrum: eigenvalue d-2k with multiplicity mult(k)
    spec = [(d-2*k, mult[k]) for k in range(d+1) if mult[k] > 0]
    print(f"   METAGRAPH SPECTRUM (eigenvalue: mult): " +
          ", ".join(f"{ev}:{m}" for ev, m in spec))
    return Pm1

# ground truth from the verified block spectra (r-block-spectra-antipodal.py, n<=6)
EXPECT = {
    4: [1,0,1,1,1,0,0],
    5: [1,0,1,1,4,1,3,0,1,0,0],
    6: [1,0,1,1,5,5,10,8,12,6,4,2,1,0,0,0],
}

if __name__ == "__main__":
    print("="*74)
    print(" Signed cycle index P_n(x) = per-level metagraph multiplicities (HYP-3540)")
    print("="*74)
    sc = {}
    for n in range(3, 9):     # n=8 reachable (8!=40320); enumeration would need 2^28
        mult_check = signed_cycle_index(n)
        if n in EXPECT:
            ok = (mult_check == EXPECT[n])
            print(f"\n [cross-check n={n} vs verified block spectrum: "
                  f"{'MATCH' if ok else 'MISMATCH '+str(mult_check)+' vs '+str(EXPECT[n])}]")
        sc[n] = report(n)
    print("\n" + "="*74)
    print(" SELF-CONVERSE COUNT SC(n) = P_n(-1) = antipodal Euler number:")
    print("   n=3..8:", [sc[n] for n in range(3, 9)])
    print("   (CLAUDE.md / mac-mini gave SC=2,2,8,12 for n=3..6; this EXTENDS to n=7,8)")
    print("   cross-check via NS-merged 0,1,2,22,184 (CLAUDE.md n=3..7): SC = A000568 - 2*NS")
    for n, ns in zip(range(3, 8), [0, 1, 2, 22, 184]):
        print(f"     n={n}: A000568-2*{ns} = {A000568[n]-2*ns}  (P_n(-1)={sc[n]})  "
              f"{'OK' if A000568[n]-2*ns==sc[n] else 'MISMATCH'}")
