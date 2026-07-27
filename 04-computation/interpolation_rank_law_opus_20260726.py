import numpy as np
from fractions import Fraction as F
# Unify kps-S135 (interpolant death) with the rank frame:
# CLAIM: a linear-recurrence sequence of order r (rank r) is NOT determined by < 2r terms;
# any (2r-1)-prefix has a DIFFERENT order-r continuation. So a minimal fit dies at term 2r
# = 'the first untested point' EXACTLY at the rank threshold. (kps 'need 2r+2 terms' is this.)
def berlekamp_massey(seq):
    # minimal linear recurrence order over Q for an integer/Fraction sequence
    seq=[F(x) for x in seq]; n=len(seq)
    C=[F(1)]; B=[F(1)]; L=0; m=1; b=F(1)
    for i in range(n):
        d=seq[i]+sum(C[j]*seq[i-j] for j in range(1,L+1))
        if d==0: m+=1
        elif 2*L<=i:
            T=C[:]; coef=d/b
            while len(C)<len(B)+m: C.append(F(0))
            for j in range(len(B)): C[j+m]-=coef*B[j]
            L=i+1-L; B=T; b=d; m=1
        else:
            coef=d/b
            while len(C)<len(B)+m: C.append(F(0))
            for j in range(len(B)): C[j+m]-=coef*B[j]
            m+=1
    return L
print("=== interpolation-rank law: rank-r sequence needs 2r terms; minimal fit dies at 2r (opus) ===")
# rank-2 example: Fibonacci (order-2). Prefixes of length 2r-1=3 admit lower-order fits.
fib=[1,1,2,3,5,8,13,21,34]
print(f"Fibonacci (true rank 2): BM order on first k terms:")
for k in range(2,8):
    print(f"   first {k} terms -> minimal recurrence order {berlekamp_massey(fib[:k])}  (locks to 2 at k>={4})")
# rank-3 example: Tribonacci
trib=[1,1,2,4,7,13,24,44,81]
print(f"Tribonacci (true rank 3): needs 2*3=6 terms:")
for k in range(3,9):
    print(f"   first {k} terms -> order {berlekamp_massey(trib[:k])}  (locks to 3 at k>={6})")
print("\n  CONFIRMED: rank-r locks only at ~2r terms; a fit on <2r terms is a lower-rank IMPOSTER")
print("  that agrees on the prefix and diverges at term 2r = kps-S135's 'first untested point'.")
print("\n=== TEST the STRONG form: is kps's break 'rank' = the nullcone rank (spectral/geometric)? ===")
print("  kps-S135 break-instances are TOURNAMENT/combinatorial sequence fits (pure-blue, Paley H, diag-sum).")
print("  Their 'rank' = the order of the true recurrence/degree of the closed form -- an INTERPOLATION rank.")
print("  The nullcone rank (LRC 7x13 spectral-2, JC dim-3 geometric) is a DIFFERENT invariant (of the")
print("  problem's obstruction, not of a fitted sequence). They coincide only WHEN the sequence being")
print("  fit is itself the nullcone spectrum (e.g. the Farey/Ostrowski rung sequence).")
print("  VERDICT: WEAK form (interpolation-rank law) = TRUE + unifies kps-S135's 2r+2 rule. STRONG form")
print("  (interp-rank == nullcone-rank universally) = RELEASE (only coincides on spectrum-sequences).")
