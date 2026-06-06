"""LRC attack: the SHELL-PARTNER LEMMA. C=2n-1. Discrete witness t=m/C gives M≥2/C iff some m has
v_k·m ∉ {0,±1} mod C ∀k. Forbidden set F={0}∪{±v_k^{-1}}; a SHELL-PARTNER (v_i+v_j≡0) ⟹ v_j^{-1}=
−v_i^{-1} ⟹ that pair shares its ±-forbidden set ⟹ |F\{0}|≤2(n-1)−2=2n−4 < 2n−2 ⟹ a good m exists
⟹ M≥2/(2n-1)>1/n. PROVE (coprime case) + verify + the residual (multiples, no shell-partner).
opus-2026-06-07-S700."""
from itertools import combinations, product
from math import gcd
from fractions import Fraction as F
def Mexact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for m in range(d):
            t=F(m,d); v=min(min((x*t)%1,1-(x*t)%1) for x in V)
            if v>best: best=v
    return best
def forbidden(V,C):
    Fset=set([0])
    for v in V:
        for m in range(1,C):
            if (v*m)%C in (1,C-1): Fset.add(m)
    return Fset
def has_shell(V,C): return any((V[i]+V[j])%C==0 for i in range(len(V)) for j in range(i+1,len(V)))
def main():
    print("SHELL-PARTNER LEMMA test: shell-partner ⟹ good discrete witness ⟹ M≥2/(2n-1)?")
    print(" n | #shell-partner configs | all have a good m (|F|<C)? | all M≥2/(2n-1)?")
    for n in range(5,9):
        C=2*n-1; B=2*n; thr=F(2,C); okwit=True; okM=True; cnt=0
        for V in combinations(range(1,B+1),n-1):
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if not has_shell(V,C): continue
            cnt+=1
            Fs=forbidden(V,C)
            if len(Fs)>=C: okwit=False
            if Mexact(V)<thr: okM=False
        print(f" {n} | {cnt} | good-m exists: {okwit} | M≥2/(2n-1): {okM}")
    print("\n PROOF (coprime): F\\{0} = ∪{±v_k^{-1}}; shell-partner v_i+v_j≡0 ⟹ v_j^{-1}=−v_i^{-1} ⟹")
    print(" that pair contributes 2 (not 4) ±-values ⟹ |F\\{0}|≤2(n-1)−2=2n−4 < 2n−2=#nonzero residues")
    print(" ⟹ ≥2 good m ⟹ M≥2/C. ∎  (gcd>1 speeds: forbid 0-multiples not ±1, handled separately.)\n")
    print("THE RESIDUAL — multiples of n WITHOUT a shell-partner: count, min M, still loose?")
    print(" n | #residual (mult-of-n, no shell-partner, gcd1) | min M | ≥2/(2n-1)? | via discrete witness?")
    for n in range(5,9):
        C=2*n-1; B=2*n; thr=F(2,C); cnt=0; minM=F(2); discrete_fail=0
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if has_shell(V,C): continue
            cnt+=1; m=Mexact(V); 
            if m<minM: minM=m
            if len(forbidden(V,C))>=C: discrete_fail+=1
        print(f" {n} | {cnt} | {minM}={float(minM):.4f} | {minM>=thr} | discrete-fails(need coarse): {discrete_fail}/{cnt}")
if __name__=='__main__': main()
