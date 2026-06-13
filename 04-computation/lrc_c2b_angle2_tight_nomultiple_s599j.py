"""C2b ANGLE 2: the clean reformulation. C' <=> 'a speed ≡0 mod n ⟹ loose (M>1/n)', i.e.
TIGHT (M=1/n) ⟹ NO multiple of n. Enumerate ALL tight configs over a window, check none has a
multiple of n. ANGLE 1b: characterize when the discrete 2n-1 witness works via SHELL collisions
(v_i+v_j≡0 mod 2n-1) shrinking the forbidden inverse-set. opus-2026-06-03-S599j."""
from itertools import combinations
from fractions import Fraction as F
from math import gcd
def nrm(a,m): a%=m; return min(a,m-a)
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for mm in range(d):
            t=F(mm,d); val=min(min((v*t)%1,1-(v*t)%1) for v in V)
            if val>best: best=val
            if best>=F(1,len(V)+1) and False: pass
    return best
def main():
    print("ANGLE 2 — TIGHT (M=1/n) ⟹ NO multiple of n ?  (the C' reformulation)")
    print(" n | #gcd1 configs | #tight (M=1/n) | #tight WITH a multiple of n (want 0)")
    for n in range(3,9):
        B=2*n; thr=F(1,n); tot=0; tight=0; tight_mult=0; egs=[]
        for V in combinations(range(1,B+1),n-1):
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            tot+=1; mm=M_exact(V)
            if mm==thr:
                tight+=1
                if any(v%n==0 for v in V): tight_mult+=1; egs.append(V)
        print(f" {n:2d} | {tot:6d} | {tight:4d} | {tight_mult}   {egs[:3]}")
    print("\nANGLE 1b — discrete 2n-1 witness works  <=>  shell collisions present?")
    print(" (forbidden j = {±v_i^{-1}}; collide iff v_i ≡ -v_k i.e. shell-partners, THM-401)")
    print(" n | #C2b | discrete-works & has-shell-pair | discrete-works & NO-shell-pair | fails&has | fails&no")
    for n in range(3,8):
        m=2*n-1; B=2*n
        a=b=c=d=0
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            # discrete witness works?
            works=any(min(nrm(v*j,m) for v in V)>=2 for j in range(1,m))
            # shell pair? v_i+v_k ≡0 mod m  (incl multiples thereof)
            shell=any((V[i]+V[k])%m==0 for i in range(len(V)) for k in range(i+1,len(V)))
            if works and shell: a+=1
            elif works and not shell: b+=1
            elif not works and shell: c+=1
            else: d+=1
        print(f" {n:2d} | {a+b+c+d:5d} | {a:4d} | {b:4d} | {c:4d} | {d:4d}")
if __name__=='__main__': main()
