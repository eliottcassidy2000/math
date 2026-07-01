"""
lrc_margin_reciprocity_descent_opus_20260630.py
The LRC construction margin as a FINITE RECIPROCITY COMPUTATION (HYP-3770; builds on HYP-3768 mac-mini/klein,
HYP-3769 opus). margin(n)=(n-1)/(n*Phi6)=-12 s(n,Phi6)/n^2 (Dedekind sum). Reciprocity computes s(h,k) in
O(log k) Euclidean descent -- no enumeration of the k residues. Push to large-lcm bases (infeasible enumeration).
 (A) reciprocity == direct enumeration (verify);
 (B) margin(n) via reciprocity, 2 steps (Phi6==1 mod n), exact to n=10^30; s(n,Phi6)->-1/12=zeta(-1);
 (C) large-lcm base D=lcm(1..m)~e^m: O(log D) descent;
 (D) NEGATIVE: general-rung margin is NOT a single observer Dedekind sum (only the construction rung n is).
ASCII only. See reflection the-LRC-margin-is-a-finite-reciprocity-computation-...-opus-20260630.md.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm

def s_direct(h,k):                                   # O(k) enumeration (the infeasible route)
    tot=F(0)
    for i in range(1,k):
        r=(h*i)%k; tot += (F(i,k)-F(1,2))*(F(0) if r==0 else F(r,k)-F(1,2))
    return tot
def s_recip(h,k):                                    # O(log k) reciprocity descent (the finite route)
    h%=k; tot=F(0); sign=1; steps=0; hh,kk=h,k
    while hh!=0:
        steps+=1; tot += sign*(F(-1,4)+(F(hh,kk)+F(kk,hh)+F(1,hh*kk))/12); sign=-sign; hh,kk=kk%hh,hh
    return tot,steps
def nextprime(m):
    x=m+1
    while any(x%p==0 for p in range(2,int(x**0.5)+1)): x+=1
    return x

def A_verify():
    print("(A) reciprocity == direct enumeration:")
    for h,k in [(2,5),(3,7),(8,21),(14,183),(23,100)]:
        sr,st=s_recip(h,k); print(f"    s({h},{k}): recip={sr} ({st} steps) direct={s_direct(h,k)} match={sr==s_direct(h,k)}")
def B_large_n():
    print("\n(B) construction margin = -12 s(n,Phi6)/n^2 via reciprocity (2 steps, Phi6==1 mod n); large n:")
    for n in [14, 10**6, 10**18, 10**30]:
        P=n*n-n+1; sr,st=s_recip(n,P); margin=-12*sr/(n*n)
        print(f"    n={n if n<100 else '10^'+str(len(str(n))-1)}: margin={float(margin):.3e} EXACT in {st} steps; enum={P} residues; s(n,Phi6)~{float(sr):.6f} (->-1/12=zeta(-1))")
def C_lcm_base():
    print("\n(C) large-lcm base D=lcm(1..m)~e^m: O(log D) descent, enumeration infeasible:")
    for m in [20,50,80]:
        D=reduce(lcm,range(1,m+1)); a=nextprime(m); sr,st=s_recip(a,D)
        print(f"    m={m}: D~10^{len(str(D))-1} ({D.bit_length()} bits), enum={D:.1e} terms INFEASIBLE; s({a},D) in {st} steps")
def D_general_rung():
    print("\n(D) NEGATIVE: is the general-rung margin a single observer Dedekind sum? n=14, D=13k+1:")
    n=14
    for k in [2,3,7,14]:
        D=k*(n-1)+1; target=-F(n*(k-1),12*D)
        hits=[a for a in range(1,D) if gcd(a,D)==1 and s_direct(a,D)==target]
        print(f"    k={k:>2} D={D:>3}: margin={F(k-1,n*D)} -> observer a with s(a,D)=target: {hits[:4] if hits else 'NONE'}  (a=n works: {n in hits})")
    print("    => only the construction rung k=n (D=Phi6) is a single observer Dedekind sum; covering-min a(n) is not.")

if __name__=="__main__":
    A_verify(); B_large_n(); C_lcm_base(); D_general_rung()
    print("\n=> LRC construction margin = finite O(log) reciprocity at ANY base (enumeration infeasible);")
    print("   residual(HYP-3769)=CF tail and Dedekind sum(HYP-3768)=CF descent are ONE Euclidean descent.")
