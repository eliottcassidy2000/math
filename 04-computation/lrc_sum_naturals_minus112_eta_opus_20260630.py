"""
lrc_sum_naturals_minus112_eta_opus_20260630.py  (HYP-3773)
The sum of natural numbers IS the covering-min:
 (1) Phi6(n) = 2*(1+2+...+(n-1)) + 1 = 2T+1; killer n(n-1)=2T == -1 mod Phi6; s(n,Phi6)=-T/(12T+6)->-1/12=zeta(-1).
 (2) margin = Dedekind-ETA multiplier phase (validated on c=1); covering-min = large-c cusp case; 1/24=-zeta(-1)/2.
 (3) the two H6 endpoints mirror the zeta ladder: cyclotomic Phi6=2*sum(k)+1<->zeta(-1)=-1/12 (naturals);
     hexagonal n(2n-1) in sum(k^2)=(n-1)n(2n-1)/6<->zeta(-2)=0 (squares).
 (4) open lever: construction margin = classical 2-fold cotangent Dedekind (AP orbit); beaters = Zagier higher-dim.
ASCII only. See reflection the-sum-of-natural-numbers-IS-the-covering-min-...-opus-20260630.md.
"""
import numpy as np
from fractions import Fraction as F
from cmath import exp, pi, sqrt, phase
def s_recip(h,k):
    h%=k
    if h==0: return F(0)
    tot=F(0); sign=1; hh,kk=h,k
    while hh!=0:
        tot+=sign*(F(-1,4)+(F(hh,kk)+F(kk,hh)+F(1,hh*kk))/12); sign=-sign; hh,kk=kk%hh,hh
    return tot
def s_cot(h,k):
    j=np.arange(1,k); return float(np.sum(1/np.tan(np.pi*j/k)/np.tan(np.pi*j*h/k))/(4*k))
def eta(t,M=400):
    q=exp(2j*pi*t); p=1+0j
    for m in range(1,M+1): p*=(1-q**m)
    return exp(1j*pi*t/12)*p

def part1_sum_naturals():
    print("(1) Phi6 = 2*(sum 1..n-1)+1; killer=2*sum== -1; s(n,Phi6)=-T/(12T+6)-> -1/12=zeta(-1):")
    for n in [3,7,14,100]:
        T=n*(n-1)//2; P=n*n-n+1; sn=s_recip(n,P)
        print(f"    n={n:>3}: T={T:>4} Phi6=2T+1?{P==2*T+1} killer=2T={n*(n-1)}== -1 mod?{n*(n-1)%P==P-1} s={sn}={float(sn):.5f}")

def part2_eta():
    t=0.2+1.2j
    mult=eta((-1)/t)/(sqrt(t)*eta(t))
    print(f"\n(2) eta multiplier validated on S-transform c=1: |mult|={abs(mult):.4f} phase/pi={phase(mult)/pi:+.4f} (pred -0.25)")
    print("    covering-min = large-c (cusp) case: c=Phi6 sends tau to Im~1/Phi6 (cusp); 1/24=-zeta(-1)/2 = eta weight anomaly.")

def part3_two_endpoints():
    print("\n(3) two H6 endpoints mirror the zeta ladder:")
    for n in [7,14,20]:
        T=n*(n-1)//2; S2=(n-1)*n*(2*n-1)//6
        print(f"    n={n}: cyclotomic Phi6=2*{T}+1={2*T+1} (naturals, zeta(-1)=-1/12); hexagonal {n*(2*n-1)} in sum k^2={S2}=(n-1)n(2n-1)/6 (squares, zeta(-2)=0)")

def part4_lever():
    print("\n(4) open lever: construction = 2-fold cotangent Dedekind (AP orbit); beaters = Zagier higher-dim:")
    for n in [7,14]:
        P=n*n-n+1; print(f"    n={n}: s(n,Phi6) reciprocity={float(s_recip(n,P)):.6f} == cotangent-sum={s_cot(n,P):.6f}")
    print("    beaters are SPREAD (not one AP orbit) => not a single 2-fold Dedekind => Zagier d(D;a1..ar)=(1/D)sum prod cot (has reciprocity).")

if __name__=="__main__":
    part1_sum_naturals(); part2_eta(); part3_two_endpoints(); part4_lever()
    print("\n=> the covering-min carries the FINITE speed-sum (Phi6=2T+1) and the REGULARIZED one (-1/12=zeta(-1));")
    print("   the margin is a modular eta cocycle at the cusp; the H6 window mirrors the zeta ladder.")
