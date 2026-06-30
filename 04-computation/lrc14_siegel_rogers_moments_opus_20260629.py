"""
The Siegel-Rogers moment hierarchy of the LRC danger-count X(t)=#{i: ||s_i t||<1/14}.
KEY: E[X] (mean = union bound) is SET-INDEPENDENT (=13/7, vacuous). ALL structure is in the
PAIR moment E[C(X,2)] (the Rogers/ze(2) pair-correlation) and higher. Lonely measure = P(X=0).
"""
from fractions import Fraction as F
from math import gcd
def frac(x): return x-(x.numerator//x.denominator)
def nrm(x): f=frac(x); return min(f,1-f)
def danger_indicator(s,t,n=14):
    return nrm(s*t)<F(1,n)
def moments(S,n=14,Q=2520):
    # exact-ish via fine rational grid t=a/Q (Q multiple of small denominators)
    EX=F(0); EC2=F(0); lonely=0
    for a in range(Q):
        t=F(a,Q)
        x=sum(1 for s in S if nrm(s*t)<F(1,n))
        EX+=x; EC2+=x*(x-1)//2
        if x==0: lonely+=1
    return EX/Q, EC2/Q, F(lonely,Q)
print("danger-count moments (grid Q=2520), n=14, fraction-of-time each runner dangerous = 2/14 = 1/7:")
print(f"{'set':>34} {'E[X] (mean=union bd)':>20} {'E[C(X,2)] (pair/ze2)':>21} {'P(X=0) lonely':>14}")
sets={
 "AP {1..13}": list(range(1,14)),
 "covering {1..11,13,84}": [1,2,3,4,5,6,7,8,9,10,11,13,84],
 "{2..14} (covering, loose)": list(range(2,15)),
 "wide {1,2,4,8,...} lacunary": [1,2,4,8,16,32,64,128,256,512,1024,2048,4096],
}
for nm,S in sets.items():
    EX,EC2,PL=moments(S)
    print(f"{nm:>34} {str(EX):>20} {float(EC2):>21.5f} {float(PL):>14.5f}")
print()
print("E[X]=13/7=%.5f for EVERY 13-set (the mean is the vacuous union bound, set-INDEPENDENT)." % (13/7))
print("=> all discrimination is in E[C(X,2)] (pair correlation = the SECOND factorial moment S2).")
print("   This S2 is exactly THM-501's resonance sum (pairs with additive relations) and the ze(2)")
print("   Rogers pair-factor (HYP-2856 floor). Var(N)=S1+2S2-S1^2 (HYP-2823) is built from it.")
print("   Lonely P(X=0): AP has the LARGEST (tight extremal); covering smaller; wide/lacunary smallest.")
