"""Confirm: the clean odd spectrum M_{k/(p-1)}(AP_2p)=(2k+1)/(2p) holds iff n=2p (p prime), capped at p=apex.
Verify n=2p (6,10,14,22) vs non-2p (12,16); witnesses at multiples of p."""
import math
from fractions import Fraction
def Mc_wit(S,c,Qmax):
    best=Fraction(-1); bw=None
    for q in range(1,Qmax+1):
        for a in range(0,q):
            mind=None
            for s in S:
                x=Fraction(s*a,q)-c; d=x-x.__floor__(); d=min(d,1-d)
                if mind is None or d<mind: mind=d
            if mind>best: best=mind; bw=(a,q)
    return best,bw
def isprime(p): return p>1 and all(p%d for d in range(2,int(p**.5)+1))
for n in [6,10,14,22,12,16,18]:
    AP=list(range(1,n)); m=(n-2)//2
    p=n//2
    twop = (n%2==0 and isprime(p))
    # check all clean points
    allmatch=True; wits=[]
    for k in range(0,m//2+1):
        c=Fraction(k,m); M,w=Mc_wit(AP,c,5*n)
        if M!=Fraction(2*k+1,n): allmatch=False
        if 0<k: wits.append(w[1])  # witness denominators
    print(f"n={n}={'2*'+str(p)+(' (p PRIME)' if twop else ' (p NOT prime)') if n%2==0 else 'odd'}: "
          f"clean odd spectrum (2k+1)/n at c=k/{m} holds={allmatch}; witness denoms={wits}")
print()
print("PINNED (honest, corrected):")
print("  * the ODD SPECTRUM M_{k/(p-1)}(AP_2p) = (2k+1)/(2p) holds IFF n=2p with p PRIME (verified 6,10,14,22;")
print("    FAILS for n=12=2^2*3, 16=2^4, 18=2*9). It is an LRC(2p) feature, NOT general.")
print("  * the witnesses are MULTIPLES OF p (the apex prime): e.g. n=14 -> q in {21,42}=3p,6p. The apex p")
print("    structures the optimal lonely time.")
print("  * THE CAP: the antipode c=1/2 gives M=p/(2p)=1/2; the top numerator is p=n/2 = THE APEX PRIME. So")
print("    '7 caps the n=14 spectrum' because 7=n/2 IS the apex prime of LRC(14) -- the spectrum literally")
print("    counts off the odds 1,3,5,7=p, and the cap is the apex.")
print("  * CORRECTION of my earlier reflection: the {1,3,5,7}/14 spectrum is SPECIAL to n=2p (not a general")
print("    landscape feature); between the clean c=k/(p-1) points M_c is messy.")
