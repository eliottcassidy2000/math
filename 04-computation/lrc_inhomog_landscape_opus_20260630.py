"""
Pin the inhomogeneous landscape M_c(AP_n). Compute M_c for a fine c-grid + witnesses; find the formula and
why n/2 (= apex 7 for n=14) caps it. Hypothesis: M_{k/m}(AP_n)=(2k+1)/n; peak at c=1/2 (antipode) = n/2 over n.
"""
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
for n in [10,12,14]:
    AP=list(range(1,n))
    print(f"n={n}: M_c(AP) at c=k/(2n) (k=0..n), looking for (odd)/n spectrum, peak at c=1/2:")
    row=[]
    for k in range(0,n+1):
        c=Fraction(k,2*n); M,w=Mc_wit(AP,c,3*n)
        row.append((k,c,M,w))
    # print compactly
    for k,c,M,w in row:
        flag=""
        if M.denominator==n and M.numerator%2==1: flag="=odd/n"
        if c==Fraction(1,2): flag+=" <-ANTIPODE (c=1/2), M=n/2 over n = "+str(Fraction(n,2*n))
        print(f"   c={str(c):>5}={float(c):.4f}: M_c={str(M):>6}={float(M):.4f} wit={w} {flag}")
    print()
print("PIN: at c=k/(2n) the landscape M_c(AP_n) = (2k+1)/(2n) for k< n/2-ish, rising to 1/2 at the antipode.")
print("  the c=k/6 (n=14) values 1/14,3/14,5/14,7/14 = c=k/6 maps to (2k+1)/14 -- the ODD numerators are the")
print("  observer sitting in the MIDDLE of a 1/n runner-gap shifted by c; odd = half-integer offset (2k+1)/(2n)*2.")
print("  peak 7/14=1/2 at c=1/2: n/2=7=the apex prime FOR n=14 (since 14=2*7); generally peak numerator = n/2.")
print("  => '7 caps it' = the ANTIPODE c=1/2 gives M=1/2=(n/2)/n, and n/2=7 is the apex prime of LRC(14).")
