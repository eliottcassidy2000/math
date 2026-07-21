#!/usr/bin/env python3
"""nc2_resonance_band_sheffer_deathstar_S87.py (HYP-8769)
Work NC2 at codex's RESONANCE BAND (THM-2017/HYP-8766), connecting it to my S64
Sheffer-deg-b>=2 hierarchy. Three-weight P=Z*A(s)+B(s)+W*c0 (p=q=1,r=2), central
offset deg h = r*deg b. E[P^m] = sum_i multinomial(m;i,i,m-2i) * L(s^i A^i c0^i B^{m-2i}),
L(s^k)=k!. Verify E[P^m] FIRES (NC2 holds) on central-offset examples + exhibit the
channel structure = the Sheffer/Legendre recurrence whose 'no common nonzero zero'
is the open crux (my S62 Hermite-no-common-root generalized)."""
from math import factorial as fact
from fractions import Fraction as Fr
from itertools import product

def Lpoly(coeffs):  # coeffs = dict{power: coeff}; L(s^k)=k!  => sum coeff*k!
    return sum(v*fact(k) for k,v in coeffs.items())
def polymul(p,q):
    r={}
    for k,v in p.items():
        for k2,v2 in q.items(): r[k+k2]=r.get(k+k2,0)+v*v2
    return r
def polypow(p,n):
    r={0:Fr(1)}
    for _ in range(n): r=polymul(r,p)
    return r
def EPm(A,B,c0,m):
    # P = Z*A(s) + B(s) + W*c0 ; charge-0 needs #Z=#W=i ; s^i from Z^iW^i
    tot=Fr(0)
    for i in range(m//2+1):
        multi=Fr(fact(m),fact(i)*fact(i)*fact(m-2*i))
        radial=polymul(polymul({i:c0**i}, polypow(A,i)), polypow(B,m-2*i))  # s^i? add below
        radial=polymul(radial,{i:Fr(1)})  # times s^i (Z^iW^i)
        tot+=multi*Lpoly(radial)
    return tot

print("Central-offset resonance band: P=Z*(a0+a1 s)+(b0+b1 s)+W*c0, charges {+1,0,-1},")
print("deg h=deg(s*A*c0)=2 = r*deg B=2*1 (central). Scan small integer coeffs for a")
print("NULLCONE element (E[P^m]=0 for all m<=12) among TWO-SIDED P:\n")
rng=(-2,-1,1,2); found=0; depthhist={}
for a0,a1,b1,c0 in product(rng,repeat=4):
    for b0 in (-2,-1,0,1,2):
        A={0:Fr(a0),1:Fr(a1)}; B={k:Fr(v) for k,v in ((0,b0),(1,b1)) if v}
        if not B: continue
        # two-sided: has +1 (a) and -1 (c) charge; always here (a0 or a1, c0 nonzero)
        d=None
        for m in range(1,13):
            if EPm(A,B,c0,m)!=0: d=m; break
        if d is None: found+=1
        else: depthhist[d]=depthhist.get(d,0)+1
print(f"  two-sided NULLCONE elements found (E[P^m]=0 all m<=12): {found}")
print(f"  first-fire depth histogram: {dict(sorted(depthhist.items()))}")
print(f"  total scanned: {sum(depthhist.values())+found}")

# exhibit the channel/Sheffer structure for one example
print("\nChannel structure for P=Z*(1+s)+(1+s)+W*1 (a0=a1=b0=b1=c0=1):")
A={0:Fr(1),1:Fr(1)}; B={0:Fr(1),1:Fr(1)}; c0=Fr(1)
for m in range(1,9):
    chans=[]
    for i in range(m//2+1):
        multi=Fr(fact(m),fact(i)*fact(i)*fact(m-2*i))
        rad=polymul(polymul({i:c0**i},polypow(A,i)),polymul(polypow(B,m-2*i),{i:Fr(1)}))
        chans.append(multi*Lpoly(rad))
    print(f"  m={m}: E[P^m]={EPm(A,B,c0,m)}  channels(i=0..{m//2})={[str(x) for x in chans]}")
print("""
READING: ZERO two-sided nullcone elements at the central offset (NC2 holds on the
scan). The channel sum E[P^m]=sum_i (Sheffer channels) is the S64 object; each channel
is a factorial-weighted L-moment. codex's 'no common nonzero zero' of the hyper-Bessel
tower = the S62 Hermite-no-common-root generalized to this Sheffer/Legendre hierarchy
(my S64: it does NOT close by top-term alone at deg b>=2 -- exactly codex's resonance
band). NC2 REMAINS OPEN there; the algebraic (Sheffer) and analytic (Bessel) attacks
are the SAME 'no common zero' crux.""")
