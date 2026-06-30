"""
Attack the off-cusp covering-min (building on klein HYP-3705/3706). (1) verify the rotation tower
q=n-1=n^2 mod Phi_6(n) (Singer multiplier = 3-fold hexagonal rotation). (2) construction {1..12,182}:
M=14/183, optimal witness a*, speed spread mod 183. (3) is 14/183 the covering-min? scan covering sets.
"""
import math
from fractions import Fraction
n=14; Phi6=n*n-n+1  # 183
print(f"n={n}, Phi_6(n)=n^2-n+1={Phi6}=3*61")
# (1) rotation tower mod Phi6
print("ROTATION TOWER mod Phi_6 (klein HYP-3706): zeta_6=*n (ord6), zeta_3=*(n-1)=*n^2 (ord3=Singer mult), -1=*n^3 (ord2)")
for name,mult in [("zeta_6 = *n",n),("zeta_3 = *n^2",n*n%Phi6),("*(n-1)",n-1),("-1 = *n^3",pow(n,3,Phi6))]:
    o=1; x=mult%Phi6
    while x!=1 and o<10: x=x*mult%Phi6; o+=1
    print(f"   {name:>16} = {mult%Phi6:>4} mod {Phi6}: order {o}")
print(f"   => n^2 mod Phi6 = {n*n%Phi6} = n-1 = {n-1} (the Singer multiplier = 3-fold rotation): {n*n%Phi6==n-1}")
# (2) construction witness
def M_and_witness(S,Qmax=400):
    best=Fraction(0); bw=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
            v=Fraction(m,q)
            if v>best: best=v; bw=(a,q)
    return best,bw
cons=list(range(1,13))+[182]
M,(a,q)=M_and_witness(cons)
print(f"\nCONSTRUCTION {{1..12,182}}: M={M}={float(M):.5f}, optimal witness t=a/q={a}/{q}")
spread=sorted((s*a)%q for s in cons)
print(f"   speeds*a mod {q} (sorted) = {spread}")
print(f"   folded distances from 0: min={min(min(x,q-x) for x in spread)} (=n={n}) => M=n/{q}=n/Phi_6")
# (3) is 14/183 the covering-min? scan covering sets for anything below
def is_cov(S): return all(any(x%qq==0 for x in S) for qq in range(2,15))
import itertools, random
random.seed(3)
target=Fraction(14,183); below=[]; tested=0; mn=Fraction(1)
# family {1..12, K}
for K in range(13,400):
    S=list(range(1,13))+[K]
    if not is_cov(S): continue
    M,_=M_and_witness(S); tested+=1
    if M<mn: mn=M
    if M<target: below.append((S,M))
# random covering sets
for _ in range(300):
    S=sorted(random.sample(range(1,60),13))
    if not is_cov(S): continue
    M,_=M_and_witness(S); tested+=1
    if M<mn: mn=M
    if M<target: below.append((sorted(S),M))
print(f"\nCOVERING-MIN SCAN: tested {tested} covering sets. min M found = {mn}={float(mn):.5f}")
print(f"   target 14/183={float(target):.5f}. Covering sets BELOW 14/183: {len(below)}")
for S,M in below[:8]:
    print(f"      M={M}={float(M):.5f}  S={S}")
