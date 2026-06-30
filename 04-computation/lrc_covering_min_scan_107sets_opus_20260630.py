"""Strengthen: broader covering-min scan + is the witness a=n=zeta_6 universal for the worst sets?"""
import math
from fractions import Fraction
from itertools import combinations
n=14; Phi6=183; target=Fraction(14,183)
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
def is_cov(S): return all(any(x%qq==0 for x in S) for qq in range(2,15))
# broader scan: {1..k} + completions; vary base interval and large elements
mn=Fraction(1); below=[]; tested=0; worst=[]
bases=[list(range(1,13)),list(range(1,12)),list(range(1,11)),list(range(2,13)),list(range(1,10))]
comps=[13,14,26,28,39,42,52,56,84,91,156,168,182,169,170,196,364]
for base in bases:
    need=13-len(base)
    for extra in combinations([c for c in comps if c not in base],need):
        S=sorted(set(base)|set(extra))
        if len(S)!=13 or not is_cov(S): continue
        M,w=M_and_witness(S); tested+=1
        if M<mn: mn=M
        if M<target: below.append((S,M))
        if M==target: worst.append((S,w))
print(f"BROADER SCAN: {tested} covering sets. min M = {mn}={float(mn):.5f}. target 14/183={float(target):.5f}")
print(f"   sets BELOW 14/183: {len(below)} (if 0, the construction is the covering-min in this scan)")
for S,M in below[:5]: print(f"      M={float(M):.5f} S={S}")
print(f"\n   sets ACHIEVING 14/183 (the worst/extremal): {len(worst)}; their witnesses a/q:")
wit_n=0
for S,w in worst[:12]:
    a,q=w; is_zeta6 = (q==Phi6 and a==n)
    if is_zeta6: wit_n+=1
    print(f"      a/q={a}/{q} {'= n/Phi6 = ZETA_6!' if is_zeta6 else ''}  S={S}")
print(f"   => of the extremal (M=14/183) sets, witness = n/Phi_6 (zeta_6 rotation) in {wit_n}/{len(worst[:12])} shown")
# the speeds-as-AP under zeta_6 for the construction
cons=list(range(1,13))+[182]
img=sorted((s*n)%Phi6 for s in cons)
print(f"\n   construction speeds * n mod Phi6 = {img}")
print(f"   = AP {{n,2n,...,(n-2)n}} u {{-n=169}}: step n={n}, min dist {min(min(x,Phi6-x) for x in img)} => M=n/Phi6")
print(f"   (n^2=(n-1) mod Phi6 sends the pronic (n-1)n -> (n-1)^2 = -n mod Phi6; the AP closes at -n)")
