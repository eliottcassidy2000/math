# opus-2026-07-17-S372 -- CORRECTED TEST.  Two fixes to the first run:
#  (a) the k=2 control was VACUOUS -- every residue class shown had a = 0 mod 14,
#      hence 7|a, hence delta = 0 by the S368 criterion.  Use 7-free residues.
#  (b) the normalisation was WRONG.  At k=2 the natural factor is ab because
#      m7 = ab/g^2 IS the multiplicative minimum there.  The k=3 analogue is
#      therefore delta * m7, NOT delta * abc (THM-1080: delta ~ 1/(pi^k m7)).
from math import sin, pi, gcd
import itertools
from collections import defaultdict
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)
def delta_m7(A,N):
    t=0.0; m=None
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0; pr=1
        for ni in n:
            p*=hhat(ni); pr*=abs(ni)
        t+=p
        if m is None or pr<m: m=pr
    return t,m

print("(3) CONTROL, DONE PROPERLY: k=2, 7-free residues, is delta*ab fixed by (a,b) mod 14?")
g2=defaultdict(list)
for a in range(1,60):
    for b in range(a+1,60):
        if gcd(a,b)!=1: continue
        if a%7==0 or b%7==0: continue
        g2[(a%14,b%14)].append((a,b))
shown=0
for r in sorted(g2):
    L=g2[r][:4]
    if len(L)<3: continue
    vals=[delta_m7(A,300)[0]*A[0]*A[1] for A in L]
    rel=(max(vals)-min(vals))/max(abs(v) for v in vals)
    print(f"    residues {str(r):10s} {[('%.5f'%v) for v in vals]}  RELATIVE spread {rel:.2%}")
    shown+=1
    if shown>=4: break
print("    (near-zero spread here = delta*ab IS a function of residues, i.e. THM-965)")

print()
print("(4) THE PROPERLY NORMALISED k=3 TEST: is delta*m7 fixed by (a,b,c) mod 14?")
g3=defaultdict(list)
for a in range(1,26):
    for b in range(a+1,26):
        for c in range(b+1,26):
            if gcd(a,b)!=1 or gcd(b,c)!=1 or gcd(a,c)!=1: continue
            if a%7==0 or b%7==0 or c%7==0: continue
            g3[(a%14,b%14,c%14)].append((a,b,c))
shown=0
for r in sorted(g3):
    L=g3[r][:4]
    if len(L)<3: continue
    out=[]
    for A in L:
        d,m=delta_m7(A,36)
        out.append((A,d,m,d*m))
    vals=[x[3] for x in out]
    rel=(max(vals)-min(vals))/max(abs(v) for v in vals)
    print(f"    residues {str(r):12s}")
    for A,d,m,v in out:
        print(f"        {str(A):14s} delta={d:+.6f}  m7={m:4d}   delta*m7 = {v:+.6f}")
    print(f"        RELATIVE spread {rel:.1%}")
    shown+=1
    if shown>=4: break
