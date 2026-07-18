# opus-2026-07-17-S372 -- HYP-7550: DOES A THM-965-STYLE CLOSED FORM EXIST AT k=3?
#
# WHY THIS IS THE DECISIVE TEST.  THM-965 is a FORMULA, not a computation:
#     mu(D_a n D_b) = [4ab + fold_M((a+b)%M) - fold_M((b-a)%M)] / (196ab),  M = 14g
# Equivalently  delta * ab  is determined ENTIRELY BY RESIDUES mod 14g.  That
# is exactly what let us bound it uniformly over infinitely many families and
# get the sawtooth floor.  Per-family EVALUATION is already cheap (THM-1065:
# 8191 subset measures in seconds), so evaluation alone buys us NOTHING -- the
# ledger needs a FORMULA IN THE PARAMETERS.
#
# THE TEST: for primitive triples, is  delta * a*b*c  determined by (a,b,c)
# mod 14?  If YES, a THM-965-style formula exists at k=3 and the ledger can be
# built.  If NO, the k=3 term is an irreducible higher-dimensional Dedekind sum
# with no elementary closed form, and the Bonferroni-ledger program is blocked.
from math import sin, pi, gcd
import itertools
from collections import defaultdict
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)
def delta(A,N):
    t=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0
        for ni in n: p*=hhat(ni)
        t+=p
    return t

print("(1) CONTROL: at k=2, is delta*ab determined by residues mod 14?  (g=1)")
groups=defaultdict(list)
for a in range(1,40):
    for b in range(a+1,40):
        if gcd(a,b)!=1: continue
        r=(a%14,b%14)
        groups[r].append((a,b))
for r in sorted(groups):
    L=[x for x in groups[r]][:4]
    if len(L)<3: continue
    vals=[delta(A,260)*A[0]*A[1] for A in L]
    spread=max(vals)-min(vals)
    print(f"    residues {str(r):10s} n={len(L)}  delta*ab = {['%.6f'%v for v in vals]}  spread {spread:.2e}")
    if len([1 for _ in range(1)])>0 and r==(3,5): pass
    if len(groups)>0 and list(sorted(groups)).index(r)>=3: break

print()
print("(2) THE TEST: at k=3, is delta*abc determined by residues mod 14?  (pairwise coprime)")
g3=defaultdict(list)
for a in range(1,30):
    for b in range(a+1,30):
        for c in range(b+1,30):
            if gcd(a,b)!=1 or gcd(b,c)!=1 or gcd(a,c)!=1: continue
            g3[(a%14,b%14,c%14)].append((a,b,c))
shown=0
for r in sorted(g3):
    L=g3[r][:4]
    if len(L)<3: continue
    vals=[delta(A,34)*A[0]*A[1]*A[2] for A in L]
    spread=max(vals)-min(vals); rel=spread/max(abs(v) for v in vals)
    print(f"    residues {str(r):12s} n={len(L)}")
    for A,v in zip(L,vals):
        print(f"        {str(A):14s} delta*abc = {v:12.6f}")
    print(f"        spread {spread:.4f}   RELATIVE {rel:.1%}")
    shown+=1
    if shown>=4: break
