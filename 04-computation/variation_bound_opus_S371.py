# opus-2026-07-17-S371 -- THE VARIATION BOUND: the actual sign-keeping bound.
#
# HONEST CORRECTION FIRST: the "leading character sum" C(A) = sum_r prod c(r_i)
# is IDENTICALLY ZERO for odd k -- under r -> -r it flips by (-1)^k, and the
# residue set is antipodally symmetric.  So C is NOT a discriminant between
# families; my framing of it as a leading coefficient was wrong.
#
# BUT THE CONSEQUENCE IS REAL AND USEFUL.  Because C = 0, the MEAN of T drops
# out of  delta = (1/pi^k) sum_r prod c(r_i) T(r) , leaving
#       delta = (1/pi^k) sum_r prod c(r_i) * (T(r) - Tbar).
# So the bound needs only the VARIATION of T across cosets, not its size:
#       |delta| <= (1/pi^k) sum_r |prod c(r_i)| * |T(r) - Tbar|      (VAR-BOUND)
# THE TEST: does VAR-BOUND have the uniform constant that S370's absolute
# bound lacked (2.4 -> 19.0 at k=3)?
from math import sin, pi
import itertools
from collections import defaultdict
LAM=1.0/14
def c(r): return sin(pi*(r%14)/7)
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)

def analyse(A,N):
    k=len(A)
    T=defaultdict(float); coef={}; signed=0.0; fullabs=0.0; m=None
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=k):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0; prod_n=1.0; cc=1.0; pr=1
        for ni in n:
            p*=hhat(ni); prod_n*=ni; cc*=c(ni); pr*=abs(ni)
        signed+=p; fullabs+=abs(p)
        r=tuple(ni%14 for ni in n)
        T[r]+=1.0/prod_n; coef[r]=cc
        if m is None or pr<m: m=pr
    if not T: return None
    Tbar=sum(T.values())/len(T)
    C=sum(coef[r] for r in T)
    var=sum(abs(coef[r])*abs(T[r]-Tbar) for r in T)/pi**k
    absb=sum(abs(coef[r])*abs(T[r]) for r in T)/pi**k
    return signed,fullabs,absb,var,C,m

print("(4) VAR-BOUND vs the S370 ABSOLUTE BOUND   (k=3)")
print("    family          |delta|   VAR-BOUND  COSET-ABS  FULL-ABS   var/|d|  full/|d|")
vr=[]; fr=[]
for A in [(2,3,5),(11,13,17),(31,37,41),(6,10,15),(5,7,11),(101,103,107),
          (1,2,3),(3,4,5),(13,17,23),(9,16,25),(7,11,19),(4,9,13)]:
    R=analyse(A,42)
    if R is None: continue
    sg,fa,ab,var,C,m=R
    vr.append(var/abs(sg)); fr.append(fa/abs(sg))
    print(f"    {str(A):16s} {abs(sg):.6f}  {var:.6f}   {ab:.6f}  {fa:.6f}  {var/abs(sg):7.2f}  {fa/abs(sg):7.2f}")
print()
print(f"    VAR-BOUND  ratio: min {min(vr):.2f}, max {max(vr):.2f}   <-- spread = {max(vr)/min(vr):.1f}x")
print(f"    FULL-ABS   ratio: min {min(fr):.2f}, max {max(fr):.2f}   <-- spread = {max(fr)/min(fr):.1f}x")
print()
print("    A UNIFORM constant (small spread) means the bound is usable in a")
print("    ledger; S370's failure was precisely that its spread grew with k.")
