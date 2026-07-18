# opus-2026-07-17-S371 -- THE CHARACTER SUM IS THE SIGN-KEEPING OBJECT.
#
# The residues of lattice vectors form a SUBGROUP Lbar = Lam mod 14 of
# (Z/14)^k.  Writing T(r) = Tbar + (T(r) - Tbar),
#     delta = (1/pi^k) [ Tbar * C(A)  +  sum_r (T(r)-Tbar) prod c(r_i) ],
#       C(A) = sum over r in Lbar of prod_i c(r_i),   c(r) = sin(pi r/7).
# C(A) is a FINITE GAUSS-STYLE CHARACTER SUM over an explicit subgroup -- no
# analysis, no lattice geometry, just a finite computation.  It is the leading
# coefficient, and where it VANISHES the deviation should drop by an order.
# PREDICTION: families with small |delta| relative to m7 should have small C.
from math import sin, pi
import itertools
from collections import defaultdict
LAM=1.0/14
def c(r): return sin(pi*(r%14)/7)
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)

def lbar(A, N):
    """the residue subgroup Lam mod 14, as realised by short lattice vectors."""
    S=set()
    for n in itertools.product(range(-N,N+1), repeat=len(A)):
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        S.add(tuple(ni%14 for ni in n))
    return S

def char_sum(A, N):
    """C(A) = sum over FULL-SUPPORT 7-free residues in Lbar of prod c(r_i)."""
    tot=0.0; cnt=0
    for r in lbar(A,N):
        if any(ri%7==0 for ri in r): continue
        p=1.0
        for ri in r: p*=c(ri)
        tot+=p; cnt+=1
    return tot, cnt

def delta_and_m7(A,N):
    sg=0.0; m=None
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0; pr=1
        for ni in n:
            p*=hhat(ni); pr*=abs(ni)
        sg+=p
        if m is None or pr<m: m=pr
    return sg,m

print("(2) THE CHARACTER SUM C(A) AS LEADING COEFFICIENT  (k=3)")
print("    C(A) = sum over full-support 7-free r in (Lam mod 14) of prod sin(pi r_i/7)")
print("    family           C(A)      #cosets   |delta|     m7    |delta|*m7")
for A in [(2,3,5),(11,13,17),(31,37,41),(6,10,15),(5,7,11),(101,103,107),
          (1,2,3),(3,4,5),(13,17,23),(9,16,25),(7,11,19),(4,9,13)]:
    C,cnt=char_sum(A,30)
    d,m=delta_and_m7(A,40)
    if m is None: continue
    print(f"    {str(A):16s} {C:+9.4f}  {cnt:5d}    {abs(d):.6f} {m:5d}   {abs(d)*m:.6f}")
print()
print("(3) DOES C(A) VANISH WHERE THE DEVIATION IS ANOMALOUSLY SMALL?")
print("    (6,10,15) had |delta| = 0.0004 with m7 = 60 -- the strongest")
print("    cancellation in the S370 table (ratio 11.2).  Check its C:")
for A in [(6,10,15),(2,3,5)]:
    C,cnt=char_sum(A,30)
    print(f"      {str(A):12s} C = {C:+.6f}  over {cnt} full-support cosets")
