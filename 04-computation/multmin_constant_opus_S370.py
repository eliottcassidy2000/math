# opus-2026-07-17-S370 -- IS |delta| <= C_k / m7 WITH A UNIFORM CONSTANT?
# The absolute sum CONVERGES (increments fall geometrically, not logarithmically
# -- my divergence heuristic was wrong: at height H a rank-(k-1) lattice has
# ~H^(k-2) vectors each contributing ~H^-k, so sum_H H^-2, convergent for all k).
# So the naive bound IS provable.  The question is the constant.
from math import sin, pi
import itertools, random
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)
def stats(A,N):
    sg=0.0; ab=0.0; m=None
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0; q=1.0; pr=1
        for ni in n:
            p*=hhat(ni); q*=1.0/(pi*abs(ni)); pr*=abs(ni)
        sg+=p; ab+=q
        if m is None or pr<m: m=pr
    return sg,ab,m

print("(3) IS THE CONSTANT UNIFORM?   test  |delta|*m7*pi^3  and  abssum*m7*pi^3")
print("    family           m7    |delta|    |d|*m7*pi^3    abssum*m7*pi^3")
random.seed(370)
fams=[(2,3,5),(3,4,5),(2,5,7),(3,5,8),(11,13,17),(31,37,41),(101,103,107),
      (5,7,11),(6,10,15),(1,2,3),(4,9,13),(7,11,19),(13,17,23),(9,16,25)]
worst_d=0; worst_a=0
for A in fams:
    sg,ab,m=stats(A,50)
    if m is None: continue
    cd=abs(sg)*m*pi**3; ca=ab*m*pi**3
    worst_d=max(worst_d,cd); worst_a=max(worst_a,ca)
    print(f"    {str(A):16s} {m:5d}  {abs(sg):.6f}   {cd:8.3f}      {ca:8.3f}")
print(f"    max over families: |delta| constant {worst_d:.3f}, abs-sum constant {worst_a:.3f}")

print()
print("(4) SAME AT k=4")
w4d=w4a=0
for A in [(1,2,3,5),(2,3,5,7),(11,13,17,19),(3,5,8,13),(5,9,14,23),(31,37,41,43)]:
    sg,ab,m=stats(A,22)
    if m is None: continue
    cd=abs(sg)*m*pi**4; ca=ab*m*pi**4
    w4d=max(w4d,cd); w4a=max(w4a,ca)
    print(f"    {str(A):18s} m7={m:4d}  |delta|={abs(sg):.6f}   |d|*m7*pi^4={cd:7.3f}   abs*m7*pi^4={ca:7.3f}")
print(f"    max: |delta| constant {w4d:.3f}, abs-sum constant {w4a:.3f}")
print()
print("  => a UNIFORM constant would give the provable bound |delta(S)| <= C_k/m7.")
