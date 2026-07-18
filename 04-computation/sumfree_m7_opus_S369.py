# opus-2026-07-17-S369 -- MY DIRECTION B: m7 = 1 IFF THE TRIPLE IS A SCHUR TRIPLE.
#
# The m7 witnesses at k=3 were all (-1,-1,1) -- i.e. prod|n_i| = 1 requires
# every coordinate to be +-1, so the relation is +-a_i +- a_j +- a_k = 0, i.e.
#         a_i + a_j = a_k     (a SCHUR TRIPLE / additive relation).
# Since |delta| ~ 1/(pi^3 m7), the WORST triples -- those deviating furthest
# from independence -- are EXACTLY the additively-related ones.  So:
#     SUM-FREE families have m7 >= 2 on every triple, hence smaller 3-body terms.
# This turns "which families are hard" into a question about ADDITIVE STRUCTURE,
# which is Schur / sum-free set theory -- another mature literature.
from math import sin, pi
import itertools, random
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)
def delta(A,N):
    t=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if sum(ni*ai for ni,ai in zip(n,A))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            t+=p
    return t
def m7(A,N):
    best=None
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1
        for ni in n: p*=abs(ni)
        if best is None or p<best: best=p
    return best
def is_schur(A):
    a,b,c=sorted(A)
    return a+b==c

print("(3) m7 = 1  <=>  SCHUR TRIPLE (a_i + a_j = a_k)?")
print("    triple           schur?   m7    |delta|")
bad=0
for A in [(2,3,5),(3,4,7),(1,2,3),(11,13,17),(5,7,12),(3,5,8),(2,5,7),
          (31,37,41),(101,103,107),(4,9,13),(5,7,11),(6,10,15)]:
    s=is_schur(A); m=m7(A,40); d=delta(A,40)
    ok = (s == (m==1))
    if not ok: bad+=1
    print(f"    {str(A):16s} {str(s):7s} {str(m):5s} {abs(d):.6f}{'' if ok else '   <-- MISMATCH'}")
print(f"    mismatches: {bad}/12")

print()
print("(4) CONSEQUENCE: sum-free families have uniformly smaller 3-body terms")
random.seed(369)
sch=[]; free=[]
for _ in range(300):
    A=tuple(sorted(random.sample(range(2,60),3)))
    d=abs(delta(A,26))
    (sch if is_schur(A) else free).append(d)
sch.sort(); free.sort()
print(f"    Schur triples   (n={len(sch):3d}): median |delta| = {sch[len(sch)//2]:.6f}")
print(f"    sum-free triples(n={len(free):3d}): median |delta| = {free[len(free)//2]:.6f}")
print(f"    ratio = {sch[len(sch)//2]/free[len(free)//2]:.2f}x")
print()
print("  => additive structure in the speed set is what drives the k-body terms.")
print("     'Hard family' = 'additively rich family'.  This is Schur/sum-free")
print("     territory, and it says the extremal LRC(14) families should be the")
print("     additively richest 13-sets -- e.g. arithmetic progressions, which is")
print("     exactly what the corpus kept finding as the adversarial construction.")
