"""Mechanism: #labeled reps of a class = n!/|Aut|; high-symmetry classes have FEW reps => hardest to reach by a
small fiber. Show at n=6: multiplicity in the minimal spans{1,3} fiber correlates inversely with |Aut|; and the
n=7 obstruction (regular tournament) has the largest |Aut|."""
import numpy as np, math
from itertools import combinations, permutations
def analyze(n, fixed_spans):
    al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
    perms=list(permutations(range(n)))
    def mat(x):
        A=np.zeros((n,n),np.int8)
        for e,(i,j) in enumerate(al):
            if (x>>e)&1: A[j,i]=1
            else: A[i,j]=1
        return A
    def canon_and_aut(A):
        best=None; aut=0
        for p in perms:
            B=A[np.ix_(p,p)]; key=B.tobytes()
            if best is None or key<best: best=key
        # count autos: perms fixing A
        for p in perms:
            if (A[np.ix_(p,p)]==A).all(): aut+=1
        return best,aut
    free=[e for e,(i,j) in enumerate(al) if (j-i) not in fixed_spans]
    fib=[]
    for s in range(2**len(free)):
        x=0
        for b in range(len(free)):
            if (s>>b)&1: x|=1<<free[b]
        fib.append(x)
    from collections import Counter
    mult=Counter(); autof={}
    for x in fib:
        A=mat(x); c,_=canon_and_aut(A); mult[c]+=1
    # aut per class (one rep)
    seen={}
    for x in fib:
        A=mat(x); c,a=canon_and_aut(A)
        if c not in seen: seen[c]=a
    return mult,seen,len(free)
n=6; mult,aut,k=analyze(6,{1,3})
print(f"n=6 minimal fiber (fix spans{{1,3}}, k={k}, 2^{k}={2**k} tilings -> {len(mult)} classes):")
pairs=sorted(((aut[c],mult[c]) for c in mult))
from collections import defaultdict
byaut=defaultdict(list)
for c in mult: byaut[aut[c]].append(mult[c])
for a in sorted(byaut):
    ms=byaut[a]; print(f"   |Aut|={a:2d}: {len(ms)} classes, multiplicity in fiber min={min(ms)} mean={np.mean(ms):.1f} (labeled reps=n!/|Aut|={math.factorial(6)//a})")
print("   => higher |Aut| (fewer labeled reps) => lower multiplicity in the fiber (harder to reach).")
