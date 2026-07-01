"""Coding-theory bounds on the S_n-orbit coloring, n<=6 (exact & fast).
KEY: MFAS(C)=iso-dist(C, transitive class) [transitive orbit = all n! linear orders].
R(n)=max_C MFAS = covering radius; D(n)=max-pair iso-dist = diameter; chain k>=D>=R & k>=info-floor."""
import numpy as np, math
from itertools import combinations, permutations
def setup(n):
    al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
    perms=list(permutations(range(n)))
    dest=np.zeros((len(perms),m),np.int64); flip=np.zeros((len(perms),m),np.int64)
    for pi,p in enumerate(perms):
        for e,(i,j) in enumerate(al):
            a,b=p[i],p[j]
            if a<b: dest[pi,e]=idx[(a,b)]
            else: dest[pi,e]=idx[(b,a)]; flip[pi,e]=1
    return al,m,perms,dest,flip
def orbit(x,m,dest,flip):
    O=set()
    for pi in range(dest.shape[0]):
        y=0
        for e in range(m):
            b=((x>>e)&1)^int(flip[pi,e]); y|=b<<int(dest[pi,e])
        O.add(y)
    return O
def canon(x,m,dest,flip): return min(orbit(x,m,dest,flip))
pc=lambda z:bin(z).count("1")
for n in [3,4,5,6]:
    al,m,perms,dest,flip=setup(n)
    reps={}
    for x in range(2**m):
        c=canon(x,m,dest,flip)
        if c not in reps: reps[c]=c
    classes=list(reps)
    # transitive class = the all-forward tournament x=0 (linear order); its orbit
    trans_orbit=orbit(0,m,dest,flip)
    # |Aut| and MFAS per class
    def aut(x): 
        o=orbit(x,m,dest,flip); return math.factorial(n)//len(o)
    mfas={c:min(pc(c^t) for t in trans_orbit) for c in classes}
    R=max(mfas.values()); maxaut=max(aut(c) for c in classes)
    ext_auts=sorted({aut(c) for c in classes if mfas[c]==R},reverse=True)
    # diameter D(n): max over pairs of min iso-dist (use orbit of one)
    orbits={c:orbit(c,m,dest,flip) for c in classes}
    D=0
    for i,ci in enumerate(classes):
        for cj in classes[i+1:]:
            d=min(pc(ci^b) for b in orbits[cj])
            if d>D: D=d
    floor=math.ceil(math.log2(len(classes)))
    k={3:1,4:2,5:4,6:7}[n]
    print(f"n={n}: #cls={len(classes):2d} floor={floor} R(n)={R} D(n)={D} k(n)={k} max|Aut|={maxaut} | chain k>=D>=R: {k>=D>=R} | R-extremizer|Aut|in{ext_auts}(=max? {ext_auts[0]==maxaut})")
