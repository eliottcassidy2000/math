# Multilinear expansion of H(x; p, q, e) in the bits x_i (Moebius inversion over subsets): degree and support.
from ihub import Hint
import itertools
def multilinear(m,p,q,e):
    vals={}
    for mask in range(1<<m):
        x=[(mask>>i)&1 for i in range(m)]
        vals[mask]=Hint(x,p,q,e)
    coef={}
    for S in range(1<<m):
        c=0; sub=S
        while True:
            c+=(-1)**(bin(S).count('1')-bin(sub).count('1'))*vals[sub]
            if sub==0: break
            sub=(sub-1)&S
        if c: coef[S]=c
    return coef
for (m,p,q,e) in [(7,0,0,0),(7,0,0,1),(7,2,3,0),(7,2,3,1),(7,1,1,0),(8,3,2,1)]:
    co=multilinear(m,p,q,e)
    deg=max(bin(S).count('1') for S in co)
    bydeg={}
    for S,c in co.items(): bydeg.setdefault(bin(S).count('1'),[]).append(c)
    print(f"m={m} I=[{p},{p+q}) e={e}: degree {deg}; #terms by degree", {d:len(v) for d,v in sorted(bydeg.items())},
          "; sample deg-2 coefs", sorted(set(abs(c) for c in bydeg.get(2,[])))[:12])
