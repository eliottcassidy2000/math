#!/usr/bin/env python3
"""Independent full 17-row combined arrangement census via determinant gcds."""
from itertools import product,combinations
from collections import Counter
from math import gcd,prod
from functools import reduce
import json
from check_independent import det,need
D=[(2,0),(0,2),(2,2),(2,-2),(2,1),(2,-1),(1,2),(1,-2)]
rows=list(product([1],range(-1,2),range(-1,2)))+[(0,a,b) for a,b in D]
rowg=[reduce(gcd,r,0) for r in rows]
pairg={(i,j):reduce(gcd,(abs(rows[i][a]*rows[j][b]-rows[i][b]*rows[j][a]) for a,b in combinations(range(3),2)),0) for i,j in combinations(range(17),2)}
triple={(i,j,k):abs(det([rows[i],rows[j],rows[k]])) for i,j,k in combinations(range(17),3)}
snf=Counter();lin=Counter();arith=Counter();maxfactor=0
for mask in range(1<<17):
    ids=[i for i in range(17) if mask>>i&1]
    g1=reduce(gcd,(rowg[i] for i in ids),0)
    g2=0
    for pair in combinations(ids,2):
        g2=gcd(g2,pairg[pair])
        if g2==1:break
    g3=0
    for tripleid in combinations(ids,3):
        g3=gcd(g3,triple[tripleid])
        if g3==1:break
    inv=() if not ids else (g1,) if not g2 else (g1,g2//g1) if not g3 else (g1,g2//g1,g3//g2)
    sign=(-1)**len(ids);snf[inv]+=sign
    lin[3-len(inv)]+=sign;arith[3-len(inv)]+=sign*prod(inv)
    maxfactor=max(maxfactor,max(inv,default=1))
need(dict(lin),{3:1,2:-17,1:79,0:-63},'combined ordinary polynomial')
need(dict(arith),{3:1,2:-21,1:164,0:-504},'combined arithmetic polynomial')
# Direct binomial conversion of -s^3 chi(-(1+s)/s).
from math import comb
P=[0]*4
for degree,c in arith.items():
    for k in range(degree+1):P[3-degree+k]+=(-1)**(degree+1)*c*comb(degree,k)
need(P,[1,24,209,690],'toric Poincare substitution')
need([(-1)**(3-j)*lin[j] for j in range(3,-1,-1)],[1,17,79,63],'additive Poincare substitution')
print(json.dumps({'status':'PASS','subsets':1<<17,'rows':rows,'signed_SNF':{str(k):v for k,v in sorted(snf.items()) if v},'ordinary_coefficients':dict(lin),'arithmetic_coefficients':dict(arith),'largest_Smith_factor':maxfactor,'toric_Betti':P},sort_keys=True,indent=2))
