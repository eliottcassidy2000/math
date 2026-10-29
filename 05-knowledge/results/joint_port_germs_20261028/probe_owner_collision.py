from itertools import combinations
from fractions import Fraction as F
from functools import reduce
from math import gcd
import check_lrc_germs as c
import json
seen={};found=None;cases=0
for n in range(2,5):
    for V in combinations(range(1,13),n):
        if reduce(gcd,V)!=1:continue
        for d in (F(1,3),F(1,4),F(1,5),F(1,2)):
            L,D,pieces,ends=c.germ(V,d)
            profile=tuple(sorted((b-a,tuple(V[i] for i in c.owners(V,a,d)[0]),
                                  tuple(V[i] for i in c.owners(V,b,d)[1])) for a,b in pieces))
            key=(d,profile);cases+=1
            if not pieces:continue
            if key in seen:
                W=seen[key]
                if c.components(W,[d]*len(W))!=pieces:
                    for w in range(1,25):
                        if w in set(V)|set(W):continue
                        aa=c.components(W+(w,),[d]*(len(W)+1));bb=c.components(V+(w,),[d]*(len(V)+1))
                        if bool(aa)!=bool(bb):
                            found=dict(V=W,W=V,delta=str(d),added_speed=w,
                                       first=[[str(a),str(b)] for a,b in c.components(W,[d]*len(W))],
                                       second=[[str(a),str(b)] for a,b in pieces],
                                       future_nonempty=[bool(aa),bool(bb)],
                                       profile=[[str(length),P,N] for length,P,N in profile])
                            break
            seen.setdefault(key,V)
            if found:break
        if found:break
    if found:break
print(json.dumps({'tested':cases,'witness':found},indent=2))
