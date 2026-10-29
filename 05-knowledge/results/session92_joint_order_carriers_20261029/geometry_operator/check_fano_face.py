#!/usr/bin/env python3
"""Fano nonface controls: exact sector gain propagation, no dense Gram."""
from itertools import combinations
import json
checks=0
def check(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)
lines=[tuple(i-1 for i in T) for T in combinations(range(1,8),3) if T[0]^T[1]^T[2]==0]
check(len(lines)==7,'lines')
f=[0]*8
surviving=[]
for S in range(128):
    coords=[i for i in range(7) if S>>i&1]
    values={0:1};stack=[0];balanced=True
    while stack:
        x=stack.pop()
        for i in coords:
            exponent=sum(all(x>>j&1 for j in F if j!=i) for F in lines if max(F)==i)
            # The actual phase formula has absent coordinates zero.
            g=-((-1)**exponent)
            y=x^(1<<i);proposed=g*values[x]
            if y in values:
                if values[y]!=proposed:balanced=False
            else:values[y]=proposed;stack.append(y)
    check(len(values)==1<<len(coords),'connected sector')
    predicted=not any(all(S>>i&1 for i in F) for F in lines)
    check(balanced==predicted,'gain/face agreement')
    if balanced:f[len(coords)]+=1;surviving.append(S)
check(f==[1,7,21,28,7,0,0,0],'f vector')
check(len(surviving)==64,'nullity')
check(3**7-len(surviving)==2123,'rank')
line_masks={sum(1<<i for i in F) for F in lines}
fours={S for S in surviving if S.bit_count()==4}
check(fours=={127^L for L in line_masks},'line complements')
# Dual-simplex nonzero linear forms on F2^3, evaluated at columns1..7.
simplex={sum((((a&v).bit_count()%2)<< (v-1)) for v in range(1,8)) for a in range(1,8)}
check(fours==simplex,'simplex supports')
# A linear binary subspace containing all seven units would be full F2^7.
check(all(1<<i in surviving for i in range(7)),'singletons')
L=next(iter(line_masks));i=(L&-L);pair=L^i
check(i in surviving and pair in surviving and (i^pair) not in surviving,'nonlinearity')
print(json.dumps({'checks':checks,'sectors':128,'f_vector':f[:5],'nullity':64,'rank':2123},sort_keys=True))
