from itertools import product,combinations
from collections import Counter
import json
choices=[sum(1<<i for i in c) for c in combinations(range(5),2)]
ms=[]
for rows in product(choices,repeat=5):
 if all(sum(r>>j&1 for r in rows)==2 for j in range(5)):ms.append(rows)
def orbit(rows):
 return {tuple(sum(((rows[(i-a)%5]>>((j-b)%5))&1)<<j for j in range(5)) for i in range(5)) for a in range(5) for b in range(5)}
def con(rows):return tuple(sum(((rows[-j%5]>>(-i%5))&1)<<j for j in range(5)) for i in range(5))
seen=set();classes=[]
for r in ms:
 if r in seen:continue
 o=orbit(r);seen|=o;rr=min(o);c=min(orbit(con(rr)));classes.append((rr,len(o),c))
print(json.dumps(dict(raw=len(ms),ordinary=len(classes),self=sum(r==c for r,s,c in classes),merged=len({min(r,c) for r,s,c in classes}),orbits=dict(Counter(s for r,s,c in classes))),indent=2))
