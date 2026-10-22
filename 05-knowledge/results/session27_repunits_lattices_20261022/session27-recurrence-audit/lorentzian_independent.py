from itertools import permutations,combinations
from math import comb
import json
rows=[]
for masks,H in [([150,172,216,49,98,69,11,120],657),([156,165,216,50,98,69,11,120],643)]:
 n=8
 assert all(((masks[a]>>b)&1)+((masks[b]>>a)&1)==1 for a,b in combinations(range(n),2))
 count=sum(all((masks[p[i]]>>p[i+1])&1 for i in range(7)) for p in permutations(range(8)))
 assert count==H
 padded=masks+[sum(1<<j for j in list(range(8))+list(range(i+1,25))) for i in range(8,25)]
 scores=[m.bit_count() for m in padded]
 assert scores==[4,4,4,3,3,3,3,4]+list(range(24,7,-1))
 c3=0
 for triple in combinations(range(25),3):
  c3+=all(sum((padded[a]>>b)&1 for b in triple if b!=a)==1 for a in triple)
 assert c3==20
 assert sum(scores)==300 and sum(a*a for a in scores)-4900==-40
 rows.append(dict(H8=count,c3_25=c3,norm25=-40,scores25=scores))
w=list(range(25))+[70]
r=[1,-1]+[0]*24
def pair(a,b):return sum(x*y for x,y in zip(a[:25],b[:25]))-a[25]*b[25]
assert pair(w,w)==0 and pair(r,w)==-1 and pair(r,r)==2
# Explicit vectors in w-perp and their representatives; y=2*z+2*(z,w)*r
# guarantees even integral coordinate sum and zero pairing with w.
for i in range(26):
 z=[int(j==i) for j in range(26)]
 y=[2*a+2*pair(z,w)*b for a,b in zip(z,r)]
 assert pair(y,w)==0 and (sum(y[:25])-y[25])%2==0
 a=(pair(y,y)+2*pair(r,y))//2
 root=[b+c+a*d for b,c,d in zip(r,y,w)]
 assert pair(root,root)==2 and pair(root,w)==-1
 assert pair([b-c for b,c in zip(root,r)],[b-c for b,c in zip(root,r)])==pair(y,y)
 # Changing representative by w must produce the same root.
 yp=[b+c for b,c in zip(y,w)]
 ap=(pair(yp,yp)+2*pair(r,yp))//2
 assert ap==a-1 and [b+c+ap*d for b,c,d in zip(r,yp,w)]==root
print(json.dumps(dict(status='PASS',cases=rows,quotient_representative_controls=26,hamiltonian_25_enumerated=False),indent=2))
