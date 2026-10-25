from itertools import combinations,product
from collections import Counter
from math import comb
rows=[(1,i,j) for i,j in product((-1,0,1),repeat=2)]
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
triples=[frozenset(s) for s in combinations(range(9),3) if det(*(rows[i] for i in s))==0]
if len(triples)!=8:raise ValueError(triples)
for S in combinations(range(9),4):
 if all(det(*(rows[i] for i in s))==0 for s in combinations(S,3)):raise ValueError('quadruple intersection')
expected={5:{0:28,1:76,2:22},6:{0:2,1:20,2:46,3:16},7:{2:4,3:18,4:12,5:2},8:{4:1,5:4,6:4},9:{8:1}}
for m in range(5,10):
 hist=Counter(sum(t<=set(S) for t in triples) for S in combinations(range(9),m))
 if dict(hist)!=expected[m]:raise ValueError((m,hist))
 d=m-6
 pg=sum((-1)**j*comb(m-3,j)*comb(m-1+d-2*j,m-1) for j in range(d//2+1)) if d>=0 else 0
 K2=d*d*2**(m-3)
 e=2**(m-4)*(m*m-9*m+24)
 if e!=12*(1+pg)-K2:raise ValueError(('Noether',m,pg,K2,e))
 for t in hist:
  double=comb(m,2)-3*t
  outside=3-2*m+comb(m,2)-t
  simple=2*m-2*comb(m,2)+3*t
  singular_e=2**(m-1)*outside+2**(m-2)*simple+2**(m-3)*double+2**(m-4)*t
  nodes=2**(m-4)*t
  if singular_e+nodes!=e:raise ValueError(('stratification',m,t))
 print('m',m,'triple_distribution',dict(sorted(hist.items())),'pg',pg,'K2',K2,'smooth_Euler',e,'PASS')
