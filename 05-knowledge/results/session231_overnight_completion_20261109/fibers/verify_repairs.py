import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import states,compose
from verify import literal,need
from itertools import permutations
import json
one=[(1,2,b,b,0,b,b) for b in (0,1)]
def closure(parts):
 a=[];mask=0;st=None
 for adj,cut in parts:
  old=len(a);n=len(adj)
  for i in range(old):a[i]|=((1<<n)-1)<<old
  a += [v<<old for v in adj];mask |=cut<<old
  z=states(adj)[cut];st=compose(st,z) if st is not None else z
 N=len(a);a += [1<<(N+1),mask]
 for i in range(N):
  a[i]|=1<<N
  if not(mask>>i&1):a[i]|=1<<(N+1)
 # Forward Hamilton subset DP, independent of composition formula.
 dp=[[0]*len(a) for _ in range(1<<len(a))]
 for i in range(len(a)):dp[1<<i][i]=1
 for s in range(1,1<<len(a)):
  for v in range(len(a)):
   if not dp[s][v]:continue
   opts=a[v]&~s
   while opts:
    bit=opts&-opts;opts-=bit;w=bit.bit_length()-1;dp[s|bit][w]+=dp[s][v]
 H=sum(dp[-1]);need(H==2*st[2]+st[0]-st[4],'closure formula')
 # Strong connectivity independently.
 for v in range(len(a)):
  seen=1<<v
  while True:
   old=seen
   for i in range(len(a)):
    if seen>>i&1:seen|=a[i]
   if seen==old:break
  need(seen==(1<<len(a))-1,'not strong')
 return dict(vertices=len(a),state=st,H=H)
p=lambda b:([0],b)
C=[28,41,58,48,34,1];D=[20,41,58,49,34,1];B=[26,44,57,16,34,9]
for X in [C,D]:need(states(X)==literal(X),'literal H29')
for X,Y in [(B,C),(B,D),(C,D)]:
 need(not any(all(bool(X[i]>>j&1)==bool(Y[perm[i]]>>perm[j]&1) for i in range(6) for j in range(i+1,6)) for perm in permutations(range(6))),'isomorphic')
out={'H29':[closure([p(1),(X,19),p(0)]) for X in [C,D]]}
C45=[48,49,3,7,12,28];D45=[40,37,17,6,11,28]
out['H45']=[closure([p(1),(X,m),p(0)]) for X,m in [(C45,9),(C45,3),(D45,3),(C45,6),(D45,10)]]
H5=[8,1,3,6]
out['H5']=[closure([p(1),(H5,4),p(1),p(1),p(0),p(0)]),closure([p(1),p(0),p(1),p(1),(H5,2),p(0)])]
H9=[16,1,3,7,14]
out['H9']=[closure([p(1),(H9,13),p(0),p(1),p(0)]),closure([p(1),p(1),(H9,4),p(1),p(0)]),closure([p(1),p(1),p(0),(H9,22),p(0)])]
print(json.dumps(out,sort_keys=True,indent=2))
