from explore import tournament,c3_blowup
from itertools import combinations
import json,hashlib

def check(x,msg):
 if not x:raise RuntimeError(msg)
def homogeneous(a,S):
 return all(len({bool(a[x]>>v&1) for v in S})==1 for x in range(len(a)) if x not in S)
quotients=[(k,c) for k in range(1,5) for c in range(1<<(k*(k-1)//2))]
for k in [5,6]:quotients += [(k,int(c)) for c in sorted(json.load(open('phase2_order%d.json'%k))['class_H'],key=int)]
rows=cards=triple_checks=pair_checks=0
for k,code in quotients:
 T=c3_blowup(tournament(k,code));n=len(T);triangles=[]
 for S in combinations(range(n),3):
  triple_checks+=1
  if homogeneous(T,S) and all(sum(T[v]>>w&1 for w in S if w!=v)==1 for v in S):triangles.append(S)
 check(triangles==[tuple(range(3*b,3*b+3)) for b in range(k)],'original triangle blocks recovered')
 for r in range(n):
  vs=[v for v in range(n) if v!=r];B=tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs);pairs=[]
  for S in combinations(range(n-1),2):
   pair_checks+=1
   if homogeneous(B,S):pairs.append(S)
  expected=tuple(vs.index(v) for v in range(3*(r//3),3*(r//3)+3) if v!=r)
  check(pairs==[expected],'unique surviving homogeneous pair')
  u,w=pairs[0]
  if not(B[u]>>w&1):u,w=w,u
  recovered=(B[u]&~(1<<w))|(1<<u)
  original=sum(1<<j for j,v in enumerate(vs) if T[r]>>v&1)
  check(recovered==original,'canonical labelled completion');cards+=1
 rows+=1
row={'quotient_presentations':rows,'all_labelled_quotient_orders':[1,2,3,4],'strong_class_quotient_orders':[5,6],'deletion_cards':cards,'triple_candidates':triple_checks,'pair_candidates':pair_checks}
print(json.dumps(row,sort_keys=True,separators=(',',':')))
print('PASS semantic_sha256='+hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).hexdigest())
