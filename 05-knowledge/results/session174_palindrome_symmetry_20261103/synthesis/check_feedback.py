from itertools import combinations,permutations
import json
checks=0
def check(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('feedback audit failure')
def h(A,S=None):
 S=range(len(A)) if S is None else S
 return sum(all(A[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in permutations(S))
def construct(Z):
 m=len(Z)
 return (2,((1<<(m+2))-1)^3)+tuple(1|(row<<2) for row in Z)
rows={}
for m in range(1,5):
 edges=list(combinations(range(m),2));vals=set()
 for mask in range(1<<len(edges)):
  Z=[0]*m
  for k,(i,j) in enumerate(edges):
   if mask>>k&1:Z[i]|=1<<j
   else:Z[j]|=1<<i
  allmask=(1<<m)-1
  H={s:h(Z,[i for i in range(m) if s>>i&1]) for s in range(1<<m)}
  F2=sum(H[s]*H[allmask^s] for s in H)
  check(F2>=2**m)
  literal=h(construct(Z));check(literal==F2+H[allmask]);check(literal not in (7,21));vals.add(literal)
 rows[str(m)]=sorted(vals)
check(rows=={'1':[3],'2':[5],'3':[9,15],'4':[17,27,37]})
print(json.dumps({'status':'PASS','explicit_checks':checks,'all_labelled_feedback_orders':[1,2,3,4],'values':rows,'method':'independent literal permutations and subset two-slot covers, no imports from author'},indent=2))
