import json,hashlib
from probe import *
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def substitute(rows,v,module):
 ex=[i for i in range(len(rows)) if i!=v]; a=ex+[v]*len(module);m=len(ex)
 return tuple(sum(((module[i-m]>>(j-m)&1) if u==v and w==v else (rows[u]>>w&1))<<j for j,w in enumerate(a) if i!=j) for i,u in enumerate(a))
def nonempty_covers(rows,K):
 n=len(rows);H=[]
 for S in range(1<<n):
  ex=[i for i in range(n) if S>>i&1]
  H.append(sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(ex)) if ex else 0)
 dp=[0]*(1<<n);dp[0]=1;out=[]
 for k in range(1,K+1):
  nd=[0]*(1<<n)
  for S in range(1<<n):
   A=S
   while A:nd[S]+=H[A]*dp[S^A];A=(A-1)&S
  dp=nd;out.append(dp[-1])
 return out
recs=json.load(open('matches.json'));controls=0
for rec in recs:
 for side in ['a','b']:
  rows=rec[side];v=rec['v'+side];P=polynomials(rows,v)
  for module in [(0,),(2,0),(6,4,0),(2,4,1)]:
   c=nonempty_covers(module,len(P[0]));pred=[sum(x*y for x,y in zip(z,c)) for z in P]
   actual=profile(substitute(rows,v,module))[:8]
   check(pred==list(actual),('module identity',rec['base'],side,module,pred,actual));controls+=1
# Exact arithmetic of adjusted-base differences with re-assigned residue labels.
checks=0
for ja in range(9):
 for jb in range(9):
  for q in range(1,102,2):
   sa,sb=q*(1+3*ja),q*(2+3*jb);ra,rb=sa%3,sb%3
   if not ra or not rb:continue
   diff=-2*q*(ja-jb)+2*((sa-ra)//3-(sb-rb)//3)
   check(3*diff==2*(-q-ra+rb),'gap identity')
   check((diff==0)==(q==1),'maximality');checks+=1
print(json.dumps({'status':'PASS','literal_module_profile_checks':controls,'arithmetic_controls':checks,'marked_matches':len(recs),'matches_sha256':hashlib.sha256(open('matches.json','rb').read()).hexdigest()},indent=2))
