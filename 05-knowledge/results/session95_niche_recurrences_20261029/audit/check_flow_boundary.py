"""Independent literal nonnegative flow audit of C3[T2,T2,T2]."""
from collections import Counter
import json
D=6;checks=0

def check(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def pmul(a,b):
 c=Counter()
 for i,x in a.items():
  for j,y in b.items():
   if i+j<=D:c[i+j]+=x*y
 return c

def part2(t):return [(t-b,b) for b in range(t//2+1)]
def boundary(lam,mu):
 out=Counter()
 for i in range(lam[0]-lam[1]+1):
  for j in range(mu[0]-mu[1]+1):
   power=(mu[0]-j)-(lam[0]-i)
   if 0<=power<=D:out[power]+=1
 return out
pred=Counter();fixpred=Counter()
for t in range(D//3+1):
 ps=part2(t);M={(a,b):boundary(a,b) for a in ps for b in ps}
 for a in ps:
  for j,c in M[a,a].items():
   if 3*t+3*j<=D:fixpred[3*t+3*j]+=c
  for b in ps:
   for c in ps:
    for power,v in pmul(pmul(M[a,b],M[b,c]),M[c,a]).items():
     if power+3*t<=D:pred[power+3*t]+=v
# Vertex numbering (block,internal). Each internal block has0->1.
edges=[(2*b,2*b+1) for b in range(3)]
for b in range(3):
 for i in range(2):
  for j in range(2):edges.append((2*b+i,2*((b+1)%3)+j))
rotation=lambda x: (x+2)%6
rot_edges=[edges.index((rotation(a),rotation(b))) for a,b in edges]
raw=Counter();fixed=Counter();balance=[0]*6;weights=[0]*len(edges)
def visit(k,left,total):
 if k==len(edges):
  if any(balance):return
  raw[total]+=1
  if all(weights[i]==weights[rot_edges[i]] for i in range(len(edges))):fixed[total]+=1
  return
 a,b=edges[k]
 for w in range(left+1):
  weights[k]=w;balance[a]+=w;balance[b]-=w
  visit(k+1,left-w,total+w)
  balance[a]-=w;balance[b]+=w
visit(0,D,0)
for j in range(D+1):
 check(raw[j]==pred[j],('raw',j,raw[j],pred[j]))
 check(fixed[j]==fixpred[j],('fixed',j,fixed[j],fixpred[j]))
 check((raw[j]+2*fixed[j])%3==0,('Burnside',j))
check(raw[3]==8,'scalar port-loss control')
# Same size and scalar circulation series, different full symmetric boundary.
# s_(3)=h3; s_(2,1) has multiplicity1 on exponent type(2,1),2 on(1,1,1).
def weak(total,n):
 if n==1:yield (total,);return
 for j in range(total+1):
  for tail in weak(total-j,n-1):yield (j,)+tail
s3={a:1 for a in weak(3,5)}
s21={a:(2 if max(a)==1 else 1) for a in weak(3,5) if max(a)<=2}
edges5=[(1,0),(2,0),(3,0),(4,0),(2,1),(3,1),(4,1),(3,2),(2,4),(4,3)]
kernel=Counter()
for weights in weak(3,len(edges5)):
 divergence=[0]*5
 for w,(a,b) in zip(weights,edges5):divergence[a]+=w;divergence[b]-=w
 kernel[tuple(divergence)]+=1
def response(A,B):
 return sum(w*ca*B.get(tuple(ai+bi for ai,bi in zip(a,div)),0) for div,w in kernel.items() for a,ca in A.items())
check(response(s3,s21)==878,'same scalar distinct boundary forward')
check(response(s21,s3)==874,'same scalar distinct boundary reverse')
print(json.dumps({'status':'PASS','checks':checks,'degree_cutoff':D,'raw':[raw[j] for j in range(D+1)],'fixed_root_rotation':[fixed[j] for j in range(D+1)],'producer_imports':False},sort_keys=True))
