import sys,itertools,json,functools
from pathlib import Path
from fractions import Fraction
P=Path(__file__).resolve().parent.parent;sys.path.insert(0,str(P))
import compiler,counting
@functools.lru_cache(None)
def direct(a):
 n=len(a);full=(1<<n)-1;H=[0]*(1<<n);E=[0]*(1<<n);S=[0]*(1<<n);H[0]=1
 for m in range(1,full+1):
  vs=[i for i in range(n) if m>>i&1]
  for p in itertools.permutations(vs):
   if all(a[x]>>y&1 for x,y in zip(p,p[1:])):
    H[m]+=1;E[m]+=p[-1]==0;S[m]+=p[0]==0
 C=[[int(m==0) for m in range(full+1)],H]
 for k in range(2,4):
  C.append([sum(C[-1][j]*H[m^j] for j in range(m+1) if j&m==j) for m in range(full+1)])
 return tuple((C[k][full],sum(E[m]*C[k-1][full^m] for m in range(full+1)),sum(S[m]*C[k][full^m] for m in range(full+1)),C[k][full^1]) for k in range(1,4)),S[full]
def nine(a):
 x,s=direct(a);h,e,R,d=x[0];C,L,T,b=x[1];return h,s,e,L,R,d,C,b,T
checks=0
for n in range(1,5):
 for a in counting.labelled(n):
  p=nine(a)
  for g in (1,2):
   actual=nine(counting.gate(a,g));pred=compiler.move(p,'K'+str(g))
   if n==1 and g==1:pred=pred[:1]+(1,)+pred[2:]
   if actual!=pred:raise RuntimeError(('gate mismatch',a,g,actual,pred))
   checks+=1
# Verify arbitrary one-port contexts, including all 3-vertex contexts, for k=1,2,3.
def join(a,g):
 n=len(a);m=len(g);b=list(a)+[0]*(m-1)
 for j in range(1,m):
  for i in range(1,n):b[i]|=1<<(n+j-1)
 for i in range(m):
  for j in range(m):
   if g[i]>>j&1:b[0 if i==0 else n+i-1]|=1<<(0 if j==0 else n+j-1)
 return counting.relabel(tuple(b),[n]+list(range(n))+list(range(n+1,n+m-1)))
def transfer(X,Y):
 z=[[Fraction(v) for v in x+y] for x,y in zip(X,Y)]
 for i in range(4):
  j=next(j for j in range(i,4) if z[j][i]);z[i],z[j]=z[j],z[i];d=z[i][i];z[i]=[v/d for v in z[i]]
  for j in range(4):
   if j!=i:
    d=z[j][i];z[j]=[v-d*w for v,w in zip(z[j],z[i])]
 return [row[4:] for row in z]
basis=[(0,),(2,0),(4,1,2),(6,0,2)];context_checks=0
for g in counting.labelled(3):
 for k in range(3):
  X=[list(direct(a)[0][k]) for a in basis];Y=[list(direct(join(a,g))[0][k]) for a in basis];M=transfer(X,Y)
  for a in counting.labelled(3):
   x=direct(a)[0][k];actual=direct(join(a,g))[0][k];pred=tuple(sum(x[j]*M[j][i] for j in range(4)) for i in range(4))
   if actual!=pred:raise RuntimeError(('context mismatch',a,g,k,actual,pred))
   context_checks+=1
# An unmerged word tree supplies an independent control for shortest-state merging.
bound=139;cap=bound.bit_length()+2;stack=[(2,compiler.SEED,())];mins={};support=set();words=0
while stack:
 n,p,w=stack.pop();words+=1
 if p not in mins or mins[p]>n:mins[p]=n
 if p[1]==2:support.add(compiler.base(p))
 for op in ('O','Z','J','K1','K2') if p[1]==1 else ('O','Z','J'):
  nn=n+compiler.COST[op]
  if nn<=cap:
   q=compiler.move(p,op)
   if compiler.base(q)<=bound:stack.append((nn,q,w+(op,)))
info,ws=compiler.inventory(bound)
if not info['complete'] or set(ws)!=support or info['states']!=len(mins):raise RuntimeError(('merge mismatch',info,len(mins),support,set(ws)))
print(json.dumps(dict(status='PASS',independent_permutation_gate_checks=checks,arbitrary_one_port_context_checks=context_checks,unmerged_words=words,merged_states=len(mins),bound=bound,support=sorted(support)),indent=2))
