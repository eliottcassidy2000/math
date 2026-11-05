"""Independent raw-coordinate orbit and exhaustive literal small graph audit."""
from itertools import permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def step(s,g):
 D,E,U,d=s
 if g=='O':return D,0,2*U,d
 if g=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

cap=50001;bank={(1,0,2,1)};queue=list(bank);vals={3};proj={(1,1,1)}
for s in queue:
 D,E,U,d=s;vals.add(D+2*U);proj.add(s[:3])
 for g in 'OZJ':
  z=step(s,g)
  if z[0]+2*z[2]<=cap and z not in bank:bank.add(z);queue.append(z)
assert len(bank)+1==908226 and len(proj)==788690 and len(vals)==24958
gaps=[n for n in range(3,cap+1,2)if n not in vals]
expected=[7,11,13,21,23,41,43,47,55,67,77,83,149,163,173,231,233,247,251,277,299,317,391,421,433,443,449,503,519,533,569,583,617,637,821,911,1067,1339,1861,1903,2647,3187]
assert gaps==expected
del bank,queue,proj,vals

def inv(T):
 n=len(T);F=(1<<n)-1;paths={M:[]for M in range(1<<n)}
 for k in range(n+1):
  for p in permutations(range(n),k):
   if all(T[a]>>b&1 for a,b in zip(p,p[1:])):paths[sum(1<<a for a in p)].append(p)
 H=len(paths[F]);s=[sum(p[0]==r for p in paths[F])for r in range(n)]
 E=[sum(p[-1]==r for p in paths[F])for r in range(n)]
 U=[sum(len(paths[F^M])*sum(bool(p)and p[0]==r for p in paths[M])for M in paths)for r in range(n)]
 d=[len(paths[F^(1<<r)])for r in range(n)]
 C=sum(bool(T[p[-1]]&1)for p in paths[F]if p[0]==0)
 return H,s,E,U,d,C

def induced(T,vertices):return [sum(((T[a]>>b)&1)<<j for j,b in enumerate(vertices))for a in vertices]

def append(T,c,g):
 n=len(T);A=T.copy()+[0]
 for v in range(n):
  if g!='O' and v==c:A[n]|=1<<v
  else:A[v]|=1<<n
 return A,n if g=='J'else c

def mapped_equal(T,X,order):
 assert len(T)==len(X)==len(order)
 return all(bool(X[i]>>j&1)==bool(T[order[i]]>>order[j]&1)for i in range(len(X))for j in range(len(X)))

def inverse(T,r):
 n=len(T);targets=[v for v in range(n)if T[r]>>v&1];assert len(targets)==1;c=targets[0]
 F=((1<<n)-1)^(1<<r);zs=[v for v in range(n)if v!=r and not T[v]&F];assert len(zs)==1;z=zs[0]
 assert c!=z;vertices=[v for v in range(n)if v not in[r,z]];W=induced(T,vertices);a=vertices.index(c)
 H,s,E,U,d,C=inv(W);assert s[a]==1
 X,b=append(*append(W,a,'O'),'J');assert mapped_equal(T,X,vertices+[z,r])
 return W,a

sources=balanced=c1=unbalanced=residuals=0
for n in range(1,6):
 edges=[(i,j)for i in range(n)for j in range(i+1,n)]
 for mask in range(1<<len(edges)):
  T=[0]*n
  for k,(i,j)in enumerate(edges):
   if mask>>k&1:T[i]|=1<<j
   else:T[j]|=1<<i
  H,s,E,U,d,C=inv(T)
  for r in range(n):
   if s[r]!=1:continue
   sources+=1
   if n>1 and U[r]==H:balanced+=1;inverse(T,r)
  if C!=1:continue
  c1+=1;roots=[r for r in range(n)if s[r]==1 and U[r]==H]
  if roots:continue
  unbalanced+=1
  rs=[r for r in range(n)if T[r].bit_count()==1];assert len(rs)==1;r=rs[0];c=T[r].bit_length()-1
  kept=set(range(n))-{r};deleted=[]
  while True:
   F=sum(1<<v for v in kept);zs=[v for v in kept if v!=c and T[v]&F==(1<<c)]
   if not zs:break
   assert len(zs)==1;z=zs[0];deleted.append(z);kept.remove(z)
  F=sum(1<<v for v in kept);sinks=[v for v in kept if not T[v]&F];assert len(sinks)==1
  z=sinks[0];kept.remove(z);vertices=sorted(kept);assert len(vertices)>=2 and deleted and c in kept
  W=induced(T,vertices);a=vertices.index(c);DD,ss,EE,UU,dd,CC=inv(W);assert ss[a]==1
  X,b=append(W,a,'O')
  for _ in deleted:X,b=append(X,b,'Z')
  X,b=append(X,b,'J');order=vertices+[z]+list(reversed(deleted))+[r]
  assert mapped_equal(T,X,order)
  k=len(deleted);assert H==DD+2*((1<<(k+1))-1)*UU[a]-((1<<(k+1))-2)*dd[a]
  residuals+=1
out={'status':'PASS','cap':cap,'source_profiles':908226,'projected_source_profiles':788690,'balanced_values':24958,'balanced_gaps':gaps,'literal_source_roots_through_order5':sources,'literal_balanced_roots':balanced,'literal_C1_graphs':c1,'unbalanced_C1_graphs':unbalanced,'canonical_residual_reconstructions':residuals}
(HERE/'BALANCED_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
