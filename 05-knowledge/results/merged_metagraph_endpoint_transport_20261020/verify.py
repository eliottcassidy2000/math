"""Independent exact controls for endpoint transport and complement leakage.
Run from this directory: python3 -u verify.py
"""
import json,itertools as it,collections as co
from fractions import Fraction as F
from pathlib import Path
from probe import adj,converse,canonical,c3,tiles

def paths(a):
 n=len(a)
 return [p for p in it.permutations(range(n)) if all(a[x]>>y&1 for x,y in zip(p,p[1:]))]
def aut(a):
 n=len(a)
 return sum(all(((a[x]>>y)&1)==((a[p[x]]>>p[y])&1) for x in range(n) for y in range(n)) for p in it.permutations(range(n)))
def blue(a,p):
 n=len(a)
 return all(((a[p[i]]>>p[j])&1)==((a[p[n-1-j]]>>p[n-1-i])&1) for i in range(n) for j in range(i))
def scorekey(s):
 n=len(s);return min(tuple(sorted(s)),tuple(sorted(n-1-v for v in s)))
def predicted_score(s,a,b):
 out=[len(s)-1-v for v in s];out[a]+=1;out[b]-=1;return scorekey(out)
def sccs(a):
 n=len(a);reach=[r|(1<<v) for v,r in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 unseen=set(range(n));comps=[]
 while unseen:
  u=min(unseen);cc={v for v in unseen if reach[u]>>v&1 and reach[v]>>u&1};comps.append(cc);unseen-=cc
 # SCC tournament order: source first.
 comps.sort(key=lambda cc: -sum((a[next(iter(cc))]>>v)&1 for v in range(n) if v not in cc))
 return comps

def regular_ends(a):
 cc=sccs(a)
 return all(len({sum((a[u]>>v)&1 for v in c) for u in c})==1 for c in [cc[0],cc[-1]])
def matrix(d):return [[a+b for a,b in zip(x,y)] for x,y in zip(d['Ablue'],d['Ablack'])]
def rank_mod(M,p=1000000007):
 M=[[int(x)%p for x in row] for row in M];r=0
 for j in range(len(M[0])):
  i=next((i for i in range(r,len(M)) if M[i][j]),None)
  if i is None:continue
  M[r],M[i]=M[i],M[r];z=pow(M[r][j],-1,p);M[r]=[x*z%p for x in M[r]]
  for i in range(len(M)):
   if i!=r and M[i][j]:z=M[i][j];M[i]=[(a-z*b)%p for a,b in zip(M[i],M[r])]
  r+=1
 return r

def audit(d):
 n=d['n'];q=d['q'];D=d['D'];A=matrix(d);h=d['c3'];fail=[];variances=[];rigid=[];endpoint_inventories=[]
 for u,t in enumerate(d['reps']):
  aa=adj(n,t);ss=[r.bit_count() for r in aa];pp=paths(aa);autn=aut(aa);k=1 if d['selfconverse'][u] else 2
  assert len(pp)*k==D[u]*autn
  raw=co.Counter();rawblue=co.Counter();target=co.Counter();targetblue=co.Counter(); ep=co.Counter()
  for p in pp:
   g=ss[p[0]]-ss[p[-1]]-1;raw[g]+=k;target[predicted_score(ss,p[0],p[-1])]+=k;ep[(ss[p[0]],ss[p[-1]])]+=1
   if blue(aa,p):rawblue[g]+=k;targetblue[predicted_score(ss,p[0],p[-1])]+=k
  for expected,M,label in [(raw,A,'all_current'),(rawblue,d['Ablue'],'blue_current')]:
   actual=co.Counter()
   for v,a in enumerate(M[u]):actual[h[v]-h[u]]+=a*autn
   assert +actual==+expected,(n,u,label,actual,expected)
  for expected,M,label in [(target,A,'all_score'),(targetblue,d['Ablue'],'blue_score')]:
   actual=co.Counter()
   for v,a in enumerate(M[u]):actual[scorekey(d['scores'][v])]+=a*autn
   assert +actual==+expected,(n,u,label,actual,expected)
  total=sum(raw.values());mean=F(sum(g*z for g,z in raw.items()),total);var=sum(F(z,autn)*(g-mean)**2 for g,z in raw.items());variances.append(var)
  assert (var==0)==regular_ends(aa),(n,u,var,sccs(aa))
  if var==0:rigid.append(u)
  endpoint_inventories.append({'u':u,'rep':t,'aut':autn,'kappa':k,'HP':len(pp),'endpoint_score_counts':[[a,b,z] for (a,b),z in sorted(ep.items())],'current_counts':[[g,z//autn] for g,z in sorted(raw.items())],'blue_current_counts':[[g,z//autn] for g,z in sorted(rawblue.items())],'variance':str(var),'regular_end_sccs':regular_ends(aa)})
 if n>=4:
  assert not any(A[u][v] for u in rigid for v in rigid)
  for M in [d['Ablue'],d['Ablack']]:assert 2*sum(sum(M[u]) for u in rigid)<=sum(map(sum,M))
 # Rational leakage modulo a good prime supplies a lower-rank certificate.
 import math
 prime=1000000007
 Lmod=[[D[i]*(i==j) for j in range(q)] for i in range(q)]
 for k in range(q):
  nz=[(i,A[i][k]) for i in range(q) if A[i][k]];inv=pow(D[k],-1,prime)
  for i,ai in nz:
   for j,aj in nz:Lmod[i][j]=(Lmod[i][j]-ai*aj*inv)%prime
 # Graph components and bipartite count.
 seen=set();cc=0;bc=0
 for seed in range(q):
  if seed in seen:continue
  cc+=1;side={seed:0};stack=[seed];ok=True;seen.add(seed)
  while stack:
   u=stack.pop()
   for v,a in enumerate(A[u]):
    if not a:continue
    if v not in side:side[v]=1-side[u];seen.add(v);stack.append(v)
    elif side[v]==side[u]:ok=False
  bc+=ok
 rank=rank_mod(Lmod);assert rank==q-cc-bc
 scalar=sum(D[i]*h[i]*h[i]-F(sum(A[i][j]*h[j] for j in range(q))**2,D[i]) for i in range(q));assert scalar==sum(variances)
 energy=sum(h[i]*(D[i]*(i==j)-A[i][j])*h[j] for i in range(q) for j in range(q));assert energy==2**d['m']*F(n-1,4)
 out={'n':n,'q':q,'rigid_triangle_rows':len(rigid),'rigid_half_edge_mass':sum(D[u] for u in rigid),'rigid_blue_half_edge_mass':sum(sum(d['Ablue'][u]) for u in rigid),'rigid_black_half_edge_mass':sum(sum(d['Ablack'][u]) for u in rigid),'nonrigid_triangle_rows':q-len(rigid),'component_count':cc,'bipartite_components':bc,'leakage_rank':rank,'triangle_leakage':str(scalar),'triangle_energy':str(energy),'endpoint_inventories':endpoint_inventories,'failures':0}
 Path(f'audit_n{n}.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='endpoint_inventories'},flush=True)

if __name__=='__main__':
 for n in range(3,8):audit(json.load(open(f'n{n}.json')))
