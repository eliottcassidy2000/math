#!/usr/bin/env python3
"""Exact replay of a finite two-retained-port, one-state digit obstruction."""
import itertools,json,pathlib,random
from fractions import Fraction as Q
from collections import Counter
HERE=pathlib.Path(__file__).resolve().parent

def patterns():
 pairs=list(itertools.product(range(-1,2),repeat=2));out=[]
 for n in range(5):
  for p in itertools.combinations_with_replacement(pairs,n):
   if sum(a==-1 for a,b in p)>2:continue
   if all(sum((a==u)+(b==u)for a,b in p)<=1 or ((u,u)in p and sum(a==u or b==u for a,b in p)==1)for u in range(2)):out.append(p)
 return out
P=patterns();IX={p:i for i,p in enumerate(P)};N=len(P);OPS=list(itertools.product(range(4),range(3),range(2)))

def intro(states,u,A):
 out=Counter()
 for parts,count in states.items():
  left=[-1]+[i for i,(a,b)in enumerate(parts)if b==-1 or b not in A]
  right=[-1]+[i for i,(a,b)in enumerate(parts)if a!=-1 and a in A]
  for i in left:
   for j in right:
    if i==j!=-1:continue
    new=[part for k,part in enumerate(parts)if k!=i and k!=j]
    new.append((parts[i][0]if i!=-1 else u,parts[j][1]if j!=-1 else u))
    out[tuple(sorted(new))]+=count
 return out

def forget(states,u):
 out=Counter()
 for p,c in states.items():
  new=tuple(sorted((-1 if a==u else a,-1 if b==u else b)for a,b in p))
  if sum(a==-1 for a,b in new)<=2:out[new]+=c
 return out

def read(states,r):
 h=s=L=0
 for p,c in states.items():
  if len(p)==1:
   h+=c
   if p[0][0]==r:s+=c
  if len(p)in(1,2)and any(b==r for a,b in p):L+=c
 return h+L-s,s

def step(states,op):
 mask,drop,swap=op;states=forget(intro(states,2,{u for u in range(2)if mask>>u&1}),drop)
 ports=[u for u in range(3)if u!=drop]
 if swap:ports.reverse()
 mp={-1:-1,ports[0]:0,ports[1]:1}
 return Counter({tuple(sorted((mp[a],mp[b])for a,b in p)):c for p,c in states.items()})

def vector(states):return [states.get(p,0)for p in P]
M=[[[ (IX[q],c)for q,c in step({p:1},op).items()]for p in P]for op in OPS]

def row_times(v,m):return [sum(v[i]*c for i,c in col)for col in m]
def col_times(m,v):
 w=[0]*N
 for j,col in enumerate(m):
  for i,c in col:w[i]+=c*v[j]
 return w

def close(states,ci):
 if ci<2:return read(forget(states,1-ci),ci)
 states=intro(states,2,{u for u in range(2)if (ci-2)>>u&1})
 return read(forget(forget(states,0),1),2)
C=[[[close({p:1},ci)[j]for p in P]for j in (0,1)]for ci in range(6)]
BAD=[i for i,p in enumerate(P)if all(a==b and a>=0 for a,b in p)]
BULK=IX[((-1,-1),(0,0),(1,1))]

def add(v,w,a=1,b=1):return [a*x+b*y for x,y in zip(v,w)]
def unit(i):return [int(j==i)for j in range(N)]

def constraints(ci,a,b):
 L=add(C[ci][0],C[ci][1],3,2)
 r0=add(row_times(L,M[a]),L,1,-2)
 delta=add(row_times(L,M[b]),row_times(L,M[a]),1,-1)
 H=[]
 for root in (r0,delta):
  rows=[root]
  for depth in range(1,6):
   rows=[row_times(v,M[o])for v in rows for o in (a,b)]
   H.extend(add(v,root,1,-1)for v in rows)
 H.extend(unit(i)for i in BAD)
 G=add(unit(BULK),delta,6,-1)
 return H,G,delta

def determinant(a):
 a=[r[:]for r in a];sgn=1;prev=1
 for k in range(len(a)-1):
  p=next((i for i in range(k,len(a))if a[i][k]),None)
  if p is None:return 0
  if p!=k:a[k],a[p]=a[p],a[k];sgn=-sgn
  pivot=a[k][k]
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):
    val=a[i][j]*pivot-a[i][k]*a[k][j]
    if val%prev:raise RuntimeError('Bareiss division')
    a[i][j]=val//prev
   a[i][k]=0
  prev=pivot
 return sgn*a[-1][-1]

def eigen2():
 ans=[];dets=[]
 for m in M:
  dense=[[0]*N for _ in range(N)]
  for j,col in enumerate(m):
   for i,c in col:dense[i][j]+=c
  for i in range(N):dense[i][i]-=2
  d=determinant(dense);dets.append(d)
  if not d:ans.append(len(dets)-1)
 return ans,dets

def literal_controls():
 from literal_forest import forest,read as literal_read
 rng=random.Random(20346);checks=0;closing=0
 for n in (3,4,5):
  for trial in range(3):
   rows=[0]*n
   for i in range(n):
    for j in range(i+1,n):
     if rng.randrange(2):rows[i]|=1<<j
     else:rows[j]|=1<<i
   state,_=forest(rows,(0,1));v=vector(state)
   if v[BULK]<1 or any(v[i]for i in BAD):raise RuntimeError('physical constraints')
   for oi,(mask,drop,swap)in enumerate(OPS):
    out=rows[:]+[0]
    for i in range(n):
     if i<2 and mask>>i&1:out[n]|=1<<i
     else:out[i]|=1<<n
    ports=[0,1,n];ports.pop(drop)
    if swap:ports.reverse()
    actual,_=forest(out,ports)
    if col_times(M[oi],v)!=vector(actual):raise RuntimeError(('transition',n,oi))
    checks+=1
   for ci in range(6):
    if ci<2:
     actual,_=forest(rows,(ci,));rr=literal_read(actual,0);expected=(rr[3],rr[1])
    else:
     out=rows[:]+[0];mask=ci-2
     for i in range(n):
      if i<2 and mask>>i&1:out[n]|=1<<i
      else:out[i]|=1<<n
     actual,_=forest(out,(n,));rr=literal_read(actual,0);expected=(rr[3],rr[1])
    got=tuple(sum(c*x for c,x in zip(row,v))for row in C[ci])
    if got!=expected:raise RuntimeError(('closure',n,ci,got,expected))
    closing+=1
 return checks,closing

def main():
 valid,dets=eigen2();cert=json.loads((HERE/'CERTIFICATE.json').read_text());seen=set()
 for case in cert['cases']:
  ci,a,b=case['ci'],case['a'],case['b'];key=ci,a,b
  if key in seen or ci not in range(6)or a not in valid or b not in valid:raise RuntimeError(('invalid case',key))
  seen.add(key);H,G,delta=constraints(ci,a,b);dual={int(i):Q(v)for i,v in case['dual']}
  if any(i<0 or i>len(H)for i in dual):raise RuntimeError('dual index')
  z=dual.get(len(H),Q(0))
  if z<0:raise RuntimeError('negative cone multiplier')
  slack=[sum(v*H[i][j]for i,v in dual.items()if i<len(H))-z*G[j]-delta[j]for j in range(N)]
  if any(v<0 for v in slack):raise RuntimeError(('negative exact slack',key))
 expected=set(itertools.product(range(6),valid,valid))
 if seen!=expected:raise RuntimeError(('incomplete',len(seen),len(expected)))
 tc,cc=literal_controls()
 print(json.dumps(dict(status='PASS',dimension=N,ops=len(M),closures=len(C),eigen2_moves=valid,nonzero_eigen_determinants={str(i):d for i,d in enumerate(dets)if d},exact_Farkas_cases=len(seen),maximum_constraint_word_length=5,impossible_coordinates=BAD,positive_bulk_coordinate=BULK,literal_transition_checks=tc,literal_closure_checks=cc),indent=2))
if __name__=='__main__':main()
