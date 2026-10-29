#!/usr/bin/env python3
"""Exact C3-wreath centralizer ports, LCA invariants and compact transporters."""
from itertools import product,combinations
import random,json
COUNT=0
def need(x,msg):
 global COUNT
 COUNT+=1
 if not x:raise AssertionError(msg)
def digits(v,k):
 z=[0]*k
 for i in range(k-1,-1,-1):z[i]=v%3;v//=3
 return tuple(z)
def integer(z):
 v=0
 for d in z:v=3*v+d
 return v
def neg(z):return tuple((-d)%3 for d in z)
def arc(x,y):
 return next(((b-a)%3==1 for a,b in zip(x,y) if a!=b),False)
def signature(xs):
 out=[]
 for i in range(len(xs)):
  for j in range(i):
   a,b=xs[j],xs[i];l=next((t for t,(x,y) in enumerate(zip(a,b)) if x!=y),len(a))
   out.append((l,int(arc(a,b))))
 return tuple(out)
def canonical(xs):
 if not xs or not xs[0]:return tuple(xs)
 d=xs[0][0];result=[None]*len(xs)
 for c in range(3):
  ids=[i for i,x in enumerate(xs) if x[0]==c]
  sub=canonical([xs[i][1:] for i in ids])
  for i,z in zip(ids,sub):result[i]=((c-d)%3,)+z
 return tuple(result)
def transporter(xs,ys,prefix=(),out=None):
 if out is None:out={}
 if not xs or not xs[0]:return out
 shift=(ys[0][0]-xs[0][0])%3
 need(all((x[0]+shift)%3==y[0] for x,y in zip(xs,ys)),'partial cyclic rotation')
 if shift:out[prefix]=shift
 for c in range(3):
  ids=[i for i,x in enumerate(xs) if x[0]==c]
  transporter([xs[i][1:] for i in ids],[ys[i][1:] for i in ids],prefix+(c,),out)
 return out
def apply(z,portrait):return tuple((x+portrait.get(z[:i],0))%3 for i,x in enumerate(z))
def split(v,k):
 z=digits(v,k)
 a=next((i for i,x in enumerate(z) if x),None)
 if a is None:return (-1,0),()
 j=k-a-1;side=z[a];suffix=z[a+1:]
 return (j,side),suffix if side==1 else neg(suffix)
def join(tag,z,k):
 j,side=tag
 return 0 if j<0 else integer((0,)*(k-j-1)+(side,)+(z if side==1 else neg(z)))
def ckey(vs,k):
 data=[split(v,k) for v in vs];out=[None]*len(vs)
 for j in range(-1,k):
  ids=[i for i,(tag,z) in enumerate(data) if tag[0]==j]
  vals=canonical([data[i][1] for i in ids])
  for i,z in zip(ids,vals):out[i]=(data[i][0],z)
 return tuple(out)
def cinvariant(vs,k):
 data=[split(v,k) for v in vs]
 return tuple(tag for tag,z in data),tuple(signature([z for tag,z in data if tag[0]==j]) for j in range(k))
def ctransport(xs,ys,k):
 dx=[split(x,k) for x in xs];dy=[split(y,k) for y in ys]
 need([t for t,z in dx]==[t for t,z in dy],'immutable tags')
 portraits={}
 for j in range(k):
  ids=[i for i,(t,z) in enumerate(dx) if t[0]==j]
  portraits[j]=transporter([dx[i][1] for i in ids],[dy[i][1] for i in ids])
 return tuple(join(t,apply(z,portraits.get(t[0],{})),k) for t,z in (split(v,k) for v in range(3**k)))
def agroup(k):
 if k==0:return [(0,)]
 child=agroup(k-1);n=3**(k-1)
 return [tuple(((v//n+a)%3)*n+gs[v//n][v%n] for v in range(3*n)) for a in range(3) for gs in product(child,repeat=3)]
def cgroup(k):
 if k==0:return [(0,)]
 C=cgroup(k-1);A=agroup(k-1);n=3**(k-1);theta=[integer(neg(digits(v,k-1))) for v in range(n)]
 return [tuple(c)+tuple(n+a[v] for v in range(n))+tuple(2*n+theta[a[theta[v]]] for v in range(n)) for c in C for a in A]
def main():
 rows=[]
 for k in range(4):
  N=3**k;G=cgroup(k);V=list(range(N));theta=tuple(integer(neg(digits(v,k))) for v in V)
  need(len(G)==3**((N-2*k-1)//4),'centralizer order')
  for g in G:
   need(sorted(g)==V,'permutation');need(all(g[theta[v]]==theta[g[v]] for v in V),'commutes with alpha')
  counts=[]
  for r in range(1,4):
   seen={};keys={}
   for vs in product(V,repeat=r):
    key=ckey(vs,k);inv=cinvariant(vs,k)
    if key not in seen:
     orbit={tuple(g[v] for v in vs) for g in G};seen[key]=orbit
     need(all(ckey(z,k)==key for z in orbit),'canonical invariant')
    need(vs in seen[key],'complete tuple orbit')
    need(inv not in keys or keys[inv]==key,'LCA and arc complete')
    keys[inv]=key
   counts.append(len(seen))
  need(counts[0]==2*k+1,'one port count');need(counts[1]==8*k*k+1,'two port count')
  base=[3**j+3*r for j in range(1,k) for r in range(3**(j-1))]
  expected=0 if k==0 else (3**(k-1)-1)//2
  need(len(base)==expected,'base formula')
  need(sum(all(g[v]==v for v in base) for g in G)==1,'explicit base')
  if k<=3 and expected:
   for subset in combinations(V,expected-1):need(sum(all(g[v]==v for v in subset) for g in G)>1,'no smaller base')
  rows.append(dict(k=k,N=N,centralizer=len(G),tuple_orbits=counts,base=base))
 # Exact mirrored-port hostile.
 need([split(x,2)[0] for x in (3,6)]==[split(x,2)[0] for x in (3,8)],'same marginal tags')
 need(arc(digits(3,2),digits(6,2))==arc(digits(3,2),digits(8,2)),'same original arc')
 need(ckey((3,6),2)!=ckey((3,8),2),'different joint orbit')
 rng=random.Random(93001)
 for k in range(1,6):
  N=3**k
  for trial in range(12):
   xs=tuple(rng.randrange(N) for _ in range(8));portraits={}
   for j in range(k):
    portraits[j]={prefix:rng.randrange(3) for depth in range(j) for prefix in product(range(3),repeat=depth)}
   g=tuple(join(t,apply(z,portraits.get(t[0],{})),k) for t,z in (split(v,k) for v in range(N)))
   ys=tuple(g[x] for x in xs)
   need(cinvariant(xs,k)==cinvariant(ys,k),'large tuple invariant')
   h=ctransport(xs,ys,k)
   need(tuple(h[x] for x in xs)==ys,'compact transporter recovers marks')
   need(sorted(h)==list(range(N)),'expanded permutation')
   for v in range(N):need(h[integer(neg(digits(v,k)))]==integer(neg(digits(h[v],k))),'large centralizer')
   if trial==0:
    for v in range(N):
     for w in range(v):need(arc(digits(v,k),digits(w,k))==arc(digits(h[v],k),digits(h[w],k)),'large tournament automorphism')
 print(json.dumps({'status':'PASS','checks':COUNT,'rows':rows,'large_compact_transporter_max_order':243},indent=2))
if __name__=='__main__':main()
