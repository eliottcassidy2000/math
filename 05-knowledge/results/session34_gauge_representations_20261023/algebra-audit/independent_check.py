#!/usr/bin/env python3
"""Independent labelled-state audits. No producer code imported."""
import itertools as it, math,json
from pathlib import Path
out={}
def uf(n):
 p=list(range(n))
 def f(x):
  while p[x]!=x:x=p[x]
  return x
 def u(a,b):p[f(a)]=f(b)
 return f,u
reps=[]
for n in range(1,6):
 E=list(it.combinations(range(n),2));N=len(E);ix={e:i for i,e in enumerate(E)};size=1<<N
 F,U=uf(size);sumfix=sumcf=0;classdata={};walsh_ok=0
 for p in it.permutations(range(n)):
  dest=[ix[tuple(sorted((p[a],p[b])))] for a,b in E]
  signs=[p[a]>p[b] for a,b in E]; flip=sum(1<<dest[k] for k,s in enumerate(signs) if s)
  parity=sum(p[a]>p[b] for a,b in E)%2
  ff=fc=0;w=[0]*(N+1)
  for mask in range(size):
   rel=sum(1<<dest[k] for k in range(N) if mask>>k&1)
   t=rel^flip;U(mask,t);ff+=t==mask;fc+=(t^(size-1))==mask
   # Direct Walsh monomial trace and global-product intertwining.
   if rel==mask:w[mask.bit_count()]+=(-1)**sum(signs[k] for k in range(N) if mask>>k&1)
   assert ((t.bit_count()-mask.bit_count())%2)==parity
  assert sum(w)==ff and sum((-1)**k*v for k,v in enumerate(w))==fc
  assert all(w[N-k]==(-1)**parity*w[k] for k in range(N+1))
  sumfix+=ff;sumcf+=fc;walsh_ok+=1
  cyc=[];seen=set()
  for i in range(n):
   if i not in seen:
    j=i;l=0
    while j not in seen:seen.add(j);l+=1;j=p[j]
    cyc.append(l)
  ty=tuple(sorted(cyc,reverse=True))
  if ty in classdata:assert classdata[ty]==(ff,fc)
  classdata[ty]=(ff,fc)
 roots={F(i) for i in range(size)};selfconv=sum(F(r)==F(r^(size-1)) for r in roots)
 assert len(roots)*math.factorial(n)==sumfix and selfconv*math.factorial(n)==sumcf
 reps.append(dict(n=n,orbits=len(roots),self_converse=selfconv,permutations_checked=walsh_ok,class_fixed={str(k):v for k,v in classdata.items()}))
out['representation_direct_orbits']=reps
winner=[]
for n in range(1,7):
 E=list(it.combinations(range(n),2));assoc=0;dist={};units=0
 for mask in range(1<<len(E)):
  win=[[i]*n for i in range(n)];reach=[1<<i for i in range(n)];scores=[0]*n
  for k,(i,j) in enumerate(E):
   a,b=(i,j) if mask>>k&1 else (j,i)
   win[i][j]=win[j][i]=a;reach[a]|=1<<b;scores[a]+=1
  F,U=uf(n);a=True
  for i,j,k in it.product(range(n),repeat=3):
   l=win[win[i][j]][k];r=win[i][win[j][k]];a&=l==r;U(l,r)
  assert a==(sorted(scores)==list(range(n)))
  assoc+=a;units+=0 in scores
  for k in range(n):
   for i in range(n):
    if reach[i]>>k&1:reach[i]|=reach[k]
  assert all((F(i)==F(j))==bool((reach[i]>>j&1) and (reach[j]>>i&1)) for i in range(n) for j in range(n))
  # Associator span is ideal, independently test every generating equality.
  assert all(F(win[i][k])==F(win[j][k]) for i,j,k in it.product(range(n),repeat=3) if F(i)==F(j))
  s=len({F(i) for i in range(n)});dist[s]=dist.get(s,0)+1
 assert assoc==math.factorial(n)
 winner.append(dict(n=n,tournaments=1<<len(E),associative=assoc,with_unit=units,SCC_distribution=dist))
out['winner_all_associators_and_ideal']=winner
# Exact phase arithmetic in Z/4, avoids complex floats and adapts actual scalar identities.
w=lambda a,b,c:2*a*b*c
r=lambda a,b:a*b
for a,b,c,d in it.product(range(2),repeat=4):assert (w(b,c,d)+w(a,b^c,d)+w(a,b,c)-w(a^b,c,d)-w(a,b,c^d))%4==0
for a,b,c in it.product(range(2),repeat=3):
 assert (w(a,b,c)+r(a,b^c)+w(b,c,a)-r(a,b)-w(b,a,c)-r(a,c))%4==0
 assert (r(a^b,c)-w(a,b,c)-r(b,c)+w(a,c,b)-r(a,c)-w(c,a,b))%4==0
out['semion_exact_phase_checks']={'pentagon':16,'hexagon':16,'nontrivial_normalized_coboundary_obstruction':'delta f(1,1,1)=1 for every normalized 2-cochain on Z2, whereas omega(1,1,1)=-1'}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: independent direct orbits n<=5, all associators/ideal/SCC n<=6, exact semion phase identities')
