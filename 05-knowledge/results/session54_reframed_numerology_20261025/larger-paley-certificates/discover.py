#!/usr/bin/env python3
import random,json,pathlib,time
ROOT=pathlib.Path(__file__).resolve().parent
for n in (19,23,31):
 rng=random.Random(540000+n);Q=sorted({x*x%n for x in range(1,n)});E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)};target=len(E)-2*n+1
 def find_cycle():
  budget=[200000]
  def dfs(P,left):
   budget[0]-=1
   if budget[0]<0:return None
   if not left:return P if (P[-1],0) in EI else None
   choices=[v for v in left if (P[-1],v) in EI];rng.shuffle(choices)
   for v in choices:
    z=dfs(P+[v],left-{v})
    if z:return z
  return dfs([0],set(range(1,n)))
 def find_path():
  budget=[200000]
  def dfs(P,left):
   budget[0]-=1
   if budget[0]<0:return None
   if not left:return P+[1] if (P[-1],1) in EI else None
   choices=[v for v in left if (P[-1],v) in EI];rng.shuffle(choices)
   for v in choices:
    z=dfs(P+[v],left-{v})
    if z:return z
  return dfs([0],set(range(2,n)))
 def vector(C):
  row=[0]*len(E)
  for e in zip(C,C[1:]+C[:1]):row[EI[e]]=1
  return row
 base=None;basis={};selected=[];seeds=[];start=time.time();attempts=0
 for trial in range(5):
  C=find_cycle();attempts+=1
  if C is None:continue
  seeds.append(C)
  if base is None:base=C;basevec=vector(C)
  for a in Q:
   for b in range(n):
    cyc=[(a*x+b)%n for x in C];vec=vector(cyc);w=[(x-y)%101 for x,y in zip(vec,basevec)]
    for p,row in sorted(basis.items()):
     f=w[p]
     if f:
      for j in range(p,len(w)):w[j]=(w[j]-f*row[j])%101
    p=next((j for j,x in enumerate(w) if x),None)
    if p is not None:
     inv=pow(w[p],-1,101);basis[p]=[(x*inv)%101 for x in w];selected.append(cyc)
    if len(basis)==target:break
   if len(basis)==target:break
  print(n,'trial',trial+1,'rank',len(basis),'target',target,'elapsed',round(time.time()-start,2),flush=True)
  if len(basis)==target:break
 path=find_path()
 out={'n':n,'modulus':101,'base_cycle':base,'independent_cycles':selected,'pivot_columns':sorted(basis),'seed_cycles':seeds,'endpoint_path':path,'attained_rank':len(basis),'target_rank':target,'seed':540000+n}
 (ROOT/f'certificate-{n}.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print('SAVED',n,'endpoint',bool(path),'full',len(basis)==target,flush=True)
