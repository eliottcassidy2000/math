import importlib.util,json,itertools
from pathlib import Path
p=Path(__file__).resolve().parents[2]/'session251_creative_completion_20261111/ternary/compiler.py'
s=importlib.util.spec_from_file_location('c',p);c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
cores=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
names=['K1','C3','S4','S5a','S5b']
def delete(r,v):
 ids=[i for i in range(len(r)) if i!=v]
 return tuple(sum(((r[i]>>j)&1)<<k for k,j in enumerate(ids)) for i in ids)
for name,r in zip(names,cores):
 print(name,[(v,c.hp_subsets(delete(r,v))[-1]) for v in range(len(r))])
for name,r in zip(names,cores):
 cycles=[]
 for k in (3,5):
  for ids in itertools.permutations(range(len(r)),k):
   if ids[0]!=min(ids):continue
   if all((r[ids[i]]>>ids[(i+1)%k])&1 for i in range(k)):cycles.append(ids)
 print(name,'odd cycles',cycles)
def mat(sig,s):
 d,kp,km,ip,im,op,om,jp,jm=sig
 k,i,o,j=(kp,ip,op,jp) if s==1 else (km,im,om,jm)
 return [[d,o,j],[0,d+k,i],[0,0,d]]
def mul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def join(a,b):
 n=len(a);m=len(b)
 return tuple(r|(((1<<m)-1)<<n) for r in a)+tuple(r<<n for r in b)
checks=0
for a in cores[:2]:
 for b in cores[:2]:
  for oa in range(1<<len(a)):
   for ob in range(1<<len(b)):
    for s in [-1,1]:
     if mul(mat(c.signature(a,oa),s),mat(c.signature(b,ob),s))!=mat(c.signature(join(a,b),oa|(ob<<len(a))),s):raise AssertionError((a,b,oa,ob,s))
     checks+=1
for n in range(4096):
 val=1+2*sum((n//(1<<(i+2)))+bool((n>>(i+1))&1) for i in range(n.bit_length()) if not ((n>>i)&1))
 if val!=c.g(n):raise AssertionError(n)
print('matrix composition checks',checks,'digit-sum checks',4096)
for n in range(1,14):
 r=tuple(sum(1<<j for j in range(n) if j==i+1 or j<i-1) for i in range(n))
 print('stem',n,c.hp_subsets(r)[-1])
