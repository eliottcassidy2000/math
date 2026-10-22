"""Independent all-position partition audit against the proposed reduced-core engine.
Reference never imports producer partitions, kernels, canonicalization, or propagators.
Engine coefficients are numerical Cauchy extractions, not exact-rational assertions.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import factorial
from collections import Counter,defaultdict
import numpy as np,json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
@lru_cache(None)
def block_partitions(mask):
 if not mask:return ((),)
 first=mask&-mask;rest=mask^first;out=[];sub=rest
 while True:
  block=first|sub
  out.extend((block,)+p for p in block_partitions(rest^sub))
  if sub==0:break
  sub=(sub-1)&rest
 return tuple(out)
def compositions(L):
 if L==0:yield ();return
 for a in range(2,L+1,2):
  for c in compositions(L-a):yield(a,)+c

def reference(S,L,scale):
 n=len(S);indices={};hcache={};answer=defaultdict(Fraction);partition_count=0
 def hom(edges,V):
  sign=1;key=[]
  for a,b in edges:
   if a>b:a,b=b,a;sign=-sign
   key.append((a,b))
  key=(V,tuple(sorted(key)))
  if key not in hcache:
   if V not in indices:indices[V]=np.array(list(product(range(n),repeat=V)),dtype=np.int8).T
   ix=indices[V];val=np.ones(ix.shape[1],dtype=np.int64)
   for a,b in key[1]:val*=S[ix[a],ix[b]]
   hcache[key]=int(val.sum())
  return sign*hcache[key]
 for runs in compositions(L):
  edges=[];pos=0
  for length in runs:
   edges.extend((pos+j,pos+j+1) for j in range(length));pos+=length+1
  for blocks in block_partitions((1<<pos)-1):
   V=len(blocks)
   if V>L+1:continue
   labels=[0]*pos;mu=1
   for j,block in enumerate(blocks):
    size=block.bit_count();mu*=(-1)**(size-1)*factorial(size-1)
    for k in range(pos):
     if block>>k&1:labels[k]=j
   es=[(labels[a],labels[b]) for a,b in edges]
   if any(a==b for a,b in es):continue
   adj=[set()for _ in blocks];degree=[0]*V
   for a,b in es:adj[a].add(b);adj[b].add(a);degree[a]+=1;degree[b]+=1
   reached={0};front=[0]
   while front:
    a=front.pop()
    for b in adj[a]-reached:reached.add(b);front.append(b)
   if len(reached)!=V:continue
   partition_count+=1
   ell=degree.count(1);e=L-V
   if e+ell not in (0,1,2):continue
   answer[e,ell]+=Fraction(mu*hom(es,V),factorial(len(runs))*scale**L)
 return answer,partition_count,len(hcache)

def main():
 sys.path.insert(0,str(HERE.parent/'session30-decorated-engine'))
 import decorated
 controls=[('weighted_skew3',np.array([[0,2,-1],[-2,0,3],[1,-3,0]],dtype=np.int64),7),('transitive4',np.array([[0,1,1,1],[-1,0,1,1],[-1,-1,0,1],[-1,-1,-1,0]],dtype=np.int64),4)]
 out=[];radius=.35;K=32
 for name,S,scale in controls:
  M=S/scale;values=defaultdict(list)
  for k in range(K):
   z=radius*np.exp(2j*np.pi*k/K);at=defaultdict(complex)
   for w in (1,2):
    for d,value in decorated.evaluate(M,w,z,detail=True):at[d['excess'],d['leaves']]+=value
   for key,value in at.items():values[key].append(value)
  for L in (2,4,6):
   t=time.time();ref,pc,hc=reference(S,L,scale);rows=[]
   for key,vals in values.items():
    observed=sum(v*np.exp(-2j*np.pi*k*L/K) for k,v in enumerate(vals))/K/radius**L
    exact=ref[key];err=abs(observed-float(exact));assert err<2e-8,(name,L,key,observed,exact,err)
    rows.append({'excess':key[0],'leaves':key[1],'exact':str(exact),'engine_real':observed.real,'engine_imag':observed.imag,'absolute_error':err})
   row={'control':name,'L':L,'connected_partitions':pc,'hom_cache':hc,'seconds':time.time()-t,'rows':rows};out.append(row)
   print(name,L,'PASS',round(row['seconds'],3),flush=True)
 (HERE/'partition-engine-results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
