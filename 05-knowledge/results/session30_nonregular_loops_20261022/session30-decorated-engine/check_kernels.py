"""Exact direct-subdivision kernel versus resummed descriptor engine.
Independent local half-edge matcher adapted from theory lane core_kernel.py.
The matcher acts on all expanded edges, not the resummed macro-run engine.
"""
import itertools,json,sys
from pathlib import Path
from decorated import decorated_cores
from loop_engine import decompositions_fast
from exact_kernel import kernel

def compositions(total,k):
 if k==1:
  if total>=1:yield(total,)
  return
 for a in range(1,total-k+2):
  for rest in compositions(total-a,k-1):yield(a,)+rest

def coeff(desc,l):
 def power(s,a):return (1+(-1)**a) if s=='*' else s**a
 if desc[0]=='m1':return power(desc[1],l)
 ans=0
 for a in range(1,l):
  for b in range(1,l-a+1):
   c=l-a-b
   z=1 if c==0 else 2*(-1)**(c//2) if c%2==0 else 0
   ans-=power(desc[1],a)*power(desc[2],b)*z
 return ans

def expanded(V,edges,lens):
 out=[];fresh=V
 for (a,b),l in zip(edges,lens):
  prev=a
  for _ in range(l-1):
   out.append((prev,fresh));prev=fresh;fresh+=1
  out.append((prev,b))
 return out

def check(maxL=6):
 cnt=0;by={}
 for w in (1,2):
  for d in decorated_cores(w):
   dec=decompositions_fast(d['V'],d['edges']);c=0
   for L in range(len(d['edges']),maxL+1):
    for lens in compositions(L,len(d['edges'])):
     direct=kernel(expanded(d['V'],d['edges'],lens))
     expected=sum(mu*__import__('math').prod(coeff(p,l) for p,l in zip(props,lens)) for props,mu in dec.items())
     if direct!=expected:raise RuntimeError((d,lens,direct,expected))
     cnt+=1;c+=1
   by[str((w,d['V'],d['edges']))]=c
 return {'status':'PASS','max_total_arcs':maxL,'exact_length_vectors':cnt,'cores':len(by),'per_core':by}
if __name__=='__main__':print(json.dumps(check(int(sys.argv[1]) if len(sys.argv)>1 else 6),indent=2))
