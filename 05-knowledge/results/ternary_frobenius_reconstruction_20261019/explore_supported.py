import sys,json
from itertools import permutations
import numpy as np
from pathlib import Path
from probe_norm import matrix,basis_mod,rref_mod,norm_monomial

def all2(n):
 a,o=matrix(n)
 js=[j for j,orb in enumerate(o) if min(orb[0])==0]
 return a[:,js],[o[j] for j in js]
def alt(t):
 a,b,c=t
 return {(a,b,c):1,(b,c,a):1,(c,a,b):1,(a,c,b):-1,(c,b,a):-1,(b,a,c):-1}
def candidate(p):
 k=(p-1)//2
 ts=[((1,3,3*p-1),-k),((1,p+1,2*p+1),k),((1,p+2,2*p),-1),((2,p,2*p+1),k),((2,p+1,2*p),k),((3,p,2*p),-k),((p,p+1,p+2),-1)]
 d={}
 for t,c in ts:
  for u,e in alt(t).items(): d[u]=(d.get(u,0)+c*e)%p
 return {t:c for t,c in d.items() if c}
def norm(d,p):
 out={}
 for t,c in d.items():
  for u,e in norm_monomial(t).items(): out[u]=(out.get(u,0)+c*e)%p
 return {u:c for u,c in out.items() if c}
if __name__=='__main__':
 res=[]
 for p in (5,7,11,13,17,19,23,29,31):
  f=candidate(p); nf=norm(f,p)
  print(p,'sparse candidate degree',3*p+3,'terms',len(f),'norm terms',len(nf)); res.append({'p':p,'candidate_degree':3*p+3,'candidate_terms':len(f),'norm_terms':len(nf)})
 for p in (5,7,11,13):
  for n in range(1,38,2):
   a,o=all2(n); dim=len(o)-len(rref_mod(a,p)[1]); k=0; pw=p
   while pw<=n: k+=1; pw*=p
   if dim!=k: raise ValueError((p,n,dim,k))
 print('supported odd dimension matches floor(log_p(n)) for 4 primes, n=1,3,...,37')
 (Path(__file__).resolve().parent / 'supported_validation.json').write_text(json.dumps(res,indent=2)+'\n')
