from math import comb
import numpy as np
import json
from pathlib import Path

def cyclic(t): return [t,(t[1],t[2],t[0]),(t[2],t[0],t[1])]
def orbits(n,interior=False):
 seen=set(); out=[]
 for a in range(n+1):
  for b in range(n-a+1):
   t=(a,b,n-a-b)
   if interior and min(t)==0: continue
   if t in seen: continue
   o=sorted(set(cyclic(t))); seen.update(o); out.append(o)
 return out

def norm_monomial(t):
 a,b,c=t; d={t:1}
 for k in range(b+1):
  u=(k,a+b-k,c); d[u]=d.get(u,0)+comb(b,k)*(-1)**(a+b-k)
 for k in range(a+1):
  u=(k+b,a-k,c); d[u]=d.get(u,0)+comb(a,k)*(-1)**(k+b)
 return {t:c for t,c in d.items() if c}

def matrix(n,interior=False):
 mons=[(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]; idx={t:i for i,t in enumerate(mons)}
 orb=orbits(n,interior); out=np.zeros((len(mons),len(orb)),dtype=np.int64)
 for j,o in enumerate(orb):
  for t in o:
   for u,c in norm_monomial(t).items(): out[idx[u],j]+=c
 return out,orb

def rref_mod(a,p):
 a=a.copy()%p; m,n=a.shape; piv=[]; r=0
 for c in range(n):
  inds=np.flatnonzero(a[r:,c])
  if not len(inds): continue
  j=r+int(inds[0]); a[[r,j]]=a[[j,r]]
  a[r]=a[r]*pow(int(a[r,c]),-1,p)%p
  inds=np.flatnonzero(a[:,c]); inds=inds[inds!=r]
  a[inds]=(a[inds]-a[inds,c,None]*a[r])%p
  piv.append(c); r+=1
  if r==m: break
 return a,piv

def basis_mod(a,p):
 b,piv=rref_mod(a,p); n=a.shape[1]; free=[j for j in range(n) if j not in piv]; out=[]
 for j in free:
  v=np.zeros(n,dtype=np.int64); v[j]=1
  for r,c in enumerate(piv): v[c]=-b[r,j]%p
  out.append(v)
 return out

def poly(v,orb,p): return {','.join(map(str,t)):int(v[j]) for j,o in enumerate(orb) if v[j]%p for t in o}

if __name__=='__main__':
 ps=[2,3,5,7,11,13]; res=[]; first={}
 for n in range(0,37):
  a,o=matrix(n); b,oi=matrix(n,True)
  row={'degree':n}
  for p in ps:
   r=len(rref_mod(a,p)[1]); bi=basis_mod(b,p); row[str(p)]={'W':len(o)-r,'K':len(bi)}
   if bi and p not in first: first[p]={'degree':n,'kernel_dim':len(bi),'witness':poly(bi[0],oi,p)}
  print(row,flush=True); res.append(row)
 (Path(__file__).resolve().parent / 'modular_census.json').write_text(json.dumps({'census':res,'first_interior':first},indent=2)+'\n')
