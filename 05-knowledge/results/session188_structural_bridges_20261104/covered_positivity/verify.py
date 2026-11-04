"""Exact formal coefficient controls for the cubic positivity identity. Stdlib only."""
from collections import defaultdict
from itertools import product
from pathlib import Path
import hashlib,json,random


def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def add(P,indices,c,m):
 if c and min(indices)>0 and max(indices)<=m:P[tuple(sorted(indices))]+=c

def identity_coefficients(m):
 L=defaultdict(int);R=defaultdict(int)
 for k in range(m):
  z=3*(k+1)
  add(L,(k+2,k+1,k+1),z,m);add(L,(k+1,k+1,k+1),2*z,m)
  add(L,(k,k+2,k+2),-z,m);add(L,(k,k+1,k+2),-2*z,m)
 add(R,(1,1,1),1,m)
 for j in range(1,m+2):
  # j(a-c)^2(a+3b+2c)
  for u,v,s in ((j-1,j-1,1),(j-1,j+1,-2),(j+1,j+1,1)):
   for w,c in ((j-1,1),(j,3),(j+1,2)):add(R,(u,v,w),j*s*c,m)
  # (b-c)^2((j+2)b+(2j+1)c)
  for u,v,s in ((j,j,1),(j,j+1,-2),(j+1,j+1,1)):
   for w,c in ((j,j+2),(j+1,2*j+1)):add(R,(u,v,w),s*c,m)
 left={c:v for c,v in L.items() if v};right={c:v for c,v in R.items() if v}
 need(left==right,'formal cubic coefficient identity')
 return {'order':m,'nonzero_monomials':len(left)}

def values(v):
 m=len(v);q=[0]+list(v)+[0,0]
 P=sum((k+1)*(q[k+2]+2*q[k+1])*(q[k+1]**2-q[k]*q[k+2]) for k in range(m))
 RHS=q[1]**3+sum(j*(q[j-1]-q[j+1])**2*(q[j-1]+3*q[j]+2*q[j+1])+(q[j]-q[j+1])**2*((j+2)*q[j]+(2*j+1)*q[j+1]) for j in range(1,m+2))
 need(3*P==RHS,'numeric SOS identity')
 return P

def main():
 formal=[identity_coefficients(m) for m in range(1,21)];positive=zero=0
 for m in range(1,8):
  for v in product(range(5),repeat=m):
   p=values(v)
   if not any(v):need(p==0,'zero profile');zero+=1;continue
   positive+=1;need(p>0,'nonpositive finite profile')
   last=max(i for i,x in enumerate(v) if x)+1
   need(3*p>=v[0]**3+(2*last+3)*v[last-1]**3,'explicit boundary lower bound')
 rng=random.Random(188820)
 for m in range(1,31):
  for _ in range(50):values([rng.randrange(-20,21) for _ in range(m)])
 root=Path(__file__).resolve().parent
 print(json.dumps({'status':'PASS','formal_coefficient_identities':formal,'nonzero_nonnegative_sequences':positive,'zero_sequences':zero,'signed_identity_controls':1500,'boundary_controls':{'single':[values([1]),values([0,0,1])],'gapped':[values([1,0,2]),values([0,3,0,1])]},'source_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
