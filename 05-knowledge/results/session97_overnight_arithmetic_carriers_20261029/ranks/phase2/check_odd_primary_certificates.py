"""Finite-reduction E(Fp)->E(Fp)[ell] certificates, avoiding an even-index ledger."""
import json,math
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent

def add(A,B,a,p):
 if A is None:return B
 if B is None:return A
 a1,a2,a3,a4,a6=a;x,y=A;u,v=B
 if x==u:
  if (y+v+a1*x+a3)%p==0:return None
  m=(3*x*x+2*a2*x+a4-a1*y)*pow((2*y+a1*x+a3)%p,-1,p)%p
 else:m=(v-y)*pow((u-x)%p,-1,p)%p
 nu=(y-m*x)%p;z=(m*m+a1*m-a2-x-u)%p;return(z,(-(m+a1)*z-nu-a3)%p)
def mul(k,A,a,p):
 B=None
 while k:
  if k&1:B=add(B,A,a,p)
  A=add(A,A,a,p);k>>=1
 return B
def prime(p):return p>1 and all(p%d for d in range(2,math.isqrt(p)+1))
def check(c,ell):
 a=list(map(int,c['ainvs']));pts=[tuple(map(Q,q))for q in c['points']];basis={};rows=[];exclude=None
 for p in range(5,10000):
  if p==ell or not prime(p)or int(c['discriminant'])%p==0:continue
  a1,a2,a3,a4,a6=a;b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
  squares={x*x%p for x in range(1,p)};N=p+1
  for x in range(p):
   v=(4*x**3+b2*x*x+2*b4*x+b6)%p;N+=0 if v==0 else(1 if v in squares else -1)
  if N%ell:
   if exclude is None:exclude={'p':p,'order':N}
   continue
  if N%(ell*ell)==0 or any(q.denominator%p==0 for xy in pts for q in xy):continue
  image=[mul(N//ell,tuple(q.numerator*pow(q.denominator,-1,p)%p for q in xy),a,p)for xy in pts]
  gen=next((z for z in image if z is not None),None)
  if gen is None:continue
  log={mul(k,gen,a,p):k for k in range(ell)}
  if len(log)!=ell or mul(ell,gen,a,p)is not None:raise ValueError('torsion failure')
  v=[log[z]for z in image];w=v[:]
  for pivot,b in sorted(basis.items(),reverse=True):
   if w[pivot]:
    k=w[pivot];w=[(x-k*y)%ell for x,y in zip(w,b)]
  pivot=next((j for j in range(len(w)-1,-1,-1)if w[j]),None)
  if pivot is not None:
   scale=pow(w[pivot],-1,ell);basis[pivot]=[(z*scale)%ell for z in w];rows.append({'p':p,'order':N,'generator':gen,'row':v})
  if len(basis)==len(pts)and exclude:break
 return {'curve':c['id'],'ell':ell,'rank':len(basis),'torsion_exclusion':exclude,'rows':rows}
if __name__=='__main__':
 cs=json.loads((P/'highrank_inputs.json').read_text())['curves'];out=[check(c,3)for c in cs if c['id']in[3567,3568]];print(json.dumps(out,indent=2))
