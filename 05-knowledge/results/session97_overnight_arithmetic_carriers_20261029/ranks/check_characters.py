"""Exact local-character witness audit. No floating-point ranks, no rank upper bound."""
import json, math
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def require(ok):
 if not ok: raise ValueError("Exact audit failed")

def inv(a):
 a1,a2,a3,a4,a6=a;b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 return b2,b4,b6,b2*b2-24*b4,-b2**3+36*b2*b4-216*b6,(-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6)

def prime(n):return n>1 and all(n%d for d in range(2,math.isqrt(n)+1))
def check(c):
 a=list(map(Q,c['ainvs']));a1,a2,a3,a4,a6=a;pts=[tuple(map(Q,p))for p in c['points']];b2,b4,b6,c4,c6,disc=inv(a)
 require(disc==int(c['discriminant']));require(1728*disc==c4**3-c6**2)
 require(all(y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6 for x,y in pts))
 basis={};rows=[];noroot=None
 for p in range(3,20000,2):
  if not prime(p) or disc.numerator%p==0 or any(x.denominator%p==0 for x,y in pts):continue
  roots=[r for r in range(p)if (4*r**3+int(b2)*r*r+2*int(b4)*r+int(b6))%p==0]
  if not roots and noroot is None:noroot=p
  for r in roots:
   v=[(x.numerator*pow(x.denominator,-1,p)-r)%p for x,y in pts]
   if 0 in v:continue # special-value extension unnecessary: omit this whole row
   row=sum((pow(t,(p-1)//2,p)==p-1)<<i for i,t in enumerate(v));w=row
   while w:
    k=w.bit_length()-1
    if k in basis:w^=basis[k]
    else:basis[k]=w;rows.append({'p':p,'root':r,'bits':row});break
  if len(basis)==len(pts) and noroot:break
 return {'curve_id':c['id'],'points_on_curve':len(pts),'character_rank':len(basis),'two_torsion_exclusion_prime':noroot,'independent_rows':rows,'claimed_lower_bound':c['rank_lower_bound']}
if __name__=='__main__':
 out=[check(json.loads((P/f'curve{i}.json').read_text()))for i in (273,302,4239)]
 print(json.dumps(out,indent=2))
 require(all(c['character_rank']==c['points_on_curve']for c in out))
