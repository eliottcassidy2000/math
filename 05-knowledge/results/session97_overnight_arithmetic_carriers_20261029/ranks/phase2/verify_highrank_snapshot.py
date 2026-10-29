"""Independent row verifier, explicit failures survive python -O."""
import json,math,copy
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parent
checks=0
def must(b):
 global checks
 checks+=1
 if not b:raise ValueError('Certificate rejected')
def verify(c,z):
 a1,a2,a3,a4,a6=map(int,c['ainvs']);b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
 pts=[tuple(map(F,p))for p in c['points']];basis={}
 for x,y in pts:must(y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6)
 for row in z['independent_rows']:
  p,r=row['p'],row['root'];must(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)));must(int(c['discriminant'])%p!=0)
  must((4*r**3+b2*r*r+2*b4*r+b6)%p==0)
  bits=0
  for k,(x,y) in enumerate(pts):
   must(x.denominator%p!=0);v=(x.numerator*pow(x.denominator,-1,p)-r)%p;must(v!=0)
   bits|=(pow(v,(p-1)//2,p)==p-1)<<k
  must(bits==row['bits']);w=bits
  for pivot in sorted(basis,reverse=True):
   if w>>pivot&1:w^=basis[pivot]
  must(w!=0);basis[w.bit_length()-1]=w
 must(len(basis)==z['character_rank'])
 p=z['two_torsion_exclusion_prime'];must(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)));must(int(c['discriminant'])%p!=0)
 must(all((4*x**3+b2*x*x+2*b4*x+b6)%p for x in range(p)))
cs={c['id']:c for c in json.loads((P/'highrank_inputs.json').read_text())['curves']}
zs=json.loads((P/'all_highrank_character_certificates.json').read_text())
for z in zs:verify(cs[z['curve_id']],z)
print(json.dumps({'checks':checks,'curves':len(zs),'full_pointsets':sum(z['character_rank']==len(cs[z['curve_id']]['points']) for z in zs),'partial':[{k:z[k]for k in ['curve_id','character_rank','claimed_lower_bound']}for z in zs if z['character_rank']<len(cs[z['curve_id']]['points'])]},sort_keys=True))
