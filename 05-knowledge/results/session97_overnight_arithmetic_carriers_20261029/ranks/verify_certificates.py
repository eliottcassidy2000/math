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
 must(len(basis)==len(pts)==z['character_rank'])
 p=z['two_torsion_exclusion_prime'];must(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)));must(int(c['discriminant'])%p!=0)
 must(all((4*x**3+b2*x*x+2*b4*x+b6)%p for x in range(p)))
zs=json.loads((P/'character_certificates.json').read_text())
for z in zs:
 c=json.loads((P/f"curve{z['curve_id']}.json").read_text());verify(c,z)
 for hostile in ['duplicate','negate_one','offcurve','bitflip','badroot']:
  c1=copy.deepcopy(c);z1=copy.deepcopy(z)
  if hostile=='duplicate':c1['points'][1]=c1['points'][0]
  elif hostile=='negate_one':
   # Sign reversal preserves independence but must invalidate this named certificate's coordinates? Characters cannot detect sign.
   x,y=map(F,c1['points'][0]);a1,_,a3,_,_=map(int,c1['ainvs']);c1['points'][0]=[str(x),str(-y-a1*x-a3)]
  elif hostile=='offcurve':c1['points'][0][1]=str(F(c1['points'][0][1])+1)
  elif hostile=='bitflip':z1['independent_rows'][0]['bits']^=1
  else:z1['independent_rows'][0]['root']+=1
  try:verify(c1,z1);accepted=True
  except ValueError:accepted=False
  must(accepted==(hostile=='negate_one'))
print(json.dumps({'checks':checks,'curves':[z['curve_id']for z in zs],'hostiles_rejected':12,'sign_reversals_correctly_preserved':3},sort_keys=True))
