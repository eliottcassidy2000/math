"""Exact independent arithmetic checks; failures survive python -O. SymPy 1.14."""
import copy, hashlib, itertools, json, math
from fractions import Fraction as Q
from functools import reduce
from pathlib import Path
import sympy as s
D=Path(__file__).resolve().parent
checks=0

def need(v):
 global checks
 checks+=1
 if not v: raise ValueError('halving certificate rejected')

def quartic(A,B,C,u):
 return [Q(1),-4*u,-2*B-4*A*u,-8*C-4*B*u,B*B-4*A*C-4*C*u]

def ev(cs,x):
 y=0
 for c in cs:y=y*x+c
 return y

def doubling(A,B,C,P):
 x,y=P
 need(y*y==x**3+A*x*x+B*x+C and y!=0)
 m=(3*x*x+2*A*x+B)/(2*y);u=m*m-A-2*x
 return u,m*(x-u)-y

def verify(c):
 curve=c['curve'];a1,a2,a3,a4,a6=map(Q,curve['ainvs']);A=a2+a1*a1/4;B=a4+a1*a3/2;C=a6+a3*a3/4
 pts=[tuple(map(Q,p))for p in curve['points']];need(len(pts)==30)
 disc=-16*(4*A*A*A*C-A*A*B*B-18*A*B*C+4*B*B*B+27*C*C)
 need(disc==Q(curve['discriminant']))
 for x,y in pts:need((y+(a1*x+a3)/2)**2==x**3+A*x*x+B*x+C)
 rows=c['character_rows'];columns=[0]*30
 for j,row in enumerate(rows):
  p,r=row['p'],row['root'];need(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)))
  red=lambda q:(q.numerator*pow(q.denominator,-1,p))%p
  need(red(disc)!=0); aa,bb,cc=map(red,(A,B,C));need((r**3+aa*r*r+bb*r+cc)%p==0)
  squares={v*v%p for v in range(1,p)};bits=0
  for i,(x,y) in enumerate(pts):
   v=(red(x)-r)%p;need(v!=0);bit=int(v not in squares);bits|=bit<<i;columns[i]|=bit<<j
  need(bits==row['bits'])
 need(len(c['records'])==30)
 z=s.Symbol('z')
 for i,record in enumerate(c['records']):
  need(record['index_zero_based']==i)
  f=s.Poly.from_list(list(map(s.Rational,quartic(A,B,C,pts[i][0]))),z,domain=s.QQ)
  product=s.Poly(s.Rational(record['unit']),z,domain=s.QQ)
  for factor in record['factors']:
   g=s.Poly.from_list(list(map(s.Rational,factor['coefficients'])),z,domain=s.QQ)
   need(g.degree()>1 and g.is_irreducible);product*=g**factor['exponent']
  need(product==f)
  row=record['nonzero_character'];need(row in rows and bool(row['bits']>>i&1))
 counts={}
 for w in range(1,5):
  count=0
  for I in itertools.combinations(range(30),w):
   need(reduce(int.__xor__,(columns[i] for i in I))!=0);count+=1
  counts[w]=count
 return {'id':curve['id'],'no_individual_halves':30,'irreducible_quartics':30,'excluded_support_counts':counts}

def run():
 packet=json.loads((D/'halving_certificate.json').read_text())
 # Source hashes are mandatory when the original minimal source bank is adjacent.
 src=D.parent/'ranks'/'phase2'
 for name,digest in packet['sources'].items():
  if (src/name).exists():need(hashlib.sha256((src/name).read_bytes()).hexdigest()==digest)
 results=[verify(c)for c in packet['curves']]
 # Positive halving control, including the sign forgotten by the quartic.
 A=B=Q(0);C=Q(-2);half=(Q(3),Q(5));target=doubling(A,B,C,half)
 need(target==(Q(129,100),Q(-383,1000)))
 need(ev(quartic(A,B,C,target[0]),half[0])==0)
 need(doubling(A,B,C,(half[0],-half[1]))==(target[0],-target[1]))
 need(doubling(A,B,C,(half[0],-half[1]))!=target)
 rejected=0
 for mode in ('point','character','factor','index'):
  c=copy.deepcopy(packet['curves'][0])
  if mode=='point':c['curve']['points'][0][1]=str(Q(c['curve']['points'][0][1])+1)
  if mode=='character':c['character_rows'][0]['bits']^=1
  if mode=='factor':c['records'][0]['factors'][0]['coefficients'][-1]=str(Q(c['records'][0]['factors'][0]['coefficients'][-1])+1)
  if mode=='index':c['records'][0]['index_zero_based']=1
  try:verify(c)
  except ValueError:rejected+=1
  else:raise ValueError('hostile accepted')
 need(rejected==4)
 print(json.dumps({'curves':results,'hostiles_rejected':rejected,'checks':checks,'scope':'Individual halves and signed distinct-point supports <=4 only; no saturation conclusion'},sort_keys=True))
if __name__=='__main__':run()
