"""Independent exact polarization and canonical-partition cover audit."""
from pathlib import Path
from itertools import product
from fractions import Fraction
from math import factorial,comb
import json,hashlib
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'profile_cone_bounded'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def plus(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def times(*ps):
 out={(0,0,0):1}
 for p in ps:
  new={}
  for e,c in out.items():
   for f,d in p.items():
    g=tuple(a+b for a,b in zip(e,f));new[g]=new.get(g,0)+c*d
  out={e:c for e,c in new.items() if c}
 return out
def scale(p,c):return {e:c*v for e,v in p.items()}
ONE={(0,0,0):1};X={(1,0,0):1};Y={(0,1,0):1};Z={(0,0,1):1};P=times(X,Y,Z)
E=plus(X,Y,Z,times(X,Y),times(X,Z),times(Y,Z),scale(P,3))
K=plus(times(X,X,Y,Y),times(X,X,Z,Z),times(Y,Y,Z,Z),times(X,X,Y),times(X,X,Z),times(X,Y,Y),times(X,Z,Z),times(Y,Y,Z),times(Y,Z,Z),scale(times(P,plus(X,Y,Z)),-1),scale(P,-6))
def ttq(n):return [0]+[sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1)) for k in range(1,n+1)]
def trilinear(num,h,qs):
 out=0
 for k in range(min(map(len,qs))):
  for ex,c in num.items():
   t=c*comb(k+h,h)
   for i,q in enumerate(qs):t*=q[k+ex[i]] if k+ex[i]<len(q) else 0
   out+=t
 return Fraction(h*out,3)
def poly(small,big,h):
 A=ttq(small)+[0]*(big-small);B=ttq(big);num=K
 for _ in range(h-1):num=times(num,E)
 out=[Fraction(0)]*4
 for choices in product(range(2),repeat=3):out[sum(choices)]+=trilinear(num,h,[A if b else B for b in choices])
 need(all(x.denominator==1 for x in out),'polarized coefficient not integral');return [int(x) for x in out],A,B,num
def main():
 d=json.loads((SOURCE/'receipt.json').read_text());records=[];certs=0
 for rec in d['records']:
  n=rec['n'];orders=rec['generator_orders'];L=factorial(n)-1;top=ttq(n);vertices=[]
  for r in orders[:-1]:vertices.append([top[k]+L*(ttq(r)[k] if k<=r else 0) for k in range(n+1)])
  vertices.append(top);num=K;coords=[];expected={c['h']:c for c in rec['coordinates']}
  for h in range(1,3*n+3):
   num={e:c for e,c in num.items() if max(e)<=n}
   if h>=2:
    coeff=[]
    for t in __import__('itertools').combinations_with_replacement(range(len(vertices)),3):
     multiplicity=6
     for i in set(t):multiplicity//=factorial(t.count(i))
     value=multiplicity*trilinear(num,h,[vertices[i] for i in t]);need(value.denominator==1 and value<=0,'nonpositive extreme-profile mixed coefficient');coeff.append([list(t),int(value)]);certs+=1
    need(coeff==expected[h]['coefficients'],'complete independent barycentric coefficient table');coords.append({'h':h,'coefficients':coeff})
   num=times(num,E)
  need(not {e:c for e,c in num.items() if max(e)<=n},'vanishing high-coordinate tail')
  records.append({'n':n,'generator_orders':orders,'budget':L,'simplex_vertex_profiles':vertices,'coordinates':coords})
 need(certs==d['integer_coefficient_certificates'],'complete certificate count')
 return {'status':'PASS','method':'Direct trilinear polarization of simplex vertex profiles, independent of lower-weight homogenization expansion.','coordinate_certificates':sum(len(r['coordinates']) for r in records),'integer_coefficients':certs,'records':records,'source_pins':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','receipt.json','receipt.optimized.json']},'scope':'Exact continuous real-weight simplex for n1..7 and every h>=2 only; no all-order or reconstruction claim.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
