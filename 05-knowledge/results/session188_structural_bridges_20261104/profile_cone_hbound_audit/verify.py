"""Independent exact polarization and canonical-partition cover audit."""
from pathlib import Path
from itertools import product
from fractions import Fraction
from math import factorial,comb
import json,hashlib
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'profile_cone_hbound'
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
 d=json.loads((SOURCE/'receipt.json').read_text());co,A,B,num=poly(5,11,9);need(co==d['coefficient_vector'],'independent eight-way cubic polarization')
 w=35000000;q=[b+w*a for a,b in zip(A,B)];N=factorial(11);f=lambda t:sum(x*t**j for j,x in enumerate(co))
 need(q==d['q'] and f(w)==d['positive_contrast']>0,'profile and positive value');need(q[1]<=N and w%2==0,'Hamilton bound/even lower coefficient')
 need(q[0]==0 and all(x>0 for x in q[1:]) and all(q[k]%factorial(k)==0 for k in range(1,12)),'positivity and divisibility');need(q[11]==N and q[10]==5*N,'top coordinates')
 boundmargin=co[1]-3*abs(co[3])*(N-1)**2;need(boundmargin>0 and co[2]>0,'positive derivative lower bound')
 t=d['first_even_positive_weight_in_this_ray'];need(t==34722116 and t%2==0 and f(t)>0>=f(t-2),'first even positive threshold by monotonicity')
 # Inverse binomial transform, independently in closed form rather than recurrence.
 a=[sum((-1)**(j-k)*comb(10-k,j-k)*q[k+1] for k in range(j+1)) for j in range(11)]
 need(a==d['descent_coefficients'] and a==a[::-1] and min(a)<0 and sum(a)==N,'exact impossible descent distribution')
 violations=[k for k in range(1,12) if q[k]>N*comb(10,k-1)];need(violations==[2,3,4],'stronger cover bounds fail')
 return {'status':'PASS','independent_coefficient_vector':co,'profile':q,'value':f(w),'derivative_positive_lower_bound':boundmargin,'threshold':t,'threshold_adjacent_values':[f(t-2),f(t)],'descent_coefficients':a,'cover_bound_violations':violations,'source_pins':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','receipt.json','receipt.optimized.json']},'scope':'Hamilton-bounded artificial cone obstruction, explicitly not tournament-realizable.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
