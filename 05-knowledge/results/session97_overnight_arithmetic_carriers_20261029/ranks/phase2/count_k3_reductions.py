"""Good-reduction point-count obstruction to Q-isomorphism, not geometric identity."""
import json
import sympy as S
from named_family import family
T=S.symbols('T');polys={};disc={}
for i in [273,302,4239]:
 a=list(map(lambda x:S.Poly(x,T),family(i,T)[0]));a1,a2,a3,a4,a6=a;b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 polys[i]=a;disc[i]=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
out=[]
for p in list(S.primerange(5,300)):
 good=True
 for D in disc.values():
  f=S.Poly(D.as_expr(),T,modulus=p)
  if f.degree()!=24 or S.gcd(f,f.diff()).degree()!=0:good=False;break
 if not good:continue
 squares={x*x%p for x in range(1,p)};row={'p':p,'counts':{}}
 for i,a in polys.items():
  total=0
  for t in range(p+1):
   aa=[int(f.eval(t))%p for f in a]if t<p else[int(f.nth(w))%p for f,w in zip(a,[2,4,6,8,12])]
   a1,a2,a3,a4,a6=aa;b2=(a1*a1+4*a2)%p;b4=(2*a4+a1*a3)%p;b6=(a3*a3+4*a6)%p
   trace=0
   for x in range(p):
    f=(4*x**3+b2*x*x+2*b4*x+b6)%p;trace+=0 if f==0 else (1 if f in squares else -1)
   total+=p+1+trace
  row['counts'][str(i)]=total
 out.append(row)
 if len(out)==5:break
print(json.dumps(out,indent=2))
