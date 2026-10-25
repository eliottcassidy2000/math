from itertools import product
from fractions import Fraction
import json
from pathlib import Path
def require(v,msg):
 if not v:raise RuntimeError(msg)
expected={3:[16,64,64,256,256,0],5:[36,88,88,320,324,4],7:[64,160,160,384,448,64]}
out=[]
for p in [3,5,7]:
 r=lambda x:sum(t*t%p==x%p for t in range(p))
 n=plus=minus=full=0;even=odd=Fraction(0)
 for t in product(range(p),repeat=5):
  if not any(t) or next(x for x in t if x)!=1:continue
  c,a,b,d,e=t
  if (a*a+b*b-2*c*c)%p or (d*d+e*e-2*c*c)%p:continue
  hp=r(a*a+d*d-c*c)*r(3*c*c-a*a-d*d)
  hm=r(c*c+a*a-d*d)*r(c*c-a*a+d*d)
  n+=1;plus+=hp;minus+=hm;full+=hp*hm
  even+=Fraction(hp+hm,2)**2;odd+=Fraction(hp-hm,2)**2
 require([n,plus,minus,full,even,odd]==expected[p],'projective direct root census')
 out.append({'p':p,'counts':[n,plus,minus,full],'even_energy':str(even),'odd_energy':str(odd)})
try:require(False,'probe')
except RuntimeError:pass
else:raise RuntimeError('guard disabled')
print(json.dumps(out,sort_keys=True))
