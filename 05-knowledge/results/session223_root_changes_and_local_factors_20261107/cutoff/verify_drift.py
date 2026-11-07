from fractions import Fraction as Q
from itertools import product

def move(s,c):
 D,E,U,d=s
 return {'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[c]
def run(s,w):
 for c in w:s=move(s,c)
 return s
def F(s):return s[0]-s[1]+s[2]
checks=0; loops=0
for k in range(2,9):
 for p in map(''.join,product('OZ',repeat=k)):
  if len(set(p))<2:continue
  b=2**k
  v=run((0,int(p[-1]=='Z'),0,1),p)
  u=Q(-v[2],b-1)
  delta=run((0,int(p[-1]=='Z'),u,1),p)[0]
  if not delta>0:raise RuntimeError((p,delta))
  loops+=1
  for a in ('','Z','OZJ','ZZJOJ'):
   s=run((1,0,2,1),a)
   for suffix in ('','J','JJ','OZJZO','JOJZZ'):
    vals=[F(run(run(s,p*m),suffix)) for m in range(1,5)]
    coeff=F(run((1,0,0,0),suffix))
    # Linearity means the signed perturbation is legal to evaluate, not physical.
    if coeff<=0:raise RuntimeError(('suffix',suffix,coeff))
    for i in (0,1):
     got=vals[i+2]-(b+1)*vals[i+1]+b*vals[i]
     want=(1-b)*delta*s[3]*coeff
     if got!=want:raise RuntimeError((p,a,suffix,got,want))
     checks+=1
print({'mixed_periods':loops,'exact_second_difference_checks':checks,'status':'PASS'})
