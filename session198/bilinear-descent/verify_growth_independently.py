"""Independent matrix reconstruction from the inherited scalar and response laws."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'dependencies'))
import bilinear,c3
basis=[tuple(int(i==j) for i in range(4))for j in range(4)]
I=[[int(i==j)for j in range(4)]for i in range(4)]
def old(x):
 a,f,e,p=x
 return a+2*e+p,e,a+f+2*e+2*p,e+p
def new(x):
 D,E,U,d=x
 return D-E-d,U-D+E-d,E,d-E
def op(k,x):
 D,E,U,d=x
 if k=='O':return D,0,2*U,d
 if k=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D
def trans(A):return list(map(list,zip(*A)))
def mul(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
def ratios(A,B):return min(Fraction(A[i][j],B[i][j])for i,j in product(range(4),repeat=2))
G={k:trans([new(op(k,old(e)))for e in basis])for k in 'OZJ'}
P={k:[[(bilinear.value if k in ('N','M')else c3.value)(k,old(e),old(f))for f in basis]for e in basis]for k in ('N','M','A','B','C','T','J0','J1','J2','J3')}
mins={}
for k,A in P.items():
 for length in (1,2,3):
  q=[]
  for word in product('OZJ',repeat=length):
   W=I
   for letter in word:W=mul(G[letter],W)
   q.extend((ratios(mul(trans(W),A),A),ratios(mul(A,W),A)))
  m=min(q);mins[k+str(length)]=str(m)
  if k in ('N','M') and length==1 and m<Fraction(4,3):raise RuntimeError('C2 one-step')
  if k in ('N','M') and length==2 and m<=2:raise RuntimeError('C2 two-step')
  if k not in ('N','M') and length==3 and m<2:raise RuntimeError('C3 three-step')
out=dict(status='PASS',derived_generators=G,minimum_ratios=mins)
Path(__file__).with_name('INDEPENDENT_GROWTH_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
