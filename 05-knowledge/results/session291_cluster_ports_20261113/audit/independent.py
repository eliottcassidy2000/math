#!/usr/bin/env python3
"""Independent arithmetic controls; no lane imports."""
from pathlib import Path
import gzip,json,itertools,math
import sympy as s
P=Path(__file__).resolve().parents[1]
R=P.parent/'session290_core_descent_20261113/descent/profiles.json.gz'
def require(x):
 if not x:raise RuntimeError('independent control failed')
r=json.loads(gzip.decompress(R.read_bytes()))['even_records']
b=[1,2,3,4,6,8,9,12]
M=s.Matrix([[dict(t['U']).get(i,0) for t in r] for i in b]);D=s.Matrix.hstack(*[M[:,2*j+1]-M[:,2*j] for j in range(12)])
require(M.rank()==8 and D.rank()==7)
require((s.Matrix([b])*D)==s.zeros(1,12))
g=0
for cols in itertools.combinations(range(12),7):g=math.gcd(g,int(D[1:,:][:,list(cols)].det(method='domain-ge')))
require(g==373248)
x,y=s.symbols('x y');mono=s.Matrix([[1,x,y,x*x,x*y,x**3,y*y,x*x*y]])
pol=list(mono*D)
for f in pol:require(s.rem(f,x-2,x).subs(y,3)==0)
require(s.expand(pol[7])==-485*x*x+336*x*y+531*x-336*y-130)
require(s.expand(pol[5])==-810*x**3+540*x*x*y-150*x*x+926*x*y+389*x-516*y*y-356*y-22)
profiles={tuple(int(z) for z in D*s.Matrix(v)) for v in itertools.product((0,1),repeat=12)}
require(len(profiles)==4096)
for m in (2,3,4,89,90,128):
 count=len({sum(z*a**m for z,a in zip(v,b)) for v in profiles})
 require(count==(4014 if m==2 else 4096))
for m in range(2,401):
 vals=[int(f.subs({x:2**m,y:3**m})) for f in pol]
 require(math.gcd(*vals)==math.gcd(2**m-2,3**m-3))
require(36288*8**90<9**90)
print(json.dumps({'status':'PASS','rank_M':8,'rank_D':7,'independent_minor_gcd':g,'binary_profiles':4096,'gcd_lengths':[2,400],'scalar_lengths':[2,3,4,89,90,128]}))
