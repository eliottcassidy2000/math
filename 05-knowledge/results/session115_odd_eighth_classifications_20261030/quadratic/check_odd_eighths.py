#!/usr/bin/env python3
from itertools import product
from collections import Counter
from math import isqrt
import json
checks=0
def gate(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(str(msg))
def dot(x,y):return (x&y).bit_count()%2
def walsh(f):return [sum(1 if f[x]^dot(u,x)==0 else -1 for x in range(len(f))) for u in range(len(f))]
def normalize(f):
 return tuple(f[x]^sum(f[1<<j]*((x>>j)&1) for j in range(3))%2 for x in range(8))
def phase(A,B,C,D):return (0,0,0,A,0,B,D,C^D)
def bits(f):return (f[3],f[5],f[6]^f[7],f[6])
def make_table(z):
 F=[[0]*8 for _ in range(8)]
 for a in range(1,8):F[a][a]=1
 lines=[(1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7),(3,5,6)]
 for (a,b,c),t in zip(lines,z):
  for x,y in [(a,b),(b,c),(c,a)]:F[x][y]=t;F[y][x]=t^1
 return F
base=make_table((0,0,0,0,1,1,1));rows=[];profiles=Counter()
for A,B,C,D in product(range(2),repeat=4):
 s=phase(A,B,C,D);W=walsh(s);kap=A^B^C
 gate(sum(s)%2==kap,('cubic',A,B,C,D))
 if kap:
  v=(1^D)+2*(1^B)+4*(1^A);eps=int(v==0)
  r=sum((int(v==1<<i)^eps)<<i for i in range(3))
  gate(all(s[x]==(eps^dot(r,x)^int(x==v)) for x in range(8)),('defect',s,v))
  gate(sorted(map(abs,W))==[2]*7+[6],('cubic spectrum',s,W))
  want=Counter({6:1,2:4,-2:3}) if v else Counter({-6:1,2:7})
  gate(Counter(W)==want,('signed profile',s,W))
 elif any(s):gate(Counter(W)==Counter({0:4,4:3,-4:1}),('quadratic',s,W))
 else:gate(Counter(W)==Counter({8:1,0:7}),('zero',s,W))
 z=(A,B,C,D,1^B^C^D,1^A^C^D,1^A^B^D);F=make_table(z)
 for x,y in product(range(8),repeat=2):gate(F[x][y]==base[x][y]^s[x]^s[y]^s[x^y],('coboundary',z,x,y))
 for t in range(1,8):
  for a,b in product(range(8),repeat=2):gate(F[a][b]^F[a^t][b^t]^F[a][b^t]^F[a^t][b]==1,('rectangle',z,t,a,b))
 m=1;sig=1
 for p,bit in zip((19,2,59,29,5,89,7),z):
  if bit:m*=p;sig*=p+1
 gate(sig%30==0 and isqrt(sig//30)**2==sig//30,('norm arithmetic',m,sig))
 gate(kap==(z[0]^z[1]^z[2]),('selection tag',m))
 profile=tuple(sorted(Counter(W).items()));profiles[str(profile)]+=1
 rows.append({'ABCD':[A,B,C,D],'kappa':kap,'walsh':W,'arithmetic_source':m})
# Every s(0)=0 gauge has the same profile as its canonical representative.
for n in range(128):
 f=(0,)+tuple((n>>i)&1 for i in range(7));nf=normalize(f)
 gate(sorted(walsh(f))==sorted(walsh(nf)),('normalization',n))
 gate(nf==phase(*bits(nf)),('coordinates',n))
# The entire Boolean cube at dimension3, not just normalized gauges.
odd=0;defects=0
for n in range(256):
 f=tuple((n>>x)&1 for x in range(8));W=walsh(f)
 odd+=sum(f)%2;isdef=6 in map(abs,W);defects+=isdef
 gate(isdef==bool(sum(f)%2),('complete odd sector',n))
# Literal GL3 action on normalized relative-gauge classes.
GL=[]
for a,b,c in product(range(1,8),repeat=3):
 imgs=[((a if x&1 else 0)^(b if x&2 else 0)^(c if x&4 else 0)) for x in range(8)]
 if len(set(imgs))==8:GL.append(imgs)
gate(len(GL)==168,'GL order')
orbits=set()
for key in product(range(2),repeat=4):
 f=phase(*key);orbit=set()
 for g in GL:
  transformed=normalize(tuple(f[g[x]] for x in range(8)));orbit.add(bits(transformed))
  gate(sum(transformed)%2==sum(f)%2,('degree sector invariant',key))
 orbits.add(tuple(sorted(orbit)))
gate(sorted(map(len,orbits))==[1,1,7,7],('orbit classification',orbits))
# Uniform single-defect formula: all parameters at d3, selected but unbounded-shape controls d4..7.
uniform=0
for d in range(3,8):
 N=1<<d
 vs=range(N) if d==3 else sorted({0,1,2,N//2,N-1,N//3})
 rs=range(N) if d==3 else sorted({0,1,N//2,N-1})
 for v,r,eps in product(vs,rs,range(2)):
  f=tuple(eps^dot(r,x)^int(x==v) for x in range(N));W=walsh(f)
  want=[(-1)**eps*(N*int(u==r)-2*((-1)**dot(u^r,v))) for u in range(N)]
  gate(W==want,('uniform',d,v,r,eps));uniform+=1
# Odd weight alone no longer identifies a one-defect affine function.
f=tuple(int(x in (0,1,2)) for x in range(16));W=walsh(f)
gate(sum(f)==3 and max(map(abs,W))==10,'dimension4 hostile')
# At d2, all four Walsh magnitudes can equal N-2, destroying unique exceptional frequency.
f=(0,0,0,1);gate(sorted(map(abs,walsh(f)))==[2]*4,'dimension2 boundary')
print(json.dumps({'status':'PASS','checks':checks,'normalized_tables':16,'all_zero_constant_gauges':128,'all_dimension3_boolean_functions':256,'odd_functions':odd,'single_defect_functions':defects,'GL3_order':len(GL),'GL3_relative_gauge_orbit_sizes':sorted(map(len,orbits)),'uniform_literal_controls':uniform,'walsh_profiles':dict(profiles),'tables':rows},indent=2))
