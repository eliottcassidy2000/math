#!/usr/bin/env python3
"""Independent exact checkpoint2 controls; infinite conclusions have written audit."""
from fractions import Fraction as F
from itertools import combinations, product
import hashlib,json,pathlib,collections
BASE=pathlib.Path(__file__).resolve().parents[2]
N=0

def check(p,why):
 global N
 N+=1
 if not p: raise RuntimeError(why)

def xy(k,j):return F(2*j+1,2**(k+1)),F(1,2**k)
def coshd(p,q):
 x,y=p;u,v=q
 return 1+((x-u)**2+(y-v)**2)/(2*y*v)
V={(k,j):xy(k,j) for k in range(7) for j in range(2**k)}
E={( (k,j),(k+1,2*j+s)) for k in range(6) for j in range(2**k) for s in (0,1)}
pair_count=0
for a,b in combinations(V,2):
 d=coshd(V[a],V[b]);is_edge=(a,b) in E or (b,a) in E
 check(d>=F(21,16),'sharp separated metric')
 check((d==F(21,16))==is_edge,'induced fixed-distance tree')
 pair_count+=1
rects=[]
for a,b in E:
 x,y=V[a];u,v=V[b]
 center=(u*u+v*v-x*x-y*y)/(2*(u-x))
 sign=1 if u>x else -1
 check(center==x-sign*F(11,8)*y,'circle center from generic endpoint equation')
 radius2=(x-center)**2+y*y
 check(radius2==F(185,64)*y*y and radius2==(u-center)**2+v*v,'circle radius')
 check((x-center)*(u-center)>0,'same monotone semicircle side')
 check(coshd(V[a],V[b])==F(21,16),'equal hyperbolic edges')
 rects.append((a,b,min(x,u),max(x,u),v,y))
edge_pairs=0
for A,B in combinations(rects,2):
 a,b,l,r,lo,hi=A;c,d,s,t,bottom,top=B
 if max(lo,bottom)<min(hi,top):
  if a==c:
   check(max(l,s)==min(r,t)==V[a][0],'siblings intervals meet only at parent abscissa')
  else:check(r<s or t<l,'same band distinct parents separated')
 elif max(lo,bottom)==min(hi,top):
  # Shared band boundary can contain only endpoints of these monotone arcs.
  h=max(lo,bottom)
  ea={p for p in (a,b) if V[p][1]==h};eb={p for p in (c,d) if V[p][1]==h}
  check(all(V[p]!=V[q] or p==q for p in ea for q in eb),'boundary intersection only genuine shared vertex')
 edge_pairs+=1
for k in range(1,65):
 a=(k,2**(k-1)-1);b=(k,2**(k-1))
 check(coshd(xy(*a),xy(*b))==F(3,2),'fixed ambient distance across root split')
 check(a[1]<2**(k-1)<=b[1],'tree distance exactly twice depth')
 check(xy(*b)[0]-xy(*a)[0]==F(1,2**k),'Euclidean separation disappears')

# Construct a small exact rational induced tree using upward arcs, finite
# exclusions, and nested intervals, independently of curvature coordinates.
def horizontal(t):return 2*t/(1+t*t)
def parameter_inside(lo,hi):
 a,b=F(-1,2),F(1,2)
 for _ in range(1000):
  c=(a+b)/2;v=horizontal(c)
  if lo<v<hi:return c
  if v<=lo:a=c
  else:b=c
 raise RuntimeError('rational parameter bisection exhausted')
positions={0:(F(0),F(0))};intervals={0:(F(-1,8),F(1,8))};parents={};max_candidates=0
for p in range(15):
 x,y=positions[p];lo,hi=intervals[p];margin=min(x-lo,hi-x,F(1,4))
 for sign in (-1,1):
  child=2*p+(1 if sign<0 else 2)
  left,right=((x-F(3,4)*margin,x-F(1,2)*margin) if sign<0 else (x+F(1,2)*margin,x+F(3,4)*margin))
  dlo,dhi=left-x,right-x
  ta=parameter_inside(dlo,(3*dlo+dhi)/4)
  tb=parameter_inside((dlo+3*dhi)/4,dhi)
  check(ta<tb,'nonempty rational arc parameter interval')
  for m in range(2,3*len(positions)+5):
   t=ta+(tb-ta)/m;dx=horizontal(t);dy=(1-t*t)/(1+t*t)
   z=(x+dx,y+dy)
   if any(z==old or (q!=p and (z[0]-old[0])**2+(z[1]-old[1])**2==1) for q,old in positions.items()):continue
   max_candidates=max(max_candidates,m-1);break
  else:raise RuntimeError('finite avoidance candidates exhausted')
  check(left<z[0]<right and lo<left<right<hi,'nested intervals')
  check(dx*dx+dy*dy==1 and dy>0 and dy*dy>F(15,16),'unit upward increment')
  positions[child]=z;intervals[child]=(left,right);parents[child]=p
for a,b in combinations(positions,2):
 x,y=positions[a];u,v=positions[b]
 check(((x-u)**2+(y-v)**2==1)==(parents.get(a)==b or parents.get(b)==a),'induced rational finite tree')
 check((x,y)!=(u,v),'rational tree injection')

# Ring enumeration, exact integer inverse radius (no floating-point ceilings),
# and one-step conjugacy are separate from local lattice edges.
dirs=((-1,1),(-1,0),(0,-1),(1,-1),(1,0),(0,1))
points=[(0,0)]
for r in range(1,101):
 x,y=r,0
 for dx,dy in dirs:
  for _ in range(r):
   check(max(abs(x),abs(y),abs(x+y))==r,'ring metric radius')
   points.append((x,y));x+=dx;y+=dy
 check((x,y)==(r,0),'ring closure')
index={p:i+1 for i,p in enumerate(points)}
check(len(points)==len(index)==30301,'ring bijection finite exhaustion')
for n,p in enumerate(points,1):
 r=max(abs(p[0]),abs(p[1]),abs(sum(p)))
 check(n<=1+3*r*(r+1) and (r==0 or n>1+3*(r-1)*r),'exact ceiling formula characterization')
for n in range(1,10001):
 successor=n//2 if n%2==0 else 3*n+1
 check(index[points[successor-1]]==successor,'exact clock conjugacy')
 if n%2==0:
  p=points[n-1];q=points[n//2-1]
  d=max(abs(p[0]-q[0]),abs(p[1]-q[1]),abs(sum(p)-sum(q)))
  rp=max(abs(p[0]),abs(p[1]),abs(sum(p)));rq=max(abs(q[0]),abs(q[1]),abs(sum(q)))
  check(d>=rp-rq,'radial jump lower bound')

# Small exact jet controls in dimension three, independent of producer's 2D test.
jet_tests=0
subsets=list(range(8))
for modulus in (2,3,5):
 for a in range(1,modulus):
  delta=min(a,modulus-a)
  for sign in (-1,0,1):
   coeff={s:modulus*sign*((s%3)-1) for s in subsets};coeff[7]=a+modulus*sign
   for x in product((F(-1),F(-1,2),F(0),F(1,2),F(1)),repeat=3):
    energy=F(0)
    for deriv in subsets:
     value=F(0)
     for s,c in coeff.items():
      if s&deriv==deriv:
       term=F(c)
       for i in range(3):
        if (s^deriv)&(1<<i):term*=x[i]
       value+=term
     energy+=value*value
    floor=F(delta*delta)
    for v in x:floor*=1+v*v
    check(energy>=floor,'3D full arithmetic jet')
    jet_tests+=1
pins=json.loads((pathlib.Path(__file__).parent/'SOURCE_PINS.json').read_text())
for entry in pins['local_packets']:
 b=(BASE/entry['path']).read_bytes()
 check(len(b)==entry['bytes'] and hashlib.sha256(b).hexdigest()==entry['sha256'],'source pin '+entry['path'])
print(json.dumps({'status':'PASS','checks':N,'hyperbolic_vertices':len(V),'hyperbolic_vertex_pairs':pair_count,'hyperbolic_edges':len(E),'hyperbolic_edge_pairs':edge_pairs,'non_quasi_isometry_depths':64,'proper_unit_rational_control_vertices':len(positions),'maximum_circle_candidates_used':max_candidates,'ring_sites':len(points),'exact_clock_tests':10000,'three_variable_jet_tests':jet_tests,'source_pins':len(pins['local_packets']),'scope':'Independent finite exact controls supplement written proofs; not Collatz convergence, full Collatz hyperbolic embedding, quasi-isometry, or full inherited theorem audit'},indent=2,sort_keys=True))
