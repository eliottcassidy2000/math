from itertools import product
import json

def grid(C,x,y):return (C+x,C-x-y,C+y,C-x+y,C,C+x-y,C-y,C+x+y,C-x)
def rot(g):return tuple(g[(2-j)*3+i] for i in range(3) for j in range(3))
def refl(g):return tuple(g[i*3+2-j] for i in range(3) for j in range(3))
cases=0
for x in range(2,31):
 for y in range(1,x):
  if x==2*y:continue
  vals=set(grid(100,x,y));assert len(vals)==9
  layouts={grid(100,u-100,v-100) for u,v in product(vals,repeat=2) if set(grid(100,u-100,v-100))==vals}
  g=grid(100,x,y);orb=set()
  for _ in range(4):orb|={g,refl(g)};g=rot(g)
  assert len(layouts)==8 and layouts==orb
  cases+=1
# All 3-by-3 magic layouts of the fixed vector plane F_3^2.
F=list(product(range(3),repeat=2))
def add(a,b):return tuple((u+v)%3 for u,v in zip(a,b))
def neg(a):return tuple(-u%3 for u in a)
def g3(c,x,y):
 return tuple(add(c,add(tuple(i*u%3 for u in x),tuple(j*v%3 for v in y))) for i,j in [(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)])
ls={g3(c,x,y) for c,x,y in product(F,repeat=3) if (x[0]*y[1]-x[1]*y[0])%3}
assert len(ls)==432
assert all(len(set(g))==9 for g in ls)
print(json.dumps({'real_parameter_cases_checked':cases,'layouts_each':8,'char3_plane_layouts':len(ls),'char3_D4_orbits':len(ls)//8},indent=2))
