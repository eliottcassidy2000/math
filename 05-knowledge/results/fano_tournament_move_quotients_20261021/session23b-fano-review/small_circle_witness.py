"""Explicit real-sphere selected-circle realization, not projective lines."""
from fractions import Fraction as F
from itertools import combinations
import json

def check(q,s):
    if not q: raise RuntimeError(s)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def minus(x,y):return tuple(a-b for a,b in zip(x,y))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
# Inverse stereographic projection of (t,t^2), t=1,...,7.
points={t:(F(2*t,t**4+t*t+1),F(2*t*t,t**4+t*t+1),F(t**4+t*t-1,t**4+t*t+1)) for t in range(1,8)}
check(all(dot(p,p)==1 for p in points.values()),'unit sphere')
check(len(set(points.values()))==7,'distinct points')
planes={}
for T in combinations(points,3):
    x,y,z=(points[t] for t in T); normal=cross(minus(y,x),minus(z,x)); offset=dot(normal,x)
    check(dot(normal,normal)>0,'noncollinear triple')
    check(all(dot(normal,p)!=offset for t,p in points.items() if t not in T),'no fourth point on selected circle')
    check(offset!=0,'all 35 triple circles are small, not great')
    check(offset*offset<dot(normal,normal),'positive circle radius')
    if T[0]^T[1]^T[2]==0:planes[T]=(normal,offset)
check(len(planes)==7,'seven designated Fano circles')
print(json.dumps({'status':'PASS','model':'Seven labelled unit-sphere points with seven designated small circles; not a linear or projective Fano representation.','points':{str(t):list(map(str,p)) for t,p in points.items()},'selected_circle_planes':{str(t):{'normal':list(map(str,n)),'offset':str(d)} for t,(n,d) in planes.items()},'all_35_triple_planes_have_no_fourth_point':True,'all_35_triple_circles_are_small_not_great':True,'extra_unlabelled_circle_intersections_not_part_of_incidence_structure':True},indent=2,sort_keys=True))
