"""Independent five exact algebra controls for P5; no assertions or floats."""
import json
import sympy as s

phi=(1+s.sqrt(5))/2
q=2-phi
D=[(0,0),(1,0),(q,q),(0,1)]
K=[(1,1),(1,0),(q,q),(0,1)]
def area(p):
    return s.simplify(sum(p[i][0]*p[(i+1)%len(p)][1]-
                         p[i][1]*p[(i+1)%len(p)][0]
                         for i in range(len(p)))/2)

tests=[area(D)==q,
       s.simplify(area(K)-(q-1))==0,
       s.simplify((1-q)/q-phi)==0,
       s.simplify(phi*(1-phi)+1)==0,
       s.simplify(s.sqrt(phi)**2-phi)==0]
if not all(tests):
    raise RuntimeError(tests)
print(json.dumps({'status':'PASS','checks':len(tests),
    'dart_oriented_oblique_area':str(area(D)),
    'kite_oriented_oblique_area':str(area(K)),
    'physical_area_ratio':str(s.simplify(-area(K)/area(D))),
    'scope':'Exact template-area controls; field theorem audited by proof.'},indent=2))
