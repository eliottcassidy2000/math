#!/usr/bin/env python3
import json,sympy as s
z,t,c=s.symbols('z t c');types=[[(0,1),(4,7)],[(1,1),(3,4),(5,3)],[(2,3),(4,4),(6,1)],[(3,7),(7,1)]]
rows=[]
for w,c3 in zip(types,[-42,-6,6,42]):
 B=s.diag(*[1+2*t*k for k,n in w])-s.ones(len(w),1)*s.Matrix([[n for k,n in w]])/8
 p=B.charpoly(z).as_expr()
 trial=7*t-7*t*t+c3*t**3+c*t**4
 series=s.series(p.subs(z,trial),t,0,5).removeO().expand()
 if any(series.coeff(t,k)!=0 for k in range(4)):raise ValueError('lower terms')
 roots=s.solve(series.coeff(t,4),c)
 if len(roots)!=1:raise ValueError('nonunique fourth coefficient')
 rows.append(int(roots[0]))
if rows!=[-203,-11,-11,-203]:raise ValueError('fourth order mismatch')
print(json.dumps({'status':'PASS','fourth_coefficients':rows,'line_model_gap':'36 h^3/J^2 + 192 h^4/J^3 + O(h^5/J^4)','scope':'exact implicit-polynomial coefficients; fixed h>0, J/h -> infinity'},indent=2,sort_keys=True))
