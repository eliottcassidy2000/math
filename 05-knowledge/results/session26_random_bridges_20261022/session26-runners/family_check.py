from probe import curve
from fractions import Fraction as F
from collections import defaultdict
import json
checks=0
for N in range(2,82,2):
 _,actual=curve((1,N));m=N//2
 expected={F(0):-4+F(2,N)}
 expected.update({F(k,N+1):2+F(2,N) for k in range(1,m+1)})
 expected.update({F(k,N-1):-2+F(2,N) for k in range(1,m)})
 if actual!=dict(sorted(expected.items())):raise ValueError((N,actual,expected))
 if len(actual)-1!=N-1:raise ValueError('order')
 checks+=1
 # Independent sample: survival from the closed positive-part sum, and a
 # rational partition of actual original inequalities at threshold delta.
 for delta in [F(1,14),F(1,5),F(1,3),F(m,N+1),F(1,2)]:
  predicted=sum(c*max(h-delta,0) for h,c in expected.items())
  walls=sorted({F(0),F(1)}|{(F(k)+sgn*delta)/a for a in (1,N) for k in range(a+1) for sgn in (-1,1) if 0<=(F(k)+sgn*delta)/a<=1})
  def dist(x):return min(x-x.__floor__(),1-(x-x.__floor__()))
  direct=sum((r-l for l,r in zip(walls,walls[1:]) if min(dist(a*(r+l)/2) for a in (1,N))>=delta),F(0))
  if direct!=predicted:raise ValueError(('mass',N,delta,direct,predicted))
  checks+=1
print(json.dumps({'family':'(1,N), all even N=2,...,80','family_formula_checks':40,'independent_fixed_threshold_wall_checks':200,'checks':checks,'unbounded_minimal_order':'N-1','status':'FINITE-EXACT verification of independently proved formula'},indent=2))
