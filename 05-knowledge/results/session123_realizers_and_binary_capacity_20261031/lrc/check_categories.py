from fractions import Fraction as F
from itertools import product
import json

def need(x,m):
 if not x:raise RuntimeError(m)
def norm(x):x=x%1;return min(x,1-x)
count=0;U=15
for A,B,C in product(range(16),repeat=3):
 imp=(U^A)|B
 need(((C&A)&~B==0)==(C&~imp==0),'residual adjunction');count+=1
need(bool(1) and bool(2) and not bool(1&2),'existential nonmonoidal hostile')
r=F(1,4);points=[F(i,4) for i in range(5)];sets=[]
for phases in [(F(0),F(1,2)),(F(1,4),F(3,4))]:
 S=[t for t in points if all(norm(t+p)>=r for p in phases)];sets.append(S)
need(set(sets[0]).isdisjoint(sets[1]),'phase loss hostile')
# Local section formula and winding arithmetic are controls, not the topology proof.
local=0
for m in range(2,20):
 for phi in [F(i,17) for i in range(-20,21)]:
  t=(F(1,2)-phi)/m;need((m*t+phi)%1==F(1,2),'local section');local+=1
 for d in range(-20,21):need(m*d+1!=0,'integer degree obstruction')
print(json.dumps({'status':'PASS','finite_residual_adjoint_controls':count,'existential_nonmonoidal_hostile':True,'same_budget_disjoint_witnesses':[[str(t) for t in s] for s in sets],'local_section_controls':local,'degree_integer_controls':18*41,'scope':'finite implementation controls; continuous no-section theorem rests on the explicit winding proof'},indent=2,sort_keys=True))
