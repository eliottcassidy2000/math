from itertools import product
from math import prod
import json

def need(b,msg):
    if not b: raise RuntimeError(msg)
def vp(x,p):
    if not x: return None
    j=0
    while x%p==0: j+=1; x//=p
    return j
cases=kernel=transfers=0
for p in (2,3,5):
 for S in range(1,8):
  for j in range(S+2):
   tau=j*(j+1)//2
   for r in range(1,S+1):
    q=prod(p**r-p**i for i in range(1,j+1))
    need(q==0 if r<=j else vp(q,p)==tau,'Newton diagonal valuation')
    if r>j: need(tau==sum(vp(p**r-p**i,p) for i in range(1,j+1)),'sum of LCA depths')
    kernel+=1
   for K in range(1,5):
    scale=p**max(K-tau,0)
    # Kernel columns for j>=S vanish exactly, others use their precision.
    for a in range(S):
     for Sp in range(1,S-a+1):
      for delta in range(3):
       for r in range(1,Sp+1):
        q=prod(p**(a+r)-p**i for i in range(1,j+1))
        need((scale*p**(r*delta)*q)%(p**K)==0,'kernel under horizon-compatible substitution')
        transfers+=1
for p in (2,3):
 for K in range(1,4):
  modulus=p**K
  for h in range(1,4):
   for S in range(1,5):
    observed=set()
    for coeff in product(range(modulus),repeat=h):
     observed.add(tuple(sum(c*p**(r*d) for d,c in enumerate(coeff))%modulus for r in range(1,S+1)))
    expected=p**sum(max(K-j*(j+1)//2,0) for j in range(min(h,S)))
    need(len(observed)==expected,'full coefficient-lattice image cardinality')
    cases+=1
print(json.dumps({'status':'PASS','prime_nodes':[2,3,5],'diagonal_lca_checks':kernel,'kernel_transition_evaluations':transfers,'exhaustive_image_cases':cases,'scope':'Independent observation-module and LCA bridge, not a replay of full clone recurrence.'},sort_keys=True,indent=2))
