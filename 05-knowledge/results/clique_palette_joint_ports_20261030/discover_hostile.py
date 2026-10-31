"""Find a small hostile to exporting Paley palette inequalities."""
from random import Random
from itertools import combinations
from scipy.optimize import linprog
import json
rng=Random(1314)
for trial in range(10000):
    rows=[0]*7
    for i,j in combinations(range(7),2):
        if rng.randrange(2): rows[i] |= 1 << j
        else: rows[j] |= 1 << i
    cyc=[sum(1 << i for i in t) for t in combinations(range(7),3)
         if all(sum(rows[i] >> j & 1 for j in t)==1 for i in t)]
    good=[s for s in range(1,128) if not any(s & c == c for c in cyc)]
    if max(s.bit_count() for s in good)!=4 or any((127^s) in good for s in good):
        continue
    lp=linprog([1]*len(good), A_ub=[[-(s >> i & 1) for s in good] for i in range(7)],
               b_ub=[-1]*7,bounds=(0,None),method='highs')
    from fractions import Fraction as F
    print(json.dumps({'trial':trial,'rows':rows,'cyclic_triples':cyc,'fractional_discovery':lp.fun,
       'positive_weights':[(s,str(F(float(x)).limit_denominator(1000))) for s,x in zip(good,lp.x) if x>1e-8],
       'dual':[str(F(float(-x)).limit_denominator(1000)) for x in lp.ineqlin.marginals]},indent=2))
    break
