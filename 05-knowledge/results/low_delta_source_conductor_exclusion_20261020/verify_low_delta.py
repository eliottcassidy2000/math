#!/usr/bin/env python3
"""Finite exact controls; mathematical proof is in the accompanying note."""
import json
from pathlib import Path
import sympy as s

y, u, v, t, z, a, b = s.symbols('y u v t z a b')
f, g = y**3, y**4+y
R = v**3-u*(u+1)**3
checks = {}
def check(name, condition):
    if not bool(condition):
        raise RuntimeError(name)
    checks[name] = True

check('sharp_equation', s.expand(R.subs({u:f,v:g})) == 0)
check('sharp_birational_identity', s.expand(g-y*(f+1)) == 0)
check('sharp_immersion', s.gcd(s.diff(f,y),s.diff(g,y)) == 1)
c = s.cancel(s.diff(R,v).subs({u:f,v:g}) / s.diff(f,y))
check('sharp_conductor', s.expand(c-(y**3+1)**2) == 0)
U = s.cancel((t**3-y**3)/(t-y))
V = s.cancel((t**4+t-y**4-y)/(t-y))
check('sharp_resultant_conductor', s.expand(s.resultant(U,V,t)-c) == 0)
gaps = [n for n in range(18) if not any(3*i+4*j == n for i in range(7) for j in range(6))]
check('sharp_semigroup_gaps', gaps == [1,2,5])
check('sharp_delta_degree', s.degree(c,y) == 2*len(gaps))
# Root-supplied three-node sharp boundary.
f2, g2 = y**3-y, y**4
R2 = v**3-2*v**2+v-4*u**2*v-u**4
c2 = (y*y-1)*(y**4+1)
check('three_node_equation', s.expand(R2.subs({u:f2,v:g2})) == 0)
check('three_node_immersion', s.gcd(s.diff(f2,y),s.diff(g2,y)) == 1)
check('three_node_conductor', s.expand(s.cancel(s.diff(R2,v).subs({u:f2,v:g2})/s.diff(f2,y))-c2) == 0)
check('three_node_conductor_squarefree', s.gcd(c2,s.diff(c2,y)) == 1)
U2=s.cancel((f2.subs(y,t)-f2)/(t-y))
V2=s.cancel((g2.subs(y,t)-g2)/(t-y))
check('three_node_resultant', s.expand(s.resultant(U2,V2,t)-c2) == 0)
sumvar,prodvar=s.symbols('sumvar prodvar')
# In symmetric pair coordinates U2=s²−p−1, V2=s(s²−2p).
check('three_node_pair_sums', s.expand((sumvar*(sumvar**2-2*prodvar)).subs(prodvar,sumvar**2-1)-sumvar*(2-sumvar**2)) == 0)
check('three_node_pair_discriminants', all(s.simplify(sv**2-4*(sv**2-1)) != 0 for sv in (0,s.sqrt(2),-s.sqrt(2))))
check('three_node_disjoint_pair_product', s.expand((t*t-s.sqrt(2)*t+1)*(t*t+s.sqrt(2)*t+1)-(t**4+1)) == 0)
D=s.diff(f2,y).subs(y,t)*s.diff(g2,y)-s.diff(f2,y)*s.diff(g2,y).subs(y,t)
check('three_node_tangent_formula', s.expand(D-4*(y-t)*(3*(y*t)**2-(y*y+y*t+t*t))) == 0)
for d in range(1,5):
    p = s.prod(z-j for j in range(1,d+1))
    F=t**2
    G=t*p.subs(z,t**2)
    eq=v**2-u*p.subs(z,u)**2
    C=s.cancel(s.diff(eq,v).subs({u:F,v:G}) / s.diff(F,t))
    check(f'quadratic_equation_d{d}', s.expand(eq.subs({u:F,v:G})) == 0)
    check(f'quadratic_conductor_d{d}', s.expand(C-p.subs(z,t**2)) == 0)
    check(f'quadratic_immersion_d{d}', s.gcd(s.diff(F,t),s.diff(G,t)) == 1)
    check(f'quadratic_conductor_degree_d{d}', s.degree(C,t) == 2*d)
# All numerical semigroups whose positive gaps lie in 1..5, with genus <=3.
from itertools import combinations
semigroups=[]
for delta in range(4):
    for gaps_tuple in combinations(range(1,6),delta):
        gs=set(gaps_tuple)
        if all(i+j not in gs for i in range(7) for j in range(7) if i not in gs and j not in gs):
            F=max(gs,default=-1)
            if F==2*delta-1:
                semigroups.append({'delta':delta,'gaps':list(gaps_tuple)})
check('symmetric_genus_at_most_three', semigroups == [
    {'delta':0,'gaps':[]}, {'delta':1,'gaps':[1]},
    {'delta':2,'gaps':[1,3]}, {'delta':3,'gaps':[1,2,5]},
    {'delta':3,'gaps':[1,3,5]}])
result={'status':'FINITE-EXACT controls, not the general proof','checks':len(checks),
        'passed':checks,'sharp_example':{'f':str(f),'g':str(g),'equation':str(R),
        'conductor':str(s.factor(c)),'delta':3,'degree_gaps':gaps},
        'three_node_example':{'f':str(f2),'g':str(g2),'equation':str(R2),'conductor':str(c2),'delta':3,'pair_sums':['0','sqrt(2)','-sqrt(2)']},
        'symmetric_semigroups_genus_at_most_three':semigroups}
print(json.dumps(result,indent=2,sort_keys=True))
