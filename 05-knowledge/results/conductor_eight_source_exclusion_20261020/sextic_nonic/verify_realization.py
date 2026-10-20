#!/usr/bin/env python3
"""Exact symbolic checks for the (9,6,2) polynomial realization."""
import json
from pathlib import Path
import sympy as s

out=Path(__file__).resolve().parent
t,z,v,U,V=s.symbols('t z v U V')
f=t**6+4*t**2
p=z**4+6*z**2+6
q=t*p.subs(z,t**2)
checks={}
checks['degrees_6_9']=(s.degree(f,t),s.degree(q,t))==(6,9)
checks['characteristic_identity']=s.expand(q**2-f**3-8*f)==4*t**2
checks['nonzero_p_at_zero']=p.subs(z,0)==6
checks['squarefree_p']=s.gcd(p,s.diff(p,z))==1
checks['immersion_original_pair']=s.gcd(s.diff(f,t),s.diff(q,t))==1
checks['conductor_degree_8']=s.degree(p.subs(z,t**2),t)==8
checks['defining_equation']=s.expand((V**2-U**3-8*U)**3+64*(V**2-U**3-8*U)-64*U).subs({U:f,V:q}).expand()==0
fv=t**6+v*t**2
pv=z**4+s.Rational(3,2)*v*z**2+s.Rational(3,8)*v**2
qv=t*pv.subs(z,t**2)
checks['family_identity']=s.expand(qv**2-fv**3-v**3*fv/8-v**4*t**2/64)==0
checks['family_critical_value']=s.rem(pv,4*z**2+3*v,z)==-3*v**2/16
checks['family_nonzero_zero_value']=pv.subs(z,0)==3*v**2/8
# Finite semigroup check through a bound exceeding the conductor. The full
# assertion is proved by the even/odd direct-sum decomposition in the note.
semigroup={2*a+9*b for a in range(40) for b in range(10)}
checks['gaps_through_40']=[n for n in range(41) if n not in semigroup]==[1,3,5,7]
result={'all_passed':all(checks.values()),'checks':checks,'gap_count':4,'conductor_degree':8,'coefficient_field':'Q, hence every characteristic-zero field'}
if not result['all_passed']:
    raise RuntimeError(result)
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
