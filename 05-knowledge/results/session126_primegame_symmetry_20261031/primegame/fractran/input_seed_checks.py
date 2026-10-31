#!/usr/bin/env python3
"""Finite exact input tests, distinct from numeric labels on other objects."""
from hashlib import sha256
from pathlib import Path
import json
from verify_primegame import step_integer, power_two_exponent

ROOT=Path(__file__).resolve().parent

def orbit(seed, cap=100000, target=3):
    seen={seed:0}; n=seed; prefix=[seed]; outputs=[]
    for t in range(1,cap+1):
        n,j=step_integer(n)
        if len(prefix)<40: prefix.append(n)
        p=power_two_exponent(n)
        if p is not None: outputs.append({'step':t,'exponent':p})
        if n in seen:
            cycle=[]; m=n
            while not cycle or m!=n:
                cycle.append(m);m,_=step_integer(m)
            return {'seed_integer':seed,'outcome':'proved eventual cycle by repeated exact state','transient_length':seen[n], 'period':t-seen[n], 'cycle':cycle,'prefix':prefix,'later_pure_two_powers':outputs}
        seen[n]=t
        if len(outputs)==target:
            return {'seed_integer':seed,'outcome':'finite replay stopped after requested outputs, not halted','updates':t,'prefix':prefix,'later_pure_two_powers':outputs}
    return {'seed_integer':seed,'outcome':'safety bound reached; no halt/cycle claim','updates':cap,'prefix':prefix,'later_pure_two_powers':outputs}

results=[orbit(x) for x in (7,17,55,15,16,23,25,343,60)]
# Exact expected cycle witnesses and requested five-seed outcomes.
expect={7:[55,65],17:[55,65],55:[55,65],25:[1375,1625],343:[55,65]}
for r in results:
    if r['seed_integer'] in expect and r.get('cycle')!=expect[r['seed_integer']]:
        raise RuntimeError('cycle mismatch')
(ROOT/'input_seed_results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
