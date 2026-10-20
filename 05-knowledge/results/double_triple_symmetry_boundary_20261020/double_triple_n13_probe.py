#!/usr/bin/env python3
"""Selected n13 double-triple embeddings, not a full n13 subgroup census."""
from pathlib import Path
import sys
base=Path(__file__).resolve().parents[1]
for candidate in (base,base/"mixed_prime_census_extension_20261020"):
    if (candidate/"extension_census.py").exists():
        sys.path.insert(0,str(candidate));break
else:raise RuntimeError("Missing predecessor extension engine")
import extension_census as e
import json,sys
m=e.m;n=13
c5=m.cyc(n,[list(range(5))]);c3=m.cyc(n,[[5,6,7]]);d3=m.cyc(n,[[8,9,10]])
rows=[]
for name,gens in [('C15_5plus3plus3',[c5,m.compose(c3,d3)]),('C5_times_C3_squared',[c5,c3,d3])]:
 r=m.census_case(n,name,gens);rows.append(r)
 print(name,'R',r['relative_tournament_orbits'],'T',r['tau_T'],'E',r['tau_E'],'FT',r['tournament_fusion_defect'],'FE',r['even_graph_fusion_defect'],'Z',r['sign_escape_relative_orbits'],file=sys.stderr,flush=True)
print(json.dumps({'status':'PASS','scope':'two specified 5+3+3+1+1 embeddings at n13, not a subgroup census','rows':rows},indent=2,sort_keys=True))
