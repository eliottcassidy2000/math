#!/usr/bin/env python3
"""FINITE-EXACT optional probe: flip the terminal involution of every 2-element."""
# This probe reuses only our independently reviewed adjacency-matrix definitions.
# Parse definitions without running the standalone marked-check output section.
from pathlib import Path
import json
from collections import Counter
from fractions import Fraction
P=Path(__file__).resolve().parent
# Functions need a shared global namespace when exec is used.
source=(P/'independent_marked_check.py').read_text().split('j=json.loads(')[0]
namespace={'__file__':str(P/'independent_marked_check.py')};exec(source,namespace)
records=[]
for case in json.loads((P/'independent_carrier.json').read_text())['cases']:
 n=case['n'];rows=[]
 for row in case['rows']:
  s=row['state'];a=namespace['matrix'](n,s)
  for xrow in row['two_element_classes']:
   x=tuple(xrow['element']);g=x
   while namespace['order'](g)>2:g=namespace['mul'](g,g)
   t=namespace['flip'](a,g)
   candidates=[h for h in namespace['centralizer_fpf'](g) if namespace['mul'](h,x)==namespace['mul'](x,h)]
   witnesses=[h for h in candidates if namespace['relabel'](a,h)==t]
   rows.append({'state':s,'element':x,'order':xrow['order'],'weight':xrow['weight'],'fixed':bool(witnesses),'one_witness':witnesses[:1]})
 records.append({'n':n,'pairs':len(rows),'fixed':sum(r['fixed'] for r in rows),'fixed_by_order':dict(Counter(r['order'] for r in rows if r['fixed'])),'fixed_mass':str(sum(Fraction(r['weight']) for r in rows if r['fixed'])),'rows':rows})
print(json.dumps({'status':'FINITE-EXACT optional all-2-element terminal-involution flip; not a pairing theorem','cases':records},indent=2,sort_keys=True))
