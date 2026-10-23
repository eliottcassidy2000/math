#!/usr/bin/env python3
"""All positive odd global color conductances reduce to four ratios modulo8."""
import json
from pathlib import Path
from check_dyadic_network import lap,need,val2
from check_dyadic_lift import dyadic_smith
HERE=Path(__file__).resolve().parent
out=[]
for a in json.loads((HERE/'dyadic_network_receipt.json').read_text()):
 n=a['nodes'];L=lap(n,a['blue_edges'])[:-1,:-1];B=lap(n,a['paired_black_edges'])[:-1,:-1];cases=[]
 for t in (1,3,5,7):
  vals=dyadic_smith((t*L+2*B).tolist());need(max(vals,default=0)<=2)
  cases.append(dict(odd_ratio_mod8=t,valuation_histogram={str(k):vals.count(k) for k in sorted(set(vals))}))
 need(all(c['valuation_histogram']==cases[0]['valuation_histogram'] for c in cases))
 out.append(dict(n=a['n'],cases=cases,all_positive_odd_color_weights_same_two_primary_group=True))
 print(out[-1],flush=True)
(HERE/'odd_color_weights_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print('ALL CHECKS PASSED')
