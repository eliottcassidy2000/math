#!/usr/bin/env python3
"""A finite hostile for naive old-cut/complement re-marking of dyadic sources."""
import json
from pathlib import Path
from audit_power_sources import path_counts,marked,append_zero
HERE=Path(__file__).resolve().parent

def main():
 recs=[];cases=0;successes=[]
 data=json.loads((HERE/'power_sources_complete8.json').read_text())['hits'][:12]
 for rec in data:
  n=len(rec['rows']);parent,x=marked(rec['rows'],rec['mask']);h,ends=path_counts(parent);rows=[]
  for name,mask in [('none',0),('old_out',rec['mask']),('old_in',((1<<n)-1)^rec['mask']),('all_old', (1<<n)-1)]:
   for flag in (0,1):
    cut=mask|(flag<<x);child,y=marked(parent,cut);F,_=path_counts(child);g=sum(ends[v] for v in range(n+1) if cut>>v&1);H,_=path_counts(append_zero(child,y));q=H-F-g
    row={'old_deficit':rec['deficit'],'cut_rule':name,'old_mark_in_cut':flag,'state':[F,q,g,h],'dyadic':bool(q and q&(q-1)==0),'new_deficit':q-F-g};rows.append(row);cases+=1
    if row['dyadic']:successes.append(row)
  recs.append({'source_deficit':rec['deficit'],'rows':rows})
 print(json.dumps({'status':'finite_exact','actual_cut_tests':cases,'dyadic_survivors':successes,'records':recs,'scope':'These eight uniform re-markings of twelve fixed sources are a probe, not an exclusion of general re-marking. Higher cover response remains missing from the scalar state.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
