#!/usr/bin/env python3
"""Merge actual dyadic sources and audit a uniform near-power interval."""
import json
from pathlib import Path
from audit_power_sources import need,path_counts,marked,append_zero
HERE=Path(__file__).resolve().parent
FILES=['power_sources_complete8.json','power_sources.json','target_sources.json','neighbor_source.json','insertion_source.json']

def bank():
 out={}
 for name in FILES:
  for raw in json.loads((HERE/name).read_text())['hits']:
   r=dict(raw);r['provenance']=name;c=r['deficit'];q=r['state'][1];key=(q,len(r['rows']))
   if c not in out or key<(out[c]['state'][1],len(out[c]['rows'])):out[c]=r
 return out

def main():
 B=bank();missing=next(c for c in range(1,max(B)+2,2) if c not in B);D=missing-2
 chosen=[B[c] for c in range(1,missing,2)];threshold=max(r['state'][1].bit_length() for r in chosen)
 overhead=max(len(r['rows'])+1-(r['state'][1].bit_length()-1) for r in chosen)
 need((D,threshold,overhead)==(235,14,1),'current uniform bank boundary')
 controls=[]
 for r in chosen:
  a=r['state'][1].bit_length()-1
  for k in [threshold,threshold+1]:
   rows,x=marked(r['rows'],r['mask'])
   for _ in range(k-a):rows=append_zero(rows,x)
   actual,_=path_counts(rows)
   need(actual==(1<<k)-r['deficit'],'actual common-scale interval count')
   need(len(rows)<=k+overhead,'uniform physical vertex budget')
   controls.append({'deficit':r['deficit'],'k':k,'raw_H':actual,'vertices':len(rows)})
 data={'status':'passed','first_missing_deficit':missing,'complete_odd_deficits':[1,D],'source_count':len(chosen),'exponent_threshold':threshold,'vertex_overhead':overhead,'sources':chosen,'common_scale_controls':controls,'scope':'For every k>=14, all118 odd counts in[2^k-235,2^k-1] are realized on at most k+1 vertices. Deficit237 is missing from this bank, not proved globally impossible.'}
 print(json.dumps(data,indent=2,sort_keys=True))
if __name__=='__main__':main()
