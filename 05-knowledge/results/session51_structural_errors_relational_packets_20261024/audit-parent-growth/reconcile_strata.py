#!/usr/bin/env python3
"""Reconcile the independent atlas and frozen structural sector witnesses."""
from pathlib import Path
from collections import Counter
import json
B=Path(__file__).resolve().parent

def cells(r):
 n=len(r);partition=[list(range(n))]
 while True:
  new=[]
  for block in partition:
   groups={}
   for v in block:
    sig=tuple(sum((r[w]>>v)&1 for w in other) for other in partition)
    groups.setdefault(sig,[]).append(v)
   new.extend(groups[s] for s in sorted(groups))
  if len(new)==len(partition):return partition
  partition=new
j=json.loads((B/'structural-audit/receipt.json').read_text());strong10=Counter()
for witness in j['representatives']:
 partition=cells(witness['rows']);size={v:len(cell) for cell in partition for v in cell};actual=0
 for orbit in witness['vertex_orbits']:
  d=size[orbit[0]]-len(orbit);actual+=d
  if d:strong10[size[orbit[0]]]+=d
 if actual!=witness['error']:raise RuntimeError('witness error mismatch')
a=json.loads((B/'independent-O2.json').read_text());strong9=Counter(dict(a['strong_order_nine_error_by_stable_cell_size']));e10=strong10.copy()
for size,value in strong9.items():e10[size]+=2*value
e10[7]+=6*j['scc_transfer']['inherited_strong_error_7']
scan=json.loads((B/'incoming-O2.json').read_text());hist=dict(scan['stable_root_histogram']);corrected={size:size*count-e10[size] for size,count in hist.items()}
if sum(e10.values())!=5508 or sum(corrected.values())!=10*scan['prediction']:raise RuntimeError('stratified trace mismatch')
result={'status':'PASS','strong_ten_error_by_stable_cell_size':sorted(strong10.items()),'independent_strong_nine_error_by_stable_cell_size':sorted(strong9.items()),'total_ten_error_by_stable_cell_size':sorted(e10.items()),'root_orbit_mass_by_stable_cell_size':sorted(corrected.items()),'total_error':sum(e10.values()),'corrected_trace':sum(corrected.values()),'prediction':sum(corrected.values())//10,'dependency':'Strong-seven error 14 from the frozen independently replayed structural packet; SCC multipliers 2 and 6 from A(z)^2.'}
print(json.dumps(result,indent=2))
