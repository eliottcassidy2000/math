#!/usr/bin/env python3
"""Classify the NEW bounded small-parent holes by explicit H9 repairs."""
import json
from pathlib import Path
from compiler import signature,generate,construct,hp_subsets
root=Path(__file__).parent
prior=json.loads((root/'explore-results.json').read_text())
holes=prior['targets_requiring_parent9_in_this_bank'];bound=prior['bound']
cores=[(2,4,9,19,7),(2,4,25,3,11)]
repairs={h:[[],[]] for h in holes}
for ci,core in enumerate(cores):
 for O in range(1<<5):
  sig=signature(core,O)
  for val,a,b in generate(sig,max(holes)):
   if val in repairs:repairs[val][ci].append((a.bit_length()+b.bit_length(),a+b,a,b,O,sig))
chosen={};cut_hist=[{},{}];missing=[[],[]]
for h,rs in repairs.items():
 entries=[]
 for ci,arr in enumerate(rs):
  if not arr:missing[ci].append(h);entries.append(None);continue
  score,total,a,b,O,sig=min(arr)
  child=construct(cores[ci],O,a,b)
  if hp_subsets(child)[-1]!=h:raise AssertionError(h)
  key=str(O);cut_hist[ci][key]=cut_hist[ci].get(key,0)+1
  entries.append({'A':a,'B':b,'O':O,'signature':sig,'child_rows':child,'number_parameter_solutions':len(arr)})
 chosen[str(h)]={'binary':bin(h)[2:],'repairs':entries}
levels={str(k):[h for h in holes if h<2**k] for k in range(8,bound.bit_length()+1)}
# Explicit hostile to replacing the sparse holes by their doubling orbit.
# Test common small affine descendants, not fit a model to the finite data.
affine={str(c):{'hits':[(h,2*h+c) for h in holes if 2*h+c in holes], 'misses':[(h,2*h+c) for h in holes if 1<2*h+c<=bound and 2*h+c not in holes]} for c in (-3,-1,1,3)}
report={'bound':bound,'scope':'Complete small-parent1/3/5 hole list within bound, using prior exact compiler. Every displayed positive repair directly recounted. No all-height regular-language assertion.', 'power_two_prefix_holes':levels,'holes_missed_by_each_single_H9_core':missing,'minimal_order_repair_cut_histograms':cut_hist,'repairs':chosen,'small_affine_hole_tests':affine}
(root/'hole-repairs.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(json.dumps({'missing_per_core':missing,'histograms':cut_hist,'prefix_counts':{k:len(v) for k,v in levels.items()}},sort_keys=True,indent=2))
