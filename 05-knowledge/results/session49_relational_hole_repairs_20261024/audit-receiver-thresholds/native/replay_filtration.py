#!/usr/bin/env python3
"""Replay the unchanged earlier independent filtration audit on this closed input set."""
from pathlib import Path
import importlib.util,json,sys
P=Path(__file__).resolve().parent;sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('prior_independent_filtration',P/'inputs/prior-filtration/audit.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Snapshot:
 def __truediv__(self,name):
  if name!='05-knowledge/results':raise RuntimeError('unexpected original lookup')
  return P/'inputs'
f,labels=m.independent_filtration(Snapshot(),P/'inputs/prior-filtration')
saved=json.loads((P/'inputs/prior-filtration/guard-completeness.normal.json').read_text())
if json.dumps(f,sort_keys=True)!=json.dumps(saved['filtration'],sort_keys=True):raise RuntimeError('independent filtration replay mismatch')
print(json.dumps({'status':'PASS','checks':m.COUNT,'paired_charge_order':f['paired_charge_order'],'filtration':'unchanged earlier independent audit re-executed','full_critical_group_exponent_claimed':False},sort_keys=True))
