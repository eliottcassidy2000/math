#!/usr/bin/env python3
"""Re-execute the earlier independent full filtration, not producer mathematics."""
import argparse,importlib.util,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
spec=importlib.util.spec_from_file_location('prior_independent_audit',HERE/'inputs/guard-theorem/audit.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class SnapshotRepo:
 def __truediv__(self,name):
  if name!='05-knowledge/results':raise RuntimeError('Unapproved original source lookup')
  return HERE/'inputs/original'
filtration,labels=m.independent_filtration(SnapshotRepo(),HERE/'inputs/guard-theorem')
saved=json.loads((HERE/'inputs/guard-theorem/normal.json').read_text())['filtration']
if json.dumps(filtration,sort_keys=True)!=json.dumps(saved,sort_keys=True):raise RuntimeError('Filtration replay mismatch')
r={'status':'PASS: replay of previously independent full original 2-primary filtration and marked harmonic completeness','checks':m.COUNT,'filtration':filtration,'scope':'2-primary group and all paired dyadic characters only; full critical-group exponent not asserted.'}
a.out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':m.COUNT,'output':str(a.out)}))
