#!/usr/bin/env python3
"""Replay only this session's new construction engines and frozen receipts."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys
p=argparse.ArgumentParser();p.add_argument('--optimized',action='store_true');args=p.parse_args()
root=Path(__file__).resolve().parent
jobs=[
 ('geometry/check_padding.py',[],['geometry/check_padding.out'],'geometry/check_padding.out'),
 ('gems/local_image.py',[],['gems/receipt.json'],'gems/receipt.json'),
 ('cycles/explore.py',[],['cycles/results.json','cycles/check.out'],'cycles/check.out'),
 ('ternary/compiler.py',[],['ternary/results.json'],'ternary/results.json'),
 ('ternary/explore.py',['--bound','262143'],['ternary/explore-results.json'],'ternary/explore-results.json'),
 ('ternary/repair_holes.py',[],['ternary/hole-repairs.json'],None),
]
def digest(data):return hashlib.sha256(data).hexdigest()
checks=[]
for script,extra,receipts,text_receipt in jobs:
 before={name:(root/name).read_bytes() for name in receipts}
 cmd=[sys.executable]+(['-O'] if args.optimized else [])+['-B',str(root/script),*extra]
 result=subprocess.run(cmd,check=True,capture_output=True)
 for name,data in before.items():
  if (root/name).read_bytes()!=data:raise RuntimeError('Receipt changed: '+name)
 if text_receipt and result.stdout!=before[text_receipt]:raise RuntimeError('Output differs: '+text_receipt)
 checks.append({'script':script,'receipts':{n:digest(b) for n,b in before.items()},'passed':True})
print(json.dumps({'scope':'Six new engines; no inherited large census rerun.','checks':checks},indent=2,sort_keys=True))
