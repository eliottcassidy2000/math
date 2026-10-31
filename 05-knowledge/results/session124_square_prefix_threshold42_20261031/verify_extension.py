#!/usr/bin/env python3
"""Replay threshold42 extension against its immutable first-checkpoint dependency."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
from core42 import create_core,DEFAULT_CORE,DEFAULT_BASE,HERE,threshold
p=argparse.ArgumentParser();p.add_argument('--core-dir',type=Path,default=DEFAULT_CORE);p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);a=p.parse_args()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
manifest=json.loads((HERE/'MANIFEST.json').read_text())
for name,h in manifest['proof_files'].items():require(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,'extension proof hash:'+name)
deps=json.loads((HERE/'DEPENDENCIES.json').read_text())
for name,h in deps['frozen_core_files'].items():require(hashlib.sha256((a.core_dir/name).read_bytes()).hexdigest()==h,'frozen dependency:'+name)
require(hashlib.sha256((HERE/'SMALL_BASES.json').read_bytes()).hexdigest()==deps['small_base_sha256'],'small base hash')
core=create_core(a.core_dir,a.base_dir);base=core.verify_bases()
# Replay complete universal certificate catalogs through frozen checkers.
paired_checks=direct_checks=0
for parity in [0,1]:
 for d in range(24,73):
  for folder,name,checker in [('paired',f'B_p{parity}_d{d}.json','verify_paired.py'),('cycles',f'cycle_p{parity}_d{d}.json','verify_direct.py')]:
   cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(a.core_dir/checker),str(a.core_dir/folder/name)]
   r=json.loads(subprocess.check_output(cmd));require(r['status']=='passed','universal replay')
   if folder=='paired':paired_checks+=r['checks']
   else:direct_checks+=r['checks']
controls=[create_core(a.core_dir,a.base_dir).cycle(n) for n in [2082,2083,2620,2621,102042,102043,5000082]]
# Check cutoff and prefix at BOTH sides of the transition; 2081 is only weakB.
cutoff=[]
for n in [2081,2082,2083,2620,2621,102042,102043]:
 c=create_core(a.core_dir,a.base_dir);pair=c.build_B(n);depth=0
 while threshold(depth+1)<=n:depth+=1
 c.old.check_pair(pair,c.original.prefix(depth));cutoff.append({'n':n,'depth':depth,'checks':c.checks})
hostile=0
for kind in ['low_B','low_cycle','bad_prefix','bad_parity']:
 c=create_core(a.core_dir,a.base_dir)
 try:
  if kind=='low_B':c.build_B(41)
  elif kind=='low_cycle':c.cycle(2081)
  else:
   pair=[list(x) for x in c.small['51']]
   if kind=='bad_prefix':pair[0][0],pair[0][1]=pair[0][1],pair[0][0]
   else:pair[1][14],pair[1][15]=pair[1][15],pair[1][14]
   c.old.check_pair(pair)
 except RuntimeError:hostile+=1
require(hostile==4,'hostile controls')
print(json.dumps({'hostile_controls':hostile,'status':'passed','base_report':base,'paired_states':98,'direct_states':98,'paired_checks':paired_checks,'direct_checks':direct_checks,'cycle_controls':controls,'cutoff_controls':cutoff,'thresholds':[threshold(i) for i in range(4)]},indent=2))
