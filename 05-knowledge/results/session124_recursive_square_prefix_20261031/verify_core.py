#!/usr/bin/env python3
"""Portable full proof replay using base data and the two new tables only."""
import argparse,copy,hashlib,json,subprocess,sys
from pathlib import Path
from core import Core,HERE,DEFAULT_BASE
from verify_direct import verify as verify_direct
p=argparse.ArgumentParser();p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);a=p.parse_args()
checks=0
def require(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
manifest=json.loads((HERE/'MANIFEST.json').read_text())
for name,h in manifest['proof_files'].items():require(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,'proof hash:'+name)
require({x.name for x in(HERE/'paired').glob('*.json')}=={f'B_p{p}_d{d}.json' for p in[0,1] for d in range(24,73)},'paired catalog')
require({x.name for x in(HERE/'cycles').glob('*.json')}=={f'cycle_p{p}_d{d}.json' for p in[0,1] for d in range(24,73)},'cycle catalog')
pairs={};pair_checks=0;cycle_checks=0
for par in[0,1]:
 for d in range(24,73):
  path=HERE/'paired'/f'B_p{par}_d{d}.json';b=json.loads(path.read_text());require((b['parity'],b['residue'])==(par,d),'paired filename/state');pairs[par,d]=b
  cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(HERE/'verify_paired.py'),str(path)]
  r=json.loads(subprocess.check_output(cmd));require(r['status']=='passed','paired state');pair_checks+=r['checks']
  c=json.loads((HERE/'cycles'/f'cycle_p{par}_d{d}.json').read_text());require((c['parity'],c['residue'])==(par,d),'cycle filename/state');r=verify_direct(c);cycle_checks+=r['checks']
for d in range(24,73):
 b=copy.deepcopy(pairs[0,d]);b['parity']=1;b['residue']=96-d;b['words'].reverse()
 for w in b['words']:
  for atom in w:
   if atom[0] in('block','tail'):atom[2]=1-atom[2]
 require(b==pairs[1,96-d],'paired exact dual')
base_report=Core(a.base_dir).verify_bases()
controls=[Core(a.base_dir).cycle(N) for N in[2621,2622,128453,128454,6294221]]
print(json.dumps({'status':'passed','optimization_level':sys.flags.optimize,'paired_states':98,'cycle_states':98,'paired_dualities':49,'paired_universal_checks':pair_checks,'direct_universal_checks':cycle_checks,'catalog_and_hash_checks':checks,'base_report':base_report,'cycle_controls':controls,'old_transition_tables_used':False},indent=2))
