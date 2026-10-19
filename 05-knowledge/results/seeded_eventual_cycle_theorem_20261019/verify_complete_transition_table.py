#!/usr/bin/env python3
"""Replay the frozen 294-state typed transition table; no discovery dependencies."""
import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
checks=0
def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
manifest=json.loads((ROOT/'certified_pair_transition_manifest.json').read_text())
checker=ROOT/'verify_pair_transition.py'
check(hashlib.sha256(checker.read_bytes()).hexdigest()==manifest['checker_sha256'],'checker checksum')
check(hashlib.sha256((ROOT/manifest['concrete_base_validation_file']).read_bytes()).hexdigest()==manifest['concrete_base_validation_sha256'],'concrete base checksum')
expected={(f,p,d) for f in ['B','Aplus','Aminus'] for p in (0,1) for d in range(24,73)}
seen=set();table={};total=0;split=0
for row in manifest['entries']:
 state=(row['family'],row['parity'],row['residue']);check(state in expected and state not in seen,'state coverage');seen.add(state)
 path=ROOT/'certified-pair-transitions'/row['path'];raw=path.read_bytes();check(hashlib.sha256(raw).hexdigest()==row['sha256'],'witness checksum')
 a=json.loads(raw);check((a['family'],a['parity'],a['residue'])==state,'state metadata');table[state]=a
 kind='single-prefix-split' if a.get('split_prefix_blocks') else 'whole-block'
 check(a['certificate_type']==row['certificate_type']==kind,'certificate type');split+=kind=='single-prefix-split'
 result=subprocess.run([sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(checker),str(path)],capture_output=True,text=True)
 check(result.returncode==0,f'{state}: {result.stderr}')
 report=json.loads(result.stdout);check(report['status']=='passed','individual replay');total+=report['checks']
check(seen==expected,'complete transition table')
for f in ['B','Aplus','Aminus']:
 for d in range(24,73):
  a=table[f,0,d];b=table[f,1,96-d];words=json.loads(json.dumps(list(reversed(a['words']))))
  for word in words:
   for atom in word:
    if atom[0] in ('block','tail'):atom[2]=1-atom[2]
  check(words==b['words'] and a['split_prefix_blocks']==b['split_prefix_blocks'],'exact duality')
print(json.dumps({'status':'passed','states':len(seen),'input_parity0_states':147,'whole_block_states':len(seen)-split,'single_prefix_split_states':split,'table_checks':checks,'individual_checks':total,'manifest_sha256':hashlib.sha256((ROOT/'certified_pair_transition_manifest.json').read_bytes()).hexdigest()},indent=2,sort_keys=True))
