#!/usr/bin/env python3
"""Independent synthesis consistency checks, not a replacement for packet proofs."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent
note=(ROOT/'NOTE.md').read_text()
links=re.findall(r'\]\(([^)]+)\)',note)
local=[];external=[]
for target in links:
 if target.startswith('../'):
  p=(ROOT/target.split('#')[0]).resolve()
  if not p.is_file():raise RuntimeError('Missing sibling target: '+target)
  local.append(target)
 elif target.startswith('https://'):
  external.append(target)
 else:raise RuntimeError('Unexpected link form: '+target)
pairs=[]
for folder in ['gates','circuits','information','spectrum','resource']:
 directory=ROOT.parent/folder
 optimized=sorted(directory.glob('*.optimized.json'))
 if not optimized:raise RuntimeError('No paired receipts for '+folder)
 for op in optimized:
  normal=op.with_name(op.name.replace('.optimized.json','.json'))
  if folder=='circuits' and op.name=='reviewer_check.optimized.json':normal=directory/'certificate.json'
  if json.loads(normal.read_text())!=json.loads(op.read_text()):raise RuntimeError('Receipt mismatch '+str(op))
  pairs.append({'normal':folder+'/'+normal.name,'optimized':folder+'/'+op.name})
expected=['at least two bits on each side','H=7 and H=21 are globally impossible','not an n-bit lower bound','All those bases are isomorphic','O((n+1)k!(k+1)3^(k+1))','No general B2 lower bound','sum_j L^(q_j) <= K L^d']
for phrase in expected:
 if phrase not in note:raise RuntimeError('Scope/correction marker missing: '+phrase)
print(json.dumps({'status':'PASS','local_link_occurrences':len(local),'distinct_local_targets':len(set(local)),'external_source_links_recorded_not_network_verified':external,'matching_receipt_pairs':pairs,'receipt_pair_count':len(pairs),'scope':'Five-packet links, complete paired JSON receipt consistency, and key scope markers. Mathematical review is documented separately.'},indent=2,sort_keys=True))
