#!/usr/bin/env python3
"""Validate this audit; optionally compare all reviewed publication sources."""
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--root-dir',type=Path)
args=parser.parse_args();count=0
for name,expected in json.loads((HERE/'MANIFEST.json').read_text())['files'].items():
 p=(HERE/name).resolve()
 if not p.is_relative_to(HERE):raise RuntimeError('Manifest escapes audit')
 if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise RuntimeError('Hash mismatch: '+name)
 count+=1
reviewed=json.loads((HERE/'REVIEWED_INPUTS.json').read_text())['documents']
for row in reviewed:
 if 'snapshot' in row:
  b=(HERE/row['snapshot']).read_bytes()
  if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError('Snapshot mismatch')
if args.root_dir:
 root=args.root_dir.resolve()
 aliases={'session47-root/final/':'','session47-helper-resilience/final-v3/':'helper-resilience-v3/','session47-audit-helper-resilience/third-audit/':'audit-helper-resilience-v3/','session47-audit-suspension-obstruction/':'audit-five-to-seven/','session47-audit-fresh570/':'fresh570-audit/'}
 for row in reviewed:
  name=row['source_path'];p=root/name
  if not p.exists():
   translated=name.removeprefix('session47-')
   for old,new in aliases.items():
    if name.startswith(old):translated=new+name[len(old):];break
   p=root/translated
  b=p.read_bytes()
  if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError('Reviewed source changed: '+str(p))
  count+=1
if (HERE/'normal.json').read_bytes()!=(HERE/'optimized.json').read_bytes():raise RuntimeError('Control modes differ')
if (HERE/'manifests.normal.json').read_bytes()!=(HERE/'manifests.optimized.json').read_bytes():raise RuntimeError('Manifest modes differ')
print(json.dumps({'status':'PASS','hashes_checked':count,'reviewed_external_sources_checked':len(reviewed) if args.root_dir else 0},sort_keys=True))
