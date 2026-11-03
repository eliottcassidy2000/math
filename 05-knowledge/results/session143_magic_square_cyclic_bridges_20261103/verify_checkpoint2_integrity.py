#!/usr/bin/env python3
"""Bind independent audit conclusions to the exact published producer files."""
import ast,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent

def need(ok,message):
 if not ok:raise RuntimeError(message)

def checked(path,expected):
 f=root/path;need(f.is_file() and not f.is_symlink(),'missing audited file '+path)
 data=f.read_bytes()
 if 'bytes' in expected:need(len(data)==expected['bytes'],'audited byte length '+path)
 need(hashlib.sha256(data).hexdigest()==expected['sha256'],'audited SHA256 '+path)
 return data
proofs=json.loads((root/'independent/reviewed_proof_pins.json').read_text())['files']
for row in proofs:checked(row['path'],row)
records=json.loads((root/'independent/lane_replay.json').read_text())['records'];manifests=0;lane_files=0
for row in records:
 name=row.get('manifest',row.get('script'));data=checked(row['lane']+'/'+name,row)
 if 'manifest' in row:
  m=json.loads(data);files=m.get('files',m.get('artifacts'))
  need(len(files)==row['files_hashed'],'audited manifest universe')
  for f in files:checked(row['lane']+'/'+f['path'],f);lane_files+=1
  manifests+=1
puzzle=json.loads((root/'independent_puzzle/reviewed_source_pins.json').read_text())
for row in puzzle['files']:checked(puzzle['producer_bundle']+'/'+row['path'],row)
scripts=json.loads((root/'independent/assert_audit.json').read_text())['scripts']
for row in scripts:
 need(not any(isinstance(x,ast.Assert) for x in ast.walk(ast.parse((root/row['path']).read_text()))),'optimization-sensitive assert '+row['path'])
print(json.dumps({'status':'PASS','audited_proofs':len(proofs),'audited_producer_manifests':manifests,'lane_files_bound_to_audit':lane_files,'optimization_safe_certificate_scripts':len(scripts),'audited_puzzle_source_files':len(puzzle['files']),'scope':'Exact content binding of published proof and producer files to independent review, not a new mathematical theorem or external reconstruction.'},indent=2,sort_keys=True))
