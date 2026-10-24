#!/usr/bin/env python3
"""Check complete outer closure, including every nested candidate payload."""
from pathlib import Path, PurePosixPath
from hashlib import sha256
import json
root=Path(__file__).resolve().parent
manifest=root/'AUDIT_MANIFEST.json'
rows=json.loads(manifest.read_text())['files']
listed=set()
for row in rows:
 name=row['path']; rel=PurePosixPath(name)
 if rel.is_absolute() or '..' in rel.parts or name in listed: raise RuntimeError('unsafe or duplicate manifest path '+name)
 listed.add(name); path=root/name
 if not path.is_file() or path.is_symlink(): raise RuntimeError('missing/nonregular payload '+name)
 data=path.read_bytes()
 if len(data)!=row['bytes'] or sha256(data).hexdigest()!=row['sha256']: raise RuntimeError('payload mismatch '+name)
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=manifest}
if actual!=listed: raise RuntimeError('incomplete closure: '+repr({'unlisted':sorted(actual-listed),'absent':sorted(listed-actual)}))
inner=root/'actual-manifest-copy/MANIFEST.json'
if sha256(inner.read_bytes()).hexdigest()!='5cf13f7a96cee50502f81f0e9e74624a6454f5702a76aa72996d5a34075d6102': raise RuntimeError('wrong final producer manifest')
inner_rows=json.loads(inner.read_text())['files']
inner_actual={str(p.relative_to(inner.parent)) for p in inner.parent.rglob('*') if p.is_file()}
if inner_actual!={'MANIFEST.json'}|{r['path'] for r in inner_rows}: raise RuntimeError('inner closure differs from actual producer manifest')
print('PASS:',len(rows),'outer payloads; all 14 final producer payloads; complete recursive file closure')
