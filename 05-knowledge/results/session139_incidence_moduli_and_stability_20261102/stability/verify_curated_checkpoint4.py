#!/usr/bin/env python3
"""Check frozen producer selection and independent evidence pins."""
from pathlib import Path
import hashlib,json
here=Path(__file__).resolve().parent
root=here.parent
manifest=json.loads((here/'CURATED_CHECKPOINT4.json').read_text())
count=0
for section in ['producer_files','independent_evidence']:
    for item in manifest[section]:
        rel=Path(item['path'])
        if rel.is_absolute() or '..' in rel.parts:
            raise RuntimeError(('unsafe path',item['path']))
        path=root/rel
        data=path.read_bytes()
        digest=hashlib.sha256(data).hexdigest()
        if digest!=item['sha256'] or len(data)!=item['bytes']:
            raise RuntimeError(('pin mismatch',item['path'],digest,len(data)))
        count+=1
print(json.dumps({'status':'PASS','pinned_files':count,'scope':manifest['scope']},sort_keys=True,indent=2))
