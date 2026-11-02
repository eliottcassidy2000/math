#!/usr/bin/env python3
"""Verify the compact n>=68 selection; checks remain active under -O."""
import argparse
import hashlib
import json
from pathlib import Path

here=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--manifest',type=Path,default=here/'CURATED_CHECKPOINT5.json')
p.add_argument('--root',type=Path,default=here.parent)
args=p.parse_args()
m=json.loads(args.manifest.read_text())
seen=set()
for item in m['files']:
    rel=Path(item['path'])
    if rel.is_absolute() or '..' in rel.parts or str(rel) in seen:
        raise RuntimeError(('unsafe or duplicate selection',item['path']))
    seen.add(str(rel))
    data=(args.root/rel).read_bytes()
    digest=hashlib.sha256(data).hexdigest()
    if digest!=item['sha256'] or len(data)!=item['bytes']:
        raise RuntimeError(('pin mismatch',item['path'],digest,len(data)))
for first,second,status in m['receipt_pairs']:
    one=(args.root/first).read_bytes()
    two=(args.root/second).read_bytes()
    if one!=two or json.loads(one)['status']!=status:
        raise RuntimeError(('receipt mode/status mismatch',first,second))
print(json.dumps({'status':'PASS','pinned_files':len(seen),
                  'concurrent_source_pins_recorded':len(m['concurrent_sources']),
                  'scope':m['scope']},sort_keys=True,indent=2))
