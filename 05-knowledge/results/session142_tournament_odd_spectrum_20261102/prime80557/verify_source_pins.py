#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
records=json.loads((P/'sources.json').read_text())
for r in records:
 data=(P/r['local_snapshot']).read_bytes()
 if len(data)!=r['bytes'] or hashlib.sha256(data).hexdigest()!=r['sha256'] or hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()!=r['blob']:
  raise ValueError(r['local_snapshot'])
print(json.dumps({'status':'passed','source_count':len(records),'source_bytes':sum(r['bytes'] for r in records)},sort_keys=True))
