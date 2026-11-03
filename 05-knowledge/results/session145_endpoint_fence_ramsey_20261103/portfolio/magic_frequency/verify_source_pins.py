#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
pins=json.loads((root/'SOURCE_PINS.json').read_text())
for pin in pins['sources']:
    data=(root/pin['snapshot']).read_bytes()
    actual={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'git_blob_sha1':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()}
    for key,value in actual.items():
        if pin[key]!=value:raise ValueError('source pin mismatch: '+pin['snapshot']+' '+key)
print(json.dumps({'status':'PASS','revision':pins['revision'],'source_count':len(pins['sources'])},sort_keys=True))
