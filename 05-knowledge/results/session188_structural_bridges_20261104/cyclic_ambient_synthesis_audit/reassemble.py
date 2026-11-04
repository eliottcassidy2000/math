"""Restore complete independent cyclic ambient synthesis records losslessly."""
import base64,gzip,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
meta=json.loads((root/'REASSEMBLY.json').read_text())
def local(name):
    if Path(name).name!=name:raise RuntimeError('nonlocal package filename')
    return root/name
ready=[];seen=set()
for group in meta['groups']:
    chunks=[]
    for item in group['parts']:
        raw=local(item['file']).read_bytes()
        if len(raw)!=item['size'] or hashlib.sha256(raw).hexdigest()!=item['sha256']:raise RuntimeError('chunk hash mismatch '+item['file'])
        chunks.append(b''.join(raw.split()))
    data=gzip.decompress(base64.b64decode(b''.join(chunks),validate=True))
    for item in group['outputs']:
        path=local(item['file'])
        if item['file'] in seen:raise RuntimeError('duplicate output filename')
        seen.add(item['file'])
        if len(data)!=item['size'] or hashlib.sha256(data).hexdigest()!=item['sha256']:raise RuntimeError('restored hash mismatch '+item['file'])
        if path.exists() and path.read_bytes()!=data:raise RuntimeError('refusing to replace different existing output '+item['file'])
        ready.append((path,data))
# No output is written until all groups and existing outputs have passed.
for path,data in ready:path.write_bytes(data)
print(json.dumps({'status':'PASS','outputs':[{'file':p.name,'size':len(b),'sha256':hashlib.sha256(b).hexdigest()} for p,b in ready]},indent=2))
