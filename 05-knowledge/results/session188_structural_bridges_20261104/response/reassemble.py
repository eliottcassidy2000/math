"""Restore both audited response receipts from their lossless chunk package."""
import base64,gzip,hashlib,json,pathlib
root=pathlib.Path(__file__).resolve().parent
meta=json.loads((root/'REASSEMBLY.json').read_text())
def checked_path(name):
    if pathlib.Path(name).name!=name: raise RuntimeError('nonlocal package filename')
    return root/name
chunks=[]
for item in meta['parts']:
    raw=checked_path(item['file']).read_bytes()
    if len(raw)!=item['size'] or hashlib.sha256(raw).hexdigest()!=item['sha256']:
        raise RuntimeError('chunk checksum mismatch: '+item['file'])
    chunks.append(b''.join(raw.split()))
data=gzip.decompress(base64.b64decode(b''.join(chunks),validate=True))
for item in meta['outputs']:
    if len(data)!=item['size'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
        raise RuntimeError('receipt checksum mismatch: '+item['file'])
    p=checked_path(item['file'])
    if p.exists() and p.read_bytes()!=data: raise RuntimeError('refusing to replace different existing file: '+item['file'])
for item in meta['outputs']: checked_path(item['file']).write_bytes(data)
print(json.dumps({'status':'PASS','restored':meta['outputs']},indent=2))
