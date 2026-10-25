import hashlib,json,pathlib
root=pathlib.Path(__file__).resolve().parent
for name,digest in json.loads((root/'MANIFEST.json').read_text()).items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Hash mismatch: '+name)
print('Manifest verified')
