import pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent
entries=json.loads((p/'MANIFEST.json').read_text())
actual={str(f.relative_to(p)) for f in p.rglob('*') if f.is_file()}
expected={x['path'] for x in entries}|{'MANIFEST.json'}
if actual!=expected:raise RuntimeError('closed inventory mismatch')
for x in entries:
 if hashlib.sha256((p/x['path']).read_bytes()).hexdigest()!=x['sha256']:raise RuntimeError('pin mismatch '+x['path'])
print('PASS: closed inventory and hashes')
