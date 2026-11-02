#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
BASE=Path(__file__).resolve().parent
manifest=json.loads((BASE/'MANIFEST.json').read_text())
actual={p.relative_to(BASE).as_posix() for p in BASE.rglob('*') if p.is_file()}
if actual!=set(manifest)|{'MANIFEST.json'}:raise RuntimeError('unlisted or missing audit file')
for p,h in manifest.items():
 if hashlib.sha256((BASE/p).read_bytes()).hexdigest()!=h:raise RuntimeError('audit hash mismatch: '+p)
for p in json.loads((BASE/'SOURCE_PINS.json').read_text())['git_sources']:
 if hashlib.sha256((BASE/p['snapshot']).read_bytes()).hexdigest()!=p['sha256']:raise RuntimeError('source pin mismatch')
for p in json.loads((BASE/'REVIEWED_PRODUCERS.json').read_text()):
 if hashlib.sha256((BASE/p['snapshot']).read_bytes()).hexdigest()!=p['proof_sha256']:raise RuntimeError('reviewed proof mismatch')
with tempfile.TemporaryDirectory() as d:
 for flags,name in [([], 'normal.json'),(['-O'],'optimized.json')]:
  output=subprocess.check_output([sys.executable,'-B']+flags+[str(BASE/'check.py')],cwd=d)
  if output!=(BASE/name).read_bytes():raise RuntimeError('replay mismatch '+name)
print(json.dumps({'status':'PASS','manifest_files':len(manifest),'replays':2,'portable':True},sort_keys=True))
