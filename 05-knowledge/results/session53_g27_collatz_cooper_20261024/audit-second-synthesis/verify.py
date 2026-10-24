from pathlib import Path
import hashlib,json,subprocess,sys
p=Path(__file__).resolve().parent
m=json.loads((p/'MANIFEST.json').read_text())
actual={f.relative_to(p).as_posix() for f in p.rglob('*') if f.is_file()}
if actual!=set(m)|{'MANIFEST.json'}:raise RuntimeError('closed recursive membership mismatch')
for f,h in m.items():
 if hashlib.sha256((p/f).read_bytes()).hexdigest()!=h:raise RuntimeError('hash mismatch '+f)
for mode,fixture in [([], 'normal.json'),(['-O'],'optimized.json')]:
 out=subprocess.check_output([sys.executable,*mode,str(p/'check.py')],cwd='/tmp')
 if out!=(p/fixture).read_bytes():raise RuntimeError('fixture mismatch '+fixture)
print('PASS: closed source-backed summary audit; both replay modes match')
