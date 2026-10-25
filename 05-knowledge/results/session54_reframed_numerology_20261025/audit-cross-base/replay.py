from pathlib import Path
import hashlib,json,subprocess,sys
root=Path(__file__).resolve().parent
m=json.loads((root/'manifest.json').read_text())
if {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}!=set(m)|{'manifest.json'}:raise RuntimeError('closed membership')
for name,digest in m.items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('digest '+name)
for opt,fixture in [([], 'check.out'),(['-O'],'check-O.out')]:
    if subprocess.check_output([sys.executable,*opt,str(root/'check_independent.py')],cwd='/tmp')!=(root/fixture).read_bytes():raise RuntimeError('independent fixture '+fixture)
print('PASS audit closure and independent normal/-O exact fixtures; packet.out records producer isolation/provenance run')
