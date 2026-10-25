from pathlib import Path
import hashlib,json,subprocess,sys,tempfile,shutil
R=Path(__file__).resolve().parent
m=json.loads((R/'manifest.json').read_text())
if {str(p.relative_to(R)) for p in R.rglob('*') if p.is_file()}!=set(m)|{'manifest.json'}:raise RuntimeError('audit exact membership')
for name,digest in m.items():
    if hashlib.sha256((R/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('audit digest '+name)
if hashlib.sha256((R/'frozen-target/manifest.json').read_bytes()).hexdigest()!='5d9519992f6ff49e9e8765a7d50928634179f890492586a60fa8818aa76b841d':raise RuntimeError('target pin')
for opt,fixture in [([], 'check.out'),(['-O'],'check-O.out')]:
    got=subprocess.check_output([sys.executable,*opt,str(R/'check_independent.py')],cwd='/tmp')
    if got!=(R/fixture).read_bytes():raise RuntimeError('independent fixture '+fixture)
with tempfile.TemporaryDirectory() as t:
    dst=Path(t)/'packet';shutil.copytree(R/'frozen-target',dst)
    lines=[]
    for opt in ([],['-O']):
        for program in ('replay.py','hostile.py'):
            lines.append(('normal' if not opt else '-O')+' '+program+': '+subprocess.check_output([sys.executable,*opt,str(dst/program)],cwd=t).decode().strip())
    lines.append('PASS actual isolated manifest-only normal/-O replay and hostile rejection')
    if '\n'.join(lines)+'\n'!=(R/'packet.out').read_text():raise RuntimeError('packet fixture')
print('PASS audit closure, independent normal/-O, isolated producer replay and hostiles')
