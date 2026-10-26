"""Verify frozen files and replay three exact paths in normal and optimized mode."""
import hashlib
from pathlib import Path
import subprocess
import sys

p=Path(__file__).resolve().parent
listed=set()
for line in (p/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
    expected,name=line.split('  ',1);listed.add(name)
    if hashlib.sha256((p/name).read_bytes()).hexdigest()!=expected:raise RuntimeError('hash mismatch: '+name)
actual={f.name for f in p.iterdir() if f.is_file() and f.name!='SHA256SUMS'}
if actual!=listed:raise RuntimeError('manifest coverage mismatch')
print('Manifest: all '+str(len(listed))+' files match',flush=True)
for script,names in [('check.py',['normal.json','optimized.json']),
                     ('audit.py',['audit-normal.json','audit-optimized.json']),
                     ('arithmetic.py',['arithmetic-normal.json','arithmetic-optimized.json'])]:
    for flags,name in [([],names[0]),(['-O'],names[1])]:
        run=subprocess.run([sys.executable,'-B',*flags,str(p/script)],check=True,capture_output=True,text=True,encoding='utf-8')
        if run.stdout.replace('\r\n','\n')!=(p/name).read_text(encoding='utf-8'):raise RuntimeError('replay differs: '+name)
    if (p/names[0]).read_bytes()!=(p/names[1]).read_bytes():raise RuntimeError('normal/optimized mismatch')
    print(script+': frozen normal and optimized replay PASS',flush=True)
