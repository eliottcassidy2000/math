"""Replay exact receipts normally and with -O, and verify the local seal."""
from pathlib import Path
import hashlib,subprocess,sys,json
p=Path(__file__).resolve().parent
listed=set()
for line in (p/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
    digest,name=line.split('  ',1);listed.add(name)
    if hashlib.sha256((p/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(name)
if listed!={f.name for f in p.iterdir() if f.is_file() and f.name!='SHA256SUMS'}:raise RuntimeError('manifest coverage')
print('Manifest PASS',flush=True)
root=p.parents[2]
for item in json.loads((p/'SOURCES.json').read_text(encoding='utf-8'))['verification_inputs']:
    data=(root/item['path']).read_bytes().replace(b'\r\n',b'\n')
    if hashlib.sha256(data).hexdigest()!=item['sha256']:raise RuntimeError('input '+item['path'])
print('Pinned verification inputs PASS',flush=True)
for script in ['check','audit']:
    for flags,label in [([], 'normal'),(['-O'],'optimized')]:
        run=subprocess.run([sys.executable,'-B',*flags,str(p/(script+'.py'))],check=True,capture_output=True,text=True,encoding='utf-8')
        if run.stdout!=(p/(script+'-'+label+'.json')).read_text(encoding='utf-8'):raise RuntimeError('receipt '+script+' '+label)
    print(script+': normal and optimized replay PASS',flush=True)
