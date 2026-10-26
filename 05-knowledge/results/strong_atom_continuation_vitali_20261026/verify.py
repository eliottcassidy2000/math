"""Check the delivery manifest, frozen response atlases, and all exact receipts."""
from pathlib import Path
import hashlib,subprocess,sys
p=Path(__file__).resolve().parent
listed=set()
for line in (p/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
    digest,name=line.split('  ',1);listed.add(name)
    if hashlib.sha256((p/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('hash mismatch: '+name)
actual={f.name for f in p.iterdir() if f.is_file() and f.name!='SHA256SUMS'}
if actual!=listed:raise RuntimeError('manifest coverage')
print('Manifest: '+str(len(listed))+' files match',flush=True)
for script in ['atlas','controls','audit']:
    for flags,label in [([], 'normal'),(['-O'],'optimized')]:
        args=[sys.executable,'-B',*flags,str(p/(script+'.py'))]
        if script=='atlas':args+=['--verify-atlases']
        run=subprocess.run(args,check=True,capture_output=True,text=True,encoding='utf-8')
        if run.stdout.replace('\r\n','\n')!=(p/(script+'-'+label+'.json')).read_text(encoding='utf-8'):
            raise RuntimeError('frozen replay mismatch: '+script+' '+label)
    if (p/(script+'-normal.json')).read_bytes()!=(p/(script+'-optimized.json')).read_bytes():
        raise RuntimeError('normal/optimized bytes differ: '+script)
    print(script+': normal/optimized frozen replay PASS',flush=True)
