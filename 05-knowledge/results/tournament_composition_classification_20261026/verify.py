"""Verify frozen packet bytes, then optionally replay without changing them."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf-8', newline='\n')

def need(test, message):
    if not test:
        raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args()
    entries=[]
    for line in (ROOT/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest,name=line.split('  ',1)
        path=ROOT/name
        need(path.is_file(), 'Missing '+name)
        need(hashlib.sha256(path.read_bytes()).hexdigest()==digest, 'Hash mismatch '+name)
        entries.append(name)
    actual={str(p.relative_to(ROOT)).replace('\\','/') for p in ROOT.rglob('*')
            if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS'}
    need(set(entries)==actual, 'Manifest coverage differs')
    print('Manifest: '+str(len(entries))+' files PASS',flush=True)
    if not args.replay:
        return
    runs=[('random_recovery.py','random_sample.json',False),
          ('check_connections.py','connections-normal.json',False),
          ('fibre/check.py','fibre/normal.json',False),
          ('arithmetic/check_arithmetic.py','arithmetic/controls-normal.json',True),
          ('arithmetic/audit_fibre.py','arithmetic/fibre-audit-normal.json',True),
          ('metagraph/check.py','metagraph/normal.json',False)]
    for script,expected,file_output in runs:
        frozen=(ROOT/expected).read_bytes()
        for optimized in (False,True):
            with tempfile.TemporaryDirectory(prefix='tournament-classification-') as temp:
                output=Path(temp)/'receipt.json'
                cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(ROOT/script)]
                if file_output:
                    cmd+=['--out',str(output)]
                result=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=True)
                received=output.read_bytes() if file_output else result.stdout
                need(received==frozen, 'Replay differs: '+script+(' -O' if optimized else ''))
        print(script+': normal and optimized replay PASS',flush=True)
    pairs=[('connections-normal.json','connections-optimized.json'),
           ('fibre/normal.json','fibre/optimized.json'),
           ('arithmetic/controls-normal.json','arithmetic/controls-optimized.json'),
           ('arithmetic/fibre-audit-normal.json','arithmetic/fibre-audit-optimized.json'),
           ('metagraph/normal.json','metagraph/optimized.json')]
    for a,b in pairs:
        need((ROOT/a).read_bytes()==(ROOT/b).read_bytes(),'Frozen receipt pair differs: '+a)
    need((ROOT/'fibre/arithmetic-crossaudit.json').read_bytes()==
         (ROOT/'arithmetic/controls-normal.json').read_bytes(), 'Cross-audit receipt differs')
    print('All frozen receipt pairs and independent cross-audit PASS')

if __name__=='__main__':
    main()
