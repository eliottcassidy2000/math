import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    args=parser.parse_args();root=Path(__file__).resolve().parent
    data=json.loads((root/'MANIFEST.json').read_text())
    for entry in data['files']:
        relative=Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise RuntimeError('unsafe manifest path')
        path=root/relative
        if path.is_symlink() or not path.is_file():raise RuntimeError('not regular file')
        raw=path.read_bytes();raw.decode('utf-8')
        if len(raw)!=entry['bytes'] or hashlib.sha256(raw).hexdigest()!=entry['sha256']:
            raise RuntimeError('content mismatch: '+str(relative))
    expected={f['path'] for f in data['files']}|{'MANIFEST.json'}
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() or p.is_symlink()}
    if actual!=expected:raise RuntimeError('extra or missing payloads')
    results=[]
    if args.full:
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)
            for name in expected:
                dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(root/name,dest)
            for flags in ([],['-O']):
                raw=subprocess.check_output([sys.executable,'-B',*flags,'check.py'],cwd=target,text=True)
                results.append(json.loads(raw))
            if results[0]!=results[1]:raise RuntimeError('optimization mismatch')
            subprocess.run([sys.executable,'-B','verify_packet.py'],cwd=target,check=True,stdout=subprocess.DEVNULL)
    print(json.dumps({'status':'PASS','files':len(data['files']),'replays':results},sort_keys=True))

if __name__=='__main__':main()
