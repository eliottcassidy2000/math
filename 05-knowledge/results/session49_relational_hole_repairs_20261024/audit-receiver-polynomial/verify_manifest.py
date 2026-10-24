#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import argparse,json,shutil,tempfile,subprocess,sys
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');ap.add_argument('--relocate',action='store_true');args=ap.parse_args()
files={p.relative_to(P).as_posix():{'bytes':p.stat().st_size,'sha256':sha256(p.read_bytes()).hexdigest()} for p in sorted(P.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
target=P/'MANIFEST.json'
if args.write:target.write_text(json.dumps(files,sort_keys=True,indent=2)+'\n')
else:
    expected=json.loads(target.read_text())
    if files!=expected:raise RuntimeError('Manifest mismatch: '+','.join(p for p in sorted(set(files)|set(expected)) if files.get(p)!=expected.get(p)))
print(json.dumps({'status':'PASS','manifest_files':len(files),'manifest_sha256':sha256(target.read_bytes()).hexdigest()},sort_keys=True),flush=True)
if args.relocate:
    with tempfile.TemporaryDirectory(prefix='receiver-threshold-independent-manifest-') as d:
        q=Path(d)
        for name in list(files)+['MANIFEST.json']:
            dst=q/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P/name,dst)
        subprocess.run([sys.executable,'-B',str(q/'reproduce.py')],cwd=q,check=True)
    print(json.dumps({'status':'PASS','manifest_only_relocation':True,'normal_and_optimized_full_suite':True},sort_keys=True))
