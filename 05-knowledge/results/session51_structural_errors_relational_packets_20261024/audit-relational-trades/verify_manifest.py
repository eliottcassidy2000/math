#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import argparse, json, shutil, subprocess, sys, tempfile
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--relocate',action='store_true');args=parser.parse_args()
files={p.relative_to(P).as_posix():{'bytes':p.stat().st_size,'sha256':sha256(p.read_bytes()).hexdigest()}
       for p in sorted(P.rglob('*')) if p.is_file() and p.relative_to(P).as_posix()!='MANIFEST.json' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
target=P/'MANIFEST.json'
if args.write:target.write_text(json.dumps(files,sort_keys=True,indent=2)+'\n')
elif files!=json.loads(target.read_text()):raise RuntimeError('closed manifest mismatch')
print(json.dumps({'status':'PASS','manifest_files':len(files),'manifest_sha256':sha256(target.read_bytes()).hexdigest()},sort_keys=True),flush=True)
if args.relocate:
    with tempfile.TemporaryDirectory(prefix='session51-independent-relational-audit-') as directory:
        destination=Path(directory)
        for name in list(files)+['MANIFEST.json']:
            file=destination/name;file.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P/name,file)
        subprocess.run([sys.executable,'-B',str(destination/'reproduce.py')],cwd=destination,check=True)
    print(json.dumps({'status':'PASS','manifest_only_relocation':True,'normal_and_optimized_full_suite':True},sort_keys=True))
