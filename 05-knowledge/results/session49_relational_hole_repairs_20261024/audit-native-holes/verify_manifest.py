#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import argparse,json,shutil,subprocess,sys,tempfile
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--relocate',action='store_true');parser.add_argument('--out',type=Path);args=parser.parse_args()
manifest=json.loads((P/'MANIFEST.json').read_text())
for name,digest in manifest.items():
    rel=Path(name)
    if rel.is_absolute() or '..' in rel.parts:raise RuntimeError('unsafe manifest path')
    if sha256((P/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('manifest mismatch: '+name)
report={'status':'PASS','manifest_payload_files':len(manifest)}
if args.relocate:
    with tempfile.TemporaryDirectory(prefix='session49-native-audit-closed-') as temporary:
        dest=Path(temporary)
        for name in manifest:
            target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P/name,target)
        shutil.copyfile(P/'MANIFEST.json',dest/'MANIFEST.json')
        actual={str(p.relative_to(dest)) for p in dest.rglob('*') if p.is_file()}
        if actual!=set(manifest)|{'MANIFEST.json'}:raise RuntimeError('copy was not manifest-only')
        check=subprocess.run([sys.executable,str(dest/'verify_manifest.py')],cwd=dest,check=True,capture_output=True,text=True)
        replay=subprocess.run([sys.executable,str(dest/'reproduce.py')],cwd=dest,check=True,capture_output=True,text=True)
        report['relocation']='actual fresh manifest-only copy';report['closed_manifest']=json.loads(check.stdout);report['closed_reproduction']=json.loads(replay.stdout)
if args.out:args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,sort_keys=True))
