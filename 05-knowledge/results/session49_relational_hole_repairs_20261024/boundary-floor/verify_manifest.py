#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import argparse,json,shutil,subprocess,sys,tempfile
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--relocate',action='store_true');parser.add_argument('--out',type=Path);args=parser.parse_args()
manifest=json.loads((P/'MANIFEST.json').read_text())
for path,digest in manifest.items():
 if Path(path).is_absolute() or '..' in Path(path).parts:raise RuntimeError('unsafe manifest path')
 if sha256((P/path).read_bytes()).hexdigest()!=digest:raise RuntimeError('hash mismatch '+path)
report={'status':'PASS','manifest_payload_files':len(manifest)}
if args.relocate:
 with tempfile.TemporaryDirectory(prefix='session49-boundary-floor-closed-') as td:
  dest=Path(td)
  for path in manifest:
   q=dest/path;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P/path,q)
  shutil.copyfile(P/'MANIFEST.json',dest/'MANIFEST.json')
  if {str(p.relative_to(dest)) for p in dest.rglob('*') if p.is_file()}!=set(manifest)|{'MANIFEST.json'}:raise RuntimeError('not a manifest-only copy')
  normal=subprocess.run([sys.executable,str(dest/'check.py')],cwd=dest,check=True,capture_output=True,text=True)
  optimized=subprocess.run([sys.executable,'-O',str(dest/'check.py')],cwd=dest,check=True,capture_output=True,text=True)
  if normal.stdout!=optimized.stdout:raise RuntimeError('normal/-O differ')
  report.update({'relocation':'actual manifest-only fresh copy','normal_optimized_equal':True,'closed_check':json.loads(normal.stdout)})
if args.out:args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,sort_keys=True))
