#!/usr/bin/env python3
"""Portable fail-closed audit replay, including manifest-only relocation controls."""
import pathlib,json,hashlib,subprocess,sys,tempfile,os
ROOT=pathlib.Path(__file__).resolve().parent
PIN='0855158ced20c231100eeb5d8faf9094337e2a28dd26ad74befb437997814f26'
def ck(p,m):
 if not p:raise RuntimeError(m)
def validate(root):
 m=json.loads((root/'MANIFEST.json').read_text());files=m['files']
 for name in files:
  p=pathlib.PurePosixPath(name)
  ck(not p.is_absolute() and '..' not in p.parts and str(p)==name,'nonrelative path')
 actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
 ck(actual==set(files)|{'MANIFEST.json'},'recursive file-set mismatch')
 for name,digest in files.items():
  ck(hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,'digest mismatch: '+name)
 ck(hashlib.sha256((root/'producer/MANIFEST.json').read_bytes()).hexdigest()==PIN,'producer pin')
 return m

def main():
 m=validate(ROOT);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 for script,expected in [('audit.py','independent.json'),('closure.py','closure.json'),('producer/replay.py','replay.log')]:
  for flags in (['-B'],['-B','-O']):
   args=[sys.executable,*flags,str(ROOT/script)]+(['--self-test'] if script.startswith('producer/') else [])
   p=subprocess.run(args,cwd='/',env=env,capture_output=True,check=True)
   ck(p.stdout==(ROOT/expected).read_bytes(),'replay mismatch: '+script)
 if '--self-test' in sys.argv:
  with tempfile.TemporaryDirectory() as td:
   q=pathlib.Path(td)/'relocated';q.mkdir()
   for name in ['MANIFEST.json',*m['files']]:
    dest=q/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((ROOT/name).read_bytes())
   # Actually execute the isolated manifest-only copy, not just its validator.
   for flags in (['-B'],['-B','-O']):
    p=subprocess.run([sys.executable,*flags,str(q/'replay.py')],cwd='/',env=env,capture_output=True,check=True)
    ck(p.stdout==b'PASS: portable audit replay and recursive manifest closure\n','relocation output')
   target=q/'independent.json';original=target.read_bytes()
   for hostile in ('corrupt','missing','nested-extra'):
    if hostile=='corrupt':target.write_bytes(original+b' ')
    if hostile=='missing':target.unlink()
    if hostile=='nested-extra':
     (q/'unexpected').mkdir();(q/'unexpected'/'extra').write_text('hostile')
    for flags in (['-B'],['-B','-O']):
     p=subprocess.run([sys.executable,*flags,str(q/'replay.py')],cwd='/',env=env,capture_output=True)
     ck(p.returncode!=0 and (b'digest mismatch' in p.stderr or b'file-set mismatch' in p.stderr),'hostile accepted: '+hostile)
    if hostile in ('corrupt','missing'):target.write_bytes(original)
    else:(q/'unexpected'/'extra').unlink();(q/'unexpected').rmdir()
   validate(q)
  print('PASS: manifest-only relocation normal/-O; corrupt/missing/nested-extra rejected normal/-O')
 else:print('PASS: portable audit replay and recursive manifest closure')
if __name__=='__main__':main()
